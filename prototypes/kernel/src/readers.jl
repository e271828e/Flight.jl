# The read-selector family and the compiled reader (§14.4, §14.7): the closed
# set of deferred reads every reader of the model addresses through, the
# declared read set they are labeled in, and the compiled gather that runs one
# such set against an executor.
#
# The reader is `apply!`'s **gather twin** (§14.4): one primitive family run in
# both directions over the same layout tables. What the condition algebra
# resolves to a plan of baked destinations, a read set resolves to a tuple of
# baked sources — the same resolve-once/execute-many shape, the same §13.1
# collecting register on the way in, and the same rule that all string work is
# a function of the *shape* and none of it survives into the read.
#
# This file sits above the data plane because the selectors are its vocabulary
# too: an output binding's `reads` (§11.2, bindings.jl) names the three table
# members of the family declared here. What it needs from the condition algebra
# — the `x`-offset walk and the "is this path a level of the build at all"
# predicate — it calls at resolution time, which is long after conditions.jl
# has been read.

# --- the selector family (§14.4), closed --------------------------------------

"""
§14.4's read-selector family, closed: `get_state(path, field[, i])`,
`get_deriv(path, field[, i])`, `get_output(path, field[, i])`,
`get_input(face)` and `get_face(name)` — one address space for every reader of
the model.

The names carry a deliberate `get_` prefix. A selector is a *deferred read*: a
value describing the read the compiled gather will perform, inert until it is
resolved against a source. The prefix names that action, and it keeps five
short common nouns out of the namespace user declarations share with domain
code.

The family splits by **source**, not by client (§14.4's source rule). The
store selectors `get_state`/`get_deriv` resolve only against live stores, which
here means an `Executor` — the one live-store holder this prototype has; the
table selectors `get_output`/`get_input`/`get_face` resolve against a table
source, an executor's own signal table or a published snapshot. A
snapshot-bound reader naming a store selector is therefore refused at attach
(`ReadBindingUnresolved`, bindings.jl), by source rather than by client.

`get_output` is the *inspection* register — a component's own declared output
port, addressed by path — and `get_face` the *integration* register: a
root-exported output face, named, curated, meaning-stable under substitution
(§11.2). `get_input` reads a root input back, the source cell it is. Only cells
and stores are addressable: there is no selector for a value a component
computes without declaring it, and the remedy is the same at every register —
the component exports it (§5.2, §8.3).

`i` is the optional component index (§14.10): the read is `v[i]`, so a vector
leaf yields named scalars. Absent, the whole value is read.
"""
struct GetState
    path::String
    field::Symbol
    i::Union{Nothing,Int}
end

struct GetDeriv
    path::String
    field::Symbol
    i::Union{Nothing,Int}
end

struct GetOutput
    path::String
    name::Symbol
    i::Union{Nothing,Int}
end

struct GetInput
    face::Symbol
end

struct GetFace
    name::Symbol
end

const ReadSelector = Union{GetState,GetDeriv,GetOutput,GetInput,GetFace}
const StoreSelector = Union{GetState,GetDeriv}

_index_arg(::Nothing) = nothing
_index_arg(i::Integer) = Int(i)
_index_arg(i) = throw(BuildError(
    ArgumentInvalid(call = :selector, reason = :index_not_integer, value = i)))

get_state(path::AbstractString, field::Union{Symbol,AbstractString}, i = nothing) =
    GetState(String(path), Symbol(field), _index_arg(i))
get_deriv(path::AbstractString, field::Union{Symbol,AbstractString}, i = nothing) =
    GetDeriv(String(path), Symbol(field), _index_arg(i))
get_output(path::AbstractString, name::Union{Symbol,AbstractString}, i = nothing) =
    GetOutput(String(path), Symbol(name), _index_arg(i))
get_input(face::Union{Symbol,AbstractString}) = GetInput(Symbol(face))
get_face(name::Union{Symbol,AbstractString}) = GetFace(Symbol(name))

# The selector as authored, for the diagnostics: a refusal names the read the
# way its author wrote it, which is what makes a collected list readable.
_ipart(i) = i === nothing ? "" : ", $i"
_spell(s::GetState) = "get_state(\"$(s.path)\", :$(s.field)$(_ipart(s.i)))"
_spell(s::GetDeriv) = "get_deriv(\"$(s.path)\", :$(s.field)$(_ipart(s.i)))"
_spell(s::GetOutput) = "get_output(\"$(s.path)\", :$(s.name)$(_ipart(s.i)))"
_spell(s::GetInput) = "get_input(:$(s.face))"
_spell(s::GetFace) = "get_face(:$(s.name))"

# --- the declared read set (§14.7) ---------------------------------------------

"""
The declared read set (§14.7): the labeled selectors in one type, so that a
bare NamedTuple of selectors reaching a service is refused with a directive
rather than a `MethodError`, exactly as `combine` refuses one in the condition
algebra (§14.2).
"""
struct Reads{NT<:NamedTuple}
    sels::NT
end

"""
    reads(; label = selector, …)

The read set a service declares: the labels are the names the gathered
NamedTuple carries, and the names a trim problem's residual function
destructures (§14.7). Order is the declared side's own, as everywhere at an
author↔framework NamedTuple seam (§9.5).
"""
reads(; sels...) = _reads(NamedTuple(sels))

function _reads(nt::NamedTuple)
    for (label, s) in pairs(nt)
        s isa ReadSelector || throw(BuildError(ReadSetMisuse(
            observed = typeof(s), reason = :not_a_selector, label = label,
            in_hand = Symbol[nameof(typeof(v)) for v in values(nt) if v isa ReadSelector])))
    end
    Reads(nt)
end

# --- the compiled reader (§14.4) ------------------------------------------------
# One entry per selector, its leaf type and its index in the entry's *type*, so
# the gather is a tuple walk the compiler unrolls: no dictionary, no address
# arithmetic and no branch survives resolution. The four entry kinds are the
# four homes a read can come from — the flat state buffer, the derivative
# buffer beside it, a discrete component's own store, and the signal table.

struct StateRead{P,I}
    off::Int
    i::I
end

struct DerivRead{P,I}
    off::Int
    i::I
end

struct StoreRead{S,F,I}
    ci::Int
    i::I
end

struct CellRead{A,I}
    addr::A
    i::I
end

_take(v, ::Nothing) = v
_take(v, i::Int) = v[i]

@inline _read(r::StateRead{P}, ex::Executor) where {P} =
    _take(reconstruct(P, ex.xbuf, r.off), r.i)
@inline _read(r::DerivRead{P}, ex::Executor) where {P} =
    _take(reconstruct(P, ex.ẋbuf, r.off), r.i)
# The `s` stores are held by component index in a `Vector{Any}` — one store
# type per component type, not per model — so the baked store type is what
# keeps the read inferable. The assertion goes on the *reference*: asserting
# the dereferenced value instead leaves the `[]` a dynamic call, which boxes.
@inline _read(r::StoreRead{S,F}, ex::Executor) where {S,F} =
    _take(getfield((ex.sstores[r.ci]::Base.RefValue{S})[], F), r.i)
@inline _read(r::CellRead, ex::Executor) = _take(gather(ex.store, r.addr), r.i)

"""
One compiled read set (§14.4): the labels as a type parameter, the resolved
entries as a tuple carrying their leaf types. `gather(reader, executor)` is the
gather twin of `apply!` — a stack-only NamedTuple per evaluation, type-stable
and allocation-free for scalar and `SVector` leaves, which is what lets a
service's per-iteration read cost nothing beyond the sweep it follows.

A reader is compiled at one activation and resolves against an executor at that
activation: store selectors reach live stores (§14.4's source rule), and the
executor is the only live-store holder here. That activation is `T`, the first
type parameter — the entries bake offsets, store types and cell addresses that
are one activation's, so the pairing is dispatch and a mismatch is the
internal-invariant refusal build.jl raises (§14.4) rather than a silent read of
another cell's slot.
"""
struct Reader{T,L,E<:Tuple}
    entries::E
end

Reader{T,L}(entries::E) where {T,L,E<:Tuple} = Reader{T,L,E}(entries)

@inline gather(r::Reader{T,L}, ex::Executor{T}) where {T,L} =
    NamedTuple{L}(map(e -> _read(e, ex), r.entries))

gather(::Reader{T}, ::Executor{S}) where {T,S} = _activation_mismatch("reader", T, S)

# --- resolution (§14.4, §13.1) --------------------------------------------------

"""
    _compile_reads(rs::Reads, b::Build, T = Float64) → Reader

Resolve a declared read set against a build and compile it, validating every
selector in §13.1's collecting register — full list, violations collected, one
`BuildError`. Schema is the authority on *may you read this, at what type*, and
the activation's layout supplies the source: an `xbuf` offset for a continuous
state field, the `ẋbuf` offset beside it for its derivative, a component index
for a discrete `s`, a cell address for a port, a root input or a root-exported
face.

The standalone entry point is internal, because a read set is only ever
compiled inside a client: `trim!` splices the collected list into its own throw
beside the problem's `TrimProblemInvalid` values (§14.8), and a device binding
raises `ReadBindingUnresolved` at attach instead (§11.2). That is why the
collecting pass is factored apart from the throw. The violations are
`TapResolution` values either way — the kind is the read's own, and it is the
*site* that decides where they surface (§13.2, Appendix C).
"""
function _compile_reads(rs::Reads, b::Build, ::Type{T} = Float64) where {T}
    reader, viol = _resolve_reads(rs, b, T)
    isempty(viol) || throw(BuildError(viol))
    reader
end

# A bare NamedTuple of selectors is the §14.2 misuse in the read register: the
# same slip, the same directive, and not a `MethodError`.
_compile_reads(other, ::Build, ::Type = Float64) = throw(BuildError(
    ReadSetMisuse(observed = typeof(other), reason = :not_a_read_set)))

"""
The collecting half of `_compile_reads`, shared with the services that own their
own setup diagnostic (§14.8): returns the compiled reader and the violation
list, the reader being `nothing` when anything failed.
"""
function _resolve_reads(rs::Reads, b::Build, ::Type{T}) where {T}
    act = activation(b, T)
    offs = _x_offsets(act.decls, b.tiers)
    viol = Diagnostic[]
    entries = Any[]
    for (label, s) in pairs(rs.sels)
        e = _resolve_selector(s, label, b, act, offs, viol)
        e === nothing || push!(entries, e)
    end
    (isempty(viol) ? Reader{T,keys(rs.sels)}(Tuple(entries)) : nothing, viol)
end

# The component a path-addressed selector names. The treatment is
# `_component`'s, one register over: the offender named plainly, an assembly
# discriminated from a path that is nothing at all (candidate lists are absent
# here, README).
function _read_component(s, label::Symbol, flat::Flat, viol::Vector{Diagnostic})
    i = findfirst(==(s.path), flat.paths)
    i === nothing || return i
    # An empty path names the root, which is a level of every build — the
    # prefix test cannot see that, the root's segment being no segment at all.
    push!(viol, _rviol(label, s, (isempty(s.path) || _addresses_level(flat, s.path)) ?
                          :assembly_path : :unknown_path))
    nothing
end

# One `TapResolution` off a selector: the label and the selector as authored are
# what makes a collected list readable, and the tap set, path and index come off
# the selector's own kind (§14.10's payload); each arm adds what it observed.
_rviol(label::Symbol, s, reason::Symbol; kw...) =
    TapResolution(; label = label, selector = _spell(s), reason = reason, tap = _tap(s),
                  path = _selpath(s), index = _selindex(s), kw...)

_tap(::Union{GetState,GetDeriv}) = :x
_tap(::Union{GetOutput,GetFace}) = :y
_tap(::GetInput) = :u

_selpath(s::Union{GetState,GetDeriv,GetOutput}) = s.path
_selpath(::Union{GetInput,GetFace}) = ""
_selindex(s::Union{GetState,GetDeriv,GetOutput}) = s.i
_selindex(::Union{GetInput,GetFace}) = nothing

# `i` is checked against the resolved leaf's declared type in exactly one
# respect: `getindex` has to mean something there. A scalar leaf is refused;
# nothing further is checked, the index being the author's own coordinate
# choice over a value whose length the schema does not fix everywhere.
function _check_index(s, label::Symbol, ::Type{P}, viol::Vector{Diagnostic}) where {P}
    (s.i === nothing || !(P <: Real)) && return true
    push!(viol, _rviol(label, s, :scalar_index; declared = P))
    false
end

_declares(label::Symbol, s, declares::Symbol, declared::NamedTuple) =
    _rviol(label, s, :undeclared; declares = declares, field = _field(s),
           candidates = collect(keys(declared)))

_field(s::Union{GetState,GetDeriv}) = s.field
_field(s::GetOutput) = s.name

function _resolve_selector(s::GetState, label::Symbol, b::Build, act::Activation,
                       offs::Vector{Int}, viol::Vector{Diagnostic})
    ci = _read_component(s, label, b.flat, viol)
    ci === nothing && return nothing
    d, t = act.decls[ci], b.tiers[ci]
    declared = state_decls(d, t)
    haskey(declared, s.field) ||
        (push!(viol, _declares(label, s, :state_field, declared)); return nothing)
    P = typeof(declared[s.field])
    _check_index(s, label, P, viol) || return nothing
    t === CONTINUOUS ?
        StateRead{P,typeof(s.i)}(offs[ci] + _leaf_offset(d.x, s.field), s.i) :
        StoreRead{typeof(d.s),s.field,typeof(s.i)}(ci, s.i)
end

function _resolve_selector(s::GetDeriv, label::Symbol, b::Build, act::Activation,
                       offs::Vector{Int}, viol::Vector{Diagnostic})
    ci = _read_component(s, label, b.flat, viol)
    ci === nothing && return nothing
    d, t = act.decls[ci], b.tiers[ci]
    if t !== CONTINUOUS
        push!(viol, _rviol(label, s, :discrete_deriv; field = s.field))
        return nothing
    end
    haskey(d.x, s.field) ||
        (push!(viol, _declares(label, s, :state_field, d.x)); return nothing)
    P = typeof(d.x[s.field])
    _check_index(s, label, P, viol) || return nothing
    # `ẋ` has `x`'s shape at the activation scalar (§7.1), so the derivative of
    # a state field sits at the state field's own offset in the other buffer.
    DerivRead{P,typeof(s.i)}(offs[ci] + _leaf_offset(d.x, s.field), s.i)
end

function _resolve_selector(s::GetOutput, label::Symbol, b::Build, act::Activation,
                       ::Vector{Int}, viol::Vector{Diagnostic})
    ci = _read_component(s, label, b.flat, viol)
    ci === nothing && return nothing
    d = act.decls[ci]
    haskey(d.outs, s.name) ||
        (push!(viol, _declares(label, s, :output_port, d.outs)); return nothing)
    _check_index(s, label, d.outs[s.name], viol) || return nothing
    addr = act.layout.addr[(s.path, s.name)]
    CellRead{typeof(addr),typeof(s.i)}(addr, s.i)
end

function _resolve_selector(s::GetInput, label::Symbol, b::Build, act::Activation,
                       ::Vector{Int}, viol::Vector{Diagnostic})
    if !(s.face in b.flat.root_inputs)
        push!(viol, _rviol(label, s, :unknown_root_input; field = s.face,
                           candidates = b.flat.root_inputs))
        return nothing
    end
    addr = act.layout.addr[("", s.face)]
    CellRead{typeof(addr),Nothing}(addr, nothing)
end

function _resolve_selector(s::GetFace, label::Symbol, b::Build, act::Activation,
                       ::Vector{Int}, viol::Vector{Diagnostic})
    exported = Symbol[f for ((p, f), _) in b.flat.out_faces if isempty(p)]
    if !(s.name in exported)
        push!(viol, s.name in b.flat.root_inputs ?
                    _rviol(label, s, :root_input_not_face; field = s.name) :
                    _rviol(label, s, :unknown_output_face; field = s.name,
                           candidates = exported))
        return nothing
    end
    addr = act.layout.addr[("", s.name)]
    CellRead{typeof(addr),Nothing}(addr, nothing)
end
