# The condition algebra and its resolution (§14.1–§14.3, §14.6): the inert lazy
# tree the fragment functions build, the collecting pass that validates it
# against a `Build`, and the plan §14.4's dynamic-walk register executes.
#
# Two properties carry the section. Composition performs no path arithmetic and
# no validation — `at` stores a prefix and never applies it — so a fragment
# function is stack-only construction and a misaddressed node fails at
# resolution, where the build is in hand. And resolution finishes before
# anything is written, which is what makes `init!` all-or-nothing (§14.6).

# --- the four node kinds (§14.2, §14.6) ---------------------------------------
# Every node is isbits except the prefix strings: a `Fragment`'s payloads are
# the author's own values, and the combinators hold their operands in tuples.

"""
Self-vocabulary payloads at one authoring level (§14.2): `x`, `s` and `m`
fields of the component addressed there, and `inputs` naming faces of *that
level's* contract. No paths — addressing children is exclusively `at`'s job.
"""
struct Fragment{X,S,M,L}
    x::X
    s::S
    m::M
    inputs::L
end

"""`at(prefix, node)`: scoping. Stores the prefix, never applies it (§14.2)."""
struct Scoped{N}
    prefix::String
    node::N
end

"""
`combine(nodes...)`: symmetric collection. Order is diagnostics only, and a
duplicate leaf is an error at resolution with both provenance chains (§14.2).
"""
struct Combined{T<:Tuple}
    nodes::T
end

"""
`override(base, patches...)`: ordered layering, the fourth node kind (§14.6).
The patch wins on a shared leaf and provenance keeps both sources; collisions
*within* one layer remain errors, and layering is variadic.
"""
struct Override{T<:Tuple}
    layers::T
end

const ConditionNode = Union{Fragment,Scoped,Combined,Override}

"""
    fragment(; x = (;), s = (;), m = (;), inputs = (;))

The leaf constructor of the condition tree (§14.2, Appendix B). Payloads are
NamedTuples in the authoring level's own vocabulary; a condition speaks state
(`x`, `s`), modes (`m`) and root inputs, never outputs and never workspace
(§14.1).
"""
function fragment(; x = (;), s = (;), m = (;), inputs = (;))
    for (name, p) in ((:x, x), (:s, s), (:m, m), (:inputs, inputs))
        p isa NamedTuple || throw(BuildError(ConditionNodeMisuse(
            observed = typeof(p), reason = :fragment_payload, payload = name)))
    end
    Fragment(x, s, m, inputs)
end

"""
    at(prefix, node)

Scope `node` under `prefix` (§14.2). Nothing is concatenated here: the prefix
is stored, and flattening at resolution is the one place path strings ever
join (§14.3).
"""
at(prefix::AbstractString, node::ConditionNode) = Scoped(String(prefix), node)
at(::AbstractString, other) = _node_misuse(other, ())

"""
    combine(nodes...)

Symmetric collection of sibling nodes (§14.2, D-204). Collision-intolerant by
design: two entries on one leaf are a resolution error naming both provenance
chains, and the layering spelling is `override`. It is not `Base.merge` and
does not extend it — `Base.merge` is last-writer-wins on NamedTuples, the
exact semantics D-065 rejects here.
"""
combine(nodes::ConditionNode...) = Combined(nodes)

# The blend of a node with a bare NamedTuple, both orders: an error method,
# raised at composition time — before any resolution pass or provenance chain
# exists, which is why it carries its own kind rather than a
# `ConditionResolution` sub-kind (§14.2).
combine(a::ConditionNode, b::NamedTuple) = _node_misuse(b, (nameof(typeof(a)),))
combine(a::NamedTuple, b::ConditionNode) = _node_misuse(a, (nameof(typeof(b)),))
combine(nodes...) = _misuse_in(nodes)

"""
    override(base, patches...)

Ordered layering (§14.6): a leaf present in several layers takes the last
layer's value, with provenance recording both sources. The use case is
"baseline plus tweaks" — the collision with `combine`'s duplicate-leaf error
*is* the intent — and trim commits `override(baseline, solution)`.

An input entry's leaf is the *root input* it resolves to, not the face it was
written with, so a full-coverage baseline authored at the root layers cleanly
under a patch a component's own `condition` method ships against its local
face.
"""
override(base::ConditionNode, patches::ConditionNode...) = Override((base, patches...))
override(layers...) = _misuse_in(layers)

_misuse_in(args::Tuple) = _node_misuse(args[findfirst(n -> !(n isa ConditionNode), args)],
                                       Tuple(nameof(typeof(n)) for n in args if n isa ConditionNode))

_node_misuse(v, in_hand) = throw(BuildError(ConditionNodeMisuse(
    observed = typeof(v), in_hand = Symbol[in_hand...])))

# --- flattening (§14.3) --------------------------------------------------------
# The only place path strings are ever concatenated: a recursion with a path
# accumulator, carrying the provenance chain beside it — the `at` prefixes and
# the payload position, which is what a collision diagnostic reports.
#
# An input entry's *leaf* is the root input its face resolves to, not the face
# it was written with, so the export chain is walked here, at the same moment
# the paths join, and the resolved face becomes the entry's key. Two spellings
# of one root input are then one leaf everywhere it matters: `override` layers
# them and `combine` collides on them, which is what §14.6's central use case —
# a root-level baseline under a component-local patch — requires.
#
# The recursion carries a third accumulator, the entry's **tree position**
# (§14.3): the `getfield`/`getindex` step tuple from the root node down to the
# authored value. The dynamic register ignores it — it bakes the value itself —
# and the specialized register lifts it to a `Getter{P}` lens, which is what
# lets one compiled plan be applied to every later tree of the same shape.

struct CEntry
    path::String
    store::Symbol        # :x | :s | :m | :input
    field::Symbol
    value::Any
    prov::String
    face::Union{Nothing,Symbol}   # input entries: the root input the chain lands on
    pos::Tuple                    # the tree position: the step tuple to this value
end

_key(e::CEntry) = e.face === nothing ? (e.path, e.store, e.field) : ("", :input, e.face)
_step(prov::String, s::String) = isempty(prov) ? s : prov * " → " * s

function _flat(n::Fragment, path::String, prov::String, pos::Tuple,
               flat::Flat, viol::Vector{Diagnostic})
    out = CEntry[]
    for (store, name, payload) in ((:x, :x, n.x), (:s, :s, n.s),
                                   (:m, :m, n.m), (:input, :inputs, n.inputs))
        for (field, v) in pairs(payload)
            e = CEntry(path, store, field, v, _step(prov, "fragment($name).$field"),
                       nothing, (pos..., name, field))
            store === :input &&
                (e = CEntry(e.path, e.store, e.field, e.value, e.prov,
                            _root_input(flat, e, viol), e.pos))
            push!(out, e)
        end
    end
    out
end

_flat(n::Scoped, path::String, prov::String, pos::Tuple,
      flat::Flat, viol::Vector{Diagnostic}) =
    _flat(n.node, _join(path, n.prefix), _step(prov, "at(\"$(n.prefix)\")"),
          (pos..., :node), flat, viol)

_flat(n::Combined, path::String, prov::String, pos::Tuple,
      flat::Flat, viol::Vector{Diagnostic}) =
    reduce(vcat, (_flat(k, path, _step(prov, "combine[$i]"), (pos..., :nodes, i), flat, viol)
                  for (i, k) in enumerate(n.nodes)); init = CEntry[])

# Layering (§14.6): each layer is flattened and checked on its own — a
# within-layer collision is still an error — and then folded onto the
# accumulator, the patch replacing the leaf it overrode and inheriting its
# provenance beside its own.
function _flat(n::Override, path::String, prov::String, pos::Tuple,
               flat::Flat, viol::Vector{Diagnostic})
    acc = CEntry[]
    for (i, layer) in enumerate(n.layers)
        label = i == 1 ? "override[base]" : "override[patch $(i - 1)]"
        es = _flat(layer, path, _step(prov, label), (pos..., :layers, i), flat, viol)
        _check_duplicates!(es, viol)
        for e in es
            j = findfirst(a -> _key(a) == _key(e), acc)
            j === nothing ? push!(acc, e) :
                (acc[j] = CEntry(e.path, e.store, e.field, e.value,
                                 "$(e.prov) (overrode $(acc[j].prov))", e.face, e.pos))
        end
    end
    acc
end

# `combine`'s one collision rule, collected rather than thrown (§13.1): both
# provenance chains and the directive naming the layering combinator.
function _check_duplicates!(es::Vector{CEntry}, viol::Vector{Diagnostic})
    seen = Dict{Tuple{String,Symbol,Symbol},CEntry}()
    for e in es
        k = _key(e)
        if haskey(seen, k)
            push!(viol, DuplicateConditionLeaf(path = e.path, store = e.store,
                                               field = e.field, face = e.face,
                                               provenance = [seen[k].prov, e.prov]))
        else
            seen[k] = e
        end
    end
    nothing
end

# --- the compiled plan (§14.3) --------------------------------------------------

"""
What a resolved condition compiles to (§14.3): per continuous-state leaf its
`xbuf` offset, per discrete or mode store the whole value the write installs —
`merge(defaults, overlay)`, the genuine last-wins NamedTuple merge baked here
(§14.1's fork, not the condition combinator) — and per root input its compiled
cell address. Every value has been through its leaf's converter already, so
application is a walk with no decisions left in it.

`faces` is the plan's root-input coverage, in entry order: the resolution-time
operand of §14.6's totality check.

`T` is the activation the plan was resolved at: the offsets, addresses and
converters baked here are that activation's, so `apply!` pairs plan and
executor by dispatch and a mismatch is the internal-invariant refusal build.jl
raises (§14.4).
"""
struct ConditionPlan{T}
    xs::Vector{Tuple{Int,Any}}             # (xbuf offset, value)
    stores::Vector{Tuple{Symbol,Int,Any}}  # (:s | :m, component index, whole store value)
    inputs::Vector{Tuple{Symbol,Any,Any}}  # (root face, cell address, value)
    faces::Vector{Symbol}
end

# --- resolution (§14.3) ---------------------------------------------------------

"""
    resolve_condition(node, b::Build, T = Float64) → ConditionPlan

Flatten the condition tree, validate every entry against `b` in §13.1's
collecting register — full list, violations collected, one `BuildError` — and
compile what survives to a plan.

The checks are §14.3's: the path resolves to a component, the field is
declared in that component's `init_x`/`init_s`/`init_m`, the value converts to
the declared leaf type, an input face resolves through the export chain to a
root input, and no leaf is written twice — an input entry's leaf being the
*root input* it resolves to, so two spellings of one root input are one leaf.
Schema is the
authority on *may you write this, at what type*; the activation's layout
supplies the destination.
"""
function resolve_condition(node::ConditionNode, b::Build, ::Type{T} = Float64) where {T}
    resolved, viol, act = _resolve_entries(node, b, T)
    _report_violations(viol)

    xs = Tuple{Int,Any}[]
    inputs = Tuple{Symbol,Any,Any}[]
    faces = Symbol[]
    overlays = Dict{Tuple{Symbol,Int},Vector{Pair{Symbol,Any}}}()
    for r in resolved
        e = r.e
        if e.store === :input
            push!(inputs, (e.face, r.dest, r.v))
            push!(faces, e.face)
        elseif e.store === :x
            push!(xs, (r.dest, r.v))
        else
            push!(get!(() -> Pair{Symbol,Any}[], overlays, (e.store, r.dest)), e.field => r.v)
        end
    end

    # Overlay partiality baked now (§14.3): the store's write is one whole
    # value, `merge(defaults, overlay)`, so application decides nothing.
    stores = Tuple{Symbol,Int,Any}[]
    for (store, ci, defaults) in _store_bases(b, act)
        ov = get(overlays, (store, ci), nothing)
        ov === nothing && continue
        push!(stores, (store, ci, convert(typeof(defaults), merge(defaults, (; ov...)))))
    end
    ConditionPlan{T}(xs, stores, inputs, faces)
end

# --- the collecting pass, shared by both registers (§14.3, §13.1) ---------------

"""
One entry that survived §14.3's checks, beside everything either register needs
to bake from it. The two registers differ in *what* they bake — the dynamic one
takes `v`, the specialized one lifts `e.pos` to a lens and keeps `L` as the
converter — and in nothing else, which is why the checks have one implementation.
"""
struct Resolved
    e::CEntry
    dest::Any   # :x → the `xbuf` offset; :s, :m → the component index; :input → the cell address
    L::Any      # the destination leaf type at this activation — §14.3's converter
    v::Any      # the authored value, through that converter
end

# §14.3's list, run once: the path resolves to a component, the field is
# declared in that component's `init_x`/`init_s`/`init_m`, the value converts to
# the declared leaf type, an input face resolves through the export chain to a
# root input, and no leaf is written twice. Violations are collected (§13.1) and
# handed back with the survivors, so a service that owns its own setup
# diagnostic can fold the same list into its kind.
function _resolve_entries(node::ConditionNode, b::Build, ::Type{T}) where {T}
    flat, tiers = b.flat, b.tiers
    act = activation(b, T)
    decls, layout = act.decls, act.layout
    viol = Diagnostic[]
    entries = _flat(node, "", "", (), flat, viol)
    _check_duplicates!(entries, viol)

    x_offs = _x_offsets(decls, tiers)
    out = Resolved[]
    for e in entries
        if e.store === :input
            e.face === nothing && continue     # the chain was reported at flattening
            addr = layout.addr[("", e.face)]
            P = _port_type(addr)
            (ok, v) = _convert(P, e.value)
            ok || (push!(viol, _unconvertible(e, e.value, P, T)); continue)
            push!(out, Resolved(e, addr, P, v))
            continue
        end
        ci = _component(flat, e, viol)
        ci === nothing && continue
        c, tier, d = flat.comps[ci], tiers[ci], decls[ci]
        declared = e.store === :x ? d.x : e.store === :s ? d.s : init_m(c)
        if isempty(declared)
            push!(viol, _no_store(e, tier))
            continue
        end
        haskey(declared, e.field) ||
            (push!(viol, _undeclared(e, c, tier, declared, T)); continue)
        L = typeof(declared[e.field])
        (ok, v) = _convert(L, e.value)
        ok || (push!(viol, _unconvertible(e, e.value, L, T)); continue)
        push!(out, Resolved(e, e.store === :x ? x_offs[ci] + _leaf_offset(d.x, e.field) : ci,
                            L, v))
    end
    (out, viol, act)
end

# The merge bases, in one order both registers walk: per component, the discrete
# store's declared defaults and then the mode store's (§14.3's fork).
_store_bases(b::Build, act::Activation) =
    [(store, ci, store === :s ? act.decls[ci].s : init_m(b.flat.comps[ci]))
     for ci in eachindex(b.flat.comps) for store in (:s, :m)]

# Anything that is not a node reaching a service entry point is the §14.2
# misuse, not a `MethodError`: the directive is the same one `combine` prints,
# and a bare NamedTuple handed to `init!` is exactly the slip it addresses.
resolve_condition(other, ::Build, ::Type = Float64) = _node_misuse(other, ())

# The `x` destinations, per component then per field: the flat buffer's layout
# is the declaration walk, so the plan bakes the offsets the state accessor
# recomputes (§7.1).
function _x_offsets(decls::Vector{Decls}, tiers::Vector{Tier})
    offs, n = Int[], 0
    for (d, t) in zip(decls, tiers)
        push!(offs, n)
        t === CONTINUOUS && (n += nleaves(typeof(d.x)))
    end
    offs
end

function _leaf_offset(x::NamedTuple, field::Symbol)
    off = 0
    for (k, v) in pairs(x)
        k === field && return off
        off += nleaves(typeof(v))
    end
    off
end

_convert(::Type{P}, v) where {P} =
    try
        (true, convert(P, v))
    catch
        (false, nothing)
    end

# One `ConditionResolution` off an entry: the leaf coordinates and the
# provenance chain are the entry's own, and each arm adds what it observed.
_cviol(e::CEntry, reason::Symbol; kw...) =
    ConditionResolution(; path = e.path, store = e.store, field = e.field, face = e.face,
                        reason = reason, provenance = e.prov, kw...)

# The component a non-input entry addresses. Assemblies are virtual for
# execution (§10.5) and own no state, so an `at` prefix stopping at one has
# nothing to write — and saying so beats "no such path".
function _component(flat::Flat, e::CEntry, viol::Vector{Diagnostic})
    i = findfirst(==(e.path), flat.paths)
    i === nothing || return i
    push!(viol, _cviol(e, any(startswith(p, e.path * "/") for p in flat.paths) ?
                          :assembly_path : :unknown_path))
    nothing
end

# An `inputs` payload names a face of the authoring level's contract, whatever
# that level's class; resolution walks the export chain to the root input and
# errors if the face never surfaces (§14.2). The chain is the `Build`'s own
# input-side face graph, total under one-level routing (§9.2, D-207): a face's
# producer is either a root input or an internal port, and a component-fed face
# reaches none — writing it would be meaningless, because the first sweep
# overwrites it.
function _root_input(flat::Flat, e::CEntry, viol::Vector{Diagnostic})
    if isempty(e.path)
        e.field in flat.root_inputs && return e.field
        push!(viol, _cviol(e, :unexported_face; candidates = flat.root_inputs))
        return nothing
    end
    k = findfirst(p -> first(p) === (e.path, e.field), flat.in_faces)
    if k === nothing
        here = [f for ((p, f), _) in flat.in_faces if p == e.path]
        push!(viol, isempty(here) && !_addresses_level(flat, e.path) ?
                    _cviol(e, :unknown_path) :
                    _cviol(e, :no_input_face; candidates = here))
        return nothing
    end
    (path, port) = last(flat.in_faces[k])
    isempty(path) && return port
    push!(viol, _cviol(e, :internally_wired; producer = (path, port)))
    nothing
end

# Does the path name a level of this build at all — a component, or an assembly
# some component sits under? Assemblies leave no row of their own in the flat
# list, so an empty face list alone does not tell a bare typo from an assembly
# that declares no input face.
_addresses_level(flat::Flat, path::String) =
    any(p == path || startswith(p, path * "/") for p in flat.paths)

# The store a condition names has to exist on the component at all: `x` is the
# continuous tier's state and `s` the discrete one's, disjoint by construction
# (D-195), and `m` is continuous-only (§3.2).
_no_store(e::CEntry, tier::Tier) =
    _cviol(e, :no_store; tier = Symbol(tier_word(tier)))

# An undeclared field, discriminated against the component's other name
# families: a condition specifies state, modes and root inputs — never outputs,
# which are derived data, and never workspace (§14.1).
function _undeclared(e::CEntry, c, tier::Tier, declared::NamedTuple, ::Type{T}) where {T}
    role = haskey(declared_at(output_types, c, tier), e.field) ? :output_port :
           haskey(declared_at(input_types, c, tier), e.field) ? :input_face :
           (_declares_workspace(c, tier) &&
            haskey(_declared_workspace(c, tier, T), e.field)) ? :workspace : nothing
    _cviol(e, :undeclared_field; candidates = collect(keys(declared)), role = role)
end

_declared_workspace(c, tier::Tier, ::Type{T}) where {T} =
    tier === CONTINUOUS ? workspace(c, T) : workspace(c)

# The one refusal §14.3's converter table cannot bake around. Its second clause
# is the non-nominal case: at a seeded activation the leaves a decision descends
# into are the ones the activation retyped, and a *frozen* discrete `s` (§9.4,
# D-166) or a leaf pinned `Float64` by its own declaration is not one of them —
# so the value cannot be carried and there is nowhere to put its partials.
function _unconvertible(e::CEntry, v, ::Type{P}, ::Type{T}) where {P,T}
    _cviol(e, :unconvertible; declared = P, observed = typeof(v), value = v,
           activation = _seeded_into_pinned(typeof(v), P, T) ? T : nothing)
end

_seeded_into_pinned(::Type{V}, ::Type{P}, ::Type{T}) where {V,P,T} =
    T !== Float64 && T in leaf_types(V) && !(T in leaf_types(P))

# §13.1's collecting register: the full list, every violation, one throw.
function _report_violations(viol::Vector{Diagnostic})
    isempty(viol) && return nothing
    throw(BuildError(viol))
end

# --- root-input totality (§14.6) ------------------------------------------------

"""
§14.6's precondition of starting, checked by the service and *not* a property
of conditions, which are legitimately partial: an application establishing a
complete world over virgin stores covers every root input. Coverage is a
plan-level fact — both operands are resolution-time data — so the check is one
comparison and runs before any evaluation, not merely before any write. A
shortfall is one collected, declaration-ordered `UninitializedInputs` naming
every uncovered face (D-068).

The services path contains no call to `probe_value`: a root input gets a
condition value or the application errors, and there is no third branch. A
fabricated zero is a fine probe input and a terrible flight condition.
"""
function assert_total(plan::ConditionPlan, flat::Flat, op::Symbol)
    covered = Set(plan.faces)
    uncovered = [f for f in flat.root_inputs if !(f in covered)]
    isempty(uncovered) && return nothing
    throw(BuildError(UninitializedInputs(op = op, faces = uncovered)))
end

# --- the dynamic-walk application register (§14.4) ------------------------------

"""
    apply!(ex::Executor, plan)
    apply!(sim, plan)

§14.4's dynamic walk: execute the validated entry list by runtime dispatch per
write — microseconds total, allocation permitted, the stopped-sim path never
having been under §7.5's zero-alloc regime — with no per-shape codegen, so
fifty structurally different scripted conditions cost fifty walks rather than
fifty compiles. Which register a service uses is internal, never user-facing:
the specialized `apply!` below is the other register over the same checks, for
the services that hold one shape fixed and vary its values.
"""
function apply!(ex::Executor{T}, plan::ConditionPlan{T}) where {T}
    for (off, v) in plan.xs
        flatten!(ex.xbuf, off, v)
    end
    for (store, ci, v) in plan.stores
        (store === :s ? ex.sstores : ex.mstores)[ci][] = v
    end
    for (_, addr, v) in plan.inputs
        scatter!(ex.store, addr, v)
    end
    nothing
end

apply!(::Executor{S}, ::ConditionPlan{T}) where {S,T} =
    _activation_mismatch("plan", T, S)

apply!(sim::Simulation, plan::ConditionPlan) = apply!(sim.exec, plan)

# --- the specialized application register (§14.3, §14.4, D-066) -----------------
# The other register over the same checks. The dynamic walk bakes *values*, so a
# plan is good for one tree; this one bakes *lenses*, so a plan compiled from a
# tree's shape applies to every later tree of that shape. That is the trade
# §14.4 states: ~10–50 ms of codegen once per shape, against a per-iteration
# write with no strings, no dispatch and no allocation — the shape being fixed
# and the values varying is exactly what an iterating service does.
#
# All string work, addressing and validation are functions of the shape, so all
# of it happens here, once. What survives into `apply!` is a tuple of baked
# writes and a tuple of prefix compares.

"""
The lens (§14.3, glossary): a condition entry's tree position lifted to a type
parameter, callable on any tree of the shape it was compiled from. Navigation is
generated from `P`, so the access is a chain of static `getfield`/`getindex`
steps the compiler folds into offsets — the authored value reached with no
search and no runtime fact consulted.
"""
struct Getter{P} end

@generated function (::Getter{P})(tree) where {P}
    ex = :tree
    for step in P
        ex = step isa Symbol ? :(getfield($ex, $(QuoteNode(step)))) : :(getindex($ex, $step))
    end
    quote
        $(Expr(:meta, :inline))
        $ex
    end
end

"""
One leaf's compiled read half: the lens that finds the authored value in the
tree, and `L`, the destination leaf type at this activation — which *is*
§14.3's converter, selected once at resolution and never consulted again.

Both of §14.3's cases are this one call. A leaf already at the activation's
scalar takes the type's own methods, partials flowing through untouched; a
plain `Float64` leaf against a seeded activation takes the zero-partial
embedding `convert` already performs, which is semantically exact for a value
held at the operating point and in no other case.
"""
struct Authored{P,L} end

@inline (::Authored{P,L})(tree) where {P,L} = convert(L, Getter{P}()(tree))

# The three destinations, one write each: the flat state buffer at a baked
# offset, a component's own store as one whole value, and a root input's cell.

struct XWrite{A}
    authored::A
    off::Int
end

struct InputWrite{A,D}
    authored::A
    addr::D
end

# The `s`/`m` fork (§14.3): the store's write is one whole value,
# `merge(defaults, overlay)`, with the base baked and the overlay's fields read
# through their lenses. `K` is `:s` or `:m` and `S` the store's value type, both
# in the type — `S` because the stores are held in a `Vector{Any}` and the
# assertion has to go on the *reference*: asserting the dereferenced value
# instead leaves the `[]` a dynamic call, which boxes.
struct StoreWrite{K,S,F,A<:Tuple}
    ci::Int
    defaults::S
    authored::A
end

StoreWrite{K,S,F}(ci::Int, defaults::S, authored::A) where {K,S,F,A<:Tuple} =
    StoreWrite{K,S,F,A}(ci, defaults, authored)

@inline _write!(w::XWrite, ex::Executor, tree) =
    flatten!(ex.xbuf, w.off, w.authored(tree))

@inline _write!(w::InputWrite, ex::Executor, tree) =
    scatter!(ex.store, w.addr, w.authored(tree))

@inline function _write!(w::StoreWrite{K,S,F}, ex::Executor, tree) where {K,S,F}
    ov = NamedTuple{F}(map(a -> a(tree), w.authored))
    ((K === :s ? ex.sstores : ex.mstores)[w.ci]::Base.RefValue{S})[] =
        convert(S, merge(w.defaults, ov))
    nothing
end

# One `Scoped` node's prefix: the tree type carries the nesting, every field
# name and every leaf type, but a prefix is a runtime `String` field, so the
# plan records the one it was compiled from and `apply!` closes the remainder
# with a `===` compare. On `String` that compare is *content* equality rather
# than pointer identity — `egal` is specialized for it — so what the sweep tests
# is that the prefix still spells the same path, whether it was written as a
# literal or computed afresh per iteration (`at("gear/$i", …)`); the cost is a
# length check and a short `memcmp`, and in the all-literal case the compiler
# folds it away.
struct Prefix{P}
    expected::String
end

"""
What a condition *shape* compiles to (§14.3, §14.4): the tree type it was
compiled from in the plan's own type, and its writes as tuples, so `apply!`
unrolls into the same machine operations an in-place write would be.

Application is `apply!(ex, plan, tree)` for any tree of that shape. Nothing is
decided there: the destinations, the converters and the merge bases were all
settled at compile time, and what is left is the fold-away shape check plus the
writes.

`T`, the leading parameter, is the activation the plan was compiled at — the
same identity `Reader` carries, for the same reason: the lenses' converters and
the writes' destinations are one activation's, and pairing them with another's
executor is refused as an internal-invariant violation (build.jl, §14.4) rather
than a wrong slot.
"""
struct SpecializedPlan{T,NT,XS<:Tuple,ST<:Tuple,IN<:Tuple,PF<:Tuple}
    xs::XS
    stores::ST
    inputs::IN
    prefixes::PF
end

SpecializedPlan{T,NT}(xs::XS, stores::ST, inputs::IN, prefixes::PF) where {T,NT,XS,ST,IN,PF} =
    SpecializedPlan{T,NT,XS,ST,IN,PF}(xs, stores, inputs, prefixes)

"""
    compile_plan(node, b::Build, T = Float64) → SpecializedPlan

Compile a condition tree's **shape** into the specialized register's plan,
running exactly the checks `resolve_condition` runs — one implementation, in
§13.1's collecting register, the two registers differing only in what they bake
(§14.3).

The tree handed here is a genuine tree, and its values matter to exactly one
check: convertibility, which is where a decision variable authored into a
frozen or pinned leaf is refused. Everything else the compile reads is shape —
paths, field names, tree positions, the destination leaf types at this
activation.
"""
function compile_plan(node::ConditionNode, b::Build, ::Type{T} = Float64) where {T}
    resolved, viol, act = _resolve_entries(node, b, T)
    _report_violations(viol)

    xs, inputs = Any[], Any[]
    overlays = Dict{Tuple{Symbol,Int},Vector{Resolved}}()
    for r in resolved
        e = r.e
        if e.store === :input
            push!(inputs, InputWrite(Authored{e.pos,r.L}(), r.dest))
        elseif e.store === :x
            push!(xs, XWrite(Authored{e.pos,r.L}(), r.dest))
        else
            push!(get!(() -> Resolved[], overlays, (e.store, r.dest)), r)
        end
    end

    stores = Any[]
    for (store, ci, defaults) in _store_bases(b, act)
        ov = get(overlays, (store, ci), nothing)
        ov === nothing && continue
        push!(stores, StoreWrite{store,typeof(defaults),Tuple(r.e.field for r in ov)}(
            ci, defaults, Tuple(Authored{r.e.pos,r.L}() for r in ov)))
    end

    SpecializedPlan{T,typeof(node)}(Tuple(xs), Tuple(stores), Tuple(inputs),
                                    Tuple(Prefix{p}(v) for (p, v) in _scoped_prefixes(node)))
end

compile_plan(other, ::Build, ::Type = Float64) = _node_misuse(other, ())

# Every `Scoped` node's prefix field, positioned, in tree order. A separate walk
# from `_flat`'s: a prefix belongs to a *node*, not to a leaf, and a scope
# wrapping no payload at all still has a prefix that can drift.
function _scoped_prefixes(node::ConditionNode)
    out = Tuple{Tuple,String}[]
    _scoped!(node, (), out)
    out
end

_scoped!(::Fragment, ::Tuple, ::Vector) = nothing
_scoped!(n::Scoped, pos::Tuple, out::Vector) =
    (push!(out, ((pos..., :prefix), n.prefix)); _scoped!(n.node, (pos..., :node), out))
_scoped!(n::Combined, pos::Tuple, out::Vector) =
    for (i, k) in enumerate(n.nodes)
        _scoped!(k, (pos..., :nodes, i), out)
    end
_scoped!(n::Override, pos::Tuple, out::Vector) =
    for (i, k) in enumerate(n.layers)
        _scoped!(k, (pos..., :layers, i), out)
    end

"""
    apply!(ex::Executor{T}, plan::SpecializedPlan{T,NT}, tree::NT)

§14.4's specialized walk: write every leaf of `tree` through its baked lens and
converter — `x` leaves flattened at their offsets, each `s` and `m` store as one
whole `merge(defaults, overlay)` value, root inputs scattered — with no
allocation, no string work and no dispatch left in the path.

**The shape check is folded.** The tree type is proven by dispatch: `NT` carries
the full nesting, every field name and every leaf type, so a tree of another
shape does not match this method and reaches the fallback below. The `===` sweep
over the `Scoped` prefixes closes the remainder, those being runtime fields the
type cannot carry. This is §9.5's mechanism transferred: on conformant code the
compiler decides it and deletes it, and drift is a structured error rather than
silent corruption.

The sweep runs before any write, so a refused application leaves the executor
exactly as it found it.
"""
function apply!(ex::Executor{T}, plan::SpecializedPlan{T,NT}, tree::NT) where {T,NT}
    _sweep_prefixes(plan.prefixes, tree)
    _writes!(plan.xs, ex, tree)
    _writes!(plan.stores, ex, tree)
    _writes!(plan.inputs, ex, tree)
    nothing
end

apply!(::Executor{T}, ::SpecializedPlan{T,NT}, tree) where {T,NT} =
    _shape_drift(NT, typeof(tree))

apply!(::Executor{S}, ::SpecializedPlan{T}, tree) where {S,T} =
    _activation_mismatch("plan", T, S)

@inline _writes!(::Tuple{}, ::Executor, tree) = nothing
@inline function _writes!(ws::Tuple, ex::Executor, tree)
    _write!(first(ws), ex, tree)
    _writes!(Base.tail(ws), ex, tree)
end

@inline _sweep_prefixes(::Tuple{}, tree) = nothing
@inline function _sweep_prefixes(ps::Tuple, tree)
    _compare(first(ps), tree)
    _sweep_prefixes(Base.tail(ps), tree)
end

@inline function _compare(p::Prefix{P}, tree) where {P}
    observed = Getter{P}()(tree)
    observed === p.expected || _prefix_drift(P, p.expected, observed)
    nothing
end

@noinline _shape_drift(::Type{NT}, ::Type{O}) where {NT,O} = throw(BuildError(
    ConditionShapeDrift(reason = :tree_type, compiled = NT, observed = O)))

@noinline _prefix_drift(P::Tuple, expected::String, observed::String) = throw(BuildError(
    ConditionShapeDrift(reason = :prefix, compiled = expected, observed = observed,
                        position = P)))

# --- capture: the gather twin of the application register (§14.1, §14.10) -------

"""
    capture(sim) → (condition, t)

Read the current committed stores **and root inputs** back as a condition value
(§14.1, glossary): every component's `x` field by field, every discrete `s`,
every mode store and every root input, as one `combine` of per-component
fragments plus the root-input fragment. The pair is capture → tweak → apply:
`init!(sim, c; t₀ = t)` re-establishes exactly the world that was read, and a
`trim!` resumes from it as `trim!(sim, problem; baseline = c, t₀ = t)`
(§14.8).

Root-input coverage is what makes the captured condition **total**, hence
re-applicable under §14.6 — a capture is by construction the one condition that
never needs a baseline under it. The condition is time-free: `t` rides beside
it, because time is not a store of any component (§14.5), and it is passed back
as the `t₀` argument.

Legal in `initialized` and `stopped`, the states whose stores are committed and
boundary-consistent. A `built` simulation has no such stores yet — boundary
zero has never run — and `running` and `errored` are the ordinary service
refusals (§14, §11.3, §13.6); all four are one `ServiceLifecycle`.

Re-applying reproduces the captured world bit for bit, with one caveat that is
boundary zero's rather than capture's: the re-application runs the sequence
(§14.5), so `project` and any guard already holding in the captured state fire
again there — which is exactly what makes a warm restart a *fresh run from
these values* rather than a resumption.
"""
function capture(sim::Simulation{T}) where {T}
    lc = lifecycle(sim)
    lc in (:initialized, :stopped) || throw(BuildError(ServiceLifecycle(
        op = :capture, status = lc, legal = [:initialized, :stopped])))
    ex, flat, tiers = sim.exec, sim.build.flat, sim.build.tiers
    act = activation(sim.build, T)
    offs = _x_offsets(act.decls, tiers)
    nodes = ConditionNode[]
    for ci in eachindex(flat.comps)
        d = act.decls[ci]
        payload = NamedTuple()
        tiers[ci] === CONTINUOUS && !isempty(d.x) &&
            (payload = merge(payload, (x = _capture_x(d.x, ex.xbuf, offs[ci]),)))
        ex.sstores[ci] === nothing || (payload = merge(payload, (s = ex.sstores[ci][],)))
        ex.mstores[ci] === nothing || (payload = merge(payload, (m = ex.mstores[ci][],)))
        isempty(payload) && continue
        push!(nodes, at(flat.paths[ci], fragment(; payload...)))
    end
    isempty(flat.root_inputs) || push!(nodes, fragment(inputs =
        NamedTuple{Tuple(flat.root_inputs)}(Tuple(gather(ex.store, act.layout.addr[("", f)])
                                                  for f in flat.root_inputs))))
    (combine(nodes...), ex.clock.t)
end

# One component's `x` payload, field by field: the flat buffer is read through
# the same declaration walk the plan writes through, so what comes back is
# exactly what goes in (§7.1).
function _capture_x(x::NamedTuple, buf::Vector, base::Int)
    off, vals = base, Any[]
    for v in values(x)
        push!(vals, reconstruct(typeof(v), buf, off))
        off += nleaves(typeof(v))
    end
    NamedTuple{keys(x)}(Tuple(vals))
end
