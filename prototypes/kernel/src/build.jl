# The build pipeline (§9.2, §9.3, §5.5). Plain printable data in, a compiled
# executor out. Three jobs, in order, because each needs the previous one's
# answer:
#
#   1. probe stage 1 — needs no wiring at all (an `h_x` bundle carries no `u`),
#      which is what makes it the fixed point the rest of the build stands on;
#   2. derive the feedthrough graph and schedule stage 2 topologically, since a
#      stage-1 port carries no dependence and therefore breaks would-be loops;
#   3. probe stage 2 in that order (real upstream values available by
#      construction) and `f` last, against the now-complete table.

struct ComponentSpec
    path::Symbol
    comp::Any
    conns::Vector{Pair{Symbol,Tuple{Symbol,Symbol}}}   # face => (producer path, port)
end

ComponentSpec(path, comp) = ComponentSpec(path, comp, Pair{Symbol,Tuple{Symbol,Symbol}}[])

struct ModelSpec
    comps::Vector{ComponentSpec}
end

function index_of(spec::ModelSpec, path::Symbol)
    i = findfirst(c -> c.path === path, spec.comps)
    i === nothing && throw(BuildError("no component at path $path"))
    i
end

struct BuildError <: Exception
    msg::String
end
Base.showerror(io::IO, e::BuildError) = print(io, "BuildError: ", e.msg)

# --- 1. declarations at the activation scalar ---------------------------------

struct Decls
    x::NamedTuple          # initial state values, walked to the activation scalar
    ins::NamedTuple        # face => type, evaluated at the activation scalar
    outs::NamedTuple       # port => type, evaluated at the activation scalar
end

# The two registers split by D-166's criterion: `init_x` is by value at nominal
# `Float64` and its types are *walked*; `input_types`/`output_types` are
# functions of the activation scalar and are *evaluated*. There is no
# output-side leaf walk — the cell types at an activation are literally what the
# declaration returns at that `T`.
function declarations(cs::ComponentSpec, ::Type{T}) where {T}
    c = cs.comp
    _reject_one_arg(cs, output_types, "output_types")
    _reject_one_arg(cs, input_types, "input_types")
    Decls(retype_value(T, init_x(c)), input_types(c, T), output_types(c, T))
end

# The whole-signature forgotten-`T`. D-166 calls this class extinct because the
# tier is read off the declaration shape; on a single-tier prototype the `::Any`
# fallback would swallow it instead, and the component would declare nothing.
function _reject_one_arg(cs::ComponentSpec, fn, name)
    hasmethod(fn, Tuple{typeof(cs.comp)}) &&
        throw(BuildError("$(cs.path): `$name` is declared in the one-argument (discrete-tier) " *
                         "form; a continuous component declares " *
                         "`$name(::C, ::Type{T}) where {T <: Real}` (D-166/D-167)"))
end

# --- 2. probing ---------------------------------------------------------------

"""
Probe `h_x` for every component: no wiring is needed, so this runs first and
tells the rest of the build which ports are stage-1 — the ones that carry no
dependence on inputs and therefore break loops (§5.3).
"""
function probe_stage1(spec::ModelSpec, decls::Vector{Decls}, ::Type{T}) where {T}
    map(zip(spec.comps, decls)) do (cs, d)
        has_stage(h_x, cs.comp) || return NamedTuple()
        bn = bundle_names(h_x, cs.comp, ())
        y = h_x(cs.comp, _bundle_values(bn, d, NamedTuple(), NamedTuple(), T))
        y isa NamedTuple ||
            throw(BuildError("$(cs.path): h_x must return a NamedTuple of port values, got $(typeof(y))"))
        _check_ports(cs.path, "h_x", y, d.outs, T)
        _embed_ports(y, d.outs, T)
    end
end

function _bundle_values(bn, d::Decls, u, y_x, ::Type{T}; y = NamedTuple()) where {T}
    vals = map(bn) do n
        n === :x   ? d.x :
        n === :u   ? u :
        n === :y_x ? y_x :
        n === :y   ? y :
        n === :t   ? zero(T) : error("no probe source for $n")
    end
    NamedTuple{bn}(vals)
end

function _check_ports(path, stage, y::NamedTuple, outs::NamedTuple, ::Type{T}) where {T}
    for (name, v) in pairs(y)
        haskey(outs, name) ||
            throw(BuildError("$path: $stage returns `$name`, which `output_types` does not declare"))
        _accepts(outs[name], typeof(v), T) ||
            throw(BuildError("$path: $stage returns `$name`::$(typeof(v)), declared $(outs[name])" *
                             _pin_hint(outs[name], typeof(v), T)))
    end
end

# --- embed-accept (D-166) -----------------------------------------------------
# A declared-`T` leaf accepts exactly two arrivals: the activation scalar, and a
# `Float64` **embedded as a zero-partial** — which is what keeps the
# constant-branch idiom (`flow > 0 ? f(x) : 0.0`) legal as written at a `Dual`
# activation. Every other leaf — a deliberately pinned `Float64`, an `Int`, a
# `Bool` — is matched exactly, so an observed `Dual` at a pinned leaf is an
# error with a hint rather than a silent narrowing.

"""Is a value of type `V` a lawful arrival at a cell declared `P`, at activation `T`?"""
function _accepts(::Type{P}, ::Type{V}, ::Type{T}) where {P,V,T}
    P === V && return true
    if P <: Real
        return V <: Real && P === T && (V === T || V === Float64)
    elseif P <: StaticArray
        return V <: StaticArray && size(P) === size(V) && _accepts(eltype(P), eltype(V), T)
    else
        P.name === V.name || return false
        fieldcount(P) === fieldcount(V) || return false
        return all(_accepts(FP, FV, T) for (FP, FV) in zip(fieldtypes(P), fieldtypes(V)))
    end
end

# The one honest cause of a `Dual` at a leaf declared `Float64`, per D-166.
_pin_hint(::Type{P}, ::Type{V}, ::Type{T}) where {P,V,T} =
    T === Float64 || P === T ? "" :
    " — if this leaf participates in differentiation, declare it `T`"

"""
What the cell will hold: an accepted `Float64` arrival stored into the
activation buffer *is* a zero-partial, so the probe hands downstream the
embedded value rather than the literal the stage returned. Without this the
probe's product types diverge from the ones the runtime gather produces.
"""
_embed(::Type{P}, v, ::Type{T}) where {P,T} =
    typeof(v) === P ? v : reconstruct(P, T[T(l) for l in _leaf_values(v)], 0)

_embed_ports(y::NamedTuple, outs::NamedTuple, ::Type{T}) where {T} =
    NamedTuple{keys(y)}(map(n -> _embed(outs[n], y[n], T), keys(y)))

# --- 3. the feedthrough graph and the stage-2 schedule -------------------------

"""
Edges run producer → consumer for every consumed **stage-2** port; consuming a
stage-1 port adds no edge, which is the whole structural payoff of the split.
Returns a topological order over component indices, or reports the cycle.
"""
function schedule_stage2(spec::ModelSpec, stage1::Vector)
    n = length(spec.comps)
    deps = [Int[] for _ in 1:n]
    for (ci, cs) in enumerate(spec.comps)
        has_stage(h_xu, cs.comp) || continue
        for (_, (ppath, pport)) in cs.conns
            pi = index_of(spec, ppath)
            haskey(stage1[pi], pport) && continue        # stage-1 port: no dependence
            push!(deps[ci], pi)
        end
    end

    order = Int[]
    ready = [ci for ci in 1:n if isempty(deps[ci])]
    remaining = Set(1:n)
    while !isempty(ready)
        ci = popfirst!(ready)
        push!(order, ci)
        delete!(remaining, ci)
        for cj in collect(remaining)
            if ci in deps[cj]
                deps[cj] = filter(!=(ci), deps[cj])
                isempty(deps[cj]) && push!(ready, cj)
            end
        end
    end

    if !isempty(remaining)
        cycle = sort!(collect(remaining))
        names = join(("$(spec.comps[ci].path)" for ci in cycle), " → ")
        throw(BuildError("algebraic loop through stage-2 ports: $names — break it with a " *
                         "stage-1 (`h_x`) port, which carries no input dependence (§5.4/§5.5)"))
    end
    order
end

# --- 4. cell layout -----------------------------------------------------------
# Every declared port gets a cell; so does every unwired input face, which is a
# root slot — the one terminal with no producer, its initial value synthesized
# by `probe_value` (§9.3).

struct Layout
    addr::Dict{Tuple{Symbol,Symbol},Any}     # (path, port|face) => CellAddr
    slots::Vector{Tuple{Symbol,Symbol,Any}}  # root slots, with their probe values
    nleaf::Int
end

function cell_layout(spec::ModelSpec, decls::Vector{Decls}, ::Type{T}) where {T}
    addr = Dict{Tuple{Symbol,Symbol},Any}()
    slots = Tuple{Symbol,Symbol,Any}[]
    off = 0
    for (cs, d) in zip(spec.comps, decls)
        for (port, P) in pairs(d.outs)
            _check_cell_eltype(cs.path, "port", port, P, T)
            addr[(cs.path, port)] = CellAddr{P}(off)
            off += nleaves(P)
        end
    end
    for (cs, d) in zip(spec.comps, decls)
        wired = Set(first(p) for p in cs.conns)
        for (face, P) in pairs(d.ins)
            face in wired && continue
            _check_cell_eltype(cs.path, "root slot", face, P, T)
            addr[(cs.path, face)] = CellAddr{P}(off)
            off += nleaves(P)
            push!(slots, (cs.path, face, probe_value(P)))
        end
    end
    Layout(addr, slots, off)
end

# Every cell leaf must be the activation scalar, because this increment builds
# the per-eltype store (D-162) at exactly one eltype: a continuous, all-walking
# model needs no other. So a pinned `Float64` leaf (D-166's schema-visible
# freezing) under a non-nominal activation, and a discrete `Int`/`Bool` leaf,
# both have nowhere of their own to live. Rejecting is the point: `flatten!`
# would otherwise convert the value into the `T` buffer and the declared pin
# would be silently promoted to a zero-partial `Dual`. Lifted in increment 3,
# where the discrete tier's cells force the second store anyway.
function _check_cell_eltype(path, kind, name, ::Type{P}, ::Type{T}) where {P,T}
    for L in leaf_types(P)
        L === T && continue
        throw(BuildError("$path: $kind `$name` declares a $L leaf at activation $T — " *
                         "pinned and off-scalar leaves need their own store, which this " *
                         "increment does not build (D-166, §9.7)"))
    end
end

"""Address of the cell feeding `face` on `cs`: the wired producer's, or its own root slot."""
function input_addr(spec::ModelSpec, layout::Layout, cs::ComponentSpec, face::Symbol)
    for (f, (ppath, pport)) in cs.conns
        f === face || continue
        haskey(layout.addr, (ppath, pport)) ||
            throw(BuildError("$(cs.path).$face is wired to $ppath.$pport, which declares no such port"))
        return layout.addr[(ppath, pport)]
    end
    layout.addr[(cs.path, face)]
end

# --- 5. the whole build -------------------------------------------------------

function build(spec::ModelSpec, ::Type{T} = Float64; chunk_size::Int = 16) where {T}
    decls = [declarations(cs, T) for cs in spec.comps]
    stage1 = probe_stage1(spec, decls, T)
    order = schedule_stage2(spec, stage1)
    layout = cell_layout(spec, decls, T)

    store = CellStore(zeros(T, layout.nleaf))
    xbuf = T[]
    x_offs = Int[]
    for d in decls
        push!(x_offs, length(xbuf))
        append!(xbuf, T[T(l) for l in _leaf_values(d.x)])
    end
    ẋbuf = zeros(T, length(xbuf))
    clock = Clock(zero(T))

    # root slots hold their synthesized values until a writer replaces them
    for (path, face, v) in layout.slots
        scatter!(store, layout.addr[(path, face)], v)
    end

    addr_group(path, names) =
        NamedTuple{tuple(names...)}(tuple((layout.addr[(path, n)] for n in names)...))
    in_group(cs, d) =
        NamedTuple{tuple(keys(d.ins)...)}(
            tuple((input_addr(spec, layout, cs, face) for face in keys(d.ins))...))

    hx_entries, hxu_entries, rhs_entries = Any[], Any[], Any[]
    products = NamedTuple[NamedTuple() for _ in spec.comps]   # probe products, per component

    for (ci, cs) in enumerate(spec.comps)
        d, s1 = decls[ci], stage1[ci]
        products[ci] = s1
        has_stage(h_x, cs.comp) || continue
        bn = bundle_names(h_x, cs.comp, ())
        push!(hx_entries, StageEntry{typeof(h_x),typeof(cs.comp),typeof(d.x),bn,
                                     typeof(NamedTuple()),typeof(NamedTuple()),
                                     typeof(addr_group(cs.path, keys(s1))),typeof(clock)}(
            h_x, cs.comp, NamedTuple(), NamedTuple(), addr_group(cs.path, keys(s1)), x_offs[ci], clock))
    end

    for ci in order
        cs, d, s1 = spec.comps[ci], decls[ci], stage1[ci]
        has_stage(h_xu, cs.comp) || continue
        bn = bundle_names(h_xu, cs.comp, tuple(keys(s1)...))
        u = NamedTuple{tuple(keys(d.ins)...)}(tuple(
            (_probe_input(spec, layout, products, cs, face, d.ins[face], T) for face in keys(d.ins))...))
        y2 = h_xu(cs.comp, _bundle_values(bn, d, u, s1, T))
        y2 isa NamedTuple ||
            throw(BuildError("$(cs.path): h_xu must return a NamedTuple, got $(typeof(y2))"))
        _check_ports(cs.path, "h_xu", y2, d.outs, T)
        isempty(intersect(keys(s1), keys(y2))) ||
            throw(BuildError("$(cs.path): $(join(intersect(keys(s1), keys(y2)), ", ")) produced by two stages"))
        y2 = _embed_ports(y2, d.outs, T)
        products[ci] = merge(s1, y2)

        yx_g, out_g, in_g = addr_group(cs.path, keys(s1)), addr_group(cs.path, keys(y2)), in_group(cs, d)
        push!(hxu_entries, StageEntry{typeof(h_xu),typeof(cs.comp),typeof(d.x),bn,
                                      typeof(in_g),typeof(yx_g),typeof(out_g),typeof(clock)}(
            h_xu, cs.comp, in_g, yx_g, out_g, x_offs[ci], clock))
    end

    # Completeness of the declaration set (§8.2), for every component and not
    # only the stateful ones: an unproduced port would otherwise own a cell that
    # no stage ever writes, and read as a silent zero forever.
    for (ci, cs) in enumerate(spec.comps)
        missing_ports = setdiff(keys(decls[ci].outs), keys(products[ci]))
        isempty(missing_ports) ||
            throw(BuildError("$(cs.path): declared port(s) $(join(missing_ports, ", ")) produced by no stage"))
    end

    for (ci, cs) in enumerate(spec.comps)
        d = decls[ci]
        isempty(d.x) && continue
        has_stage(f, cs.comp) ||
            throw(BuildError("$(cs.path) declares `init_x` but defines no `f`"))
        bn = bundle_names(f, cs.comp, tuple(keys(stage1[ci])...))
        u = NamedTuple{tuple(keys(d.ins)...)}(tuple(
            (_probe_input(spec, layout, products, cs, face, d.ins[face], T) for face in keys(d.ins))...))
        ẋ = f(cs.comp, _bundle_values(bn, d, u, stage1[ci], T; y = products[ci]))
        _check_derivative(cs.path, ẋ, d.x)

        y_g, in_g = addr_group(cs.path, keys(d.outs)), in_group(cs, d)
        push!(rhs_entries, RHSEntry{typeof(cs.comp),typeof(d.x),bn,typeof(in_g),typeof(y_g),typeof(clock)}(
            cs.comp, in_g, y_g, x_offs[ci], clock))
    end

    bodies = (sweep_hx = chunked_body(hx_entries, store, xbuf, ẋbuf, clock; chunk_size),
              sweep_hxu = chunked_body(hxu_entries, store, xbuf, ẋbuf, clock; chunk_size),
              rhs = chunked_body(rhs_entries, store, xbuf, ẋbuf, clock; chunk_size),
              ticks = chunked_body(Any[], store, xbuf, ẋbuf, clock; chunk_size))

    Simulation(store, xbuf, ẋbuf, clock, bodies, layout, spec, T)
end

# A probed input value: the upstream product where wired, the synthesized slot
# value at a root. Stage-2 probing runs in topological order, so upstream
# products exist by construction.
#
# The entry is a *bound*, read permissively (D-167): a `T` entry is tolerant of
# both lawful arrivals, a pinned `Float64` entry demands a frozen one. The value
# passed on is the producer's, unembedded — the consumer gathers the producer's
# cell at runtime, so the cell's type is what its bundle carries.
function _probe_input(spec, layout, products, cs, face, P, ::Type{T}) where {T}
    for (f, (ppath, pport)) in cs.conns
        f === face || continue
        pi = index_of(spec, ppath)
        haskey(products[pi], pport) ||
            throw(BuildError("$(cs.path).$face reads $ppath.$pport, which no stage of $ppath produces"))
        v = products[pi][pport]
        _accepts(P, typeof(v), T) ||
            throw(BuildError("$(cs.path).$face declared $P, wired to $ppath.$pport::$(typeof(v))" *
                             _pin_hint(P, typeof(v), T)))
        return v
    end
    probe_value(P)
end

# §7.1: `Ẋ` has exactly `X`'s shape at the activation scalar. Checked
# structurally here so the runtime `flatten!` into the derivative block is safe.
function _check_derivative(path, ẋ, x::NamedTuple)
    ẋ isa NamedTuple ||
        throw(BuildError("$path: f must return a NamedTuple shaped like `init_x`, got $(typeof(ẋ))"))
    keys(ẋ) === keys(x) ||
        throw(BuildError("$path: f returns fields $(keys(ẋ)), state has $(keys(x)) — " *
                         "derivative completeness is structural (§7.1)"))
    for k in keys(x)
        nleaves(typeof(ẋ[k])) == nleaves(typeof(x[k])) ||
            throw(BuildError("$path: derivative field `$k` is $(typeof(ẋ[k])), state field is $(typeof(x[k]))"))
    end
end
