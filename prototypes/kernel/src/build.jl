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
#
# Its input is the flattened tree (`src/assembly.jl`): primitives by absolute
# path, one resolved producer per input. Assemblies are virtual from here on.

# --- 1. declarations at the activation scalar ---------------------------------

struct Decls
    x::NamedTuple          # initial state values, walked to the activation scalar
    ins::NamedTuple        # face => type, evaluated at the activation scalar
    outs::NamedTuple       # port => type, evaluated at the activation scalar
end

# The two registers split by D-166's criterion: `init_x` is by value and, on the
# continuous tier, its types are *walked*; `input_types`/`output_types` are
# functions of the activation scalar there and are *evaluated*. There is no
# output-side leaf walk — the cell types at an activation are literally what the
# declaration returns at that `T`. On the discrete tier the plain forms declare
# the pinned world, so nothing walks and nothing is evaluated at `T`.
function declarations(c, t::Tier, ::Type{T}) where {T}
    t === CONTINUOUS ?
        Decls(retype_value(T, init_x(c)), input_types(c, T), output_types(c, T)) :
        Decls(init_x(c), declared_at(input_types, c, t), declared_at(output_types, c, t))
end

# --- tier classification (§8.2) -----------------------------------------------
# Tier is read off the declaration shape, never announced. For a **stateful**
# leaf the update law carries it — `f` continuous, `g` discrete, the stage names
# being shared (D-173). A **stateless** leaf has no update law, so its contract
# arities decide, `output_types` being mandatory and therefore the decider: the
# two-argument forms declare cells at the activation scalar, the plain forms the
# pinned discrete world (D-166/D-167). Every other tier-implying declaration
# must then agree, and disagreement names the offending one.
#
# The classifier sees primitives only: a component that declares nothing at all
# has no *class* to read, which §8.5 settles before this runs.

"""The tier the primitive at `path` announces, or a `BuildError` naming what disagrees."""
function classify_tier(path::String, c)
    votes = Tuple{Symbol,Tier}[]
    has_stage(f, c) && push!(votes, (:f, CONTINUOUS))
    has_stage(g, c) && push!(votes, (:g, DISCRETE))
    !isempty(init_m(c)) && push!(votes, (:init_m, CONTINUOUS))
    for (name, fn) in ((:output_types, output_types), (:input_types, input_types),
                       (:workspace, workspace))
        _declares(fn, c, Type{Float64}) && push!(votes, (name, CONTINUOUS))
        _declares(fn, c) && push!(votes, (name, DISCRETE))
    end

    # The decider, by §8.2's two cases.
    if !isempty(init_x(c))
        i = findfirst(v -> first(v) === :f || first(v) === :g, votes)
        i === nothing &&
            throw(BuildError("`$path` declares `init_x` but defines neither `f` nor `g` — " *
                             "a store needs its update (§8.2)"))
    else
        i = findfirst(v -> first(v) === :output_types, votes)
        i === nothing &&
            throw(BuildError("`$path` declares no `output_types` and owns no state — there is " *
                             "nothing for the tier to be read off (§8.2)"))
    end

    t = last(votes[i])
    for (name, vt) in votes
        vt === t ||
            throw(BuildError("`$path`: `$name` is declared in the $(tier_word(vt))-tier " *
                             "form, but this component's other declarations announce the " *
                             "$(tier_word(t)) tier (§8.2)"))
    end
    t
end

# --- 2. probing ---------------------------------------------------------------

"""
Probe `h_x` for every component: no wiring is needed, so this runs first and
tells the rest of the build which ports are stage-1 — the ones that carry no
dependence on inputs and therefore break loops (§5.3).
"""
function probe_stage1(flat::Flat, decls::Vector{Decls}, tiers::Vector{Tier},
                      wss::Vector, mstores::Vector, Δt::Float64, ::Type{T}) where {T}
    map(zip(flat.paths, flat.comps, decls, tiers, wss, mstores)) do (path, c, d, t, ws, m)
        has_stage(h_x, c) || return NamedTuple()
        bn = bundle_names(h_x, c, t, ())
        y = h_x(c, _bundle_values(bn, d, NamedTuple(), NamedTuple(), T; ws, m, Δt))
        y isa NamedTuple ||
            throw(BuildError("`$path`: h_x must return a NamedTuple of port values, got $(typeof(y))"))
        _check_ports(path, "h_x", y, d.outs, T)
        _embed_ports(y, d.outs, T)
    end
end

function _bundle_values(bn, d::Decls, u, y_x, ::Type{T}; y = NamedTuple(), ws = nothing,
                        m = nothing, Δt = 0.0) where {T}
    vals = map(bn) do n
        n === :x   ? d.x :
        n === :m   ? m[] :
        n === :u   ? u :
        n === :y_x ? y_x :
        n === :y   ? y :
        n === :ws  ? ws :
        n === :t   ? zero(T) :
        n === :Δt  ? Δt : error("no probe source for $n")
    end
    NamedTuple{bn}(vals)
end

function _check_ports(path, stage, y::NamedTuple, outs::NamedTuple, ::Type{T}) where {T}
    for (name, v) in pairs(y)
        haskey(outs, name) ||
            throw(BuildError("`$path`: $stage returns `$name`, which `output_types` does not declare"))
        _accepts(outs[name], typeof(v), T) ||
            throw(BuildError("`$path`: $stage returns `$name`::$(typeof(v)), declared $(outs[name])" *
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
function schedule_stage2(flat::Flat, stage1::Vector)
    n = length(flat.comps)
    deps = [Int[] for _ in 1:n]
    for ci in 1:n
        has_stage(h_xu, flat.comps[ci]) || continue
        for (_, (ppath, pport)) in flat.conns[ci]
            isempty(ppath) && continue                   # a root slot: no producer to wait for
            pi = index_of(flat, ppath)
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
        names = join((flat.paths[ci] for ci in cycle), " → ")
        throw(BuildError("algebraic loop through stage-2 ports: $names — break it with a " *
                         "stage-1 (`h_x`) port, which carries no input dependence (§5.4/§5.5)"))
    end
    order
end

# --- 4. cell layout -----------------------------------------------------------
# Every declared port gets a cell; so does every root input face, the one
# terminal legitimately fed by no component, its initial value synthesized by
# `probe_value` (§6.1, §9.3, §11.3).
#
# An assembly's faces get no cells of their own. A face *is* its ultimate
# internal endpoint (§8.6), so it is entered as an alias onto that endpoint's
# address — which is what makes a face's type and tier derived rather than
# declared.

struct Layout
    addr::Dict{Tuple{String,Symbol},Any}     # (path, port|face) => CellAddr
    slots::Vector{Tuple{Symbol,Any}}         # root slots, with their probe values
    sizes::Vector{Pair{DataType,Int}}        # leaf eltype => buffer length, name-sorted
end

function cell_layout(flat::Flat, decls::Vector{Decls}, ::Type{T}) where {T}
    addr = Dict{Tuple{String,Symbol},Any}()
    slots = Tuple{Symbol,Any}[]
    offs = Dict{DataType,Int}()
    function place!(path, kind, name, P)
        L = _cell_leaf_eltype(path, kind, name, P)
        off = get(offs, L, 0)
        addr[(path, name)] = CellAddr{P}(off)
        offs[L] = off + nleaves(P)
    end
    for (path, d) in zip(flat.paths, decls)
        for (port, P) in pairs(d.outs)
            place!(path, "port", port, P)
        end
    end
    for face in flat.slots
        P = _slot_type(flat, decls, face)
        place!("", "root slot", face, P)
        push!(slots, (face, probe_value(P)))
    end
    for (alias, target) in flat.faces
        addr[alias] = addr[target]
    end
    sizes = sort!([L => n for (L, n) in offs]; by = p -> string(first(p)))
    Layout(addr, slots, sizes)
end

# A root slot's cell follows the same derivation an assembly face does: it is the
# internal endpoint's. With fan-out there are several, and the first one decides —
# a divergent sibling then meets the ordinary entry check at its own probe.
function _slot_type(flat::Flat, decls::Vector{Decls}, face::Symbol)
    for (ci, conns) in enumerate(flat.conns), (f, producer) in conns
        producer === ("", face) && return decls[ci].ins[f]
    end
    throw(BuildError("root input face `$face` routes to no input"))
end

# A cell lives whole in the buffer of its single leaf eltype, so an address is
# one offset. A mixed-leaf cell — legal in the design via pinning inside a
# declared struct (D-166) — needs multi-cursor addressing, which this increment
# deliberately does not build. Rejecting by name keeps the restriction a scope
# cut rather than a silent mislayout.
function _cell_leaf_eltype(path, kind, name, ::Type{P}) where {P}
    Ls = leaf_types(P)
    isempty(Ls) &&
        throw(BuildError("$(_at(path)): $kind `$name` declares $P, which has no leaves"))
    all(L -> L === Ls[1], Ls) ||
        throw(BuildError("$(_at(path)): $kind `$name` declares $P, a mixed-leaf cell — this " *
                         "increment lays out leaf-homogeneous cells only (multi-cursor " *
                         "addressing not built)"))
    Ls[1]
end

"""Address of the cell feeding `face`: its resolved producer's port, or a root slot."""
input_addr(layout::Layout, conns::Vector{Pair{Symbol,Tuple{String,Symbol}}}, face::Symbol) =
    layout.addr[last(conns[findfirst(p -> first(p) === face, conns)])]

# --- 5. the whole build -------------------------------------------------------

function build(root::AbstractComponent, ::Type{T} = Float64;
               chunk_size::Int = 16, Δt::Union{Nothing,Float64} = nothing) where {T}
    flat = flatten(root)
    tiers = [classify_tier(p, c) for (p, c) in zip(flat.paths, flat.comps)]
    decls = [declarations(c, t, T) for (c, t) in zip(flat.comps, tiers)]

    # `Δt` is bundle data carried in entry fields (§9.7), so the sample period
    # is fixed before the executor exists. Here one call binds deployment and
    # compilation together; in the design they are separate strata, and every
    # component sits at `D = 1, Φ = 0` until increment 5 brings the grid.
    if any(t -> t === DISCRETE, tiers) && Δt === nothing
        throw(BuildError("this model has discrete components; `build` needs the base tick " *
                         "period, as `build(root; Δt = …)` (§10.5)"))
    end
    Δt_base = something(Δt, 0.0)
    Δt_base >= 0 || throw(BuildError("Δt must be positive, got $Δt_base"))

    # Three homes for state, and no store mirrors another (§7.3): the flat
    # buffer for continuous `x`, one store per discrete `x`, one per mode set.
    # A store's *type* is shared by every instance of a component type, so
    # instances still compile to one body; only the reference varies.
    xstores = Any[t === DISCRETE && !isempty(d.x) ? Ref(d.x) : nothing
                  for (d, t) in zip(decls, tiers)]
    mstores = Any[isempty(init_m(c)) ? nothing : Ref(init_m(c)) for c in flat.comps]
    # Declaration by allocation (§7.3, D-077): sizes from the instance, eltypes
    # from the activation. Called once per activation, per component.
    wss = Any[_declares_workspace(c, t) ?
              (t === CONTINUOUS ? workspace(c, T) : workspace(c)) : nothing
              for (c, t) in zip(flat.comps, tiers)]

    stage1 = probe_stage1(flat, decls, tiers, wss, mstores, Δt_base, T)
    order = schedule_stage2(flat, stage1)
    layout = cell_layout(flat, decls, T)

    store = StoreBundle(NamedTuple{tuple((Symbol(L) for (L, _) in layout.sizes)...)}(
        tuple((CellStore(zeros(L, n)) for (L, n) in layout.sizes)...)))

    xbuf = T[]
    x_offs = Int[]
    for (d, t) in zip(decls, tiers)
        push!(x_offs, length(xbuf))
        t === CONTINUOUS && append!(xbuf, T[T(l) for l in _leaf_values(d.x)])
    end
    ẋbuf = zeros(T, length(xbuf))
    clock = Clock(zero(T))

    # root slots hold their synthesized values until a writer replaces them
    for (face, v) in layout.slots
        scatter!(store, layout.addr[("", face)], v)
    end

    addr_group(path, names) =
        NamedTuple{tuple(names...)}(tuple((layout.addr[(path, n)] for n in names)...))
    in_group(ci, d) =
        NamedTuple{tuple(keys(d.ins)...)}(
            tuple((input_addr(layout, flat.conns[ci], face) for face in keys(d.ins))...))

    hx_entries, hxu_entries, rhs_entries, tick_entries = Any[], Any[], Any[], Any[]
    hx_tiers, hxu_tiers, rhs_tiers, tick_tiers = Tier[], Tier[], Tier[], Tier[]
    products = NamedTuple[NamedTuple() for _ in flat.comps]   # probe products, per component

    # A discrete component's stages never run at a non-nominal activation: its
    # cells are frozen `Float64` constants with zero partials, holding what the
    # probe wrote, which is exactly what a tick at `t₀⁻` would have produced
    # (§7.2, §10.5; `frozen_discrete_walkthrough.md`).
    frozen(ci) = tiers[ci] === DISCRETE && T !== Float64

    # A frozen component's inputs are *synthesized*, not read from the upstream
    # product. Its stages never run at this activation, so it never gathers the
    # `Dual` cell its pinned input declaration would otherwise have to refuse —
    # the dependence through it is temporal, not instantaneous. Its own cells
    # then hold pinned constants, which downstream continuous consumers embed as
    # zero-partials. (The design carries the *nominal* activation's cell
    # contents across instead of synthesizing them: §9.4's activation story,
    # which this increment does not build.)
    in_values(ci, d) = NamedTuple{tuple(keys(d.ins)...)}(tuple(
        (frozen(ci) ? probe_value(d.ins[face]) :
         _probe_input(flat, layout, products, ci, face, d.ins[face], T)
         for face in keys(d.ins))...))

    for (ci, c) in enumerate(flat.comps)
        d, s1 = decls[ci], stage1[ci]
        products[ci] = s1
        (has_stage(h_x, c) && !frozen(ci)) || continue
        bn = bundle_names(h_x, c, tiers[ci], ())
        push!(hx_entries, StageEntry{typeof(d.x),bn}(
            h_x, c, NamedTuple(), NamedTuple(), addr_group(flat.paths[ci], keys(s1)),
            x_offs[ci], clock, xstores[ci], mstores[ci], wss[ci], Δt_base))
        push!(hx_tiers, tiers[ci])
    end

    for ci in order
        c, path, d, s1 = flat.comps[ci], flat.paths[ci], decls[ci], stage1[ci]
        has_stage(h_xu, c) || continue
        bn = bundle_names(h_xu, c, tiers[ci], tuple(keys(s1)...))
        u = in_values(ci, d)
        y2 = h_xu(c, _bundle_values(bn, d, u, s1, T; ws = wss[ci], m = mstores[ci],
                                    Δt = Δt_base))
        y2 isa NamedTuple ||
            throw(BuildError("`$path`: h_xu must return a NamedTuple, got $(typeof(y2))"))
        _check_ports(path, "h_xu", y2, d.outs, T)
        isempty(intersect(keys(s1), keys(y2))) ||
            throw(BuildError("`$path`: $(join(intersect(keys(s1), keys(y2)), ", ")) produced by two stages"))
        y2 = _embed_ports(y2, d.outs, T)
        products[ci] = merge(s1, y2)
        frozen(ci) && continue

        yx_g, out_g, in_g = addr_group(path, keys(s1)), addr_group(path, keys(y2)), in_group(ci, d)
        push!(hxu_entries, StageEntry{typeof(d.x),bn}(
            h_xu, c, in_g, yx_g, out_g, x_offs[ci], clock,
            xstores[ci], mstores[ci], wss[ci], Δt_base))
        push!(hxu_tiers, tiers[ci])
    end

    # Completeness of the declaration set (§8.2), for every component and not
    # only the stateful ones: an unproduced port would otherwise own a cell that
    # no stage ever writes, and read as a silent zero forever.
    for ci in eachindex(flat.comps)
        missing_ports = setdiff(keys(decls[ci].outs), keys(products[ci]))
        isempty(missing_ports) ||
            throw(BuildError("`$(flat.paths[ci])`: declared port(s) " *
                             "$(join(missing_ports, ", ")) produced by no stage"))
    end

    # The update law, one block per tier: `f` into the flat derivative buffer,
    # `g` into the component's own store. Both read the complete fresh table.
    for (ci, c) in enumerate(flat.comps)
        path, d, t = flat.paths[ci], decls[ci], tiers[ci]
        isempty(d.x) && continue
        update = t === CONTINUOUS ? f : g
        bn = bundle_names(update, c, t, tuple(keys(stage1[ci])...))
        u = in_values(ci, d)
        vals = _bundle_values(bn, d, u, stage1[ci], T; y = products[ci], ws = wss[ci],
                              m = mstores[ci], Δt = Δt_base)
        y_g, in_g = addr_group(path, keys(d.outs)), in_group(ci, d)

        if t === CONTINUOUS
            ẋ = f(c, vals)
            _check_derivative(path, ẋ, d.x)
            push!(rhs_entries, RHSEntry{typeof(d.x),bn}(
                c, in_g, y_g, x_offs[ci], clock, mstores[ci], wss[ci]))
            push!(rhs_tiers, t)
        else
            x⁺ = g(c, vals)
            _check_update(path, x⁺, d.x)
            frozen(ci) && continue
            push!(tick_entries, UpdateEntry{bn}(
                c, in_g, y_g, clock, xstores[ci], wss[ci], Δt_base))
            push!(tick_tiers, t)
        end
    end

    body(es, ts) = chunked_body(es, ts, store, xbuf, ẋbuf, clock; chunk_size)
    bodies = (sweep_hx = body(hx_entries, hx_tiers),
              sweep_hxu = body(hxu_entries, hxu_tiers),
              rhs = body(rhs_entries, rhs_tiers),
              ticks = body(tick_entries, tick_tiers))

    Simulation(store, xbuf, ẋbuf, clock, bodies, layout, flat, Δt_base, xstores, mstores, T)
end

# A probed input value: the producer's product, or the synthesized value of the
# root slot the obligation chain ends at. Stage-2 probing runs in topological
# order, so upstream products exist by construction.
#
# The entry is a *bound*, read permissively (D-167): a `T` entry is tolerant of
# both lawful arrivals, a pinned `Float64` entry demands a frozen one. The value
# passed on is the producer's, unembedded — the consumer gathers the producer's
# cell at runtime, so the cell's type is what its bundle carries.
function _probe_input(flat::Flat, layout::Layout, products, ci, face, P, ::Type{T}) where {T}
    path = flat.paths[ci]
    (ppath, pport) = last(flat.conns[ci][findfirst(p -> first(p) === face, flat.conns[ci])])
    v = if isempty(ppath)
        last(layout.slots[findfirst(s -> first(s) === pport, layout.slots)])
    else
        products[index_of(flat, ppath)][pport]
    end
    _accepts(P, typeof(v), T) ||
        throw(BuildError("`$path`.$face declared $P, fed from " *
                         "$(isempty(ppath) ? "root slot `$pport`" : "`$ppath`.$pport")" *
                         "::$(typeof(v))" * _pin_hint(P, typeof(v), T)))
    v
end

# §7.1: `Ẋ` has exactly `X`'s shape at the activation scalar. Checked
# structurally here so the runtime `flatten!` into the derivative block is safe.
function _check_derivative(path, ẋ, x::NamedTuple)
    ẋ isa NamedTuple ||
        throw(BuildError("`$path`: f must return a NamedTuple shaped like `init_x`, got $(typeof(ẋ))"))
    keys(ẋ) === keys(x) ||
        throw(BuildError("`$path`: f returns fields $(keys(ẋ)), state has $(keys(x)) — " *
                         "derivative completeness is structural (§7.1)"))
    for k in keys(x)
        nleaves(typeof(ẋ[k])) == nleaves(typeof(x[k])) ||
            throw(BuildError("`$path`: derivative field `$k` is $(typeof(ẋ[k])), state field is $(typeof(x[k]))"))
    end
end

# §7.3: a discrete store is overwritten wholesale with what `g` returns, so the
# successor must be the store's own type exactly. The discrete world is pinned —
# no walk, no embedding — which makes the store assignment type-stable and the
# ban on arithmetic over stores enforceable by construction.
function _check_update(path, x⁺, x::NamedTuple)
    x⁺ isa NamedTuple ||
        throw(BuildError("`$path`: g must return a NamedTuple shaped like `init_x`, got $(typeof(x⁺))"))
    typeof(x⁺) === typeof(x) ||
        throw(BuildError("`$path`: g returns $(typeof(x⁺)), state store is $(typeof(x)) — a " *
                         "discrete successor is the store's own type exactly (§7.3)"))
end
