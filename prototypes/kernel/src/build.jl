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

The `Δt` the discrete bundles carry is a fabricated, probe-scoped placeholder
(§9.3): `Δt` in seconds does not exist until `Simulation` binds `Δt_base`, and
the probe checks types, not physics.
"""
function probe_stage1(flat::Flat, decls::Vector{Decls}, tiers::Vector{Tier},
                      wss::Vector, mstores::Vector, ::Type{T}) where {T}
    map(zip(flat.paths, flat.comps, decls, tiers, wss, mstores)) do (path, c, d, t, ws, m)
        has_stage(h_x, c) || return NamedTuple()
        bn = bundle_names(h_x, c, t, ())
        y = h_x(c, _bundle_values(bn, d, NamedTuple(), NamedTuple(), T; ws, m, Δt = 1.0))
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

# --- 5. the Build artifact (§9.2) -----------------------------------------------

"""
The deployment-free product of the build pipeline (§9.2): everything the strata
settle before `Δt_base` exists. The schedule it carries is anchor-relative —
`flat.triples` against `flat.anchors` — because final divisors for anchored
entries do not exist until `Δt_base` binds; one `Build` backs any number of
`Simulation`s, each materializing its own stores and buffers, so nothing
writable lives here. The probe products do: they are what a cell holds until
first written (§10.5).
"""
struct Build{T}
    flat::Flat
    tiers::Vector{Tier}
    decls::Vector{Decls}
    stage1::Vector{NamedTuple}     # stage-1 probe products, per component
    products::Vector{NamedTuple}   # complete probe products, per component
    order::Vector{Int}
    layout::Layout
end

"""
    build(root, T = Float64) → Build

The strata (§9.1), deployment-free: flatten, classify, probe, schedule, lay
out. Nothing here needs `Δt_base`, `h` or `n` — those are `Simulation`'s.
"""
function build(root::AbstractComponent, ::Type{T} = Float64) where {T}
    flat = flatten(root)
    tiers = [classify_tier(p, c) for (p, c) in zip(flat.paths, flat.comps)]
    decls = [declarations(c, t, T) for (c, t) in zip(flat.comps, tiers)]

    # Probe-scoped mode stores and workspaces (§9.3): the probes need `m` and
    # `ws` to build bundles, and everything these hold is garbage once the
    # build finishes — each `Simulation` materializes its own.
    mstores = Any[isempty(init_m(c)) ? nothing : Ref(init_m(c)) for c in flat.comps]
    wss = _workspaces(flat, tiers, T)

    stage1 = probe_stage1(flat, decls, tiers, wss, mstores, T)
    order = schedule_stage2(flat, stage1)
    layout = cell_layout(flat, decls, T)
    products = probe_stage2(flat, decls, tiers, stage1, order, layout, wss, mstores, T)
    Build{T}(flat, tiers, decls, collect(stage1), products, order, layout)
end

# Declaration by allocation (§7.3, D-077): sizes from the instance, eltypes from
# the activation. Called once per probe and once per `Simulation`.
_workspaces(flat::Flat, tiers::Vector{Tier}, ::Type{T}) where {T} =
    Any[_declares_workspace(c, t) ?
        (t === CONTINUOUS ? workspace(c, T) : workspace(c)) : nothing
        for (c, t) in zip(flat.comps, tiers)]

# A discrete component's stages never run at a non-nominal activation: its
# cells are frozen `Float64` constants with zero partials, holding what the
# probe wrote, which is exactly what a tick at `t₀⁻` would have produced
# (§7.2, §10.5; `frozen_discrete_walkthrough.md`).
_frozen(tiers::Vector{Tier}, ci::Int, ::Type{T}) where {T} =
    tiers[ci] === DISCRETE && T !== Float64

"""
Probe stage 2 in topological order — real upstream values available by
construction — and the update laws last, checking every return against the
declaration (§9.3); then the declaration-completeness check, for every
component. The discrete bundles' `Δt` is the probe-scoped placeholder of
`probe_stage1`. Returns the complete probe products.
"""
function probe_stage2(flat::Flat, decls::Vector{Decls}, tiers::Vector{Tier},
                      stage1, order::Vector{Int}, layout::Layout,
                      wss::Vector, mstores::Vector, ::Type{T}) where {T}
    products = NamedTuple[s1 for s1 in stage1]

    # A frozen component's inputs are *synthesized*, not read from the upstream
    # product. Its stages never run at this activation, so it never gathers the
    # `Dual` cell its pinned input declaration would otherwise have to refuse —
    # the dependence through it is temporal, not instantaneous. Its own cells
    # then hold pinned constants, which downstream continuous consumers embed as
    # zero-partials. (The design carries the *nominal* activation's cell
    # contents across instead of synthesizing them: §9.4's activation story,
    # which this increment does not build.)
    in_values(ci, d) = NamedTuple{tuple(keys(d.ins)...)}(tuple(
        (_frozen(tiers, ci, T) ? probe_value(d.ins[face]) :
         _probe_input(flat, layout, products, ci, face, d.ins[face], T)
         for face in keys(d.ins))...))

    for ci in order
        c, path, d, s1 = flat.comps[ci], flat.paths[ci], decls[ci], stage1[ci]
        has_stage(h_xu, c) || continue
        bn = bundle_names(h_xu, c, tiers[ci], tuple(keys(s1)...))
        u = in_values(ci, d)
        y2 = h_xu(c, _bundle_values(bn, d, u, s1, T; ws = wss[ci], m = mstores[ci],
                                    Δt = 1.0))
        y2 isa NamedTuple ||
            throw(BuildError("`$path`: h_xu must return a NamedTuple, got $(typeof(y2))"))
        _check_ports(path, "h_xu", y2, d.outs, T)
        isempty(intersect(keys(s1), keys(y2))) ||
            throw(BuildError("`$path`: $(join(intersect(keys(s1), keys(y2)), ", ")) produced by two stages"))
        products[ci] = merge(s1, _embed_ports(y2, d.outs, T))
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

    # The update laws, probed against the now-complete table: `f` for shape,
    # `g` for the store's own type.
    for (ci, c) in enumerate(flat.comps)
        path, d, t = flat.paths[ci], decls[ci], tiers[ci]
        isempty(d.x) && continue
        update = t === CONTINUOUS ? f : g
        bn = bundle_names(update, c, t, tuple(keys(stage1[ci])...))
        vals = _bundle_values(bn, d, in_values(ci, d), stage1[ci], T; y = products[ci],
                              ws = wss[ci], m = mstores[ci], Δt = 1.0)
        t === CONTINUOUS ? _check_derivative(path, f(c, vals), d.x) :
                           _check_update(path, g(c, vals), d.x)
    end
    products
end

# --- 6. deployment binding (§9.1) -----------------------------------------------
# Everything below post-dates the strata: it exists per `Simulation`, not per
# `Build`. Grid arithmetic is exact — GCD over `Rational{Int}` — and floats are
# refused at the door.

_exact(name, v::Rational{Int}) = v
_exact(name, v::Integer) = Rational{Int}(v)
_exact(name, v::Period) = v.T
_exact(name, v::AbstractFloat) =
    throw(BuildError("`$name` must be exact — a Rational (`$name = 1//100`) or a " *
                     "`Period`/`Hz` value: grid derivation is GCD arithmetic, ill-defined " *
                     "over floats (§9.1, §10.5)"))
_exact(name, v) =
    throw(BuildError("`$name` must be a Rational or a `Period`/`Hz` value, got $(typeof(v))"))

_as_int(r::Rational) = denominator(r) == 1 ? Int(numerator(r)) : nothing

"""
Deployment binding (§9.1): `Δt_base` from one of three cross-validated sources —
the explicit keyword, the `n·h` product (default `n = 1`), or GCD derivation
over the constraint pool, requested as `Δt_base = :derive` and permitted only
with every discrete component anchored. Resolution is one exact division pair
per anchor and one multiply-add per component. Returns the bound deployment:
`h`, `n`, `Δt_base`, the per-component `(D, Φ, Δt)` columns, and the bound
schedule (§9.2's printable artifact, as plain data).
"""
function bind_schedule(b::Build, h, n, Δt_base)
    h === nothing &&
        throw(BuildError("deployment needs the continuous step: `Simulation(…; h = 1//100)` " *
                         "— a domain rate is not a framework default (§9.1)"))
    h_r = _exact("h", h)
    h_r > 0 || throw(BuildError("h must be positive, got $h_r"))
    n === nothing || n ≥ 1 || throw(BuildError("n must be an integer ≥ 1, got $n (§9.1)"))

    anchors, prov, triples = b.flat.anchors, b.flat.aprov, b.flat.triples
    # The constraint pool: every anchor's period and every nonzero offset (§9.1).
    pool = vcat([Tk for (Tk, _) in anchors], [τk for (_, τk) in anchors if τk != 0])

    if Δt_base === :derive
        unanchored = [b.flat.paths[ci] for ci in eachindex(b.tiers)
                      if b.tiers[ci] === DISCRETE && triples[ci][1] == 0]
        isempty(unanchored) ||
            throw(BuildError("Δt_base cannot be derived: `$(join(unanchored, "`, `"))` " *
                             "is/are unanchored, with period `m·Δt_base` — an anchor edit " *
                             "anywhere in the tree would silently rescale it. Declare the " *
                             "base tick period instead: `Δt_base = …`, or `n = …` (§9.1)"))
        isempty(pool) &&
            throw(BuildError("Δt_base cannot be derived: no anchor declares a constraint " *
                             "to derive it from (§9.1)"))
        Δt_r = reduce(gcd, pool)                     # the coarsest admissible value
    elseif Δt_base !== nothing
        Δt_r = _exact("Δt_base", Δt_base)
    else
        Δt_r = something(n, 1) * h_r                 # the default path (§15.4)
    end

    n_i = _as_int(Δt_r / h_r)
    (n_i === nothing || n_i < 1) &&
        throw(BuildError("harmonic grid: Δt_base = $Δt_r is not an integer multiple of " *
                         "h = $h_r (Δt_base = n·h, n ≥ 1, §10.5)"))
    n === nothing || n == n_i ||
        throw(BuildError("Δt_base = $Δt_r disagrees with n = $n: Δt_base/h = $n_i (§9.1)"))

    # Per anchor, one exact division pair; anchor 0 is the base grid itself.
    admissible() = "an admissible Δt_base divides gcd(pool) = $(reduce(gcd, pool))"
    Dk, Φk = [1], [0]
    for (k, (Tk, τk)) in enumerate(anchors)
        D = _as_int(Tk / Δt_r)
        D === nothing &&
            throw(BuildError("$(prov[k]): period $Tk is not an integer multiple of " *
                             "Δt_base = $Δt_r — $(admissible()) (§9.1)"))
        Φ = _as_int(τk / Δt_r)
        Φ === nothing &&
            throw(BuildError("$(prov[k]): offset $τk does not land on the base grid at " *
                             "Δt_base = $Δt_r — $(admissible()) (§9.1)"))
        push!(Dk, D); push!(Φk, Φ)
    end

    # Per component, one multiply-add; the canonical residue 0 ≤ Φ < D survives
    # composition (§10.5), which is what the gate's truncated rem relies on.
    Δtb = Float64(Δt_r)
    D_c, Φ_c, Δt_c = Int[], Int[], Float64[]
    sched = @NamedTuple{path::String, D::Int, Φ::Int, Δt::Float64}[]
    for ci in eachindex(b.tiers)
        (a, m, c) = triples[ci]
        if b.tiers[ci] === DISCRETE
            D, Φ = m * Dk[a + 1], Φk[a + 1] + c * Dk[a + 1]
            push!(D_c, D); push!(Φ_c, Φ); push!(Δt_c, D * Δtb)
            push!(sched, (path = b.flat.paths[ci], D = D, Φ = Φ, Δt = D * Δtb))
        else
            push!(D_c, 1); push!(Φ_c, 0); push!(Δt_c, 0.0)
        end
    end
    (h = Float64(h_r), n = n_i, Δt_base = Δtb, sched = sched, D = D_c, Φ = Φ_c, Δt = Δt_c)
end

# --- 7. entry compilation, per deployment ---------------------------------------
# What moves behind deployment binding: `Δt`, `D` and `Φ` are entry data (§9.7),
# so the executor cannot exist until `Δt_base` does. Each call materializes its
# own stores and buffers — the `Build` stays immutable and backs many
# `Simulation`s (§9.2). No user stage runs here: every check already ran at the
# probe, and what compiles is the checked shape.

function compile(b::Build{T}, D_c::Vector{Int}, Φ_c::Vector{Int}, Δt_c::Vector{Float64};
                 chunk_size::Int = 16) where {T}
    flat, tiers, decls, layout = b.flat, b.tiers, b.decls, b.layout

    # Three homes for state, and no store mirrors another (§7.3): the flat
    # buffer for continuous `x`, one store per discrete `x`, one per mode set.
    # A store's *type* is shared by every instance of a component type, so
    # instances still compile to one body; only the reference varies.
    xstores = Any[t === DISCRETE && !isempty(d.x) ? Ref(d.x) : nothing
                  for (d, t) in zip(decls, tiers)]
    mstores = Any[isempty(init_m(c)) ? nothing : Ref(init_m(c)) for c in flat.comps]
    wss = _workspaces(flat, tiers, T)

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

    addr_group(path, names) =
        NamedTuple{tuple(names...)}(tuple((layout.addr[(path, n)] for n in names)...))
    in_group(ci, d) =
        NamedTuple{tuple(keys(d.ins)...)}(
            tuple((input_addr(layout, flat.conns[ci], face) for face in keys(d.ins))...))

    # A cell holds what the build probe populated until a sweep first writes it
    # (§10.5): an offset component's pre-first-tick reads and a frozen
    # component's pinned cells are both this seed.
    for (ci, path) in enumerate(flat.paths)
        isempty(b.products[ci]) ||
            scatter_group!(store, addr_group(path, keys(b.products[ci])), b.products[ci])
    end
    # Root slots hold their synthesized values until a writer replaces them.
    for (face, v) in layout.slots
        scatter!(store, layout.addr[("", face)], v)
    end

    frozen(ci) = _frozen(tiers, ci, T)
    gate(ci) = tiers[ci] === DISCRETE ? (D_c[ci], Φ_c[ci]) : nothing
    y2keys(ci) = keys(b.products[ci])[length(keys(b.stage1[ci]))+1:end]

    hx_entries, hxu_entries, rhs_entries, tick_entries = Any[], Any[], Any[], Any[]
    hx_gates, hxu_gates, rhs_gates, tick_gates = Any[], Any[], Any[], Any[]

    for (ci, c) in enumerate(flat.comps)
        (has_stage(h_x, c) && !frozen(ci)) || continue
        d, s1 = decls[ci], b.stage1[ci]
        bn = bundle_names(h_x, c, tiers[ci], ())
        push!(hx_entries, StageEntry{typeof(d.x),bn}(
            h_x, c, NamedTuple(), NamedTuple(), addr_group(flat.paths[ci], keys(s1)),
            x_offs[ci], clock, xstores[ci], mstores[ci], wss[ci], Δt_c[ci]))
        push!(hx_gates, gate(ci))
    end

    for ci in b.order
        c, path, d, s1 = flat.comps[ci], flat.paths[ci], decls[ci], b.stage1[ci]
        (has_stage(h_xu, c) && !frozen(ci)) || continue
        bn = bundle_names(h_xu, c, tiers[ci], tuple(keys(s1)...))
        push!(hxu_entries, StageEntry{typeof(d.x),bn}(
            h_xu, c, in_group(ci, d), addr_group(path, keys(s1)),
            addr_group(path, y2keys(ci)), x_offs[ci], clock,
            xstores[ci], mstores[ci], wss[ci], Δt_c[ci]))
        push!(hxu_gates, gate(ci))
    end

    # The update law, one block per tier: `f` into the flat derivative buffer,
    # `g` into the component's own store. Both read the complete fresh table.
    for (ci, c) in enumerate(flat.comps)
        path, d, t = flat.paths[ci], decls[ci], tiers[ci]
        (isempty(d.x) || frozen(ci)) && continue
        update = t === CONTINUOUS ? f : g
        bn = bundle_names(update, c, t, tuple(keys(b.stage1[ci])...))
        y_g, in_g = addr_group(path, keys(d.outs)), in_group(ci, d)
        if t === CONTINUOUS
            push!(rhs_entries, RHSEntry{typeof(d.x),bn}(
                c, in_g, y_g, x_offs[ci], clock, mstores[ci], wss[ci]))
            push!(rhs_gates, nothing)
        else
            push!(tick_entries, UpdateEntry{bn}(
                c, in_g, y_g, clock, xstores[ci], wss[ci], Δt_c[ci]))
            push!(tick_gates, gate(ci))
        end
    end

    body(es, gs) = chunked_body(es, gs, store, xbuf, ẋbuf, clock; chunk_size)
    bodies = (sweep_hx = body(hx_entries, hx_gates),
              sweep_hxu = body(hxu_entries, hxu_gates),
              rhs = body(rhs_entries, rhs_gates),
              ticks = body(tick_entries, tick_gates))

    (store = store, xbuf = xbuf, ẋbuf = ẋbuf, clock = clock, bodies = bodies,
     xstores = xstores, mstores = mstores)
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
