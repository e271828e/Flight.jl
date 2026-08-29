# The build pipeline (§9.2, §9.3, §5.5). Plain printable data in, a compiled
# executor out. Three jobs, in order, because each needs the previous one's
# answer:
#
#   1. probe stage 1 — needs no wiring at all (an `h_x`/`h_s` bundle carries no `u`),
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
    x::NamedTuple          # continuous state, walked to the activation scalar
    s::NamedTuple          # discrete state, pinned wholesale (D-195)
    ins::NamedTuple        # face => type, evaluated at the activation scalar
    outs::NamedTuple       # port => type, evaluated at the activation scalar
end

"""The tier's own state register: exactly one of the two is ever populated."""
state_decls(d::Decls, t::Tier) = t === CONTINUOUS ? d.x : d.s

# The two registers split by D-166's criterion: `init_x` is by value and its
# types are *walked*; `input_types`/`output_types` are functions of the
# activation scalar on the continuous tier and are *evaluated*. There is no
# output-side leaf walk — the cell types at an activation are literally what the
# declaration returns at that `T`. On the discrete tier the plain forms declare
# the pinned world, so `init_s` does not walk and nothing is evaluated at `T`.
function declarations(c, t::Tier, ::Type{T}) where {T}
    t === CONTINUOUS ?
        Decls(retype_value(T, init_x(c)), NamedTuple(), input_types(c, T), output_types(c, T)) :
        Decls(NamedTuple(), init_s(c),
              declared_at(input_types, c, t), declared_at(output_types, c, t))
end

# --- tier classification (§8.2) -----------------------------------------------
# Tier is read off the declaration shape, never announced. For a **stateful**
# leaf the whole name family carries it — `init_x`/`h_x`/`h_xu`/`f` continuous
# against `init_s`/`h_s`/`h_su`/`g` discrete, the two disjoint (D-195), with the
# update law the decider. A **stateless** leaf has no update law, so its
# contract arities decide, `output_types` being mandatory and therefore the
# decider: the two-argument forms declare cells at the activation scalar, the
# plain forms the pinned discrete world (D-166/D-167). Every other tier-implying
# declaration must then agree, and disagreement names the offending one —
# including the wrong-letter cases the split state letters make visible, an
# `h_x` or an `init_x` on a leaf whose update law is `g` and the converse.
#
# The classifier sees primitives only: a component that declares nothing at all
# has no *class* to read, which §8.5 settles before this runs.

"""The tier the primitive at `path` announces, or a `BuildError` naming what disagrees."""
function classify_tier(path::String, c)
    votes = Tuple{Symbol,Tier}[]
    has_stage(f, c) && push!(votes, (:f, CONTINUOUS))
    has_stage(g, c) && push!(votes, (:g, DISCRETE))
    !isempty(init_x(c)) && push!(votes, (:init_x, CONTINUOUS))
    !isempty(init_s(c)) && push!(votes, (:init_s, DISCRETE))
    has_stage(h_x, c) && push!(votes, (:h_x, CONTINUOUS))
    has_stage(h_xu, c) && push!(votes, (:h_xu, CONTINUOUS))
    has_stage(h_s, c) && push!(votes, (:h_s, DISCRETE))
    has_stage(h_su, c) && push!(votes, (:h_su, DISCRETE))
    !isempty(init_m(c)) && push!(votes, (:init_m, CONTINUOUS))
    !isempty(events(c)) && push!(votes, (:events, CONTINUOUS))
    for (name, fn) in ((:output_types, output_types), (:input_types, input_types),
                       (:workspace, workspace))
        _declares(fn, c, Type{Float64}) && push!(votes, (name, CONTINUOUS))
        _declares(fn, c) && push!(votes, (name, DISCRETE))
    end

    # The decider, by §8.2's two cases.
    state = !isempty(init_x(c)) ? :init_x : !isempty(init_s(c)) ? :init_s : nothing
    if state !== nothing
        i = findfirst(v -> first(v) === :f || first(v) === :g, votes)
        i === nothing &&
            throw(BuildError(StoreWithoutUpdate(path = path, store = state)))
    else
        i = findfirst(v -> first(v) === :output_types, votes)
        i === nothing &&
            throw(BuildError(TierUnreadable(path = path,
                                            declarations = Symbol[first(v) for v in votes])))
    end

    # The vote loop collects (§13.1): a leaf written half in each tier's spelling
    # names every declaration that disagrees, not the first one found.
    t = last(votes[i])
    viol = Diagnostic[]
    for (name, vt) in votes
        vt === t ||
            push!(viol, DeclarationOnWrongTier(path = path, declaration = name,
                                               reason = :tier_form,
                                               found = Symbol(tier_word(vt)),
                                               announced = Symbol(tier_word(t))))
    end
    isempty(viol) || throw(BuildError(viol))
    t
end

# --- 2. probing ---------------------------------------------------------------

"""
Probe the stage-1 function — `h_x` continuous, `h_s` discrete (D-195) — for
every component in this activation's executable set: no wiring
is needed, so this runs first and tells the rest of the build which ports are
stage-1 — the ones that carry no dependence on inputs and therefore break loops
(§5.3). A frozen component's stages are outside the set (§9.4), so its stage-1
product is the nominal activation's, carried across from `carry`.

The `Δt` the discrete bundles carry is a fabricated, probe-scoped placeholder
(§9.3): `Δt` in seconds does not exist until `Simulation` binds `Δt_base`, and
the probe checks types, not physics.
"""
function probe_stage1(flat::Flat, decls::Vector{Decls}, tiers::Vector{Tier},
                      wss::Vector, mstores::Vector, carry, ::Type{T}) where {T}
    map(eachindex(flat.comps)) do ci
        path, c, d = flat.paths[ci], flat.comps[ci], decls[ci]
        _frozen(tiers, ci, T) && return carry.stage1[ci]
        h1 = stage1_of(tiers[ci])
        has_stage(h1, c) || return NamedTuple()
        stage = String(nameof(h1))
        bn = bundle_names(h1, c, tiers[ci], ())
        y = h1(c, _bundle_values(bn, d, NamedTuple(), NamedTuple(), T;
                                 ws = wss[ci], m = mstores[ci], Δt = 1.0))
        y isa NamedTuple ||
            throw(BuildError(ConformanceFailure(path = path, what = stage,
                                                reason = :return_type, shape = :ports,
                                                observed = typeof(y))))
        _check_ports(path, stage, y, d.outs, T)
        _embed_ports(y, d.outs, T)
    end
end

function _bundle_values(bn, d::Decls, u, y1, ::Type{T}; y = NamedTuple(), ws = nothing,
                        m = nothing, Δt = 0.0) where {T}
    vals = map(bn) do n
        n === :x   ? d.x :
        n === :s   ? d.s :
        n === :m   ? m[] :
        n === :u   ? u :
        n === :y_x ? y1 :
        n === :y_s ? y1 :
        n === :y   ? y :
        n === :ws  ? ws :
        n === :t   ? zero(T) :
        n === :Δt  ? Δt : throw(InternalInvariant("no probe source for $n"))
    end
    NamedTuple{bn}(vals)
end

# The per-port loop collects (§13.1): a stage returning three wrong types names
# all three.
function _check_ports(path, stage, y::NamedTuple, outs::NamedTuple, ::Type{T}) where {T}
    viol = Diagnostic[]
    for (name, v) in pairs(y)
        if !haskey(outs, name)
            push!(viol, UndeclaredReturnField(path = path, stage = stage, name = name,
                                              candidates = collect(keys(outs))))
            continue                               # nothing declared to compare against
        end
        _accepts(outs[name], typeof(v), T) ||
            push!(viol, ConformanceFailure(path = path, what = stage, reason = :field_type,
                                           shape = :ports, field = name,
                                           observed = typeof(v), declared = outs[name],
                                           activation = T))
    end
    isempty(viol) || throw(BuildError(viol))
    nothing
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
function schedule_stage2(flat::Flat, tiers::Vector{Tier}, stage1::Vector)
    n = length(flat.comps)
    deps = [Int[] for _ in 1:n]
    for ci in 1:n
        has_stage(stage2_of(tiers[ci]), flat.comps[ci]) || continue
        for (_, (ppath, pport)) in flat.conns[ci]
            isempty(ppath) && continue                   # a root input: no producer to wait for
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
        throw(BuildError(AlgebraicCycle(members = String[flat.paths[ci] for ci in cycle])))
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
    root_inputs::Vector{Tuple{Symbol,Any}}   # root inputs, with their probe values
    sizes::Vector{Pair{DataType,Int}}        # leaf eltype => buffer length, name-sorted
end

function cell_layout(flat::Flat, decls::Vector{Decls}, ::Type{T}) where {T}
    addr = Dict{Tuple{String,Symbol},Any}()
    root_inputs = Tuple{Symbol,Any}[]
    offs = Dict{DataType,Int}()
    # Placement collects (§13.1): every leafless declaration in the model is
    # named, and the barrier throws before the alias pass, which would otherwise
    # look up an address placement never made.
    viol = Diagnostic[]
    function place!(path, site::Symbol, name, ::Type{P}) where {P}
        lts = leaf_types(P)
        if isempty(lts)
            push!(viol, IllegalPortType(path = path, site = site, name = name, declared = P))
            return false
        end
        Ls = leaf_eltypes(P)
        addr[(path, name)] = CellAddr{P,length(Ls)}(Tuple(get(offs, L, 0) for L in Ls))
        for L in Ls
            offs[L] = get(offs, L, 0) + count(==(L), lts)
        end
        true
    end
    for (path, d) in zip(flat.paths, decls)
        for (port, P) in pairs(d.outs)
            place!(path, :port, port, P)
        end
    end
    for face in flat.root_inputs
        P = _root_input_type(flat, decls, face)
        place!("", :root_input, face, P) || continue
        push!(root_inputs, (face, probe_value(P)))
    end
    isempty(viol) || throw(BuildError(viol))
    for (alias, target) in flat.out_faces
        addr[alias] = addr[target]
    end
    sizes = sort!([L => n for (L, n) in offs]; by = p -> string(first(p)))
    Layout(addr, root_inputs, sizes)
end

# A root input's cell follows the same derivation an assembly face does: it is the
# internal endpoint's. With fan-out there are several, and the first one decides —
# a divergent sibling then meets the ordinary entry check at its own probe.
function _root_input_type(flat::Flat, decls::Vector{Decls}, face::Symbol)
    for (ci, conns) in enumerate(flat.conns), (f, producer) in conns
        producer === ("", face) && return decls[ci].ins[f]
    end
    throw(InternalInvariant("root input face `$face` routes to no input"))
end

"""Address of the cell feeding `face`: its resolved producer's port, or a root input."""
input_addr(layout::Layout, conns::Vector{Pair{Symbol,Tuple{String,Symbol}}}, face::Symbol) =
    layout.addr[last(conns[findfirst(p -> first(p) === face, conns)])]

# --- 5. the Build artifact and its activations (§9.2, §9.4) ---------------------

"""
One activation: Stratum C's products at a concrete scalar `T` (§9.1) —
declarations evaluated at `T`, the probe chain run over exactly the function
set this activation can execute (§9.4), cells laid out. Immutable once
constructed; the probe products are what a cell holds until first written
(§10.5).
"""
struct Activation{T}
    decls::Vector{Decls}
    stage1::Vector{NamedTuple}     # stage-1 probe products, per component
    products::Vector{NamedTuple}   # complete probe products, per component
    layout::Layout
end

"""
The deployment-free product of the build pipeline (§9.2): everything the strata
settle before `Δt_base` exists. Structure and schedule are `T`-independent by
construction (§9.1); the nominal `Float64` activation runs at build, and other
activations re-run Stratum C only, at first request, cached on the `Build`
(§9.4). The schedule it carries is anchor-relative — `flat.triples` against
`flat.anchors` — because final divisors for anchored entries do not exist until
`Δt_base` binds; one `Build` backs any number of `Simulation`s, each
materializing its own stores and buffers, so nothing writable lives here.
"""
struct Build
    flat::Flat
    tiers::Vector{Tier}
    order::Vector{Int}
    nominal::Activation{Float64}
    policies::Vector{NamedTuple}   # per component: event name => :boundary | :localized (§10.4)
    cache::Dict{DataType,Any}      # non-nominal activations, lazily materialized
end

"""
    build(root; activations = ()) → Build

Strata A and B plus the nominal activation (§9.1): flatten, classify, probe at
`Float64`, schedule, lay out. Nothing here needs `Δt_base`, `h` or `n` — those
are `Simulation`'s. `activations` is §9.4's opt-in exhaustive mode: each listed
scalar's activation is materialized eagerly instead of at first request.
"""
function build(root::AbstractComponent; activations::Tuple = ())
    flat = flatten(root)
    tiers = [classify_tier(p, c) for (p, c) in zip(flat.paths, flat.comps)]
    _check_event_declarations(flat)
    nominal, order = _stratum_c(flat, tiers, nothing, nothing, Float64)
    policies = probe_events(flat, tiers, nominal)
    b = Build(flat, tiers, order, nominal, policies, Dict{DataType,Any}())
    for A in activations
        activation(b, A)
    end
    b
end

# An event needs both halves (§8.2): a guard or handler with no method for the
# component type is caught by method lookup at declaration-reading time, rather
# than as a `MethodError` at the first firing — an event firing only in a corner
# of the envelope would otherwise hide the omission indefinitely.
function _check_event_declarations(flat::Flat)
    # One barrier for the whole pass (§13.1): every malformed entry in the model
    # is named, not the first one the walk reaches.
    viol = Diagnostic[]
    for (path, c) in zip(flat.paths, flat.comps)
        for (name, ev) in pairs(events(c))
            if !(ev isa Event)
                push!(viol, EventHalfMissing(path = path, event = name,
                                             reason = :not_an_event, found = typeof(ev)))
                continue                           # neither half exists to look up
            end
            for (half, fn) in ((:guard, ev.guard), (:handler, ev.handler))
                hasmethod(fn, Tuple{typeof(c),NamedTuple}) ||
                    push!(viol, EventHalfMissing(path = path, event = name, reason = half,
                                                 found = typeof(c)))
            end
        end
    end
    isempty(viol) || throw(BuildError(viol))
    nothing
end

"""
    activation(b, T) → Activation{T}

The activation at `T`: the nominal one directly, any other from the cache or by
a Stratum-C re-run at first request (§9.4). An activation is a pure function of
the build and the concrete scalar type, so caching is invisible; the spec's
torn-state guarantee is normative for concurrent first requests, which this
single-threaded prototype meets by having none.
"""
function activation(b::Build, ::Type{T}) where {T}
    T === Float64 && return b.nominal
    get!(b.cache, T) do
        first(_stratum_c(b.flat, b.tiers, b.order, b.nominal, T))
    end::Activation{T}
end

# Stratum C at `T`, parametric in the scalar (§9.1): declarations evaluated,
# probe chain run, cells laid out. At the nominal activation `carry` is
# `nothing`, every stage is probed (§9.3's probe-everything scope), and Stratum
# B's schedule falls out between the two probe passes — the stage-1 run at
# `Float64` serves classification and nominal products alike. A non-nominal
# activation receives both: the schedule is `T`-independent, and the frozen
# components' products are carried across from `carry` rather than probed,
# their stages being outside this activation's executable set (§9.4).
function _stratum_c(flat::Flat, tiers::Vector{Tier}, order, carry, ::Type{T}) where {T}
    decls = [declarations(c, t, T) for (c, t) in zip(flat.comps, tiers)]

    # Probe-scoped mode stores and workspaces (§9.3): the probes need `m` and
    # `ws` to build bundles, and everything these hold is garbage once the
    # build finishes — each `Simulation` materializes its own.
    mstores = Any[isempty(init_m(c)) ? nothing : Ref(init_m(c)) for c in flat.comps]
    wss = _workspaces(flat, tiers, T)

    stage1 = probe_stage1(flat, decls, tiers, wss, mstores, carry, T)
    order === nothing && (order = schedule_stage2(flat, tiers, stage1))
    layout = cell_layout(flat, decls, T)
    products = probe_stage2(flat, decls, tiers, stage1, order, layout, wss, mstores, carry, T)
    Activation{T}(decls, collect(stage1), products, layout), order
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
component. Frozen components are not probed: their complete products come from
`carry`, the nominal activation (§9.4). The discrete bundles' `Δt` is the
probe-scoped placeholder of `probe_stage1`. Returns the complete probe
products.
"""
function probe_stage2(flat::Flat, decls::Vector{Decls}, tiers::Vector{Tier},
                      stage1, order::Vector{Int}, layout::Layout,
                      wss::Vector, mstores::Vector, carry, ::Type{T}) where {T}
    products = NamedTuple[s1 for s1 in stage1]

    # A frozen component's stages never run at this activation, so its complete
    # product is the *nominal* activation's, carried across (§9.4): its cells
    # hold what a tick at `t₀⁻` computed from real nominal inputs — pinned
    # `Float64` constants, which downstream continuous consumers embed as
    # zero-partials. The dependence through it is temporal, not instantaneous
    # (`frozen_discrete_walkthrough.md`), so nothing here needs to gather the
    # `Dual` cell its pinned input declaration would otherwise have to refuse.
    for ci in eachindex(flat.comps)
        _frozen(tiers, ci, T) && (products[ci] = carry.products[ci])
    end

    in_values(ci, d) = NamedTuple{tuple(keys(d.ins)...)}(tuple(
        (_probe_input(flat, layout, products, ci, face, d.ins[face], T)
         for face in keys(d.ins))...))

    for ci in order
        c, path, d, s1 = flat.comps[ci], flat.paths[ci], decls[ci], stage1[ci]
        h2 = stage2_of(tiers[ci])
        (has_stage(h2, c) && !_frozen(tiers, ci, T)) || continue
        stage = String(nameof(h2))
        bn = bundle_names(h2, c, tiers[ci], tuple(keys(s1)...))
        u = in_values(ci, d)
        y2 = h2(c, _bundle_values(bn, d, u, s1, T; ws = wss[ci], m = mstores[ci],
                                  Δt = 1.0))
        y2 isa NamedTuple ||
            throw(BuildError(ConformanceFailure(path = path, what = stage,
                                                reason = :return_type, shape = :namedtuple,
                                                observed = typeof(y2))))
        _check_ports(path, stage, y2, d.outs, T)
        isempty(intersect(keys(s1), keys(y2))) ||
            throw(BuildError(ProducedByTwoStages(path = path,
                                                 ports = collect(intersect(keys(s1),
                                                                           keys(y2))))))
        products[ci] = merge(s1, _embed_ports(y2, d.outs, T))
    end

    # Completeness of the declaration set (§8.2), for every component and not
    # only the stateful ones: an unproduced port would otherwise own a cell that
    # no stage ever writes, and read as a silent zero forever. The pass collects
    # (§13.1): every component with an unproduced port is named, and the barrier
    # throws once for the whole model — as the two below it do.
    viol = Diagnostic[]
    for ci in eachindex(flat.comps)
        missing_ports = setdiff(keys(decls[ci].outs), keys(products[ci]))
        isempty(missing_ports) ||
            push!(viol, DeclaredNotProduced(path = flat.paths[ci],
                                            ports = collect(missing_ports),
                                            products = collect(keys(products[ci]))))
    end
    isempty(viol) || throw(BuildError(viol))

    # The update laws, probed against the now-complete table: `f` for shape,
    # `g` for the store's own type. A frozen component's `g` is outside the
    # executable set like its output stages (§9.4).
    empty!(viol)
    for (ci, c) in enumerate(flat.comps)
        path, d, t = flat.paths[ci], decls[ci], tiers[ci]
        (isempty(state_decls(d, t)) || _frozen(tiers, ci, T)) && continue
        update = update_of(t)
        bn = bundle_names(update, c, t, tuple(keys(stage1[ci])...))
        vals = _bundle_values(bn, d, in_values(ci, d), stage1[ci], T; y = products[ci],
                              ws = wss[ci], m = mstores[ci], Δt = 1.0)
        append!(viol, t === CONTINUOUS ? _check_derivative(path, f(c, vals), d.x) :
                                         _check_update(path, g(c, vals), d.s))
    end
    isempty(viol) || throw(BuildError(viol))

    # `project`, probed at every activation it runs at — its result is written
    # back to the buffer wholesale at both schedule positions (§5.3), so the
    # check holds it *complete* against `X`'s own shape at `T` (§9.3).
    empty!(viol)
    for (ci, c) in enumerate(flat.comps)
        has_stage(project, c) || continue
        path, d = flat.paths[ci], decls[ci]
        if tiers[ci] !== CONTINUOUS
            push!(viol, DeclarationOnWrongTier(path = path, declaration = :project,
                                               reason = :continuous_only))
            continue                               # no manifold to run it against
        end
        if isempty(d.x)
            push!(viol, DeclarationOnWrongTier(path = path, declaration = :project,
                                               reason = :no_manifold))
            continue
        end
        append!(viol, _check_state_write(path, "project", project(c, d.x), d.x))
    end
    isempty(viol) || throw(BuildError(viol))
    products
end

# The complete state write-back (§9.3): the same predicate for `project`'s
# return and a handler's `x` key, both written to the flat buffer wholesale.
#
# It returns its violation list rather than throwing, so its two callers — the
# `project` pass and the handler check — can put it under their own barrier
# (§13.1). The two shape checks are sequential: neither later one is meaningful
# once an earlier one fails.
function _check_state_write(path, what, x⁺, x::NamedTuple)
    x⁺ isa NamedTuple ||
        return Diagnostic[ConformanceFailure(path = path, what = what,
                                             reason = :return_type, shape = :state,
                                             observed = typeof(x⁺))]
    keys(x⁺) === keys(x) ||
        return Diagnostic[ConformanceFailure(path = path, what = what, reason = :field_set,
                                             shape = :state,
                                             observed_fields = collect(keys(x⁺)),
                                             declared_fields = collect(keys(x)))]
    viol = Diagnostic[]
    for k in keys(x)
        nleaves(typeof(x⁺[k])) == nleaves(typeof(x[k])) ||
            push!(viol, ConformanceFailure(path = path, what = what, reason = :field_type,
                                           shape = :state, field = k,
                                           observed = typeof(x⁺[k]),
                                           declared = typeof(x[k])))
    end
    viol
end

"""
Probe the event system, at the nominal activation only (§9.3, D-052): guards
and handlers never run at a non-nominal activation (§9.4), so nothing here is
parametric in the scalar. Each guard runs against real probed values and its
return type *is* the detection policy (§10.4, D-179) — `Bool` boundary-detected,
the nominal scalar localized, anything else an error naming both admissible
forms. Each handler runs once and its return is held to the §5.2 return law,
key by key. Returns the per-component policy register.
"""
function probe_events(flat::Flat, tiers::Vector{Tier}, act::Activation{Float64})
    decls, layout, products = act.decls, act.layout, act.products
    mstores = Any[isempty(init_m(c)) ? nothing : Ref(init_m(c)) for c in flat.comps]
    wss = _workspaces(flat, tiers, Float64)
    map(eachindex(flat.comps)) do ci
        c, path, d = flat.comps[ci], flat.paths[ci], decls[ci]
        evs = events(c)
        isempty(evs) && return NamedTuple()
        bn = event_bundle_names(c)
        u = NamedTuple{tuple(keys(d.ins)...)}(tuple(
            (_probe_input(flat, layout, products, ci, face, d.ins[face], Float64)
             for face in keys(d.ins))...))
        vals = _bundle_values(bn, d, u, NamedTuple(), Float64; y = products[ci],
                              ws = wss[ci], m = mstores[ci])
        NamedTuple{tuple(keys(evs)...)}(map(tuple(keys(evs)...)) do name
            σ = evs[name].guard(c, vals)
            policy = σ isa Bool ? :boundary :
                     σ isa Float64 ? :localized :
                     throw(BuildError(GuardForm(path = path, event = name,
                                                observed = typeof(σ))))
            _check_handler(path, name, evs[name].handler(c, vals), d, c)
            policy
        end)
    end
end

# The handler return law (§5.2, §9.3): a key is present iff the store exists on
# the component and the handler updates it. Key set first — an unknown key, or a
# key naming a store the component does not declare, names the stores that
# exist. Then per key: `x` complete (the flat buffer is written back wholesale),
# `m` a names-subset with matching types (per-field stores merge naturally).
#
# The key loops collect under one barrier per handler (§13.1): a handler naming
# three stores it does not own names all three.
function _check_handler(path, name, ret, d::Decls, c)
    what = "event `$name`'s handler"
    ret isa NamedTuple ||
        throw(BuildError(ConformanceFailure(path = path, what = what, reason = :return_type,
                                            shape = :stores, observed = typeof(ret))))
    m₀ = init_m(c)
    stores = Symbol[]
    isempty(d.x) || push!(stores, :x)
    isempty(m₀) || push!(stores, :m)
    viol = Diagnostic[]
    for k in keys(ret)
        k in stores ||
            push!(viol, HandlerReturnKey(path = path, event = name, key = k, stores = stores))
    end
    # Per key, but only for a store the component actually declares: the key set
    # is the outer fact, and holding a write to `x` against an empty state would
    # report the same omission twice in different words.
    haskey(ret, :x) && :x in stores &&
        append!(viol, _check_state_write(path, "$what `x`", ret.x, d.x))
    if haskey(ret, :m) && :m in stores
        if !(ret.m isa NamedTuple)
            push!(viol, ConformanceFailure(path = path, what = "$what `m`",
                                           reason = :return_type, shape = :mode,
                                           observed = typeof(ret.m)))
        else
            for k in keys(ret.m)
                if !haskey(m₀, k)
                    push!(viol, ConformanceFailure(path = path, what = what,
                                                   reason = :field_set, shape = :mode,
                                                   field = k,
                                                   declared_fields = collect(keys(m₀))))
                elseif typeof(ret.m[k]) !== typeof(m₀[k])
                    push!(viol, ConformanceFailure(path = path, what = what,
                                                   reason = :field_type, shape = :mode,
                                                   field = k, observed = typeof(ret.m[k]),
                                                   declared = typeof(m₀[k])))
                end
            end
        end
    end
    isempty(viol) || throw(BuildError(viol))
    nothing
end

# --- 6. deployment binding (§9.1) -----------------------------------------------
# Everything below post-dates the strata: it exists per `Simulation`, not per
# `Build`. Grid arithmetic is exact — GCD over `Rational{Int}` — and floats are
# refused at the door.

_exact(name::Symbol, v::Rational{Int}) = v
_exact(name::Symbol, v::Integer) = Rational{Int}(v)
_exact(name::Symbol, v::Period) = v.T
_exact(name::Symbol, v::AbstractFloat) =
    throw(BuildError(DeploymentInvalid(parameter = name, reason = :inexact, value = v)))
_exact(name::Symbol, v) =
    throw(BuildError(DeploymentInvalid(parameter = name, reason = :not_a_quantity,
                                       value = typeof(v))))

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
        throw(BuildError(DeploymentInvalid(parameter = :h, reason = :missing)))
    h_r = _exact(:h, h)
    h_r > 0 ||
        throw(BuildError(DeploymentInvalid(parameter = :h, reason = :range, value = h_r)))
    n === nothing || n ≥ 1 ||
        throw(BuildError(DeploymentInvalid(parameter = :n, reason = :range, value = n)))

    anchors, prov, triples = b.flat.anchors, b.flat.aprov, b.flat.triples
    # The constraint pool: every anchor's period and every nonzero offset (§9.1).
    pool = vcat([Tk for (Tk, _) in anchors], [τk for (_, τk) in anchors if τk != 0])

    if Δt_base === :derive
        unanchored = [b.flat.paths[ci] for ci in eachindex(b.tiers)
                      if b.tiers[ci] === DISCRETE && triples[ci][1] == 0]
        isempty(unanchored) ||
            throw(BuildError(DeploymentInvalid(parameter = :Δt_base, reason = :unanchored,
                                               paths = unanchored)))
        isempty(pool) &&
            throw(BuildError(DeploymentInvalid(parameter = :Δt_base,
                                               reason = :no_constraint)))
        Δt_r = reduce(gcd, pool)                     # the coarsest admissible value
    elseif Δt_base !== nothing
        Δt_r = _exact(:Δt_base, Δt_base)
    else
        Δt_r = something(n, 1) * h_r                 # the default path (§15.4)
    end

    n_i = _as_int(Δt_r / h_r)
    (n_i === nothing || n_i < 1) &&
        throw(BuildError(DeploymentInvalid(parameter = :Δt_base, reason = :not_harmonic,
                                           value = Δt_r, related = h_r)))
    n === nothing || n == n_i ||
        throw(BuildError(DeploymentInvalid(parameter = :Δt_base, reason = :disagrees_with_n,
                                           value = Δt_r, related = n, quotient = n_i)))

    # Per anchor, one exact division pair; anchor 0 is the base grid itself. The
    # loop collects (§13.1): every anchor the chosen base grid cannot express is
    # named, so the coarsest admissible value is chosen against the whole list.
    adm = isempty(pool) ? nothing : reduce(gcd, pool)
    viol = Diagnostic[]
    Dk, Φk = [1], [0]
    for (k, (Tk, τk)) in enumerate(anchors)
        D = _as_int(Tk / Δt_r)
        D === nothing &&
            push!(viol, DeploymentInvalid(parameter = :Δt_base, reason = :anchor_period,
                                          value = Tk, related = Δt_r, provenance = prov[k],
                                          admissible = adm))
        Φ = _as_int(τk / Δt_r)
        Φ === nothing &&
            push!(viol, DeploymentInvalid(parameter = :Δt_base, reason = :anchor_offset,
                                          value = τk, related = Δt_r, provenance = prov[k],
                                          admissible = adm))
        push!(Dk, something(D, 1)); push!(Φk, something(Φ, 0))
    end
    isempty(viol) || throw(BuildError(viol))

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

"""
The declared defaults, established into the three state homes (§7.3, §8.2):
`x` into the flat buffer by the declaration walk, `s` and `m` into the
component stores they were allocated for. The stores' *types* are fixed by
their allocation; this only ever writes values into them.

Construction and `init!` share this one path, which is what makes an
application "fresh run from the `init_*` defaults, with these overrides"
(D-063) rather than an overlay on whatever the last trajectory left behind.
"""
function establish_defaults!(xbuf::Vector{T}, sstores::Vector, mstores::Vector,
                             comps::Vector, decls::Vector{Decls},
                             tiers::Vector{Tier}) where {T}
    off = 0
    for (ci, (d, t)) in enumerate(zip(decls, tiers))
        if t === CONTINUOUS
            for l in _leaf_values(d.x)
                off += 1
                xbuf[off] = T(l)
            end
        end
        sstores[ci] === nothing || (sstores[ci][] = d.s)
        mstores[ci] === nothing || (mstores[ci][] = init_m(comps[ci]))
    end
    nothing
end

"""
One executor (§9.7, glossary): the compiled execution form of the schedule over
its own buffers, at one scalar — the entries' concretely-typed tuples closed
into the phase bodies, beside the store set they read and write.

Buffers are never cached, because every buffer set has exactly one owner
(§9.2): a `Simulation` owns its nominal executor, and every service invocation
instantiates its own from the same cached layouts.
"""
struct Executor{T,S,B,CL,EV}
    act::Activation{T}     # the layout and probe products it was materialized from
    store::S               # the signal table: cells and root inputs
    xbuf::Vector{T}        # continuous state, the flat buffer (§7.1)
    ẋbuf::Vector{T}        # its derivative — integrator scratch (§7.5)
    sstores::Vector{Any}   # discrete state stores, by component index
    mstores::Vector{Any}   # mode stores, by component index
    clock::CL
    bodies::B              # the phase bodies, closed over the buffers above
    events::EV             # the compiled event set, likewise (nominal only)
end

"""
The refusal every product compiled *against one activation* owes the executor
it is handed. A reader's `xbuf` offsets and store types and a plan's cell
addresses and converters are one activation's (§9.4), so running one against
another's buffer set would read and write the wrong slot in silence rather than
fail. The scalar rides in each product's own type and the executor carries it
too, so the pairing is dispatch — the fallback method, never a runtime test on
the hot path.

It carries no kind name, because it names no user-facing failure: §14.4 makes
the pairing a framework invariant the services uphold, neither plans nor
readers being user values, so reaching here is an internal assertion firing.
"""
@noinline _activation_mismatch(what::String, ::Type{T}, ::Type{S}) where {T,S} =
    throw(InternalInvariant(
        "this $what was compiled at $T and the executor " *
        "handed to it runs at $S — a service paired its products wrongly (§14.4). The " *
        "offsets, store types and cell addresses it bakes belong to one activation (§9.4), " *
        "and against another's buffer set they would read and write the wrong slot in " *
        "silence; one product is compiled per activation (§9.2)"))

function compile(b::Build, act::Activation{T}, D_c::Vector{Int}, Φ_c::Vector{Int},
                 Δt_c::Vector{Float64}; chunk_size::Int = 16) where {T}
    flat, tiers, decls, layout = b.flat, b.tiers, act.decls, act.layout

    # Three homes for state, and no store mirrors another (§7.3): the flat
    # buffer for continuous `x`, one store per discrete `s`, one per mode set.
    # A store's *type* is shared by every instance of a component type, so
    # instances still compile to one body; only the reference varies.
    sstores = Any[t === DISCRETE && !isempty(d.s) ? Ref(d.s) : nothing
                  for (d, t) in zip(decls, tiers)]
    mstores = Any[isempty(init_m(c)) ? nothing : Ref(init_m(c)) for c in flat.comps]
    wss = _workspaces(flat, tiers, T)

    store = StoreBundle(NamedTuple{tuple((_cell_key(L) for (L, _) in layout.sizes)...)}(
        tuple((CellStore(zeros(L, n)) for (L, n) in layout.sizes)...)))

    x_offs, nx = Int[], 0
    for (d, t) in zip(decls, tiers)
        push!(x_offs, nx)
        t === CONTINUOUS && (nx += nleaves(typeof(d.x)))
    end
    xbuf = zeros(T, nx)
    establish_defaults!(xbuf, sstores, mstores, flat.comps, decls, tiers)
    ẋbuf = zeros(T, nx)
    clock = Clock(zero(T))

    addr_group(path, names) =
        NamedTuple{tuple(names...)}(tuple((layout.addr[(path, n)] for n in names)...))
    in_group(ci, d) =
        NamedTuple{tuple(keys(d.ins)...)}(
            tuple((input_addr(layout, flat.conns[ci], face) for face in keys(d.ins))...))

    # A cell holds what the build probe populated until a sweep first writes it
    # (§10.5): this is the table's pre-`init!` content, and a frozen
    # component's pinned cells are this seed for the whole run (§9.4).
    for (ci, path) in enumerate(flat.paths)
        isempty(act.products[ci]) ||
            scatter_group!(store, addr_group(path, keys(act.products[ci])), act.products[ci])
    end
    # Root inputs hold their synthesized values until a writer replaces them.
    for (face, v) in layout.root_inputs
        scatter!(store, layout.addr[("", face)], v)
    end

    frozen(ci) = _frozen(tiers, ci, T)
    gate(ci) = tiers[ci] === DISCRETE ? (D_c[ci], Φ_c[ci]) : nothing
    y2keys(ci) = keys(act.products[ci])[length(keys(act.stage1[ci]))+1:end]

    stage1_entries, stage2_entries, rhs_entries, tick_entries = Any[], Any[], Any[], Any[]
    stage1_gates, stage2_gates, rhs_gates, tick_gates = Any[], Any[], Any[], Any[]

    for (ci, c) in enumerate(flat.comps)
        h1 = stage1_of(tiers[ci])
        (has_stage(h1, c) && !frozen(ci)) || continue
        d, s1 = decls[ci], act.stage1[ci]
        bn = bundle_names(h1, c, tiers[ci], ())
        push!(stage1_entries, StageEntry{typeof(d.x),bn}(
            h1, c, NamedTuple(), NamedTuple(), addr_group(flat.paths[ci], keys(s1)),
            x_offs[ci], clock, sstores[ci], mstores[ci], wss[ci], Δt_c[ci]))
        push!(stage1_gates, gate(ci))
    end

    for ci in b.order
        c, path, d, s1 = flat.comps[ci], flat.paths[ci], decls[ci], act.stage1[ci]
        h2 = stage2_of(tiers[ci])
        (has_stage(h2, c) && !frozen(ci)) || continue
        bn = bundle_names(h2, c, tiers[ci], tuple(keys(s1)...))
        push!(stage2_entries, StageEntry{typeof(d.x),bn}(
            h2, c, in_group(ci, d), addr_group(path, keys(s1)),
            addr_group(path, y2keys(ci)), x_offs[ci], clock,
            sstores[ci], mstores[ci], wss[ci], Δt_c[ci]))
        push!(stage2_gates, gate(ci))
    end

    # The update law, one block per tier: `f` into the flat derivative buffer,
    # `g` into the component's own store. Both read the complete fresh table.
    for (ci, c) in enumerate(flat.comps)
        path, d, t = flat.paths[ci], decls[ci], tiers[ci]
        (isempty(state_decls(d, t)) || frozen(ci)) && continue
        update = update_of(t)
        bn = bundle_names(update, c, t, tuple(keys(act.stage1[ci])...))
        y_g, in_g = addr_group(path, keys(d.outs)), in_group(ci, d)
        if t === CONTINUOUS
            push!(rhs_entries, RHSEntry{typeof(d.x),bn}(
                c, in_g, y_g, x_offs[ci], clock, mstores[ci], wss[ci]))
            push!(rhs_gates, nothing)
        else
            push!(tick_entries, UpdateEntry{bn}(
                c, in_g, y_g, clock, sstores[ci], wss[ci], Δt_c[ci]))
            push!(tick_gates, gate(ci))
        end
    end

    # The event system compiles at the nominal activation only: guards and
    # handlers are outside every other activation's executable set (§9.4,
    # D-052), so the event phase there is the bare sweep. Entries carry a global
    # index into the register vectors, in executor component order then
    # declaration order within a component (§10.6), and each carries its
    # `Build.policies` verdict into the compiled mask — which is what the frame
    # loop's trigger check reads (§10.4).
    ev_entries, ev_owner = Any[], Int[]
    ev_names = Tuple{String,Symbol}[]
    ev_localized = Bool[]
    if T === Float64
        for (ci, c) in enumerate(flat.comps)
            evs = events(c)
            isempty(evs) && continue
            d = decls[ci]
            bn = event_bundle_names(c)
            pj = has_stage(project, c) ? project : nothing
            for name in keys(evs)
                push!(ev_entries, EventEntry{typeof(d.x),bn}(
                    evs[name].guard, evs[name].handler, pj, c, length(ev_entries) + 1,
                    in_group(ci, d), addr_group(flat.paths[ci], keys(d.outs)),
                    x_offs[ci], clock, mstores[ci], wss[ci]))
                push!(ev_owner, ci)
                push!(ev_names, (flat.paths[ci], name))
                push!(ev_localized, b.policies[ci][name] === :localized)
            end
        end
    end

    # Projection runs at every activation — it is continuous machinery, inside
    # every executable set — between the integrate's state write and its decode
    # (§5.3).
    proj_entries = Any[ProjectEntry{typeof(decls[ci].x)}(c, x_offs[ci])
                       for (ci, c) in enumerate(flat.comps)
                       if tiers[ci] === CONTINUOUS && has_stage(project, c)]

    body(es, gs) = chunked_body(es, gs, store, xbuf, ẋbuf, clock; chunk_size)
    bodies = (sweep_1 = body(stage1_entries, stage1_gates),
              sweep_2 = body(stage2_entries, stage2_gates),
              rhs = body(rhs_entries, rhs_gates),
              ticks = body(tick_entries, tick_gates))

    evset = EventSet(ev_entries, proj_entries, store, xbuf, ev_owner, ev_names,
                     ev_localized, length(flat.comps))
    Executor(act, store, xbuf, ẋbuf, sstores, mstores, clock, bodies, evset)
end

# A probed input value: the producer's product, or the synthesized value of the
# root input the obligation chain ends at. Stage-2 probing runs in topological
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
        last(layout.root_inputs[findfirst(s -> first(s) === pport, layout.root_inputs)])
    else
        products[index_of(flat, ppath)][pport]
    end
    _accepts(P, typeof(v), T) ||
        throw(BuildError(WireTypeMismatch(path = path, face = face, declared = P,
                                          producer_path = ppath, producer_port = pport,
                                          observed = typeof(v), activation = T)))
    v
end

# §7.1: `Ẋ` has exactly `X`'s shape at the activation scalar. Checked
# structurally here so the runtime `flatten!` into the derivative block is safe.
function _check_derivative(path, ẋ, x::NamedTuple)
    ẋ isa NamedTuple ||
        return Diagnostic[ConformanceFailure(path = path, what = "f", reason = :return_type,
                                             shape = :init_x, observed = typeof(ẋ))]
    keys(ẋ) === keys(x) ||
        return Diagnostic[ConformanceFailure(path = path, what = "f", reason = :field_set,
                                             shape = :init_x,
                                             observed_fields = collect(keys(ẋ)),
                                             declared_fields = collect(keys(x)))]
    viol = Diagnostic[]
    for k in keys(x)
        nleaves(typeof(ẋ[k])) == nleaves(typeof(x[k])) ||
            push!(viol, ConformanceFailure(path = path, what = "f", reason = :field_type,
                                           shape = :init_x, field = k,
                                           observed = typeof(ẋ[k]), declared = typeof(x[k])))
    end
    viol
end

# §7.3: a discrete store is overwritten wholesale with what `g` returns, so the
# successor must be the store's own type exactly. The discrete world is pinned —
# no walk, no embedding — which makes the store assignment type-stable and the
# ban on arithmetic over stores enforceable by construction.
function _check_update(path, s⁺, s::NamedTuple)
    s⁺ isa NamedTuple ||
        return Diagnostic[ConformanceFailure(path = path, what = "g", reason = :return_type,
                                             shape = :init_s, observed = typeof(s⁺))]
    typeof(s⁺) === typeof(s) ||
        return Diagnostic[ConformanceFailure(path = path, what = "g", reason = :field_set,
                                             shape = :init_s, observed = typeof(s⁺),
                                             declared = typeof(s))]
    Diagnostic[]
end
