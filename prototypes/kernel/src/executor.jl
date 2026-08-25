# The compiled executor (§9.7): schedule entries carrying what selects code in
# type parameters and what is plain data in fields, gathered into chunked
# statically-typed tuples behind non-inlined barriers, traversed by a
# compile-time-unrolled walk.

# --- entries ------------------------------------------------------------------
# Three kinds, by where the product goes: a stage entry writes cells (both
# tiers — one entry type carries either tier's stage function, whose names are
# disjoint, D-195), an RHS entry writes the flat `ẋ` buffer, an update entry
# writes its own discrete state store. All build their §5.2 bundle the same way,
# from `BN` — the bundle name set the law fixed at build time.
#
# What selects code sits in type parameters, what varies per instance in fields
# (§9.7): the state store is a `Ref` whose *type* is shared by every instance of
# a component type, so two instances still compile to one body.

struct StageEntry{F,Comp,XT,BN,IA<:NamedTuple,YA<:NamedTuple,OA<:NamedTuple,CL,SS,MS,WS}
    fn::F
    comp::Comp
    inputs::IA      # input face => cell address (the wiring's name binding)
    y1::YA          # own stage-1 port => cell address (`y_x`, `y_s` on the discrete tier)
    outs::OA        # port this entry writes => cell address
    x_off::Int      # continuous state offset into the flat buffer
    clock::CL
    sstore::SS      # discrete state store, or nothing on the continuous tier
    mstore::MS      # mode store, or nothing
    ws::WS          # workspace, or nothing
    Δt::Float64     # sample period; unused on the continuous tier
end

struct RHSEntry{Comp,XT,BN,IA<:NamedTuple,YA<:NamedTuple,CL,MS,WS}
    comp::Comp
    inputs::IA
    y::YA           # every own port — `f` reads the complete fresh table (§5.3)
    x_off::Int
    clock::CL
    mstore::MS
    ws::WS
end

struct UpdateEntry{Comp,BN,IA<:NamedTuple,YA<:NamedTuple,CL,SS,WS}
    comp::Comp
    inputs::IA
    y::YA           # every own port — `g` reads the complete fresh table too
    clock::CL
    sstore::SS      # written by this entry, and by nothing else
    ws::WS
    Δt::Float64
end

# Outer constructors: only `XT`/`BN` cannot be inferred from the arguments.
StageEntry{XT,BN}(fn, comp, inputs, y1, outs, x_off, clock, sstore, mstore, ws, Δt) where {XT,BN} =
    StageEntry{typeof(fn),typeof(comp),XT,BN,typeof(inputs),typeof(y1),typeof(outs),
               typeof(clock),typeof(sstore),typeof(mstore),typeof(ws)}(
        fn, comp, inputs, y1, outs, x_off, clock, sstore, mstore, ws, Δt)

RHSEntry{XT,BN}(comp, inputs, y, x_off, clock, mstore, ws) where {XT,BN} =
    RHSEntry{typeof(comp),XT,BN,typeof(inputs),typeof(y),typeof(clock),
             typeof(mstore),typeof(ws)}(comp, inputs, y, x_off, clock, mstore, ws)

UpdateEntry{BN}(comp, inputs, y, clock, sstore, ws, Δt) where {BN} =
    UpdateEntry{typeof(comp),BN,typeof(inputs),typeof(y),typeof(clock),
                typeof(sstore),typeof(ws)}(comp, inputs, y, clock, sstore, ws, Δt)

# One bundle-expression builder, three @generated entry points. Absent names are
# absent, never `nothing`-filled: a body destructuring what it does not own
# fails at the destructuring (§5.2). The two state letters name the two homes
# outright (D-195): `x` reconstructs from the flat buffer, `s` reads the
# component's own store, and no name selects a home by tier at compile time
# because the tier already picked the name.
function _bundle_expr(BN, XT)
    args = map(BN) do n
        n === :x   ? :(reconstruct($XT, xbuf, e.x_off)) :
        n === :s   ? :(e.sstore[]) :
        n === :m   ? :(e.mstore[]) :
        n === :u   ? :(gather_group(e.inputs, store)) :
        n === :y_x ? :(gather_group(e.y1, store)) :
        n === :y_s ? :(gather_group(e.y1, store)) :
        n === :y   ? :(gather_group(e.y, store)) :
        n === :ws  ? :(e.ws) :
        n === :t   ? :(e.clock.t) :
        n === :Δt  ? :(e.Δt) :
        error("no source for bundle field $n")
    end
    :(NamedTuple{$BN}(($(args...),)))
end

@generated function make_bundle(e::StageEntry{F,Comp,XT,BN}, store,
                                xbuf) where {F,Comp,XT,BN}
    quote
        $(Expr(:meta, :inline))
        $(_bundle_expr(BN, XT))
    end
end

@generated function make_bundle(e::RHSEntry{Comp,XT,BN}, store, xbuf) where {Comp,XT,BN}
    quote
        $(Expr(:meta, :inline))
        $(_bundle_expr(BN, XT))
    end
end

@generated function make_bundle(e::UpdateEntry{Comp,BN}, store, xbuf) where {Comp,BN}
    quote
        $(Expr(:meta, :inline))
        $(_bundle_expr(BN, Nothing))
    end
end

@inline function run!(e::StageEntry, store, xbuf, ẋbuf)
    y = e.fn(e.comp, make_bundle(e, store, xbuf))
    scatter_group!(store, e.outs, y)
end

@inline function run!(e::RHSEntry, store, xbuf, ẋbuf)
    ẋ = f(e.comp, make_bundle(e, store, xbuf))
    flatten!(ẋbuf, e.x_off, ẋ)      # shape conformance established at probe time
    nothing
end

# The jump map: `g` reads the fresh table and writes only its own store, which
# is what makes the update block order-free with disjoint writes (§9.7).
@inline function run!(e::UpdateEntry, store, xbuf, ẋbuf)
    e.sstore[] = g(e.comp, make_bundle(e, store, xbuf))
    nothing
end

# --- the event machinery (§10.6, §5.3) ------------------------------------------
# Not sweep entries: guards and handlers are driven by the boundary iteration in
# `sim.jl`, against per-event registers, so their entries live in their own
# compiled set. One entry per declared event, carrying both halves plus its
# component's `project` (or `nothing`), and a global index into the register
# vectors — global order is executor component order, then declaration order
# within a component, which is what makes the §13.4-style dispatch order
# deterministic. Bundles are built exactly like every other entry's, from the
# event name set the law fixed at build time (`x, m, u, y, ws, t`, §5.2).

struct EventEntry{G,H,P,Comp,XT,BN,IA<:NamedTuple,YA<:NamedTuple,CL,MS,WS}
    guard::G
    handler::H
    proj::P         # the component's `project`, or nothing
    comp::Comp
    idx::Int        # global event index into the register vectors
    inputs::IA
    y::YA           # every own port — guards and handlers read the complete fresh table
    x_off::Int
    clock::CL
    mstore::MS
    ws::WS
end

EventEntry{XT,BN}(guard, handler, proj, comp, idx, inputs, y, x_off, clock,
                  mstore, ws) where {XT,BN} =
    EventEntry{typeof(guard),typeof(handler),typeof(proj),typeof(comp),XT,BN,
               typeof(inputs),typeof(y),typeof(clock),typeof(mstore),typeof(ws)}(
        guard, handler, proj, comp, idx, inputs, y, x_off, clock, mstore, ws)

@generated function make_bundle(e::EventEntry{G,H,P,Comp,XT,BN}, store,
                                xbuf) where {G,H,P,Comp,XT,BN}
    quote
        $(Expr(:meta, :inline))
        $(_bundle_expr(BN, XT))
    end
end

"""
Projection between a state write and its decode (§5.3): reconstruct, project,
write back wholesale — the completeness the probe established is what makes the
wholesale write safe by construction.
"""
struct ProjectEntry{Comp,XT}
    comp::Comp
    x_off::Int
end
ProjectEntry{XT}(comp, x_off) where {XT} = ProjectEntry{typeof(comp),XT}(comp, x_off)

@inline function run_project!(e::ProjectEntry{Comp,XT}, xbuf) where {Comp,XT}
    flatten!(xbuf, e.x_off, project(e.comp, reconstruct(XT, xbuf, e.x_off)))
    nothing
end

"""
The per-`Simulation` compiled event set: the entry tuples plus the §10.6
registers. The three normative registers — prior, last-observed sample, firing
count — are detection bookkeeping, not model memory: plain vectors indexed by
the global event index, in no state store, reconstructed deterministically.
`now`, `fire` and `comp_fired` are the iteration's round-scoped scratch, and
`warned` implements "at most once per event per boundary" for the
`FiringBudget` degradation.
"""
struct EventSet{E<:Tuple,P<:Tuple,S,X}
    entries::E
    projects::P
    store::S
    xbuf::X
    owner::Vector{Int}                   # component index per event
    names::Vector{Tuple{String,Symbol}}  # (path, event name), for the degradation warning
    now::Vector{Bool}
    prior::Vector{Bool}
    last::Vector{Bool}
    fire::Vector{Bool}
    count::Vector{Int}
    warned::Vector{Bool}
    comp_fired::Vector{Bool}             # per component, round-scoped
end

function EventSet(entries::Vector, projects::Vector, store, xbuf,
                  owner::Vector{Int}, names::Vector{Tuple{String,Symbol}}, ncomps::Int)
    n = length(entries)
    EventSet(tuple(entries...), tuple(projects...), store, xbuf, owner, names,
             fill(false, n), fill(false, n), fill(false, n), fill(false, n),
             zeros(Int, n), fill(false, n), fill(false, ncomps))
end

# The three walks the iteration drives, each the compile-time-unrolled tuple
# recursion of the phase bodies. Guard evaluation writes each predicate sample
# into `now` by global index; the fire walk runs `handler → project` for
# exactly the masked entries, latching the returned stores — `x` into the flat
# buffer, `m` merged into the mode store, per the return law's iff shape (§5.2).

@noinline _projects!(es::EventSet) = _proj_walk(es.projects, es.xbuf)
@inline _proj_walk(::Tuple{}, xbuf) = nothing
@inline function _proj_walk(t::Tuple, xbuf)
    run_project!(t[1], xbuf)
    _proj_walk(Base.tail(t), xbuf)
end

@noinline _guards!(es::EventSet) = _guard_walk(es.entries, es.store, es.xbuf, es.now)
@inline _guard_walk(::Tuple{}, store, xbuf, now) = nothing
@inline function _guard_walk(t::Tuple, store, xbuf, now)
    e = t[1]
    now[e.idx] = _holding(e.guard(e.comp, make_bundle(e, store, xbuf)))
    _guard_walk(Base.tail(t), store, xbuf, now)
end

@noinline _fire!(es::EventSet) = _fire_walk(es.entries, es.store, es.xbuf, es.fire)
@inline _fire_walk(::Tuple{}, store, xbuf, fire) = nothing
@inline function _fire_walk(t::Tuple, store, xbuf, fire)
    e = t[1]
    if fire[e.idx]
        _latch!(e, e.handler(e.comp, make_bundle(e, store, xbuf)), xbuf)
        _fire_project!(e, xbuf)
    end
    _fire_walk(Base.tail(t), store, xbuf, fire)
end

@inline function _latch!(e::EventEntry, ret::NamedTuple, xbuf)
    haskey(ret, :x) && flatten!(xbuf, e.x_off, ret.x)
    haskey(ret, :m) && (e.mstore[] = merge(e.mstore[], ret.m))
    nothing
end

@inline function _fire_project!(e::EventEntry{G,H,P,Comp,XT}, xbuf) where {G,H,P,Comp,XT}
    P === Nothing && return nothing
    flatten!(xbuf, e.x_off, e.proj(e.comp, reconstruct(XT, xbuf, e.x_off)))
    nothing
end

# --- the gate (§10.5, D-185) ----------------------------------------------------
# The boundary sweep walks the full list with *discrete* entries gated by
# `(idx − Φ) % D == 0`. The gate is a wrapper only discrete entries wear, so a
# continuous entry pays nothing at a boundary, and the interior walk — compiled
# from continuous entries alone — never meets an index at all: an empty due set
# is arity selection, never a sentinel index failing every gate (D-185).

struct Gated{E}
    e::E
    D::Int
    Φ::Int
end

@inline run_at!(e, store, xbuf, ẋbuf, tick::Int) = run!(e, store, xbuf, ẋbuf)

@inline function run_at!(g::Gated, store, xbuf, ẋbuf, tick::Int)
    # Under the canonical residue 0 ≤ Φ < D, truncated rem is never 0 on the
    # negative pre-first-tick differences, so one subtraction and one remainder
    # are the whole admission test — and boundary zero's "everything with Φ = 0"
    # is this same gate at index 0, implemented by nothing (§10.5).
    (tick - g.Φ) % g.D == 0 && run!(g.e, store, xbuf, ẋbuf)
    nothing
end

# --- the walk -----------------------------------------------------------------

struct Chunk{E<:Tuple,S,X,CL}
    entries::E
    store::S
    xbuf::X
    ẋbuf::X
    clock::CL
end

@noinline (c::Chunk)() = _walk(c.entries, c.store, c.xbuf, c.ẋbuf)
@noinline (c::Chunk)(tick::Int) = _walk_at(c.entries, c.store, c.xbuf, c.ẋbuf, tick)

@inline _walk(::Tuple{}, store, xbuf, ẋbuf) = nothing
@inline function _walk(t::Tuple, store, xbuf, ẋbuf)
    run!(t[1], store, xbuf, ẋbuf)
    _walk(Base.tail(t), store, xbuf, ẋbuf)
end

@inline _walk_at(::Tuple{}, store, xbuf, ẋbuf, tick::Int) = nothing
@inline function _walk_at(t::Tuple, store, xbuf, ẋbuf, tick::Int)
    run_at!(t[1], store, xbuf, ẋbuf, tick)
    _walk_at(Base.tail(t), store, xbuf, ẋbuf, tick)
end

"""
A phase body: **two chunk tuples compiled from one entry list** (§9.7, §10.5).

The zero-arg call is the *interior* variant, walking continuous entries only —
what RK stage evaluations and guard trial evaluations run, and what
`@ballocated(body()) == 0` measures. The ZOH therefore holds mid-step **by
construction**: discrete entries are not gated out at runtime, they are absent
from the compiled walk, and the hot path carries no gating test at all.

The one-arg call is the *boundary* variant, walking the full list with each
discrete entry gated by `(tick - Φ) % D` against the index it takes (§10.5,
D-185).
"""
struct PhaseBody{I<:Tuple,B<:Tuple}
    interior::I
    boundary::B
end

@inline (b::PhaseBody)() = _walkchunks(b.interior)
@inline (b::PhaseBody)(tick::Int) = _walkchunks(b.boundary, tick)

@inline _walkchunks(::Tuple{}) = nothing
@inline function _walkchunks(t::Tuple)
    t[1]()
    _walkchunks(Base.tail(t))
end

@inline _walkchunks(::Tuple{}, tick::Int) = nothing
@inline function _walkchunks(t::Tuple, tick::Int)
    t[1](tick)
    _walkchunks(Base.tail(t), tick)
end

# Construction is type-opaque: entries are built into untyped buffers and
# splatted once per chunk; the compiled tuple's only consumer is the walk.
# `gates` runs parallel to `entries`: `nothing` for a continuous entry, the
# compiled `(D, Φ)` pair for a discrete one — which is also what selects the
# interior subset, discreteness being a build-time fact (§10.5, D-147).
function chunked_body(entries::Vector, gates::Vector, store, xbuf, ẋbuf, clock;
                      chunk_size::Int = 16)
    chunks(es) = tuple((Chunk(tuple(es[lo:min(lo + chunk_size - 1, length(es))]...),
                              store, xbuf, ẋbuf, clock)
                        for lo in 1:chunk_size:length(es))...)
    PhaseBody(chunks([e for (e, gt) in zip(entries, gates) if gt === nothing]),
              chunks([gt === nothing ? e : Gated(e, gt[1], gt[2])
                      for (e, gt) in zip(entries, gates)]))
end
