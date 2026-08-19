# The compiled executor (§9.7): schedule entries carrying what selects code in
# type parameters and what is plain data in fields, gathered into chunked
# statically-typed tuples behind non-inlined barriers, traversed by a
# compile-time-unrolled walk.

# --- entries ------------------------------------------------------------------
# Three kinds, by where the product goes: a stage entry writes cells (both
# tiers — the stage names are shared, D-173), an RHS entry writes the flat `ẋ`
# buffer, an update entry writes its own discrete state store. All build their
# §5.2 bundle the same way, from `BN` — the bundle name set the law fixed at
# build time.
#
# What selects code sits in type parameters, what varies per instance in fields
# (§9.7): the state store is a `Ref` whose *type* is shared by every instance of
# a component type, so two instances still compile to one body.

struct StageEntry{F,Comp,XT,BN,IA<:NamedTuple,YA<:NamedTuple,OA<:NamedTuple,CL,XS,MS,WS}
    fn::F
    comp::Comp
    inputs::IA      # input face => cell address (the wiring's name binding)
    yx::YA          # own stage-1 port => cell address
    outs::OA        # port this entry writes => cell address
    x_off::Int      # continuous state offset into the flat buffer
    clock::CL
    xstore::XS      # discrete state store, or nothing on the continuous tier
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

struct UpdateEntry{Comp,BN,IA<:NamedTuple,YA<:NamedTuple,CL,XS,WS}
    comp::Comp
    inputs::IA
    y::YA           # every own port — `g` reads the complete fresh table too
    clock::CL
    xstore::XS      # written by this entry, and by nothing else
    ws::WS
    Δt::Float64
end

# Outer constructors: only `XT`/`BN` cannot be inferred from the arguments.
StageEntry{XT,BN}(fn, comp, inputs, yx, outs, x_off, clock, xstore, mstore, ws, Δt) where {XT,BN} =
    StageEntry{typeof(fn),typeof(comp),XT,BN,typeof(inputs),typeof(yx),typeof(outs),
               typeof(clock),typeof(xstore),typeof(mstore),typeof(ws)}(
        fn, comp, inputs, yx, outs, x_off, clock, xstore, mstore, ws, Δt)

RHSEntry{XT,BN}(comp, inputs, y, x_off, clock, mstore, ws) where {XT,BN} =
    RHSEntry{typeof(comp),XT,BN,typeof(inputs),typeof(y),typeof(clock),
             typeof(mstore),typeof(ws)}(comp, inputs, y, x_off, clock, mstore, ws)

UpdateEntry{BN}(comp, inputs, y, clock, xstore, ws, Δt) where {BN} =
    UpdateEntry{typeof(comp),BN,typeof(inputs),typeof(y),typeof(clock),
                typeof(xstore),typeof(ws)}(comp, inputs, y, clock, xstore, ws, Δt)

# One bundle-expression builder, three @generated entry points. Absent names are
# absent, never `nothing`-filled: a body destructuring what it does not own
# fails at the destructuring (§5.2). `x` reads from the flat buffer on the
# continuous tier and from the component's own store on the discrete one — the
# shared state letter (D-173), two homes.
function _bundle_expr(BN, XT, XS)
    args = map(BN) do n
        n === :x   ? (XS === Nothing ? :(reconstruct($XT, xbuf, e.x_off)) : :(e.xstore[])) :
        n === :m   ? :(e.mstore[]) :
        n === :u   ? :(gather_group(e.inputs, store)) :
        n === :y_x ? :(gather_group(e.yx, store)) :
        n === :y   ? :(gather_group(e.y, store)) :
        n === :ws  ? :(e.ws) :
        n === :t   ? :(e.clock.t) :
        n === :Δt  ? :(e.Δt) :
        error("no source for bundle field $n")
    end
    :(NamedTuple{$BN}(($(args...),)))
end

@generated function make_bundle(e::StageEntry{F,Comp,XT,BN,IA,YA,OA,CL,XS}, store,
                                xbuf) where {F,Comp,XT,BN,IA,YA,OA,CL,XS}
    quote
        $(Expr(:meta, :inline))
        $(_bundle_expr(BN, XT, XS))
    end
end

@generated function make_bundle(e::RHSEntry{Comp,XT,BN}, store, xbuf) where {Comp,XT,BN}
    quote
        $(Expr(:meta, :inline))
        $(_bundle_expr(BN, XT, Nothing))
    end
end

@generated function make_bundle(e::UpdateEntry{Comp,BN}, store, xbuf) where {Comp,BN}
    quote
        $(Expr(:meta, :inline))
        $(_bundle_expr(BN, Nothing, Ref))
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
    e.xstore[] = g(e.comp, make_bundle(e, store, xbuf))
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
