# The compiled executor (§12.7): schedule entries carrying what selects code in
# type parameters and what is plain data in fields, gathered into chunked
# statically-typed tuples behind non-inlined barriers, traversed by a
# compile-time-unrolled walk.

# --- entries ------------------------------------------------------------------
# Two kinds on the continuous tier: a stage entry writes cells, an RHS entry
# writes the flat `ẋ` buffer. Both build their §5.2 bundle the same way, from
# `BN` — the bundle name set the law fixed at build time.

struct StageEntry{F,Comp,XT,BN,IA<:NamedTuple,YA<:NamedTuple,OA<:NamedTuple,CL}
    fn::F
    comp::Comp
    inputs::IA      # input face => cell address (the wiring's name binding)
    yx::YA          # own stage-1 port => cell address
    outs::OA        # port this entry writes => cell address
    x_off::Int
    clock::CL
end

struct RHSEntry{Comp,XT,BN,IA<:NamedTuple,YA<:NamedTuple,CL}
    comp::Comp
    inputs::IA
    y::YA           # every own port — `f` reads the complete fresh table (§5.3)
    x_off::Int
    clock::CL
end

# One bundle-expression builder, two @generated entry points. Absent names are
# absent, never `nothing`-filled: a body destructuring what it does not own
# fails at the destructuring (§5.2).
function _bundle_expr(BN, XT)
    args = map(BN) do n
        n === :x   ? :(reconstruct($XT, xbuf, e.x_off)) :
        n === :u   ? :(gather_group(e.inputs, store)) :
        n === :y_x ? :(gather_group(e.yx, store)) :
        n === :y   ? :(gather_group(e.y, store)) :
        n === :t   ? :(e.clock.t) :
        error("no source for bundle field $n")
    end
    :(NamedTuple{$BN}(($(args...),)))
end

@generated function make_bundle(e::StageEntry{F,Comp,XT,BN}, store, xbuf) where {F,Comp,XT,BN}
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

@inline function run!(e::StageEntry, store, xbuf, ẋbuf)
    y = e.fn(e.comp, make_bundle(e, store, xbuf))
    scatter_group!(store, e.outs, y)
end

@inline function run!(e::RHSEntry, store, xbuf, ẋbuf)
    ẋ = f(e.comp, make_bundle(e, store, xbuf))
    flatten!(ẋbuf, e.x_off, ẋ)      # shape conformance established at probe time
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

@inline _walk(::Tuple{}, store, xbuf, ẋbuf) = nothing
@inline function _walk(t::Tuple, store, xbuf, ẋbuf)
    run!(t[1], store, xbuf, ẋbuf)
    _walk(Base.tail(t), store, xbuf, ẋbuf)
end

"""
A phase body. Both arities off one entry list (§12.7): the zero-arg call is the
*interior* variant, which is what RK stage evaluations and guard probes run and
what `@ballocated(body()) == 0` measures; the one-arg call is the *boundary*
variant, which gates discrete entries by modulo against the tick index.
Increment 2 has no discrete entries, so the two coincide — the arity exists so
increment 3 fills the boundary variant rather than restructuring this.
"""
struct PhaseBody{C<:Tuple}
    chunks::C
end

@inline (b::PhaseBody)() = _walkchunks(b.chunks)
@inline (b::PhaseBody)(tick::Int) = _walkchunks(b.chunks)

@inline _walkchunks(::Tuple{}) = nothing
@inline function _walkchunks(t::Tuple)
    t[1]()
    _walkchunks(Base.tail(t))
end

# Construction is type-opaque: entries are built into untyped buffers and
# splatted once per chunk; the compiled tuple's only consumer is the walk.
function chunked_body(entries::Vector, store, xbuf, ẋbuf, clock; chunk_size::Int = 16)
    chunks = Any[]
    for lo in 1:chunk_size:length(entries)
        hi = min(lo + chunk_size - 1, length(entries))
        push!(chunks, Chunk(tuple(entries[lo:hi]...), store, xbuf, ẋbuf, clock))
    end
    PhaseBody(tuple(chunks...))
end
