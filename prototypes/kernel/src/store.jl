# The signal table: per-eltype homogeneous cell stores with build-time offsets
# (D-162, §9.7). Cells are flattened by the §7.1 leaf walk into one
# contiguous buffer per element type; an address carries the port type as its
# parameter and the offset as an ordinary field, so instances of one component
# type share one compiled body.

struct CellAddr{P}
    off::Int
end

struct CellStore{T}
    buf::Vector{T}
end

@inline gather(s::CellStore, a::CellAddr{P}) where {P} = reconstruct(P, s.buf, a.off)
@inline scatter!(s::CellStore, a::CellAddr{P}, v) where {P} = flatten!(s.buf, a.off, v)

"""
The store bundle: one homogeneous `CellStore` per element type present in the
model's cells, keyed by the eltype's name — the *plural* in §9.7's "per-eltype
stores". One concrete bundle type per model, so `Chunk`'s store parameter stays
a single type: chunk-type count, not model size, is what bounds the compile
curve D-162 measured. Selection is static — the address's port type names its
leaf eltype at compile time, so a gather touches exactly one buffer with no
runtime lookup — and an address's offset is relative to its own eltype's buffer.
"""
struct StoreBundle{NT<:NamedTuple}
    stores::NT
end

@generated function gather(b::StoreBundle, a::CellAddr{P}) where {P}
    key = QuoteNode(Symbol(leaf_eltype(P)))
    quote
        $(Expr(:meta, :inline))
        gather(getfield(b.stores, $key), a)
    end
end

@generated function scatter!(b::StoreBundle, a::CellAddr{P}, v) where {P}
    key = QuoteNode(Symbol(leaf_eltype(P)))
    quote
        $(Expr(:meta, :inline))
        scatter!(getfield(b.stores, $key), a, v)
    end
end

# Address groups are NamedTuples — `face => address` — so the name binding the
# wiring established is carried in the type and the gather reassembles the
# bundle the author destructures.

@generated function gather_group(addrs::NamedTuple{Ns}, store) where {Ns}
    args = [:(gather(store, addrs[$i])) for i in 1:length(Ns)]
    quote
        $(Expr(:meta, :inline))
        NamedTuple{$Ns}(($(args...),))
    end
end

@generated function scatter_group!(store, addrs::NamedTuple{Ns}, y::NamedTuple) where {Ns}
    stmts = [:(scatter!(store, addrs[$i], getfield(y, $(QuoteNode(Ns[i]))))) for i in 1:length(Ns)]
    quote
        $(Expr(:meta, :inline))
        $(stmts...)
        nothing
    end
end

"""
The clock, in its own mutable cell so the zero-arg phase bodies can close over
it. `t` is a bundle field for every stage (§5.2's bundle law) and varies within
a step — RK stages evaluate at internal times.
"""
mutable struct Clock{T}
    t::T
end
