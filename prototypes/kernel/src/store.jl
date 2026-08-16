# The signal table: per-eltype homogeneous cell stores with build-time offsets
# (D-162, §9.7). Cells are flattened by the §7.1 leaf walk into one
# contiguous buffer per element type; an address carries the port type as its
# parameter and the offset as an ordinary field, so instances of one component
# type share one compiled body.
#
# Increment 2 is continuous-only, so one buffer at the activation scalar covers
# the whole table. The discrete tier's `Int`/`Bool` cells are the same
# construction at another eltype.

struct CellAddr{P}
    off::Int
end

struct CellStore{T}
    buf::Vector{T}
end

@inline gather(s::CellStore, a::CellAddr{P}) where {P} = reconstruct(P, s.buf, a.off)
@inline scatter!(s::CellStore, a::CellAddr{P}, v) where {P} = flatten!(s.buf, a.off, v)

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
