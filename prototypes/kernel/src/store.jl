# The signal table: per-eltype homogeneous cell stores with build-time offsets
# (D-162, §9.7). Cells are flattened by the §7.1 leaf walk into one contiguous
# buffer per element type; a mixed-leaf cell (§7.2's per-leaf table, D-166
# pinning) simply spans several of them. An address carries the port type as
# its parameter and one cursor per distinct leaf eltype as an `NTuple` *field*,
# so instances of one component type share one compiled body; `K` is a pure
# function of `P` (bench-confirmed 2026-08-20: the C2M point in
# `../cellstore_bench` keeps D-162's flat curve), and the homogeneous cell is
# the `K = 1` case of the same representation.

struct CellAddr{P,K}
    offs::NTuple{K,Int}
end

struct CellStore{T}
    buf::Vector{T}
end

"""
The store bundle: one homogeneous `CellStore` per element type present in the
model's cells, keyed by the eltype's name — the *plural* in §9.7's "per-eltype
stores". One concrete bundle type per model, so `Chunk`'s store parameter stays
a single type: chunk-type count, not model size, is what bounds the compile
curve D-162 measured. Selection is static — the address's port type names its
leaf eltypes at compile time, so a gather binds exactly the buffers `P`'s
leaves live in, with no runtime lookup — and each of an address's offsets is
relative to its own eltype's buffer, `leaf_eltypes` fixing the order.
"""
struct StoreBundle{NT<:NamedTuple}
    stores::NT
end

# The bundle's field name for one leaf eltype. `Symbol(L)` will not do: `show`
# for a type abbreviates the module prefix when the type is visible from the
# *printing* module, and the module a `@generated` body prints from is not the
# one `compile` prints from — so a scalar tagged by a type outside `Base` (a
# service's own seeding tag, §9.4) keys the bundle one way at construction and
# is looked up another way at the gather. Printing with no module context at all
# fixes one fully-qualified spelling for every site.
_cell_key(::Type{L}) where {L} = Symbol(sprint(show, L; context = :module => nothing))

@generated function gather(b::StoreBundle, a::CellAddr{P,K}) where {P,K}
    Ls = leaf_eltypes(P)
    binds = [:($(Symbol(:buf, k)) = getfield(b.stores, $(QuoteNode(_cell_key(L)))).buf)
             for (k, L) in enumerate(Ls)]
    expr = _mreconstruct_expr(P, Ls, zeros(Int, K))
    quote
        $(Expr(:meta, :inline))
        offs = a.offs
        $(binds...)
        $expr
    end
end

@generated function scatter!(b::StoreBundle, a::CellAddr{P,K}, v) where {P,K}
    Ls = leaf_eltypes(P)
    binds = [:($(Symbol(:buf, k)) = getfield(b.stores, $(QuoteNode(_cell_key(L)))).buf)
             for (k, L) in enumerate(Ls)]
    block = _mflatten_expr(P, :v, Ls, zeros(Int, K))
    quote
        $(Expr(:meta, :inline))
        offs = a.offs
        $(binds...)
        $block
        nothing
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
a step — RK stages evaluate at internal times. `step` counts completed
continuous steps since `t₀`; every `n`-th step boundary is a base tick (§10.5),
and no entry reads it — it is the loop's, not the bundle's. `t₀` anchors the
indexed grid: frame tops are `t₀ + k·h`, computed from the index and never
accumulated, so a remainder step's float arithmetic cannot drift the grid
(§10.4).
"""
mutable struct Clock{T}
    t::T
    step::Int
    t₀::T
end
Clock(t) = Clock(t, 0, t)
