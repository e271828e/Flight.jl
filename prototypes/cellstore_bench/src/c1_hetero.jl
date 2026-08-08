# C1 — heterogeneous store, type-domain indices.
#
# One store enumerating every cell in its own type: a tuple of concretely-typed
# mutable cells. An address is a *singleton whose type carries the cell index*,
# so a gather is `getfield` at a compile-time constant followed by a load. No
# offsets exist at runtime; addressing is entirely in the type domain.
#
# The structural consequence, which is what gate 2 prices: two instances of the
# same component type have different cell indices, hence different address-token
# types, hence different `Entry` types — so the executor compiles a separate
# body per instance, and the store type itself grows with the model and appears
# in every one of those signatures.
#
# Variant not built: a per-model generated mutable struct with inline typed
# fields (better locality, one nominal type per model instead of a long tuple
# type). It differs from this form only on the gate-3 tiebreakers — the
# type-domain addressing that gate 2 measures is identical — so it is worth
# building only if C1 survives gate 2 and loses on runtime alone.

struct C1 <: Candidate end

struct CellIdx{I} end

struct C1Store{C<:Tuple}
    cells::C
end

@inline gather(s::C1Store, ::CellIdx{I}) where {I} = getfield(s.cells, I)[]

@inline function scatter!(s::C1Store, ::CellIdx{I}, v) where {I}
    getfield(s.cells, I)[] = v
    nothing
end

function build_store(::Type{C1}, spec::ModelSpec, ::Type{T}) where {T}
    inv = cell_inventory(spec, T)
    cells = Any[Ref(zero_value(P, T)) for (_, _, P) in inv]
    index = Dict{Tuple{Int,Symbol},Int}(
        (ci, port) => i for (i, (ci, port, _)) in enumerate(inv))
    C1Store(tuple(cells...)), index
end

cell_addr(::Type{C1}, layout, ci::Int, port::Symbol) = CellIdx{layout[(ci, port)]}()
