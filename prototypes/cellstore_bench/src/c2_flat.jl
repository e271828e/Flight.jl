# C2 — per-eltype homogeneous stores, field offsets.
#
# Cells are flattened by the §7.1 leaf walk into one contiguous buffer per
# element type — the same construction the state buffer already uses, pointed at
# the signal table. An address carries the port *type* as a parameter (so the
# gather stays statically dispatched and the reconstruction unrolls) and the
# offset as an ordinary `Int` *field*.
#
# The structural consequence, and the whole claim under test: two instances of
# the same component type have the same address-token types and differ only in
# field values, so they share one `Entry` type and one compiled body. The store
# type is `C2Store{T}` whatever the model.
#
# The risk this measures rather than argues: an offset in an immutable field of
# an object the compiler can see as constant is const-propagable, which would
# respecialize the shared body back into one body per instance. `invoke_sweep`'s
# `@constprop :none` barrier keeps the measurement on the real-loop side of
# that, where the simulation is a runtime value.
#
# Only one buffer appears here because the interior continuous sweep is entirely
# at the activation scalar. Discrete-side stores (`Int`, `Bool`, mode tags) are
# the same construction at another eltype and do not enter this measurement.

struct C2 <: Candidate end

struct FlatAddr{P}
    off::Int
end

struct C2Store{T}
    buf::Vector{T}
end

@inline gather(s::C2Store, a::FlatAddr{P}) where {P} = reconstruct(P, s.buf, a.off)

@inline scatter!(s::C2Store, a::FlatAddr{P}, v) where {P} = flatten!(s.buf, a.off, v)

function build_store(::Type{C2}, spec::ModelSpec, ::Type{T}) where {T}
    layout = Dict{Tuple{Int,Symbol},Tuple{Int,Type}}()
    off = 0
    for (ci, port, P) in cell_inventory(spec, T)
        layout[(ci, port)] = (off, P)
        off += nleaves(P)
    end
    C2Store(zeros(T, off)), layout
end

function cell_addr(::Type{C2}, layout, ci::Int, port::Symbol)
    off, P = layout[(ci, port)]
    FlatAddr{P}(off)
end
