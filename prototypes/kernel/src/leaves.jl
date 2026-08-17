# Leaf walk over the closed value vocabulary (spec §7.1): real scalars, static
# arrays, and isbits structs whose fields are drawn from the same vocabulary.
#
# Shared by both candidates: the continuous state buffer is flat by decision
# (§7.1), so *some* flatten/reconstruct machinery is needed either way. C2 then
# reuses it for cells; C1 uses it only for state.

using StaticArrays

"""
    nleaves(P)

Number of activation-scalar leaves a value of type `P` occupies in a flat
buffer. Layout-time only — never on an evaluation path.
"""
nleaves(::Type{<:Real}) = 1
nleaves(::Type{P}) where {P<:StaticArray} = length(P) * nleaves(eltype(P))
nleaves(::Type{P}) where {P} = sum(nleaves, fieldtypes(P); init = 0)

"""
    leaf_types(P)

The element type of each leaf a value of type `P` occupies, in flat order.
Layout-time only — the shape counterpart of `nleaves`, used to check a declared
cell against the store it has to live in.
"""
leaf_types(::Type{P}) where {P<:Real} = Type[P]
leaf_types(::Type{P}) where {P<:StaticArray} = repeat(leaf_types(eltype(P)), length(P))
leaf_types(::Type{P}) where {P} = reduce(vcat, (leaf_types(FT) for FT in fieldtypes(P)); init = Type[])

"""
    leaf_eltype(P)

The element type shared by every leaf of `P`. The compile-time selector the
store bundle dispatches on; meaningful only for leaf-homogeneous `P`, which
layout has already checked via `leaf_types`.
"""
leaf_eltype(::Type{P}) where {P<:Real} = P
leaf_eltype(::Type{P}) where {P<:StaticArray} = leaf_eltype(eltype(P))
leaf_eltype(::Type{P}) where {P} = leaf_eltype(fieldtype(P, 1))

# --- expression builders (compile time) --------------------------------------

# Returns (expr, next_base): `expr` reconstructs a `P` from `buf` starting at
# `off + base + 1`, with all indices static relative to `off`.
function _reconstruct_expr(::Type{P}, base::Int) where {P}
    if P <: Real
        return :(@inbounds buf[off + $(base + 1)]), base + 1
    elseif P <: StaticArray
        args = Expr[]
        b = base
        for _ in 1:length(P)
            e, b = _reconstruct_expr(eltype(P), b)
            push!(args, e)
        end
        return Expr(:call, P, args...), b
    else
        args = Expr[]
        b = base
        for FT in fieldtypes(P)
            e, b = _reconstruct_expr(FT, b)
            push!(args, e)
        end
        # NamedTuples take their fields as one tuple; everything else positionally
        P <: NamedTuple && return Expr(:call, P, Expr(:tuple, args...)), b
        return Expr(:call, P, args...), b
    end
end

# Returns (block, next_base): statements storing the leaves of the value
# denoted by expression `v` into `buf` starting at `off + base + 1`.
function _flatten_expr(::Type{P}, v, base::Int) where {P}
    stmts = Expr[]
    if P <: Real
        push!(stmts, :(@inbounds buf[off + $(base + 1)] = $v))
        return Expr(:block, stmts...), base + 1
    elseif P <: StaticArray
        b = base
        for i in 1:length(P)
            blk, b = _flatten_expr(eltype(P), :(@inbounds $v[$i]), b)
            push!(stmts, blk)
        end
        return Expr(:block, stmts...), b
    else
        b = base
        for (i, FT) in enumerate(fieldtypes(P))
            blk, b = _flatten_expr(FT, :(getfield($v, $i)), b)
            push!(stmts, blk)
        end
        return Expr(:block, stmts...), b
    end
end

# --- evaluation-path entry points --------------------------------------------

"""
    reconstruct(P, buf, off)

Materialize a `P` from `nleaves(P)` consecutive entries of `buf` starting at
`off + 1`. Fully unrolled; register-level for isbits `P`.
"""
@generated function reconstruct(::Type{P}, buf::AbstractVector, off::Int) where {P}
    expr, _ = _reconstruct_expr(P, 0)
    quote
        $(Expr(:meta, :inline))
        $expr
    end
end

"""Zero-valued `P` at scalar `T`. Build time only (cell initialization)."""
zero_value(::Type{P}, ::Type{T}) where {P,T} = reconstruct(P, zeros(T, nleaves(P)), 0)

# --- the activation walk (§7.2) ----------------------------------------------
# Declarations are written at concrete `Float64`; cell and view types per
# activation come from this walk over them, never from inference through user
# code. `Float64` leaves and `Float64` type parameters follow the activation
# scalar; everything else (`Int`, `Bool`, sizes) is pinned and passes through.

"""
    retype(T, P)

`P` with every `Float64` position replaced by `T`. Build time only.
"""
retype(::Type{T}, ::Type{Float64}) where {T} = T
function retype(::Type{T}, ::Type{P}) where {T,P}
    isconcretetype(P) && isempty(P.parameters) && return P
    isempty(P.parameters) && return P
    P.name.wrapper{(p isa Type ? retype(T, p) : p for p in P.parameters)...}
end

"""Value counterpart: the same value with its leaves converted to `T`."""
function retype_value(::Type{T}, v) where {T}
    P = retype(T, typeof(v))
    reconstruct(P, T[T(l) for l in _leaf_values(v)], 0)
end

_leaf_values(v::Real) = (v,)
_leaf_values(v::StaticArray) = Iterators.flatten(map(_leaf_values, Tuple(v)))
_leaf_values(v::NamedTuple) = Iterators.flatten(map(_leaf_values, values(v)))
_leaf_values(v::Tuple) = Iterators.flatten(map(_leaf_values, v))
_leaf_values(v) = Iterators.flatten(map(_leaf_values,
    ntuple(i -> getfield(v, i), fieldcount(typeof(v)))))

"""
    flatten!(buf, off, v)

Store the leaves of `v` into `buf` starting at `off + 1`. Returns `nothing`.
"""
@generated function flatten!(buf::AbstractVector, off::Int, v::P) where {P}
    block, _ = _flatten_expr(P, :v, 0)
    quote
        $(Expr(:meta, :inline))
        $block
        nothing
    end
end
