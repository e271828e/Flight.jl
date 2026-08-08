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
