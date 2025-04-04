module Types

using StaticArrays, StructArrays, StructTypes

using ..GUI

export shf
export Ranged, saturation, linear_scaling

################################################################################
################################################################################

function shf(s::S) where {S}
    for f in fieldnames(S)
        println("$f: $(getfield(s,f))")
    end
end

################################################################################
################################ Ranged ########################################

#needs unit tests
struct Ranged{T<:Real, Min, Max}
    val::T
    function Ranged(val::T, min_val::T, max_val::T) where {T <: Real}
        new{T, min_val, max_val}(min(max(val, min_val), max_val))
    end
end

Ranged(val::T, vmin::Real, vmax::Real) where {T} = Ranged(val, T(vmin), T(vmax))
Ranged{T}(x::Ranged) where {T} = convert(Ranged{T}, x)
Ranged{T,Min,Max}(x::Ranged) where {T,Min,Max} = convert(Ranged{T,Min,Max}, x)
Ranged{T,Min,Max}(x::Real) where {T,Min,Max} = Ranged(x, Min, Max)
(T::Type{<:Real})(x::Ranged) = convert(T, x)

Base.typemin(::Type{Ranged{T,Min,Max}}) where {T, Min, Max} = Min
Base.typemin(::T) where {T <: Ranged} = typemin(T)
Base.typemax(::Type{Ranged{T,Min,Max}}) where {T, Min, Max} = Max
Base.typemax(::T) where {T <: Ranged} = typemax(T)
Base.convert(::Type{T}, x::Ranged) where {T<:Real} = T(x.val)
Base.convert(::Type{Ranged{T,Min,Max}}, x::Real) where {T, Min, Max} = Ranged(T(x), Min, Max)

function Base.convert(::Type{Ranged{T1,Min,Max}}, x::Ranged{T2}) where {T1, T2, Min, Max}
    Ranged(T1(x.val), Min, Max)
end

#if the conversion target does not specify bounds, take them from the source
function Base.convert(::Type{Ranged{T1}}, x::Ranged{T2,Min,Max}) where {T1, T2, Min, Max}
    Ranged(T1(x.val), Min, Max)
end

function Base.promote_rule(::Type{Ranged{T1,Min,Max}}, ::Type{T2}) where {T1, T2, Min, Max}
    Ranged{promote_type(T1,T2),Min,Max}
end

#basic addition and subtraction
Base.:+(x::Ranged{T1,Min,Max}, y::Real) where {T1, Min, Max} = Ranged(x.val + y, Min, Max)
Base.:-(x::Ranged{T1,Min,Max}, y::Real) where {T1, Min, Max} = Ranged(x.val - y, Min, Max)
Base.:-(x::Ranged{T1,Min,Max}) where {T1, Min, Max} = Ranged(-x.val, Min, Max)

#bounds must be identical, since there is no easy way of deciding whose bounds
#should win
Base.:+(x::Ranged{T1,Min,Max}, y::Ranged{T2,Min,Max}) where {T1,T2,Min,Max} = Ranged(x.val + y.val, Min, Max)
Base.:-(x::Ranged{T1,Min,Max}, y::Ranged{T2,Min,Max}) where {T1,T2,Min,Max} = Ranged(x.val - y.val, Min, Max)

#basic equality
Base.:(==)(x::Ranged{T1}, y::Real) where {T1} = (==)(promote(x.val, y)...)

saturation(x::Ranged)::Int64 = (x == typemax(x)) - (x == typemin(x))

function linear_scaling(u::Ranged, range::AbstractRange)
    linear_scaling(u, (range[1], range[end]))
end

function linear_scaling(u::Ranged{T, UMin, UMax}, range::NTuple{2,Real}) where {T, UMin, UMax}
    @assert UMin != UMax
    return range[1] + (range[2] - range[1])/(UMax - UMin) * (T(u) - UMin)
end


function GUI.display_bar(label::String, source::Ranged{T,Min,Max}, args...; kwargs...) where {T<:AbstractFloat,Min,Max}
    display_bar(label, Float64(source), Min, Max, args...; kwargs...)
end

function GUI.safe_slider(label::String, source::Ranged{T,Min,Max}, args...; kwargs...) where {T<:AbstractFloat,Min,Max}
    safe_slider(label, Float64(source), Min, Max, args...; kwargs...)
end

function GUI.safe_input(label::String, source::Ranged{T,Min,Max}, args...; kwargs...) where {T<:AbstractFloat,Min,Max}
    safe_input(label, Float64(source), args...; kwargs...)
end

#enable JSON3 parsing
StructTypes.StructType(::Type{<:Ranged}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{Ranged{T,Min,Max}}) where {T,Min,Max} = T
StructTypes.lower(x::Ranged) = x.val
StructTypes.construct(::Type{Ranged{T,Min,Max}}, x::Real) where {T, Min, Max} = Ranged(x,Min,Max)

#Example 1: reading a numeric value into a Ranged field of a mutable struct
# @kwdef mutable struct MyMutableStruct
#     a::Ranged{Float64, 0.0, 1.0} = Ranged(0.5, 0.0, 0.1)
#     b::Bool = false
# end
# u = MyMutableStruct()
# JSON3.read!("""{"a": 0.1209}""", u)

#Example 2: constructing an immutable struct containing a Ranged field
# @kwdef struct MyStruct
#     a::Ranged{Float64, 0.0, 1.0} = Ranged(0.5, 0.0, 0.1)
#     b::Bool = false
# end
# StructTypes.StructType(::Type{MyStruct}) = StructTypes.Struct()
# y = JSON3.read("""{"a": 0.1209, "b": true}""", MyStruct)


# function test()
#     a = Ranged(1, 0, 2)
#     b = Ranged(2.0, 0, 2)
#     A = fill(a, 100)
#     B = fill(b, 100)
#     C = copy(B)

#     C .= A .+ B #no allocations
# end



end #module