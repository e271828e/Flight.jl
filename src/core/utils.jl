module Utils

using StructArrays

export Ranged, linear_scaling
export TimeHistory
export pwf

#print with propertynames
function pwf(s)
    for f in propertynames(s)
        println("$f: $(getproperty(s,f))")
    end
end

################################################################################
################################ Ranged ########################################

#needs some unit tests
struct Ranged{T<:Real, Min, Max}
    val::T
    function Ranged(val::T, min_val::Real, max_val::Real) where {T}
        new{T, min_val, max_val}(min(max(val, min_val), max_val))
    end
end

Ranged{T}(x::Ranged) where {T} = convert(Ranged{T}, x)
Ranged{T,Min,Max}(x::Ranged) where {T,Min,Max} = convert(Ranged{T,Min,Max}, x)
(T::Type{<:Real})(x::Ranged) = convert(T, x)

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
Base.:(==)(x::Ranged{T1}, y::Real) where {T1, T2} = (==)(promote(x.val, y)...)

function linear_scaling(u::Ranged{T, UMin, UMax}, range::NTuple{2,Real}) where {T, UMin, UMax}
    @assert UMin != UMax
    return range[1] + (range[2] - range[1])/(UMax - UMin) * (T(u) - UMin)
end

function test()

    a = Ranged(1, 0, 2)
    b = Ranged(2.0, 0, 2)
    A = fill(a, 100)
    B = fill(b, 100)
    C = copy(B)

    #no allocations
    C .= A .+ B

end

################################################################################
############################ TimeHistory #######################################

mutable struct TimeHistory{V, T <: AbstractVector{Float64}, D <: AbstractVector{V}}
    _t::T
    _data::D
    function TimeHistory(t::T, data::D) where {T, D <: AbstractVector{V}} where {V}
        @assert length(t) == length(data)
        new{V, T, D}(t, data)
    end
end

TimeHistory(t::Real, data) = TimeHistory([Float64(t)], [data])

function TimeHistory(t::AbstractVector, M::Matrix)
    #each Matrix column interpreted as one Vector value
    TimeHistory(t, [M[:, i] for i in 1:size(M,2)])
end

Base.length(th::TimeHistory) = length(th._t)

function Base.getproperty(th::TimeHistory, s::Symbol)
    t = getfield(th, :_t)
    y = getfield(th, :_data)
    if s === :_t
        return t
    elseif s === :_data
        return y
    else
        return TimeHistory(t, getproperty(StructArray(y), s))
    end
end

Base.getindex(th::TimeHistory, i) = TimeHistory(th._t[i], th._data[i])
Base.view(th::TimeHistory, i) = TimeHistory(view(th._t, i), view(th._data, i))

#for inspection
get_child_names(::T) where {T <: TimeHistory} = get_child_names(T)
get_child_names(::Type{<:TimeHistory{V}}) where {V} = fieldnames(V)

#could be rewritten as @generated to avoid allocation if needed
function get_scalar_components(th::TimeHistory{<:AbstractVector{T}}) where {T<:Real}
    [TimeHistory(th._t, y) for y in th._data |> StructArray |> StructArrays.components]
end

end #module