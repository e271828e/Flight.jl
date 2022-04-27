module Misc

using DataStructures: OrderedDict

export Bounded, linear_scaling
export pwf

#convenience
Base.NamedTuple(dict::OrderedDict{Symbol, T} where {T}) = NamedTuple{Tuple(keys(dict))}(values(dict))

#print with fieldnames
function pwf(s)
    for f in fieldnames(typeof(s))
        println("$f: $(getfield(s,f))")
    end
end

function linear_scaling(u::Bounded{T, UMin, UMax}, range::NTuple{2,Real}) where {T, UMin, UMax}
    @assert UMin != UMax
    return range[1] + (range[2] - range[1])/(UMax - UMin) * (T(u) - UMin)
end


################################### Bounded ####################################
#needs some unit tests
struct Bounded{T<:Real, Min, Max}
    val::T
    function Bounded(val::T, min_val::Real, max_val::Real) where {T}
        new{T, min_val, max_val}(min(max(val, min_val), max_val))
    end
end

Bounded{T}(x::Bounded) where {T} = convert(Bounded{T}, x)
Bounded{T,Min,Max}(x::Bounded) where {T,Min,Max} = convert(Bounded{T,Min,Max}, x)
(T::Type{<:Real})(x::Bounded) = convert(T, x)

Base.convert(::Type{T}, x::Bounded) where {T<:Real} = T(x.val)
Base.convert(::Type{Bounded{T,Min,Max}}, x::Real) where {T, Min, Max} = Bounded(T(x), Min, Max)

function Base.convert(::Type{Bounded{T1,Min,Max}}, x::Bounded{T2}) where {T1, T2, Min, Max}
    Bounded(T1(x.val), Min, Max)
end

function Base.convert(::Type{Bounded{T1}}, x::Bounded{T2,Min,Max}) where {T1, T2, Min, Max}
    #since the conversion target does not specify bounds, take them from the conversion source
    Bounded(T1(x.val), Min, Max)
end

function Base.promote_rule(::Type{Bounded{T1,Min,Max}}, ::Type{T2}) where {T1, T2, Min, Max}
    Bounded{promote_type(T1,T2),Min,Max}
end

#basic addition and subtraction
Base.:+(x::Bounded{T1,Min,Max}, y::Real) where {T1, Min, Max} = Bounded(x.val + y, Min, Max)
Base.:-(x::Bounded{T1,Min,Max}, y::Real) where {T1, Min, Max} = Bounded(x.val - y, Min, Max)
Base.:-(x::Bounded{T1,Min,Max}) where {T1, Min, Max} = Bounded(-x.val, Min, Max)

#both Bounded must have identical underlying types and bounds, otherwise there
#is no easy way of deciding whose bounds should win
Base.:+(x::B, y::B) where {B <: Bounded{T,Min,Max}} where {T, Min, Max} = Bounded(x.val + y.val, Min, Max)
Base.:-(x::B, y::B) where {B <: Bounded{T,Min,Max}} where {T, Min, Max} = Bounded(x.val - y.val, Min, Max)

#basic equality
Base.:(==)(x::Bounded{T1}, y::Real) where {T1, T2} = (==)(promote(x.val, y)...)


end#module