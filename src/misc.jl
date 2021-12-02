module Misc

export Bounded

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

#addition and subtraction between Bounded and Reals not required by now

end#module