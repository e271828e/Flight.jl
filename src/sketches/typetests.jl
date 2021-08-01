

struct Type1{D <: Int}
    a::Vector{T} where {T<:Real}
    b::Dict{Symbol,D}
end