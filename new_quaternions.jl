"""
Fast, lightweight quaternion module.
"""
# module NewQuaternions

using StaticArrays: SVector, MVector
using LinearAlgebra

export Quaternion, Quat, UnitQuat

######################## Quat #############################

abstract type NormConstraint end
struct FreeNorm <: NormConstraint end
struct UnitNorm <: NormConstraint end

struct Quaternion{N <: NormConstraint} <: AbstractVector{Float64}
    _sv::SVector{4, Float64}
    function Quaternion{N}(v::AbstractVector, normalization::Bool) where {N}
        sv = SVector{4,Float64}(v) #faster normalization method
        return normalization ? new{N}(normalize(sv)) : new{N}(sv)
    end
end

const FreeQuat = Quaternion{FreeNorm}
const UnitQuat = Quaternion{UnitNorm}

Quaternion{FreeNorm}(v::AbstractVector; normalization::Bool = false) = FreeQuat(v, normalization)
Quaternion{UnitNorm}(v::AbstractVector; normalization::Bool = true) = UnitQuat(v, normalization)

Quaternion{N}(s::Real) where {N} = convert(Quaternion{N}, s)

function Quaternion{FreeNorm}(; real = 0.0, imag = zeros(SVector{3}))
    FreeQuat(vcat(real, SVector{3,Float64}(imag)))
end

function Quaternion{UnitNorm}(; real = 0.0, imag = zeros(SVector{3}), normalization::Bool = true)
    sv = vcat(real, SVector{3,Float64}(imag))
    any(sv .!= 0) || throw(ArgumentError("UnitQuat requires a non-zero real or imaginary part"))
    UnitQuat(vcat(real, SVector{3,Float64}(imag)); normalization)
end

Base.convert(::Type{Quaternion{N}}, s::Real) where {N} = Quaternion{N}(SVector{4,Float64}(s, 0, 0, 0))
Base.convert(::Type{Quaternion{N}}, v::AbstractVector) where {N} = Quaternion{N}(v)

Base.promote_rule(::Type{<:Quaternion}, ::Type{S}) where {S<:Real} = FreeQuat
Base.promote_rule(::Type{Quaternion{N}}, ::Type{<:AbstractVector{<:Real}}) where {N} = Quaternion{N}
Base.promote_rule(::Type{UnitQuat}, ::Type{FreeQuat}) = FreeQuat

#### AbstractVector interface ####
Base.size(::Quaternion) = (4,)
Base.length(::Quaternion) = 4
Base.firstindex(::Quaternion) = 1
Base.lastindex(::Quaternion) = 4
Base.eltype(::Quaternion) = Float64

#show the specific type when printing
Base.show(io::IO, q::Quaternion) = print(io, "$(typeof(q))($(q[:]))")
Base.show(io::IO, ::MIME"text/plain", q::Quaternion) = print(io, "$(typeof(q))($(q[:]))")

#for retrieving real and imaginary parts
Base.copy(q::Quaternion{N}) where {N} = Quaternion{N}(copy(getfield(q, :_sv)))
Base.getindex(q::Quaternion, i) = getfield(q, :_sv)[i]
Base.getindex(q::Quaternion, s::Symbol) = getindex(q, Val(s))
Base.getindex(q::Quaternion, ::Val{:real}) = getfield(q, :_sv)[1]
Base.getindex(q::Quaternion, ::Val{:imag}) = SVector{3, Float64}(@view getfield(q, :_sv)[2:4])
Base.getproperty(q::Quaternion, s::Symbol) = getindex(q, Val(s))

LinearAlgebra.norm(q::Quaternion) = norm(getfield(q, :_sv)) #uses StaticArrays implementation
LinearAlgebra.normalize(f::FreeQuat) = FreeQuat(normalize(getfield(f, :_sv)))
LinearAlgebra.normalize(u::UnitQuat) = UnitQuat(getfield(u, :_sv))
@inline norm_sqr(q::Quaternion) = (data = getfield(q,:_sv); sum(data.*data))

#### Adjoint & Inverse
@inline Base.adjoint(q::Quaternion) = conj(q)

@inline Base.conj(q::FreeQuat) = FreeQuat(vcat(q.real, -q.imag))
@inline Base.inv(q::FreeQuat) = FreeQuat(getfield(q', :_sv) / norm_sqr(q))

@inline Base.conj(u::UnitQuat) = UnitQuat(vcat(u.real, -u.imag), false)
@inline Base.inv(u::UnitQuat) = u'

#### Operators
Base.:+(q::FreeQuat) = q
Base.:-(q::FreeQuat) = FreeQuat(-getfield(q, :_sv))
Base.:+(u::UnitQuat) = u
Base.:-(u::UnitQuat) = UnitQuat(-getfield(u, :_q), false)

Base.:(==)(q1::Quaternion, q2::Quaternion) = getfield(q1,:_sv) == getfield(q2,:_sv)
Base.:(≈)(q1::Quaternion, q2::Quaternion) = getfield(q1,:_sv) ≈ getfield(q2,:_sv)

Base.:+(q1::FreeQuat, q2::FreeQuat) = FreeQuat(getfield(q1,:_sv) + getfield(q2,:_sv))
Base.:+(q1::Quaternion, q2::Quaternion) = +(promote(q1,q2)...)
Base.:+(q::Quaternion, a::Real) = +(promote(q, a)...)
Base.:+(a::Real, q::Quaternion) = +(promote(a, q)...)

Base.:-(q1::FreeQuat, q2::FreeQuat) = FreeQuat(getfield(q1,:_sv) + getfield(q2,:_sv))
Base.:-(q1::Quaternion, q2::Quaternion) = -(promote(q1,q2)...)
Base.:-(q::Quaternion, a::Real) = -(promote(q, a)...)
Base.:-(a::Real, q::Quaternion) = -(promote(a, q)...)

Base.:*(q1::Quaternion{N1}, q2::Quaternion{N2}) where {N1,N2} = *(promote(q1,q2)...)

@inline function Base.:*(q1::Quaternion{N}, q2::Quaternion{N}) where {N}

    p_re = q1.real * q2.real - q1.imag ⋅ q2.imag
    p_im = q1.real * q2.imag + q2.real * q1.imag + q1.imag × q2.imag

    #bypass normalization for UnitQuat!
    Quaternion{N}(vcat(p_re, p_im), false)
end

#multiplying by a scalar must not yield a UnitQuat
Base.:*(q::Quaternion, a::Real) = a * q
Base.:*(a::Real, q::Quaternion) = FreeQuat(a * getfield(q, :_sv))

Base.:/(q1::Quaternion{N}, q2::Quaternion{N}) where {N} = Quaternion{N}(q1 * inv(q2))
Base.:/(q1::Quaternion{N1}, q2::Quaternion{N2}) where {N1,N2} = /(promote(q1,q2)...)
Base.:/(q::Quaternion, a::Real) = FreeQuat(getfield(q, :_sv)/ a)
Base.:/(a::Real, q::Quaternion) = /(promote(a, q)...)

Base.:\(q1::Quaternion{N}, q2::Quaternion{N}) where {N} = Quaternion{N}(inv(q1) * q2) #!= /(q2, q1) == q2 * inv(q1)
Base.:\(q1::Quaternion{N1}, q2::Quaternion{N2}) where {N1,N2} = \(promote(q1,q2)...)
Base.:\(q::Quaternion, a::Real) = \(promote(q, a)...)
Base.:\(a::Real, q::Quaternion) = q / a

#end