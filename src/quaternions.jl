module Quaternions

using StaticArrays: MVector #this form of using does not allow method extension
using LinearAlgebra #this form does, with LinearAlgebra.method

export AbstractQuat, Quat, UnitQuat

######################## AbstractQuat #############################

abstract type AbstractQuat <: AbstractVector{Float64} end

#indexing and iterable interfaces; see https://docs.julialang.org/en/v1/manual/interfaces/
Base.size(::AbstractQuat) = (4,)
Base.length(::AbstractQuat) = 4
Base.firstindex(::AbstractQuat) = 1
Base.lastindex(::AbstractQuat) = 4
Base.getindex(::AbstractQuat, i) = error("AbstractQuat: getindex not implemented")
Base.setindex!(::AbstractQuat, v, i) = error("AbstractQuat: setindex! not implemented")
#not needed, inherited from AbstractVector
# Base.iterate(q::AbstractQuat, state = 1) = (state > 4 ? nothing : (q[state], state + 1))
Base.eltype(::AbstractQuat) = Float64 #helps with allocation efficiency

#display functions
Base.show(io::IO, ::MIME"text/plain", q::AbstractQuat) = print(io, "$(typeof(q)): $(q[:])")
Base.show(io::IO, q::AbstractQuat) = print(io, "$(typeof(q)): $(q[:])")

#real and imaginary parts
Base.propertynames(::Type{AbstractQuat}) = (:real, :imag)
function Base.getproperty(q::AbstractQuat, s::Symbol)
    if s == :real
        return q[1]
    elseif s == :imag
        return q[2:4]
    else
        error("No property $s defined for AbstractQuat types")
    end
end

function Base.setproperty!(q::AbstractQuat, s::Symbol, v)
    if s == :real
        q[1] = v
    elseif s == :imag
        q[2:4] = v
    else
        error("No property $s defined for Quat")
    end
end

#norm
norm_sqr(q::AbstractQuat) = sum(abs2.(q))
LinearAlgebra.norm(q::AbstractQuat) = norm(q[:]) #uses StaticArrays implementation


######################## Quat #############################

#notes:
#the implicit, automatically generated inner constructor has the form:
#Quat(input) = new(data).

#whenever typeof(input) != QData, new will call convert(QData, input). as long
#as QData provides a convert method that can handle typeof(input), then we do
#not need to handle it explicitly in an outer constructor.

# in this case, QData is a StaticVector: MVector{Float64,4}, and it can handle
#basically any AbstractVector{T} where {T<:Real} of length 4. this includes any
#AbstractQuat subtype. but not a real scalar

#Julia does NOT dispatch on keyword arguments:
# https://discourse.julialang.org/t/keyword-argument-constructor-breaks-incomplete-constructor/34198/3
# https://docs.julialang.org/en/v1/manual/methods/#Note-on-Optional-and-keyword-Arguments-1
# this means that we cannot define methods with keyword arguments of the same
# types but different names and expect Julia to choose the right one

QData = MVector{4, Float64}

struct Quat <: AbstractQuat
    data::QData
    Quat(data) = new(copy(data))
end

#outer constructors
Quat(s::Real) = Quat([s, 0, 0, 0])
Quat() = Quat(zeros(4))
Quat(; real = 0.0, imag = zeros(3)) = Quat([real, imag...])

Base.copy(q::Quat) = Quat(copy(getfield(q, :data)))
Base.getindex(q::Quat, i) = (getfield(q, :data)[i])
Base.setindex!(q::Quat, v, i) = (getfield(q, :data)[i] = v)

#more efficient than the general normalization for AbstractVector
LinearAlgebra.normalize!(q::Quat) = (normalize!(getfield(q, :data)); return q)
LinearAlgebra.normalize(q::Quat) = Quat(normalize(getfield(q, :data)))

Base.promote_rule(::Type{Quat}, ::Type{S}) where {S<:Real} = Quat
Base.convert(::Type{Quat}, a::Real) = Quat(a)
Base.convert(::Type{Quat}, a::AbstractVector) = Quat(a)

#### Adjoint & Inverse
Base.conj(q::Quat)= Quat([q.real, -q.imag...])
Base.adjoint(q::Quat) = conj(q)
Base.inv(q::Quat) = Quat(q'[:] / norm_sqr(q))

#### Operators
Base.:+(q::Quat) = q
Base.:-(q::Quat) = Quat(-q[:])

#(==) is inherited from AbstractVector, but will return true for any
#AbstractVector as long as it matches q[:], to avoid it we need to define these:
Base.:(==)(q1::AbstractVector, q2::Quat) = false
Base.:(==)(q1::Quat, q2::AbstractVector) = false
Base.:(==)(q1::Quat, q2::Quat) = q1[:] == q2[:]

Base.:+(q1::Quat, q2::Quat) = Quat(q1[:] + q2[:])
Base.:+(q::Quat, a::Real) = +(promote(q, a)...)
Base.:+(a::Real, q::Quat) = +(promote(a, q)...)

Base.:-(q1::Quat, q2::Quat) = Quat(q1[:] - q2[:])
Base.:-(q::Quat, a::Real) = -(promote(q, a)...)
Base.:-(a::Real, q::Quat) = -(promote(a, q)...)

function Base.:*(q1::Quat, q2::Quat)
    p_real = q1.real * q2.real - dot(q1.imag, q2.imag)
    p_imag = q1.real * q2.imag + q2.real * q1.imag + cross(q1.imag, q2.imag)
    Quat([p_real, p_imag...])
end
Base.:*(q::Quat, a::Real) = a * q
Base.:*(a::Real, q::Quat) = Quat(a * q[:])

Base.:/(q1::Quat, q2::Quat) = q1 * inv(q2)
Base.:/(q::Quat, a::Real) = Quat(q[:] / a)
Base.:/(a::Real, q::Quat) = /(promote(a, q)...)

Base.:\(q1::Quat, q2::Quat) = inv(q1) * q2 #!= /(q2, q1) == q2 * inv(q1)
Base.:\(q::Quat, a::Real) = \(promote(q, a)...)
Base.:\(a::Real, q::Quat) = q / a


######################## UnitQuat #############################

struct UnitQuat <: AbstractQuat
    quat::Quat
    #restrict inner constructor to the declared field types
    function UnitQuat(quat::AbstractVector; enforce_norm::Bool = true)
        # println("Normalization $enforce_norm")
        quat_copy = copy(quat)
        enforce_norm && normalize!(quat_copy)
        return new(quat_copy) #tries to convert quat to Quat using Quat convert methods
    end
end

#outer constructors
UnitQuat(s::Real) = UnitQuat([s, 0, 0, 0])
UnitQuat() = UnitQuat([1, 0, 0, 0], enforce_norm = false)
UnitQuat(; real::Real, imag::AbstractVector{T}) where {T<:Real} = UnitQuat([real, imag...])

#bypass normalization on copy
Base.copy(u::UnitQuat) = UnitQuat(copy(getfield(u, :quat)), enforce_norm = false)
Base.getindex(u::UnitQuat, i) = (getfield(u, :quat)[i])
Base.setindex!(::UnitQuat, v, i) = error(
    "UnitQuat: Directly setting components not allowed, cast to Quat first")
Base.setproperty!(::UnitQuat, ::Symbol, v) = error(
    "UnitQuat: Directly setting real and imaginary parts not allowed, cast to Quat first")

LinearAlgebra.normalize!(u::UnitQuat) = (normalize!(getfield(u, :quat)), return u)
LinearAlgebra.normalize(u::UnitQuat) = UnitQuat(normalize(getfield(u, :quat)))

Base.promote_rule(::Type{UnitQuat}, ::Type{Quat}) = Quat
Base.convert(::Type{UnitQuat}, a::AbstractVector) = UnitQuat(a)

#### Adjoint & Inverse
Base.conj(u::UnitQuat)= UnitQuat([u.real, -u.imag...])
Base.adjoint(u::UnitQuat) = conj(u)
Base.inv(u::UnitQuat) = u'

#### Operators
Base.:+(u::UnitQuat) = u
Base.:-(u::UnitQuat) = UnitQuat(-u[:])

#(==) is inherited from AbstractVector, but will return true for any
#AbstractVector as long as it matches u[:], to avoid it we need to define these:
Base.:(==)(v1::AbstractVector, u2::UnitQuat) = false
Base.:(==)(u1::UnitQuat, v2::AbstractVector) = false
Base.:(==)(u1::UnitQuat, q2::Quat) = ==(promote(u1, q2)...)
Base.:(==)(q1::Quat, u2::UnitQuat) = ==(promote(q1, u2)...)
Base.:(==)(u1::UnitQuat, u2::UnitQuat) = u1[:] == u2[:]

Base.:*(u1::UnitQuat, u2::UnitQuat) = UnitQuat(getfield(u1, :quat) * getfield(u2, :quat), enforce_norm = false)
Base.:*(u::UnitQuat, q::Quat) = *(promote(u, q)...)
Base.:*(q::Quat, u::UnitQuat) = *(promote(q, u)...)

Base.:/(u1::UnitQuat, u2::UnitQuat) = u1 * inv(u2)
Base.:/(u::UnitQuat, q::Quat) = /(promote(u, q)...)
Base.:/(q::Quat, u::UnitQuat) = /(promote(q, u)...)

Base.:\(u1::UnitQuat, u2::UnitQuat) = inv(u1) * u2 #!= /(u2, u1) == u2 * inv(u1)
Base.:\(u::UnitQuat, q::Quat) = \(promote(u, q)...)
Base.:\(q::Quat, u::UnitQuat) = \(promote(q, u)...)

end #module