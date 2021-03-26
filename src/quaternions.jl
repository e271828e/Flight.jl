
#TODO: write missing unit tests for operators, read up on testing facilities
#TODO: prettify functions with subindices and symbols
# q₁.re * q₂.re - q₁.im · q.im
# q₁⊙q₂
# × q1.imag
#TODO: generate getproperty with code generation from a simple macro that receives
# properties and their associated expressions and expands into an if sequence
# properties = Dict(:real => argname[1], :imag => argname[2:4])

module Quaternions

using StaticArrays: MVector #this form of using does not allow method extension
using LinearAlgebra #this form does, with LinearAlgebra.method

export Quat, UnitQuat, AbstractQuat
export Quat16, Quat32, Quat64
export UnitQuat16, UnitQuat32, UnitQuat64

######################## AbstractQuat #############################

abstract type AbstractQuat{T<:AbstractFloat} <: AbstractVector{T} end

#indexing and iterable interfaces; see https://docs.julialang.org/en/v1/manual/interfaces/
Base.size(::AbstractQuat) = (4,)
Base.length(::AbstractQuat) = 4
Base.firstindex(::AbstractQuat) = 1
Base.lastindex(::AbstractQuat) = 4
Base.getindex(::AbstractQuat, i) = error("AbstractQuat: getindex not implemented")
Base.setindex!(::AbstractQuat, v, i) = error("AbstractQuat: setindex! not implemented")
#not needed, inherited from AbstractVector
# Base.iterate(q::AbstractQuat, state = 1) = (state > 4 ? nothing : (q[state], state + 1))
Base.eltype(::AbstractQuat{T}) where {T} = T #helps with allocation efficiency

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
function norm_sqr(q::AbstractQuat)
    println("Called custom norm squared")
    sum(abs2.(q))
end
LinearAlgebra.norm(q::AbstractQuat) = norm(q[:]) #uses StaticArrays implementation

#not needed, inherited from AbstractVector
# Base.:+(q::AbstractQuat) = q
# Base.:-(q::AbstractQuat) = (q[:] = -q[:]; return q)
# Base.:(==)(q1::T, q2::T) where {T<:AbstractQuat} = (q1[:] == q2[:])
# Base.:(==)(q1::AbstractQuat, q2::AbstractQuat) = ==(promote(q1, q2)...)

######################## Quat #############################

QData{T} = MVector{4, T}

struct Quat{T} <: AbstractQuat{T}
    data::QData{T}
    #new() implicitly tries to convert() its inputs to the declared field types,
    #maybe it's clearer to do it explicitly upstream. thus, require the inner
    #constructor to accept only arguments of the exact declared type and let the
    #outer constructors handle conversion explicitly.
    Quat{T}(data::QData{T}) where {T} = new{T}(data)
end

Quat16 = Quat{Float16}
Quat32 = Quat{Float32}
Quat64 = Quat{Float64}

#main outer constructor. it would be generated automatically if no explicit
#inner constructor were provided
Quat(data::QData{T}) where {T} = Quat{T}(data)

#explicit type parameter
Quat{T}(v::AbstractVector) where {T} = Quat{T}(QData{T}(v))
Quat{T}(q::AbstractQuat) where {T} = Quat{T}(QData{T}(q[:]))
Quat{T}(s::Real) where {T} = Quat{T}(QData{T}(s, 0, 0, 0))

#inferred type parameter
Quat(x::Union{AbstractVector{T}, AbstractQuat{T}, T}) where {T<:AbstractFloat} = Quat{T}(x)

#real and imaginary kwargs
function Quat{T}(; kwargs...) where {T}
    # println("Called kwargs const")
    kwargs = Dict(kwargs)
    if haskey(kwargs, :real) && !haskey(kwargs, :imag)
        data = QData{T}(kwargs[:real], 0, 0, 0)
    elseif !haskey(kwargs, :real) && haskey(kwargs, :imag)
        data = QData{T}(0, kwargs[:imag]...)
    elseif haskey(kwargs, :real) && haskey(kwargs, :imag)
        data = QData{T}(kwargs[:real], kwargs[:imag]...)
    else
        data = QData{T}(zeros(T,4))
    end
    return Quat{T}(data)
end
#if no type parameter is provided with keyword constructor, default to Float64
Quat(; kwargs...) = Quat{Float64}(; kwargs...) #semicolon is essential here!

#basics
Base.getindex(q::Quat, i) = (getfield(q, :data)[i])
Base.setindex!(q::Quat, v, i) = (getfield(q, :data)[i] = v)

#### Promotion
Base.promote_rule(::Type{Quat{T}}, ::Type{Quat{S}}) where {T, S} = Quat{promote_type(T,S)}
#prioritize quaternion type parameter
Base.promote_rule(::Type{Quat{T}}, ::Type{S}) where {T, S<:Real} = Quat{T}
#prioritize better precision between quaternion type parameter and real subtype
# Base.promote_rule(::Type{Quat{T}}, ::Type{S}) where {T, S<:Real} = Quat{promote_type(T,S)}

#### Conversion
Base.convert(::Type{Quat{T}}, q::Quat) where {T} = Quat{T}(q)
Base.convert(::Type{Quat{T}}, a::Real) where {T} = Quat{T}(a)

#### Functions
Base.copy(q::T) where{T<:Quat} = T(q)
Base.conj(q::T) where {T<:Quat} = T([q.real, -q.imag...])
Base.adjoint(q::T) where {T<:Quat} = Base.conj(q)
Base.inv(q::T) where {T<:Quat} = T(q'[:] / norm_sqr(q))
LinearAlgebra.normalize!(q::Quat) = (normalize!(getfield(q, :data)); return q)

#### Operators
Base.:+(q1::T, q2::T) where {T<:Quat} = T(q1[:] + q2[:])
Base.:+(q1::Quat, q2::Quat) = +(promote(q1, q2)...)
Base.:+(q::Quat, a::Real) = +(promote(q, a)...)
Base.:+(a::Real, q::Quat) = +(promote(a, q)...)

Base.:-(q1::T, q2::T) where {T<:Quat} = T(q1[:] - q2[:])
Base.:-(q1::Quat, q2::Quat) = -(promote(q1, q2)...)
Base.:-(q::Quat, a::Real) = -(promote(q, a)...)
Base.:-(a::Real, q::Quat) = -(promote(a, q)...)

function Base.:*(q1::T, q2::T) where {T<:Quat}
    p_real = q1.real * q2.real - dot(q1.imag, q2.imag)
    p_imag = q1.real * q2.imag + q2.real * q1.imag + cross(q1.imag, q2.imag)
    T([p_real, p_imag...])
end
Base.:*(q1::Quat, q2::Quat) = *(promote(q1, q2)...)
Base.:*(a::Real, q::Quat) = *(promote(a, q)...)
Base.:*(q::Quat, a::Real) = *(promote(q, a)...)

Base.:/(q1::T, q2::T) where {T<:Quat} = q1 * inv(q2)
Base.:/(q1::Quat, q2::Quat) = /(promote(q1, q2)...)
Base.:/(a::Real, q::Quat) = /(promote(a, q)...)
Base.:/(q::Quat, a::Real) = /(promote(q, a)...)

Base.:\(q1::Quat, q2::Quat) = /(q2, q1)
Base.:\(a::Real, q::Quat) = /(q, a)
Base.:\(q::Quat, a::Real) = /(a, q)

#multiplication and division by scalar could be implemented more efficiently
#without promotion to Quat, but that makes the outcome when Number != eltype(q)
#harder to control

######################## UnitQuat #############################

struct UnitQuat{T} <: AbstractQuat{T}
    quat::Quat{T}
    #restrict inner constructor to the declared field types
    function UnitQuat{T}(quat::Quat{T}; enforce_norm::Bool = true) where {T}
        enforce_norm && normalize!(quat)
        return new{T}(quat)
    end
end

UnitQuat16 = UnitQuat{Float16}
UnitQuat32 = UnitQuat{Float32}
UnitQuat64 = UnitQuat{Float64}

#explicit type parameter
UnitQuat{T}(v::AbstractVector) where {T} = UnitQuat{T}(Quat{T}(v))
UnitQuat{T}(q::AbstractQuat) where {T} = UnitQuat{T}(Quat{T}(q))
UnitQuat{T}(args...; kwargs...) where {T} = UnitQuat(Quat{T}(args...; kwargs...))

#implicit type parameter; let the Quat constructor infer it...
UnitQuat(args...; kwargs...) = UnitQuat(Quat(args...; kwargs...))
#... then dispatch to the UnitQuat inner constructor using the inferred parameter
UnitQuat(quat::Quat{T}) where {T} = UnitQuat{T}(quat)

#handle the zero-argument case explicitly, because Quat defaults to the null
#quaternion, for which normalization yields NaNs
UnitQuat{T}(::Vararg{Any,0}) where {T} = UnitQuat{T}(Quat{T}(1.0), enforce_norm = false)
UnitQuat(::Vararg{Any,0}) = UnitQuat{Float64}(Quat{Float64}(1.0), enforce_norm = false)

#basics
#bypass normalization on copy
Base.getindex(u::UnitQuat, i) = (getfield(u, :quat)[i])
Base.setindex!(::UnitQuat, v, i) = error(
    "UnitQuat: Directly setting components not allowed, cast to Quat first")
Base.setproperty!(::UnitQuat, ::Symbol, v) = error(
    "UnitQuat: Directly setting real and imaginary parts not allowed, cast to Quat first")

Base.promote_rule(::Type{UnitQuat{T}}, ::Type{UnitQuat{S}}) where {T, S} = UnitQuat{promote_type(T,S)}
Base.promote_rule(::Type{UnitQuat{T}}, ::Type{Quat{S}}) where {T, S} = Quat{promote_type(T,S)}
# another option, which may downcast UnitQuat's type parameter to that of Quat:
#Base.promote_rule(::Type{UnitQuat{T}}, ::Type{Quat{S}}) where {T, S} = Quat{S}

#the inner constructor normalizes when downcasting or upcasting from another UnitQuat
Base.convert(::Type{UnitQuat{T}}, u::UnitQuat) where {T} = UnitQuat{T}(u)
Base.convert(::Type{Quat{T}}, u::UnitQuat{S}) where {T, S}  = Quat{T}(UnitQuat{T}(u))

#functions
Base.copy(u::UnitQuat{T}) where{T} = UnitQuat{T}(Quat{T}(u), enforce_norm = false)
Base.adjoint(u::UnitQuat{T}) where {T} = UnitQuat{T}(Quat{T}([u.real, -u.imag...]))
Base.inv(u::UnitQuat) = u'
LinearAlgebra.normalize!(u::UnitQuat) = (normalize!(getfield(u, :quat)))

#operators
Base.:*(u1::T, u2::T) where {T<:UnitQuat} = T(getfield(u1, :quat) * getfield(u2, :quat))
Base.:*(u1::UnitQuat, u2::UnitQuat) = *(promote(u1, u2)...)
Base.:*(u::UnitQuat, q::Quat) = *(promote(u, q)...)
Base.:*(q::Quat, u::UnitQuat) = *(promote(q, u)...)

Base.:/(u1::T, u2::T) where {T<:UnitQuat} = u1 * inv(u2)
Base.:/(u1::UnitQuat, u2::UnitQuat) = /(promote(u1, u2)...)
Base.:/(u::UnitQuat, q::Quat) = /(promote(u, q)...)
Base.:/(q::Quat, u::UnitQuat) = /(promote(q, u)...)

Base.:\(u1::UnitQuat, u2::UnitQuat) = /(u2, u1)
Base.:\(u::UnitQuat, q::Quat) = /(q, u)
Base.:\(q::Quat, u::UnitQuat) = /(u, q)

end #module