module Quaternions

using StaticArrays: SVector, MVector #this form of using does not allow method extension
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
# Base.iterate(q::AbstractQuat, state = 1) = (state > 4 ? nothing : (q[state], state + 1))
Base.eltype(::AbstractQuat) = Float64
Base.show(io::IO, q::AbstractQuat) = print(io, "$(typeof(q))($(q[:]))")
Base.show(io::IO, ::MIME"text/plain", q::AbstractQuat) = print(io, "$(typeof(q))($(q[:]))")

#real and imaginary parts
Base.getindex(q::AbstractQuat, s::Symbol) = getindex(q, Val(s))
Base.getproperty(q::AbstractQuat, s::Symbol) = getindex(q, Val(s))
Base.setproperty!(q::AbstractQuat, s::Symbol, v) = setindex!(q, v, Val(s))

######################## Quat #############################

struct Quat <: AbstractQuat
    _sv::SVector{4, Float64}
    #whenever typeof(input) != QData, new will call convert(QData, input). as long
    #as QData provides a convert method that can handle typeof(input), then we do
    #not need to handle it explicitly in an outer constructor.
end

#outer constructors
Quat(s::Real) = Quat(SVector{4,Float64}(s, 0, 0, 0))
Quat(; real = 0.0, imag = zeros(3)) = Quat(SVector{4,Float64}(real, imag[1], imag[2], imag[3]))
Quat(q::AbstractQuat) = Quat(q[:])

Base.copy(q::Quat) = Quat(copy(getfield(q, :_sv)))
Base.getindex(q::Quat, i) = getfield(q, :_sv)[i]
Base.getindex(q::Quat, ::Val{:real}) = getfield(q, :_sv)[1]
Base.getindex(q::Quat, ::Val{:imag}) = SVector{3, Float64}(@view getfield(q, :_sv)[2:4])
#to avoid allocation, use @view instead of:
# Base.getindex(q::Quat, ::Val{:imag}) = SVector{3, Float64}(q[2:4])

LinearAlgebra.norm(q::Quat) = norm(getfield(q, :_sv)) #uses StaticArrays implementation
LinearAlgebra.normalize(q::Quat) = Quat(normalize(getfield(q, :_sv)))
norm_sqr(q::Quat) = (data = getfield(q,:_sv); sum(data.*data))

Base.promote_rule(::Type{Quat}, ::Type{S}) where {S<:Real} = Quat
Base.convert(::Type{Quat}, a::Real) = Quat(a)
Base.convert(::Type{Quat}, v::Union{AbstractQuat, AbstractVector}) = Quat(v[:])
Base.convert(::Type{Quat}, q::Quat) = q #if already a Quat, don't do anything

#### Adjoint & Inverse
Base.conj(q::Quat)= Quat([q.real, -q.imag...])
Base.adjoint(q::Quat) = conj(q)
Base.inv(q::Quat) = Quat(getfield(q', :_sv) / norm_sqr(q))

#### Operators
Base.:+(q::Quat) = q
Base.:-(q::Quat) = Quat(-getfield(q, :_sv))

Base.:(==)(q1::Quat, q2::Quat) = getfield(q1,:_sv) == getfield(q2,:_sv)
Base.:(≈)(q1::Quat, q2::Quat) = getfield(q1,:_sv) ≈ getfield(q2,:_sv)

Base.:+(q1::Quat, q2::Quat) = Quat(getfield(q1,:_sv) + getfield(q2,:_sv))
Base.:+(q::Quat, a::Real) = +(promote(q, a)...)
Base.:+(a::Real, q::Quat) = +(promote(a, q)...)

Base.:-(q1::Quat, q2::Quat) = Quat(getfield(q1,:_sv) - getfield(q2,:_sv))
Base.:-(q::Quat, a::Real) = -(promote(q, a)...)
Base.:-(a::Real, q::Quat) = -(promote(a, q)...)

function Base.:*(q1::Quat, q2::Quat)

    q1_re = q1[Val(:real)]; q2_re = q2[Val(:real)]
    q1_im = q1[Val(:imag)]; q2_im = q2[Val(:imag)]

    p_re = q1_re * q2_re - q1_im ⋅ q2_im
    p_im = q1_re * q2_im + q2_re * q1_im + q1_im × q2_im

    # Quat([p_re, p_im...]) #much slower!
    Quat(vcat(p_re, p_im))
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
    _q::Quat
    function UnitQuat(input::AbstractVector; normalization::Bool = true)
        v = SVector{4,Float64}(input)
        return normalization ? new(normalize(v)) : new(v)
    end
end

#outer constructors
UnitQuat(::Real) = UnitQuat(SVector{4,Float64}(1.0, 0, 0, 0), normalization = false)
function UnitQuat(; real::Union{Real, Nothing} = nothing,
                    imag::Union{AbstractVector{T} where {T<:Real}, Nothing} = nothing,
                    normalization::Bool = true)
    if imag === nothing
        return UnitQuat(1)
    elseif real === nothing
        return UnitQuat([0, imag...], normalization = normalization)
    else
        return UnitQuat([real, imag...], normalization = normalization)
    end
end
UnitQuat(q::AbstractQuat) = UnitQuat(q[:])

#bypass normalization on copy
Base.copy(u::UnitQuat) = UnitQuat(copy(getfield(u, :_q)), normalization = false) #saves normalization
Base.getindex(u::UnitQuat, i) = (getfield(u, :_q)[i])
Base.getindex(u::UnitQuat, ::Val{:real}) = getindex(getfield(u, :_q), Val(:real))
Base.getindex(u::UnitQuat, ::Val{:imag}) = getindex(getfield(u, :_q), Val(:imag))

LinearAlgebra.norm(u::UnitQuat) = norm(getfield(u, :_q)) #uses StaticArrays implementation
LinearAlgebra.normalize(u::UnitQuat) = UnitQuat(normalize(getfield(u, :_q)), normalization = false)
# LinearAlgebra.normalize!(u::UnitQuat) = (setfield!(u, :_q, normalize(getfield(u, :_q))); return u)

Base.promote_rule(::Type{UnitQuat}, ::Type{Quat}) = Quat
Base.convert(::Type{UnitQuat}, u::UnitQuat) = u #do not normalize on convert
Base.convert(::Type{UnitQuat}, v::Union{AbstractQuat, AbstractVector}) = UnitQuat(v[:])

#### Adjoint & Inverse
Base.conj(u::UnitQuat)= UnitQuat([u.real, -u.imag...], normalization = false)
Base.adjoint(u::UnitQuat) = conj(u)
Base.inv(u::UnitQuat) = u'

#### Operators
Base.:+(u::UnitQuat) = u
Base.:-(u::UnitQuat) = UnitQuat(-getfield(u, :_q), normalization = false)

Base.:(==)(u1::UnitQuat, q2::Quat) = ==(promote(u1, q2)...)
Base.:(==)(q1::Quat, u2::UnitQuat) = ==(promote(q1, u2)...)
Base.:(==)(u1::UnitQuat, u2::UnitQuat) = getfield(u1,:_q) == getfield(u2,:_q)

Base.:(≈)(u1::UnitQuat, q2::Quat) = ≈(promote(u1, q2)...)
Base.:(≈)(q1::Quat, u2::UnitQuat) = ≈(promote(q1, u2)...)
Base.:(≈)(u1::UnitQuat, u2::UnitQuat) = getfield(u1,:_q) ≈ getfield(u2,:_q)

Base.:*(u1::UnitQuat, u2::UnitQuat) = UnitQuat(getfield(u1, :_q) * getfield(u2, :_q), normalization = false)
Base.:*(u::UnitQuat, q::Quat) = *(promote(u, q)...)
Base.:*(q::Quat, u::UnitQuat) = *(promote(q, u)...)

Base.:/(u1::UnitQuat, u2::UnitQuat) = u1 * inv(u2)
Base.:/(u::UnitQuat, q::Quat) = /(promote(u, q)...)
Base.:/(q::Quat, u::UnitQuat) = /(promote(q, u)...)

Base.:\(u1::UnitQuat, u2::UnitQuat) = inv(u1) * u2 #!= /(u2, u1) == u2 * inv(u1)
Base.:\(u::UnitQuat, q::Quat) = \(promote(u, q)...)
Base.:\(q::Quat, u::UnitQuat) = \(promote(q, u)...)

end #module