module Attitude

using StaticArrays: SVector
using LinearAlgebra
using ..Quaternions: UnitQuat, Quat

export Rotation, RQuat, RAxAng, REuler, invert, compose, transform, dt

const ε_small = 1e-8 #threshold for small angle approximation
const ε_null = 1e-10 #threshold for null rotation
const half_π = π/2

######################### Rotation ###########################

#implements a passive (alias) rotation, that is, a rotation describing the
#relative attitude between two reference frames.

abstract type Rotation end


############################# RQuat ###############################

mutable struct RQuat <: Rotation
    _quat::UnitQuat
end

RQuat(r::Rotation) = convert(RQuat, r)
RQuat() = RQuat(UnitQuat())

Base.:(==)(q1::RQuat, q2::RQuat) = (q1._quat == q2._quat || q1._quat == -q2._quat)
Base.:(≈)(q1::RQuat, q2::RQuat) = (q1._quat ≈ q2._quat || q1._quat ≈ -q2._quat)
Base.:∘(q1::RQuat, q2::RQuat) = RQuat(q1._quat * q2._quat)
Base.adjoint(q::RQuat) = RQuat(q._quat')

function Base.:*(q_ab::RQuat, v_b::AbstractVector{T} where {T<:Real})
    q = q_ab._quat; q_re = q.real; q_im = q.imag
    v_a = v_b + 2q_im × (q_re * v_b + q_im × v_b)
    return v_a
end

LinearAlgebra.norm(q::RQuat) = norm(q._quat)
LinearAlgebra.normalize(q::RQuat) = RQuat(normalize(q._quat))
LinearAlgebra.normalize!(q::RQuat) = (q._quat = normalize(q._quat))

dt(q_ab::RQuat, ω_ab_b::AbstractVector{T} where {T<:Real}) = 0.5 * (q_ab._quat * Quat(imag=ω_ab_b))

#require each Rotation subtype to implement conversions to and from RQuat
Base.convert(::Type{RQuat}, r::R) where {R<:Rotation} = error("Implement RQuat to $R conversion")
Base.convert(::Type{R}, r::RQuat) where {R<:Rotation} = error("Implement $R to RQuat conversion")
#this enables a two-step conversion via RQuat as a fallback
Base.convert(::Type{R}, r::Rotation) where {R<:Rotation} = convert(R, RQuat(r))
#trivial conversions
Base.convert(::Type{R}, r::R) where {R<:Rotation} = r
Base.convert(::Type{RQuat}, r::RQuat) = r

#unless the representation defines its own methods, fall back to RQuat
Base.adjoint(r::R) where {R<:Rotation} = convert(R, RQuat(r)')
Base.:(≈)(r1::Rotation, r2::Rotation) = RQuat(r1) ≈ RQuat(r2)
Base.:(==)(r1::Rotation, r2::Rotation) = RQuat(r1) == RQuat(r2)
Base.:*(r::Rotation, v::AbstractVector{T} where {T<:Real}) = RQuat(r) * v
Base.:∘(r1::Rotation, r2::Rotation) = RQuat(r1) ∘ RQuat(r2)

#only allow composition of Rotations
Base.:∘(r::Rotation, x::Any) = error("$(typeof(r)) ∘ $(typeof(x)) composition not allowed")
Base.:∘(x::Any, r::Rotation) = error("$(typeof(x)) ∘ $(typeof(r)) composition not allowed")


############################# RAxAng ###############################

const RAxis = SVector{3, Float64}

struct RAxAng <: Rotation
    axis::RAxis
    angle::Float64
    function RAxAng(axis::AbstractVector{T} where {T<:Real}, angle::Real; normalization::Bool = true)
        return normalization ? new(normalize(axis), angle) : new(axis, angle)
    end
end

RAxAng(r::Rotation) = convert(RAxAng, r)
RAxAng(input::Tuple{Union{Nothing, AbstractVector{T} where T<:Real}, Real}) = RAxAng(input...)
RAxAng(::Nothing, ::Real; normalization::Bool = false) = RAxAng()
RAxAng() = RAxAng(RAxis([1, 0, 0]), 0, normalization = false)

function Base.convert(::Type{RAxAng}, q::RQuat)
    q_re = q._quat.real
    q_im = q._quat.imag
    norm_im = norm(q_im)
    μ = 2atan(norm_im, q_re)
    u = (norm_im > ε_null ? q_im / norm_im : nothing)
    return RAxAng(u, μ, normalization = false)
end

function Base.convert(::Type{RQuat}, r::RAxAng)
    u = r.axis
    μ = r.angle
    RQuat(UnitQuat(real = cos(0.5μ), imag = u*sin(0.5μ), normalization = false))
end

#the only ad-hoc operator we can easily define for axis-angle
Base.adjoint(r::RAxAng) = RAxAng(r.axis, -r.angle)


############################# REuler #############################

struct REuler <: Rotation
    ψ::Float64
    θ::Float64
    φ::Float64
    function REuler(ψ, θ, φ)
        @assert abs(ψ<=π) "Heading must be within [-π, π]"
        @assert abs(θ<=half_π) "Inclination must be within [-π/2, π/2]"
        @assert abs(φ<=π) "Bank must be within [-π, π]"
        new(ψ, θ, φ)
    end
end

REuler(r::Rotation) = convert(REuler, r)
REuler(input::Tuple{Real, Real, Real}) = REuler(input...)
REuler() = REuler(0, 0, 0)

function Base.getproperty(r::REuler, s::Symbol)
    if s == :heading
        return r.ψ
    elseif s == :inclination
        return r.θ
    elseif s == :bank
        return r.φ
    else
        return getfield(r, s)
    end
end

#MAKE SURE ENFORCE NORM IS ONLY CALLED WHEN NECESSARY!

#RAxAng, RVec and Euler should all promote to Quat.

#conversions to provide:
# RMat, RAxAng, RVec and Euler to and from RQuat
# RMat to and from Quat
# with that, can provide two-step conversion to and from RMat to RAxAng, RVec and
# Euler. maybe Euler could be direct.


#put println statements in each conversion method to know when it's been called

#a RMat can be constructed from a 3x3 Matrix. in that case, orthonormality is
#always enforced by applying QR factorization. if it comes from a conversion, it
#is not necessary

#create normalize! and normalize methods for RMat

#i can admit * of rotation matrices without converting to RQuat, and also
#transforming vectors! only use promote when absolutely required!
#define promote_rule(R, q) = q (always prefer q)

#go to Flight folder
#enter package manager
#activate .
#do not do activate Flight. in that case, the test and using commands do not
#work. they expect to be run from Flight's parent folder. why??

#must fail when passed a 1D array of size other than 3
#how do we do this? by defining StaticArrays types


end