module Attitude

using StaticArrays: SVector
using LinearAlgebra
using ..Quaternions: UnitQuat, Quat

export Rotation, RQuat, RAxAng, REuler, Rx, Ry, Rz, invert, compose, transform, dt

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

Base.:(==)(r1::RQuat, r2::RQuat) = (r1._quat == r2._quat || r1._quat == -r2._quat)
Base.:(≈)(r1::RQuat, r2::RQuat) = (r1._quat ≈ r2._quat || r1._quat ≈ -r2._quat)
Base.:∘(r1::RQuat, r2::RQuat) = RQuat(r1._quat * r2._quat)
Base.adjoint(r::RQuat) = RQuat(r._quat')

function Base.:*(r_ab::RQuat, v_b::AbstractVector{T} where {T<:Real})
    q = r_ab._quat; q_re = q.real; q_im = q.imag
    v_a = v_b + 2q_im × (q_re * v_b + q_im × v_b)
    return v_a
end

LinearAlgebra.norm(r::RQuat) = norm(r._quat)
LinearAlgebra.normalize(r::RQuat) = RQuat(normalize(r._quat))
LinearAlgebra.normalize!(r::RQuat) = (r._quat = normalize(r._quat))

dt(r_ab::RQuat, ω_ab_b::AbstractVector{T} where {T<:Real}) = 0.5 * (r_ab._quat * Quat(imag=ω_ab_b))

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
Base.:*(r::Rotation, v::AbstractVector{T} where {T<:Real}) = RQuat(r) * v
Base.:∘(r1::Rotation, r2::Rotation) = RQuat(r1) ∘ RQuat(r2)
#cannot define absolute equality between two Rotation subtypes other than RQuat,
#because it requires promotion to RQuat, which itself is inaccurate

#only allow composition of Rotations
Base.:∘(r::Rotation, x::Any) = error("$(typeof(r)) ∘ $(typeof(x)) composition not allowed")
Base.:∘(x::Any, r::Rotation) = error("$(typeof(x)) ∘ $(typeof(r)) composition not allowed")


############################# RAxAng ###############################

struct RAxAng <: Rotation
    axis::SVector{3, Float64}
    angle::Float64
    function RAxAng(axis::AbstractVector{T} where {T<:Real}, angle::Real; normalization::Bool = true)
        return normalization ? new(normalize(axis), angle) : new(axis, angle)
    end
end

RAxAng(r::Rotation) = convert(RAxAng, r)
RAxAng(input::Tuple{Union{Nothing, AbstractVector{T} where T<:Real}, Real}) = RAxAng(input...)
RAxAng(::Nothing, ::Real; normalization::Bool = false) = RAxAng()

Rz(ψ::Real) = RAxAng([0,0,1], ψ, normalization = false)
Ry(θ::Real) = RAxAng([0,1,0], θ, normalization = false)
Rx(φ::Real) = RAxAng([1,0,0], φ, normalization = false)

RAxAng() = Rx(0)

function Base.convert(::Type{RAxAng}, r::RQuat)
    q_re = r._quat.real
    q_im = r._quat.imag
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
    ψ::Float64 #heading
    θ::Float64 #inclination
    φ::Float64 #bank
    function REuler(ψ, θ, φ)
        @assert abs(ψ<=π) "Heading must be within [-π, π]"
        @assert abs(θ<=half_π) "Inclination must be within [-π/2, π/2]"
        @assert abs(φ<=π) "Bank must be within [-π, π]"
        new(ψ, θ, φ)
    end
end

REuler(r::Rotation) = convert(REuler, r)
REuler(input::Tuple{Real, Real, Real}) = REuler(input...)
REuler(; ψ = 0, θ = 0, φ = 0) = REuler(ψ, θ, φ)

function Base.convert(::Type{REuler}, r::RQuat)
        q = r._quat
        q_sq = q[:] .* q[:]

        ψ = atan( 2*(q[1]*q[4] + q[2]*q[3]), 1 - 2*(q_sq[3] + q_sq[4]))
        θ = asin( clamp(2*(q[1]*q[3] - q[2]*q[4]), -1, 1) )
        φ = atan( 2*(q[1]*q[2] + q[3]*q[4]), 1 - 2*(q_sq[2] + q_sq[3]))

        return REuler(ψ, θ, φ)
end

function Base.convert(::Type{RQuat}, r::REuler)
    Rz(r.ψ) ∘ Ry(r.θ) ∘ Rx(r.φ)
end


end