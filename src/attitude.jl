module Attitude

using StaticArrays: SVector, SMatrix
using LinearAlgebra
using Flight.Quaternions: UnitQuat, Quat
#if we didn't want to rely on Flight using Quaternions, we could do simply:
#using ..Quaternions: UnitQuat, Quat
#however, assumes a specific hierarchy

export Rotation, RQuat, RAxAng, REuler, RMatrix, Rx, Ry, Rz, dt

const ε_null = 1e-10 #threshold for null rotation
const half_π = π/2

function skew(v::AbstractVector{T} where T<:Real)
    v = SVector{3,Float64}(v)
    vx = SMatrix{3,3,Float64}([
        0      -v[3]    v[2];
        v[3]    0      -v[1];
       -v[2]    v[1]    0])
    return vx
end

######################### Rotation ###########################

#implements a passive (alias) rotation, that is, a rotation describing the
#relative attitude between two reference frames.

#TODO: Add display methods

abstract type Rotation end


############################# RQuat ###############################

struct RQuat <: Rotation
    _quat::UnitQuat
end

RQuat(r::Rotation) = convert(RQuat, r)
RQuat() = RQuat(UnitQuat())

Base.getindex(r::RQuat, i) = getindex(r._quat, i)
Base.length(r::RQuat) = length(r._quat)
Base.iterate(r::RQuat, state = 1) = iterate(r._quat, state)
Base.show(io::IO, r::RQuat) = print(io, "$(typeof(r))($(r[:]))")

Base.:(==)(r1::RQuat, r2::RQuat) = (r1._quat == r2._quat || r1._quat == -r2._quat)
Base.:(≈)(r1::RQuat, r2::RQuat) = (r1._quat ≈ r2._quat || r1._quat ≈ -r2._quat)
Base.:∘(r1::RQuat, r2::RQuat) = RQuat(r1._quat * r2._quat)
Base.adjoint(r::RQuat) = RQuat(r._quat')
function Base.:*(r_ab::RQuat, v_b::AbstractVector{T} where {T<:Real})::SVector{3,Float64}
    #this conversion yields a threefold speed gain
    v_b = SVector{3, Float64}(v_b)
    # v_b = SVector{3, Float64}(v_b) #this causes type instability between SV3
    # and Vector, but does not hurt performance
    q = r_ab._quat; q_re = q.real; q_im = q.imag
    v_a = v_b + 2q_im × (q_re * v_b + q_im × v_b)
    return v_a
end

LinearAlgebra.norm(r::RQuat) = norm(r._quat)
LinearAlgebra.normalize(r::RQuat) = RQuat(normalize(r._quat))
# LinearAlgebra.normalize!(r::RQuat) = (r._quat = normalize(r._quat))

dt(r_ab::RQuat, ω_ab_b::AbstractVector{T} where {T<:Real}) = 0.5 * (r_ab._quat * Quat(imag=ω_ab_b))

#require each Rotation subtype to implement conversions to and from RQuat
Base.convert(::Type{RQuat}, r::R) where {R<:Rotation} = error("Implement $R to RQuat conversion")
Base.convert(::Type{R}, r::RQuat) where {R<:Rotation} = error("Implement RQuat to $R conversion")
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


######################### RMatrix #######################

struct RMatrix <: Rotation
    _mat::SMatrix{3, 3, Float64}
    function RMatrix(input::AbstractArray{T, 2} where {T<:Real}; normalization::Bool = true)
        return normalization ? new(qr(input).Q) : new(input)
    end
end

RMatrix(r::Rotation) = convert(RMatrix, r)
RMatrix() = RMatrix(SMatrix{3,3,Float64}(I), normalization = false)

Base.size(r::RMatrix) = size(getfield(r,:_mat))
Base.getindex(r::RMatrix, i...) = getindex(getfield(r,:_mat), i...)

Base.convert(::Type{RMatrix}, r::RMatrix) = r

function Base.convert(::Type{RMatrix}, r::RQuat)

    q = normalize(r._quat) #cheap
    q_sq = q[:] .* q[:]
    dq12 = 2*q[1]*q[2]; dq13 = 2*q[1]*q[3]; dq14 = 2*q[1]*q[4];
    dq23 = 2*q[2]*q[3]; dq24 = 2*q[2]*q[4]; dq34 = 2*q[3]*q[4];

    RMatrix([
    1 - 2*(q_sq[3]+q_sq[4])  dq23 - dq14                dq24 + dq13;
    dq23 + dq14              1 - 2*(q_sq[2]+q_sq[4])    dq34 - dq12;
    dq24 - dq13              dq34 + dq12                1 - 2*(q_sq[2]+q_sq[3])
    ], normalization = false)

end

function Base.convert(::Type{RQuat}, r::RMatrix)

    #find the maximum amongst the absolute values of quaternion components
    R = r._mat
    tr_R = tr(R)
    i_max = findmax([tr_R; diag(R)])[2]

    if i_max == 1
        return RQuat([
            1 + tr_R,
            R[3,2] - R[2,3],
            R[1,3] - R[3,1],
            R[2,1] - R[1,2]
        ])
    elseif i_max == 2
        return RQuat([
            R[3,2] - R[2,3],
            1 + 2R[1,1] - tr_R,
            R[1,2] + R[2,1],
            R[3,1] + R[1,3]
        ])
    elseif i_max == 3
        return RQuat([
            R[1,3] - R[3,1],
            R[1,2] + R[2,1],
            1 + 2R[2,2] - tr_R,
            R[2,3] + R[3,2]
        ])
    else #i_max == 4
        return RQuat([
            R[2,1] - R[1,2],
            R[3,1] + R[1,3],
            R[2,3] + R[3,2],
            1 + 2R[3,3] - tr_R
        ])
    end
end

#a rotation matrix can define its own equality, composition and inversion
#operators without relying on RQuat
Base.:(==)(r1::RMatrix, r2::RMatrix) = (r1._mat == r2._mat)
Base.:(≈)(r1::RMatrix, r2::RMatrix) = r1._mat ≈ r2._mat
Base.:*(r::RMatrix, v::AbstractVector{T} where {T<:Real}) = r._mat * v
Base.:∘(r1::RMatrix, r2::RMatrix) = RMatrix(r1._mat * r2._mat, normalization = false)
Base.:*(r1::RMatrix, r2::RMatrix) = r1 ∘ r2 #allow overload of * for composition
Base.adjoint(r::RMatrix) = RMatrix(r._mat', normalization = false)

#using StaticArrays, a QR factorization takes ~25 ns. it is basically free!
# assembling the RMatrix object adds 240 ns.
#normalizing a RQuat takes ~15ns.
LinearAlgebra.normalize(r::RMatrix) = RMatrix(r._mat, normalization = true)
# LinearAlgebra.normalize!(r::RMatrix) = (r._mat = qr(r._mat).Q; return r)
LinearAlgebra.det(r::RMatrix) = det(r._mat)

function dt(r_ab::RMatrix, ω_ab_b::AbstractVector{T} where {T<:Real})
    Ω_ab_b = skew(ω_ab_b)
    return r_ab._mat * Ω_ab_b
end


############################# RAxAng ###############################

struct RAxAng <: Rotation
    axis::SVector{3, Float64}
    angle::Float64
    function RAxAng(axis::AbstractVector{T} where {T<:Real}, angle::Real; normalization::Bool = true)
        axis = SVector{3,Float64}(axis) #normalization will be faster for an SVector
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