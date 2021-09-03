"""
Lightweight 3D attitude representation module.

Defines the abstract type `Rotation`, and implements four concrete subtypes:
- Unit quaternion (`RQuat`)
- Rotation matrix (`RMatrix`)
- Axis-angle (`RAxAng`)
- Euler angles (`REuler`)

All of them support the essential attitude operations used in this package:
- Equality (`==`) and approximate equality (`≈`)
- Inversion (`r'`)
- Composition (`r1 ∘ r2`)
- Vector coordinate transformation (`r * v == r(v)`)
- Conversion to and from other representations (`r1 = R1(r2)`)

`RQuat` operates as a primary representation. It implements all of the
operations natively, and therefore it can be used for promotion whenever a
specific operation cannot be directly performed on its original operand(s).

This allows for a common interface; any `Rotation` subtype need only implement
conversions to and from `RQuat` to be fully compatible with all other subtypes
and operations. However, for efficiency reasons, some direct methods are
provided for the implemented subtypes.
"""
module Attitude

using StaticArrays
using StructArrays
using LinearAlgebra
using Flight.Quaternions
using Flight.Plotting
#using ..Quaternions: UnitQuat, Quat #also works (but relies on folder hierarchy)

export Rotation, RQuat, RAxAng, REuler, RMatrix, Rx, Ry, Rz, dt

const half_π = π/2

"""
Compute the skew-symmetric matrix corresponding to a 3-element vector.
"""
function skew(v::AbstractVector)

    @SMatrix [0 -v[3] v[2];
            v[3] 0 -v[1];
            -v[2] v[1] 0]

end

######################### Rotation ###########################

"""
Generic 3D rotation descriptor
"""
abstract type Rotation end

"Function call notation for coordinate transformations"
(r::Rotation)(v::AbstractVector{<:Real}) = r * v

############################# RQuat ###############################

"Unit quaternion representation"
struct RQuat <: Rotation
    _u::UnitQuat #no normalization is performed for an UnitQuat input
end

#generic AbstractVectors are normalized by default, this can be overridden with
#the optional keyword argument
function RQuat(v::AbstractVector{<:Real}; normalization::Bool = true)
    RQuat(UnitQuat(SVector{4,Float64}(v), normalization = normalization))
end
RQuat() = RQuat(UnitQuat(1.0))
RQuat(r::Rotation) = convert(RQuat, r)

Base.getindex(r::RQuat, i) = getindex(r._u, i)
Base.length(r::RQuat) = length(r._u)
Base.size(r::RQuat) = size(r._u)
# Base.iterate(r::RQuat, state = 1) = iterate(r._u, state)

# direct iteration unsupported because it allocates. doing x.=r[:] instead of
#x.=r gets the underlying SVector first by forwarding getindex, then iterates on
#the StateVector itself. this avoids allocations.

Base.show(io::IO, r::RQuat) = print(io, "$(typeof(r))($(r[:]))")

#Strict equality (accounts for the unit quaternions' double cover of SO(3))
Base.:(==)(r1::RQuat, r2::RQuat) = (r1._u == r2._u || r1._u == -r2._u)
# Approximate equality (accounts for the unit quaternions' double cover of SO(3))
Base.:(≈)(r1::RQuat, r2::RQuat) = (r1._u ≈ r2._u || r1._u ≈ -r2._u)
# Composition
Base.:∘(r1::RQuat, r2::RQuat) = RQuat(r1._u * r2._u)
# Inversion
Base.adjoint(r::RQuat) = RQuat(r._u', normalization = false)

# Coordinate transformation (function call notation also supported)
function Base.:*(r_ab::RQuat, v_b_in::AbstractVector{<:Real})::SVector{3,Float64}
    v_b = SVector{3, Float64}(v_b_in) #casting to SVector yields a 3x speed gain
    q = r_ab._u; q_re = q.real; q_im = q.imag
    v_a = v_b + 2q_im × (q_re * v_b + q_im × v_b)
    return v_a
end

LinearAlgebra.norm(r::RQuat) = norm(r._u)
LinearAlgebra.normalize(r::RQuat) = RQuat(normalize(r._u))

"""
    dt(r_ab::RQuat, ω_ab_b::AbstractVector{<:Real})

Unit quaternion time derivative

If `r_ab` represents the rotation from axes εa to axes εb, and `ω_ab_b` is the
angular velocity of εb with respect to εa projected in εb, then `ṙ_ab ==
dt(r_ab, ω_ab_b)`
"""
dt(r_ab::RQuat, ω_ab_b::AbstractVector{<:Real}) = 0.5 * (r_ab._u * Quat(imag=ω_ab_b))

##### RQuat fallbacks #####

#require each Rotation subtype to implement conversions to and from RQuat
Base.convert(::Type{RQuat}, r::R) where {R<:Rotation} = error("Implement $R to RQuat conversion")
Base.convert(::Type{R}, r::RQuat) where {R<:Rotation} = error("Implement RQuat to $R conversion")

#this enables a two-step conversion via RQuat, but requires defining trivial
#conversions to avoid stack overflow
Base.convert(::Type{R}, r::Rotation) where {R<:Rotation} = convert(R, RQuat(r))
Base.convert(::Type{R}, r::R) where {R<:Rotation} = r
Base.convert(::Type{RQuat}, r::RQuat) = r

#unless the representation defines its own methods, fall back to RQuat;
#promote rules don't work here, because promoting two inputs of the same subtype
#leaves them unchanged, and here we need them converted to RQuat
Base.adjoint(r::R) where {R<:Rotation} = convert(R, RQuat(r)')
Base.:*(r::Rotation, v::AbstractVector{<:Real}) = RQuat(r) * v
Base.:(≈)(r1::Rotation, r2::Rotation) = RQuat(r1) ≈ RQuat(r2)
Base.:∘(r1::Rotation, r2::Rotation) = RQuat(r1) ∘ RQuat(r2)

#absolute equality is not defined in general, because comparing two different
#subtypes requires promotion to RQuat, whose outcome is affected by floating
#point precision

######################### RMatrix #######################

"Rotation matrix representation"
struct RMatrix <: Rotation
    _mat::SMatrix{3, 3, Float64, 9}
    function RMatrix(input::AbstractArray{<:Real, 2}; normalization::Bool = true)
        return normalization ? new(qr(input).Q) : new(input)
    end
end

#normalize by default
RMatrix() = RMatrix(SMatrix{3,3,Float64}(I), normalization = false)
RMatrix(r::Rotation) = convert(RMatrix, r)

Base.size(r::RMatrix) = size(getfield(r,:_mat))
Base.getindex(r::RMatrix, i...) = getindex(getfield(r,:_mat), i...)

Base.convert(::Type{RMatrix}, r::RMatrix) = r

function Base.convert(::Type{RMatrix}, r::RQuat)

    q = normalize(r._u) #cheap
    q_sq = q[:] .* q[:]
    dq12 = 2*q[1]*q[2]; dq13 = 2*q[1]*q[3]; dq14 = 2*q[1]*q[4];
    dq23 = 2*q[2]*q[3]; dq24 = 2*q[2]*q[4]; dq34 = 2*q[3]*q[4];

    M = @SMatrix [
    1 - 2*(q_sq[3]+q_sq[4])  dq23 - dq14                dq24 + dq13;
    dq23 + dq14              1 - 2*(q_sq[2]+q_sq[4])    dq34 - dq12;
    dq24 - dq13              dq34 + dq12                1 - 2*(q_sq[2]+q_sq[3])
    ]

    RMatrix(M, normalization = false)

end

function Base.convert(::Type{RQuat}, r::RMatrix)
    R = r._mat
    tr_R = tr(R)
    i_max = findmax(vcat(tr_R, diag(R)))[2]

    #these 4-element vectors are all proportional to the the desired unit
    #quaternion's components. to minimize numerical errors, we choose the one
    #with the largest norm (indicated by i_max), and construct the quaternion by
    #normalizing it
    if i_max == 1
        v = SVector{4,Float64}(
            1 + tr_R,
            R[3,2] - R[2,3],
            R[1,3] - R[3,1],
            R[2,1] - R[1,2]
        )
    elseif i_max == 2
        v = SVector{4,Float64}(
            R[3,2] - R[2,3],
            1 + 2R[1,1] - tr_R,
            R[1,2] + R[2,1],
            R[3,1] + R[1,3]
        )
    elseif i_max == 3
        v = SVector{4,Float64}(
            R[1,3] - R[3,1],
            R[1,2] + R[2,1],
            1 + 2R[2,2] - tr_R,
            R[2,3] + R[3,2]
        )
    else #i_max == 4
        v = SVector{4, Float64}(
            R[2,1] - R[1,2],
            R[3,1] + R[1,3],
            R[2,3] + R[3,2],
            1 + 2R[3,3] - tr_R
        )
    end

    return RQuat(v, normalization = true)

end

#RMatrix also implements all operators natively
Base.:(==)(r1::RMatrix, r2::RMatrix) = (r1._mat == r2._mat)
Base.:(≈)(r1::RMatrix, r2::RMatrix) = r1._mat ≈ r2._mat

# Coordinate transformation
Base.:*(r::RMatrix, v::AbstractVector{<:Real}) = r._mat * v

# Composition
Base.:∘(r1::RMatrix, r2::RMatrix) = RMatrix(r1._mat * r2._mat, normalization = false)
Base.:*(r1::RMatrix, r2::RMatrix) = r1 ∘ r2 #allow overload of * for composition

# Inversion
Base.adjoint(r::RMatrix) = RMatrix(r._mat', normalization = false)

#let the constructor handle it
LinearAlgebra.normalize(r::RMatrix) = RMatrix(r._mat, normalization = true)

LinearAlgebra.det(r::RMatrix) = det(r._mat)

"""
    dt(r_ab::RMatrix, ω_ab_b::AbstractVector{<:Real})

Rotation matrix time derivative

If `R_ab` represents the rotation from axes εa to axes εb, and `ω_ab_b` is the
angular velocity of εb with respect to εa projected in εb, then `Ṙ_ab ==
dt(R_ab, ω_ab_b)`
"""
function dt(R_ab::RMatrix, ω_ab_b::AbstractVector{<:Real})
    Ω_ab_b = skew(ω_ab_b)
    return R_ab._mat * Ω_ab_b
end


############################# RAxAng ###############################

"Axis-angle representation"
struct RAxAng <: Rotation
    axis::SVector{3, Float64}
    angle::Float64
    function RAxAng(axis::AbstractVector{<:Real}, angle::Real; normalization::Bool = true)
        axis = SVector{3,Float64}(axis) #normalization will be faster for an SVector
        return normalization ? new(normalize(axis), angle) : new(axis, angle)
    end
end

RAxAng(r::Rotation) = convert(RAxAng, r)
RAxAng(::Nothing, ::Real; normalization::Bool = false) = RAxAng()
RAxAng(input::Tuple{Union{Nothing, AbstractVector{<:Real}}, Real}) = RAxAng(input...)

#shorthand constructors for rotations around X, Y, Z
Rz(ψ::Real = 0.0) = RAxAng(SVector{3}(0.0, 0.0, 1.0), ψ, normalization = false)
Ry(θ::Real = 0.0) = RAxAng(SVector{3}(0.0, 1.0, 0.0), θ, normalization = false)
Rx(φ::Real = 0.0) = RAxAng(SVector{3}(1.0, 0.0, 0.0), φ, normalization = false)

RAxAng() = Rx(0)

#conversions to and from RQuat
function Base.convert(::Type{RAxAng}, r::RQuat)
    q_re = r._u.real
    q_im = r._u.imag
    norm_im = norm(q_im)
    μ = 2atan(norm_im, q_re)
    u = (norm_im > 0 ? q_im / norm_im : nothing)
    return RAxAng(u, μ, normalization = false)
end

function Base.convert(::Type{RQuat}, r::RAxAng)
    u = r.axis
    μ = r.angle
    RQuat(UnitQuat(real = cos(0.5μ), imag = u*sin(0.5μ), normalization = false))
end

# Inversion (the only native operation we can easily define for axis-angle)
Base.adjoint(r::RAxAng) = RAxAng(r.axis, -r.angle)


############################# REuler #############################

"Euler angle representation (convention is ZYX)"
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

#conversions to and from RQuat
function Base.convert(::Type{REuler}, r::RQuat)
        q = r._u
        q_sq = q[:] .* q[:]

        ψ = atan( 2*(q[1]*q[4] + q[2]*q[3]), 1 - 2*(q_sq[3] + q_sq[4]))
        θ = asin( clamp(2*(q[1]*q[3] - q[2]*q[4]), -1, 1) )
        φ = atan( 2*(q[1]*q[2] + q[3]*q[4]), 1 - 2*(q_sq[2] + q_sq[3]))

        return REuler(ψ, θ, φ)
end

function Base.convert(::Type{RQuat}, r::REuler)
    Rz(r.ψ) ∘ Ry(r.θ) ∘ Rx(r.φ)
end

########################## Plotting #################################

#unless a more specialized method is defined, a TimeHistory{Rotations} is
#converted to REuler for plotting
@recipe function plot_rotation(th::TimeHistory{<:AbstractVector{<:Rotation}})

    v_euler = Vector{REuler}(undef, length(th.data))
    for i in 1:length(v_euler)
        v_euler[i] = REuler(th.data[i])
    end
    sa = StructArray(v_euler)
    data = hcat(sa.ψ, sa.θ, sa.φ)./π #plot as π factors

    label --> ["Heading" "Inclination" "Bank"]
    yguide --> [L"$\psi \ (\pi \ rad)$" L"$\theta \ (\pi \ rad)$" L"$\phi \ (\pi \ rad)$"]
    th_split --> :h #custom TimeHistory attribute

    return TimeHistory(th.t, data)

end


end