"""
Lightweight 3D attitude representation module.

Defines the abstract type `Abstract3DRotation`, and implements four concrete
subtypes:
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

This allows for a common interface; any `Abstract3DRotation` subtype need only
implement conversions to and from `RQuat` to be fully compatible with all other
subtypes and operations. For efficiency, some direct methods are provided for
the implemented subtypes.
"""
module Attitude

using LinearAlgebra, StaticArrays, ComponentArrays, StructArrays
using Plots, LaTeXStrings, DataStructures

using FlightCore

using ..Quaternions

export Abstract3DRotation, RQuat, RAxAng, RVec, REuler, RMatrix, Rx, Ry, Rz
export azimuth, inclination, wrap_to_π

const half_π = π/2

"""
Compute the skew-symmetric matrix corresponding to a 3-element vector.
"""
function v2skew(v::AbstractVector)

    @assert length(v) == 3

    @SMatrix [0 -v[3] v[2];
            v[3] 0 -v[1];
            -v[2] v[1] 0]

end

######################### Abstract3DRotation ###########################

"""
Continuous 3D rotation descriptor
"""
abstract type Abstract3DRotation end

"Function call notation for coordinate transformations"
(r::Abstract3DRotation)(v::AbstractVector{<:Real}) = r * v

############################# RQuat ###############################

"Unit quaternion representation"
struct RQuat <: Abstract3DRotation
    _u::UnitQuat #no normalization is performed for an UnitQuat input
end

#generic AbstractVectors are normalized by default
function RQuat(v::AbstractVector{<:Real}; normalization::Bool = true)
    RQuat(UnitQuat(SVector{4,Float64}(v); normalization))
end
RQuat() = RQuat(UnitQuat(1.0), normalization = false)
RQuat(r::Abstract3DRotation) = convert(RQuat, r)

Base.getindex(r::RQuat, i) = getindex(r._u, i)
Base.length(r::RQuat) = length(r._u)
Base.size(r::RQuat) = size(r._u)
# Base.iterate(r::RQuat, state = 1) = iterate(getfield(getfield(r._u, :_q), :_sv)[:], state)

# direct iteration unsupported because it allocates. doing x.=r[:] instead of
#x.=r gets the underlying SVector first by forwarding getindex, then iterates on
#the StateVector itself. this avoids allocations.

Base.show(io::IO, r::RQuat) = print(io, "$(typeof(r))($(r[:]))")

#Strict equality (accounts for the unit quaternions' double cover of SO(3))
Base.:(==)(r1::RQuat, r2::RQuat) = (r1._u == r2._u || r1._u == -r2._u)
# Approximate equality (accounts for the unit quaternions' double cover of SO(3))
Base.:(≈)(r1::RQuat, r2::RQuat; kwargs...) = (≈(r1._u, r2._u; kwargs...) || ≈(r1._u, -(r2._u); kwargs...))
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
Quaternions.UnitQuat(r::RQuat) = r._u

"""
    dt(r_ab::RQuat, ω_ab_b::AbstractVector{<:Real})

Rotation quaternion time derivative

If `r_ab` represents the rotation from axes εa to axes εb, and `ω_ab_b` is the
angular velocity of εb with respect to εa projected in εb, then `ṙ_ab ==
dt(r_ab, ω_ab_b)`
"""
dt(r_ab::RQuat, ω_ab_b::AbstractVector{<:Real}) = 0.5 * (r_ab._u * FreeQuat(imag=ω_ab_b))

"""
    ω(r_ab::RQuat, ṙ_ab::AbstractVector{<:Real})

Angular velocity from rotation quaternion and its time derivative

If `r_ab` is a `RQuat` representing the rotation from axes εa to axes εb, and
`̇r_ab` is its time derivative, `ω_ab_b = ω(r_ab, ̇r_ab)` is the angular
velocity of εb with respect to εa projected in εb axes
"""
ω(r_ab::RQuat, ṙ_ab::AbstractVector{<:Real}) = 2 * (UnitQuat(r_ab)' * FreeQuat(ṙ_ab)).imag

##### RQuat fallbacks #####

#require each Abstract3DRotation subtype to implement conversions to and from RQuat
Base.convert(::Type{RQuat}, r::R) where {R<:Abstract3DRotation} = error("Implement $R to RQuat conversion")
Base.convert(::Type{R}, r::RQuat) where {R<:Abstract3DRotation} = error("Implement RQuat to $R conversion")

#this enables a two-step conversion via RQuat, but requires defining trivial
#conversions to avoid infinite recursion
Base.convert(::Type{R}, r::Abstract3DRotation) where {R<:Abstract3DRotation} = convert(R, RQuat(r))
Base.convert(::Type{R}, r::R) where {R<:Abstract3DRotation} = r
Base.convert(::Type{RQuat}, r::RQuat) = r

#unless the representation defines its own methods, fall back to RQuat;
#promote rules don't work here, because promoting two inputs of the same subtype
#leaves them unchanged, and here we need them converted to RQuat
Base.adjoint(r::R) where {R<:Abstract3DRotation} = convert(R, RQuat(r)')
Base.:*(r::Abstract3DRotation, v::AbstractVector{<:Real}) = RQuat(r) * v
Base.:(≈)(r1::Abstract3DRotation, r2::Abstract3DRotation; kwargs...) = ≈(RQuat(r1), RQuat(r2); kwargs...)
Base.:∘(r1::Abstract3DRotation, r2::Abstract3DRotation) = RQuat(r1) ∘ RQuat(r2)

#absolute equality is not defined in general, because comparing two different
#subtypes requires promotion to RQuat, whose outcome is affected by floating
#point precision

######################### RMatrix #######################

"Rotation matrix representation"
struct RMatrix <: Abstract3DRotation
    _mat::SMatrix{3, 3, Float64, 9}
    function RMatrix(input::AbstractArray{<:Real, 2}; normalization::Bool = true)
        sm = SMatrix{3,3}(input)
        return normalization ? new(qr(sm).Q) : new(sm)
    end
end

#normalize by default
RMatrix() = RMatrix(SMatrix{3,3,Float64}(I), normalization = false)
RMatrix(r::Abstract3DRotation) = convert(RMatrix, r)

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
Base.:*(r::RMatrix, m::AbstractMatrix{<:Real}) = r._mat * m
Base.:*(m::AbstractMatrix{<:Real}, r::RMatrix) = m * r._mat

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
    Ω_ab_b = v2skew(ω_ab_b)
    return R_ab._mat * Ω_ab_b
end


############################# RAxAng ###############################

"Axis-angle representation"
struct RAxAng <: Abstract3DRotation
    axis::SVector{3, Float64}
    angle::Float64
    function RAxAng(axis::AbstractVector{<:Real}, angle::Real; normalization::Bool = true)
        axis = SVector{3,Float64}(axis) #normalization will be faster for an SVector
        return normalization ? new(normalize(axis), angle) : new(axis, angle)
    end
end

RAxAng(r::Abstract3DRotation) = convert(RAxAng, r)
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

# Inversion
Base.adjoint(r::RAxAng) = RAxAng(r.axis, -r.angle)

############################# RVec ###############################

"Rotation vector representation"
@kwdef struct RVec <: Abstract3DRotation
    _data::SVector{3, Float64} = @SVector[0.0, 0.0, 0.0]
end

RVec(r::Abstract3DRotation) = convert(RVec, r)

Base.getindex(ρ::RVec, i) = getindex(ρ._data, i)
Base.length(ρ::RVec) = length(ρ._data)
Base.size(ρ::RVec) = size(ρ._data)

Base.:+(ρ1::RVec, ρ2::RVec) = RVec(ρ1._data + ρ2._data)
Base.:-(ρ1::RVec, ρ2::RVec) = RVec(ρ1._data - ρ2._data)
Base.:*(a::Real, ρ::RVec) = RVec(a * ρ._data)
Base.:*(ρ::RVec, a::Real) = a * ρ

#conversions to and from RQuat
function Base.convert(::Type{RVec}, r::RQuat)
    q_re = r._u.real
    q_im = r._u.imag
    norm_im = norm(q_im)
    if norm_im > 0
        μ = 2atan(norm_im, q_re)
        u = q_im / norm_im
        ρ = RVec(μ * u)
    else
        ρ = RVec()
    end
    return ρ
end

LinearAlgebra.norm(ρ::RVec) = norm(ρ._data)

function Base.convert(::Type{RQuat}, ρ::RVec)
    μ = norm(ρ._data)
    if μ > 0
        u = ρ._data/μ
        r = RQuat(UnitQuat(real = cos(0.5μ), imag = u*sin(0.5μ), normalization = false))
    else
        r = RQuat()
    end
    return r
end

# Inversion
Base.adjoint(ρ::RVec) = RVec(-ρ._data)

############################# REuler #############################

"""
Euler angle representation. Rotation order is ZYX, so angles follow the ordering
ψ (heading), θ (inclination), φ (bank).
"""
struct REuler <: Abstract3DRotation
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

REuler(r::Abstract3DRotation) = convert(REuler, r)
REuler(input::Tuple{Real, Real, Real}) = REuler(input...)
REuler(; ψ = 0, θ = 0, φ = 0) = REuler(ψ, θ, φ)
function REuler(v::AbstractVector{<:Real})
    @assert length(v) == 3
    REuler(v[1], v[2], v[3])
end

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

function Base.convert(::Type{RMatrix}, r::REuler)

    cψ = cos(r.ψ); sψ = sin(r.ψ)
    cθ = cos(r.θ); sθ = sin(r.θ)
    cφ = cos(r.φ); sφ = sin(r.φ)

    M = @SMatrix [
    cψ * cθ         -sψ * cφ + cψ * sθ * sφ      sψ * sφ + cψ * sθ * cφ;
    sψ * cθ          cψ * cφ + sψ * sθ * sφ     -cψ * sφ + sψ * sθ * cφ;
    -sθ              cθ * sφ                     cθ * cφ
    ]

    RMatrix(M, normalization = false)

end

function Base.convert(::Type{REuler}, r::RMatrix)

    R = r._mat

    ψ = atan(R[2,1], R[1,1])
    θ = -asin( clamp(R[3,1], -1, 1) )
    φ = atan(R[3,2], R[3,3])

    REuler(ψ, θ, φ)

end


"""
    dt(e_ab::REuler, ω_ab_b::AbstractVector{<:Real})

Time derivative of Euler angles

If `e_ab` represents the rotation from axes εa to axes εb, and `ω_ab_b` is the
angular velocity of εb with respect to εa projected in εb, then `̇e_ab ==
dt(e_ab, ω_ab_b) is an array containing the time derivative of each angle in
e_ab`
"""
function dt(e_ab::REuler, ω_ab_b::AbstractVector{<:Real})

    sin_φ = sin(e_ab.φ); cos_φ = cos(e_ab.φ);
    tan_θ = tan(e_ab.θ); sec_θ = sec(e_ab.θ)

    M = @SMatrix [
                    0   sin_φ * sec_θ   cos_φ * sec_θ;
                    0   cos_φ           -sin_φ;
                    1   sin_φ * tan_θ   cos_φ * tan_θ
    ]

    ė_ab = M * SVector{3}(ω_ab_b)
    return ComponentVector(ė_ab, Axis(:ψ, :θ, :φ))
end

"""
    ω(e_ab::REuler, ė_ab::AbstractVector{<:Real})

Angular velocity vector from a set of Euler angles and its time derivative

If `e_ab` is a set of Euler angles representing the rotation from axes εa to
axes εb, and `̇e_ab` is its time derivative, `ω_ab_b = ω(e_ab, ̇e_ab)` is the
angular velocity of εb with respect to εa projected in εb axes
"""
function ω(e_ab::REuler, ė_ab::AbstractVector{<:Real})

    sin_θ = sin(e_ab.θ); cos_θ = cos(e_ab.θ);
    sin_φ = sin(e_ab.φ); cos_φ = cos(e_ab.φ);

    M = @SMatrix [
                  -sin_θ           0       1;
                  cos_θ * sin_φ    cos_φ   0;
                  cos_θ * cos_φ    -sin_φ  0;
    ]

    ω_ab_b = M * SVector{3}(ė_ab)
    return ω_ab_b

end

azimuth(v::AbstractVector{<:Real}) = atan(v[2], v[1])
inclination(v::AbstractVector{<:Real}) = atan(-v[3], √(v[1]^2 + v[2]^2))
wrap_to_π(x) = x + 2π*floor((π-x)/(2π))


################################# Plotting #####################################

#if no specific method available, convert to REuler for plotting
@recipe function f(ts::TimeSeries{<:Abstract3DRotation}; rot_ref = "", rot_target = "")

    return TimeSeries(ts._t, [REuler(v) for v in ts._data])

end

@recipe function f(ts::TimeSeries{<:REuler}; rot_ref = "", rot_target = "")

    label --> ["Heading" "Inclination" "Bank"]
    yguide --> hcat(L"$\psi_{%$rot_ref %$rot_target} \ (deg)$",
                    L"$\theta_{%$rot_ref %$rot_target} \ (deg)$",
                    L"$\phi_{%$rot_ref %$rot_target} \ (deg)$")
    ts_split --> :v #custom TimeSeries attribute
    link --> :none

    data = rad2deg.(hcat(ts.ψ._data, ts.θ._data, ts.φ._data)') #
    return TimeSeries(ts._t, data)

end


end