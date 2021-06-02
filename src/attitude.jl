module Attitude

# https://discourse.julialang.org/t/writing-functions-for-types-defined-in-another-module/31895
# https://discourse.julialang.org/t/what-is-the-preferred-way-to-use-multiple-files/30687/6

using StaticArrays: SVector
using LinearAlgebra
using ..Quaternions: UnitQuat, Quat

export Rotation, RQuat, AxAng, invert, compose, transform, dt

#########################

#implements a passive (alias) rotation, that is, a rotation describing the
#relative attitude between two reference frames.

abstract type Rotation end

Base.adjoint(r::Rotation) = invert(r)
Base.:∘(r1::Rotation, r2::Rotation) = compose(r1, r2)
Base.:∘(r::Rotation, x::Any) = error("$(typeof(r)) ∘ $(typeof(x)) composition not allowed")
Base.:∘(x::Any, r::Rotation) = error("$(typeof(x)) ∘ $(typeof(r)) composition not allowed")
Base.:*(r::Rotation, v::AbstractVector{T} where {T<:Real})  = transform(r, v)

############################# RQuat ###############################

mutable struct RQuat <: Rotation
    _quat::UnitQuat
    #if input is already a UnitQuat, new will not call convert, and therefore a
    #reference instead of a copy will be assigned to field quat. however, this
    #is no concern, since UnitQuats are virtually immutable (except for
    #normalize!), and there is garbage collection
end

RQuat() = RQuat(UnitQuat())

Base.:(==)(q1::RQuat, q2::RQuat) = (q1._quat == q2._quat || q1._quat == -q2._quat)
Base.:(≈)(q1::RQuat, q2::RQuat) = (q1._quat ≈ q2._quat || q1._quat ≈ -q2._quat)

LinearAlgebra.norm(q::RQuat) = norm(q._quat) #uses StaticArrays implementation
LinearAlgebra.normalize(q::RQuat) = RQuat(normalize(q._quat))
LinearAlgebra.normalize!(q::RQuat) = (q._quat = normalize(q._quat))

invert(q::RQuat) = RQuat(q._quat')

compose(q1::RQuat, q2::RQuat) = RQuat(q1._quat * q2._quat)

function transform(q_ab::RQuat, v_b::AbstractVector{T} where {T<:Real})
    q_ab = q_ab._quat; q_ab_imag = q_ab.imag
    v_a = v_b + 2 * cross(q_ab_imag, q_ab.real * v_b + cross(q_ab_imag, v_b))
    return v_a
end

dt(q_ab::RQuat, ω_ab_b::AbstractVector{T} where {T<:Real}) = 0.5 * (q_ab._quat * Quat(imag=ω_ab_b))


############################# RQuat ###############################

const RAxis = SVector{3, Float64}

struct AxAng <: Rotation

    axis::RAxis
    angle::Float64

    function AxAng(axis::AbstractVector{T} where {T<:Real}, angle::Real; enforce_norm::Bool = true)
        return enforce_norm ? new(normalize(axis), angle) : new(axis, angle)
    end

end

AxAng(input::Tuple{AbstractVector{S}, T} where {S, T <: Real}) = AxAng(input[1], input[2])
AxAng() = AxAng(RAxis([1, 0, 0]), 0, enforce_norm = false)

#MAKE SURE ENFORCE NORM IS ONLY CALLED WHEN NECESSARY!

# declare RMat, RVec, RQuat, AxAng, Euler types and work with them, providing
# convert() methods between them as required. in general, only to and from
# quaternion will be required.

#AxAng, RVec and Euler should all promote to Quat.

#conversions to provide:
# RMat, AxAng, RVec and Euler to and from RQuat
# RMat to and from Quat
# with that, can provide two-step conversion to and from RMat to AxAng, RVec and
# Euler. maybe Euler could be direct.


#RQuat is the center. every RotationType subtype must provide a:
#convert(::Type{RotationType}, input::RQuat)
#convert(::Type{RQuat}, input::RotationType)

#for example:
# convert(::Type{AxAng}, input::RQuat) = ...

#with this, any conversion is possible in two steps as:
# convert(::Type{AxAng}, input::Euler) = convert(AxAng, convert(RQuat, input))
#more generally:
# convert(::Type{T}, input::Rotation) where {T<:Rotation}= convert(T, convert(RQuat, input))
#in the specific case in which input is already a RQuat, the conversion output
#is itself, it will be bypassed

#in some cases, we can override the convert method and replace it with a direct
#conversion. for example, from axis-angle to rotation vector it is
#straightforward, makes little sense to go through quaternion. or maybe i could
#provide direct conversions to and from RMat and Euler

#put println statements in each conversion method to know when it's been called

# struct AxAng <: Rotation
#     axis::Float64
#     angle::SVector{3, Float64}
# end

#Rotation is an abstract type, of which these are concrete subtypes. not all
#operations are defined for every one of them.

#how do i make RVec(input) work for input<:Rotation? Simply make the constructor
#convert input to RQuat, then from RQuat to RVec
#convert(RQuat, input)
#convert(RVec, input)
#in fact, since we know that we will support conversions to Quat for all
#Rotation subtypes, we can simply do:
#RVec(input::Rotation) = RVec(RQ)
#Euler(1, 2, 3)

#a RMat can be constructed from a 3x3 Matrix. in that case, orthonormality is
#always enforced by applying QR factorization. if it comes from a conversion, it
#is not necessary

#create normalize! and normalize methods for RMat

#i can admit * of rotation matrices without converting to RQuat, and also
#transforming vectors! only use promote when absolutely required!
#define promote_rule(R, q) = q (always prefer q)

#Euler can be a struct with fields upsi, utheta, uphi

#replace argument types with StaticArrays to fix lengths. also, allow different
#Real subtypes

#go to Flight folder
#enter package manager
#activate .
#do not do activate Flight. in that case, the test and using commands do not
#work. they expect to be run from Flight's parent folder. why??

#must fail when passed a 1D array of size other than 3
#how do we do this? by defining StaticArrays types


end