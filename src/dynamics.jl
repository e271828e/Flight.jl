module Dynamics

using StaticArrays
using LinearAlgebra
using UnPack
using ComponentArrays

using Flight.WGS84
using Flight.Attitude
using Flight.Kinematics
import Flight.System: Y

export Acc, AccY
export MassData, WrenchCV, Wrench, Frame
export inertia_wrench, gravity_wrench, f_dyn!


################ Acceleration Output Vector ##############

struct Acc end

const AccYTemplate = ComponentVector(
    α_eb_b = zeros(3), α_ib_b = zeros(3), a_eOb_b = zeros(3), a_iOb_b = zeros(3))

const AccY{T, D} = ComponentVector{T, D, typeof(getaxes(AccYTemplate))} where {T, D}
Y(::Acc) = similar(AccYTemplate)


################ MassData ########################

Base.@kwdef struct MassData
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end


################# Wrench ########################

Base.@kwdef struct Wrench
    F::SVector{3,Float64} = zeros(SVector{3})
    M::SVector{3,Float64} = zeros(SVector{3})
 end

####################### Frame ###############

"""
#Specifies a local Frame fc(Oc, Ɛc) relative to the airframe reference
frame fb(Ob, Ɛb) by:
#a) the position vector of the local frame origin Oc relative to the reference
#frame origin Ob, projected in the reference frame axes
# b) the attitude of the local frame axes relative to the reference
#frame axes, given by rotation quaternion q_bc
"""
Base.@kwdef struct Frame
    r_ObOc_b::SVector{3,Float64} = zeros(SVector{3})
    q_bc::RQuat = RQuat()
end

"""
Translate a Wrench specified on a local Frame fc(Oc, εc) to the
airframe reference frame fb(Ob, εb) given the relative Frame
specification f_bc
"""

function Base.:*(f_bc::Frame, wr_Oc_c::Wrench)

    F_Oc_c = wr_Oc_c.F
    M_Oc_c = wr_Oc_c.M

    #project on the reference axes
    F_Oc_b = f_bc.q_bc * F_Oc_c
    M_Oc_b = f_bc.q_bc * M_Oc_c

    #translate them to airframe origin
    F_Ob_b = F_Oc_b
    M_Ob_b = M_Oc_b + f_bc.r_ObOc_b × F_Oc_b
    Wrench(F = F_Ob_b, M = M_Ob_b) #wr_Ob_b

end

################## f_dyn! and helper functions ####################

function inertia_wrench(mass::MassData, y_vel::VelY, h_rot_b::AbstractVector{<:Real})

    @unpack m, J_Ob_b, r_ObG_b = mass

    ω_ie_b = SVector{3,Float64}(y_vel.ω_ie_b)
    ω_eb_b = SVector{3,Float64}(y_vel.ω_eb_b)
    ω_ib_b = SVector{3,Float64}(y_vel.ω_ib_b)
    v_eOb_b = SVector{3,Float64}(y_vel.v_eOb_b)

    #additional angular momentum due to the angular velocity of the rotating
    #elements with respect to the airframe
    h_rot_b = SVector{3,Float64}(h_rot_b)

    #angular momentum of the overall airframe as a rigid body
    h_rbd_b = J_Ob_b * ω_ib_b

    #total angular momentum
    h_all_b = h_rbd_b + h_rot_b

    #exact
    a_1_b = (ω_eb_b + 2 * ω_ie_b) × v_eOb_b
    F_in_Ob_b = -m * (a_1_b + ω_ib_b × (ω_ib_b × r_ObG_b) + r_ObG_b × (ω_eb_b × ω_ie_b ))
    M_in_Ob_b = - ( J_Ob_b * (ω_ie_b × ω_eb_b) + ω_ib_b × h_all_b + m * r_ObG_b × a_1_b)

    Wrench(F = F_in_Ob_b, M = M_in_Ob_b)

end

function gravity_wrench(mass::MassData, y_pos::PosY)

    #strictly, the gravity vector should be evaluated at G, with its direction
    #given by the z-axis of LTF(G). however, since g(G) ≈ g(Ob) and LTF(G) ≈
    #LTF(Ob), we can instead evaluate g at Ob, assuming its direction given by
    #LTF(Ob), and then apply it at G.
    q_el = RQuat(y_pos.q_el, normalization = false)
    q_lb = RQuat(y_pos.q_lb, normalization = false)
    Ob = WGS84Pos(NVector(q_el), y_pos.h)

    g_G_l = gravity(Ob)

    #the resultant consists of the force of gravity acting on G along the local
    #vertical and a null torque
    F_G_l = mass.m * g_G_l
    M_G_l = zeros(SVector{3})
    wr_G_l = Wrench(F = F_G_l, M = M_G_l)

    #with the previous assumption, the transformation from body frame to local
    #gravity frame is given by the translation r_ObG_b and the (passive)
    #rotation from b to LTF(Ob) (instead of LTF(G)), which is given by pos.l_b'
    wr_Oc_c = wr_G_l
    f_bc = Frame(r_ObOc_b = mass.r_ObG_b, q_bc = q_lb')
    return f_bc * wr_Oc_c #wr_Ob_b

end

function f_dyn!(y_acc::AccY, ẋ_vel::VelX, wr_ext_Ob_b::Wrench,
    h_rot_b::AbstractVector{<:Real}, mass::MassData, y_kin::KinY)

    #wr_ext_Ob_b: External wrench on the airframe due to aircraft components

    #h_rot_b: Additional angular momentum due to the angular velocity of the
    #rotating aircraft components with respect to the airframe. these are
    #computed individually by each component relative to its center of mass and
    #then summed

    #wr_ext_Ob_b and h_rot_b, as well as mass data, are produced by aircraft
    #components, so they must be computed by the aircraft's x_dot method. and,
    #since y_kin is needed by those components, it must be called from the
    #aircraft's kinematic state vector

    #clearly, r_ObG_b cannot be arbitrarily large, because J_Ob_b is larger than
    #J_G_b (Steiner). therefore, at some point J_G_b would become zero (or at
    #least singular)!

    @unpack m, J_Ob_b, r_ObG_b = mass

    ω_eb_b = SVector{3,Float64}(y_kin.vel.ω_eb_b)
    ω_ie_b = SVector{3,Float64}(y_kin.vel.ω_ie_b)
    v_eOb_b = SVector{3,Float64}(y_kin.vel.v_eOb_b)
    q_el = RQuat(y_kin.pos.q_el, normalization = false)
    q_eb = RQuat(y_kin.pos.q_eb, normalization = false)
    Ob = WGS84Pos(NVector(q_el), y_kin.pos.h)

    #preallocating is faster than directly concatenating the blocks
    A = Array{Float64}(undef, (6,6))

    r_ObG_b_sk = Attitude.skew(r_ObG_b)
    A[1:3, 1:3] .= J_Ob_b
    A[1:3, 4:6] .= m * r_ObG_b_sk
    A[4:6, 1:3] .= -m * r_ObG_b_sk
    A[4:6, 4:6] .= m * SMatrix{3,3,Float64}(I)

    A = SMatrix{6,6}(A)

    wr_g_Ob_b = gravity_wrench(mass, y_kin.pos)
    wr_in_Ob_b = inertia_wrench(mass, y_kin.vel, h_rot_b)
    wr_Ob_b = wr_ext_Ob_b + wr_g_Ob_b + wr_in_Ob_b
    b = SVector{6}([wr_Ob_b.M ; wr_Ob_b.F])

    ẋ_vel .= A\b

    #update y_acc
    v̇_eOb_b = SVector{3}(ẋ_vel.v_eOb_b)
    r_eO_e = rECEF(Ob)
    r_eO_b = q_eb' * r_eO_e

    α_eb_b = SVector{3}(ẋ_vel.ω_eb_b) #α_eb_b == ω_eb_b_dot
    α_ib_b = α_eb_b - ω_eb_b × ω_ie_b
    a_eOb_b = v̇_eOb_b + ω_eb_b × v_eOb_b
    a_iOb_b = v̇_eOb_b + (ω_eb_b + 2ω_ie_b) × v_eOb_b + ω_ie_b × (ω_ie_b × r_eO_b)

    y_acc.α_eb_b .= α_eb_b #α_eb_b == ω_eb_b_dot
    y_acc.α_ib_b .= α_ib_b
    y_acc.a_eOb_b .= a_eOb_b
    y_acc.a_iOb_b .= a_iOb_b

    return nothing

end

end #module