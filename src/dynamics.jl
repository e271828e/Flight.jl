module Dynamics

using StaticArrays
using LinearAlgebra
using UnPack

using Flight.WGS84
using Flight.Attitude
using Flight.Airframe
using Flight.Kinematics
using Flight.System

export AccY
export MassData
export f_dyn!


############################## MassData #################################

"""
Groups the mass properties required by the aircraft's dynamic equations

These are:
- `m`: Total aircraft mass, including that of the rotating elements.
- `J_Ob_b`: Overall aircraft inertia tensor with respect to the airframe origin
  Ob, projected on the airframe axes εb.
- `r_ObG_b`: Position vector from the airframe origin Ob to the aircraft's
  center of mass G, projected on the airframe axes εb.

Considerations on the role of mass properties in the aircraft dynamic equations:
- The inertia tensor `J_Ob_b` must include the contributions of any rotating
  elements on the aircraft. As long as rotating elements have axial symmetry
  around their axes of rotation, their contributions to `J_Ob_b` will be
  strictly constant.
- Aircraft dynamics and kinematics are formulated on the airframe origin Ob
  instead of the aircraft's center of mass G. This allows for any of the
  aircraft's mass properties to change, either gradually (for example, due to
  fuel consumption) or suddenly (due to a payload release), without causing
  discontinuities in the kinematic state vector.
"""
Base.@kwdef struct MassData
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end


################## Dynamic Equations and helper functions ####################

"""
    inertia_wrench(mass::MassData, y_vel::VelY, hr_b::AbstractVector{<:Real})

Compute the equivalent `Wrench` arising from inertia terms in the dynamic
equations

The resulting `Wrench` is defined on the airframe's reference frame.

# Arguments:
- `mass::MassData`: Current aircraft mass properties
- `y_vel::VelY`: Velocity outputs
- `hr_b::AbstractVector{<:Real}`: Additional angular momentum due to the
  angular velocity of any rotating elements with respect to the airframe,
  projected on the airframe axes

"""
function inertia_wrench(mass::MassData, y_vel::VelY, hr_b::AbstractVector{<:Real})

    @unpack m, J_Ob_b, r_ObG_b = mass
    @unpack ω_ie_b, ω_eb_b, ω_ib_b, v_eOb_b = y_vel

    #angular momentum of the overall airframe as a rigid body
    h_rbd_b = J_Ob_b * ω_ib_b

    #total angular momentum
    h_all_b = h_rbd_b + SVector{3,Float64}(hr_b)

    #inertia terms (exact... overkill, but very cheap anyway)
    a_1_b = (ω_eb_b + 2 * ω_ie_b) × v_eOb_b
    F_in_Ob_b = -m * (a_1_b + ω_ib_b × (ω_ib_b × r_ObG_b) + r_ObG_b × (ω_eb_b × ω_ie_b ))
    M_in_Ob_b = - ( J_Ob_b * (ω_ie_b × ω_eb_b) + ω_ib_b × h_all_b + m * r_ObG_b × a_1_b)

    return Wrench(F = F_in_Ob_b, M = M_in_Ob_b)

end

function gravity_wrench(mass::MassData, y_pos::PosY)

    #gravity can be viewed as an entity acting on a local frame with its origin
    #at G and its axes aligned with the local tangent frame

    #strictly, the gravity vector should be evaluated at G, with its direction
    #given by the z-axis of LTF(G). however, since g(G) ≈ g(Ob) and LTF(G) ≈
    #LTF(Ob), we can instead evaluate g at Ob, assuming its direction given by
    #LTF(Ob), and then apply it at G.
    g_G_l = g_Ob_l = g_l(y_pos.Ob_nvh)

    #the resultant consists of the gravity force acting on G along the local
    #vertical and a null torque
    F_G_l = mass.m * g_G_l
    M_G_l = zeros(SVector{3})
    wr_G_l = Wrench(F = F_G_l, M = M_G_l)

    #with the previous assumption, the transformation from body frame to local
    #gravity frame is given by the translation r_ObG_b and the (passive)
    #rotation from b to LTF(Ob) (instead of LTF(G)), which is given by pos.l_b'
    wr_c = wr_G_l
    f_bc = FrameSpec(r_ObOc_b = mass.r_ObG_b, q_bc = y_pos.q_lb')
    return f_bc(wr_c) #wr_b

end


###################### Acceleration Outputs #####################

struct Acc <: AbstractComponent end

Base.@kwdef struct AccY <: AbstractY{Acc}
    α_eb_b::SVector{3,Float64}
    α_ib_b::SVector{3,Float64}
    a_eOb_b::SVector{3,Float64}
    a_iOb_b::SVector{3,Float64}
    f_Ob_b::SVector{3,Float64} #specific force
end

function f_dyn!(ẋ_vel::VelX, wr_ext_b::Wrench, hr_b::AbstractVector{<:Real},
    mass::MassData, y_kin::KinY)

    #wr_ext_b: Total external wrench on the airframe

    #hr_b: Additional angular momentum due to the angular velocity of the
    #rotating aircraft components with respect to the airframe. these are
    #computed individually by each component relative to its center of mass and
    #then summed

    #wr_ext_b and hr_b, as well as mass data, are produced by aircraft
    #components, so they must be computed by the aircraft's x_dot method. and,
    #since y_kin is needed by those components, it must be called from the
    #aircraft's kinematic state vector

    #clearly, r_ObG_b cannot be arbitrarily large, because J_Ob_b is larger than
    #J_G_b (Steiner). therefore, at some point J_G_b would become zero (or at
    #least singular)!

    @unpack m, J_Ob_b, r_ObG_b = mass
    @unpack q_lb, q_el, q_eb, Ob_nvh = y_kin.pos
    @unpack ω_eb_b, ω_ie_b, v_eOb_b = y_kin.vel

    r_ObG_b_sk = Attitude.skew(r_ObG_b)
    A11 = J_Ob_b
    A12 = m * r_ObG_b_sk
    A21 = -m * r_ObG_b_sk
    A22 = m * SMatrix{3,3,Float64}(I)

    A = vcat(hcat(A11, A12), hcat(A21, A22))

    wr_g_b = gravity_wrench(mass, y_kin.pos)
    wr_in_b = inertia_wrench(mass, y_kin.vel, hr_b)
    wr_b = wr_ext_b + wr_g_b + wr_in_b
    b = SVector{6}(vcat(wr_b.M, wr_b.F))

    ẋ_vel .= A\b

    # update ẋ_vel
    v̇_eOb_b = SVector{3}(ẋ_vel.v_eOb_b)
    r_eOb_e = CartECEF(Ob_nvh).data
    r_eOb_b = q_eb'(r_eOb_e)

    #angular accelerations
    α_eb_b = SVector{3}(ẋ_vel.ω_eb_b) #α_eb_b == ω_eb_b_dot
    α_ib_b = α_eb_b - ω_eb_b × ω_ie_b

    #linear accelerations and specific force
    a_eOb_b = v̇_eOb_b + ω_eb_b × v_eOb_b
    a_iOb_b = v̇_eOb_b + (ω_eb_b + 2ω_ie_b) × v_eOb_b + ω_ie_b × (ω_ie_b × r_eOb_b)

    g_Ob_b = q_lb'(g_l(Ob_nvh))
    G_Ob_b = g_Ob_b + ω_ie_b × (ω_ie_b × r_eOb_b)
    f_Ob_b = a_iOb_b - G_Ob_b

    return AccY(α_eb_b, α_ib_b, a_eOb_b, a_iOb_b, f_Ob_b)

end

end #module