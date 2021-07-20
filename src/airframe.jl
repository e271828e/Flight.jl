module Airframe

# using Base: Float64
using StaticArrays: SVector, SMatrix
using LinearAlgebra
using UnPack

using Flight.LBV
using Flight.WGS84
using Flight.Attitude
using Flight.Kinematics

export XAirframe, Wrench, FrameTransform, MassData
export inertia_wrench, gravity_wrench

@define_node XAirframe (kin = XKinWGS84,)

Base.@kwdef struct Wrench
    T::SVector{3,Float64} = zeros(SVector{3})
    F::SVector{3,Float64} = zeros(SVector{3})
end
Base.show(io::IO, wr::Wrench) = print(io, "Wrench(T = $(wr.T), F = $(wr.F))")


#defines the transform f_bc from the airframe reference frame Fb(Ob, Ɛb) to a
#local component frame Fc(Oc, Ɛc) by:
#a) the position vector of the local frame origin Oc relative to the reference
#frame origin Ob, projected in the reference frame axes
# b) the attitude of the local frame axes relative to the reference
#frame axes, given by rotation b_c

Base.@kwdef struct FrameTransform
    r_ObOc_b::SVector{3,Float64} = zeros(SVector{3})
    q_bc::RQuat = RQuat()
end

function Base.:*(f_bc::FrameTransform, wr_Oc_c::Wrench)

    #translates a wrench specified on a local frame f2(O2, ε2) to a
    #reference frame f1(O1, ε1) given the frame transform from 1 to 2

    T_Oc_c = wr_Oc_c.T
    F_Oc_c = wr_Oc_c.F

    #project on the reference axes
    T_Oc_b = f_bc.q_bc * T_Oc_c
    F_Oc_b = f_bc.q_bc * F_Oc_c

    #translate them to airframe origin
    F_Ob_b = F_Oc_b
    T_Ob_b = T_Oc_b + f_bc.r_ObOc_b × F_Oc_b

    wr_Ob_b = Wrench(T = T_Ob_b, F = F_Ob_b)

    return wr_Ob_b

end

Base.@kwdef struct MassData
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end

function inertia_wrench(mass::MassData, vel::VelDataWGS84, h_add_b::AbstractVector{T} where {T<:Real})

    #h_add_b: angular momentum due to rotating airframe components (computed using
    #their angular velocity wrt the airframe, not the inertial frame)

    @unpack m, J_Ob_b, r_ObG_b = mass
    @unpack ω_ie_b, ω_eb_b, ω_ib_b, v_eOb_b = vel #these are already SVectors

    h_add_b = SVector{3,Float64}(h_add_b)

    #angular momentum of the overall airframe as a rigid body
    h_rbd_b = J_Ob_b * ω_ib_b

    #total angular momentum
    h_all_b = h_rbd_b + h_add_b

    #exact
    a_1_b = (ω_eb_b + 2 * ω_ie_b) × v_eOb_b
    T_in_Ob_b = - ( J_Ob_b * (ω_ie_b × ω_eb_b) + ω_ib_b × h_all_b + m * r_ObG_b × a_1_b)
    F_in_Ob_b = -m * (a_1_b + ω_ib_b × (ω_ib_b × r_ObG_b) + r_ObG_b × (ω_eb_b × ω_ie_b ))

    Wrench(T = T_in_Ob_b, F = F_in_Ob_b)

end

    # q_bl = q_be * q_el
function gravity_wrench(mass::MassData, pos::PosDataWGS84)

    #strictly, the gravity vector should be evaluated at G, with its direction
    #given by the z-axis of LTF(G). however, since g(G) ≈ g(Ob) and LTF(G) ≈
    #LTF(Ob), we can instead evaluate g at Ob, assuming its direction given by
    #LTF(Ob), and then apply it at G.
    g_G_l = gravity(pos.Ob)

    #the resultant consists of the force of gravity acting on G along the local
    #vertical and a null torque
    F_G_l = mass.m * g_G_l
    T_G_l = zeros(SVector{3})
    wr_G_l = Wrench(T = T_G_l, F = F_G_l)

    #with the previous assumption, the transformation from body frame to local
    #gravity frame is given by the translation r_ObG_b and the (passive)
    #rotation from b to LTF(Ob) (instead of LTF(G)), which is given by pos.l_b'
    f_bc = FrameTransform(r_ObOc_b = mass.r_ObG_b, q_bc = pos.l_b')
    wr_Oc_c = wr_G_l
    wr_Ob_b = f_bc * wr_Oc_c

    return wr_Ob_b

end


function x_vel_dot(wr_Ob_b::Wrench, mass::MassData, kin::PVDataWGS84)::XVel
    return XVel()
end

# function x_vel_dot(wr_Ob_b::Wrench, mass::MassData, kin::KinDataFlat)::XVel
#     return XVel()
# end



end #module