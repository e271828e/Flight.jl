module Dynamics

using StaticArrays: SVector, SMatrix
using LinearAlgebra
using UnPack

using Flight.LBV
using Flight.WGS84
using Flight.Attitude
using Flight.Kinematics

export XAirframe, Wrench, FrameTransform, MassData
export v2skew, inertia_wrench, gravity_wrench, x_vel_dot

"""
Computes the skew-symmetric matrix corresponding to 3-element vector v.
"""
# function v2skew(v::AbstractVector{T} where {T<:Real})
function v2skew(v::AbstractVector{T}) where {T<:Real}
    #much slower, each indexing operation yields an allocation
    # [0. -v[3] v[2]; v[3] 0. -v[1]; -v[2] v[1] 0.]
    M = zeros(T, 3, 3)
                    M[1,2] = -v[3];  M[1,3] = v[2]
    M[2,1] = v[3];                   M[2,3] = -v[1]
    M[3,1] = -v[2]; M[3,2] = v[1]

    SMatrix{3,3}(M)
end


Base.@kwdef struct Wrench
    T::SVector{3,Float64} = zeros(SVector{3})
    F::SVector{3,Float64} = zeros(SVector{3})
end
Base.show(io::IO, wr::Wrench) = print(io, "Wrench(T = $(wr.T), F = $(wr.F))")
Base.:+(wr1::Wrench, wr2::Wrench) = Wrench(T = wr1.T + wr2.T, F = wr1.F + wr2.F)

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
    f_bc = FrameTransform(r_ObOc_b = mass.r_ObG_b, q_bc = pos.q_lb')
    wr_Oc_c = wr_G_l
    wr_Ob_b = f_bc * wr_Oc_c

    return wr_Ob_b

end

#maybe define abstract type PVData
#then PVDataWGS84 <: PVData. then we can define different gravity_wrench and
#inertia_wrench methods for PVDataWGS84 and PVDataFlatEarth

#maybe consider splitting

function x_vel_dot(wr_ext_Ob_b::Wrench, h_ext_b::AbstractVector{T} where {T<:Real}, mass::MassData, pv::PVDataWGS84)

    #wr_ext_Ob_b: Wrench due to aircraft components

    #h_ext_b: Additional angular momentum due to rotating aircraft components
    #(computed using their angular velocity wrt the airframe, not the inertial
    #frame)

    #wr_ext_Ob_b and h_ext_b, as well as mass data, are produced by aircraft
    #components, so they must be computed by the aircraft's x_dot method. and,
    #since pv is needed by those components, it must be called from the
    #aircraft's kinematic state vector

    #clearly, r_ObG_b cannot be arbitrarily large, because J_Ob_b is larger than
    #J_G_b (Steiner). therefore, at some point J_G_b would become zero (or at
    #least singular)!

    @unpack m, J_Ob_b, r_ObG_b = mass

    #preallocating is faster than directly concatenating the blocks
    A = Array{Float64}(undef, (6,6))

    r_ObG_b_sk = v2skew(r_ObG_b)
    A[1:3, 1:3] .= J_Ob_b
    A[1:3, 4:6] .= m * r_ObG_b_sk
    A[4:6, 1:3] .= -m * r_ObG_b_sk
    A[4:6, 4:6] .= m * SMatrix{3,3,Float64}(I)

    A = SMatrix{6,6}(A)

    wr_g_Ob_b = gravity_wrench(mass, pv.pos)
    wr_in_Ob_b = inertia_wrench(mass, pv.vel, h_ext_b)
    wr_Ob_b = wr_ext_Ob_b + wr_g_Ob_b + wr_in_Ob_b
    b = [wr_Ob_b.T ; wr_Ob_b.F]

    XVel(A\b)

end


end #module