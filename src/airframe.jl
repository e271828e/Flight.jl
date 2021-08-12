module Airframe

using StaticArrays: SVector, SMatrix
using LinearAlgebra
using UnPack

using Flight.WGS84
using Flight.Attitude
using Flight.Kinematics

export Wrench, MassData, ComponentFrame, ComponentGroup
export v2skew, inertia_wrench, gravity_wrench, f_vel!

Base.@kwdef struct MassData
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end

Base.@kwdef struct Wrench
    F::SVector{3,Float64} = zeros(SVector{3})
    M::SVector{3,Float64} = zeros(SVector{3})
end
Base.show(io::IO, wr::Wrench) = print(io, "Wrench(F = $(wr.F), M = $(wr.M))")
Base.:+(wr1::Wrench, wr2::Wrench) = Wrench(F = wr1.F + wr2.F, M = wr1.M + wr2.M)


"""
#Defines a local component frame Fc(Oc, Ɛc) related to the airframe reference
frame Fb(Ob, Ɛb) by:
#a) the position vector of the local frame origin Oc relative to the reference
#frame origin Ob, projected in the reference frame axes
# b) the attitude of the local frame axes relative to the reference
#frame axes, given by rotation b_c
"""
Base.@kwdef struct ComponentFrame
    r_ObOc_b::SVector{3,Float64} = zeros(SVector{3})
    q_bc::RQuat = RQuat()
end

"""
translates a wrench specified on a local frame f2(O2, ε2) to a
reference frame f1(O1, ε1) given the frame transform from 1 to 2
"""
function Base.:*(f_bc::ComponentFrame, wr_Oc_c::Wrench)

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


function inertia_wrench(mass::MassData, vel::VelY, h_ext_b::AbstractVector{<:Real})

    @unpack m, J_Ob_b, r_ObG_b = mass
    @unpack ω_ie_b, ω_eb_b, ω_ib_b, v_eOb_b = vel #these are already SVectors

    h_ext_b = SVector{3,Float64}(h_ext_b)

    #angular momentum of the overall airframe as a rigid body
    h_rbd_b = J_Ob_b * ω_ib_b

    #total angular momentum
    h_all_b = h_rbd_b + h_ext_b

    #exact
    a_1_b = (ω_eb_b + 2 * ω_ie_b) × v_eOb_b
    F_in_Ob_b = -m * (a_1_b + ω_ib_b × (ω_ib_b × r_ObG_b) + r_ObG_b × (ω_eb_b × ω_ie_b ))
    M_in_Ob_b = - ( J_Ob_b * (ω_ie_b × ω_eb_b) + ω_ib_b × h_all_b + m * r_ObG_b × a_1_b)

    Wrench(F = F_in_Ob_b, M = M_in_Ob_b)

end

function gravity_wrench(mass::MassData, pos::PosY)

    #strictly, the gravity vector should be evaluated at G, with its direction
    #given by the z-axis of LTF(G). however, since g(G) ≈ g(Ob) and LTF(G) ≈
    #LTF(Ob), we can instead evaluate g at Ob, assuming its direction given by
    #LTF(Ob), and then apply it at G.
    g_G_l = gravity(pos.Ob)

    #the resultant consists of the force of gravity acting on G along the local
    #vertical and a null torque
    F_G_l = mass.m * g_G_l
    M_G_l = zeros(SVector{3})
    wr_G_l = Wrench(F = F_G_l, M = M_G_l)

    #with the previous assumption, the transformation from body frame to local
    #gravity frame is given by the translation r_ObG_b and the (passive)
    #rotation from b to LTF(Ob) (instead of LTF(G)), which is given by pos.l_b'
    wr_Oc_c = wr_G_l
    f_bc = ComponentFrame(r_ObOc_b = mass.r_ObG_b, q_bc = pos.q_lb')
    return f_bc * wr_Oc_c #wr_Ob_b

end

function f_vel!(ẋ_vel, wr_ext_Ob_b::Wrench, h_ext_b::AbstractVector{<:Real},
                mass::MassData, pv::PosVelY)

    #wr_ext_Ob_b: External Wrench on the airframe due to aircraft components

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
    @unpack ω_eb_b, ω_ie_b, v_eOb_b = pv.vel
    @unpack Ob, q_eb = pv.pos

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
    b = SVector{6}([wr_Ob_b.M ; wr_Ob_b.F])

    ẋ_vel .= A\b

    α_eb_b = SVector{3}(ẋ_vel.ω_eb_b) #α_eb_b == ω_eb_b_dot
    α_ib_b = α_eb_b - ω_eb_b × ω_ie_b

    v̇_eOb_b = SVector{3}(ẋ_vel.v_eOb_b)
    r_eO_e = rECEF(Ob)
    r_eO_b = q_eb' * r_eO_e
    a_eOb_b = v̇_eOb_b + ω_eb_b × v_eOb_b
    a_iOb_b = v̇_eOb_b + (ω_eb_b + 2ω_ie_b) × v_eOb_b + ω_ie_b × (ω_ie_b × r_eO_b)

    AccelY(α_eb_b, α_ib_b, a_eOb_b, a_iOb_b)

end


"""
Computes the skew-symmetric matrix corresponding to 3-element vector v.
"""
function v2skew(v::AbstractVector{T}) where {T<:Real}
    #much slower, each indexing operation yields an allocation
    # [0. -v[3] v[2]; v[3] 0. -v[1]; -v[2] v[1] 0.]
    M = zeros(T, 3, 3)
                    M[1,2] = -v[3];  M[1,3] = v[2]
    M[2,1] = v[3];                   M[2,3] = -v[1]
    M[3,1] = -v[2]; M[3,2] = v[1]

    SMatrix{3,3}(M)
end


end #module