module Dynamics

using StaticArrays
using LinearAlgebra
using UnPack
using ComponentArrays

using Flight.WGS84
using Flight.Rotations
using Flight.Kinematics

export AccY
export MassData, Wrench, FrameSpec
export get_wr_Ob_b, get_h_Gc_b, inertia_wrench, gravity_wrench, translate, f_dyn!


############################# FrameSpec ###############################

"""
Specifies a reference frame `fc(Oc, Ɛc)` relative to another `fb(Ob, Ɛb)`

Frame `fc(Oc, Ɛc)` is specified by:
- The position vector from fb's origin Ob to fc's origin Oc, projected on fb's
  axes εb (`r_ObOc_b`)
- The rotation quaternion from fb's axes εb to fc's axes εc (`q_bc`)
"""
Base.@kwdef struct FrameSpec
    r_ObOc_b::SVector{3,Float64} = zeros(SVector{3})
    q_bc::RQuat = RQuat()
end


############################## Wrench #################################

"""
Force and torque combination defined on a concrete reference frame

A `Wrench` is defined on reference frame fc(Oc,εc) when it is applied at its
origin Oc and projected on its axes εc
"""
Base.@kwdef struct Wrench
    F::SVector{3,Float64} = zeros(SVector{3})
    M::SVector{3,Float64} = zeros(SVector{3})
 end


 """
    Base.:+(wr1::Wrench, wr2::Wrench)

Add two compatible `Wrench` instances.

`Wrench` addition should only be performed between compatible `Wrench`
instances, i.e., those defined in the same reference frame
"""
Base.:+(wr1::Wrench, wr2::Wrench) = Wrench(F = wr1.F + wr2.F, M = wr1.M + wr2.M)


"""
    translate(f_bc::FrameSpec, wr_Oc_c::Wrench)

Translate a Wrench from one reference frame to another.

If `f_bc` is a `FrameSpec` specifying frame fc(Oc, εc) relative to fb(Ob, εb),
and `wr_Oc_c` is a `Wrench` defined on fc, then `wr_Ob_b = translate(f_bc,
wr_Oc_c)` is the equivalent `Wrench` defined on fb.

An alternative function call notation is available for this method:
`f_bc(wr_Oc_c) == translate(f_bc, wr_Oc_c)`
"""
function translate(f_bc::FrameSpec, wr_Oc_c::Wrench)

    @unpack q_bc, r_ObOc_b = f_bc
    F_Oc_c = wr_Oc_c.F
    M_Oc_c = wr_Oc_c.M

    #project onto airframe axes
    F_Oc_b = q_bc(F_Oc_c)
    M_Oc_b = q_bc(M_Oc_c)

    #translate to airframe origin
    F_Ob_b = F_Oc_b
    M_Ob_b = M_Oc_b + r_ObOc_b × F_Oc_b

    return Wrench(F = F_Ob_b, M = M_Ob_b) #wr_Ob_b

end

(f_bc::FrameSpec)(wr_Oc_c::Wrench) = translate(f_bc, wr_Oc_c)


############################## MassData #################################

"""
Groups the mass properties required by the aircraft's dynamic equations

These are:
- `m`: Total aircraft mass, including that of the rotating elements.
- `J_Ob_b`: Overall aircraft inertia tensor with respect to the airframe origin
  Ob, projected on the airframe axes εb.
- `r_ObG_b`: Position vector from the airframe origin Ob to the aircraft's
  center of mass G, projected on the airframe axes εb.

Some considerations regarding the role of mass properties on the aircraft
dynamic equations:
- All rotating elements on the aircraft are assumed to have axial symmetry.
  Under this assumption, it can be shown that the contribution of each rotating
  element to `J_Ob_b` is due to the lever arm from the airframe origin Ob to the
  center of mass of that rotating element (see tensor version of Steiner's
  theorem). These contributions must be included in `J_Ob_b`.
- Aircraft dynamics and kinematics are formulated on the airframe origin Ob
  instead of the aircraft's center of mass G. This allows for any of these mass
  properties to change, either gradually (for example, due to fuel consumption)
  or suddenly (due to a payload release), without causing discontinuities in the
  kinematic state vector.
"""
Base.@kwdef struct MassData
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end


##################
#every airframe component must output a struct that implements these methods for
#retrieving wr_Ob_b and h_Gc_b
function get_wr_Ob_b(::T, args...) where {T}
    error("Method get_wr_Ob_b not implemented for type $T or incorrect call signature")
end

function get_h_Gc_b(::T, args...) where {T}
    error("Method get_h_Gc_b not implemented for type $T or incorrect call signature")
end

#these automate wr_Ob_b and h_Gc_b retrieval for airframe component groups,
#whose output structs are gathered in a NamedTuple (for example, HybridSystems
#built from SystemDescriptorGroups)
@inline @generated function get_wr_Ob_b(y::NamedTuple{L}) where {L}

    ex = Expr(:block)
    push!(ex.args, :(wr = Wrench())) #allocate a zero wrench

    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            wr += get_wr_Ob_b(y[$label])
        end
        push!(ex.args, ex_ss)
    end
    return ex
end

@inline @generated function get_h_Gc_b(y::NamedTuple{L}) where {L}

    ex = Expr(:block)
    push!(ex.args, :(h = SVector(0., 0., 0.))) #allocate

    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            h += get_h_Gc_b(y[$label])
        end
        push!(ex.args, ex_ss)
    end
    return ex
end


################## Dynamic Equations and helper functions ####################

"""
    inertia_wrench(mass::MassData, y_vel::VelY, h_rot_b::AbstractVector{<:Real})

Compute the equivalent `Wrench` arising from inertia terms in the dynamic
equations

The resulting `Wrench` is defined on the airframe's reference frame.

# Arguments:
- `mass::MassData`: Current aircraft mass properties
- `y_vel::VelY`: Velocity outputs
- `h_rot_b::AbstractVector{<:Real}`: Additional angular momentum due to the
  angular velocity of any rotating elements with respect to the airframe,
  projected on the airframe axes

"""
function inertia_wrench(mass::MassData, y_vel::VelY, h_rot_b::AbstractVector{<:Real})

    @unpack m, J_Ob_b, r_ObG_b = mass
    @unpack ω_ie_b, ω_eb_b, ω_ib_b, v_eOb_b = y_vel

    #angular momentum of the overall airframe as a rigid body
    h_rbd_b = J_Ob_b * ω_ib_b

    #total angular momentum
    h_all_b = h_rbd_b + SVector{3,Float64}(h_rot_b)

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
    Ob = y_pos.Ob

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
    f_bc = FrameSpec(r_ObOc_b = mass.r_ObG_b, q_bc = y_pos.q_lb')
    return f_bc(wr_Oc_c) #wr_Ob_b

end


###################### Acceleration Outputs #####################

Base.@kwdef struct AccY
    α_eb_b::SVector{3,Float64}
    α_ib_b::SVector{3,Float64}
    a_eOb_b::SVector{3,Float64}
    a_iOb_b::SVector{3,Float64}
end

function f_dyn!(ẋ_vel::VelX, wr_ext_Ob_b::Wrench, h_rot_b::AbstractVector{<:Real},
    mass::MassData, y_kin::KinY)

    #wr_ext_Ob_b: Total external wrench on the airframe

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
    @unpack q_el, q_eb, Ob = y_kin.pos
    @unpack ω_eb_b, ω_ie_b, v_eOb_b = y_kin.vel

    r_ObG_b_sk = Rotations.skew(r_ObG_b)
    A11 = J_Ob_b
    A12 = m * r_ObG_b_sk
    A21 = -m * r_ObG_b_sk
    A22 = m * SMatrix{3,3,Float64}(I)

    A = vcat(hcat(A11, A12), hcat(A21, A22))

    wr_g_Ob_b = gravity_wrench(mass, y_kin.pos)
    wr_in_Ob_b = inertia_wrench(mass, y_kin.vel, h_rot_b)
    wr_Ob_b = wr_ext_Ob_b + wr_g_Ob_b + wr_in_Ob_b
    b = SVector{6}(vcat(wr_Ob_b.M, wr_Ob_b.F))

    ẋ_vel .= A\b

    # update ẋ_vel
    v̇_eOb_b = SVector{3}(ẋ_vel.v_eOb_b)
    r_eO_e = rECEF(Ob)
    r_eO_b = q_eb'(r_eO_e)

    α_eb_b = SVector{3}(ẋ_vel.ω_eb_b) #α_eb_b == ω_eb_b_dot
    α_ib_b = α_eb_b - ω_eb_b × ω_ie_b
    a_eOb_b = v̇_eOb_b + ω_eb_b × v_eOb_b
    a_iOb_b = v̇_eOb_b + (ω_eb_b + 2ω_ie_b) × v_eOb_b + ω_ie_b × (ω_ie_b × r_eO_b)

    return AccY(α_eb_b, α_ib_b, a_eOb_b, a_iOb_b)

end

end #module