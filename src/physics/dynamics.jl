module Dynamics

using StaticArrays
using LinearAlgebra
using UnPack
using Plots

using Flight.Utils
using Flight.Systems
using Flight.Plotting

using Flight.Attitude
using Flight.Geodesy
using Flight.Kinematics

import Flight.Plotting: make_plots

export FrameTransform, transform
export Wrench
export AbstractMassDistribution, PointMass, RigidBody, MassProperties
export HasMass, HasNoMass, get_mp_b
export GetsExternalWrench, GetsNoExternalWrench, get_wr_b
export HasAngularMomentum, HasNoAngularMomentum, get_hr_b
export DynData, f_dyn!



############################# FrameTransform ###############################

"""
Specifies a reference frame `fc(Oc, Ɛc)` relative to another `fb(Ob, Ɛb)`

Frame `fc(Oc, Ɛc)` is defined by:
- The position vector from fb's origin Ob to fc's origin Oc, projected on fb's
  axes εb (`r_ObOc_b`)
- The rotation quaternion from fb's axes εb to fc's axes εc (`q_bc`)
"""
Base.@kwdef struct FrameTransform
    r::SVector{3,Float64} = zeros(SVector{3})
    q::RQuat = RQuat()
end


"""
    Base.:∘(t_bc::FrameTransform, t_cd::FrameTransform)

Concatenate two `FrameTransform` instances.
"""
function Base.:∘(t_bc::FrameTransform, t_cd::FrameTransform)

    r_ObOc_b = t_bc.r; q_bc = t_bc.q
    r_OcOd_c = t_cd.r; q_cd = t_cd.q

    r_ObOd_b = r_ObOc_b + q_bc(r_OcOd_c)
    q_bd = q_bc ∘ q_cd

    t_bd = FrameTransform(r_ObOd_b, q_bd)

    return t_bd

end


"""
    Base.:adjoint(t_bc::FrameTransform)

Get the reciprocal `FrameTransform`.
"""
function Base.:adjoint(t_bc::FrameTransform)

    r_ObOc_b = t_bc.r; q_bc = t_bc.q

    q_cb = q_bc'
    r_OcOb_c = q_cb(-r_ObOc_b)
    t_cb = FrameTransform(r_OcOb_c, q_cb)

    return t_cb
end


"""
    transform(t_bc::FrameTransform, r_OcP_c::AbstractVector{<:Real})

Given the (3D) position vector of some point P in reference frame fc, and the
`FrameTransform` from fb to fc, compute the position vector of P in fb
"""
function transform(t_bc::FrameTransform, r_OcP_c::AbstractVector{<:Real})

    r_ObOc_b = t_bc.r; q_bc = t_bc.q

    r_ObP_b = r_ObOc_b + q_bc(SVector{3,Float64}(r_OcP_c))
    return r_ObP_b
end


"""Alternative function call notation: `t_bc(p_c) == transform(t_bc, p_c)`
"""
(t_bc::FrameTransform)(x) = transform(t_bc, x)



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
    transform(t_bc::FrameTransform, wr_c::Wrench)

Transform a `Wrench` from one reference frame to another.

If `t_bc` is a `FrameTransform` specifying frame fc(Oc, εc) relative to fb(Ob, εb),
and `wr_c` is a `Wrench` defined on fc, then `wr_b = transform(t_bc,
wr_c)` is the equivalent `Wrench` defined on fb.
"""

function transform(t_bc::FrameTransform, wr_c::Wrench)

    F_Oc_c = wr_c.F
    M_Oc_c = wr_c.M

    #project onto vehicle axes
    F_Oc_b = t_bc.q(F_Oc_c)
    M_Oc_b = t_bc.q(M_Oc_c)

    #transform to vehicle origin
    F_Ob_b = F_Oc_b
    M_Ob_b = M_Oc_b + t_bc.r × F_Oc_b

    return Wrench(F = F_Ob_b, M = M_Ob_b) #wr_b

end


############################## MassDistributions #################################

abstract type AbstractMassDistribution end

struct PointMass <: AbstractMassDistribution
    m::Float64
end

"""Defines a rigid body mass distribution.

A `RigidBody` is defined in a local reference frame fc, whose origin Oc is
located at the body's center of mass G, and whose axes εc are arbitrarily chosen
to express the inertia tensor conveniently (typically they will be the body's
principal axes of inertia)
- m: Total mass
- J: Inertia tensor computed with respect to Oc and projected in εc
"""
struct RigidBody <: AbstractMassDistribution
    m::Float64
    J::SMatrix{3,3,Float64,9}
end


"""
Groups the mass properties of a component expressed on a specific reference
frame fb(Ob,εb).

These are:
- `m`: Total mass, including that of rotating elements.
- `J_O`: Overall inertia tensor with respect to the origin Ob, projected on
  axes εb, more explicitly, J_Ob_b.
- `r_OG`: Position vector from the origin Ob to the component's center of
  mass G, projected on axes εb, more explicitly r_ObG_b.

Notes:
- The inertia tensor `J_O` must include the contributions of any rotating
  elements. As long as rotating elements have axial symmetry around their axes
  of rotation, their contributions to `J_O` will be constant.
"""
Base.@kwdef struct MassProperties
    m::Float64 = 0.0
    J_O::SMatrix{3, 3, Float64, 9} = zeros(SMatrix{3,3,Float64,9})
    r_OG::SVector{3, Float64} = zeros(SVector{3})
end


"""
Compute the `MassProperties` of a `PointMass` located at point P in reference
frame fb, given the position vector of P with respect to Ob projected in axes εb
"""
function MassProperties(p::PointMass, r_ObP_b::AbstractVector{<:Real})
    J_Ob_b = -p.m * Attitude.skew(r_ObP_b)^2
    MassProperties(p.m, J_Ob_b, r_ObP_b)
end

"""
Compute the `MassProperties` of a `PointMass` located at the origin Oc of
reference frame fc, given the `FrameTransform` from fb to fc
"""
MassProperties(p::PointMass, t_bc::FrameTransform) = MassProperties(p, t_bc.r)

"""
Return the `MassProperties` of a `RigidBody` C expressed in its own reference
frame fc.
"""
MassProperties(c::RigidBody) = MassProperties(c.m, c.J, zeros(SVector{3}))

"""
    MassProperties(c::RigidBody, t_bc::FrameTransform)

Compute the `MassProperties` of a `RigidBody` C in reference frame fb, given
the `FrameTransform` t_bc from fb to C's local reference frame fc.
"""
function MassProperties(c::RigidBody, t_bc::FrameTransform)

    m = c.m
    J_G_c = c.J

    q_bc = t_bc.q
    #tensor rotation is somewhat expensive, so if it's just a translation skip it
    if q_bc != RQuat()
        R_bc = RMatrix(q_bc)
        J_G_b = R_bc * J_G_c * R_bc'
    else
        J_G_b = J_G_c
    end

    r_ObG_b = t_bc.r
    J_Ob_b = J_G_b - m * Attitude.skew(r_ObG_b)^2

    return MassProperties(m, J_Ob_b, r_ObG_b) #p_b

end

 """
    Base.:+(p1::MassProperties, p2::MassProperties)

Aggregate the `MassProperties` of two components expressed in a common reference
frame. If neither component has mass, return a null MassProperties() instance
"""
function Base.:+(p1::MassProperties, p2::MassProperties)
    m = p1.m + p2.m
    if m > 0
        return MassProperties( m = m,
        r_OG = 1/m * (p1.m * p1.r_OG + p2.m * p2.r_OG),
        J_O = p1.J_O + p2.J_O)
    else
        return MassProperties()
    end
end

"""
    transform(t_bc::FrameTransform, p_c::MassProperties)

Transform a `MassProperties` instance from one reference frame fc to another fb.

If `t_bc` is a `FrameTransform` specifying frame fc(Oc, εc) relative to fb(Ob,
εb), and `p_c` is a `MassProperties` defined with respect to fc, then `p_b =
transform(t_bc, p_c)` is the equivalent `MassProperties` defined with respect to
fb.
"""
function transform(t_bc::FrameTransform, p_c::MassProperties)

    r_ObOc_b = t_bc.r
    q_bc = t_bc.q

    m = p_c.m
    J_Oc_c = p_c.J_O
    r_OcG_c = p_c.r_OG

    r_OcG_b = q_bc(r_OcG_c)
    r_ObG_b = r_ObOc_b + r_OcG_b

    J_G_c = J_Oc_c + m * Attitude.skew(r_OcG_c)^2

    #tensor rotation is somewhat expensive, so if it's just a translation skip it
    if q_bc != RQuat()
        R_bc = RMatrix(q_bc)
        J_G_b = R_bc * J_G_c * R_bc'
    else
        J_G_b = J_G_c
    end

    J_Ob_b = J_G_b - m * Attitude.skew(r_ObG_b)^2

    return MassProperties(m, J_Ob_b, r_ObG_b) #p_b

end

############################ MassTrait #################################

abstract type MassTrait end
struct HasMass <: MassTrait end
struct HasNoMass <: MassTrait end

"""
Notes:
- When get_mp_b is called on a System, the returned MassProperties instance must
  be expressed in the System's parent reference frame.
- At the root of the component hierarchy we have the vehicle. The vehicle is
  its own parent, so the MassProperties it returns must be expressed in its own
  reference frame: total vehicle mass, position vector from the vehicle origin
  Ob to the vehicle center of mass G expressed in vehicle axes, and inertia
  tensor of the vehicle with respect to its origin, expressed in vehicle axes.
  These are the properties expected by the dynamics equations.
- Aircraft dynamics and kinematics are formulated on the vehicle origin Ob
  instead of the overall aircraft's center of mass G. This allows for any of the
  overall aircraft's mass properties to change, either gradually (for example, due to
  fuel consumption) or suddenly (due to a payload release), without having to
  worry about discontinuities in the kinematic state vector.
"""

MassTrait(::S) where {S<:System} = error(
    "Please extend Dynamics.MassTrait for $S")

get_mp_b(sys::System) = get_mp_b(MassTrait(sys), sys)

get_mp_b(::HasNoMass, sys::System) = MassProperties()

#default implementation for a System with the HasMass trait tries to compute
#the aggregate mass properties for all its the subsystems
@inline @generated function (get_mp_b(::HasMass, sys::System{T, X, Y, U, D, P, S})
    where {T<:SystemDescriptor, X, Y, U, D, P, S})

    ex = Expr(:block)

    if isempty(fieldnames(S))
        push!(ex.args,
            :(error("System{$(T)} has the HasMass trait and no subsystems, "*
                "but it does not extend the get_mp_b method")))
    else
        push!(ex.args, :(p = MassProperties()))
        for label in fieldnames(S)
            push!(ex.args,
                :(p += get_mp_b(sys.subsystems[$(QuoteNode(label))])))
        end
    end
    return ex

end
###################### WrenchTrait ##########################

abstract type WrenchTrait end
struct GetsExternalWrench <: WrenchTrait end
struct GetsNoExternalWrench <: WrenchTrait end

#prevents the trait system from failing silently when wrongly extended
WrenchTrait(::S) where {S<:System} = error(
    "Please extend Components.WrenchTrait for $S")

get_wr_b(sys::System) = get_wr_b(WrenchTrait(sys), sys)

get_wr_b(::GetsNoExternalWrench, sys::System) = Wrench()

#default implementation for a System with the GetsExternalWrench trait, tries
#to sum all the Wrenches from its individual components. override as required
@inline @generated function (get_wr_b(::GetsExternalWrench, sys::System{T, X, Y, U, D, P, S})
    where {T<:SystemDescriptor, X, Y, U, D, P, S})

    # Core.print("Generated function called")

    ex = Expr(:block)

    if isempty(fieldnames(S))
        push!(ex.args,
            :(error("System{$(T)} has the GetsExternalWrench trait and no subsystems, "*
                "but it does not extend the get_wr_b method")))
    else
        push!(ex.args, :(wr = Wrench())) #initialize a zero wrench
        for label in fieldnames(S)
            push!(ex.args,
                :(wr += get_wr_b(sys.subsystems[$(QuoteNode(label))])))
        end
    end

    return ex

end

###################### AngularMomentumTrait ##########################

#accounts only for additional angular momentum due to rotating components
abstract type AngularMomentumTrait end
struct HasAngularMomentum <: AngularMomentumTrait end
struct HasNoAngularMomentum <: AngularMomentumTrait end

#prevents the trait system from failing silently when wrongly extended
AngularMomentumTrait(::S) where {S<:System} = error(
    "Please extend Dynamics.AngularMomentumTrait for $S")

get_hr_b(sys::System) = get_hr_b(AngularMomentumTrait(sys), sys)

get_hr_b(::HasNoAngularMomentum, sys::System) = zeros(SVector{3})

#default implementation for a System with the HasAngularMomentum trait, tries to
#sum the angular momentum from its individual components. override as required
@inline @generated function (get_hr_b(::HasAngularMomentum, sys::System{T, X, Y, U, D, P, S})
    where {T<:SystemDescriptor, X, Y, U, D, P, S})

    # Core.print("Generated function called")
    ex = Expr(:block)

    if isempty(fieldnames(S))
        push!(ex.args,
            :(error("System{$(T)} has the HasAngularMomentum trait and no subsystems, "*
                "but it does not extend the get_hr_b method")))
    else
        push!(ex.args, :(h = SVector(0., 0., 0.))) #initialize
        for label in fieldnames(S)
            push!(ex.args,
                :(h += get_hr_b(sys.subsystems[$(QuoteNode(label))])))
        end
    end

    return ex

end

########################### Dynamic Equations ################################

"""
    inertia_wrench(mass::MassProperties, vel::VelData, hr_b::AbstractVector{<:Real})

Compute the equivalent `Wrench` arising from inertia terms in the dynamic
equations

The resulting `Wrench` is defined on the vehicle's reference frame fb.

# Arguments:
- `mp_b::MassProperties`: Current aircraft mass properties in frame fb
- `vel::VelData`: Velocity outputs
- `hr_b::AbstractVector{<:Real}`: Additional angular momentum due to the
  angular velocity of any rotating elements with respect to the vehicle,
  projected on vehicle axes

"""
function inertia_wrench(mp_b::MassProperties, kin::Kinematics.Common, hr_b::AbstractVector{<:Real})

    @unpack ω_ie_b, ω_eb_b, ω_ib_b, v_eOb_b = kin

    m = mp_b.m; J_Ob_b = mp_b.J_O; r_ObG_b = mp_b.r_OG

    #angular momentum of the vehicle as a rigid body (excluding rotating
    #components)
    h_rbd_b = J_Ob_b * ω_ib_b

    #total angular momentum
    h_all_b = h_rbd_b + SVector{3,Float64}(hr_b)

    #inertia terms (exact... overkill, but very cheap anyway)
    a_1_b = (ω_eb_b + 2 * ω_ie_b) × v_eOb_b
    F_in_Ob_b = -m * (a_1_b + ω_ib_b × (ω_ib_b × r_ObG_b) + r_ObG_b × (ω_eb_b × ω_ie_b ))
    M_in_Ob_b = - ( J_Ob_b * (ω_ie_b × ω_eb_b) + ω_ib_b × h_all_b + m * r_ObG_b × a_1_b)

    return Wrench(F = F_in_Ob_b, M = M_in_Ob_b)

end

function gravity_wrench(mp_b::MassProperties, kin::Kinematics.Common)

    #gravity can be viewed as an entity acting on a local frame with its origin
    #at G and its axes aligned with the local tangent frame

    #strictly, the gravity vector should be evaluated at G, with its direction
    #given by the z-axis of LTF(G). however, since g(G) ≈ g(Ob) and LTF(G) ≈
    #LTF(Ob), we can instead evaluate g at Ob, assuming its direction given by
    #LTF(Ob), and then apply it at G.
    @unpack n_e, h_e, q_nb = kin

    Ob = GeographicLocation(n_e, h_e)
    g_G_n = g_Ob_n = g_n(Ob)

    #the resultant consists of the gravity force acting on G along the local
    #vertical and a null torque
    F_G_n = mp_b.m * g_G_n
    M_G_n = zeros(SVector{3})
    wr_G_n = Wrench(F = F_G_n, M = M_G_n)

    #with the previous assumption, the transformation from body frame to local
    #gravity frame is given by the translation r_ObG_b and the (passive)
    #rotation from b to LTF(Ob) (instead of LTF(G)), which is given by q_lb'
    wr_c = wr_G_n
    r_ObG_b = mp_b.r_OG
    t_bc = FrameTransform(r = r_ObG_b, q = q_nb')
    return t_bc(wr_c) #wr_b

end


###################### Acceleration Outputs #####################

Base.@kwdef struct DynDataIn
    wr_g_b::Wrench = Wrench()
    wr_in_b::Wrench = Wrench()
    wr_ext_b::Wrench = Wrench()
    hr_b::SVector{3,Float64} = zeros(SVector{3})
end

Base.@kwdef struct DynDataOut
    α_eb_b::SVector{3,Float64} = zeros(SVector{3})
    α_ib_b::SVector{3,Float64} = zeros(SVector{3})
    a_eOb_b::SVector{3,Float64} = zeros(SVector{3})
    a_eOb_n::SVector{3,Float64} = zeros(SVector{3})
    a_iOb_b::SVector{3,Float64} = zeros(SVector{3})
    f_Ob_b::SVector{3,Float64} = zeros(SVector{3}) #specific force (g) maybe two y axes??
end

Base.@kwdef struct DynData
    input::DynDataIn = DynDataIn()
    output::DynDataOut = DynDataOut()
end


function f_dyn!(ẋ_vel::Kinematics.XVel, kin::Kinematics.Common, mp_b::MassProperties,
    wr_ext_b::Wrench, hr_b::AbstractVector{<:Real})

    #wr_ext_b: Total external wrench on the vehicle

    #hr_b: Additional angular momentum due to the angular velocity of the
    #rotating aircraft components with respect to the vehicle. these are
    #computed individually by each component relative to its center of mass and
    #then summed

    #wr_ext_b and hr_b, as well as mass data, are produced by aircraft
    #components, so they must be computed by the aircraft's x_dot method. and,
    #since kin is needed by those components, it must be called from the
    #aircraft's kinematic state vector

    #clearly, r_ObG_b cannot be arbitrarily large, because J_Ob_b is larger than
    #J_G_b (Steiner). therefore, at some point J_G_b would become zero (or at
    #least singular)!

    @unpack q_eb, q_nb, n_e, h_e, ω_eb_b, ω_ie_b, v_eOb_b = kin

    m = mp_b.m; J_Ob_b = mp_b.J_O; r_ObG_b = mp_b.r_OG

    r_ObG_b_sk = Attitude.skew(r_ObG_b)
    A11 = J_Ob_b
    A12 = m * r_ObG_b_sk
    A21 = -m * r_ObG_b_sk
    A22 = m * SMatrix{3,3,Float64}(I)

    A = vcat(hcat(A11, A12), hcat(A21, A22))

    wr_g_b = gravity_wrench(mp_b, kin)
    wr_in_b = inertia_wrench(mp_b, kin, hr_b)
    wr_b = wr_ext_b + wr_g_b + wr_in_b
    b = SVector{6}(vcat(wr_b.M, wr_b.F))

    # update ẋ_vel
    ẋ_vel .= A\b

    #compute outputs
    Ob = GeographicLocation(n_e, h_e)
    r_eOb_e = CartesianLocation(Ob)[:]
    r_eOb_b = q_eb'(r_eOb_e)
    v̇_eOb_b = SVector{3}(ẋ_vel.v_eOb_b)

    #angular accelerations
    α_eb_b = SVector{3}(ẋ_vel.ω_eb_b) #α_eb_b == ω_eb_b_dot
    α_ib_b = α_eb_b - ω_eb_b × ω_ie_b

    #linear accelerations and specific force
    a_eOb_b = v̇_eOb_b + ω_eb_b × v_eOb_b
    a_eOb_n = q_nb(a_eOb_b)
    a_iOb_b = v̇_eOb_b + (ω_eb_b + 2ω_ie_b) × v_eOb_b + ω_ie_b × (ω_ie_b × r_eOb_b)

    g_Ob_b = q_nb'(g_n(Ob))
    G_Ob_b = g_Ob_b + ω_ie_b × (ω_ie_b × r_eOb_b)
    f_Ob_b = a_iOb_b - G_Ob_b

    data_in = DynDataIn(wr_g_b, wr_in_b, wr_ext_b, hr_b)
    data_out = DynDataOut(α_eb_b, α_ib_b, a_eOb_b, a_eOb_n, a_iOb_b, f_Ob_b)

    return DynData(data_in, data_out)

end

################################# Dynamics #####################################

@recipe function f(th::TimeHistory{<:Wrench}; wr_frame = "", wr_source = "")

    layout := (1, 2)
    seriestype --> :path

    @series begin
        subplot := 1
        title --> "Force"
        yguide --> L"$F_{O%$wr_frame \ (%$wr_source)}^{%$wr_frame} \ (N)$"
        th_split --> :none
        th.F
    end

    @series begin
        subplot := 2
        title --> "Torque"
        yguide --> L"$M_{O%$wr_frame \ (%$wr_source)}^{%$wr_frame} \ (N \ m)$"
        th_split --> :none
        th.M
    end

end

function make_plots(th::TimeHistory{<:DynData}; kwargs...)

    return OrderedDict(
        :input => make_plots(th.input; kwargs...),
        :output => make_plots(th.output; kwargs...)
    )

end

function make_plots(th::TimeHistory{<:DynDataIn}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    pd[:wr_g_b] = plot(th.wr_g_b;
        plot_title = "Gravity Wrench [Vehicle Axes]",
        wr_source = "g", wr_frame = "b",
        kwargs...)

    pd[:wr_in_b] = plot(th.wr_in_b;
        plot_title = "Inertia Wrench [Vehicle Axes]",
        wr_source = "in", wr_frame = "b",
        kwargs...)

    pd[:wr_ext_b] = plot(th.wr_ext_b;
        plot_title = "External Wrench [Vehicle Axes]",
        wr_source = "ext", wr_frame = "b",
        kwargs...)

    pd[:hr_b] = plot(th.hr_b;
        plot_title = "Angular Momentum from Rotating Components [Vehicle Axes]",
        ylabel = hcat(
            L"$h_{Ob \ (r)}^{x_b} \ (kg \ m^2 / s)$",
            L"$h_{Ob \ (r)}^{y_b} \ (kg \ m^2 / s)$",
            L"$h_{Ob \ (r)}^{z_b} \ (kg \ m^2 / s)$"),
        th_split = :h, link = :none,
        kwargs...)

    return pd

end

function make_plots(th::TimeHistory{<:DynDataOut}; kwargs...)

    #standard gravity for specific force normalization
    g₀ = 9.80665

    pd = OrderedDict{Symbol, Plots.Plot}()

    pd[:α_eb_b] = plot(th.α_eb_b;
        plot_title = "Angular Acceleration (Vehicle/ECEF) [Vehicle Axes]",
        ylabel = hcat(
            L"$\alpha_{eb}^{x_b} \ (rad/s^2)$",
            L"$\alpha_{eb}^{y_b} \ (rad/s^2)$",
            L"$\alpha_{eb}^{z_b} \ (rad/s^2)$"),
        th_split = :h,
        kwargs...)

    pd[:a_eOb_b] = plot(th.a_eOb_b;
        plot_title = "Linear Acceleration (Vehicle/ECEF) [Vehicle Axes]",
        ylabel = hcat(
            L"$a_{eb}^{x_b} \ (m/s^{2})$",
            L"$a_{eb}^{y_b} \ (m/s^{2})$",
            L"$a_{eb}^{z_b} \ (m/s^{2})$"),
        th_split = :h,
        kwargs...)

    pd[:a_eOb_n] = plot(th.a_eOb_n;
        plot_title = "Linear Acceleration (Vehicle/ECEF) [NED Axes]",
        ylabel = hcat(
            L"$a_{eb}^{N} \ (m/s^{2})$",
            L"$a_{eb}^{E} \ (m/s^{2})$",
            L"$a_{eb}^{D} \ (m/s^{2})$"),
        th_split = :h, link = :none,
        kwargs...)

    pd[:f_Ob_b] = plot(TimeHistory(th._t, th.f_Ob_b._data / g₀);
        plot_title = "Specific Force [Vehicle Axes]",
        ylabel = hcat(
            L"$f_{Ob}^{x_b} \ (g)$",
            L"$f_{Ob}^{y_b} \ (g)$",
            L"$f_{Ob}^{z_b} \ (g)$"),
        th_split = :h,
        kwargs...)

    return pd

end

end #module