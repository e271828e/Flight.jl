module Dynamics

using StaticArrays, LinearAlgebra, UnPack

using Flight.FlightCore

using ..Attitude
using ..Geodesy
using ..Kinematics

export FrameTransform, transform
export Wrench
export AbstractMassDistribution, PointDistribution, RigidBodyDistribution, MassProperties
export RigidBodyData, RigidBodyDynamics, RigidBodyDynamicsY
export get_mp_b, get_wr_b, get_hr_b

#standard gravity for specific force normalization
const g₀ = 9.80665


################################################################################
############################# FrameTransform ###################################

"""
Specifies a reference frame `fc(Oc, Ɛc)` relative to another `fb(Ob, Ɛb)`

Frame `fc(Oc, Ɛc)` is defined by:
- The position vector from fb's origin Ob to fc's origin Oc, projected on fb's
  axes εb (`r_ObOc_b`)
- The rotation quaternion from fb's axes εb to fc's axes εc (`q_bc`)
"""
@kwdef struct FrameTransform
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



################################################################################
################################## Wrench ######################################

"""
Force and torque combination defined on a concrete reference frame

A `Wrench` is defined on reference frame fc(Oc,εc) when it is applied at its
origin Oc and projected on its axes εc
"""
@kwdef struct Wrench
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


################################################################################
############################## MassDistributions ###############################

abstract type AbstractMassDistribution end

struct PointDistribution <: AbstractMassDistribution
    m::Float64
end

"""Defines a rigid body mass distribution.

A `RigidBodyDistribution` is defined in a local reference frame fc, whose origin Oc is
located at the body's center of mass G, and whose axes εc are arbitrarily chosen
to express the inertia tensor conveniently (typically they will be the body's
principal axes of inertia)
- m: Total mass
- J: Inertia tensor computed with respect to Oc and projected in εc
"""
@kwdef struct RigidBodyDistribution <: AbstractMassDistribution
    m::Float64 = 1.0
    J::SMatrix{3,3,Float64,9} = diagm(ones(SVector{3}))
end


################################################################################
################################ MassProperties ################################


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
@kwdef struct MassProperties
    m::Float64 = 0.0
    J_O::SMatrix{3, 3, Float64, 9} = zeros(SMatrix{3,3,Float64,9})
    r_OG::SVector{3, Float64} = zeros(SVector{3})
end


"""
Compute the `MassProperties` of a `PointDistribution` located at point P in reference
frame fb, given the position vector of P with respect to Ob projected in axes εb
"""
function MassProperties(p::PointDistribution, r_ObP_b::AbstractVector{<:Real})
    J_Ob_b = -p.m * Attitude.v2skew(r_ObP_b)^2
    MassProperties(p.m, J_Ob_b, r_ObP_b)
end

"""
Compute the `MassProperties` of a `PointDistribution` located at the origin Oc of
reference frame fc, given the `FrameTransform` from fb to fc
"""
MassProperties(p::PointDistribution, t_bc::FrameTransform) = MassProperties(p, t_bc.r)

"""
Return the `MassProperties` of a `RigidBodyDistribution` C expressed in its own reference
frame fc.
"""
MassProperties(c::RigidBodyDistribution) = MassProperties(c.m, c.J, zeros(SVector{3}))

"""
    MassProperties(c::RigidBodyDistribution, t_bc::FrameTransform)

Compute the `MassProperties` of a `RigidBodyDistribution` C in reference frame fb, given
the `FrameTransform` t_bc from fb to C's local reference frame fc.
"""
function MassProperties(c::RigidBodyDistribution, t_bc::FrameTransform)

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
    J_Ob_b = J_G_b - m * Attitude.v2skew(r_ObG_b)^2

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

    J_G_c = J_Oc_c + m * Attitude.v2skew(r_OcG_c)^2

    #tensor rotation is somewhat expensive, so if it's just a translation skip it
    if q_bc != RQuat()
        R_bc = RMatrix(q_bc)
        J_G_b = R_bc * J_G_c * R_bc'
    else
        J_G_b = J_G_c
    end

    J_Ob_b = J_G_b - m * Attitude.v2skew(r_ObG_b)^2

    return MassProperties(m, J_Ob_b, r_ObG_b) #p_b

end


################################################################################
############################### RigidBodyData ##################################

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


#mp_Ob: Mass properties of the System, translated to the vehicle frame

#wr_ext_Ob: Total external wrench contributed to the System, resolved in the
#vehicle frame

#hr_b: Internal angular momentum due to the angular velocity of the rotating
#System's components with respect to the vehicle. Computed individually by
#each component relative to its center of mass and then summed

#note: for a given J_Ob_b, r_ObG_b cannot be arbitrarily large, because
#moments of inertia are minimum at G. therefore, at some point J_G_b would
#become non-positive definite

@kwdef struct RigidBodyData
    mp_Ob::MassProperties = MassProperties(RigidBodyDistribution())
    wr_ext_Ob::Wrench = Wrench()
    hr_b::SVector{3,Float64} = zeros(SVector{3})
end

function RigidBodyData(sys::System)
    mp_Ob = get_mp_b(sys)
    wr_ext_Ob = get_wr_b(sys)
    hr_b = get_hr_b(sys)
    RigidBodyData(mp_Ob, wr_ext_Ob, hr_b)
end

#generated functions are needed here because type inference does not work
#throughout the whole System hierarchy

# default implementation tries to compute the aggregate mass properties for all
# its subsystems. override if possible to reduce compilation time
@inline @generated function (get_mp_b(sys::System{T, X, Y, U, S, P, B})
    where {T<:SystemDefinition, X, Y, U, S, P, B})

    ex = Expr(:block)

    if isempty(fieldnames(B))
        push!(ex.args,
            :(error("System{$(T)} is a leaf System "*
                "but it does not extend the get_mp_b method")))
    else
        push!(ex.args, :(p = MassProperties()))
        for label in fieldnames(B)
            push!(ex.args,
                :(p += get_mp_b(sys.subsystems[$(QuoteNode(label))])))
        end
    end
    return ex

end


#default implementation tries to sum the angular momentum from its individual
#components. override if possible to reduce compilation time
@inline @generated function (get_hr_b(sys::System{T, X, Y, U, S, P, B})
    where {T<:SystemDefinition, X, Y, U, S, P, B})

    # Core.print("Generated function called")
    ex = Expr(:block)

    if isempty(fieldnames(B))
        push!(ex.args,
            :(error("System{$(T)} is a leaf System "*
                "but it does not extend the get_hr_b method")))
    else
        push!(ex.args, :(h = SVector(0., 0., 0.))) #initialize
        for label in fieldnames(B)
            push!(ex.args,
                :(h += get_hr_b(sys.subsystems[$(QuoteNode(label))])))
        end
    end

    return ex

end

#default implementation tries to sum all the Wrenches from its individual
#components. override if possible to reduce compilation time
@inline @generated function (get_wr_b(sys::System{T, X, Y, U, S, P, B})
    where {T<:SystemDefinition, X, Y, U, S, P, B})

    # Core.print("Generated function called")

    ex = Expr(:block)

    if isempty(fieldnames(B))
        push!(ex.args,
            :(error("System{$(T)} is a leaf System, "*
                "but it does not extend the get_wr_b method")))
    else
        push!(ex.args, :(wr = Wrench())) #initialize a zero wrench
        for label in fieldnames(B)
            push!(ex.args,
                :(wr += get_wr_b(sys.subsystems[$(QuoteNode(label))])))
        end
    end

    return ex

end


################################################################################
############################### Dynamics #######################################

struct RigidBodyDynamics <: SystemDefinition end

#all magnitudes resolved in body axes unless otherwise noted
@kwdef struct RigidBodyDynamicsY
    mp_Ob::MassProperties = MassProperties() #aircraft mass properties at Ob
    mp_Gb::MassProperties = MassProperties() #aircraft mass properties at G
    wr_g_Ob::Wrench = Wrench() #gravity wrench at Ob
    wr_in_Ob::Wrench = Wrench() #inertia wrench at Ob
    wr_ext_Ob::Wrench = Wrench() #externally applied wrench at Ob
    wr_net_Ob::Wrench = Wrench() #net wrench at Ob
    wr_net_Gb::Wrench = Wrench() #net wrench at G
    hr_b::SVector{3,Float64} = zeros(SVector{3}) #angular momentum from rotating components
    α_eb_b::SVector{3,Float64} = zeros(SVector{3}) #ECEF-to-body angular acceleration
    α_ib_b::SVector{3,Float64} = zeros(SVector{3}) #ECI-to-body angular acceleration
    v̇_eOb_b::SVector{3,Float64} = zeros(SVector{3}) #time derivative of ECEF-relative velocity
    a_eOb_b::SVector{3,Float64} = zeros(SVector{3}) #ECEF-relative acceleration of Ob
    a_eOb_n::SVector{3,Float64} = zeros(SVector{3}) #ECEF-relative acceleration of Ob, NED axes
    a_iOb_b::SVector{3,Float64} = zeros(SVector{3}) #ECI-relative acceleration of Ob
    a_iG_b::SVector{3,Float64} = zeros(SVector{3}) #ECI-relative acceleration of G
    f_Gb_b::SVector{3,Float64} = zeros(SVector{3}) #specific force at G
end

Systems.X(::RigidBodyDynamics) = zero(Kinematics.XVelTemplate)
Systems.Y(::RigidBodyDynamics) = RigidBodyDynamicsY()


function Systems.f_ode!(sys::System{RigidBodyDynamics}, kin_data::KinData, rb_data::RigidBodyData)

    @unpack q_eb, q_nb, n_e, h_e, r_eOb_e, ω_eb_b, ω_ie_b, ω_ib_b, v_eOb_b = kin_data
    @unpack mp_Ob, wr_ext_Ob, hr_b = rb_data
    @unpack ẋ = sys

    m = mp_Ob.m; J_Ob_b = mp_Ob.J_O; r_ObGb_b = mp_Ob.r_OG

    ######################## Gravity Wrench ####################################

    #compute geographic position of Gb
    r_ObGb_e = q_eb(r_ObGb_b)
    r_eGb_e = r_eOb_e + r_ObGb_e
    Gb = Cartesian(r_eGb_e)

    #create an auxiliary frame c with axes parallel to the local NED frame at
    #Gb to compute gravity vector at Gb
    q_ec = ltf(Gb)
    q_bc = q_eb' ∘ q_ec
    g_Gb_c = SVector{3,Float64}(0, 0, gravity(Gb))
    g_Gb_b = q_bc(g_Gb_c)

    #the resultant from gravity on the vehicle at its center of gravity Gb
    #consists of gravity force, plus a null torque
    F_Gb_b = m * g_Gb_b
    M_Gb_b = zeros(SVector{3})
    wr_g_Gb = Wrench(F = F_Gb_b, M = M_Gb_b)

    #define pure translation from Ob to Gb
    t_ObGb = FrameTransform(r = r_ObGb_b) #Gb to Ob

    #translate gravity wrench to Ob
    wr_g_Ob = t_ObGb(wr_g_Gb)

    ########################## Inertia Wrench ##################################

    #angular momentum of the vehicle as a rigid body (excluding rotating
    #components)
    h_rbd_b = J_Ob_b * ω_ib_b

    #total angular momentum
    h_all_b = h_rbd_b + SVector{3,Float64}(hr_b)

    #inertia terms (exact version... maybe overkill, but very cheap)
    a_1_b = (ω_eb_b + 2 * ω_ie_b) × v_eOb_b
    F_in_Ob_b = -m * (a_1_b + ω_ib_b × (ω_ib_b × r_ObGb_b) + r_ObGb_b × (ω_eb_b × ω_ie_b ))
    M_in_Ob_b = - ( J_Ob_b * (ω_ie_b × ω_eb_b) + ω_ib_b × h_all_b + m * r_ObGb_b × a_1_b)

    wr_in_Ob = Wrench(F = F_in_Ob_b, M = M_in_Ob_b)

    ########################### Dynamic Equations ##############################

    r_ObGb_b_sk = Attitude.v2skew(r_ObGb_b)
    A11 = J_Ob_b
    A12 = m * r_ObGb_b_sk
    A21 = -m * r_ObGb_b_sk
    A22 = m * SMatrix{3,3,Float64}(I)

    A = vcat(hcat(A11, A12), hcat(A21, A22)) #mass matrix

    wr_net_Ob = wr_ext_Ob + wr_g_Ob + wr_in_Ob
    b = SVector{6}(vcat(wr_net_Ob.M, wr_net_Ob.F)) #forcing vector

    ẋ .= A\b

    ########################## Additional Outputs ##############################

    #angular accelerations
    α_eb_b = SVector{3}(ẋ.ω_eb_b) #α_eb_b == ω_eb_b_dot
    α_ib_b = α_eb_b - ω_eb_b × ω_ie_b

    #linear accelerations at Ob
    v̇_eOb_b = SVector{3}(ẋ.v_eOb_b)
    r_eOb_b = q_eb'(r_eOb_e)
    a_eOb_b = v̇_eOb_b + ω_eb_b × v_eOb_b
    a_eOb_n = q_nb(a_eOb_b)
    a_iOb_b = v̇_eOb_b + (ω_eb_b + 2ω_ie_b) × v_eOb_b + ω_ie_b × (ω_ie_b × r_eOb_b)

    #define pure translation from Gb to Ob
    t_GbOb = FrameTransform(r = -r_ObGb_b) #Gb to Ob

    #translate mass properties and net wrench to Gb
    mp_Gb = t_GbOb(mp_Ob)
    wr_net_Gb = t_GbOb(wr_net_Ob)

    #gravitation at Gb, body axes
    r_eGb_b = q_eb'(r_eGb_e)
    G_Gb_b = g_Gb_b + ω_ie_b × (ω_ie_b × r_eGb_b)

    #linear acceleration and specific force at Gb
    a_iG_b = a_iOb_b + ω_ib_b × (ω_ib_b × r_ObGb_b) + α_ib_b × r_ObGb_b
    f_Gb_b = a_iG_b - G_Gb_b

    sys.y = RigidBodyDynamicsY(; mp_Ob, mp_Gb, wr_g_Ob, wr_in_Ob, wr_ext_Ob,
        wr_net_Ob, wr_net_Gb, hr_b, α_eb_b, α_ib_b, v̇_eOb_b, a_eOb_b, a_eOb_n,
        a_iOb_b, a_iG_b, f_Gb_b)

end

################################# Dynamics #####################################

@recipe function f(ts::TimeSeries{<:Wrench}; wr_frame = "", wr_source = "")

    layout := (1, 2)
    seriestype --> :path

    @series begin
        subplot := 1
        title --> "Force"
        yguide --> L"$F_{O%$wr_frame \ (%$wr_source)}^{%$wr_frame} \ (N)$"
        ts_split --> :none
        ts.F
    end

    @series begin
        subplot := 2
        title --> "Torque"
        yguide --> L"$M_{O%$wr_frame \ (%$wr_source)}^{%$wr_frame} \ (N \ m)$"
        ts_split --> :none
        ts.M
    end

end

function Plotting.make_plots(ts::TimeSeries{<:RigidBodyDynamicsY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    pd[:wr_g_Ob] = plot(ts.wr_g_Ob;
        plot_title = "Gravity Wrench at Ob [Vehicle Axes]",
        wr_source = "g", wr_frame = "b",
        kwargs...)

    pd[:wr_in_Ob] = plot(ts.wr_in_Ob;
        plot_title = "Inertia Wrench at Ob [Vehicle Axes]",
        wr_source = "in", wr_frame = "b",
        kwargs...)

    pd[:wr_ext_Ob] = plot(ts.wr_ext_Ob;
        plot_title = "External Wrench at Ob [Vehicle Axes]",
        wr_source = "ext", wr_frame = "b",
        kwargs...)

    pd[:wr_net_Ob] = plot(ts.wr_net_Ob;
        plot_title = "Total Wrench at Ob [Vehicle Axes]",
        wr_source = "ext", wr_frame = "b",
        kwargs...)

    pd[:hr_b] = plot(ts.hr_b;
        plot_title = "Angular Momentum from Rotating Components [Vehicle Axes]",
        ylabel = hcat(
            L"$h_{Ob \ (r)}^{x_b} \ (kg \ m^2 / s)$",
            L"$h_{Ob \ (r)}^{y_b} \ (kg \ m^2 / s)$",
            L"$h_{Ob \ (r)}^{z_b} \ (kg \ m^2 / s)$"),
        ts_split = :h, link = :none,
        kwargs...)

    pd[:α_eb_b] = plot(ts.α_eb_b;
        plot_title = "Angular Acceleration (Vehicle/ECEF) [Vehicle Axes]",
        ylabel = hcat(
            L"$\alpha_{eb}^{x_b} \ (rad/s^2)$",
            L"$\alpha_{eb}^{y_b} \ (rad/s^2)$",
            L"$\alpha_{eb}^{z_b} \ (rad/s^2)$"),
        ts_split = :h,
        kwargs...)

    pd[:a_eOb_b] = plot(ts.a_eOb_b;
        plot_title = "Linear Acceleration (Ob/ECEF) [Vehicle Axes]",
        ylabel = hcat(
            L"$a_{eb}^{x_b} \ (m/s^{2})$",
            L"$a_{eb}^{y_b} \ (m/s^{2})$",
            L"$a_{eb}^{z_b} \ (m/s^{2})$"),
        ts_split = :h,
        kwargs...)

    pd[:a_eOb_n] = plot(ts.a_eOb_n;
        plot_title = "Linear Acceleration (Ob/ECEF) [NED Axes]",
        ylabel = hcat(
            L"$a_{eb}^{N} \ (m/s^{2})$",
            L"$a_{eb}^{E} \ (m/s^{2})$",
            L"$a_{eb}^{D} \ (m/s^{2})$"),
        ts_split = :h, link = :none,
        kwargs...)

    pd[:f_Gb_b] = plot(TimeSeries(ts._t, ts.f_Gb_b._data / g₀);
        plot_title = "Specific Force (Center of Mass) [Vehicle Axes]",
        ylabel = hcat(
            L"$f_{G}^{x_b} \ (g)$",
            L"$f_{G}^{y_b} \ (g)$",
            L"$f_{G}^{z_b} \ (g)$"),
        ts_split = :h, link = :none,
        kwargs...)

    return pd

end

################################################################################
################################# GUI ##########################################

function GUI.draw(sys::System{RigidBodyDynamics}, p_open::Ref{Bool} = Ref(true),
                    label::String = "Rigid Body")

    @unpack mp_Ob, mp_Gb, wr_net_Ob, wr_net_Gb, hr_b, α_eb_b, a_eOb_b, f_Gb_b = sys.y

    CImGui.Begin(label, p_open)

    CImGui.Text(@sprintf("Mass: %.3f kg", mp_Ob.m))
    GUI.draw(mp_Ob.r_OG, "CG Position (O) [Body]", "m")

    #inertia tensor is given at G instead of O
    if CImGui.TreeNode("Inertia Tensor (G) [Body]")
        CImGui.Text(@sprintf("XX: %.3f kg m2", mp_Gb.J_O[1,1]))
        CImGui.Text(@sprintf("YY: %.3f kg m2", mp_Gb.J_O[2,2]))
        CImGui.Text(@sprintf("ZZ: %.3f kg m2", mp_Gb.J_O[3,3]))
        CImGui.Text(@sprintf("XY: %.3f kg m2", mp_Gb.J_O[1,2]))
        CImGui.Text(@sprintf("XZ: %.3f kg m2", mp_Gb.J_O[1,3]))
        CImGui.Text(@sprintf("YZ: %.3f kg m2", mp_Gb.J_O[2,3]))
        CImGui.TreePop()
    end

    #net force and torque are given at G instead of O
    GUI.draw(wr_net_Gb.F, "Net Force (G) [Body]", "N")
    GUI.draw(wr_net_Gb.M, "Net Torque (G) [Body]", "N*m")
    GUI.draw(hr_b, "Angular Momentum From Rotating Components [Body]", "kg*(m^2)/s")
    GUI.draw(α_eb_b, "Angular Acceleration (Body / ECEF) [Body]", "rad/(s^2)")
    GUI.draw(a_eOb_b, "Linear Acceleration (O / ECEF) [Body]", "m/(s^2)")
    GUI.draw(f_Gb_b/g₀, "Specific Force (G) [Body]", "g")

    CImGui.End()

end


end #module