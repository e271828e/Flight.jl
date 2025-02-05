module Dynamics

using StaticArrays, LinearAlgebra, UnPack

using Flight.FlightCore

using ..Attitude
using ..Geodesy
using ..Kinematics

export FrameTransform, transform
export Wrench
export AbstractMassDistribution, PointDistribution, RigidBodyDistribution, MassProperties
export RigidBodyDynamics, Accelerations, Actions
export get_mp_b, get_wr_b, get_hr_b

#standard gravity for specific force normalization
const g₀ = 9.80665


################################################################################
############################# FrameTransform ###################################

"""
Specifies a reference frame `fc(Oc, Ɛc)` relative to another `fb(Ob, Ɛb)`

Frame `fc(Oc, Ɛc)` is defined by:
- The position vector from fb's origin Ob to fc's origin Oc, projected on fb's
  axes εb (`r_bc_b`)
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

    r_bc_b = t_bc.r; q_bc = t_bc.q
    r_cd_c = t_cd.r; q_cd = t_cd.q

    r_bd_b = r_bc_b + q_bc(r_cd_c)
    q_bd = q_bc ∘ q_cd

    t_bd = FrameTransform(r_bd_b, q_bd)

    return t_bd

end

"""
    Base.:adjoint(t_bc::FrameTransform)

Get the reciprocal `FrameTransform`.
"""
function Base.:adjoint(t_bc::FrameTransform)

    r_bc_b = t_bc.r; q_bc = t_bc.q

    q_cb = q_bc'
    r_cb_c = q_cb(-r_bc_b)
    t_cb = FrameTransform(r_cb_c, q_cb)

    return t_cb
end


"""
    transform(t_bc::FrameTransform, r_cP_c::AbstractVector{<:Real})

Given the (3D) position vector of some point P in reference frame fc, and the
`FrameTransform` from fb to fc, compute the position vector of P in fb (r_bP_b)
"""
function transform(t_bc::FrameTransform, r_cP_c::AbstractVector{<:Real})

    r_bc_b = t_bc.r; q_bc = t_bc.q

    r_bP_b = r_bc_b + q_bc(SVector{3,Float64}(r_cP_c))
    return r_bP_b
end


"""Alternative function call notation: `t_bc(mp_c) == transform(t_bc, mp_c)`
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

    F_c_c = wr_c.F
    M_c_c = wr_c.M

    #project onto b axes
    F_c_b = t_bc.q(F_c_c)
    M_c_b = t_bc.q(M_c_c)

    #translate to b origin
    F_b_b = F_c_b
    M_b_b = M_c_b + t_bc.r × F_c_b

    return Wrench(F = F_b_b, M = M_b_b) #wr_b

end


################################################################################
############################## MassDistributions ###############################

abstract type AbstractMassDistribution end

struct PointDistribution <: AbstractMassDistribution
    m::Float64
end

"""Defines a rigid body mass distribution.

A `RigidBodyDistribution` is defined in a local reference frame fc, whose origin
Oc is located at the body's center of mass G, and whose axes εc are arbitrarily
chosen to express the inertia tensor conveniently (typically they will be the
body's principal axes of inertia)
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
Mass properties of a component expressed on some reference frame fb(Ob,εb):
- `m`: Total mass
- `J`: Overall inertia tensor with respect to Ob, projected on axes εb (more
  explicitly, J_Ob_b)
- `r_OG`: Position vector from the origin Ob to the component's center of mass
  G, projected on axes εb (more explicitly r_ObG_b)

Notes:
- The mass and inertia tensor must include the contributions of any rotating
  elements. As long as rotating elements have axial symmetry around their axes
  of rotation, their contributions to `J` will be constant.
"""
@kwdef struct MassProperties
    m::Float64 = 0.0
    J::SMatrix{3, 3, Float64, 9} = zeros(SMatrix{3,3,Float64,9})
    r_OG::SVector{3, Float64} = zeros(SVector{3})
end


"""
Compute the `MassProperties` of a `PointDistribution` located at point P in reference
frame fb, given the position vector of P with respect to b projected in axes εb
"""
function MassProperties(p::PointDistribution, r_bP_b::AbstractVector{<:Real})
    J_b_b = -p.m * Attitude.v2skew(r_bP_b)^2
    MassProperties(p.m, J_b_b, r_bP_b)
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

Compute the `MassProperties` of a `RigidBodyDistribution` C in reference frame
fb given the `FrameTransform` t_bc from fb to C's local reference frame fc(G,εc)
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

    #Steiner
    r_bG_b = t_bc.r
    J_b_b = J_G_b - m * Attitude.v2skew(r_bG_b)^2

    return MassProperties(m, J_b_b, r_bG_b)

end

 """
    Base.:+(p1::MassProperties, p2::MassProperties)

Aggregate the `MassProperties` of two components expressed in a common reference
frame. If neither component has mass, return a null MassProperties() instance
"""
function Base.:+(p1::MassProperties, p2::MassProperties)
    m = p1.m + p2.m
    if m > 0
        return MassProperties(
            m = m,
            J = p1.J + p2.J,
            r_OG = 1/m * (p1.m * p1.r_OG + p2.m * p2.r_OG))
    else
        return MassProperties()
    end
end

"""
    transform(t_bc::FrameTransform, mp_c::MassProperties)

Transform a `MassProperties` instance from its current frame fc to another fb.

If `t_bc` is a `FrameTransform` specifying frame fc(Oc, εc) relative to fb(Ob,
εb), and `mp_c` is a `MassProperties` defined with respect to fc, then `mp_b =
transform(t_bc, p_c)` is the equivalent `MassProperties` defined with respect to
fb.
"""
function transform(t_bc::FrameTransform, mp_c::MassProperties)

    r_bc_b = t_bc.r
    q_bc = t_bc.q

    m = mp_c.m
    J_c_c = mp_c.J
    r_cG_c = mp_c.r_OG

    #translate inertia tensor to G in c axes
    J_G_c = J_c_c + m * Attitude.v2skew(r_cG_c)^2

    #rotate inertia tensor to b axes
    if q_bc != RQuat()
        R_bc = RMatrix(q_bc)
        J_G_b = R_bc * J_G_c * R_bc'
    else
        J_G_b = J_G_c
    end

    #get position vector of G in b frame
    r_cG_b = q_bc(r_cG_c)
    r_bG_b = r_bc_b + r_cG_b

    #translate inertia tensor to b frame origin
    J_b_b = J_G_b - m * Attitude.v2skew(r_bG_b)^2

    return MassProperties(m, J_b_b, r_bG_b) #p_b

end

MassProperties(sys::System) = get_mp_b(sys)


################################################################################
############################### RigidBodyData ##################################

"""
Notes:
- When get_mp_b is called on a component, the returned MassProperties instance
  must be expressed in the System's parent reference frame.
- At the root of the component hierarchy lies an AbstractComponents System,
  whose reference frame is the vehicle frame b. Therefore, the aggregate mass
  properties for the complete component hierarchy will be expressed in the
  vehicle frame, wherein aircraft dynamics and kinematics are formulated.
- Having an arbitrary origin for the vehicle frame, fixed to the vehicle body,
  instead of the current center of mass, slightly complicates the dynamics
  equations. In exchange, it allows for the aircraft's mass properties to
  change, either gradually (for example due to fuel consumption) or suddenly
  (due to a payload release) without having to worry about discontinuities in
  the kinematic state vector.
"""

#mp_b: Mass properties of the System, translated to its parent's frame

#wr_ext_b: Total external wrench contributed to the System, resolved in the
#vehicle frame

#hr_b: Intrinsic angular momentum due to the angular velocity of the rotating
#System's components with respect to the vehicle. Computed individually by
#each component relative to its center of mass and then summed

#note: for a given J_b_b, r_bG_b cannot be arbitrarily large, because moments
#of inertia are minimum at G. therefore, at some point J_G_b would become
#non-positive definite

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
################################## Actions ###################################

#for dynamics computations, we define a special reference frame G, with origin
#at the center of mass and axes parallel to the b frame axes (q_bG = 1)

#all magnitudes resolved in body axes unless otherwise noted
#of these, only mp_b and wr_ext_b are required by dynamics equations
@kwdef struct Actions
    g_G_b::SVector{3,Float64} = zeros(SVector{3}) #gravity at G
    G_G_b::SVector{3,Float64} = zeros(SVector{3}) #gravitational attraction at G
    hr_b::SVector{3,Float64} = zeros(SVector{3}) #intrinsic angular momentum
    wr_g_b::Wrench = Wrench() #gravity wrench at b
    wr_in_b::Wrench = Wrench() #inertia wrench at b
    wr_ext_b::Wrench = Wrench() #external wrench at b
    wr_net_b::Wrench = Wrench() #net wrench at b
    wr_net_G::Wrench = Wrench() #net wrench at G
end


function Actions(sys::System;
    mp_b::MassProperties = get_mp_b(sys),
    wr_ext_b::Wrench = get_wr_b(sys),
    hr_b::SVector{3,Float64} = get_hr_b(sys),
    kin_data::KinData = KinData())

    @unpack q_eb, q_nb, n_e, h_e, r_eb_e, ω_eb_b, v_eb_b = kin_data
    m = mp_b.m; J_b_b = mp_b.J; r_bG_b = mp_b.r_OG

    ω_ie_e = SVector{3, Float64}(0, 0, ω_ie) #use WGS84 constant
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    ######################## Gravity Wrench ####################################

    #compute geographic position of frame G
    r_bG_e = q_eb(r_bG_b)
    r_eG_e = r_eb_e + r_bG_e
    G = Cartesian(r_eG_e)

    #to compute gravity vector at G, create an auxiliary frame c with origin at
    #G and axes parallel to the local NED frame
    q_ec = ltf(G)
    q_bc = q_eb' ∘ q_ec
    g_G_c = SVector{3,Float64}(0, 0, gravity(G))
    g_G_b = q_bc(g_G_c)

    r_eG_b = q_eb'(r_eG_e)
    G_G_b = g_G_b + ω_ie_b × (ω_ie_b × r_eG_b)

    #the resultant from gravity on the vehicle at its center of gravity G
    #consists of gravity force, plus a null torque
    F_G_b = m * g_G_b
    M_G_b = zeros(SVector{3})
    wr_g_G = Wrench(F = F_G_b, M = M_G_b)

    #define the transform from b to G as a pure translation
    t_bG = FrameTransform(r = r_bG_b)

    #translate gravity wrench to Ob
    wr_g_b = t_bG(wr_g_G)

    ########################## Inertia Wrench ##################################

    #angular momentum of the vehicle as a rigid body (excluding intrinsic
    #angular momentum due to rotating components)
    h_rbd_b = J_b_b * ω_ib_b

    #total angular momentum
    h_all_b = h_rbd_b + SVector{3,Float64}(hr_b)

    #inertia terms (exact version... certainly overkill, but very cheap)
    a_1_b = (ω_eb_b + 2 * ω_ie_b) × v_eb_b
    F_in_b = -m * (a_1_b + ω_ib_b × (ω_ib_b × r_bG_b) + r_bG_b × (ω_eb_b × ω_ie_b ))
    M_in_b = - ( J_b_b * (ω_ie_b × ω_eb_b) + ω_ib_b × h_all_b + m * r_bG_b × a_1_b)

    wr_in_b = Wrench(F = F_in_b, M = M_in_b)

    wr_net_b = wr_ext_b + wr_g_b + wr_in_b

    #define translation from G to b
    t_Gb = t_bG'
    wr_net_G = t_Gb(wr_net_b)

    Actions(; g_G_b, G_G_b, hr_b, wr_g_b, wr_in_b, wr_ext_b, wr_net_b, wr_net_G)

end


################################################################################
############################### Dynamics #######################################

#all magnitudes resolved in body axes unless otherwise noted
@kwdef struct Accelerations
    α_eb_b::SVector{3,Float64} = zeros(SVector{3}) #ECEF-to-body angular acceleration
    α_ib_b::SVector{3,Float64} = zeros(SVector{3}) #ECI-to-body angular acceleration
    v̇_eb_b::SVector{3,Float64} = zeros(SVector{3}) #time derivative of ECEF-relative velocity
    a_eb_b::SVector{3,Float64} = zeros(SVector{3}) #ECEF-relative acceleration of b
    a_eb_n::SVector{3,Float64} = zeros(SVector{3}) #ECEF-relative acceleration of b, NED axes
    a_ib_b::SVector{3,Float64} = zeros(SVector{3}) #ECI-relative acceleration of b
    a_iG_b::SVector{3,Float64} = zeros(SVector{3}) #ECI-relative acceleration of G
    f_G_b::SVector{3,Float64} = zeros(SVector{3}) #specific force at G
end

struct RigidBodyDynamics <: SystemDefinition end

Systems.X(::RigidBodyDynamics) = zero(Kinematics.XVelTemplate)
Systems.Y(::RigidBodyDynamics) = Accelerations()

Accelerations(sys::System{<:RigidBodyDynamics}) = sys.y

function Systems.f_ode!(sys::System{RigidBodyDynamics}, mp_b::MassProperties,
                        kin_data::KinData, actions::Actions)

    @unpack q_eb, q_nb, n_e, h_e, r_eb_e, ω_eb_b, v_eb_b = kin_data
    @unpack wr_net_b, G_G_b = actions
    @unpack ẋ = sys

    ########################### Dynamic Equations ##############################

    m = mp_b.m; J_b_b = mp_b.J; r_bG_b = mp_b.r_OG

    r_bG_b_sk = Attitude.v2skew(r_bG_b)
    A11 = J_b_b
    A12 = m * r_bG_b_sk
    A21 = -m * r_bG_b_sk
    A22 = m * SMatrix{3,3,Float64}(I)

    A = vcat(hcat(A11, A12), hcat(A21, A22)) #mass matrix
    b = SVector{6}(vcat(wr_net_b.M, wr_net_b.F)) #forcing vector

    ẋ .= A\b

    ########################## Additional Outputs ##############################

    ω_ie_e = SVector{3, Float64}(0, 0, ω_ie) #use WGS84 constant
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    #angular accelerations
    α_eb_b = SVector{3}(ẋ.ω_eb_b) #α_eb_b == ω_eb_b_dot
    α_ib_b = α_eb_b - ω_eb_b × ω_ie_b

    #linear accelerations at b
    v̇_eb_b = SVector{3}(ẋ.v_eb_b)
    r_eb_b = q_eb'(r_eb_e)
    a_eb_b = v̇_eb_b + ω_eb_b × v_eb_b
    a_eb_n = q_nb(a_eb_b)
    a_ib_b = v̇_eb_b + (ω_eb_b + 2ω_ie_b) × v_eb_b + ω_ie_b × (ω_ie_b × r_eb_b)

    #linear acceleration and specific force at G
    a_iG_b = a_ib_b + ω_ib_b × (ω_ib_b × r_bG_b) + α_ib_b × r_bG_b
    f_G_b = a_iG_b - G_G_b

    sys.y = Accelerations(; α_eb_b, α_ib_b, v̇_eb_b, a_eb_b, a_eb_n,
        a_ib_b, a_iG_b, f_G_b)

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

function Plotting.make_plots(ts::TimeSeries{<:Actions}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    pd[:wr_g_b] = plot(ts.wr_g_b;
        plot_title = "Gravity Wrench at Ob [Vehicle Axes]",
        wr_source = "g", wr_frame = "b",
        kwargs...)

    pd[:wr_in_b] = plot(ts.wr_in_b;
        plot_title = "Inertia Wrench at Ob [Vehicle Axes]",
        wr_source = "in", wr_frame = "b",
        kwargs...)

    pd[:wr_ext_b] = plot(ts.wr_ext_b;
        plot_title = "External Wrench at Ob [Vehicle Axes]",
        wr_source = "ext", wr_frame = "b",
        kwargs...)

    pd[:wr_net_b] = plot(ts.wr_net_b;
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

    return pd

end

function Plotting.make_plots(ts::TimeSeries{<:Accelerations}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    pd[:α_eb_b] = plot(ts.α_eb_b;
        plot_title = "Angular Acceleration (Vehicle/ECEF) [Vehicle Axes]",
        ylabel = hcat(
            L"$\alpha_{eb}^{x_b} \ (rad/s^2)$",
            L"$\alpha_{eb}^{y_b} \ (rad/s^2)$",
            L"$\alpha_{eb}^{z_b} \ (rad/s^2)$"),
        ts_split = :h,
        kwargs...)

    pd[:a_eb_b] = plot(ts.a_eb_b;
        plot_title = "Linear Acceleration (Vehicle/ECEF) [Vehicle Axes]",
        ylabel = hcat(
            L"$a_{eb}^{x_b} \ (m/s^{2})$",
            L"$a_{eb}^{y_b} \ (m/s^{2})$",
            L"$a_{eb}^{z_b} \ (m/s^{2})$"),
        ts_split = :h,
        kwargs...)

    pd[:a_eb_n] = plot(ts.a_eb_n;
        plot_title = "Linear Acceleration (Vehicle/ECEF) [NED Axes]",
        ylabel = hcat(
            L"$a_{eb}^{N} \ (m/s^{2})$",
            L"$a_{eb}^{E} \ (m/s^{2})$",
            L"$a_{eb}^{D} \ (m/s^{2})$"),
        ts_split = :h, link = :none,
        kwargs...)

    pd[:f_G_b] = plot(TimeSeries(ts._t, ts.f_G_b._data / g₀);
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

function GUI.draw(mp_b::MassProperties, p_open::Ref{Bool} = Ref(true),
                    label::String = "Mass Properties")

    m = mp_b.m; r_bG_b = mp_b.r_OG

    #define frame transform from G to b
    t_Gb = FrameTransform(r = -r_bG_b)

    #translate mass properties to G to get inertia tensor at G
    mp_G = t_Gb(mp_b)

    CImGui.Begin(label, p_open)

        CImGui.Text(@sprintf("Mass: %.3f kg", m))
        GUI.draw(r_bG_b, "Center of Mass Position", "m")

        if CImGui.TreeNode("Inertia Tensor [Origin]")
            CImGui.Text(@sprintf("XX: %.3f kg*m2", mp_b.J[1,1]))
            CImGui.Text(@sprintf("YY: %.3f kg*m2", mp_b.J[2,2]))
            CImGui.Text(@sprintf("ZZ: %.3f kg*m2", mp_b.J[3,3]))
            CImGui.Text(@sprintf("XY: %.3f kg*m2", mp_b.J[1,2]))
            CImGui.Text(@sprintf("XZ: %.3f kg*m2", mp_b.J[1,3]))
            CImGui.Text(@sprintf("YZ: %.3f kg*m2", mp_b.J[2,3]))
            CImGui.TreePop()
        end

        if CImGui.TreeNode("Inertia Tensor [CoM]")
            CImGui.Text(@sprintf("XX: %.3f kg*m2", mp_G.J[1,1]))
            CImGui.Text(@sprintf("YY: %.3f kg*m2", mp_G.J[2,2]))
            CImGui.Text(@sprintf("ZZ: %.3f kg*m2", mp_G.J[3,3]))
            CImGui.Text(@sprintf("XY: %.3f kg*m2", mp_G.J[1,2]))
            CImGui.Text(@sprintf("XZ: %.3f kg*m2", mp_G.J[1,3]))
            CImGui.Text(@sprintf("YZ: %.3f kg*m2", mp_G.J[2,3]))
            CImGui.TreePop()
        end

    CImGui.End()

end

function GUI.draw(wr::Wrench, label::String)

    @unpack F, M = wr

    if CImGui.TreeNode(label)
        GUI.draw(F, "Force", "N")
        GUI.draw(M, "Torque", "N*m")
        CImGui.TreePop()
    end

end


function GUI.draw(dyn::Actions, p_open::Ref{Bool} = Ref(true),
                    label::String = "Actions")

    @unpack g_G_b, G_G_b, hr_b, wr_g_b, wr_in_b, wr_ext_b, wr_net_b, wr_net_G = dyn

    CImGui.Begin(label, p_open)

        GUI.draw(g_G_b, "Gravity (Vehicle Frame)", "m/(s^2)")
        GUI.draw(G_G_b, "Gravitation (CoM)", "m/(s^2)")
        GUI.draw(wr_g_b, "Gravity Wrench (Vehicle Frame)")
        GUI.draw(wr_in_b, "Inertia Wrench (Vehicle Frame)")
        GUI.draw(wr_ext_b, "External Wrench (Vehicle Frame)")
        GUI.draw(wr_net_b, "Net Wrench (Vehicle Frame)")
        GUI.draw(wr_net_G, "Net Wrench (CoM)")
        GUI.draw(hr_b, "Intrinsic Angular Momentum", "kg*(m^2)/s")

    CImGui.End()

end

function GUI.draw(dyn::Accelerations, p_open::Ref{Bool} = Ref(true),
                    label::String = "Accelerations")


    @unpack α_eb_b, α_ib_b, v̇_eb_b, a_eb_b, a_eb_n, a_ib_b, a_iG_b, f_G_b = dyn

    CImGui.Begin(label, p_open)

    GUI.draw(α_eb_b, "Angular Acceleration (Vehicle / ECEF) [Vehicle]", "rad/(s^2)")
    GUI.draw(α_ib_b, "Angular Acceleration (Vehicle / ECI) [Vehicle]", "rad/(s^2)")
    GUI.draw(v̇_eb_b, "Velocity Time-Derivative (Vehicle / ECEF) [Vehicle]", "rad/(s^2)")
    GUI.draw(a_eb_b, "Linear Acceleration (Vehicle / ECEF) [Vehicle]", "m/(s^2)")
    GUI.draw(a_eb_n, "Linear Acceleration (Vehicle / ECEF) [NED]", "m/(s^2)")
    GUI.draw(a_ib_b, "Linear Acceleration (Vehicle / ECI) [Vehicle]", "m/(s^2)")
    GUI.draw(a_iG_b, "Linear Acceleration (CoM / ECI) [Vehicle]", "m/(s^2)")
    GUI.draw(f_G_b/g₀, "Specific Force (CoM) [Vehicle]", "g")

    CImGui.End()

end


end #module