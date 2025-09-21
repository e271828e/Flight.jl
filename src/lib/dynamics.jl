module Dynamics

using StaticArrays, LinearAlgebra, UnPack
using Plots: Plots, @recipe, @series, plot
using LaTeXStrings, DataStructures

using Flight.FlightCore

using ..Attitude
using ..Geodesy
using ..Kinematics

export FrameTransform
export Wrench
export AbstractMassDistribution, PointDistribution, RigidBodyDistribution, MassProperties
export get_mp_b, get_wr_b, get_hr_b
export VehicleDynamics, DynamicsData

#standard gravity for specific force normalization
const g₀ = 9.80665


################################################################################
############################# FrameTransform ###################################

"""
Specify a reference frame ``c`` with respect to another ``b``.

``c`` is defined by:
- ``r_{bc}^{b}``: Position vector from ``b``'s origin ``O_b`` to ``c``'s origin
  ``O_c``, projected on ``b``'s axes ``Ɛ_b``
- ``q^b_c``: Rotation quaternion from ``b``'s axes ``Ɛ_b`` to ``c``'s axes
  ``Ɛ_c``

`FrameTransform` objects are callable; the syntax `t_bc(x_c)` is equivalent to
`translate(t_bc, x_c)`.
"""
@kwdef struct FrameTransform
    r::SVector{3,Float64} = zeros(SVector{3})
    q::RQuat = RQuat()
end


"""
    translate(t_bc::FrameTransform, r_cP_c::AbstractVector{<:Real})

Translate a 3D position vector from one reference frame to another.

Given the `FrameTransform` `t_bc` from frame ``b`` to frame ``c`` and the
position vector of some point ``P`` in ``c`` (``r_{cP}^{c}``), return the
position vector of ``P`` in ``b`` (``r_{bP}^{b}``).
"""
function translate(t_bc::FrameTransform, r_cP_c::AbstractVector{<:Real})

    r_bc_b = t_bc.r; q_bc = t_bc.q

    r_bP_b = r_bc_b + q_bc(SVector{3,Float64}(r_cP_c))
    return r_bP_b
end


"""
    (t_bc::FrameTransform)(x_c::Any)

Convenience method, equivalent to `translate(t_bc, x_c)`.
"""
(t_bc::FrameTransform)(x_c::Any) = translate(t_bc, x_c)


"""
    Base.:adjoint(t_bc::FrameTransform)

Given `t_bc`, return its reciprocal `t_cb`.
"""
function Base.:adjoint(t_bc::FrameTransform)

    r_bc_b = t_bc.r; q_bc = t_bc.q

    q_cb = q_bc'
    r_cb_c = q_cb(-r_bc_b)
    t_cb = FrameTransform(r_cb_c, q_cb)

    return t_cb
end


"""
    Base.:∘(t_bc::FrameTransform, t_cd::FrameTransform)

Concatenate `FrameTransform`s `t_bc` and `t_cd` and return the resulting `t_bd`.
"""
function Base.:∘(t_bc::FrameTransform, t_cd::FrameTransform)

    r_bc_b = t_bc.r; q_bc = t_bc.q
    r_cd_c = t_cd.r; q_cd = t_cd.q

    r_bd_b = r_bc_b + q_bc(r_cd_c)
    q_bd = q_bc ∘ q_cd

    t_bd = FrameTransform(r_bd_b, q_bd)

    return t_bd

end

################################################################################
################################## Wrench ######################################

"""
Force and torque combination defined in a specific reference frame.

A `Wrench` is defined in frame ``c`` when it is applied at ``c``'s origin
``O_c`` and projected on its axes ``Ɛ_c``
"""
@kwdef struct Wrench
    F::SVector{3,Float64} = zeros(SVector{3})
    τ::SVector{3,Float64} = zeros(SVector{3})
 end


"""
    Base.:+(wr1::Wrench, wr2::Wrench)

Add two compatible `Wrench` instances.

`Wrench` addition should only be performed between compatible `Wrench`
instances, i.e., those defined in the same reference frame
"""
Base.:+(wr1::Wrench, wr2::Wrench) = Wrench(F = wr1.F + wr2.F, τ = wr1.τ + wr2.τ)


"""
    translate(t_bc::FrameTransform, wr_c::Wrench)

Translate a `Wrench` from one reference frame to another.

If `t_bc` is a `FrameTransform` specifying frame fc(Oc, εc) relative to fb(Ob, εb),
and `wr_c` is a `Wrench` defined on fc, then `wr_b = translate(t_bc,
wr_c)` is the equivalent `Wrench` defined on fb.
"""

function translate(t_bc::FrameTransform, wr_c::Wrench)

    F_c_c = wr_c.F
    τ_c_c = wr_c.τ

    #project onto b axes
    F_c_b = t_bc.q(F_c_c)
    τ_c_b = t_bc.q(τ_c_c)

    #translate to b origin
    F_b_b = F_c_b
    τ_b_b = τ_c_b + t_bc.r × F_c_b

    return Wrench(F = F_b_b, τ = τ_b_b) #wr_b

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
    J::SMatrix{3,3,Float64,9} = SMatrix{3,3,Float64}(I)
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
    translate(t_bc::FrameTransform, mp_c::MassProperties)

Translate a `MassProperties` instance from its current frame fc to another fb.

If `t_bc` is a `FrameTransform` specifying frame fc(Oc, εc) relative to fb(Ob,
εb), and `mp_c` is a `MassProperties` defined with respect to fc, then `mp_b =
translate(t_bc, p_c)` is the equivalent `MassProperties` defined with respect to
fb.
"""
function translate(t_bc::FrameTransform, mp_c::MassProperties)

    #note: for a given J_b, r_bc_b cannot be arbitrarily large, because moments
    #of inertia are minimum at G. therefore, at some point J_c would become
    #non-positive definite

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


############################ Input Gathering ###################################

#generated functions are required because type inference does not work
#throughout the whole Systems hierarchy

# mp_b: Mass properties translated to body frame.
# default implementation tries to compute the aggregate mass properties for all
# its submodels. override if possible to reduce compilation time
@inline @generated function (get_mp_b(mdl::Model{T, X, Y, U, S, P, B})
    where {T, X, Y, U, S, P, B})

    ex = Expr(:block)

    if isempty(fieldnames(B))
        push!(ex.args,
            :(error("Model{$(T)} is a leaf Model "*
                "but it does not extend the get_mp_b method")))
    else
        push!(ex.args, :(p = MassProperties()))
        for label in fieldnames(B)
            push!(ex.args,
                :(p += get_mp_b(mdl.submodels[$(QuoteNode(label))])))
        end
    end
    return ex

end


#wr_b: Total external wrench contributed, resolved in body frame.
#default implementation tries to sum the angular momentum from its individual
#components. override if possible to reduce compilation time
@inline @generated function (get_hr_b(mdl::Model{T, X, Y, U, S, P, B})
    where {T, X, Y, U, S, P, B})

    # Core.print("Generated function called")
    ex = Expr(:block)

    if isempty(fieldnames(B))
        push!(ex.args,
            :(error("Model{$(T)} is a leaf Model "*
                "but it does not extend the get_hr_b method")))
    else
        push!(ex.args, :(h = SVector(0., 0., 0.))) #initialize
        for label in fieldnames(B)
            push!(ex.args,
                :(h += get_hr_b(mdl.submodels[$(QuoteNode(label))])))
        end
    end

    return ex

end

#ho_b: Internal angular momentum relative to the airframe due to rotating
#components, resolved in the vehicle frame.
#default implementation tries to sum all the Wrenches from its individual
#components. override if possible to reduce compilation time
@inline @generated function (get_wr_b(mdl::Model{T, X, Y, U, S, P, B})
    where {T, X, Y, U, S, P, B})

    # Core.print("Generated function called")

    ex = Expr(:block)

    if isempty(fieldnames(B))
        push!(ex.args,
            :(error("Model{$(T)} is a leaf Model, "*
                "but it does not extend the get_wr_b method")))
    else
        push!(ex.args, :(wr = Wrench())) #initialize a zero wrench
        for label in fieldnames(B)
            push!(ex.args,
                :(wr += get_wr_b(mdl.submodels[$(QuoteNode(label))])))
        end
    end

    return ex

end

# MassProperties(mdl::Model) = get_mp_b(mdl)

################################################################################
########################### VehicleDynamics ####################################

struct VehicleDynamics <: ModelDefinition end

@kwdef mutable struct DynamicsU
    mp_Σ_b::MassProperties = MassProperties(RigidBodyDistribution()) #total mass properties, Body frame
    wr_Σ_b::Wrench = Wrench() #total external wrench, CoM frame
    ho_Σ_b::SVector{3,Float64} = zeros(SVector{3}) #total internal angular momentum, Body axes
    q_eb::RQuat = KinData().q_eb #Body frame attitude with respect to ECEF frame
    r_eb_e::SVector{3, Float64} = KinData().r_eb_e #Cartesian position of Body frame with respecto ECEF frame, ECEF coordinates
end

@kwdef struct DynamicsData
    wr_Σ_c::Wrench = Wrench() #total external wrench at CoM
    wr_Σ_b::Wrench = Wrench() #total external wrench at Body
    mp_Σ_c::MassProperties = MassProperties() #mass properties at CoM
    mp_Σ_b::MassProperties = MassProperties() #mass properties at Body
    ho_Σ_b::SVector{3,Float64} = zeros(SVector{3}) #internal angular momentum, Body axes
    ω̇_ec_c::SVector{3,Float64} = zeros(SVector{3}) #ECEF-to-CoM angular acceleration
    v̇_ec_c::SVector{3,Float64} = zeros(SVector{3}) #ECEF-to-CoM velocity derivative
    a_ec_c::SVector{3,Float64} = zeros(SVector{3}) #ECEF-to-CoM acceleration
    a_ic_c::SVector{3,Float64} = zeros(SVector{3}) #ECI-to-CoM acceleration
    g_c_c::SVector{3,Float64} = zeros(SVector{3}) #gravity at CoM
    γ_c_c::SVector{3,Float64} = zeros(SVector{3}) #gravitational attraction at CoM
    f_c_c::SVector{3,Float64} = zeros(SVector{3}) #specific force at CoM
    ω̇_eb_b::SVector{3,Float64} = zeros(SVector{3}) #ECEF-to-Body angular acceleration
    v̇_eb_b::SVector{3,Float64} = zeros(SVector{3}) #ECEF-to-Body velocity derivative
    α_ib_b::SVector{3,Float64} = zeros(SVector{3}) #ECI-to-Body angular acceleration
    a_eb_b::SVector{3,Float64} = zeros(SVector{3}) #ECEF-to-Body acceleration
    a_ib_b::SVector{3,Float64} = zeros(SVector{3}) #ECI-to-Body acceleration
end

Modeling.X(::VehicleDynamics) = zero(Kinematics.XVelTemplate)
Modeling.U(::VehicleDynamics) = DynamicsU()
Modeling.Y(::VehicleDynamics) = DynamicsData()

function Modeling.f_ode!(mdl::Model{VehicleDynamics})

    @unpack x, ẋ, u = mdl
    @unpack mp_Σ_b, wr_Σ_b, ho_Σ_b, q_eb, r_eb_e = u

    ω_eb_b = SVector{3, Float64}(x.ω_eb_b) #also available as kin_data.ω_eb_b
    v_eb_b = SVector{3, Float64}(x.v_eb_b) #also available as kin_data.v_eb_b

    ω_ie_e = SVector{3, Float64}(0, 0, ω_ie) #use WGS84 constant
    ω_ie_b = q_eb'(ω_ie_e)

    #frame transform from c (CoM) to b (body)
    r_bc_b = mp_Σ_b.r_OG
    t_cb = FrameTransform(r = -r_bc_b) #pure translation


    ###################### angular & linear momentum equations #################

    #translate data to frame c
    mp_Σ_c = t_cb(mp_Σ_b)
    wr_Σ_c = t_cb(wr_Σ_b)
    ho_Σ_c = ho_Σ_b

    F_Σ_c = wr_Σ_c.F; τ_Σ_c = wr_Σ_c.τ
    m_Σ = mp_Σ_c.m; J_Σ_c = mp_Σ_c.J

    ω_ec_c = ω_eb_b
    v_ec_c = v_eb_b + ω_ec_c × r_bc_b

    ω_ie_c = ω_ie_b
    ω_ic_c = ω_ie_c + ω_ec_c

    #compute geographic position of Oc
    r_bc_e = q_eb(r_bc_b)
    r_ec_e = r_eb_e + r_bc_e
    Oc = Cartesian(r_ec_e)

    #define auxiliary local-level frame l with Ol = Oc
    q_el = ltf(Oc)
    q_be = q_eb'
    q_ce = q_be
    q_cl = q_ce ∘ q_el

    #compute gravity at c
    g_c_l = SVector{3,Float64}(0, 0, gravity(Oc)) #gravity at c, l coordinates\
    g_c_c = q_cl(g_c_l) #gravity at c, c coordinates

    #solve dynamic equations at c
    hc_Σ_c = J_Σ_c * ω_ic_c + ho_Σ_c
    ω̇_ec_c = J_Σ_c \ (τ_Σ_c - J_Σ_c * (ω_ie_c × ω_ec_c) - ω_ic_c × hc_Σ_c)
    v̇_ec_c = 1/m_Σ * F_Σ_c + g_c_c - (ω_ec_c + 2ω_ie_c) × v_ec_c

    #translate derivatives back to b
    ω̇_eb_b = ω̇_ec_c
    v̇_eb_b = v̇_ec_c - ω̇_ec_c × r_bc_b


    ########################## additional outputs ##############################

    r_ec_c = q_ce(r_ec_e)
    r_eb_b = q_be(r_eb_e)

    a_ec_c = v̇_ec_c + ω_ec_c × v_ec_c
    a_ic_c = v̇_ec_c + (ω_ec_c + 2ω_ie_c) × v_ec_c + ω_ie_c × (ω_ie_c × r_ec_c)
    γ_c_c = g_c_c + ω_ie_c × (ω_ie_c × r_ec_c)
    f_c_c = a_ic_c - γ_c_c

    α_eb_b = ω̇_eb_b
    α_ib_b = α_eb_b - ω_eb_b × ω_ie_b
    a_eb_b = v̇_eb_b + ω_eb_b × v_eb_b
    a_ib_b = v̇_eb_b + (ω_eb_b + 2ω_ie_b) × v_eb_b + ω_ie_b × (ω_ie_b × r_eb_b)


    ############################## Model update ###############################

    ẋ.ω_eb_b = ω̇_eb_b
    ẋ.v_eb_b = v̇_eb_b

    mdl.y = DynamicsData(; wr_Σ_c, wr_Σ_b, mp_Σ_c, mp_Σ_b, ho_Σ_b,
        ω̇_ec_c, v̇_ec_c, a_ec_c, a_ic_c, g_c_c, γ_c_c, f_c_c,
        ω̇_eb_b, v̇_eb_b, α_ib_b, a_eb_b, a_ib_b)

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
        yguide --> L"$\tau_{O%$wr_frame \ (%$wr_source)}^{%$wr_frame} \ (N \ m)$"
        ts_split --> :none
        ts.τ
    end

end

function Plotting.make_plots(ts::TimeSeries{<:DynamicsData}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    pd[:g_c] = plot(ts.g_c_c;
        plot_title = "Gravity at CoM [CoM Axes]",
        ylabel = hcat(
            L"$g_{c}^{x_c} \ (m / s^2)$",
            L"$g_{c}^{y_c} \ (m / s^2)$",
            L"$b_{c}^{z_c} \ (m / s^2)$"),
        ts_split = :v, link = :none,
        kwargs...)

    pd[:wr_ext_c] = plot(ts.wr_Σ_c;
        plot_title = "External Wrench at CoM [CoM Axes]",
        wr_source = "ext", wr_frame = "c",
        kwargs...)

    pd[:wr_ext_b] = plot(ts.wr_Σ_b;
        plot_title = "External Wrench at Ob [Body Axes]",
        wr_source = "ext", wr_frame = "b",
        kwargs...)

    pd[:hr_b] = plot(ts.ho_Σ_b;
        plot_title = "Internal Angular Momentum [Body Axes]",
        ylabel = hcat(
            L"$h_{Ob \ (r)}^{x_b} \ (kg \ m^2 / s)$",
            L"$h_{Ob \ (r)}^{y_b} \ (kg \ m^2 / s)$",
            L"$h_{Ob \ (r)}^{z_b} \ (kg \ m^2 / s)$"),
        ts_split = :v, link = :none,
        kwargs...)

    pd[:α_eb_b] = plot(ts.ω̇_eb_b;
        plot_title = "Angular Acceleration (Body/ECEF) [Body Axes]",
        ylabel = hcat(
            L"$\alpha_{eb}^{x_b} \ (rad/s^2)$",
            L"$\alpha_{eb}^{y_b} \ (rad/s^2)$",
            L"$\alpha_{eb}^{z_b} \ (rad/s^2)$"),
        ts_split = :v,
        kwargs...)

    pd[:a_eb_b] = plot(ts.a_eb_b;
        plot_title = "Linear Acceleration (Body/ECEF) [Body Axes]",
        ylabel = hcat(
            L"$a_{eb}^{x_b} \ (m/s^{2})$",
            L"$a_{eb}^{y_b} \ (m/s^{2})$",
            L"$a_{eb}^{z_b} \ (m/s^{2})$"),
        ts_split = :v,
        kwargs...)

    pd[:f_c_c] = plot(TimeSeries(ts._t, ts.f_c_c._data / g₀);
        plot_title = "Specific Force at CoM [CoM Axes]",
        ylabel = hcat(
            L"$f_{c}^{x_c} \ (g)$",
            L"$f_{c}^{y_c} \ (g)$",
            L"$f_{c}^{z_c} \ (g)$"),
        ts_split = :v, link = :none,
        kwargs...)

    return pd

end

################################################################################
################################# GUI ##########################################

function GUI.draw(mp::MassProperties,
                    label::String = "Mass Properties")

    m = mp.m; r_bc_b = mp.r_OG

    if CImGui.TreeNode(label)
        CImGui.Text(@sprintf("Mass: %.3f kg", m))
        GUI.draw(r_bc_b, "CoM Position", "m")

        if CImGui.TreeNode("Inertia Tensor")
            CImGui.Text(@sprintf("XX: %.3f kg*m2", mp.J[1,1]))
            CImGui.Text(@sprintf("YY: %.3f kg*m2", mp.J[2,2]))
            CImGui.Text(@sprintf("ZZ: %.3f kg*m2", mp.J[3,3]))
            CImGui.Text(@sprintf("XY: %.3f kg*m2", mp.J[1,2]))
            CImGui.Text(@sprintf("XZ: %.3f kg*m2", mp.J[1,3]))
            CImGui.Text(@sprintf("YZ: %.3f kg*m2", mp.J[2,3]))
            CImGui.TreePop()
        end
        CImGui.TreePop()
    end

end

function GUI.draw(wr::Wrench, label::String)

    @unpack F, τ = wr
    if CImGui.TreeNode(label)
        GUI.draw(F, "Force", "N")
        GUI.draw(τ, "Torque", "N*m")
        CImGui.TreePop()
    end

end

GUI.draw(dyn::Model{VehicleDynamics}) = GUI.draw(dyn.y)

function GUI.draw(dyn::DynamicsData, p_open::Ref{Bool} = Ref(true),
                    label::String = "Vehicle Dynamics")

    @unpack wr_Σ_c, wr_Σ_b, mp_Σ_c, mp_Σ_b, ho_Σ_b,
            a_ec_c, a_ic_c, g_c_c, γ_c_c, f_c_c, ω̇_eb_b, α_ib_b,
            a_eb_b, a_ib_b = dyn

    CImGui.Begin(label, p_open)

        GUI.draw(mp_Σ_b, "Mass Properties (Body)")
        GUI.draw(mp_Σ_c, "Mass Properties (CoM)")
        GUI.draw(wr_Σ_b, "External Wrench (Body)")
        GUI.draw(wr_Σ_c, "External Wrench (CoM)")
        GUI.draw(ho_Σ_b, "Internal Angular Momentum (Body)", "kg*m^2/s")
        GUI.draw(ω̇_eb_b, "Angular Acceleration (Body / ECEF) [Body Axes]", "rad/(s^2)")
        GUI.draw(α_ib_b, "Angular Acceleration (Body / ECI) [Body Axes]", "rad/(s^2)")
        GUI.draw(a_ec_c, "Linear Acceleration (CoM / ECEF) [CoM Axes]", "m/(s^2)")
        GUI.draw(a_ic_c, "Linear Acceleration (CoM / ECI) [CoM Axes]", "m/(s^2)")
        GUI.draw(a_eb_b, "Linear Acceleration (Body / ECEF) [Body Axes]", "m/(s^2)")
        GUI.draw(a_ib_b, "Linear Acceleration (Body / ECI) [Body Axes]", "m/(s^2)")
        GUI.draw(g_c_c, "Gravity (CoM) [CoM Axes]", "m/s^2")
        GUI.draw(γ_c_c, "Gravitational Attraction (CoM) [CoM Axes]", "m/s^2")
        GUI.draw(f_c_c/g₀, "Specific Force (CoM) [CoM Axes]", "g")

    CImGui.End()

end


end #module