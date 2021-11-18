module Dynamics

using StaticArrays
using LinearAlgebra
using UnPack

using Flight.Geodesy
using Flight.Attitude
using Flight.Kinematics
using Flight.ModelingTools

using Flight.Plotting
import Flight.Plotting: plots

export DynData, MassData, FrameTransform, Wrench
export translate, f_dyn!



############################# FrameTransform ###############################

"""
Specifies a reference frame `fc(Oc, Ɛc)` relative to another `fb(Ob, Ɛb)`

Frame `fc(Oc, Ɛc)` is specified by:
- The position vector from fb's origin Ob to fc's origin Oc, projected on fb's
  axes εb (`r_ObOc_b`)
- The rotation quaternion from fb's axes εb to fc's axes εc (`q_bc`)
"""
Base.@kwdef struct FrameTransform
    r::SVector{3,Float64} = zeros(SVector{3})
    q::RQuat = RQuat()
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
    translate(t_bc::FrameTransform, wr_c::Wrench)

Translate a Wrench from one reference frame to another.

If `t_bc` is a `FrameTransform` specifying frame fc(Oc, εc) relative to fb(Ob, εb),
and `wr_c` is a `Wrench` defined on fc, then `wr_b = translate(t_bc,
wr_c)` is the equivalent `Wrench` defined on fb.

An alternative function call notation is provided:
`t_bc(wr_c) == translate(t_bc, wr_c)`
"""
(t_bc::FrameTransform)(wr_c::Wrench) = translate(t_bc, wr_c)

function translate(t_bc::FrameTransform, wr_c::Wrench)

    F_Oc_c = wr_c.F
    M_Oc_c = wr_c.M

    #project onto airframe axes
    F_Oc_b = t_bc.q(F_Oc_c)
    M_Oc_b = t_bc.q(M_Oc_c)

    #translate to airframe origin
    F_Ob_b = F_Oc_b
    M_Ob_b = M_Oc_b + t_bc.r × F_Oc_b

    return Wrench(F = F_Ob_b, M = M_Ob_b) #wr_b

end


############################## MassData #################################

"""
Groups the mass properties required by the aircraft's dynamic equations

These are:
- `m`: Total aircraft mass, including that of the rotating elements.
- `J_Ob_b`: Overall aircraft inertia tensor with respect to the airframe origin
  Ob, projected on the airframe axes εb.
- `r_ObG_b`: Position vector from the airframe origin Ob to the aircraft's
  center of mass G, projected on the airframe axes εb.

Notes:
- The inertia tensor `J_Ob_b` must include the contributions of any rotating
  elements on the aircraft. As long as rotating elements have axial symmetry
  around their axes of rotation, their contributions to `J_Ob_b` will be
  strictly constant.
- Aircraft dynamics and kinematics are formulated on the airframe origin Ob
  instead of the aircraft's center of mass G. This allows for any of the
  aircraft's mass properties to change, either gradually (for example, due to
  fuel consumption) or suddenly (due to a payload release), without having to
  worry about discontinuities in the kinematic state vector.
"""
Base.@kwdef struct MassData
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end


################## Dynamic Equations and helper functions ####################


"""
    inertia_wrench(mass::MassData, vel::VelData, hr_b::AbstractVector{<:Real})

Compute the equivalent `Wrench` arising from inertia terms in the dynamic
equations

The resulting `Wrench` is defined on the airframe's reference frame.

# Arguments:
- `mass::MassData`: Current aircraft mass properties
- `vel::VelData`: Velocity outputs
- `hr_b::AbstractVector{<:Real}`: Additional angular momentum due to the
  angular velocity of any rotating elements with respect to the airframe,
  projected on the airframe axes

"""
function inertia_wrench(mass::MassData, vel::VelData, hr_b::AbstractVector{<:Real})

    @unpack m, J_Ob_b, r_ObG_b = mass
    @unpack ω_ie_b, ω_eb_b, ω_ib_b, v_eOb_b = vel

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

function gravity_wrench(mass::MassData, pos::PosData)

    #gravity can be viewed as an entity acting on a local frame with its origin
    #at G and its axes aligned with the local tangent frame

    #strictly, the gravity vector should be evaluated at G, with its direction
    #given by the z-axis of LTF(G). however, since g(G) ≈ g(Ob) and LTF(G) ≈
    #LTF(Ob), we can instead evaluate g at Ob, assuming its direction given by
    #LTF(Ob), and then apply it at G.
    @unpack n_e, h_e, q_nb = pos

    Ob = Geographic(n_e, h_e)
    g_G_n = g_Ob_n = g_n(Ob)

    #the resultant consists of the gravity force acting on G along the local
    #vertical and a null torque
    F_G_n = mass.m * g_G_n
    M_G_n = zeros(SVector{3})
    wr_G_n = Wrench(F = F_G_n, M = M_G_n)

    #with the previous assumption, the transformation from body frame to local
    #gravity frame is given by the translation r_ObG_b and the (passive)
    #rotation from b to LTF(Ob) (instead of LTF(G)), which is given by pos.l_b'
    wr_c = wr_G_n
    t_bc = FrameTransform(r = mass.r_ObG_b, q = q_nb')
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


function f_dyn!(ẋ_vel::VelX, kin::KinData, mass::MassData,
    wr_ext_b::Wrench, hr_b::AbstractVector{<:Real})

    #wr_ext_b: Total external wrench on the airframe

    #hr_b: Additional angular momentum due to the angular velocity of the
    #rotating aircraft components with respect to the airframe. these are
    #computed individually by each component relative to its center of mass and
    #then summed

    #wr_ext_b and hr_b, as well as mass data, are produced by aircraft
    #components, so they must be computed by the aircraft's x_dot method. and,
    #since kin is needed by those components, it must be called from the
    #aircraft's kinematic state vector

    #clearly, r_ObG_b cannot be arbitrarily large, because J_Ob_b is larger than
    #J_G_b (Steiner). therefore, at some point J_G_b would become zero (or at
    #least singular)!

    @unpack m, J_Ob_b, r_ObG_b = mass
    @unpack q_eb, q_nb, n_e, h_e = kin.pos
    @unpack ω_eb_b, ω_ie_b, v_eOb_b = kin.vel

    r_ObG_b_sk = Attitude.skew(r_ObG_b)
    A11 = J_Ob_b
    A12 = m * r_ObG_b_sk
    A21 = -m * r_ObG_b_sk
    A22 = m * SMatrix{3,3,Float64}(I)

    A = vcat(hcat(A11, A12), hcat(A21, A22))

    wr_g_b = gravity_wrench(mass, kin.pos)
    wr_in_b = inertia_wrench(mass, kin.vel, hr_b)
    wr_b = wr_ext_b + wr_g_b + wr_in_b
    b = SVector{6}(vcat(wr_b.M, wr_b.F))

    # update ẋ_vel
    ẋ_vel .= A\b

    #compute outputs
    Ob = Geographic(n_e, h_e)
    r_eOb_e = CartECEF(Ob)[:]
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


######################## Plots ############################

@recipe function plot_wrench(th::TimeHistory{<:AbstractVector{<:Wrench}};
                             wr_frame = "", wr_source = "")

    @unpack F, M = StructArray(th.data)

    layout := (1, 2)
    seriestype --> :path

    @series begin
        subplot := 1
        title --> "Force"
        yguide --> L"$F_{O%$wr_frame \ (%$wr_source)}^{%$wr_frame} \ (N)$"
        # yguide --> L"$F \ (N)$"
        th_split --> :none
        TimeHistory(th.t, F)
    end

    @series begin
        subplot := 2
        title --> "Torque"
        yguide --> L"$M_{O%$wr_frame \ (%$wr_source)}^{%$wr_frame} \ (N \ m)$"
        th_split --> :none
        TimeHistory(th.t, M)
    end

end


function plots(t, data::AbstractVector{<:DynData}; mode, save_path, kwargs...)

    sa = StructArray(data)
    plots(t, sa.input; mode, save_path, kwargs...)
    plots(t, sa.output; mode, save_path, kwargs...)

end


function plots(t, data::AbstractVector{<:DynDataIn}; mode, save_path, kwargs...)

    @unpack wr_g_b, wr_in_b, wr_ext_b, hr_b = StructArray(data)

    pd = Dict{String, Plots.Plot}()

    pd["01_wr_g_b"] = thplot(t, wr_g_b;
        plot_title = "Gravity Wrench [Airframe]",
        wr_source = "g", wr_frame = "b",
        kwargs...)

    pd["02_wr_in_b"] = thplot(t, wr_in_b;
        plot_title = "Inertia Wrench [Airframe]",
        wr_source = "in", wr_frame = "b",
        kwargs...)

    pd["03_wr_ext_b"] = thplot(t, wr_ext_b;
        plot_title = "External Wrench [Airframe]",
        wr_source = "ext", wr_frame = "b",
        kwargs...)

    pd["04_hr_b"] = thplot(t, hr_b;
        plot_title = "Angular Momentum from Rotating Components [Airframe]",
        ylabel = hcat(
            L"$h_{Ob \ (r)}^{x_b} \ (kg \ m^2 / s)$",
            L"$h_{Ob \ (r)}^{y_b} \ (kg \ m^2 / s)$",
            L"$h_{Ob \ (r)}^{z_b} \ (kg \ m^2 / s)$"),
        th_split = :h, link = :none,
        kwargs...)

    save_plots(pd; save_path)

end


function plots(t, data::AbstractVector{<:DynDataOut}; mode, save_path, kwargs...)

    @unpack α_eb_b, a_eOb_b, a_eOb_n, f_Ob_b = StructArray(data)

    #standard gravity for specific force normalization
    g₀ = 9.80665

    pd = Dict{String, Plots.Plot}()

    pd["05_α_eb_b"] = thplot(t, α_eb_b;
        plot_title = "Angular Acceleration (Airframe/ECEF) [Airframe]",
        ylabel = hcat(
            L"$\alpha_{eb}^{x_b} \ (rad/s^2)$",
            L"$\alpha_{eb}^{y_b} \ (rad/s^2)$",
            L"$\alpha_{eb}^{z_b} \ (rad/s^2)$"),
        th_split = :h,
        kwargs...)

    pd["06_a_eOb_b"] = thplot(t, a_eOb_b;
        plot_title = "Linear Acceleration (Airframe/ECEF) [Airframe]",
        ylabel = hcat(
            L"$a_{eb}^{x_b} \ (m/s^{2})$",
            L"$a_{eb}^{y_b} \ (m/s^{2})$",
            L"$a_{eb}^{z_b} \ (m/s^{2})$"),
        th_split = :h,
        kwargs...)

    pd["07_a_eOb_n"] = thplot(t, a_eOb_n;
        plot_title = "Linear Acceleration (Airframe/ECEF) [NED]",
        ylabel = hcat(
            L"$a_{eb}^{N} \ (m/s^{2})$",
            L"$a_{eb}^{E} \ (m/s^{2})$",
            L"$a_{eb}^{D} \ (m/s^{2})$"),
        th_split = :h, link = :none,
        kwargs...)

    pd["08_f_Ob_b"] = thplot(t, f_Ob_b / g₀;
        plot_title = "Specific Force [Airframe]",
        ylabel = hcat(
            L"$f_{Ob}^{x_b} \ (g)$",
            L"$f_{Ob}^{y_b} \ (g)$",
            L"$f_{Ob}^{z_b} \ (g)$"),
        th_split = :h,
        kwargs...)

    save_plots(pd; save_path)

end

end #module