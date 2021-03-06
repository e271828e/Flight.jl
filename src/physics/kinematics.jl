module Kinematics

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack
using Plots

using Flight.Utils
using Flight.Systems
using Flight.Plotting

using Flight.Attitude
using Flight.Geodesy

import Flight.Systems: init, f_cont!, f_disc!
import Flight.Plotting: make_plots

export AbstractKinematics, ECEF, LTF, NED

########################## AbstractKinematics #############################
##############################################################################

abstract type AbstractKinematics <: SystemDescriptor end

Base.@kwdef struct Initializer
    ω_lb_b::SVector{3, Float64} = zeros(SVector{3})
    v_eOb_n::SVector{3, Float64} = zeros(SVector{3})
    q_nb::RQuat = RQuat()
    Ob::GeographicLocation{NVector,Ellipsoidal} = GeographicLocation()
    Δx::Float64 = 0.0
    Δy::Float64 = 0.0
end

struct Common
    q_nb::RQuat
    q_eb::RQuat
    n_e::NVector
    h_e::Altitude{Ellipsoidal}
    h_o::Altitude{Orthometric}
    Δxy::SVector{2,Float64}
    ω_lb_b::SVector{3,Float64}
    ω_eb_b::SVector{3,Float64}
    ω_ie_b::SVector{3,Float64}
    ω_ib_b::SVector{3,Float64}
    v_eOb_b::SVector{3,Float64}
    v_eOb_n::SVector{3,Float64}
end

Base.@kwdef struct ECEFSpecific
    q_en::RQuat
    ω_el_n::SVector{3,Float64}
end

Base.@kwdef struct LTFSpecific
    q_lb::RQuat
    q_el::RQuat
    ω_el_l::SVector{3,Float64}
end

Base.@kwdef struct NEDSpecific
    e_nb::REuler
    ϕ_λ::LatLon
    ω_nb_b::SVector{3,Float64}
    ω_en_n::SVector{3,Float64}
end

struct Y{S <: Union{ECEFSpecific, LTFSpecific, NEDSpecific}}
    common::Common
    specific::S
end

Base.getproperty(y::Y, s::Symbol) = getproperty(y, Val(s))

@generated function Base.getproperty(y::Y{T}, ::Val{S}) where {T, S}
    if S === :common || S === :specific
        return :(getfield(y, $(QuoteNode(S))))
    elseif S ∈ fieldnames(Common)
        return :(getfield(getfield(y, :common), $(QuoteNode(S))))
    elseif S ∈ fieldnames(T)
        return :(getfield(getfield(y, :specific), $(QuoteNode(S))))
    else
        error("$(typeof(Y)) has no property $S")
    end
end

Common(ic::Initializer = Initializer()) = Y(LTF(), ic).common

function init(::SystemX, kin::AbstractKinematics, ic::Initializer = Initializer())
    x = similar(x_template(kin))
    init!(x, ic)
    return x
end

function init(::SystemY, kin::AbstractKinematics, ic::Initializer = Initializer())
    return Y(kin, ic)
end

function Y(kin::AbstractKinematics, ic::Initializer = Initializer())
    x = init(SystemX(), kin, ic)
    return Y(x)
end

init!(sys::System{<:AbstractKinematics}, ic::Initializer = Initializer()) = init!(sys.x, ic)

#for dispatching
const XVelTemplate = ComponentVector(ω_eb_b = zeros(3), v_eOb_b = zeros(3))
const XVel{T, D} = ComponentVector{T, D, typeof(getaxes(XVelTemplate))} where {T, D}

function renormalize_block!(x, ε) #returns true if norm was corrected
    norm_x = norm(x)
    abs(norm_x - 1.0) > ε ? (x ./= norm_x; return true) : return false
end

@inline function get_ω_el_n(v_eOb_n::AbstractVector{<:Real}, Ob::GeographicLocation)

    (R_N, R_E) = radii(Ob)
    h_e = AltE(Ob)

    return SVector{3}(
        v_eOb_n[2] / (R_E + Float64(h_e)),
        -v_eOb_n[1] / (R_N + Float64(h_e)),
        0.0)

end

@inline function get_ω_en_n(v_eOb_n::AbstractVector{<:Real}, Ob::GeographicLocation)

    (R_N, R_E) = radii(Ob)
    h_e = AltE(Ob)
    ϕ = LatLon(Ob).ϕ

    return SVector{3}(
        v_eOb_n[2] / (R_E + Float64(h_e)),
        -v_eOb_n[1] / (R_N + Float64(h_e)),
        -v_eOb_n[2] * tan(ϕ) / (R_E + Float64(h_e))
        )
end


########################### LTF-based Kinematics #########################
##########################################################################

#Fast, robust, singularity-free (all-attitude, all-latitude) kinematic
#description, appropriate for simulation.

#The salient feature of this kinematic description is that the azimuth of the
#local tangent frame (LTF) is not slaved to the geographic North, as it is in a
#NED-based kinematic description. Instead, the vertical component of the LTF's
#transport rate is arbitrarily set to zero. This avoids polar singularities, but
#also means that the azimuth angle of the LTF with respect to the geographic
#North will fluctuate as the vehicle moves around the Earth's surface. The
#resulting LTF is sometimes known as Wander Azimuth Frame. Position is defined
#by the rotation from the ECEF axes to the LTF axes, and attitude is defined by
#the rotation from the LTF axes to the vehicle axes.

struct LTF <: AbstractKinematics end

const XPosLTFTemplate = ComponentVector(q_lb = zeros(4), q_el = zeros(4), Δx = 0.0, Δy = 0.0, h_e = 0.0)
const XLTFTemplate = ComponentVector(pos = similar(XPosLTFTemplate), vel = similar(XVelTemplate))
const XLTF{T, D} = ComponentVector{T, D, typeof(getaxes(XLTFTemplate))} where {T, D}
x_template(::LTF) = XLTFTemplate

function init!(x::XLTF, ic::Initializer = Initializer())

    @unpack q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy = ic

    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_eb_b = ω_el_b + ω_lb_b
    v_eOb_b = q_nb'(v_eOb_n)
    h_e = AltE(Ob)

    q_lb = q_nb #arbitrarily initialize ψ_nl to 1

    x.pos.q_lb .= q_lb[:]
    x.pos.q_el .= ltf(Ob)[:]
    x.pos.Δx = Δx
    x.pos.Δy = Δy
    x.pos.h_e = h_e
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

end

function Y(x::XLTF)

    x_pos = x.pos; x_vel = x.vel
    q_lb = RQuat(x_pos.q_lb, normalization = false)
    q_el = RQuat(x_pos.q_el, normalization = false)
    ω_eb_b = SVector{3}(x_vel.ω_eb_b)
    v_eOb_b = SVector{3}(x_vel.v_eOb_b)
    h_e = AltE(x_pos.h_e[1])
    Δxy = SVector(x_pos.Δx, x_pos.Δy)

    ψ_nl = get_ψ_nl(q_el)
    q_nl = Rz(ψ_nl)
    q_nb = q_nl ∘ q_lb
    q_eb = q_el ∘ q_lb

    n_e = NVector(q_el)
    h_o = AltO(h_e, n_e)

    v_eOb_n = q_nb(v_eOb_b)
    Ob = GeographicLocation(n_e, h_e)
    ω_el_n = get_ω_el_n(v_eOb_n, Ob)

    ω_el_l = q_nl'(ω_el_n)
    ω_el_b = q_lb'(ω_el_l)
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector(0., 0., ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    return Y(
        Common(q_nb, q_eb, n_e, h_e, h_o, Δxy, ω_lb_b, ω_eb_b, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n),
        LTFSpecific(; q_lb, q_el, ω_el_l)
    )

end

#only updates xpos_dot, f_dyn! performs the xvel_dot update
function f_cont!(sys::System{LTF})

    #compute and update y
    sys.y = Y(sys.x)

    @unpack q_lb, q_el, ω_lb_b, ω_el_l, v_eOb_n = sys.y

    #update ẋ_pos
    ẋ_pos = sys.ẋ.pos
    ẋ_pos.q_lb .= Attitude.dt(q_lb, ω_lb_b)
    ẋ_pos.q_el .= Attitude.dt(q_el, ω_el_l)
    ẋ_pos.Δx = v_eOb_n[1]
    ẋ_pos.Δy = v_eOb_n[2]
    ẋ_pos.h_e = -v_eOb_n[3]

end

function f_disc!(sys::System{LTF}, ε = 1e-10)
    #we need both calls executed, so | must be used here instead of ||
    x_pos = sys.x.pos
    renormalize_block!(x_pos.q_lb, ε) | renormalize_block!(x_pos.q_el, ε)
end


########################## ECEF-based Kinematics #########################
##########################################################################

#Robust, singularity-free (all-attitude, all-latitude) kinematic description,
#appropriate for simulation.

struct ECEF <: AbstractKinematics end

const XPosECEFTemplate = ComponentVector(q_eb = zeros(4), n_e = zeros(3), Δx = 0.0, Δy = 0.0, h_e = 0.0)
const XECEFTemplate = ComponentVector(pos = similar(XPosECEFTemplate), vel = similar(XVelTemplate))
const XECEF{T, D} = ComponentVector{T, D, typeof(getaxes(XECEFTemplate))} where {T, D}

x_template(::ECEF) = XECEFTemplate

function init!(x::XECEF, ic::Initializer = Initializer())

    @unpack q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy = ic

    n_e = NVector(Ob)
    h_e = AltE(Ob)

    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_eb_b = ω_el_b + ω_lb_b
    v_eOb_b = q_nb'(v_eOb_n)

    x.pos.q_eb .= q_eb[:]
    x.pos.n_e .= n_e[:]
    x.pos.Δx = Δx
    x.pos.Δy = Δy
    x.pos.h_e = h_e
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

end

function Y(x::XECEF)

    x_pos = x.pos; x_vel = x.vel
    q_eb = RQuat(x_pos.q_eb, normalization = false)
    n_e = NVector(x_pos.n_e, normalization = false)
    ω_eb_b = SVector{3}(x_vel.ω_eb_b)
    v_eOb_b = SVector{3}(x_vel.v_eOb_b)
    Δxy = SVector(x_pos.Δx, x_pos.Δy)
    h_e = AltE(x_pos.h_e[1])
    h_o = AltO(h_e, n_e)

    q_en = ltf(n_e)
    q_nb = q_en' ∘ q_eb

    Ob = GeographicLocation(n_e, h_e)
    v_eOb_n = q_nb(v_eOb_b)
    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector(0., 0., ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    return Y(
        Common(q_nb, q_eb, n_e, h_e, h_o, Δxy, ω_lb_b, ω_eb_b, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n),
        ECEFSpecific(; q_en, ω_el_n)
    )

end

#only updates xpos_dot, xvel_dot update can only be performed by f_dyn!
function f_cont!(sys::System{ECEF})

    #compute and update y
    sys.y = Y(sys.x)

    @unpack q_eb, q_en, ω_el_n, ω_eb_b, v_eOb_n = sys.y

    #update ẋ_pos
    ẋ_pos = sys.ẋ.pos
    ẋ_pos.q_eb .= Attitude.dt(q_eb, ω_eb_b)
    ẋ_pos.n_e .= q_en(ω_el_n × SVector{3,Float64}(0,0,-1))
    ẋ_pos.Δx = v_eOb_n[1]
    ẋ_pos.Δy = v_eOb_n[2]
    ẋ_pos.h_e = -v_eOb_n[3]

end

function f_disc!(sys::System{ECEF}, ε = 1e-10)
    x_pos = sys.x.pos
    #we need both calls executed, so | must be used here instead of ||
    renormalize_block!(x_pos.q_eb, ε) | renormalize_block!(x_pos.n_e, ε)
end


################################ NED Kinematics ################################
################################################################################

#slower, non-singularity free implementation. useful for analysis and control
#design

struct NED <: AbstractKinematics end

const XPosNEDTemplate = ComponentVector(e_nb = zeros(3), ϕ = 0.0, λ = 0.0, Δx = 0.0, Δy = 0.0, h_e = 0.0)
const XNEDTemplate = ComponentVector(pos = similar(XPosNEDTemplate), vel = similar(XVelTemplate))
const XNED{T, D} = ComponentVector{T, D, typeof(getaxes(XNEDTemplate))} where {T, D}
x_template(::NED) = XNEDTemplate


function init!(x::XNED, ic::Initializer = Initializer())

    @unpack q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy = ic

    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_eb_b = ω_el_b + ω_lb_b
    v_eOb_b = q_nb'(v_eOb_n)

    e_nb = REuler(q_nb)
    ϕ_λ = LatLon(Ob)
    h_e = AltE(Ob)

    x.pos.e_nb .= SVector(e_nb.ψ, e_nb.θ, e_nb.φ)
    x.pos.ϕ = ϕ_λ.ϕ
    x.pos.λ = ϕ_λ.λ
    x.pos.Δx = Δx
    x.pos.Δy = Δy
    x.pos.h_e = h_e
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

end

function Y(x::XNED)

    e_nb = REuler(x.pos.e_nb)
    ϕ_λ = LatLon(x.pos.ϕ, x.pos.λ)
    h_e = AltE(x.pos.h_e[1])
    Δxy = SVector(x.pos.Δx, x.pos.Δy)
    ω_eb_b = SVector{3}(x.vel.ω_eb_b)
    v_eOb_b = SVector{3}(x.vel.v_eOb_b)

    n_e = NVector(ϕ_λ)
    h_o = AltO(h_e, n_e)

    q_nb = RQuat(e_nb)
    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    v_eOb_n = q_nb(v_eOb_b)
    Ob = GeographicLocation(n_e, h_e)

    ω_en_n = get_ω_en_n(v_eOb_n, Ob)
    ω_en_b = q_nb'(ω_en_n)
    ω_nb_b = ω_eb_b - ω_en_b

    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector{3}(0, 0, ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    return Y(
        Common(q_nb, q_eb, n_e, h_e, h_o, Δxy, ω_lb_b, ω_eb_b, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n),
        NEDSpecific(; e_nb, ϕ_λ, ω_nb_b, ω_en_n)
    )

end

#only updates xpos_dot, xvel_dot update can only be performed by f_dyn!
function f_cont!(sys::System{NED})

    #compute and update y
    sys.y = Y(sys.x)

    @unpack e_nb, ϕ_λ, ω_nb_b, ω_en_n, v_eOb_n = sys.y

    #update ẋ_pos
    ẋ_pos = sys.ẋ.pos
    ẋ_pos.e_nb .= Attitude.dt(e_nb, ω_nb_b)
    ẋ_pos.ϕ = -ω_en_n[2] #can be verified in [Groves]
    ẋ_pos.λ = ω_en_n[1] / cos(ϕ_λ.ϕ) #can be verified in [Groves]
    ẋ_pos.Δx = v_eOb_n[1]
    ẋ_pos.Δy = v_eOb_n[2]
    ẋ_pos.h_e = -v_eOb_n[3]

end

f_disc!(sys::System{NED}) = false


################################# Plotting #####################################

# #@userplot allows defining a custom plot for a specific dataset without having
# #to create a custom type for dispatch. we just wrap the data in the userplot
# #type generated by Plots, and is received inside the recipe in its field "args"
@userplot Trajectory3D
@recipe function f(t3d::Trajectory3D)

    # https://daschw.github.io/recipes/#series_recipes

    xs, ys, zs = t3d.args
    @assert length(xs) == length(ys) == length(zs)
    n = length(xs)

    xe, ye, ze = map(extrema, (xs, ys, zs))
    x_mid, y_mid, _ = map(v -> 0.5sum(v), (xe, ye, ze))
    x_span, y_span, z_span = map(v -> v[2] - v[1], (xe, ye, ze))
    span = max(x_span, y_span, z_span)

    xl = (x_mid - 0.5span, x_mid + 0.5span)
    yl = (y_mid - 0.5span, y_mid + 0.5span)
    zl = (ze[1], ze[1] + span)

    seriestype --> :path
    xguide --> L"$\Delta x\ (m)$"
    yguide --> L"$\Delta y\ (m)$"
    zguide --> L"$h\ (m)$"
    legend --> false

    xlims --> xl
    ylims --> yl
    zlims --> zl

    yflip --> true

    @series begin
        linecolor --> :lightgray
        xs, ys, fill(zl[1], n)
    end

    @series begin
        linecolor --> :lightgray
        xs, fill(yl[1], n), zs
    end

    @series begin
        linecolor --> :lightgray
        fill(xl[1], n), ys, zs
    end

    @series begin
        linecolor --> :blue
        linewidth --> 3
        xs, ys, zs
    end

    return nothing

end

function make_plots(th::TimeHistory{<:Y}; kwargs...)

    return OrderedDict(
        :common => make_plots(th.common; kwargs...),
        # :vel => make_plots(th.vel; kwargs...)
    )

end

function make_plots(th::TimeHistory{<:Common}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    plot_level = get(kwargs, :plot_level, :full)

    #example of capturing the plot_level keyword to control which plots are generated
    if plot_level == :simplified
        return pd #nothing also works
    end

    pd[:e_nb] = plot(
        th.q_nb; #will automatically be converted to REuler for plotting
        plot_title = "Attitude (Vehicle/NED)",
        rot_ref = "n", rot_target = "b",
        kwargs...)

    #will be automatically converted to LatLon for plotting
    #remove the title added by the LatLon TH recipe
    subplot_latlon = plot(th.n_e; title = "", th_split = :v, kwargs...)

    #remove the title added by the Altitude TH recipe
    subplot_h = plot(th.h_e; title = "", kwargs...)
                plot!(th.h_o; title = "", kwargs...)

    subplot_xy = plot(
        th.Δxy;
        label = [L"$\int v_{eb}^{x_n} dt$" L"$\int v_{eb}^{y_n} dt$"],
        ylabel = [L"$\Delta x\ (m)$" L"$\Delta y \ (m)$"],
        th_split = :h, link = :none, kwargs...)

    pd[:Ob_geo] = plot(
        subplot_latlon, subplot_h;
        layout = grid(1, 2, widths = [0.67, 0.33]),
        plot_title = "Position (WGS84)",
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd[:Ob_xyh] = plot(
        subplot_xy, subplot_h;
        layout = grid(1, 2, widths = [0.67, 0.33]),
        plot_title = "Position (Local Cartesian)",
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    #when a plot is assembled from multiple subplots, the plot_titlefontsize
    #attribute no longer works, and it is titlefontisze what determines the font
    #size of the overall figure title (which normally is used for subplots).

    th_Δx, th_Δy = Utils.get_scalar_components(th.Δxy)
    xs, ys, zs = th_Δx._data, th_Δy._data, Float64.(th.h_e._data)

    pd[:Ob_t3d] = plot(
        Trajectory3D((xs, ys, zs));
        plot_title = "Trajectory (Local Cartesian, Ellipsoidal Altitude)",
        titlefontsize = 20,
        camera = (30, 45),
        kwargs...
        )

    pd[:ω_lb_b] = plot(
        th.ω_lb_b;
        plot_title = "Angular Velocity (Vehicle/LTF) [Vehicle Axes]",
        label = ["Roll Rate" "Pitch Rate" "Yaw Rate"],
        ylabel = [L"$p \ (rad/s)$" L"$q \ (rad/s)$" L"$r \ (rad/s)$"],
        th_split = :h,
        kwargs...)

    pd[:v_eOb_n] = plot(
        th.v_eOb_n;
        plot_title = "Velocity (Vehicle/ECEF) [NED Axes]",
        label = ["North" "East" "Down"],
        ylabel = [L"$v_{eb}^{N} \ (m/s)$" L"$v_{eb}^{E} \ (m/s)$" L"$v_{eb}^{D} \ (m/s)$"],
        th_split = :h,
        kwargs...)

    pd[:v_eOb_b] = plot(
        th.v_eOb_b;
        plot_title = "Velocity (Vehicle/ECEF) [Vehicle Axes]",
        ylabel = [L"$v_{eb}^{x_b} \ (m/s)$" L"$v_{eb}^{y_b} \ (m/s)$" L"$v_{eb}^{z_b} \ (m/s)$"],
        th_split = :h,
        kwargs...)

    # pd[:ω_el_n] = plot(
    #     th.ω_el_n;
    #     plot_title = "Local Tangent Frame Transport Rate (LTF/ECEF) [NED Axes]",
    #     ylabel = L"$\omega_{el}^{l} \ (rad/s)$",
    #     th_split = :h,
    #     kwargs...)

    return pd

end


end #module