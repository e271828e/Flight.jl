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

export AbstractKinematics, KinLTF, KinECEF
export VelX, PosData, VelData, KinData, KinInit

abstract type AbstractKinematics <: SystemDescriptor end

Base.@kwdef struct KinInit
    ω_lb_b::SVector{3, Float64} = zeros(SVector{3})
    v_eOb_n::SVector{3, Float64} = zeros(SVector{3})
    q_nb::RQuat = RQuat()
    Ob::GeographicLocation{NVector,Ellipsoidal} = GeographicLocation()
    Δx::Float64 = 0.0
    Δy::Float64 = 0.0
end

struct PosData
    q_nb::RQuat
    q_eb::RQuat
    e_nb::REuler
    q_en::RQuat
    n_e::NVector
    ϕ_λ::LatLon
    h_e::Altitude{Ellipsoidal}
    h_o::Altitude{Orthometric}
    Δxy::SVector{2,Float64}
end

struct VelData
    ω_eb_b::SVector{3,Float64}
    ω_el_n::SVector{3,Float64}
    ω_lb_b::SVector{3,Float64}
    ω_ie_b::SVector{3,Float64}
    ω_ib_b::SVector{3,Float64}
    v_eOb_b::SVector{3,Float64}
    v_eOb_n::SVector{3,Float64}
end

struct KinData
    pos::PosData
    vel::VelData
end

function KinData(init::KinInit)

    @unpack q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy = init

    h_e = Ob.alt
    q_el = ltf(Ob)

    #arbitrarily initialize ψ_nl to 0
    q_nl = RQuat()
    q_lb = q_nb
    q_en = q_el
    q_eb = q_el ∘ q_lb

    n_e = NVector(q_el)

    v_eOb_b = q_nb'(v_eOb_n)

    (R_N, R_E) = radii(Ob)
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + Float64(h_e)),
        -v_eOb_n[1] / (R_N + Float64(h_e)),
        0.0)
    ω_el_b = q_nb'(ω_el_n)
    ω_eb_b = ω_el_b + ω_lb_b

    ω_el_l = q_nl'(ω_el_n)
    ω_el_b = q_lb'(ω_el_l)
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector{3}(0, 0, ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    pos = PosData(q_nb, q_eb, REuler(q_nb), q_en, n_e, LatLon(n_e), h_e,
        Altitude{Orthometric}(h_e, n_e), SVector{2}(Δx, Δy))

    vel = VelData(ω_eb_b, ω_el_n, ω_lb_b, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n)

    return KinData(pos, vel)

end

KinData() = KinData(KinInit())

#every AbstractKinematics implementation must comply with the same outputs
init(::SystemY, ::AbstractKinematics) = KinData()

#for dispatching
const VelXTemplate = ComponentVector(ω_eb_b = zeros(3), v_eOb_b = zeros(3))
const VelX{T, D} = ComponentVector{T, D, typeof(getaxes(VelXTemplate))} where {T, D}

function renormalize_block!(x, ε) #returns true if norm was corrected
    norm_x = norm(x)
    abs(norm_x - 1.0) > ε ? (x ./= norm_x; return true) : return false
end

######################### LTF Kinematics #########################

struct KinLTF <: AbstractKinematics end

init(::SystemX, kin::KinLTF) = init(kin, KinInit())

function init(::KinLTF, kin_init::KinInit = KinInit())

    @unpack q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy = kin_init

    x = ComponentVector(
        pos = ComponentVector(q_lb = zeros(4), q_el = zeros(4), Δx = 0.0, Δy = 0.0, h_e = 0.0),
        vel = similar(VelXTemplate) #using similar instead of simple assignment is essential
        )

    h_e = Ob.alt
    (R_N, R_E) = radii(Ob)
    v_eOb_b = q_nb'(v_eOb_n)
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + Float64(h_e)),
        -v_eOb_n[1] / (R_N + Float64(h_e)),
        0.0)

    ω_el_b = q_nb'(ω_el_n)
    ω_eb_b = ω_el_b + ω_lb_b

    q_lb = q_nb #arbitrarily initialize ψ_nl to 1

    x.pos.q_lb .= q_lb[:]
    x.pos.q_el .= ltf(Ob)[:]
    x.pos.Δx = Δx
    x.pos.Δy = Δy
    x.pos.h_e = h_e
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

    return x

end

#this only updates xpos_dot, we need f_dyn! to perform the xvel_dot update
function f_cont!(sys::System{KinLTF})

    x = sys.x

    q_lb = RQuat(x.pos.q_lb, normalization = false)
    q_el = RQuat(x.pos.q_el, normalization = false)
    ω_eb_b = SVector{3}(x.vel.ω_eb_b)
    v_eOb_b = SVector{3}(x.vel.v_eOb_b)
    h_e = Altitude{Ellipsoidal}(x.pos.h_e[1])

    n_e = NVector(q_el)
    ψ_nl = get_ψ_nl(q_el)
    q_nl = Rz(ψ_nl)
    q_nb = q_nl ∘ q_lb
    q_eb = q_el ∘ q_lb
    q_en = q_el ∘ q_nl'

    (R_N, R_E) = radii(n_e)
    v_eOb_n = q_nb(v_eOb_b)
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + Float64(h_e)),
        -v_eOb_n[1] / (R_N + Float64(h_e)),
        0.0)

    ω_el_l = q_nl'(ω_el_n)
    ω_el_b = q_lb'(ω_el_l)
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector{3}(0, 0, ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    #update ẋ_pos
    ẋ_pos = sys.ẋ.pos
    ẋ_pos.q_lb .= Attitude.dt(q_lb, ω_lb_b)
    ẋ_pos.q_el .= Attitude.dt(q_el, ω_el_l)
    ẋ_pos.Δx = v_eOb_n[1]
    ẋ_pos.Δy = v_eOb_n[2]
    ẋ_pos.h_e = -v_eOb_n[3]

    #build and assign output
    pos = PosData(q_nb, q_eb, REuler(q_nb), q_en, n_e, LatLon(n_e), h_e,
        Altitude{Orthometric}(h_e, n_e), SVector{2}(x.pos.Δx, x.pos.Δy))

    vel = VelData(ω_eb_b, ω_el_n, ω_lb_b, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n)

    sys.y = KinData(pos, vel)

end

function f_disc!(sys::System{KinLTF}, ε = 1e-10)
    #we need both calls executed, so | must be used here instead of ||
    x_pos = sys.x.pos
    renormalize_block!(x_pos.q_lb, ε) | renormalize_block!(x_pos.q_el, ε)
end

######################### ECEF Kinematics #########################

struct KinECEF <: AbstractKinematics end

init(::SystemX, kin::KinECEF) = init(kin, KinInit())

function init(::KinECEF, kin_init::KinInit = KinInit())

    @unpack q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy = kin_init

    x = ComponentVector(
        pos = ComponentVector(q_eb = zeros(4), n_e = zeros(3), Δx = 0.0, Δy = 0.0, h_e = 0.0),
        vel = similar(VelXTemplate) #using similar rather than assignment is essential here!
        )

    n_e = Ob.l2d
    h_e = Ob.alt
    (R_N, R_E) = radii(Ob)
    v_eOb_b = q_nb'(v_eOb_n)
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + Float64(h_e)),
        -v_eOb_n[1] / (R_N + Float64(h_e)),
        0.0)

    ω_el_b = q_nb'(ω_el_n)
    ω_eb_b = ω_el_b + ω_lb_b

    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    x.pos.q_eb .= q_eb[:]
    x.pos.n_e .= n_e[:]
    x.pos.Δx = Δx
    x.pos.Δy = Δy
    x.pos.h_e = h_e
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

    return x

end

#this only updates xpos_dot, we need f_dyn! to perform the xvel_dot update
function f_cont!(sys::System{KinECEF})

    x = sys.x

    q_eb = RQuat(x.pos.q_eb, normalization = false)
    n_e = NVector(x.pos.n_e, normalization = false)
    ω_eb_b = SVector{3}(x.vel.ω_eb_b)
    v_eOb_b = SVector{3}(x.vel.v_eOb_b)
    h_e = Altitude{Ellipsoidal}(x.pos.h_e[1])

    q_en = ltf(n_e)
    q_nb = q_en' ∘ q_eb

    (R_N, R_E) = radii(n_e)
    v_eOb_n = q_nb(v_eOb_b)
    ω_el_n = SVector{3,Float64}(
        v_eOb_n[2] / (R_E + Float64(h_e)),
        -v_eOb_n[1] / (R_N + Float64(h_e)),
        0.0)

    ω_el_b = q_nb'(ω_el_n)
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector{3,Float64}(0, 0, ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    #update ẋ_pos
    ẋ_pos = sys.ẋ.pos
    ẋ_pos.q_eb .= Attitude.dt(q_eb, ω_eb_b)
    ẋ_pos.n_e .= q_en(ω_el_n × SVector{3,Float64}(0,0,-1))
    ẋ_pos.Δx = v_eOb_n[1]
    ẋ_pos.Δy = v_eOb_n[2]
    ẋ_pos.h_e = -v_eOb_n[3]

    #build output
    pos = PosData(q_nb, q_eb, REuler(q_nb), q_en, n_e, LatLon(n_e), h_e,
        Altitude{Orthometric}(h_e, n_e), SVector{2}(x.pos.Δx, x.pos.Δy))

    vel = VelData(ω_eb_b, ω_el_n, ω_lb_b, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n)

    sys.y = KinData(pos, vel)

end

function f_disc!(sys::System{KinECEF}, ε = 1e-10)
    x_pos = sys.x.pos
    #we need both calls executed, so | must be used here instead of ||
    renormalize_block!(x_pos.q_eb, ε) | renormalize_block!(x_pos.n_e, ε)
end

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

function make_plots(th::TimeHistory{<:KinData}; kwargs...)

    return OrderedDict(
        :pos => make_plots(th.pos; kwargs...),
        :vel => make_plots(th.vel; kwargs...)
    )

end

function make_plots(th::TimeHistory{<:PosData}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    plot_level = get(kwargs, :plot_level, :full)

    #example of capturing the plot_level keyword to control which plots are generated
    if plot_level == :simplified
        return pd #nothing also works
    end

    pd[:e_nb] = plot(
        th.e_nb;
        plot_title = "Attitude (Vehicle/NED)",
        rot_ref = "n", rot_target = "b",
        kwargs...)

    #remove the title added by the LatLon TH recipe
    subplot_latlon = plot(th.ϕ_λ; title = "", th_split = :v, kwargs...)

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

    return pd

end

function make_plots(th::TimeHistory{<:VelData}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    pd[:ω_lb_b] = plot(
        th.ω_lb_b;
        plot_title = "Angular Velocity (Vehicle/LTF) [Vehicle Axes]",
        label = ["Roll Rate" "Pitch Rate" "Yaw Rate"],
        ylabel = [L"$p \ (rad/s)$" L"$q \ (rad/s)$" L"$r \ (rad/s)$"],
        th_split = :h,
        kwargs...)

    pd[:ω_el_n] = plot(
        th.ω_el_n;
        plot_title = "Local Tangent Frame Transport Rate (LTF/ECEF) [NED Axes]",
        ylabel = L"$\omega_{el}^{l} \ (rad/s)$",
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

    return pd

end


end #module