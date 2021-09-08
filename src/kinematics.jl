module Kinematics

using LinearAlgebra
using StaticArrays: SVector
using ComponentArrays
using UnPack

using Flight.Geodesy
using Flight.Attitude
using Flight.System
import Flight.System: get_x0

using Flight.Plotting
import Flight.Plotting: plots

export Pos, Vel, Kin, KinInit
export PosX, PosY, VelX, VelY, KinX, KinY
export init!, f_kin!, renormalize!

struct Pos <: AbstractComponent end
struct Vel <: AbstractComponent end
struct Kin <: AbstractComponent end

Base.@kwdef struct KinInit
    ω_lb_b::SVector{3, Float64} = zeros(SVector{3})
    v_eOb_b::SVector{3, Float64} = zeros(SVector{3})
    q_nb::RQuat = RQuat()
    Ob::Geographic{NVector,Ellipsoidal} = Geographic()
    Δx::Float64 = 0.0
    Δy::Float64 = 0.0
end

const PosXTemplate = ComponentVector(q_lb = zeros(4), q_el = zeros(4), Δx = 0.0, Δy = 0.0, h_e = 0.0)
const VelXTemplate = ComponentVector(ω_eb_b = zeros(3), v_eOb_b = zeros(3))
const KinXTemplate = ComponentVector(pos = PosXTemplate, vel = VelXTemplate)

"""
Type definition for dispatching on position state vector instances"

"""
const PosX{T, D} = ComponentVector{T, D, typeof(getaxes(PosXTemplate))} where {T, D}
"Type definition for dispatching on velocity state vector instances"
const VelX{T, D} = ComponentVector{T, D, typeof(getaxes(VelXTemplate))} where {T, D}
"Type definition for dispatching on velocity state vector instances"
const KinX{T, D} = ComponentVector{T, D, typeof(getaxes(KinXTemplate))} where {T, D}

Base.@kwdef struct PosY
    q_lb::RQuat = RQuat()
    q_nb::RQuat = RQuat()
    q_eb::RQuat = RQuat()
    e_nb::REuler = REuler()
    ψ_nl::Float64 = 0.0
    q_el::RQuat = RQuat() #attitude of the local tangent frame at Ob
    n_e::NVector = NVector()
    ϕ_λ::LatLon = LatLon()
    h_e::Altitude{Ellipsoidal} = Altitude{Ellipsoidal}()
    h_o::Altitude{Orthometric} = Altitude{Orthometric}()
    Δxy::SVector{2,Float64} = zeros(2) #v_eOb_n[1:2] integrals
end

Base.@kwdef struct VelY
    ω_eb_b::SVector{3,Float64} = zeros(3)
    ω_lb_b::SVector{3,Float64} = zeros(3)
    ω_el_l::SVector{3,Float64} = zeros(3)
    ω_ie_b::SVector{3,Float64} = zeros(3)
    ω_ib_b::SVector{3,Float64} = zeros(3)
    v_eOb_b::SVector{3,Float64} = zeros(3)
    v_eOb_n::SVector{3,Float64} = zeros(3)
end

Base.@kwdef struct KinY
    pos::PosY = PosY()
    vel::VelY = VelY()
end

#Kin is not a System, so we do not really need to define get_x0 to comply with the
#System interface. however, it is convenient for testing, and to ensure the
#aircraft state has its kinematic block initialized to reasonable values
get_x0(::Kin) = get_x0(KinInit())
get_x0(init::KinInit) = (x=similar(KinXTemplate); init!(x, init); return x)

function init!(x::KinX, init::KinInit)

    @unpack q_nb, Ob, ω_lb_b, v_eOb_b, Δx, Δy = init

    h_e = Ob.alt
    (R_N, R_E) = radii(Ob)
    v_eOb_n = q_nb * v_eOb_b
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + h_e),
        -v_eOb_n[1] / (R_N + h_e),
        0.0)

    ω_el_b = q_nb' * ω_el_n
    ω_eb_b = ω_el_b + ω_lb_b

    q_lb = q_nb #arbitrarily initialize ψ_nl to -1

    x.pos.q_lb .= q_lb[:]
    x.pos.q_el .= ltf(Ob)[:]
    x.pos.Δx = Δx
    x.pos.Δy = Δy
    x.pos.h_e = h_e
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

end

function f_kin!(ẋ_pos::PosX, x::KinX)

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

    (R_N, R_E) = radii(n_e)
    v_eOb_n = q_nb(v_eOb_b)
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + h_e),
        -v_eOb_n[1] / (R_N + h_e),
        0.0)

    ω_el_l = q_nl'(ω_el_n)
    ω_el_b = q_lb'(ω_el_l)
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector{3}(0, 0, ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    #update ẋ_pos
    ẋ_pos.q_lb .= dt(q_lb, ω_lb_b)
    ẋ_pos.q_el .= dt(q_el, ω_el_l)
    ẋ_pos.Δx = v_eOb_n[1]
    ẋ_pos.Δy = v_eOb_n[2]
    ẋ_pos.h_e = -v_eOb_n[3]

    #build output
    y_pos = PosY(q_lb, q_nb, q_eb, REuler(q_nb), ψ_nl, q_el, n_e, LatLon(n_e),
                h_e, Altitude{Orthometric}(h_e, n_e), SVector{2}(x.pos.Δx, x.pos.Δy))

    y_vel = VelY(ω_eb_b, ω_lb_b, ω_el_l, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n)

    return KinY(y_pos, y_vel)

end

function renormalize!(x_kin::KinX, ε = 1e-10)
    #we need both calls executed, so | must be used here instead of ||
    renormalize_q!(x_kin.pos.q_lb, ε) | renormalize_q!(x_kin.pos.q_el, ε)
end

function renormalize_q!(x_q, ε) #returns true if norm restored, false otherwise
    norm_q = norm(x_q)
    abs(norm_q - 1.0) > ε ? (x_q ./= norm_q; return true) : return false
end


############################ Plotting ################################

function plots(t, data::AbstractVector{<:KinY}; mode, save_path, kwargs...)

    sa = StructArray(data)

    plots(t, sa.pos; mode, save_path, kwargs...)
    plots(t, sa.vel; mode, save_path, kwargs...)
end

function plots(t, data::AbstractVector{<:PosY}; mode, save_path, kwargs...)

    @unpack e_nb, ϕ_λ, h_e, h_o, Δxy = StructArray(data)

    plt_e_nb = thplot(t, e_nb;
        plot_title = "Attitude (Airframe/NED)",
        kwargs...)

    #remove the title added by the Altitude TH recipe
    splt_h = thplot(t, h_e; title = "", kwargs...)
    thplot!(t, h_o; title = "", kwargs...)

    #remove the title added by the LatLon TH recipe
    splt_latlon = thplot(t, ϕ_λ;
                         title = "",
                         th_split = :v,
                         kwargs...)

    splt_xy = thplot(t, Δxy;
        label = [L"$\int v_{eO_b}^{x_n} dt$" L"$\int v_{eO_b}^{y_n} dt$"],
        ylabel = [L"$\Delta x\ (m)$" L"$\Delta y \ (m)$"],
        th_split = :h,
        link = :none,
        kwargs...)

    #when we assemble a plot from multiple subplots, the plot_titlefontsize
    #attribute no longer works, and it is titlefontisze what determines the font
    #size of the overall figure title (which normally is used for subplots).
    #however, we can still override it specifically for this plot
    plt_geo = plot(splt_latlon, splt_h;
                layout = grid(1, 2, widths = [0.67, 0.33]),
                plot_title = "Position (WGS84)",
                kwargs...,
                titlefontsize = 20) #override titlefontsize after kwargs


    plt_xyh = plot(splt_xy, splt_h;
                layout = grid(1, 2, widths = [0.67, 0.33]),
                plot_title = "Position (Local Cartesian)",
                kwargs...,
                titlefontsize = 20)

    savefig(plt_e_nb, joinpath(save_path, "e_nb.png"))
    savefig(plt_geo, joinpath(save_path, "Ob_geo.png"))
    savefig(plt_xyh, joinpath(save_path, "Ob_xyh.png"))

    #debug mode plots:
    # wander angle

    #maybe add a Trajectory user recipe for Vectors of 3DVector so that i can
    #pass it a Vector series directly.
    #also trplot
    # Ob_xyh_voa = VectorOfArray(Ob_xyh)
    # plt_Ob_xyh_3D = plot(collect(view(Ob_xyh_voa,i,:) for i ∈ 1:3)...;
    #     camera = (30, 45))
    #aspect_ratio attribute does not work for 3d figures
    # savefig(plt_Ob_xyh_3D, joinpath(save_path, "Ob_xyh_3D.png"))
end

function plots(t, data::AbstractVector{<:VelY}; mode, save_path, kwargs...)

    @unpack v_eOb_b, v_eOb_n, ω_lb_b, ω_el_l = StructArray(data)

    plt_ω_lb_b = thplot(t, ω_lb_b;
        plot_title = "Angular Velocity (Airframe/LTF) [Airframe]",
        label = ["Roll Rate" "Pitch Rate" "Yaw Rate"],
        ylabel = [L"$p \ (rad/s)$" L"$q \ (rad/s)$" L"$r \ (rad/s)$"],
        th_split = :h,
        kwargs...)

    plt_ω_el_l = thplot(t, ω_el_l;
        plot_title = "LTF Transport Rate (LTF/ECEF) [LTF]",
        ylabel = L"$\omega_{el}^{l} \ (rad/s)$",
        th_split = :h,
        kwargs...)

    plt_v_eOb_n = thplot(t, v_eOb_n;
        plot_title = "Velocity (Airframe/ECEF) [NED]",
        label = ["North" "East" "Down"],
        ylabel = [L"$v_{eO_b}^{N} \ (m/s)$" L"$v_{eO_b}^{E} \ (m/s)$" L"$v_{eO_b}^{D} \ (m/s)$"],
        th_split = :h,
        kwargs...)

    plt_v_eOb_b = thplot(t, v_eOb_b;
        plot_title = "Velocity (Airframe/ECEF) [Airframe]",
        ylabel = [L"$v_{eO_b}^{x_b} \ (m/s)$" L"$v_{eO_b}^{y_b} \ (m/s)$" L"$v_{eO_b}^{z_b} \ (m/s)$"],
        th_split = :h,
        kwargs...)

    savefig(plt_ω_lb_b, joinpath(save_path, "ω_lb_b.png"))
    savefig(plt_ω_el_l, joinpath(save_path, "ω_el_l.png"))
    savefig(plt_v_eOb_n, joinpath(save_path, "v_eOb_n.png"))
    savefig(plt_v_eOb_b, joinpath(save_path, "v_eOb_b.png"))

end

end #module