module Kinematics

using LinearAlgebra
using StaticArrays: SVector
using ComponentArrays
using UnPack

using Flight.WGS84
using Flight.Attitude
using Flight.System
import Flight.System: x0

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
    Ob::NVectorAlt = NVectorAlt()
    Δx::Float64 = 0.0
    Δy::Float64 = 0.0
end

const PosXTemplate = ComponentVector(q_lb = zeros(4), q_el = zeros(4), Δx = 0.0, Δy = 0.0, h = 0.0)
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

Base.@kwdef struct PosY <: AbstractY{Pos}
    q_lb::RQuat
    q_nb::RQuat
    q_eb::RQuat
    e_nb::REuler
    ψ_nl::Float64
    q_el::RQuat #attitude of the local tangent frame at Ob
    Ob_nvh::NVectorAlt #position of Ob, NVectorAlt descriptor
    Ob_llh::LatLonAlt #position of Ob, LatLonAlt descriptor
    Ob_xyh::SVector{3,Float64} #velocity integral
end

Base.@kwdef struct VelY <: AbstractY{Vel}
    ω_eb_b::SVector{3,Float64}
    ω_lb_b::SVector{3,Float64}
    ω_el_l::SVector{3,Float64}
    ω_ie_b::SVector{3,Float64}
    ω_ib_b::SVector{3,Float64}
    v_eOb_b::SVector{3,Float64}
    v_eOb_n::SVector{3,Float64}
end

Base.@kwdef struct KinY <: AbstractY{Kin}
    pos::PosY
    vel::VelY
end

#Kin is not a System, so we do not really need to define x0 to comply with the
#System interface. however, it is convenient for testing, and to ensure the
#aircraft state has its kinematic block initialized to reasonable values
x0(::Kin) = x0(KinInit())
x0(init::KinInit) = (x=similar(KinXTemplate); init!(x, init); return x)

function init!(x::KinX, init::KinInit)

    @unpack q_nb, Ob, ω_lb_b, v_eOb_b, Δx, Δy = init

    h = Ob.h[1]
    (R_N, R_E) = radii(Ob)
    v_eOb_n = q_nb * v_eOb_b
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + h),
        -v_eOb_n[1] / (R_N + h),
        0.0)

    ω_el_b = q_nb' * ω_el_n
    ω_eb_b = ω_el_b + ω_lb_b

    q_lb = q_nb #arbitrarily initialize ψ_nl to -1

    x.pos.q_lb .= q_lb[:]
    x.pos.q_el .= ltf(Ob)[:]
    x.pos.Δx = Δx
    x.pos.Δy = Δy
    x.pos.h = h
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

end

function f_kin!(ẋ_pos::PosX, x::KinX)

    q_lb = RQuat(x.pos.q_lb, normalization = false)
    q_el = RQuat(x.pos.q_el, normalization = false)
    Δx = x.pos.Δx
    Δy = x.pos.Δy
    h = x.pos.h[1]
    ω_eb_b = SVector{3}(x.vel.ω_eb_b)
    v_eOb_b = SVector{3}(x.vel.v_eOb_b)

    Ob_nvh = NVectorAlt(NVector(q_el), h)
    _ψ_nl = ψ_nl(q_el)
    q_nl = Rz(_ψ_nl)
    q_nb = q_nl ∘ q_lb
    q_eb = q_el ∘ q_lb

    (R_N, R_E) = radii(Ob_nvh)
    v_eOb_n = q_nb(v_eOb_b)
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + h),
        -v_eOb_n[1] / (R_N + h),
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
    ẋ_pos.h = -v_eOb_n[3]

    #build outputs
    y_pos = PosY(q_lb, q_nb, q_eb, REuler(q_nb), _ψ_nl, q_el, Ob_nvh,
                 LatLonAlt(Ob_nvh), SVector{3}(Δx, Δy, h))
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

    @unpack e_nb, Ob_xyh = StructArray(data)
    # println(mode); println(save_path); println(kwargs)
    plt_e_nb = thplot(t, e_nb; label = ["Heading" "Inclination" "Bank"], kwargs...)
    plt_Ob_xyh = thplot(t, Ob_xyh;
        label = [L"$\int v_{eO_b}^{x_n} dt$" L"$\int v_{eO_b}^{y_n} dt$" L"$\mathrm{Altitude}$"],
        plot_title = "Local Cartesian Coordinates",
        ylabel = [L"$\Delta x\quad (m)$" L"$\Delta y \quad (m)$" L"h \quad (m)"],
        split = :h, link = :none, kwargs...)

    savefig(plt_e_nb, joinpath(save_path, "e_nb.png"))
    savefig(plt_Ob_xyh, joinpath(save_path, "Ob_xyh.png"))
end

function plots(t, data::AbstractVector{<:VelY}; mode, save_path, kwargs...)

    @unpack v_eOb_b, v_eOb_n, ω_lb_b = StructArray(data)
    plt_ω_lb_b = thplot(t, ω_lb_b; label = ["Yaw rate" "Pitch rate" "Roll rate"], kwargs...)
    savefig(plt_ω_lb_b, joinpath(save_path, "ω_lb_b.png"))
    # plot(t, pos_data; mode, save_path, kwargs...)
end

end #module