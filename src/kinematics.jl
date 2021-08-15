module Kinematics

using LinearAlgebra
using StaticArrays: SVector
using ComponentArrays
using UnPack

using Flight.WGS84
using Flight.Attitude
import Flight.System: X, Y

export Pos, Vel, Kin, KinInit
export PosX, PosY, VelX, VelY, KinX, KinY
export init!, f_pos!

abstract type KinematicStruct end

struct Pos end
struct Vel end
struct Kin end

Base.@kwdef struct KinInit <: KinematicStruct
    ω_lb_b::SVector{3, Float64} = zeros(SVector{3})
    v_eOb_b::SVector{3, Float64} = zeros(SVector{3})
    q_nb::RQuat = RQuat()
    Ob::WGS84Pos = WGS84Pos()
    Δx::Float64 = 0.0
    Δy::Float64 = 0.0
end

const PosXTemplate = ComponentVector(q_lb = zeros(4), q_el = zeros(4), Δx = 0.0, Δy = 0.0, h = 0.0)
const PosYTemplate = ComponentVector(
    q_lb = zeros(4), q_nb = zeros(4), q_eb = zeros(4),
    ψ_nl = 0, ψ_nb = 0, θ_nb = 0, φ_nb = 0,
    q_el = zeros(4), ϕ = 0, λ = 0, h = 0, Δx = 0, Δy = 0)
const VelXTemplate = ComponentVector(ω_eb_b = zeros(3), v_eOb_b = zeros(3))
const VelYTemplate = ComponentVector(
    ω_eb_b = zeros(3), ω_lb_b = zeros(3), ω_el_l = zeros(3), ω_ie_b = zeros(3),
    ω_ib_b = zeros(3), v_eOb_b = zeros(3), v_eOb_n = zeros(3))
const KinXTemplate = ComponentVector(pos = PosXTemplate, vel = VelXTemplate)
const KinYTemplate = ComponentVector(pos = PosYTemplate, vel = VelYTemplate)

const PosX{D} = ComponentVector{Float64, D, typeof(getaxes(PosXTemplate))} where {D<:AbstractVector{Float64}}
const PosY{D} = ComponentVector{Float64, D, typeof(getaxes(PosYTemplate))} where {D<:AbstractVector{Float64}}
const VelX{D} = ComponentVector{Float64, D, typeof(getaxes(VelXTemplate))} where {D<:AbstractVector{Float64}}
const VelY{D} = ComponentVector{Float64, D, typeof(getaxes(VelYTemplate))} where {D<:AbstractVector{Float64}}

const KinX{D} = ComponentVector{Float64, D, typeof(getaxes(KinXTemplate))} where {D<:AbstractVector{Float64}}
const KinY{D} = ComponentVector{Float64, D, typeof(getaxes(KinYTemplate))} where {D<:AbstractVector{Float64}}

#Kin is not a System, so we do not really need to define Y, U and D to comply
#with the System interface. however, these are convenient for testing, and to
#ensure AircraftXTemplate has its kinematic block initialized to reasonable values
X(::Kin) = X(KinInit())
X(init::KinInit) = (x=similar(KinXTemplate); init!(x, init); return x)
Y(::Pos) = similar(PosYTemplate)
Y(::Vel) = similar(VelYTemplate)
Y(::Kin) = similar(KinYTemplate)

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

function f_pos!(y::KinY, ẋ_pos::PosX, x::KinX)

    #careful here: x.pos.h, x.vel.ω_eb_b and x.vel.v_eOb_b create views (this is
    #how LBV behaves by design). to copy the data, we can extract their
    #components using slices
    q_lb = RQuat(x.pos.q_lb, normalization = false)
    q_el = RQuat(x.pos.q_el, normalization = false)
    Δx = x.pos.Δx
    Δy = x.pos.Δy
    h = x.pos.h[1]
    ω_eb_b = SVector{3}(x.vel.ω_eb_b)
    v_eOb_b = SVector{3}(x.vel.v_eOb_b)

    Ob = WGS84Pos(NVector(q_el), h)
    ψ_nl_tmp = ψ_nl(q_el)
    q_nl = Rz(ψ_nl_tmp)
    q_nb = q_nl ∘ q_lb
    q_eb = q_el ∘ q_lb
    euler_nb = REuler(q_nb)

    (R_N, R_E) = radii(Ob)
    v_eOb_n = q_nb * v_eOb_b
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + h),
        -v_eOb_n[1] / (R_N + h),
        0.0)

    ω_el_l = q_nl' * ω_el_n
    ω_el_b = q_lb' * ω_el_l
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector{3}(0, 0, ω_ie)
    ω_ie_b = q_eb' * ω_ie_e
    ω_ib_b = ω_ie_b + ω_eb_b

    #update ẋ_pos
    ẋ_pos.q_lb .= dt(q_lb, ω_lb_b)
    ẋ_pos.q_el .= dt(q_el, ω_el_l)
    ẋ_pos.Δx = v_eOb_n[1]
    ẋ_pos.Δy = v_eOb_n[2]
    ẋ_pos.h = -v_eOb_n[3]

    #update y
    y_pos = y.pos
    y_pos.q_lb .= q_lb[:]
    y_pos.q_nb .= q_nb[:]
    y_pos.q_eb .= q_eb[:]
    y_pos.ψ_nl = ψ_nl_tmp
    y_pos.ψ_nb = euler_nb.ψ
    y_pos.θ_nb = euler_nb.θ
    y_pos.φ_nb = euler_nb.φ
    y_pos.q_el .= q_el[:]
    y_pos.ϕ = lat(Ob)
    y_pos.λ = lon(Ob)
    y_pos.h = h
    y_pos.Δx = Δx
    y_pos.Δy = Δy

    y_vel = y.vel
    y_vel.ω_eb_b .= ω_eb_b
    y_vel.ω_lb_b .= ω_lb_b
    y_vel.ω_el_l .= ω_el_l
    y_vel.ω_ie_b .= ω_ie_b
    y_vel.ω_ib_b .= ω_ib_b
    y_vel.v_eOb_b .= v_eOb_b
    y_vel.v_eOb_n .= v_eOb_n

    return nothing

end

end #module