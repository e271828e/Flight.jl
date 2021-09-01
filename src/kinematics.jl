module Kinematics

using LinearAlgebra
using StaticArrays: SVector
using ComponentArrays
using UnPack

using Flight.WGS84
using Flight.Attitude
import Flight.System: x0

export Pos, Vel, Kin, KinInit
export PosX, PosY, VelX, VelY, KinX, KinY
export init!, f_kin!, renormalize!

abstract type KinematicStruct end

struct Pos end
struct Vel end
struct Kin end

Base.@kwdef struct KinInit <: KinematicStruct
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

Base.@kwdef struct PosY
    q_lb::RQuat
    q_nb::RQuat
    q_eb::RQuat
    e_nb::REuler
    ψ_nl::Float64
    q_el::RQuat
<<<<<<< HEAD
    Ob::NVectorAlt #may need to add a LatLonAlt field
=======
    Ob::NVectorAlt
    ϕ::Float64
    λ::Float64
    h::Float64
>>>>>>> 338417e (Updated WGS84 to an abstract type with subtypes NVectorAlt, LatLonAlt, CartECEF)
    Δx::Float64
    Δy::Float64
end

Base.@kwdef struct VelY
    ω_eb_b::SVector{3,Float64}
    ω_lb_b::SVector{3,Float64}
    ω_el_l::SVector{3,Float64}
    ω_ie_b::SVector{3,Float64}
    ω_ib_b::SVector{3,Float64}
    v_eOb_b::SVector{3,Float64}
    v_eOb_n::SVector{3,Float64}
end

Base.@kwdef struct KinY
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

    Ob = NVectorAlt(NVector(q_el), h)
    _ψ_nl = ψ_nl(q_el)
    q_nl = Rz(_ψ_nl)
    q_nb = q_nl ∘ q_lb
    q_eb = q_el ∘ q_lb

    (R_N, R_E) = radii(Ob)
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
    y_pos = PosY(q_lb, q_nb, q_eb, REuler(q_nb), _ψ_nl, q_el, Ob, Δx, Δy)
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


end #module