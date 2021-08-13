module Kinematics

using LinearAlgebra
using StaticArrays: SVector
using ComponentArrays
using UnPack

using Flight.WGS84
using Flight.Attitude
import Flight.System: X

export Pos, Vel, PosVel, PosVelInit
export PosX, VelX, PosVelX
export PosY, VelY, PosVelY, AccY, KinY
export init!, f_pos!

abstract type KinematicStruct end

struct Pos end
struct Vel end
struct PosVel end

const PosXTemplate = ComponentVector(q_lb = zeros(4), q_el = zeros(4), Δx = 0.0, Δy = 0.0, h = 0.0)
const VelXTemplate = ComponentVector(ω_eb_b = zeros(3), v_eOb_b = zeros(3))
const PosVelXTemplate = ComponentVector(pos = PosXTemplate, vel = VelXTemplate)

const PosX{D} = ComponentVector{Float64, D, typeof(getaxes(PosXTemplate))} where {D<:AbstractVector{Float64}}
const VelX{D} = ComponentVector{Float64, D, typeof(getaxes(VelXTemplate))} where {D<:AbstractVector{Float64}}
const PosVelX{D} = ComponentVector{Float64, D, typeof(getaxes(PosVelXTemplate))} where {D<:AbstractVector{Float64}}

Base.@kwdef struct PosVelInit <: KinematicStruct
    ω_lb_b::SVector{3, Float64} = zeros(SVector{3})
    v_eOb_b::SVector{3, Float64} = zeros(SVector{3})
    q_nb::RQuat = RQuat()
    Ob::WGS84Pos = WGS84Pos()
    Δx::Float64 = 0.0
    Δy::Float64 = 0.0
end

struct PosY <: KinematicStruct
    q_lb::RQuat
    q_nl::RQuat
    q_nb::RQuat
    q_eb::RQuat
    q_el::RQuat
    Ob::WGS84Pos
    Δx::Float64
    Δy::Float64
end

struct VelY <: KinematicStruct
    ω_eb_b::SVector{3,Float64}
    ω_lb_b::SVector{3,Float64}
    ω_el_l::SVector{3,Float64}
    ω_ie_b::SVector{3,Float64}
    ω_ib_b::SVector{3,Float64}
    v_eOb_b::SVector{3,Float64}
    v_eOb_n::SVector{3,Float64}
end

struct PosVelY <: KinematicStruct
    pos::PosY
    vel::VelY
end

struct AccY <: KinematicStruct
    α_eb_b::SVector{3,Float64}
    α_ib_b::SVector{3,Float64}
    a_eOb_b::SVector{3,Float64}
    a_iOb_b::SVector{3,Float64}
end

struct KinY <: KinematicStruct
    pos::PosY
    vel::VelY
    acc::AccY
end

X(::PosVel) = X(PosVelInit())
X(init::PosVelInit) = (x=similar(PosVelXTemplate); init!(x, init); return x)

function init!(x::PosVelX, init::PosVelInit)

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

    x.pos.q_lb .= q_lb #assignment without colon slicing courtesy of RQuat iterability
    x.pos.q_el .= ltf(Ob) #assignment without colon slicing courtesy of RQuat iterability
    x.pos.Δx = Δx
    x.pos.Δy = Δy
    x.pos.h = h
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

end

function f_pos!(ẋ_pos::PosX, x::PosVelX)::PosVelY

    #careful here: x.pos.h, x.vel.ω_eb_b and x.vel.v_eOb_b create views (this is
    #how LBV behaves by design). to copy the data, we can extract their
    #components using slices
    q_lb = RQuat(x.pos.q_lb)
    q_el = RQuat(x.pos.q_el)
    Δx = x.pos.Δx
    Δy = x.pos.Δy
    h = x.pos.h[1]
    ω_eb_b = SVector{3}(x.vel.ω_eb_b)
    v_eOb_b = SVector{3}(x.vel.v_eOb_b)

    Ob = WGS84Pos(NVector(q_el), h)
    q_nl = Rz(ψ_nl(q_el))
    q_nb = q_nl ∘ q_lb
    q_eb = q_el ∘ q_lb

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
    ẋ_pos.q_lb .= Attitude.dot(q_lb, ω_lb_b)
    ẋ_pos.q_el .= Attitude.dot(q_el, ω_el_l)
    ẋ_pos.Δx = v_eOb_n[1]
    ẋ_pos.Δy = v_eOb_n[2]
    ẋ_pos.h = -v_eOb_n[3]

    #build outputs
    pos = PosY(q_lb, q_nl, q_nb, q_eb, q_el, Ob, Δx, Δy)
    vel = VelY(ω_eb_b, ω_lb_b, ω_el_l, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n)

    return PosVelY(pos, vel)

end

function print_data_oneline(data::KinematicStruct, io::IO)
    fnames = fieldnames(typeof(data))
    print(io, typeof(data), "(")
    for f in fnames
        print(io, f, " = ", getproperty(data, f))
        f == fnames[end] ? print(io, ")") : print(io, ", ")
    end
end

function print_data_multiline(data::KinematicStruct, io::IO)
    fnames = fieldnames(typeof(data))
    print(io, typeof(data), " with fields:\n")
    for f in fnames
        print(io, "\t", f, ": ", getproperty(data, f))
        f == fnames[end] ? break : println(io)
    end
end

Base.show(io::IO, data::KinematicStruct) = print_data_oneline(data, io)
Base.show(io::IO, ::MIME"text/plain", data::KinematicStruct) = print_data_multiline(data, io)

end #module