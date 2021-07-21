module Kinematics

using LinearAlgebra
using StaticArrays: SVector
using UnPack

using Flight.LBV
using Flight.WGS84
using Flight.Attitude

export KinInit, XVel
export XPosWGS84, XKinWGS84, VelDataWGS84, PosDataWGS84, PVDataWGS84, AccDataWGS84
export x_pos_dot

#for some reason, defining and using using the type alias SV3 =
#SVector{3,Float64} yields type instability

@define_node XPosWGS84 (q_lb = LBVLeaf{4}, q_el = LBVLeaf{4}, h = LBVLeaf{1})
@define_node XVel (ω_eb_b = LBVLeaf{3}, v_eOb_b = LBVLeaf{3})
@define_node XKinWGS84 (pos = XPosWGS84, vel = XVel)

abstract type KinematicData end

function print_data_oneline(data::KinematicData, io::IO)
    fnames = fieldnames(typeof(data))
    print(io, typeof(data), "(")
    for f in fnames
        print(io, f, " = ", getproperty(data, f))
        f == fnames[end] ? print(io, ")") : print(io, ", ")
    end
end

function print_data_multiline(data::KinematicData, io::IO)
    fnames = fieldnames(typeof(data))
    print(io, typeof(data), " with fields:\n")
    for f in fnames
        print(io, "\t", f, ": ", getproperty(data, f))
        f == fnames[end] ? break : println(io)
    end
end

Base.show(io::IO, data::KinematicData) = print_data_oneline(data, io)
Base.show(io::IO, ::MIME"text/plain", data::KinematicData) = print_data_multiline(data, io)

Base.@kwdef struct KinInit <: KinematicData
    q_nb::RQuat = RQuat()
    Ob::WGS84Pos = WGS84Pos()
    ω_lb_b::SVector{3, Float64} = zeros(SVector{3})
    v_eOb_b::SVector{3, Float64} = zeros(SVector{3})
end

Base.@kwdef struct PosDataWGS84 <: KinematicData
    q_lb::RQuat
    q_nl::RQuat
    q_nb::RQuat
    q_eb::RQuat
    q_el::RQuat
    Ob::WGS84Pos
end

Base.@kwdef struct VelDataWGS84 <: KinematicData
    ω_eb_b::SVector{3,Float64}
    ω_lb_b::SVector{3,Float64}
    ω_el_l::SVector{3,Float64}
    ω_ie_b::SVector{3,Float64}
    ω_ib_b::SVector{3,Float64}
    v_eOb_b::SVector{3,Float64}
    v_eOb_n::SVector{3,Float64}
end

Base.@kwdef struct PVDataWGS84 <: KinematicData
    pos::PosDataWGS84
    vel::VelDataWGS84
end

Base.@kwdef struct AccDataWGS84 <: KinematicData
    α_eb_b::SVector{3,Float64}
    α_ib_b::SVector{3,Float64}
    a_eOb_b::SVector{3,Float64}
    a_iOb_b::SVector{3,Float64}
end


function XKinWGS84(init::KinInit)

    @unpack q_nb, Ob, ω_lb_b, v_eOb_b = init

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

    x = XKinWGS84(zeros(length(XKinWGS84))) #avoid infinite recursion
    x.pos.q_lb .= q_lb #assignment without colon slicing courtesy of RQuat iterability
    x.pos.q_el .= ltf(Ob) #assignment without colon slicing courtesy of RQuat iterability
    x.pos.h .= h
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

    return x

end

XKinWGS84() = XKinWGS84(KinInit())

function PVDataWGS84(x::XKinWGS84)

    #careful here: x.pos.h, x.vel.ω_eb_b and x.vel.v_eOb_b create views (this is
    #how LBV behaves by design). to copy the data, we can extract their
    #components using slices
    q_lb = RQuat(x.pos.q_lb)
    q_el = RQuat(x.pos.q_el)
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

    pos = PosDataWGS84(q_lb, q_nl, q_nb, q_eb, q_el, Ob)
    vel = VelDataWGS84(ω_eb_b, ω_lb_b, ω_el_l, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n)

    PVDataWGS84(pos, vel)

end

function AccDataWGS84(x_vel_dot::XVel, pv::PVDataWGS84)

    @unpack ω_eb_b, ω_ie_b, v_eOb_b = pv.vel
    @unpack Ob, q_eb = pv.pos

    ω_eb_b_dot = SVector{3}(x_vel_dot.ω_eb_b)
    v_eOb_b_dot = SVector{3}(x_vel_dot.v_eOb_b)

    r_eO_e = rECEF(Ob)
    r_eO_b = q_eb' * r_eO_e

    α_eb_b = ω_eb_b_dot
    α_ib_b = ω_eb_b_dot - ω_eb_b × ω_ie_b

    a_eOb_b = v_eOb_b_dot + ω_eb_b × v_eOb_b
    a_iOb_b = v_eOb_b_dot + (ω_eb_b + 2ω_ie_b) × v_eOb_b + ω_ie_b × (ω_ie_b × r_eO_b)

    AccDataWGS84(α_eb_b, α_ib_b, a_eOb_b, a_iOb_b)

end

function x_pos_dot(pv::PVDataWGS84)::XPosWGS84
    x_dot = XPosWGS84()
    x_dot.q_lb .= dt(pv.pos.q_lb, pv.vel.ω_lb_b)
    x_dot.q_el .= dt(pv.pos.q_el, pv.vel.ω_el_l)
    x_dot.h .= -pv.vel.v_eOb_n[3]
    return x_dot
end


end #module