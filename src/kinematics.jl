module Kinematics

using LinearAlgebra
using StaticArrays: SVector
using UnPack

using Flight.LBV
using Flight.WGS84
using Flight.Attitude

export XKinWGS84, KinInit, VelDataWGS84, PosDataWGS84, PVDataWGS84, AccDataWGS84
export x_pos_dot

#for some reason, defining and using using the type alias SV3 =
#SVector{3,Float64} yields type instability

@define_node XPosWGS84 (l_b = LBVLeaf{4}, e_l = LBVLeaf{4}, h = LBVLeaf{1})
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
    n_b::RQuat = RQuat()
    Ob::WGS84Pos = WGS84Pos()
    ω_lb_b::SVector{3, Float64} = zeros(SVector{3})
    v_eOb_b::SVector{3, Float64} = zeros(SVector{3})
end

function XKinWGS84(init::KinInit)

    @unpack n_b, Ob, ω_lb_b, v_eOb_b = init

    h = Ob.h[1]
    (R_N, R_E) = radii(Ob)
    v_eOb_n = n_b * v_eOb_b
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + h),
        -v_eOb_n[1] / (R_N + h),
        0.0)

    ω_el_b = n_b' * ω_el_n
    ω_eb_b = ω_el_b + ω_lb_b

    l_b = n_b #arbitrarily initialize ψ_nl to -1

    x = XKinWGS84(zeros(length(XKinWGS84))) #avoid infinite recursion
    x.pos.l_b .= l_b #assignment without colon slicing courtesy of RQuat iterability
    x.pos.e_l .= ltf(Ob) #assignment without colon slicing courtesy of RQuat iterability
    x.pos.h .= h
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

    return x

end

XKinWGS84() = XKinWGS84(KinInit())

Base.@kwdef struct PosDataWGS84 <: KinematicData
    l_b::RQuat
    n_l::RQuat
    n_b::RQuat
    e_b::RQuat
    e_l::RQuat
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

function PVDataWGS84(x::XKinWGS84)

    #careful here: x.pos.h, x.vel.ω_eb_b and x.vel.v_eOb_b create views (this is
    #how LBV behaves by design). to copy the data, we can extract their
    #components using slices
    l_b = RQuat(x.pos.l_b)
    e_l = RQuat(x.pos.e_l)
    h = x.pos.h[1]
    ω_eb_b = SVector{3}(x.vel.ω_eb_b)
    v_eOb_b = SVector{3}(x.vel.v_eOb_b)

    Ob = WGS84Pos(NVector(e_l), h)
    n_l = Rz(ψ_nl(e_l))
    n_b = n_l ∘ l_b
    e_b = e_l ∘ l_b

    (R_N, R_E) = radii(Ob)
    v_eOb_n = n_b * v_eOb_b
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + h),
        -v_eOb_n[1] / (R_N + h),
        0.0)

    ω_el_l = n_l' * ω_el_n
    ω_el_b = l_b' * ω_el_l
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector{3}(0, 0, ω_ie)
    ω_ie_b = e_b' * ω_ie_e
    ω_ib_b = ω_ie_b + ω_eb_b

    pos = PosDataWGS84(l_b, n_l, n_b, e_b, e_l, Ob)
    vel = VelDataWGS84(ω_eb_b, ω_lb_b, ω_el_l, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n)

    PVDataWGS84(pos, vel)

end

function x_pos_dot(pv::PVDataWGS84)::XPosWGS84
    x_dot = XPosWGS84()
    x_dot.l_b .= dt(pv.pos.l_b, pv.vel.ω_lb_b)
    x_dot.e_l .= dt(pv.pos.e_l, pv.vel.ω_el_l)
    x_dot.h .= -pv.vel.v_eOb_n[3]
    return x_dot
end

Base.@kwdef struct AccDataWGS84 <: KinematicData
    α_eb_b::SVector{3,Float64}
    α_ib_b::SVector{3,Float64}
    a_eOb_b::SVector{3,Float64}
    a_iOb_b::SVector{3,Float64}
end

function AccDataWGS84(x_vel_dot::XVel, pv::PVDataWGS84)

    @unpack ω_eb_b, ω_ie_b, v_eOb_b = pv.vel
    @unpack Ob, e_b = pv.pos

    ω_eb_b_dot = SVector{3}(x_vel_dot.ω_eb_b)
    v_eOb_b_dot = SVector{3}(x_vel_dot.v_eOb_b)

    r_eO_e = rECEF(Ob)
    r_eO_b = e_b' * r_eO_e

    α_eb_b = ω_eb_b_dot
    α_ib_b = ω_eb_b_dot - ω_eb_b × ω_ie_b

    a_eOb_b = v_eOb_b_dot + ω_eb_b × v_eOb_b
    a_iOb_b = v_eOb_b_dot + (ω_eb_b + 2ω_ie_b) × v_eOb_b + ω_ie_b × (ω_ie_b × r_eO_b)

    AccDataWGS84(α_eb_b, α_ib_b, a_eOb_b, a_iOb_b)

end


end #module