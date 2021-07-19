module Kinematics

using StaticArrays: SVector
using UnPack

using Flight.LBV
using Flight.WGS84
using Flight.Attitude

export XKinWGS84, KinDataWGS84, KinInit

Base.@kwdef struct KinInit
    n_b::RQuat = RQuat()
    Ob::WGS84Pos = WGS84Pos()
    ω_lb_b::SVector{3, Float64} = zeros(3)
    v_eOb_b::SVector{3, Float64} = zeros(3)
end

@define_node XAttWGS84 (l_b = LBVLeaf{4},)
@define_node XPosWGS84 (e_l = LBVLeaf{4}, h = LBVLeaf{1})
@define_node XVel (ω_eb_b = LBVLeaf{3}, v_eOb_b = LBVLeaf{3})
@define_node XKinWGS84 (att = XAttWGS84, pos = XPosWGS84, vel = XVel)

function XKinWGS84(init::KinInit)

    x = XKinWGS84()

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

    x.att.l_b .= l_b #direct assignment possible thanks to RQuat being iterable
    x.pos.e_l .= ltf(Ob) #ditto
    x.pos.h .= h
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

    return x

end

#since this constructor is defined specifically for LBVNode{:XKinWGS84}, it
#overrides the default no-argument constructor for LBVNode{S} with generic type
#parameter
# XKinWGS84() = XKinWGS84(KinInit())

Base.@kwdef struct AccDataWGS84
    α_eb_b::SVector{3,Float64}
    α_ib_b::SVector{3,Float64}
    a_eOb_b::SVector{3,Float64}
    a_iOb_b::SVector{3,Float64}
end

function AccDataWGS84(x_vel_dot::XVel, kin::XKinWGS84)

end

Base.@kwdef struct AttDataWGS84
    l_b::RQuat
    n_l::RQuat
    n_b::RQuat
    e_b::RQuat
end

Base.@kwdef struct PosDataWGS84
    Ob::WGS84Pos
    e_l::RQuat
end

Base.@kwdef struct VelDataWGS84
    ω_eb_b::SVector{3,Float64}
    ω_lb_b::SVector{3,Float64}
    ω_el_l::SVector{3,Float64}
    ω_ie_b::SVector{3,Float64}
    ω_ib_b::SVector{3,Float64}
    v_eOb_b::SVector{3,Float64}
    v_eOb_n::SVector{3,Float64}
end

Base.@kwdef struct KinDataWGS84
    att::AttDataWGS84
    pos::PosDataWGS84
    vel::VelDataWGS84
end

function Base.show(io::IO, data::KinDataWGS84)
    println(io, "Attitude:")
    for f in fieldnames(AttDataWGS84)
        println(io, "\t", f, ": ", getproperty(data.att, f))
    end
    println(io, "Position:")
    for f in fieldnames(PosDataWGS84)
        println(io, "\t", f, ": ", getproperty(data.pos, f))
    end
    println(io, "Velocity:")
    for f in fieldnames(VelDataWGS84)
        println(io, "\t", f, ": ", getproperty(data.vel, f))
    end
end
Base.show(io::IO, ::MIME"text/plain", data::AttDataWGS84) = show(io, data)

function KinDataWGS84(x::XKinWGS84)

    #note 1: it seems that when few operations are performed with a given vector, the
    #initial cost of allocating it as StaticVectors is not worth it

    #note 2: for some reason, using the type alias SV3 = SVector{3,Float64} as SV3{} yields type
    #instability. instead, SVector(1.0,2,3) should be used

    #careful here: x.pos.h, x.vel.ω_eb_b and x.vel.v_eOb_b create views (this is
    #how LBV behaves by design). to copy the data, we can extract their
    #components using slices
    l_b = RQuat(x.att.l_b)
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

    att = AttDataWGS84(l_b, n_l, n_b, e_b)
    pos = PosDataWGS84(Ob, e_l)
    vel = VelDataWGS84(ω_eb_b, ω_lb_b, ω_el_l, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n)

    KinDataWGS84(att, pos, vel)

end

function x_pos_dot(data::KinDataWGS84)::XPosWGS84
    return XPosWGS84()
end

function x_att_dot(data::KinDataWGS84)::XAttWGS84
    return XAttWGS84()
end


end #module