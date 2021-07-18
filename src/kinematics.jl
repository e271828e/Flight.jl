module Kinematics

using StaticArrays: SVector
using Flight.LBV
using Flight.WGS84
using Flight.Attitude

export XKinWGS84, KinDataWGS84, KinInit

SV3 = SVector{3, Float64}
# SV4 = SVector{4, Float64}

@define_node XVel (ω_eb_b = LBVLeaf{3}, v_eOb_b = LBVLeaf{3})
@define_node XPosWGS84 (q_e_w = LBVLeaf{4}, h = LBVLeaf{1})
@define_node XAttWGS84 (q_w_b = LBVLeaf{4},)

@define_node XKinWGS84 (vel = XVel, pos = XPosWGS84, att = XAttWGS84)

struct KinInit
    n_b::RQuat
    Ob::WGS84Pos
    ω_wb_b::SV3
    v_eOb_b::SV3
end
KinInit(; n_b = RQuat(), Ob = WGS84Pos(), ω_wb_b = zeros(3), v_eOb_b = zeros(3)) =
    KinInit(n_b, Ob, ω_wb_b, v_eOb_b)

function XKinWGS84(init::KinInit)
    print("OK, got $init")
end



struct VelDataWGS84
    ω_eb_b::SV3
    ω_wb_b::SV3
    ω_ew_w::SV3
    ω_ie_b::SV3
    ω_ib_b::SV3
    v_eOb_b::SV3
    v_eOb_n::SV3
end

struct PosDataWGS84
    q_e_w::RQuat
    Ob::WGS84Pos
end

struct AttDataWGS84
    q_w_b::RQuat
    q_n_w::RQuat
    q_n_b::RQuat
    q_e_b::RQuat
end

struct KinDataWGS84
    vel::VelDataWGS84
    pos::PosDataWGS84
    att::AttDataWGS84
end

struct AccDataWGS84
    α_eb_b::SV3
    α_ib_b::SV3
    a_eOb_b::SV3
    a_iOb_b::SV3
end

function x_pos_dot(data::KinDataWGS84)::XPosWGS84
    return XPosWGS84()
end

function x_att_dot(data::KinDataWGS84)::XAttWGS84
    return XAttWGS84()
end

end #module