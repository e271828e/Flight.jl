module Airframe

using Base: Float64
using StaticArrays: SVector, SMatrix
using LinearAlgebra
using Flight.LBV
using Flight.WGS84
using Flight.Attitude
using Flight.Kinematics

export XAirframe, Wrench, MassData

SV3 = SVector{3, Float64}
SM3 = SMatrix{3, 3, Float64}

@define_node XAirframe (kin = XKinWGS84,)

struct Wrench
    F::SV3
    T::SV3
end
Wrench(; F = zeros(3), T = zeros(3)) = Wrench(F, T)
Base.show(io::IO, wr::Wrench) = print(io, "Wrench(F = $(wr.F), T = $(wr.T))")

struct MassData
    m::Float64
    J_Ob_b::SM3
    r_ObG_b::SV3
end
MassData(; m = 1.0, J_Ob_b = SM3(I), r_ObG_b = zeros(3)) = MassData(m, J_Ob_b, r_ObG_b)

function x_vel_dot(wr_Ob_b::Wrench, mass::MassData, kin::KinDataWGS84)::XVel
    return XVel()
end

# function x_vel_dot(wr_Ob_b::Wrench, mass::MassData, kin::KinDataFlat)::XVel
#     return XVel()
# end



end #module