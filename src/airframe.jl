module Airframe

using Base: Float64
using StaticArrays: SVector, SMatrix
using LinearAlgebra
using Flight.LBV
using Flight.WGS84
using Flight.Attitude
using Flight.Kinematics

export XAirframe, Wrench, MassData

@define_node XAirframe (kin = XKinWGS84,)

struct Wrench
    F::SVector{3,Float64}
    T::SVector{3,Float64}
end
Wrench(; F = zeros(3), T = zeros(3)) = Wrench(F, T)
Base.show(io::IO, wr::Wrench) = print(io, "Wrench(F = $(wr.F), T = $(wr.T))")

struct MassData
    m::Float64
    J_Ob_b::SMatrix{3, 3, Float64}
    r_ObG_b::SMatrix{3, 3, Float64}
end
MassData(; m = 1.0, J_Ob_b = SM3(I), r_ObG_b = zeros(3)) = MassData(m, J_Ob_b, r_ObG_b)

function x_vel_dot(wr_Ob_b::Wrench, mass::MassData, kin::AttPosVelWGS84)::XVel
    return XVel()
end

# function x_vel_dot(wr_Ob_b::Wrench, mass::MassData, kin::KinDataFlat)::XVel
#     return XVel()
# end



end #module