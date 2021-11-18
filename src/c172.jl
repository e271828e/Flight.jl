module C172

using LinearAlgebra
using StaticArrays
using UnPack
using Unitful

using Flight.ModelingTools

using Flight.Attitude

using Flight.Airframe
using Flight.Mass
# using Flight.Airdata
# using Flight.Propulsion
# using Flight.Aerodynamics
using Flight.LandingGear

using Flight.Kinematics
using Flight.Dynamics

using Flight.Aircraft

# using Flight.Plotting
# import Flight.Plotting: plots

export C172Aircraft

struct C172ID <: AbstractAircraftID end

const C172_kin = KinLTF()

const C172_mass = ConstantMass(
    m = upreferred(2650u"lb") |> ustrip,
    J_Ob_b = upreferred.([948, 1346, 1967]u"m") |> ustrip |> diagm |> SMatrix{3,3,Float64},
    r_ObG_b = zeros(SVector{3}))

const C172_aero = NullAirframeComponent()
const C172_pwp = NullAirframeComponent()
const C172_ldg = TricycleLandingGear(

    left = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-1, -1.25, 1], q = RQuat() ),
            l_0 = 0.0,
            damper = SimpleDamper(
                k_s = 25000,
                k_d_ext = 1000,
                k_d_cmp = 1000,
                ξ_min = -1)
            ),
        braking = DirectBraking()),

    right = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-1, 1.25, 1], q = RQuat() ),
            l_0 = 0.0,
            damper = SimpleDamper(
                k_s = 25000,
                k_d_ext = 1000,
                k_d_cmp = 1000,
                ξ_min = -1)
            ),
        braking = DirectBraking()),

    center = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [2, 0, 1] , q = RQuat()),
            l_0 = 0.0),
        steering = DirectSteering())
)

const C172_srf = NullAirframeComponent()

function C172Aircraft(; kin = C172_kin, mass = C172_mass, aero = C172_aero,
                pwp = C172_pwp, ldg = C172_ldg, srf = C172_srf)
    AircraftBase(C172ID(), kin, mass, aero, pwp, ldg, srf)
end

end #module