module C172

using LinearAlgebra
using StaticArrays
using UnPack
using Unitful

using Flight.Modeling
import Flight.Modeling: init_x0, init_y0, init_u0, init_d0, f_cont!, f_disc!
using Flight.Plotting
import Flight.Plotting: plots

using Flight.Attitude

using Flight.Terrain
using Flight.Airdata
using Flight.Kinematics
using Flight.Dynamics

using Flight.Components
import Flight.Components: MassTrait, WrenchTrait, AngularMomentumTrait, get_mass_properties

using Flight.Propulsion
using Flight.Aerodynamics
using Flight.LandingGear

using Flight.Aircraft

export C172Aircraft

struct C172ID <: AbstractAircraftID end


############################## C172Powerplant ################################

Base.@kwdef struct C172Pwp <: SystemGroupDescriptor
    left::EThruster = EThruster(motor = ElectricMotor(α = CW))
    right::EThruster = EThruster(motor = ElectricMotor(α = CCW))
end

WrenchTrait(::System{<:C172Pwp}) = HasWrench()
AngularMomentumTrait(::System{<:C172Pwp}) = HasAngularMomentum()

############################## C172Airframe #################################

Base.@kwdef struct C172Airframe{ Aero <: SystemDescriptor, Pwp <: SystemDescriptor,
                        Ldg <: SystemDescriptor} <: SystemGroupDescriptor
    aero::Aero = C172_aero()
    pwp::Pwp = C172_pwp()
    ldg::Ldg = C172_ldg()
end

MassTrait(::System{<:C172Airframe}) = HasMass()
WrenchTrait(::System{<:C172Airframe}) = HasWrench()
AngularMomentumTrait(::System{<:C172Airframe}) = HasAngularMomentum()

function get_mass_properties(::System{<:C172Airframe})

    #for an aircraft implementing a fuel system here we could compute the actual
    #mass properties from the different fuel tank contents

    MassProperties(
        #upreferred(2650u"lb") |> ustrip,
        m = 1202.0197805,
        #upreferred.([948, 1346, 1967]u"m") |> ustrip |> diagm |> SMatrix{3,3,Float64},
        J_Ob_b = SA[948.0 0 0; 0 1346.0 0; 0 0 1967.0],
        r_ObG_b = zeros(SVector{3}))
end

#get_wr_b and get_hr_b use the fallback for SystemGroups, which will in turn
#call get_wr_b and get_hr_b on aero, pwp and ldg

#uses default SystemGroup f_disc! implementation

############################## C172Controls #################################

struct C172Controls <: SystemDescriptor end

Base.@kwdef mutable struct C172ControlsU
    throttle::Float64 = 0.0 #[0, 1]
    yoke_x::Float64 = 0.0 #[-1, 1], ailerons (+ bank right)
    yoke_y::Float64 = 0.0 #[-1, 1], elevator (+ pitch up)
    pedals::Float64 = 0.0 #[-1, 1], rudder and nose wheel (+ yaw right)
    brake_left::Float64 = 0.0 #[0, 1]
    brake_right::Float64 = 0.0 #[0, 1]
end

#const is essential when declaring type aliases!
struct C172ControlsY
    throttle::Float64
    yoke_x::Float64
    yoke_y::Float64
    pedals::Float64
    brake_left::Float64
    brake_right::Float64
end

init_u0(::C172Controls) = C172ControlsU()

init_y0(::C172Controls) = C172ControlsY(zeros(SVector{6})...)


###################### Continuous update functions ##########################

function assign_component_inputs!(afr::System{<:C172Airframe},
    ctl::System{<:C172Controls})

    @unpack throttle, yoke_x, yoke_y, pedals, brake_left, brake_right = ctl.u
    @unpack aero, pwp, ldg = afr.subsystems

    pwp.u.left.throttle = throttle
    pwp.u.right.throttle = throttle
    ldg.u.center.steering[] = pedals
    ldg.u.left.braking[] = brake_left
    ldg.u.right.braking[] = brake_right

    return nothing
end

function f_cont!(afr::System{<:C172Airframe}, ctl::System{<:C172Controls},
                kin::KinData, air::AirData, trn::AbstractTerrain)

    @unpack aero, pwp, ldg = afr.subsystems

    #could this go in the main aircraft f_cont!?
    assign_component_inputs!(afr, ctl)
    # f_cont!(srf, air) #update surface actuator continuous state & outputs
    f_cont!(ldg, kin, trn) #update landing gear continuous state & outputs
    f_cont!(pwp, kin, air) #update powerplant continuous state & outputs
    f_cont!(aero, air, kin, trn) #requires previous srf update
    # f_cont!(aero, air, kin, srf, trn) #requires previous srf update

    afr.y = (aero = aero.y, pwp = pwp.y, ldg = ldg.y)

end

function f_cont!(ctl::System{<:C172Controls}, ::System{<:C172Airframe},
                ::KinData, ::AirData, ::AbstractTerrain)

    #here, controls do nothing but update their output state. in a more complex
    #aircraft this could even be a complete autopilot implementation
    @unpack throttle, yoke_x, yoke_y, pedals, brake_left, brake_right = ctl.u
    return C172ControlsY(throttle, yoke_x, yoke_y, pedals, brake_left,
    brake_right)
    # return (throttle = 0.0, yoke_x = 0.0, yoke_y = 0.0, pedals = 0.0,
    #         brake_left = 0.0, brake_right = 0.0)

end

######################## Discrete Update Functions ########################

function f_disc!(::System{<:C172Controls}, ::System{<:C172Airframe})
    #this does nothing, but could also be used to implement open loop or closed
    #loop control laws, predefined maneuvers, etc. by calling an additional
    #function we could call as part of
    return false
end

function f_disc!(afr::System{<:C172Airframe}, ::System{<:C172Controls})
    #fall back to the default SystemGroup implementation, the f_disc! for the
    #C172 components don't have to deal with the C172Controls
    return f_disc!(afr)
end

#these should eventually be declared as constants:

C172_kin() = KinLTF()

C172_aero() = SimpleDrag()

C172_pwp() = C172Pwp()

function C172_ldg()

    mlg_damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)
    nlg_damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)

    left = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-1, -1.25, 1], q = RQuat() ),
            l_0 = 0.0,
            damper = mlg_damper),
        braking = DirectBraking())

    right = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-1, 1.25, 1], q = RQuat() ),
            l_0 = 0.0,
            damper = mlg_damper),
        braking = DirectBraking())

    center = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [2, 0, 1] , q = RQuat()),
            l_0 = 0.0,
            damper = nlg_damper),
        steering = DirectSteering())

    TricycleLandingGear(; left, right, center)

end

function C172Aircraft(; kin = C172_kin(), aero = C172_aero(), pwp = C172_pwp(),
                        ldg = C172_ldg())
    AircraftBase(
        C172ID(),
        kin,
        C172Airframe(aero, pwp, ldg),
        C172Controls())
end

end #module