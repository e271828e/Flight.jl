module C172

using LinearAlgebra
using StaticArrays
using UnPack
using Unitful

using Flight.Modeling
import Flight.Modeling: init_x0, init_y0, init_u0, init_d0, f_cont!, f_disc!
using Flight.Plotting
import Flight.Plotting: plots
using Flight.Misc

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
import Flight.Aircraft: assign_joystick_inputs!

using Flight.Input

export C172Aircraft, C172Pwp, C172Ldg

struct C172ID <: AbstractAircraftID end


############################## C172Powerplant ################################

struct C172Pwp <: SystemGroupDescriptor
    left::EThruster
    right::EThruster
end

WrenchTrait(::System{<:C172Pwp}) = HasWrench()
AngularMomentumTrait(::System{<:C172Pwp}) = HasAngularMomentum()

function C172Pwp()

    prop = SimpleProp(kF = 4e-3, J = 0.25)

    left = EThruster(propeller = prop, motor = ElectricMotor(α = CW))
    right = EThruster(propeller = prop, motor = ElectricMotor(α = CCW))

    C172Pwp(left, right)

end

############################ C172LandingGear ############################

struct C172Ldg{L <: LandingGearUnit, R <: LandingGearUnit,
    C <: LandingGearUnit} <: SystemGroupDescriptor
    left::L
    right::R
    center::C
end

WrenchTrait(::System{<:C172Ldg}) = HasWrench()
AngularMomentumTrait(::System{<:C172Ldg}) = HasNoAngularMomentum()

function C172Ldg()

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

    C172Ldg(left, right, center)

end

############################### C172 Aerodynamics ###########################

C172_aero() = SimpleDrag()

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
    throttle::Bounded{Float64, 0, 1} = 0.0
    yoke_x::Bounded{Float64, -1, 1} = 0.0 #ailerons (+ bank right)
    yoke_y::Bounded{Float64, -1, 1} = 0.0 #elevator (+ pitch up)
    pedals::Bounded{Float64, -1, 1} = 0.0 #rudder and nose wheel (+ yaw right)
    brake_left::Bounded{Float64, 0, 1} = 0.0 #[0, 1]
    brake_right::Bounded{Float64, 0, 1} = 0.0 #[0, 1]
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


function C172Aircraft(; kin = C172_kin(), aero = C172_aero(), pwp = C172Pwp(),
                        ldg = C172Ldg())
    AircraftBase(
        C172ID(),
        kin,
        C172Airframe(aero, pwp, ldg),
        C172Controls())
end

######################## Joystick Input Interface ########################

yoke_curve(x) = exp_axis_curve(x, strength = 0.5, deadzone = 0.05)
pedal_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.1)
brake_curve(x) = exp_axis_curve(x, strength = 0, deadzone = 0.05)

function assign_joystick_inputs!(ac::System{<:AircraftBase{C172ID}}, joystick::XBoxController)

    if get_button_change(joystick, :dpad_up) === button_released
        ac.u.throttle = Float64(ac.u.throttle) + 0.1
    elseif get_button_change(joystick, :dpad_down) === button_released
        ac.u.throttle = Float64(ac.u.throttle) - 0.1
    end

    ac.u.yoke_x = get_axis_data(joystick, :right_analog_x) |> yoke_curve
    ac.u.yoke_y = get_axis_data(joystick, :right_analog_y) |> yoke_curve
    ac.u.pedals = get_axis_data(joystick, :left_analog_x) |> pedal_curve
    ac.u.brake_left = get_axis_data(joystick, :left_trigger) |> brake_curve
    ac.u.brake_right = get_axis_data(joystick, :right_trigger) |> brake_curve
end


end #module