module C172RDirect

using UnPack
using Printf
using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic

using Flight.FlightCore
using Flight.FlightPhysics
using Flight.FlightAircraft

include("airframe.jl"); using .Airframe

export Cessna172R

################################################################################
############################### DirectControls #################################

struct DirectControls <: AbstractAvionics end

Base.@kwdef mutable struct DirectControlsU
    eng_start::Bool = false
    eng_stop::Bool = false
    throttle::Ranged{Float64, 0, 1} = 0.0
    mixture::Ranged{Float64, 0, 1} = 0.5
    aileron::Ranged{Float64, -1, 1} = 0.0
    elevator::Ranged{Float64, -1, 1} = 0.0
    rudder::Ranged{Float64, -1, 1} = 0.0
    aileron_trim::Ranged{Float64, -1, 1} = 0.0
    elevator_trim::Ranged{Float64, -1, 1} = 0.0
    rudder_trim::Ranged{Float64, -1, 1} = 0.0
    flaps::Ranged{Float64, 0, 1} = 0.0
    brake_left::Ranged{Float64, 0, 1} = 0.0
    brake_right::Ranged{Float64, 0, 1} = 0.0
end

Base.@kwdef struct DirectControlsY
    eng_start::Bool = false
    eng_stop::Bool = false
    throttle::Float64 = 0.0
    mixture::Float64 = 0.5
    aileron::Float64 = 0.0
    elevator::Float64 = 0.0
    rudder::Float64 = 0.0
    aileron_trim::Float64 = 0.0
    elevator_trim::Float64 = 0.0
    rudder_trim::Float64 = 0.0
    flaps::Float64 = 0.0
    brake_left::Float64 = 0.0
    brake_right::Float64 = 0.0
end

Systems.init(::SystemU, ::DirectControls) = DirectControlsU()
Systems.init(::SystemY, ::DirectControls) = DirectControlsY()

############################ hasta aqui ########################################

########################### Update Methods #####################################

function Systems.f_ode!(avionics::System{DirectControls}, ::System{<:Airframe},
                ::KinematicData, ::AirData, ::RigidBodyData,
                ::System{<:AbstractTerrain})

    #DirectControls has no internal dynamics, just input-output feedthrough
    @unpack eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right = avionics.u

    avionics.y = DirectControlsY(;
            eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right)

end

function Aircraft.map_controls!(airframe::System{<:Airframe}, avionics::System{DirectControls})

    @unpack eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right = avionics.u

    @pack!  airframe.act.u =
            eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right

end


############################ Joystick Mappings #################################

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign!(u::DirectControlsU,
            joystick::Joystick{XBoxControllerID}, ::DefaultMapping)

    u.aileron = get_axis_value(joystick, :right_analog_x) |> aileron_curve
    u.elevator = get_axis_value(joystick, :right_analog_y) |> elevator_curve
    u.rudder = get_axis_value(joystick, :left_analog_x) |> rudder_curve
    u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

    u.aileron_trim -= 0.01 * was_released(joystick, :dpad_left)
    u.aileron_trim += 0.01 * was_released(joystick, :dpad_right)
    u.elevator_trim += 0.01 * was_released(joystick, :dpad_down)
    u.elevator_trim -= 0.01 * was_released(joystick, :dpad_up)

    u.throttle += 0.1 * was_released(joystick, :button_Y)
    u.throttle -= 0.1 * was_released(joystick, :button_A)

    u.flaps += 0.3333 * was_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

end


################################## GUI #########################################


function GUI.draw!(sys::System{<:DirectControls}, label::String = "Cessna 172R Direct Controls")

    u = sys.u

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    u.eng_start = dynamic_button("Engine Start", 0.4); CImGui.SameLine()
    u.eng_stop = dynamic_button("Engine Stop", 0.0)
    u.throttle = safe_slider("Throttle", u.throttle, "%.6f")
    u.mixture = safe_slider("Mixture", u.mixture, "%.6f")
    u.aileron = safe_slider("Aileron", u.aileron, "%.6f")
    u.elevator = safe_slider("Elevator", u.elevator, "%.6f")
    u.rudder = safe_slider("Rudder", u.rudder, "%.6f")
    u.aileron_trim = safe_input("Aileron Trim", u.aileron_trim, 0.001, 0.1, "%.6f")
    u.elevator_trim = safe_input("Elevator Trim", u.elevator_trim, 0.001, 0.1, "%.6f")
    u.rudder_trim = safe_input("Rudder Trim", u.rudder_trim, 0.001, 0.1, "%.6f")
    u.flaps = safe_slider("Flaps", u.flaps, "%.6f")
    u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
    u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")

    CImGui.PopItemWidth()

    CImGui.End()

    GUI.draw(sys, label)

end


################################################################################
############################### Cessna172R #####################################

#here we define the basic Cessna172R, which installs the C172RDirect.DirectControls
#avionics, which provides a basic reversible direct control system
const Cessna172RDirect{K, F} = AircraftTemplate{K, F, DirectControls} where {K, F <: Airframe}
Cessna172RDirect(kinematics = LTF()) = AircraftTemplate(kinematics, C172RAirframe(), DirectControls())



end #module