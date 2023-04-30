module C172Rv1

using UnPack
using Printf
using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic

using LinearAlgebra
using StaticArrays
using ComponentArrays

using Flight.FlightCore
using Flight.FlightPhysics
using Flight.FlightAircraft

using ..Airframe

export Cessna172Rv1

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


########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:DirectControls}, ::Real,
                        ::System{<:C172RAirframe}, ::KinematicData,
                        ::RigidBodyData, ::AirData, ::TerrainData)

    #DirectControls has no internal dynamics, just input-output feedthrough
    @unpack eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right = avionics.u

    avionics.y = DirectControlsY(;
            eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right)

    return false

end

function Aircraft.map_controls!(airframe::System{<:C172RAirframe},
                                avionics::System{DirectControls})

    @unpack eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right = avionics.y

    @pack!  airframe.u.act = eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right

end


################################## GUI #########################################

function GUI.draw!(avionics::System{<:DirectControls}, airframe::System{<:C172RAirframe},
                    label::String = "Cessna 172R Direct Controls")

    u = avionics.u

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    u.eng_start = dynamic_button("Engine Start", 0.4); CImGui.SameLine()
    u.eng_stop = dynamic_button("Engine Stop", 0.0)
    u.throttle = safe_slider("Throttle", u.throttle, "%.6f")
    u.aileron = safe_slider("Aileron", u.aileron, "%.6f")
    u.elevator = safe_slider("Elevator", u.elevator, "%.6f")
    u.rudder = safe_slider("Rudder", u.rudder, "%.6f")
    u.aileron_trim = safe_input("Aileron Trim", u.aileron_trim, 0.001, 0.1, "%.6f")
    u.elevator_trim = safe_input("Elevator Trim", u.elevator_trim, 0.001, 0.1, "%.6f")
    u.rudder_trim = safe_input("Rudder Trim", u.rudder_trim, 0.001, 0.1, "%.6f")
    u.flaps = safe_slider("Flaps", u.flaps, "%.6f")
    u.mixture = safe_slider("Mixture", u.mixture, "%.6f")
    u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
    u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")

    CImGui.Text(@sprintf("Engine Speed: %.3f RPM",
                        Piston.radpersec2RPM(airframe.y.pwp.engine.ω)))

    CImGui.PopItemWidth()

    CImGui.End()

end

################################################################################
############################# Cessna172Rv1 #####################################

#Cessna172R variant with DirectControls avionics
const Cessna172Rv1{K, F} = AircraftTemplate{K, F, DirectControls} where {K, F <: C172RAirframe}
Cessna172Rv1(kinematics = LTF()) = AircraftTemplate(kinematics, C172RAirframe(), DirectControls())


############################ Joystick Mappings #################################

function IODevices.assign!(sys::System{<:Cessna172Rv1}, joystick::Joystick,
                           mapping::InputMapping)
    IODevices.assign!(sys.avionics, joystick, mapping)
end

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign!(sys::System{DirectControls},
                           joystick::XBoxController,
                           ::DefaultMapping)

    u = sys.u

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

function IODevices.assign!(sys::System{DirectControls},
                           joystick::T16000M,
                           ::DefaultMapping)

    u = sys.u

    u.throttle = get_axis_value(joystick, :throttle)
    u.aileron = get_axis_value(joystick, :stick_x) |> aileron_curve
    u.elevator = get_axis_value(joystick, :stick_y) |> elevator_curve
    u.rudder = get_axis_value(joystick, :stick_z) |> rudder_curve

    u.brake_left = is_pressed(joystick, :button_1)
    u.brake_right = is_pressed(joystick, :button_1)

    u.aileron_trim -= 2e-4 * is_pressed(joystick, :hat_left)
    u.aileron_trim += 2e-4 * is_pressed(joystick, :hat_right)
    u.elevator_trim += 2e-4 * is_pressed(joystick, :hat_down)
    u.elevator_trim -= 2e-4 * is_pressed(joystick, :hat_up)

    u.flaps += 0.3333 * was_released(joystick, :button_3)
    u.flaps -= 0.3333 * was_released(joystick, :button_2)

end


end #module