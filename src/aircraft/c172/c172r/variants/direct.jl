module C172RDirect

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI
using Flight.FlightCore.IODevices
using Flight.FlightCore.Joysticks
using Flight.FlightCore.Utils: Ranged

using Flight.FlightPhysics.Kinematics
using Flight.FlightPhysics.RigidBody
using Flight.FlightPhysics.Environment

using Flight.FlightComponents.Piston
using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.World

using ...C172
using ..C172R

export Cessna172RDirect

################################################################################
############################### DirectControls #################################

struct DirectControls <: AbstractAvionics end

@kwdef mutable struct DirectControlsU
    eng_start::Bool = false
    eng_stop::Bool = false
    throttle::Ranged{Float64, 0., 1.} = 0.0
    mixture::Ranged{Float64, 0., 1.} = 0.5
    aileron::Ranged{Float64, -1., 1.} = 0.0
    elevator::Ranged{Float64, -1., 1.} = 0.0
    rudder::Ranged{Float64, -1., 1.} = 0.0
    aileron_offset::Ranged{Float64, -1., 1.} = 0.0
    elevator_offset::Ranged{Float64, -1., 1.} = 0.0
    rudder_offset::Ranged{Float64, -1., 1.} = 0.0
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end

@kwdef struct DirectControlsY
    eng_start::Bool = false
    eng_stop::Bool = false
    throttle::Float64 = 0.0
    mixture::Float64 = 0.5
    aileron::Float64 = 0.0
    elevator::Float64 = 0.0
    rudder::Float64 = 0.0
    aileron_offset::Float64 = 0.0
    elevator_offset::Float64 = 0.0
    rudder_offset::Float64 = 0.0
    flaps::Float64 = 0.0
    brake_left::Float64 = 0.0
    brake_right::Float64 = 0.0
end

Systems.init(::SystemU, ::DirectControls) = DirectControlsU()
Systems.init(::SystemY, ::DirectControls) = DirectControlsY()


########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:DirectControls}, ::Real,
                        ::System{<:C172R.Physics},
                        ::System{<:AbstractEnvironment})

    #DirectControls has no internal dynamics, just input-output feedthrough
    @unpack eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset, flaps,
            brake_left, brake_right = avionics.u

    avionics.y = DirectControlsY(;
            eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset, flaps,
            brake_left, brake_right)

    return false

end

function Aircraft.assign!(airframe::System{<:C172.Airframe},
                                avionics::System{DirectControls})

    @unpack eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset, flaps,
            brake_left, brake_right = avionics.y

    @pack!  airframe.act.u = eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset, flaps,
            brake_left, brake_right

end


################################## GUI #########################################

function GUI.draw!(avionics::System{<:DirectControls},
                    physics::System{<:C172R.Physics},
                    label::String = "Cessna 172R Direct Controls")

    @unpack airframe = physics
    u = avionics.u

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    if airframe.y.pwp.engine.state === Piston.eng_off
        eng_start_HSV = HSV_gray
    elseif airframe.y.pwp.engine.state === Piston.eng_starting
        eng_start_HSV = HSV_amber
    else
        eng_start_HSV = HSV_green
    end
    dynamic_button("Engine Start", eng_start_HSV, 0.1, 0.2)
    u.eng_start = CImGui.IsItemActive()
    CImGui.SameLine()
    dynamic_button("Engine Stop", HSV_gray, (HSV_gray[1], HSV_gray[2], HSV_gray[3] + 0.1), (0.0, 0.8, 0.8))
    u.eng_stop = CImGui.IsItemActive()
    CImGui.SameLine()
    CImGui.Text(@sprintf("%.3f RPM", Piston.radpersec2RPM(airframe.y.pwp.engine.ω)))
    CImGui.Separator()

    u.throttle = safe_slider("Throttle", u.throttle, "%.6f")
    u.aileron = safe_slider("Aileron", u.aileron, "%.6f")
    u.elevator = safe_slider("Elevator", u.elevator, "%.6f")
    u.rudder = safe_slider("Rudder", u.rudder, "%.6f")
    u.aileron_offset = safe_input("Aileron Offset", u.aileron_offset, 0.001, 0.1, "%.6f")
    u.elevator_offset = safe_input("Elevator Offset", u.elevator_offset, 0.001, 0.1, "%.6f")
    u.rudder_offset = safe_input("Rudder Offset", u.rudder_offset, 0.001, 0.1, "%.6f")
    u.flaps = safe_slider("Flaps", u.flaps, "%.6f")
    u.mixture = safe_slider("Mixture", u.mixture, "%.6f")
    u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
    u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")


    CImGui.PopItemWidth()

    CImGui.End()

end

################################################################################
############################# Cessna172RDirect #################################

#Cessna172R with direct control Avionics
const Cessna172RDirect{K} = C172R.Template{K, DirectControls} where {K}
Cessna172RDirect(kinematics = LTF()) = C172R.Template(kinematics, DirectControls())

############################ Joystick Mappings #################################

function IODevices.assign!(sys::System{<:Cessna172RDirect}, joystick::Joystick,
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

    u.aileron_offset -= 0.01 * was_released(joystick, :dpad_left)
    u.aileron_offset += 0.01 * was_released(joystick, :dpad_right)
    u.elevator_offset += 0.01 * was_released(joystick, :dpad_down)
    u.elevator_offset -= 0.01 * was_released(joystick, :dpad_up)

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

    u.aileron_offset -= 2e-4 * is_pressed(joystick, :hat_left)
    u.aileron_offset += 2e-4 * is_pressed(joystick, :hat_right)
    u.elevator_offset += 2e-4 * is_pressed(joystick, :hat_down)
    u.elevator_offset -= 2e-4 * is_pressed(joystick, :hat_up)

    u.flaps += 0.3333 * was_released(joystick, :button_3)
    u.flaps -= 0.3333 * was_released(joystick, :button_2)

end


end #module