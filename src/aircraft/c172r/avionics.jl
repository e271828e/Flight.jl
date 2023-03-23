module C172RAvionics

using UnPack
using Printf
using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic

using Flight.FlightCore
using Flight.FlightPhysics

using Flight.FlightAircraft.Aircraft

using ..C172RAirframe

export ReversibleControls

################################################################################
########################### ReversibleControls #################################

struct ReversibleControls <: AbstractAvionics end

#elevator↑ (stick forward) -> e↑ -> δe↑ -> trailing edge down -> Cm↓ -> pitch down
#aileron↑ (stick right) -> a↑ -> δa↑ -> left trailing edge down, right up -> Cl↓ -> roll right
#rudder↑ (right pedal forward) -> r↓ -> δr↓ -> rudder trailing edge right -> Cn↑ -> yaw right
#rudder↑ (right pedal forward) -> nose wheel steering right -> yaw right
#flaps↑ -> δf↑ -> flap trailing edge down -> CL↑
Base.@kwdef mutable struct ReversibleControlsU
    throttle::Ranged{Float64, 0, 1} = 0.0
    aileron_trim::Ranged{Float64, -1, 1} = 0.0
    aileron::Ranged{Float64, -1, 1} = 0.0 #incremental command, for input devices
    elevator_trim::Ranged{Float64, -1, 1} = 0.0
    elevator::Ranged{Float64, -1, 1} = 0.0 #incremental command, for input devices
    rudder_trim::Ranged{Float64, -1, 1} = 0.0
    rudder::Ranged{Float64, -1, 1} = 0.0 #incremental command, for input devices
    brake_left::Ranged{Float64, 0, 1} = 0.0
    brake_right::Ranged{Float64, 0, 1} = 0.0
    flaps::Ranged{Float64, 0, 1} = 0.0
    mixture::Ranged{Float64, 0, 1} = 0.5
    eng_start::Bool = false
    eng_stop::Bool = false
end

Base.@kwdef struct ReversibleControlsY
    throttle::Float64 = 0.0
    aileron_trim::Float64 = 0.0
    aileron::Float64 = 0.0
    elevator_trim::Float64 = 0.0
    elevator::Float64 = 0.0
    rudder_trim::Float64 = 0.0
    rudder::Float64 = 0.0
    brake_left::Float64 = 0.0
    brake_right::Float64 = 0.0
    flaps::Float64 = 0.0
    mixture::Float64 = 0.5
    eng_start::Bool = false
    eng_stop::Bool = false
end

Systems.init(::SystemU, ::ReversibleControls) = ReversibleControlsU()
Systems.init(::SystemY, ::ReversibleControls) = ReversibleControlsY()


########################### Update Methods #####################################

function Systems.f_ode!(avionics::System{ReversibleControls}, ::System{<:Airframe},
                ::KinematicData, ::AirData, ::System{<:AbstractTerrain})

    #ReversibleControls has no internal dynamics, just input-output feedthrough
    @unpack throttle, aileron_trim, aileron, elevator_trim, elevator,
            rudder_trim, rudder, brake_left, brake_right, flaps, mixture,
            eng_start, eng_stop = avionics.u

    avionics.y = ReversibleControlsY(;
            throttle, aileron_trim, aileron, elevator_trim, elevator,
            rudder_trim, rudder, brake_left, brake_right, flaps, mixture,
            eng_start, eng_stop)

end

#no digital components or state machines in ReversibleControls
@inline Systems.f_step!(::System{ReversibleControls}, ::System{<:Airframe}, ::KinematicSystem) = false
@inline Systems.f_disc!(::System{ReversibleControls}, ::System{<:Airframe}, ::KinematicSystem, Δt) = false


function Aircraft.map_controls!(airframe::System{<:Airframe}, avionics::System{ReversibleControls})

    @unpack throttle, aileron_trim, aileron, elevator_trim, elevator,
            rudder_trim, rudder, brake_left, brake_right, flaps, mixture,
            eng_start, eng_stop = avionics.u

    @unpack aero, pwp, ldg = airframe

    pwp.u.engine.start = eng_start
    pwp.u.engine.stop = eng_stop
    pwp.u.engine.thr = throttle
    pwp.u.engine.mix = mixture
    ldg.u.nose.steering[] = (rudder_trim + rudder) #rudder↑ (right pedal forward) -> nose wheel steering right
    ldg.u.left.braking[] = brake_left
    ldg.u.right.braking[] = brake_right
    aero.u.e = (elevator_trim + elevator) #elevator↑ (stick forward) -> e↑ -> pitch down
    aero.u.a = (aileron_trim + aileron) #aileron↑ (stick right) -> a↑ -> roll right
    aero.u.r = -(rudder_trim + rudder) #rudder↑ (right pedal forward) -> r↓ -> yaw right
    aero.u.f = flaps #flaps↑ -> δf↑

    return nothing
end


############################ Joystick Mappings #################################

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign!(u::ReversibleControlsU,
            joystick::Joystick{XBoxControllerID}, ::DefaultMapping)

    u.aileron = get_axis_value(joystick, :right_analog_x) |> aileron_curve
    u.elevator = -get_axis_value(joystick, :right_analog_y) |> elevator_curve
    u.rudder = get_axis_value(joystick, :left_analog_x) |> rudder_curve
    u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

    u.aileron_trim -= 0.01 * was_released(joystick, :dpad_left)
    u.aileron_trim += 0.01 * was_released(joystick, :dpad_right)
    u.elevator_trim -= 0.01 * was_released(joystick, :dpad_down)
    u.elevator_trim += 0.01 * was_released(joystick, :dpad_up)

    u.throttle += 0.1 * was_released(joystick, :button_Y)
    u.throttle -= 0.1 * was_released(joystick, :button_A)

    u.flaps += 0.3333 * was_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

end


################################## GUI #########################################


function GUI.draw!(sys::System{<:ReversibleControls}, label::String = "Cessna 172R Reversible Controls")

    u = sys.u

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    u.eng_start = dynamic_button("Engine Start", 0.4); CImGui.SameLine()
    u.eng_stop = dynamic_button("Engine Stop", 0.0)
    u.throttle = safe_slider("Throttle", u.throttle, "%.6f")
    u.mixture = safe_slider("Mixture", u.mixture, "%.6f")
    u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
    u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")
    u.aileron = safe_slider("Aileron", u.aileron, "%.6f")
    u.elevator = safe_slider("Elevator", u.elevator, "%.6f")
    u.rudder = safe_slider("Rudder", u.rudder, "%.6f")
    u.flaps = safe_slider("Flaps", u.flaps, "%.6f")
    u.aileron_trim = safe_input("Aileron Trim", u.aileron_trim, 0.001, 0.1, "%.6f")
    u.elevator_trim = safe_input("Elevator Trim", u.elevator_trim, 0.001, 0.1, "%.6f")
    u.rudder_trim = safe_input("Rudder Trim", u.rudder_trim, 0.001, 0.1, "%.6f")

    CImGui.PopItemWidth()

    CImGui.End()

    GUI.draw(sys, label)

end

function GUI.draw(sys::System{<:ReversibleControls}, label::String = "Cessna 172R Reversible Controls")

    y = sys.y

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    @running_plot("Throttle", y.throttle, 0, 1, 0.0, 60)
    @running_plot("Mixture", y.mixture, 0, 1, 0.5, 60)
    @running_plot("Left Brake", y.brake_left, 0, 1, 0.0, 60)
    @running_plot("Right Brake", y.brake_right, 0, 1, 0.0, 60)
    @running_plot("Aileron", y.aileron, -1, 1, 0.0, 60)
    @running_plot("Elevator", y.elevator, -1, 1, 0.0, 60)
    @running_plot("Rudder", y.rudder, -1, 1, 0.0, 60)
    @running_plot("Aileron Trim", y.aileron_trim, -1, 1, 0.0, 60)
    @running_plot("Elevator Trim", y.elevator_trim, -1, 1, 0.0, 60)
    @running_plot("Rudder Trim", y.rudder_trim, -1, 1, 0.0, 60)
    @running_plot("Flap Setting", y.flaps, 0, 1, 0.0, 60)

    CImGui.PopItemWidth()

    CImGui.End()

end


end #module