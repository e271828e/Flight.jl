module C172FBW

using LinearAlgebra, StaticArrays, ComponentArrays
using UnPack, Reexport

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI
using Flight.FlightCore.IODevices
using Flight.FlightCore.Joysticks
using Flight.FlightCore.Utils: Ranged

using Flight.FlightPhysics.Kinematics
using Flight.FlightPhysics.RigidBody
using Flight.FlightPhysics.Environment

using Flight.FlightComponents.Propellers
using Flight.FlightComponents.Piston
using Flight.FlightComponents.Aircraft

using ..C172

################################################################################
################################ Powerplant ####################################

PowerPlant() = Piston.Thruster(propeller = Propeller(t_bp = FrameTransform(r = [2.055, 0, 0.833])))

################################################################################
################################## Actuator ####################################

@kwdef struct Actuator <: SystemDefinition #second order linear actuator model
    ω_n::Float64 = 10*2π #natural frequency (default 10 Hz)
    ζ::Float64 = 0.5 #damping ratio (default underdamped with minimal resonance)
    range::Tuple{Float64, Float64} = (-1.0, 1.0)
end

@kwdef struct ActuatorY
    cmd::Float64 = 0.0
    pos_free::Float64 = 0.0
    pos::Float64 = 0.0
    sat::Int64 = 0 #output saturation status
end

#with an underdamped actuator, the position state can still transiently exceed
#the intended range due to overshoot. the true actuator position should
#therefore be clamped. in the real world, this behaviour could correspond to a
#clutched output actuator, where the output position saturates beyond a given
#opposing torque (for example, if the surface's mechanical limits are hit)
Systems.init(::SystemU, act::Actuator) = Ref(Ranged(0.0, act.range[1], act.range[2]))
Systems.init(::SystemX, ::Actuator) = ComponentVector(v = 0.0, p = 0.0)
Systems.init(::SystemY, ::Actuator) = ActuatorY()

function Systems.f_ode!(sys::System{Actuator})

    @unpack ẋ, x, u, params = sys
    @unpack ω_n, ζ, range = params

    cmd = Float64(u[])
    pos_free = x.p
    pos = clamp(pos_free, range[1], range[2]) #clamped output
    sat_hi = pos_free >= range[1]
    sat_lo = pos_free <= range[2]
    sat = sat_hi - sat_lo

    ẋ.v = ω_n^2 * (cmd - x.p) - 2ζ*ω_n*x.v
    ẋ.p = x.v

    sys.y = ActuatorY(; cmd, pos_free, pos, sat)

end

################################################################################
########################### FlyByWireActuation #################################

@kwdef struct FlyByWireActuation <: C172.Actuation
    aileron_act::Actuator = Actuator(range = (-1.0, 1.0))
    elevator_act::Actuator = Actuator(range = (-1.0, 1.0))
    rudder_act::Actuator = Actuator(range = (-1.0, 1.0))
    throttle_act::Actuator = Actuator(range = (0.0, 1.0))
end

@kwdef mutable struct FlyByWireActuationU
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

@kwdef struct FlyByWireActuationY
    eng_start::Bool = false
    eng_stop::Bool = false
    mixture::Float64 = 0.5
    aileron_offset::Float64 = 0.0
    elevator_offset::Float64 = 0.0
    rudder_offset::Float64 = 0.0
    flaps::Float64 = 0.0
    brake_left::Float64 = 0.0
    brake_right::Float64 = 0.0
    throttle_act::ActuatorY = ActuatorY()
    aileron_act::ActuatorY = ActuatorY()
    elevator_act::ActuatorY = ActuatorY()
    rudder_act::ActuatorY = ActuatorY()
end

Systems.init(::SystemU, ::FlyByWireActuation) = FlyByWireActuationU()
Systems.init(::SystemY, ::FlyByWireActuation) = FlyByWireActuationY()

RigidBody.MassTrait(::System{FlyByWireActuation}) = HasNoMass()
RigidBody.AngMomTrait(::System{FlyByWireActuation}) = HasNoAngularMomentum()
RigidBody.WrenchTrait(::System{FlyByWireActuation}) = GetsNoExternalWrench()

function Systems.f_ode!(sys::System{FlyByWireActuation})

    @unpack throttle_act, aileron_act, elevator_act, rudder_act = sys

    #assign inputs to individual actuator subsystems
    throttle_act.u[] = Float64(sys.u.throttle)
    aileron_act.u[] = Float64(sys.u.aileron)
    elevator_act.u[] = Float64(sys.u.elevator)
    rudder_act.u[] = Float64(sys.u.rudder)

    #update individual actuator subsystems
    f_ode!(throttle_act)
    f_ode!(aileron_act)
    f_ode!(elevator_act)
    f_ode!(rudder_act)

    #assign direct feedthrough inputs
    @unpack eng_start, eng_stop, mixture,
            aileron_offset, elevator_offset, rudder_offset,
            flaps, brake_left, brake_right = sys.u

    sys.y = FlyByWireActuationY(; eng_start, eng_stop,
            mixture, aileron_offset, elevator_offset, rudder_offset,
            flaps, brake_left, brake_right,
            throttle_act = throttle_act.y, aileron_act = aileron_act.y,
            elevator_act = elevator_act.y, rudder_act = rudder_act.y)

end

function C172.assign!(aero::System{<:C172.Aero},
                ldg::System{<:C172.Ldg},
                pwp::System{<:Piston.Thruster},
                act::System{<:FlyByWireActuation})

    @unpack eng_start, eng_stop,
            throttle_act, aileron_act, elevator_act, rudder_act, mixture,
            aileron_offset, elevator_offset, rudder_offset,
            brake_left, brake_right, flaps = act.y

    pwp.engine.u.start = eng_start
    pwp.engine.u.stop = eng_stop
    pwp.engine.u.throttle = throttle_act.pos
    pwp.engine.u.mixture = mixture
    ldg.nose.steering.u[] = (rudder_offset + rudder_act.pos)
    ldg.left.braking.u[] = brake_left
    ldg.right.braking.u[] = brake_right
    aero.u.e = -(elevator_offset + elevator_act.pos)
    aero.u.a = (aileron_offset + aileron_act.pos)
    aero.u.r = -(rudder_offset + rudder_act.pos)
    aero.u.f = flaps

    return nothing
end


function GUI.draw(sys::System{FlyByWireActuation}, label::String = "Cessna 172R Fly-By-Wire Actuation")

    @unpack eng_start, eng_stop,
            aileron_offset, elevator_offset, rudder_offset,
            flaps, mixture, brake_left, brake_right,
            throttle_act, aileron_act, elevator_act, rudder_act = sys.y

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    CImGui.Dummy(10.0, 10.0)
    CImGui.Text("Engine Start: $(eng_start)")
    CImGui.Text("Engine Stop: $(eng_stop)")
    CImGui.Dummy(10.0, 10.0); CImGui.Separator()

     if CImGui.CollapsingHeader("Throttle")
        CImGui.Text("Throttle Command"); CImGui.SameLine(160); display_bar("", throttle_act.cmd, 0, 1)
        CImGui.Text("Throttle Position"); CImGui.SameLine(160); display_bar("", throttle_act.pos, 0, 1)
        @running_plot("Throttle Position", throttle_act.pos, 0, 1, 0.0, 120)
    end

    if CImGui.CollapsingHeader("Aileron")
        CImGui.Text("Aileron Command"); CImGui.SameLine(160); display_bar("", aileron_act.cmd, -1, 1)
        CImGui.Text("Aileron Position"); CImGui.SameLine(160); display_bar("", aileron_act.pos, -1, 1)
        @running_plot("Aileron Position", aileron_act.pos, -1, 1, 0.0, 120)
    end

    if CImGui.CollapsingHeader("Elevator")
        CImGui.Text("Elevator Command"); CImGui.SameLine(160); display_bar("", elevator_act.cmd, -1, 1)
        CImGui.Text("Elevator Position"); CImGui.SameLine(160); display_bar("", elevator_act.pos, -1, 1)
        @running_plot("Elevator Position", elevator_act.pos, -1, 1, 0.0, 120)
    end

    if CImGui.CollapsingHeader("Rudder")
        CImGui.Text("Rudder Command"); CImGui.SameLine(160); display_bar("", rudder_act.cmd, -1, 1)
        CImGui.Text("Rudder Position"); CImGui.SameLine(160); display_bar("", rudder_act.pos, -1, 1)
        @running_plot("Rudder Position", rudder_act.pos, -1, 1, 0.0, 120)
    end

    display_bar("Aileron Offset", aileron_offset, -1, 1)
    display_bar("Elevator Offset", elevator_offset, -1, 1)
    display_bar("Rudder Offset", rudder_offset, -1, 1)
    display_bar("Flaps", flaps, 0, 1)
    display_bar("Mixture", mixture, 0, 1)
    display_bar("Left Brake", brake_left, 0, 1)
    display_bar("Right Brake", brake_right, 0, 1)

    CImGui.PopItemWidth()

    CImGui.End()

end

function GUI.draw!(sys::System{FlyByWireActuation}, label::String = "Cessna 172R Mechanical Actuation")

    @unpack u, y = sys

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    CImGui.Dummy(10.0, 10.0)
    dynamic_button("Engine Start", 0.4); CImGui.SameLine()
    u.eng_start = CImGui.IsItemActive()
    dynamic_button("Engine Stop", 0.0)
    u.eng_stop = CImGui.IsItemActive()
    CImGui.Dummy(10.0, 10.0); CImGui.Separator()

    u.throttle = safe_slider("Throttle Command", u.throttle, "%.6f")
    @running_plot("Throttle Position", y.throttle_act.pos, 0, 1, 0.0, 120)
    u.aileron = safe_slider("Aileron Command", u.aileron, "%.6f")
    @running_plot("Aileron Position", y.aileron_act.pos, -1, 1, 0.0, 120)
    u.elevator = safe_slider("Elevator Command", u.elevator, "%.6f")
    @running_plot("Elevator Position", y.elevator_act.pos, -1, 1, 0.0, 120)
    u.rudder = safe_slider("Rudder Command", u.rudder, "%.6f")
    @running_plot("Rudder Position", y.rudder_act.pos, -1, 1, 0.0, 120)

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

# ################################## IODevices ###################################

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign!(sys::System{<:FlyByWireActuation},
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

    u.flaps += 0.3333 * was_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

    u.throttle += 0.1 * was_released(joystick, :button_Y)
    u.throttle -= 0.1 * was_released(joystick, :button_A)
end

function IODevices.assign!(sys::System{<:FlyByWireActuation},
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


################################################################################
################################# Template #####################################

#Cessna172 with default power plant and fly-by-wire actuation
const Template{K, V} = C172.Template{K, typeof(PowerPlant()), FlyByWireActuation, V} where {K, V}
Template(kinematics, avionics) = C172.Template(kinematics, PowerPlant(), FlyByWireActuation(), avionics)

# include(normpath("variants/base.jl")); @reexport using .C172FBWBase
# include(normpath("variants/cas.jl")); @reexport using .C172FBWCAS

end