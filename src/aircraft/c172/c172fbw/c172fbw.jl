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
################################ EMActuator ####################################

#second order system, defined in terms of damping ratio ζ and natural frequency
#ω_n. both command and output must be normalized between -1 and 1. see Ogata,
#p164
#input will be a Ranged{Float64, -1.0, 1.0}. with an underdamped system, due to
#overshoot, the position state can still transiently exceed [-1, 1]. therefore,
#we must bound the output. see PIDDiscrete
#x = ComponentVector(v = 0.0, p = 0.0)

################################################################################
############################### EMActuation ####################################

struct EMActuation <: C172.Actuation end

# @kwdef mutable struct MechanicalActuationU
#     eng_start::Bool = false
#     eng_stop::Bool = false
#     throttle::Ranged{Float64, 0., 1.} = 0.0
#     mixture::Ranged{Float64, 0., 1.} = 0.5
#     aileron::Ranged{Float64, -1., 1.} = 0.0
#     elevator::Ranged{Float64, -1., 1.} = 0.0
#     rudder::Ranged{Float64, -1., 1.} = 0.0
#     aileron_offset::Ranged{Float64, -1., 1.} = 0.0
#     elevator_offset::Ranged{Float64, -1., 1.} = 0.0
#     rudder_offset::Ranged{Float64, -1., 1.} = 0.0
#     flaps::Ranged{Float64, 0., 1.} = 0.0
#     brake_left::Ranged{Float64, 0., 1.} = 0.0
#     brake_right::Ranged{Float64, 0., 1.} = 0.0
# end

# @kwdef struct MechanicalActuationY
#     eng_start::Bool = false
#     eng_stop::Bool = false
#     throttle::Float64 = 0.0
#     mixture::Float64 = 0.5
#     aileron::Float64 = 0.0
#     elevator::Float64 = 0.0
#     rudder::Float64 = 0.0
#     aileron_offset::Float64 = 0.0
#     elevator_offset::Float64 = 0.0
#     rudder_offset::Float64 = 0.0
#     flaps::Float64 = 0.0
#     brake_left::Float64 = 0.0
#     brake_right::Float64 = 0.0
# end

# Systems.init(::SystemU, ::MechanicalActuation) = MechanicalActuationU()
# Systems.init(::SystemY, ::MechanicalActuation) = MechanicalActuationY()

# RigidBody.MassTrait(::System{MechanicalActuation}) = HasNoMass()
# RigidBody.AngMomTrait(::System{MechanicalActuation}) = HasNoAngularMomentum()
# RigidBody.WrenchTrait(::System{MechanicalActuation}) = GetsNoExternalWrench()

# function Systems.f_ode!(act::System{MechanicalActuation})

#     @unpack eng_start, eng_stop,
#             throttle, mixture, aileron, elevator, rudder,
#             aileron_offset, elevator_offset, rudder_offset,
#             flaps, brake_left, brake_right= act.u

#     act.y = MechanicalActuationY(; eng_start, eng_stop,
#             throttle, mixture, aileron, elevator, rudder,
#             aileron_offset, elevator_offset, rudder_offset,
#             flaps, brake_left, brake_right)

# end

# function C172.assign!(aero::System{<:C172.Aero},
#                 ldg::System{<:C172.Ldg},
#                 pwp::System{<:Piston.Thruster},
#                 act::System{<:MechanicalActuation})

#     @unpack eng_start, eng_stop,
#             throttle, mixture, aileron, elevator, rudder,
#             aileron_offset, elevator_offset, rudder_offset,
#             brake_left, brake_right, flaps = act.y

#     pwp.engine.u.start = eng_start
#     pwp.engine.u.stop = eng_stop
#     pwp.engine.u.throttle = throttle
#     pwp.engine.u.mixture = mixture
#     ldg.nose.steering.u[] = (rudder_offset + rudder)
#     ldg.left.braking.u[] = brake_left
#     ldg.right.braking.u[] = brake_right
#     aero.u.e = -(elevator_offset + elevator)
#     aero.u.a = (aileron_offset + aileron)
#     aero.u.r = -(rudder_offset + rudder)
#     aero.u.f = flaps

#     return nothing
# end


# function GUI.draw(sys::System{MechanicalActuation}, label::String = "Cessna 172R Mechanical Actuation")

#     y = sys.y

#     CImGui.Begin(label)

#     CImGui.PushItemWidth(-60)

#     CImGui.Text("Engine Start: $(y.eng_start)")
#     CImGui.Text("Engine Stop: $(y.eng_stop)")

#     @running_plot("Throttle", y.throttle, 0, 1, 0.0, 120)
#     display_bar("Throttle", y.throttle, 0, 1)
#     @running_plot("Aileron", y.aileron, -1, 1, 0.0, 120)
#     display_bar("Aileron", y.aileron, -1, 1)
#     @running_plot("Elevator", y.elevator, -1, 1, 0.0, 120)
#     display_bar("Elevator", y.elevator, -1, 1)
#     @running_plot("Rudder", y.rudder, -1, 1, 0.0, 120)
#     display_bar("Rudder", y.rudder, -1, 1)

#     display_bar("Aileron Offset", y.aileron_offset, -1, 1)
#     display_bar("Elevator Offset", y.elevator_offset, -1, 1)
#     display_bar("Rudder Offset", y.rudder_offset, -1, 1)
#     display_bar("Flaps", y.flaps, 0, 1)
#     display_bar("Mixture", y.mixture, 0, 1)
#     display_bar("Left Brake", y.brake_left, 0, 1)
#     display_bar("Right Brake", y.brake_right, 0, 1)

#     CImGui.PopItemWidth()

#     CImGui.End()

# end

# function GUI.draw!(sys::System{MechanicalActuation}, label::String = "Cessna 172R Mechanical Actuation")

#     u = sys.u

#     CImGui.Begin(label)

#     CImGui.PushItemWidth(-60)

#     dynamic_button("Engine Start", 0.4); CImGui.SameLine()
#     u.eng_start = CImGui.IsItemActive()
#     dynamic_button("Engine Stop", 0.0)
#     u.eng_stop = CImGui.IsItemActive()

#     u.throttle = safe_slider("Throttle", u.throttle, "%.6f")
#     @running_plot("Throttle", u.throttle, 0, 1, 0.0, 120)
#     u.aileron = safe_slider("Aileron", u.aileron, "%.6f")
#     @running_plot("Aileron", u.aileron, -1, 1, 0.0, 120)
#     u.elevator = safe_slider("Elevator", u.elevator, "%.6f")
#     @running_plot("Elevator", u.elevator, -1, 1, 0.0, 120)
#     u.rudder = safe_slider("Rudder", u.rudder, "%.6f")
#     @running_plot("Rudder", u.rudder, -1, 1, 0.0, 120)

#     u.aileron_offset = safe_input("Aileron Offset", u.aileron_offset, 0.001, 0.1, "%.6f")
#     u.elevator_offset = safe_input("Elevator Offset", u.elevator_offset, 0.001, 0.1, "%.6f")
#     u.rudder_offset = safe_input("Rudder Offset", u.rudder_offset, 0.001, 0.1, "%.6f")
#     u.flaps = safe_slider("Flaps", u.flaps, "%.6f")
#     u.mixture = safe_slider("Mixture", u.mixture, "%.6f")
#     u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
#     u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")

#     CImGui.PopItemWidth()

#     CImGui.End()

# end

# ################################## IODevices ###################################

# elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
# aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
# rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
# brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

# function IODevices.assign!(sys::System{<:MechanicalActuation},
#                            joystick::XBoxController,
#                            ::DefaultMapping)

#     u = sys.u

#     u.aileron = get_axis_value(joystick, :right_analog_x) |> aileron_curve
#     u.elevator = get_axis_value(joystick, :right_analog_y) |> elevator_curve
#     u.rudder = get_axis_value(joystick, :left_analog_x) |> rudder_curve
#     u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
#     u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

#     u.aileron_offset -= 0.01 * was_released(joystick, :dpad_left)
#     u.aileron_offset += 0.01 * was_released(joystick, :dpad_right)
#     u.elevator_offset += 0.01 * was_released(joystick, :dpad_down)
#     u.elevator_offset -= 0.01 * was_released(joystick, :dpad_up)

#     u.flaps += 0.3333 * was_released(joystick, :right_bumper)
#     u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

#     u.throttle += 0.1 * was_released(joystick, :button_Y)
#     u.throttle -= 0.1 * was_released(joystick, :button_A)
# end

# function IODevices.assign!(sys::System{<:MechanicalActuation},
#                            joystick::T16000M,
#                            ::DefaultMapping)

#     u = sys.u

#     u.throttle = get_axis_value(joystick, :throttle)
#     u.aileron = get_axis_value(joystick, :stick_x) |> aileron_curve
#     u.elevator = get_axis_value(joystick, :stick_y) |> elevator_curve
#     u.rudder = get_axis_value(joystick, :stick_z) |> rudder_curve

#     u.brake_left = is_pressed(joystick, :button_1)
#     u.brake_right = is_pressed(joystick, :button_1)

#     u.aileron_offset -= 2e-4 * is_pressed(joystick, :hat_left)
#     u.aileron_offset += 2e-4 * is_pressed(joystick, :hat_right)
#     u.elevator_offset += 2e-4 * is_pressed(joystick, :hat_down)
#     u.elevator_offset -= 2e-4 * is_pressed(joystick, :hat_up)

#     u.flaps += 0.3333 * was_released(joystick, :button_3)
#     u.flaps -= 0.3333 * was_released(joystick, :button_2)

# end


################################################################################
################################# Template #####################################

#Cessna172 with default power plant and mechanical actuation
const Template{K, V} = C172.Template{K, typeof(PowerPlant()), EMActuation, V} where {K, V}
Template(kinematics, avionics) = C172.Template(kinematics, PowerPlant(), EMActuation(), avionics)

# include(normpath("variants/base.jl")); @reexport using .C172FBWBase
# include(normpath("variants/cas.jl")); @reexport using .C172FBWCAS

end