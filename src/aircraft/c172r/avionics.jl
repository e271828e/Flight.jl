module C172RAvionics

using UnPack
using Printf
using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic

using Flight.FlightCore
using Flight.FlightPhysics
using Flight.FlightAircraft

using ..C172RAirframe

export DirectControls

################################################################################
############################### DirectControls #################################

struct DirectControls <: AbstractAvionics end

#we can simply reuse the airframe MechanicalActuation definitions, this is just
#a wrapper
const DirectControlsU = C172RAirframe.MechanicalActuationU
const DirectControlsY = C172RAirframe.MechanicalActuationY

Systems.init(::SystemU, ::DirectControls) = DirectControlsU()
Systems.init(::SystemY, ::DirectControls) = DirectControlsY()

########################### Update Methods #####################################

function Systems.f_ode!(avionics::System{DirectControls}, ::System{<:Airframe},
                ::KinematicData, ::AirData, ::System{<:AbstractTerrain})

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
            brake_left, brake_right = avionics.y

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
################################################################################

################################################################################
############################## StateMachine ####################################

################################################################################
############################### Autopilot ######################################

# struct AutopilotLogic <: Component end

Base.@kwdef struct RateCAS <: Component
    roll::PICompensator{1} = PICompensator{1}()
    pitch::PICompensator{1} = PICompensator{1}()
    yaw::PICompensator{1} = PICompensator{1}()
end

Base.@kwdef struct Autopilot <: Component
    # logic::AutopilotLogic = AutopilotLogic()
    rate::RateCAS = RateCAS()
end

# function Systems.f_ode!()

################################################################################
############################### XAvionics ######################################

struct StateMachine <: Component end

@enum FlightPhase begin
    on_ground = 0
    on_air = 1
end

@enum AutopilotState begin
    ap_disabled = 0
    ap_standby = 1
    ap_active = 2
end

Base.@kwdef mutable struct StateMachineU
    ap_enable::Bool = false
end

Base.@kwdef mutable struct StateMachineS
    flight_phase::FlightPhase = on_ground
    ap_state::AutopilotState = ap_disabled
end

Base.@kwdef struct StateMachineY
    flight_phase::FlightPhase = on_ground
    ap_state::AutopilotState = ap_disabled
end

#autopilot solo debe tener como inputs roll, pitch, yaw

Base.@kwdef struct XAvionics <: AbstractAvionics
    sm::StateMachine = StateMachine()
    ap::Autopilot = Autopilot()
end

#we can reuse DirectControlsU here, but noticing that for XAvionics aileron,
#elevator and rudder actually mean roll_input, pitch_input, yaw_input

#elevator↑ (stick forward) -> e↑ -> δe↑ -> trailing edge down -> Cm↓ -> pitch down
#aileron↑ (stick right) -> a↑ -> δa↑ -> left trailing edge down, right up -> Cl↓ -> roll right
#rudder↑ (right pedal forward) -> r↓ -> δr↓ -> rudder trailing edge right -> Cn↑ -> yaw right
#rudder↑ (right pedal forward) -> nose wheel steering right -> yaw right
#flaps↑ -> δf↑ -> flap trailing edge down -> CL↑
# Base.@kwdef mutable struct XAvionicsU
#     eng_start::Bool = false
#     eng_stop::Bool = false
#     throttle::Ranged{Float64, 0, 1} = 0.0
#     mixture::Ranged{Float64, 0, 1} = 0.5
#     roll_input::Ranged{Float64, -1, 1} = 0.0
#     pitch_input::Ranged{Float64, -1, 1} = 0.0
#     yaw_input::Ranged{Float64, -1, 1} = 0.0
#     aileron_trim::Ranged{Float64, -1, 1} = 0.0
#     elevator_trim::Ranged{Float64, -1, 1} = 0.0
#     rudder_trim::Ranged{Float64, -1, 1} = 0.0
#     flaps::Ranged{Float64, 0, 1} = 0.0
#     brake_left::Ranged{Float64, 0, 1} = 0.0
#     brake_right::Ranged{Float64, 0, 1} = 0.0
# end

# Base.@kwdef struct XAvionicsY
#     throttle::Float64 = 0.0
#     aileron_trim::Float64 = 0.0
#     aileron::Float64 = 0.0
#     elevator_trim::Float64 = 0.0
#     elevator::Float64 = 0.0
#     rudder_trim::Float64 = 0.0
#     rudder::Float64 = 0.0
#     brake_left::Float64 = 0.0
#     brake_right::Float64 = 0.0
#     flaps::Float64 = 0.0
#     mixture::Float64 = 0.5
#     eng_start::Bool = false
#     eng_stop::Bool = false
# end

# Systems.init(::SystemU, ::XAvionics) = FeedthroughActuationU()
# Systems.init(::SystemY, ::XAvionics) = (ap = init_y(ap), controls = init_y(controls))


########################### Update Methods #####################################

#fallback method?
# function Systems.f_ode!(avionics::System{XAvionics}, airframe::System{<:Airframe},
#                 kin::KinematicData, air::AirData, trn::System{<:AbstractTerrain})

#     @unpack sm, ap, act

#     #pregunta: en qué instantes deben mapearse los outputs del ap a controls? se
#     #hace en map_controls!, que va despues de esta llamada
#     f_ode!(ap)
#     f_ode!(act, airframe, kin, air, trn)
#     update_y!(avionics)
# end

# #no digital components or state machines in FeedthroughActuation
# @inline Systems.f_step!(::System{FeedthroughActuation}, ::System{<:Airframe}, ::KinematicSystem) = false
# @inline Systems.f_disc!(::System{FeedthroughActuation}, ::System{<:Airframe}, ::KinematicSystem, Δt) = false


# function Aircraft.map_controls!(airframe::System{<:Airframe}, avionics::System{XAvionics})

#     @unpack sm, ap, act = avionics

#     @unpack

#     # @unpack throttle, aileron_trim, aileron, elevator_trim, elevator,
#     #         rudder_trim, rudder, brake_left, brake_right, flaps, mixture,
#     #         eng_start, eng_stop = act.y

#     @unpack aero, pwp, ldg = airframe

#     pwp.u.engine.start = eng_start
#     pwp.u.engine.stop = eng_stop
#     pwp.u.engine.thr = throttle
#     pwp.u.engine.mix = mixture
#     ldg.u.nose.steering[] = (rudder_trim + rudder) #rudder↑ (right pedal forward) -> nose wheel steering right
#     ldg.u.left.braking[] = brake_left
#     ldg.u.right.braking[] = brake_right
#     aero.u.e = (elevator_trim + elevator) #elevator↑ (stick forward) -> e↑ -> pitch down
#     aero.u.a = (aileron_trim + aileron) #aileron↑ (stick right) -> a↑ -> roll right
#     aero.u.r = -(rudder_trim + rudder) #rudder↑ (right pedal forward) -> r↓ -> yaw right
#     aero.u.f = flaps #flaps↑ -> δf↑

#     return nothing
# end

#opcion 1:
#definimos un XAvionicsU identico a RevControlsU pero sustituyendo elevator por
#pitch_input, aileron por roll_input y rudder por yaw_input

#internamente hacemos que el origen de elevator de RevAvionicsU se determine de
#una de dos maneras, wow o no wow, y en funcion de ello le asignamos pitch_input
#o pitch_cas_output

#opcion 2:
#prescindimos de FeedthroughActuation, y consideramos solo XAvionics / AutopilotU
#mas adelante podemos meter actuators



end #module