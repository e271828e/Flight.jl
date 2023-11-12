module C172R

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, Reexport
using ControlSystems, RobustAndOptimalControl
using NLopt
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI
using Flight.FlightCore.IODevices
using Flight.FlightCore.Joysticks
using Flight.FlightCore.Utils: Ranged

using Flight.FlightPhysics.Attitude
using Flight.FlightPhysics.Geodesy
using Flight.FlightPhysics.Kinematics
using Flight.FlightPhysics.RigidBody
using Flight.FlightPhysics.Environment

using Flight.FlightComponents.Propellers
using Flight.FlightComponents.Piston
using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.Control
using Flight.FlightComponents.World

using ..C172


################################################################################
################################ Powerplant ####################################

function PowerPlant()

    prop_data = Propellers.Lookup(Propellers.Blade(), 2)

    propeller = Propeller(prop_data;
        sense = Propellers.CW, d = 2.0, J_xx = 0.3,
        t_bp = FrameTransform(r = [2.055, 0, 0.833]))

    Piston.Thruster(; propeller)

end

################################################################################
############################# Actuation ##################################

#pitch up -> Cm↑ -> trailing edge up -> δe↓ -> aero.e↓
### in order for a positive elevator actuation input (stick back) to
#produce a positive pitching moment, we need an inversion from elevator
#actuation input to actual elevator deflection (aero.e = -act.elevator). then we
#get: act.elevator↑ -> aero.e↓ -> pitch up

#roll right -> Cl↑ -> left trailing edge down, right up -> δa↑ -> aero.a↑ ->
#act.aileron↑ ### no act-aero inversion

#yaw right -> Cn↑ -> rudder trailing edge right -> δr↓ -> aero.r↓
### in order for a positive rudder actuation input (right pedal forward) to
#produce a positive yawing moment, we need an inversion from rudder actuation
#input to actual rudder deflection (aero.r = -act.rudder). then we get:
#act.rudder↑ -> aero.r↓ -> yaw right

#yaw right -> nose wheel steering right -> act.rudder↑ (right pedal forward) ### no
#act-nws inversion

#more lift -> CL↑ -> flap trailing edge down -> δf↑ -> aero.f↑ -> act.flaps↑ ### no
#act-aero inversion

#Reversible actuation system

struct Actuation <: C172.Actuation end

@kwdef mutable struct ActuationU
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

@kwdef struct ActuationY
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

Systems.init(::SystemU, ::Actuation) = ActuationU()
Systems.init(::SystemY, ::Actuation) = ActuationY()

function Systems.f_ode!(act::System{Actuation})

    @unpack eng_start, eng_stop,
            throttle, mixture, aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset,
            flaps, brake_left, brake_right= act.u

    act.y = ActuationY(; eng_start, eng_stop,
            throttle, mixture, aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset,
            flaps, brake_left, brake_right)

end

function C172.assign!(aero::System{<:C172.Aero},
                ldg::System{<:C172.Ldg},
                pwp::System{<:Piston.Thruster},
                act::System{<:Actuation})

    @unpack eng_start, eng_stop,
            throttle, mixture, aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset,
            brake_left, brake_right, flaps = act.y

    pwp.engine.u.start = eng_start
    pwp.engine.u.stop = eng_stop
    pwp.engine.u.throttle = throttle
    pwp.engine.u.mixture = mixture
    ldg.nose.steering.u[] = (rudder_offset + rudder)
    ldg.left.braking.u[] = brake_left
    ldg.right.braking.u[] = brake_right
    aero.u.e = -(elevator_offset + elevator)
    aero.u.a = (aileron_offset + aileron)
    aero.u.r = -(rudder_offset + rudder)
    aero.u.f = flaps

    return nothing
end


function GUI.draw(sys::System{Actuation}, label::String = "Cessna 172R Actuation")

    y = sys.y

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    CImGui.Text("Engine Start: $(y.eng_start)")
    CImGui.Text("Engine Stop: $(y.eng_stop)")

    @running_plot("Throttle", y.throttle, 0, 1, 0.0, 120)
    display_bar("Throttle", y.throttle, 0, 1)
    @running_plot("Aileron", y.aileron, -1, 1, 0.0, 120)
    display_bar("Aileron", y.aileron, -1, 1)
    @running_plot("Elevator", y.elevator, -1, 1, 0.0, 120)
    display_bar("Elevator", y.elevator, -1, 1)
    @running_plot("Rudder", y.rudder, -1, 1, 0.0, 120)
    display_bar("Rudder", y.rudder, -1, 1)

    display_bar("Aileron Offset", y.aileron_offset, -1, 1)
    display_bar("Elevator Offset", y.elevator_offset, -1, 1)
    display_bar("Rudder Offset", y.rudder_offset, -1, 1)
    display_bar("Flaps", y.flaps, 0, 1)
    display_bar("Mixture", y.mixture, 0, 1)
    display_bar("Left Brake", y.brake_left, 0, 1)
    display_bar("Right Brake", y.brake_right, 0, 1)

    CImGui.PopItemWidth()

    CImGui.End()

end

function GUI.draw!(sys::System{Actuation}, label::String = "Cessna 172R Actuation")

    u = sys.u

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    dynamic_button("Engine Start", 0.4); CImGui.SameLine()
    u.eng_start = CImGui.IsItemActive()
    dynamic_button("Engine Stop", 0.0)
    u.eng_stop = CImGui.IsItemActive()

    u.throttle = safe_slider("Throttle", u.throttle, "%.6f")
    @running_plot("Throttle", u.throttle, 0, 1, 0.0, 120)
    u.aileron = safe_slider("Aileron", u.aileron, "%.6f")
    @running_plot("Aileron", u.aileron, -1, 1, 0.0, 120)
    u.elevator = safe_slider("Elevator", u.elevator, "%.6f")
    @running_plot("Elevator", u.elevator, -1, 1, 0.0, 120)
    u.rudder = safe_slider("Rudder", u.rudder, "%.6f")
    @running_plot("Rudder", u.rudder, -1, 1, 0.0, 120)

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

################################## IODevices ###################################

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign!(sys::System{<:Actuation},
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

function IODevices.assign!(sys::System{<:Actuation},
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

const Airframe = C172.Airframe{typeof(PowerPlant()), Actuation}
const Physics{K} = Aircraft.Physics{K, Airframe} where {K <: AbstractKinematicDescriptor}
const Template{K, A} = Aircraft.Template{Physics{K}, A} where {K <: AbstractKinematicDescriptor, A <: AbstractAvionics}

Physics(kinematics = LTF()) = Aircraft.Physics(kinematics, C172.Airframe(PowerPlant(), Actuation()))
Template(kinematics = LTF(), avionics = NoAvionics()) = Aircraft.Template(Physics(kinematics), avionics)


############################### Trimming #######################################
################################################################################

#assigns trim state and parameters to aircraft physics, then updates aircraft physics
function C172.assign!(physics::System{<:C172R.Physics},
                env::System{<:AbstractEnvironment},
                trim_params::TrimParameters,
                trim_state::TrimState)

    @unpack EAS, β_a, x_fuel, flaps, mixture, payload = trim_params
    @unpack n_eng, α_a, throttle, aileron, elevator, rudder = trim_state
    @unpack act, pwp, aero, fuel, ldg, pld = physics.airframe

    init_kinematics!(physics, Kinematics.Initializer(trim_state, trim_params, env))

    #for trimming, control surface inputs are set to zero, and we work only with
    #their offsets
    act.u.throttle = throttle
    act.u.elevator = 0
    act.u.aileron = 0
    act.u.rudder = 0
    act.u.aileron_offset = aileron
    act.u.elevator_offset = elevator
    act.u.rudder_offset = rudder
    act.u.flaps = flaps
    act.u.mixture = mixture

    #assign payload
    @unpack m_pilot, m_copilot, m_lpass, m_rpass, m_baggage = payload
    @pack! pld.u = m_pilot, m_copilot, m_lpass, m_rpass, m_baggage

    #engine must be running
    pwp.engine.s.state = Piston.eng_running

    #set engine speed state
    ω_eng = n_eng * pwp.engine.constants.ω_rated
    pwp.x.engine.ω = ω_eng

    #engine idle compensator: as long as the engine remains at normal
    #operational speeds, well above its nominal idle speed, the idle controller
    #compensator's output will be saturated at its lower bound by proportional
    #error. its integrator will be disabled, its state will not change nor have
    #any effect on the engine. we can simply set it to zero
    pwp.x.engine.idle .= 0.0

    #engine friction compensator: with the engine running at normal operational
    #speeds, the engine's friction constraint compensator will be saturated, so
    #its integrator will be disabled and its state will not change. furthermore,
    #with the engine running friction is ignored. we can simply set it to zero.
    pwp.x.engine.frc .= 0.0

    #fuel content
    fuel.x .= Float64(x_fuel)

    aero.x.α_filt = α_a #ensures zero state derivative
    aero.x.β_filt = β_a #ensures zero state derivative

    f_ode!(physics, env)

    #check essential assumptions about airframe systems states & derivatives
    @assert !any(SVector{3}(leg.strut.wow for leg in ldg.y))
    @assert pwp.x.engine.ω > pwp.engine.constants.ω_idle
    @assert pwp.x.engine.idle[1] .== 0
    @assert pwp.x.engine.frc[1] .== 0
    @assert abs(aero.ẋ.α_filt) < 1e-10
    @assert abs(aero.ẋ.β_filt) < 1e-10

end




include(normpath("variants/base.jl")); @reexport using .C172RBase
# include(normpath("variants/cas.jl")); @reexport using .C172RCAS
# include(normpath("variants/direct.jl")); @reexport using .C172RDirect

end