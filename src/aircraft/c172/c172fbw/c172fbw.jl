module C172FBW

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

    #cache propeller lookup data to speed up aircraft instantiation. WARNING: if
    #the propeller definition or the lookup data generation methods in the
    #Propellers module are modified, the cache file must be regenerated
    # cache_file = joinpath(@__DIR__, "prop.h5")
    # if !isfile(cache_file)
    #     prop_data = Propellers.Lookup(Propellers.Blade(), 2)
    #     Propellers.save_lookup(prop_data, cache_file)
    # end
    # prop_data = Propellers.load_lookup(cache_file)

    #always generate the lookup data from scratch
    prop_data = Propellers.Lookup(Propellers.Blade(), 2)

    propeller = Propeller(prop_data;
        sense = Propellers.CW, d = 2.0, J_xx = 0.3,
        t_bp = FrameTransform(r = [2.055, 0, 0.833]))

    Piston.Thruster(; propeller)

end

################################################################################
################################## Actuator ####################################

@kwdef struct Actuator <: SystemDefinition #second order linear actuator model
    ω_n::Float64 = 5*2π #natural frequency (default: 10 Hz)
    ζ::Float64 = 0.6 #damping ratio (default: underdamped with minimal resonance)
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

    @unpack ẋ, x, u, constants = sys
    @unpack ω_n, ζ, range = constants

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
#################################### Actuation #################################

#Fly-by-wire actuation system. Throttle and aerodynamic surfaces are controlled
#via actuators, the rest of the actuation system is direct feedthrough

@kwdef struct Actuation <: C172.Actuation
    throttle_act::Actuator = Actuator(range = (0.0, 1.0))
    aileron_act::Actuator = Actuator(range = (-1.0, 1.0))
    elevator_act::Actuator = Actuator(range = (-1.0, 1.0))
    rudder_act::Actuator = Actuator(range = (-1.0, 1.0))
end

@kwdef mutable struct ActuationU
    eng_start::Bool = false
    eng_stop::Bool = false
    throttle_cmd::Ranged{Float64, 0., 1.} = 0.0
    mixture::Ranged{Float64, 0., 1.} = 0.5
    aileron_cmd::Ranged{Float64, -1., 1.} = 0.0
    elevator_cmd::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd::Ranged{Float64, -1., 1.} = 0.0
    aileron_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    elevator_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end

@kwdef struct ActuationY
    eng_start::Bool = false
    eng_stop::Bool = false
    throttle_cmd::Float64 = 0.0
    mixture::Float64 = 0.5
    aileron_cmd::Float64 = 0.0
    elevator_cmd::Float64 = 0.0
    rudder_cmd::Float64 = 0.0
    aileron_cmd_offset::Float64 = 0.0
    elevator_cmd_offset::Float64 = 0.0
    rudder_cmd_offset::Float64 = 0.0
    flaps::Float64 = 0.0
    brake_left::Float64 = 0.0
    brake_right::Float64 = 0.0
    throttle_act::ActuatorY = ActuatorY()
    aileron_act::ActuatorY = ActuatorY()
    elevator_act::ActuatorY = ActuatorY()
    rudder_act::ActuatorY = ActuatorY()
end

Systems.init(::SystemU, ::Actuation) = ActuationU()
Systems.init(::SystemY, ::Actuation) = ActuationY()

RigidBody.MassTrait(::System{Actuation}) = HasNoMass()
RigidBody.AngMomTrait(::System{Actuation}) = HasNoAngularMomentum()
RigidBody.WrenchTrait(::System{Actuation}) = GetsNoExternalWrench()

function Systems.f_ode!(sys::System{Actuation})

    @unpack throttle_act, aileron_act, elevator_act, rudder_act = sys

    @unpack eng_start, eng_stop, throttle_cmd, mixture,
            aileron_cmd, elevator_cmd, rudder_cmd,
            aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
            flaps, brake_left, brake_right = sys.u

    #assign inputs to actuator subsystems
    throttle_act.u[] = Float64(throttle_cmd)
    aileron_act.u[] = Float64(aileron_cmd + aileron_cmd_offset)
    elevator_act.u[] = Float64(elevator_cmd + elevator_cmd_offset)
    rudder_act.u[] = Float64(rudder_cmd + rudder_cmd_offset)

    #update actuator subsystems
    f_ode!(throttle_act)
    f_ode!(aileron_act)
    f_ode!(elevator_act)
    f_ode!(rudder_act)

    sys.y = ActuationY(;
            eng_start, eng_stop, throttle_cmd, mixture,
            aileron_cmd, elevator_cmd, rudder_cmd,
            aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
            flaps, brake_left, brake_right,
            throttle_act = throttle_act.y, aileron_act = aileron_act.y,
            elevator_act = elevator_act.y, rudder_act = rudder_act.y)

end

function C172.assign!(aero::System{<:C172.Aero},
                ldg::System{<:C172.Ldg},
                pwp::System{<:Piston.Thruster},
                act::System{<:Actuation})

    @unpack eng_start, eng_stop, mixture, flaps, brake_left, brake_right,
            throttle_act, aileron_act, elevator_act, rudder_act = act.y

    pwp.engine.u.start = eng_start
    pwp.engine.u.stop = eng_stop
    pwp.engine.u.throttle = throttle_act.pos
    pwp.engine.u.mixture = mixture
    ldg.nose.steering.u[] = rudder_act.pos
    ldg.left.braking.u[] = brake_left
    ldg.right.braking.u[] = brake_right
    aero.u.e = -elevator_act.pos
    aero.u.a = aileron_act.pos
    aero.u.r = -rudder_act.pos
    aero.u.f = flaps

    return nothing
end


function GUI.draw(sys::System{Actuation}, label::String = "Cessna 172R Fly-By-Wire Actuation")

    @unpack eng_start, eng_stop, throttle_cmd, mixture,
            aileron_cmd, elevator_cmd, rudder_cmd,
            aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
            flaps, brake_left, brake_right,
            throttle_act, aileron_act, elevator_act, rudder_act = sys.y

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    CImGui.Dummy(10.0, 10.0)
    CImGui.Text("Engine Start: $(eng_start)")
    CImGui.Text("Engine Stop: $(eng_stop)")
    CImGui.Dummy(10.0, 10.0);

    CImGui.Separator()

     if CImGui.CollapsingHeader("Throttle")
        CImGui.Text("Throttle Actuator Command"); CImGui.SameLine(200); display_bar("", throttle_act.cmd, 0, 1)
        CImGui.Text("Throttle Actuator Position"); CImGui.SameLine(200); display_bar("", throttle_act.pos, 0, 1)
        @running_plot("Throttle Actuator Position", throttle_act.pos, 0, 1, 0.0, 120)
        CImGui.Dummy(10.0, 10.0);
    end

    if CImGui.CollapsingHeader("Aileron")
        CImGui.Text("Aileron Command"); CImGui.SameLine(200); display_bar("", aileron_cmd, -1, 1)
        CImGui.Text("Aileron Command Offset"); CImGui.SameLine(200); display_bar("", aileron_cmd_offset, -1, 1)
        CImGui.Text("Aileron Actuator Command"); CImGui.SameLine(200); display_bar("", aileron_act.cmd, -1, 1)
        CImGui.Text("Aileron Actuator Position"); CImGui.SameLine(200); display_bar("", aileron_act.pos, -1, 1)
        @running_plot("Aileron Actuator Position", aileron_act.pos, -1, 1, 0.0, 120)
        CImGui.Dummy(10.0, 10.0);
    end

    if CImGui.CollapsingHeader("Elevator")
        CImGui.Text("Elevator Command"); CImGui.SameLine(200); display_bar("", elevator_cmd, -1, 1)
        CImGui.Text("Elevator Command Offset"); CImGui.SameLine(200); display_bar("", elevator_cmd_offset, -1, 1)
        CImGui.Text("Elevator Actuator Command"); CImGui.SameLine(200); display_bar("", elevator_act.cmd, -1, 1)
        CImGui.Text("Elevator Actuator Position"); CImGui.SameLine(200); display_bar("", elevator_act.pos, -1, 1)
        @running_plot("Elevator Actuator Position", elevator_act.pos, -1, 1, 0.0, 120)
        CImGui.Dummy(10.0, 10.0)
    end

    if CImGui.CollapsingHeader("Rudder")
        CImGui.Text("Rudder Command"); CImGui.SameLine(200); display_bar("", rudder_cmd, -1, 1)
        CImGui.Text("Rudder Command Offset"); CImGui.SameLine(200); display_bar("", rudder_cmd_offset, -1, 1)
        CImGui.Text("Rudder Actuator Command"); CImGui.SameLine(200); display_bar("", rudder_act.cmd, -1, 1)
        CImGui.Text("Rudder Actuator Position"); CImGui.SameLine(200); display_bar("", rudder_act.pos, -1, 1)
        @running_plot("Rudder Position", rudder_act.pos, -1, 1, 0.0, 120)
        CImGui.Dummy(10.0, 10.0)
    end

    CImGui.Separator()

    CImGui.Dummy(10.0, 10.0)
    display_bar("Flaps", flaps, 0, 1)
    display_bar("Mixture", mixture, 0, 1)
    display_bar("Left Brake", brake_left, 0, 1)
    display_bar("Right Brake", brake_right, 0, 1)

    CImGui.PopItemWidth()

    CImGui.End()

end

function GUI.draw!(sys::System{Actuation}, label::String = "Cessna 172R Fly By Wire Actuation")

    @unpack u, y = sys

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    CImGui.Dummy(10.0, 10.0)
    dynamic_button("Engine Start", 0.4); CImGui.SameLine()
    u.eng_start = CImGui.IsItemActive()
    dynamic_button("Engine Stop", 0.0)
    u.eng_stop = CImGui.IsItemActive()
    CImGui.Dummy(10.0, 10.0);

    CImGui.Separator()
    u.throttle_cmd = safe_slider("Throttle Command", u.throttle_cmd, "%.6f")
    CImGui.Text("Throttle Actuator Command"); CImGui.SameLine(200); display_bar("", y.throttle_act.cmd, 0, 1)
    CImGui.Text("Throttle Actuator Position"); CImGui.SameLine(200); display_bar("", y.throttle_act.pos, 0, 1)
    @running_plot("Throttle Actuator Position", y.throttle_act.pos, 0, 1, 0.0, 120)
    CImGui.Dummy(10.0, 10.0)

    u.aileron_cmd = safe_slider("Aileron Command", u.aileron_cmd, "%.6f")
    u.aileron_cmd_offset = safe_input("Aileron Command Offset", u.aileron_cmd_offset, 0.001, 0.1, "%.6f")
    CImGui.Text("Aileron Actuator Command"); CImGui.SameLine(200); display_bar("", y.aileron_act.cmd, -1, 1)
    CImGui.Text("Aileron Actuator Position"); CImGui.SameLine(200); display_bar("", y.aileron_act.pos, -1, 1)
    @running_plot("Aileron Actuator Position", y.aileron_act.pos, -1, 1, 0.0, 120)
    CImGui.Dummy(10.0, 10.0)

    u.elevator_cmd = safe_slider("Elevator Command", u.elevator_cmd, "%.6f")
    u.elevator_cmd_offset = safe_input("Elevator Command Offset", u.elevator_cmd_offset, 0.001, 0.1, "%.6f")
    CImGui.Text("Elevator Actuator Command"); CImGui.SameLine(200); display_bar("", y.elevator_act.cmd, -1, 1)
    CImGui.Text("Elevator Actuator Position"); CImGui.SameLine(200); display_bar("", y.elevator_act.pos, -1, 1)
    @running_plot("Elevator Position", y.elevator_act.pos, -1, 1, 0.0, 120)
    CImGui.Dummy(10.0, 10.0)

    u.rudder_cmd = safe_slider("Rudder Command", u.rudder_cmd, "%.6f")
    u.rudder_cmd_offset = safe_input("Rudder Command Offset", u.rudder_cmd_offset, 0.001, 0.1, "%.6f")
    CImGui.Text("Rudder Actuator Command"); CImGui.SameLine(200); display_bar("", y.rudder_act.cmd, -1, 1)
    CImGui.Text("Rudder Actuator Position"); CImGui.SameLine(200); display_bar("", y.rudder_act.pos, -1, 1)
    @running_plot("Rudder Position", y.rudder_act.pos, -1, 1, 0.0, 120)
    CImGui.Separator()

    CImGui.Dummy(10.0, 10.0)
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

function IODevices.assign!(sys::System{<:Actuation},
                           joystick::XBoxController,
                           ::DefaultMapping)

    u = sys.u

    u.aileron_cmd = get_axis_value(joystick, :right_analog_x) |> aileron_curve
    u.elevator_cmd = get_axis_value(joystick, :right_analog_y) |> elevator_curve
    u.rudder_cmd = get_axis_value(joystick, :left_analog_x) |> rudder_curve
    u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

    u.aileron_cmd_offset -= 0.01 * was_released(joystick, :dpad_left)
    u.aileron_cmd_offset += 0.01 * was_released(joystick, :dpad_right)
    u.elevator_cmd_offset += 0.01 * was_released(joystick, :dpad_down)
    u.elevator_cmd_offset -= 0.01 * was_released(joystick, :dpad_up)

    u.flaps += 0.3333 * was_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

    u.throttle_cmd += 0.1 * was_released(joystick, :button_Y)
    u.throttle_cmd -= 0.1 * was_released(joystick, :button_A)
end

function IODevices.assign!(sys::System{<:Actuation},
                           joystick::T16000M,
                           ::DefaultMapping)

    u = sys.u

    u.throttle_cmd = get_axis_value(joystick, :throttle)
    u.aileron_cmd = get_axis_value(joystick, :stick_x) |> aileron_curve
    u.elevator_cmd = get_axis_value(joystick, :stick_y) |> elevator_curve
    u.rudder_cmd = get_axis_value(joystick, :stick_z) |> rudder_curve

    u.brake_left = is_pressed(joystick, :button_1)
    u.brake_right = is_pressed(joystick, :button_1)

    u.aileron_cmd_offset -= 2e-4 * is_pressed(joystick, :hat_left)
    u.aileron_cmd_offset += 2e-4 * is_pressed(joystick, :hat_right)
    u.elevator_cmd_offset += 2e-4 * is_pressed(joystick, :hat_down)
    u.elevator_cmd_offset -= 2e-4 * is_pressed(joystick, :hat_up)

    u.flaps += 0.3333 * was_released(joystick, :button_3)
    u.flaps -= 0.3333 * was_released(joystick, :button_2)

end


################################################################################
################################# Template #####################################

const Airframe = C172.Airframe{typeof(PowerPlant()), Actuation}
const Physics{K} = Aircraft.Physics{K, Airframe}
const Template{K, A} = Aircraft.Template{Physics{K}, A} where {K, A}

Physics(kinematics = LTF()) = Aircraft.Physics(kinematics, C172.Airframe(PowerPlant(), Actuation()))
Template(kinematics = LTF(), avionics = NoAvionics()) = Aircraft.Template(Physics(kinematics), avionics)


############################### Trimming #######################################
################################################################################

#first 2 are aircraft-agnostic
@kwdef struct TrimState <: FieldVector{7, Float64}
    α_a::Float64 = 0.1 #angle of attack, aerodynamic axes
    φ_nb::Float64 = 0.0 #bank angle
    n_eng::Float64 = 0.75 #normalized engine speed (ω/ω_rated)
    throttle::Float64 = 0.47
    aileron::Float64 = 0.014
    elevator::Float64 = -0.0015
    rudder::Float64 = 0.02 #rudder↑ -> aero.u.r↓ -> right yaw
end

@kwdef struct TrimParameters <: AbstractTrimParameters
    Ob::Geographic{NVector, Ellipsoidal} = Geographic(NVector(), HOrth(1000))
    ψ_nb::Float64 = 0.0 #geographic heading
    EAS::Float64 = 50.0 #equivalent airspeed
    γ_wOb_n::Float64 = 0.0 #wind-relative flight path angle
    ψ_lb_dot::Float64 = 0.0 #LTF-relative turn rate
    θ_lb_dot::Float64 = 0.0 #LTF-relative pitch rate
    β_a::Float64 = 0.0 #sideslip angle measured in the aerodynamic reference frame
    x_fuel::Ranged{Float64, 0., 1.} = 0.5 #normalized fuel load
    mixture::Ranged{Float64, 0., 1.} = 0.5 #engine mixture control
    flaps::Ranged{Float64, 0., 1.} = 0.0 #flap setting
    payload::C172.PayloadU = C172.PayloadU()
end


function Kinematics.Initializer(trim_state::TrimState,
                                trim_params::TrimParameters,
                                env::System{<:AbstractEnvironment})

    @unpack EAS, β_a, γ_wOb_n, ψ_nb, ψ_lb_dot, θ_lb_dot, Ob = trim_params
    @unpack α_a, φ_nb = trim_state

    atm_data = AtmosphericData(env.atm, Ob)
    TAS = Atmosphere.EAS2TAS(EAS; ρ = atm_data.ISA.ρ)
    v_wOb_a = Atmosphere.get_velocity_vector(TAS, α_a, β_a)
    v_wOb_b = C172.f_ba.q(v_wOb_a) #wind-relative aircraft velocity, body frame

    θ_nb = Aircraft.θ_constraint(; v_wOb_b, γ_wOb_n, φ_nb)
    e_nb = REuler(ψ_nb, θ_nb, φ_nb)
    q_nb = RQuat(e_nb)

    e_lb = e_nb #initialize LTF arbitrarily to NED
    ė_lb = SVector(ψ_lb_dot, θ_lb_dot, 0.0)
    ω_lb_b = Attitude.ω(e_lb, ė_lb)

    loc = NVector(Ob)
    h = HEllip(Ob)

    v_wOb_n = q_nb(v_wOb_b) #wind-relative aircraft velocity, NED frame
    v_ew_n = atm_data.wind.v_ew_n
    v_eOb_n = v_ew_n + v_wOb_n

    Kinematics.Initializer(; q_nb, loc, h, ω_lb_b, v_eOb_n, Δx = 0.0, Δy = 0.0)

end

#assigns trim state and parameters to the aircraft system, then updates it
function assign!(physics::System{<:C172FBW.Physics},
                env::System{<:AbstractEnvironment},
                trim_params::TrimParameters,
                trim_state::TrimState)

    @unpack EAS, β_a, x_fuel, flaps, mixture, payload = trim_params
    @unpack n_eng, α_a, throttle, aileron, elevator, rudder = trim_state
    @unpack act, pwp, aero, fuel, ldg, pld = physics.airframe

    init_kinematics!(physics, Kinematics.Initializer(trim_state, trim_params, env))

    #for trimming, control surface inputs are set to zero, and we work only with
    #their offsets
    act.u.throttle_cmd = throttle
    act.u.elevator_cmd = 0
    act.u.aileron_cmd = 0
    act.u.rudder_cmd = 0
    act.u.aileron_cmd_offset = aileron
    act.u.elevator_cmd_offset = elevator
    act.u.rudder_cmd_offset = rudder
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

    #actuator states: in steady state every actuator's velocity state must be
    #zero, and its position state must be equal to the actuator command. the
    #actuator command is in turn equal to the surface command plus its offset,
    #which we have set to zero
    act.x.throttle_act.v = 0.0
    act.x.throttle_act.p = throttle
    act.x.aileron_act.v = 0.0
    act.x.aileron_act.p = aileron
    act.x.elevator_act.v = 0.0
    act.x.elevator_act.p = elevator
    act.x.rudder_act.v = 0.0
    act.x.rudder_act.p = rudder

    aero.x.α_filt = α_a #ensures zero state derivative
    aero.x.β_filt = β_a #ensures zero state derivative
    fuel.x .= Float64(x_fuel)

    f_ode!(physics, env)

    #check essential assumptions about airframe systems states & derivatives
    @assert !any(SVector{3}(leg.strut.wow for leg in ldg.y))
    @assert pwp.x.engine.ω > pwp.engine.constants.ω_idle
    @assert pwp.x.engine.idle[1] .== 0
    @assert pwp.x.engine.frc[1] .== 0
    @assert abs(aero.ẋ.α_filt) < 1e-10
    @assert abs(aero.ẋ.β_filt) < 1e-10

    @assert all(SVector{8,Float64}(act.ẋ) .== 0)

end

function cost(physics::System{<:C172FBW.Physics})

    @unpack ẋ, y = physics

    v_nd_dot = SVector{3}(ẋ.kinematics.vel.v_eOb_b) / norm(y.kinematics.common.v_eOb_b)
    ω_dot = SVector{3}(ẋ.kinematics.vel.ω_eb_b) #ω should already of order 1
    n_eng_dot = ẋ.airframe.pwp.engine.ω / physics.airframe.pwp.engine.constants.ω_rated

    sum(v_nd_dot.^2) + sum(ω_dot.^2) + n_eng_dot^2

end

function get_f_target(physics::System{<:C172FBW.Physics},
                      trim_params::TrimParameters,
                      env::System{<:AbstractEnvironment})

    let physics = physics, env = env, trim_params = trim_params
        function (x::TrimState)
            assign!(physics, env, trim_params, x)
            return cost(physics)
        end
    end

end

function Aircraft.trim!(physics::System{<:C172FBW.Physics},
                        trim_params::TrimParameters = TrimParameters(),
                        env::System{<:AbstractEnvironment} = System(SimpleEnvironment()))

    trim_state = TrimState() #could initial condition as an optional input

    f_target = get_f_target(physics, trim_params, env)

    #wrapper with the interface expected by NLopt
    f_opt(x::Vector{Float64}, ::Vector{Float64}) = f_target(TrimState(x))

    n = length(trim_state)
    x0 = zeros(n); lower_bounds = similar(x0); upper_bounds = similar(x0); initial_step = similar(x0)

    x0[:] .= trim_state

    lower_bounds[:] .= TrimState(
        α_a = -π/12,
        φ_nb = -π/3,
        n_eng = 0.4,
        throttle = 0,
        aileron = -1,
        elevator = -1,
        rudder = -1)

    upper_bounds[:] .= TrimState(
        α_a = physics.airframe.aero.constants.α_stall[2], #critical AoA is 0.28 < 0.36
        φ_nb = π/3,
        n_eng = 1.1,
        throttle = 1,
        aileron = 1,
        elevator = 1,
        rudder = 1)

    initial_step[:] .= 0.05 #safe value for all optimization variables

    #any of these three algorithms works
    # opt = Opt(:LN_NELDERMEAD, length(x0))
    opt = Opt(:LN_BOBYQA, length(x0))
    # opt = Opt(:GN_CRS2_LM, length(x0))
    opt.min_objective = f_opt
    opt.maxeval = 100000
    opt.stopval = 1e-14
    opt.lower_bounds = lower_bounds
    opt.upper_bounds = upper_bounds
    opt.initial_step = initial_step

    # @btime optimize($opt, $x0)

    (minf, minx, exit_flag) = optimize(opt, x0)

    success = (exit_flag === :STOPVAL_REACHED)
    if !success
        println("Warning: Trimming optimization failed with exit_flag $exit_flag")
    end
    trim_state_opt = TrimState(minx)
    assign!(physics, env, trim_params, trim_state_opt)
    return (success = success, trim_state = trim_state_opt)

end


################################################################################
############################### Linearization ##################################

#flaps and mixture are trim parameters and thus omitted from the control vector
@kwdef struct ULinear <: FieldVector{4, Float64}
    throttle_cmd::Float64 = 0.0
    aileron_cmd::Float64 = 0.0
    elevator_cmd::Float64 = 0.0
    rudder_cmd::Float64 = 0.0
end

@kwdef struct XLinear <: FieldVector{24, Float64}
    ψ::Float64 = 0.0; θ::Float64 = 0.0; φ::Float64 = 0.0; #heading, inclination, bank (body/NED)
    ϕ::Float64 = 0.0; λ::Float64 = 0.0; h::Float64 = 0.0; #latitude, longitude, ellipsoidal altitude
    p::Float64 = 0.0; q::Float64 = 0.0; r::Float64 = 0.0; #angular rates (ω_eb_b)
    v_x::Float64 = 0.0; v_y::Float64 = 0.0; v_z::Float64 = 0.0; #Ob/ECEF velocity, body axes
    α_filt::Float64 = 0.0; β_filt::Float64 = 0.0; #filtered airflow angles
    ω_eng::Float64 = 0.0; fuel::Float64 = 0.0; #engine speed, fuel fraction
    thr_v::Float64 = 0.0; thr_p::Float64 = 0.0; #throttle actuator states
    ail_v::Float64 = 0.0; ail_p::Float64 = 0.0; #aileron actuator states
    ele_v::Float64 = 0.0; ele_p::Float64 = 0.0; #elevator actuator states
    rud_v::Float64 = 0.0; rud_p::Float64 = 0.0 #rudder actuator states
end

@kwdef struct YLinear <: FieldVector{24, Float64}
    ψ::Float64 = 0.0; θ::Float64 = 0.0; φ::Float64 = 0.0; #heading, inclination, bank (body/NED)
    ϕ::Float64 = 0.0; λ::Float64 = 0.0; h::Float64 = 0.0; #latitude, longitude, ellipsoidal altitude
    p::Float64 = 0.0; q::Float64 = 0.0; r::Float64 = 0.0; #angular rates (ω_eb_b)
    EAS::Float64 = 0.0; TAS::Float64 = 0.0; #airspeed
    α::Float64 = 0.0; β::Float64 = 0.0; #airspeed
    f_x::Float64 = 0.0; f_y::Float64 = 0.0; f_z::Float64 = 0.0; #specific force at G (f_iG_b)
    v_N::Float64 = 0.0; v_E::Float64 = 0.0; v_D::Float64 = 0.0; #Ob/ECEF velocity, NED axes
    χ::Float64 = 0.0; γ::Float64 = 0.0; c::Float64 = 0.0; #track and flight path angles, climb rate
    ω_eng::Float64 = 0.0; m_fuel::Float64 = 0.0 #engine speed, fuel mass
end


function ULinear(physics::System{<:C172FBW.Physics})

    @unpack throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd = physics.airframe.act.u
    ULinear(; throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd)

end

function XLinear(x_physics::ComponentVector)

    x_kinematics = x_physics.kinematics
    x_airframe = x_physics.airframe

    @unpack ψ_nb, θ_nb, φ_nb, ϕ, λ, h_e = x_kinematics.pos
    p, q, r = x_kinematics.vel.ω_eb_b
    v_x, v_y, v_z = x_kinematics.vel.v_eOb_b
    α_filt, β_filt = x_airframe.aero
    ω_eng = x_airframe.pwp.engine.ω
    fuel = x_airframe.fuel[1]
    thr_v = x_airframe.act.throttle_act.v
    thr_p = x_airframe.act.throttle_act.p
    ail_v = x_airframe.act.aileron_act.v
    ail_p = x_airframe.act.aileron_act.p
    ele_v = x_airframe.act.elevator_act.v
    ele_p = x_airframe.act.elevator_act.p
    rud_v = x_airframe.act.rudder_act.v
    rud_p = x_airframe.act.rudder_act.p

    ψ, θ, φ, h = ψ_nb, θ_nb, φ_nb, h_e

    XLinear(;  ψ, θ, φ, ϕ, λ, h, p, q, r, v_x, v_y, v_z,
                        α_filt, β_filt, ω_eng, fuel,
                        thr_v, thr_p, ail_v, ail_p, ele_v, ele_p, rud_v, rud_p)

end

ẊLinear(physics::System{<:C172FBW.Physics}) = XLinear(physics.ẋ)
XLinear(physics::System{<:C172FBW.Physics}) = XLinear(physics.x)

function YLinear(physics::System{<:C172FBW.Physics})

    @unpack e_nb, ϕ_λ, h_e, ω_eb_b, v_eOb_n = physics.y.kinematics
    @unpack ψ, θ, φ = e_nb
    @unpack ϕ, λ = ϕ_λ

    h = h_e
    p, q, r = ω_eb_b
    v_N, v_E, v_D = v_eOb_n
    χ = Attitude.azimuth(v_eOb_n)
    γ = Attitude.inclination(v_eOb_n)
    c = -v_D
    f_x, f_y, f_z = physics.y.rigidbody.f_G_b
    EAS = physics.y.air.EAS
    TAS = physics.y.air.TAS
    α = physics.y.air.α_b
    β = physics.y.air.β_b
    ω_eng = physics.y.airframe.pwp.engine.ω
    m_fuel = physics.y.airframe.fuel.m_avail

    YLinear(; ψ, θ, φ, ϕ, λ, h, p, q, r, EAS, TAS, α, β,
            f_x, f_y, f_z, v_N, v_E, v_D, χ, γ, c, ω_eng, m_fuel)

end

function assign!(physics::System{<:C172FBW.Physics}, u::ULinear)

    @unpack throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd = u
    @pack! physics.airframe.act.u = throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd

end

function assign!(physics::System{<:C172FBW.Physics}, x::XLinear)

    @unpack ψ, θ, φ, ϕ, λ, h, p, q, r, v_x, v_y, v_z, α_filt, β_filt, ω_eng,
            fuel, thr_v, thr_p, ail_v, ail_p, ele_v, ele_p, rud_v, rud_p = x

    x_kinematics = physics.x.kinematics
    x_airframe = physics.x.airframe

    ψ_nb, θ_nb, φ_nb, h_e = ψ, θ, φ, h

    @pack! x_kinematics.pos = ψ_nb, θ_nb, φ_nb, ϕ, λ, h_e
    x_kinematics.vel.ω_eb_b .= p, q, r
    x_kinematics.vel.v_eOb_b .= v_x, v_y, v_z
    x_airframe.aero .= α_filt, β_filt
    x_airframe.pwp.engine.ω = ω_eng
    x_airframe.fuel .= fuel
    x_airframe.act.throttle_act.v = thr_v
    x_airframe.act.throttle_act.p = thr_p
    x_airframe.act.aileron_act.v = ail_v
    x_airframe.act.aileron_act.p = ail_p
    x_airframe.act.elevator_act.v = ele_v
    x_airframe.act.elevator_act.p = ele_p
    x_airframe.act.rudder_act.v = rud_v
    x_airframe.act.rudder_act.p = rud_p

end


function Aircraft.linearize!(
            physics::System{<:C172FBW.Physics{NED}},
            trim_params::TrimParameters = TrimParameters(),
            env::System{<:AbstractEnvironment} = System(SimpleEnvironment()))

    (_, trim_state) = trim!(physics, trim_params, env)

    ẋ0 = ẊLinear(physics)
    x0 = XLinear(physics)
    u0 = ULinear(physics)
    y0 = YLinear(physics)

    #f_main will not be returned for use in another scope, so we don't need to
    #capture physics and env with a let block, because they are guaranteed not
    #be reassigned within the scope of linearize!
    function f_main(x, u)

        assign!(physics, XLinear(x))
        assign!(physics, ULinear(u))
        f_ode!(physics, env)

        return (ẋ = ẊLinear(physics), y = YLinear(physics))

    end

    (A, B, C, D) = ss_matrices(f_main, x0, u0)

    #restore the System to its trimmed condition
    assign!(physics, env, trim_params, trim_state)

    #now we need to rebuild vectors and matrices for the LinearStateSpace as
    #ComponentArrays, because we want matrix components to remain labelled,
    #which cannot be achieved with FieldVectors

    x_axis = Axis(propertynames(x0))
    u_axis = Axis(propertynames(u0))
    y_axis = Axis(propertynames(y0))

    ẋ0_cv = ComponentVector(ẋ0, x_axis)
    x0_cv = ComponentVector(x0, x_axis)
    u0_cv = ComponentVector(u0, u_axis)
    y0_cv = ComponentVector(y0, y_axis)

    A_cv = ComponentMatrix(A, x_axis, x_axis)
    B_cv = ComponentMatrix(B, x_axis, u_axis)
    C_cv = ComponentMatrix(C, y_axis, x_axis)
    D_cv = ComponentMatrix(D, y_axis, u_axis)

    return LinearStateSpace(ẋ0_cv, x0_cv, u0_cv, y0_cv, A_cv, B_cv, C_cv, D_cv)

end


function ss_matrices(f_main::Function, x0::XLinear, u0::ULinear)

    f_ẋ(x, u) = f_main(x, u).ẋ
    f_y(x, u) = f_main(x, u).y

    #none of these closures will be returned for use in another scope, so we
    #don't need to capture x0 and u0 with a let block, because they are
    #guaranteed not be reassigned within the scope of ss_matrices
    f_A!(ẋ, x) = (ẋ .= f_ẋ(x, u0))
    f_B!(ẋ, u) = (ẋ .= f_ẋ(x0, u))
    f_C!(y, x) = (y .= f_y(x, u0))
    f_D!(y, u) = (y .= f_y(x0, u))

    #preallocate mutable arrays
    A = XLinear() * XLinear()' |> Matrix
    B = XLinear() * ULinear()' |> Matrix
    C = YLinear() * XLinear()' |> Matrix
    D = YLinear() * ULinear()' |> Matrix

    jacobian!(A, f_A!, Vector(x0))
    jacobian!(B, f_B!, Vector(u0))
    jacobian!(C, f_C!, Vector(x0))
    jacobian!(D, f_D!, Vector(u0))

    return (A, B, C, D)

end

function Control.LinearStateSpace(
            physics::System{<:C172FBW.Physics{NED}},
            trim_params::TrimParameters = TrimParameters(),
            env::System{<:AbstractEnvironment} = System(SimpleEnvironment());
            model::Symbol = :full)

    lm = linearize!(physics, trim_params, env)

    if model === :full
        return lm

    elseif model === :lon
        x_labels = [:q, :θ, :v_x, :v_z, :α_filt, :ω_eng, :thr_v, :thr_p, :ele_v, :ele_p]
        u_labels = [:throttle_cmd, :elevator_cmd]
        y_labels = [:q, :θ, :α, :EAS, :TAS, :f_x, :f_z, :γ, :c, :ω_eng, :v_D]
        return submodel(lm; x = x_labels, u = u_labels, y = y_labels)

    elseif model === :lat
        u_labels = [:aileron_cmd, :rudder_cmd]
        x_labels = [:p, :r, :φ, :ψ, :v_x, :v_y, :β_filt, :ail_v, :ail_p, :rud_v, :rud_p]
        y_labels = [:p, :r, :φ, :ψ, :β, :f_y, :χ]
        return submodel(lm; x = x_labels, u = u_labels, y = y_labels)

    else
        error("Valid model keyword values: :full, :lon, :lat")

    end

end

function Control.LinearStateSpace(ac::System{<:C172FBW.Template{NED}}, args...; kwargs...)
    LinearStateSpace(ac.physics, args...; kwargs...)
end

function RobustAndOptimalControl.named_ss(
            physics::System{<:C172FBW.Physics{NED}}, args...; kwargs...)

    lss = LinearStateSpace(physics, args...; kwargs...)
    x_labels, u_labels, y_labels = map(collect ∘ propertynames, (lss.x0, lss.u0, lss.y0))
    return named_ss(ss(lss), x = x_labels, u = u_labels, y = y_labels)
end

function RobustAndOptimalControl.named_ss(ac::System{<:C172FBW.Template{NED}}, args...; kwargs...)
    named_ss(ac.physics, args...; kwargs...)
end

################################################################################
################################## Variants ####################################

include(normpath("variants/base.jl")); @reexport using .C172FBWBase
include(normpath("variants/cas/cas.jl")); @reexport using .C172FBWCAS

end