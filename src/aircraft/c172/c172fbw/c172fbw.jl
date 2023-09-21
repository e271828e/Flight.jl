module C172FBW

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, Reexport
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
using Flight.FlightComponents.World

using ..C172

################################################################################
################################ Powerplant ####################################

PowerPlant() = Piston.Thruster(propeller = Propeller(t_bp = FrameTransform(r = [2.055, 0, 0.833])))

################################################################################
################################## Actuator ####################################

@kwdef struct Actuator <: SystemDefinition #second order linear actuator model
    ω_n::Float64 = 10*2π #natural frequency (default: 10 Hz)
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

#Cessna172 with default power plant, fly-by-wire actuation and any avionics
const Template{K, V} = C172.Template{K, typeof(PowerPlant()), Actuation, V} where {K, V}
Template(kinematics, avionics) = C172.Template(kinematics, PowerPlant(), Actuation(), avionics)


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

@kwdef struct TrimParameters
    Ob::Geographic{NVector, Ellipsoidal} = Geographic(NVector(), HOrth(1000))
    ψ_nb::Float64 = 0.0 #geographic heading
    TAS::Float64 = 40.0 #true airspeed
    γ_wOb_n::Float64 = 0.0 #wind-relative flight path angle
    ψ_lb_dot::Float64 = 0.0 #LTF-relative turn rate
    θ_lb_dot::Float64 = 0.0 #LTF-relative pitch rate
    β_a::Float64 = 0.0 #sideslip angle measured in the aerodynamic reference frame
    fuel::Float64 = 0.5 #fuel load, 0 to 1
    mixture::Float64 = 0.5 #engine mixture control, 0 to 1
    flaps::Float64 = 0.0 #flap setting, 0 to 1
end


function Kinematics.Initializer(trim_state::TrimState,
                                trim_params::TrimParameters,
                                env::System{<:AbstractEnvironment})

    @unpack TAS, β_a, γ_wOb_n, ψ_nb, ψ_lb_dot, θ_lb_dot, Ob = trim_params
    @unpack α_a, φ_nb = trim_state

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
    v_ew_n = AtmosphericData(env.atm, Ob).wind.v_ew_n
    v_eOb_n = v_ew_n + v_wOb_n

    Kinematics.Initializer(; q_nb, loc, h, ω_lb_b, v_eOb_n, Δx = 0.0, Δy = 0.0)

end

#assigns trim state and parameters to the aircraft system, then updates it
function assign!(ac::System{<:C172FBW.Template},
                env::System{<:AbstractEnvironment},
                trim_params::TrimParameters,
                trim_state::TrimState)

    @unpack TAS, β_a, fuel, flaps, mixture = trim_params
    @unpack n_eng, α_a, throttle, aileron, elevator, rudder = trim_state

    init_kinematics!(ac, Kinematics.Initializer(trim_state, trim_params, env))

    #for trimming, control surface inputs are set to zero, and we work only with
    #their offsets
    ac.airframe.act.u.throttle_cmd = throttle
    ac.airframe.act.u.elevator_cmd = 0
    ac.airframe.act.u.aileron_cmd = 0
    ac.airframe.act.u.rudder_cmd = 0
    ac.airframe.act.u.aileron_cmd_offset = aileron
    ac.airframe.act.u.elevator_cmd_offset = elevator
    ac.airframe.act.u.rudder_cmd_offset = rudder
    ac.airframe.act.u.flaps = flaps
    ac.airframe.act.u.mixture = mixture

    #engine must be running
    ac.airframe.pwp.engine.s.state = Piston.eng_running

    #set engine speed state
    ω_eng = n_eng * ac.airframe.pwp.engine.params.ω_rated
    ac.x.airframe.pwp.engine.ω = ω_eng

    #engine idle compensator: as long as the engine remains at normal
    #operational speeds, well above its nominal idle speed, the idle controller
    #compensator's output will be saturated at its lower bound by proportional
    #error. its integrator will be disabled, its state will not change nor have
    #any effect on the engine. we can simply set it to zero
    ac.x.airframe.pwp.engine.idle .= 0.0

    #engine friction compensator: with the engine running at normal operational
    #speeds, the engine's friction constraint compensator will be saturated, so
    #its integrator will be disabled and its state will not change. furthermore,
    #with the engine running friction is ignored. we can simply set it to zero.
    ac.x.airframe.pwp.engine.frc .= 0.0

    #actuator states: in steady state every actuator's velocity state must be
    #zero, and its position state must be equal to the actuator command. the
    #actuator command is in turn equal to the surface command plus its offset,
    #which we have set to zero
    ac.x.airframe.act.throttle_act.v = 0.0
    ac.x.airframe.act.throttle_act.p = throttle
    ac.x.airframe.act.aileron_act.v = 0.0
    ac.x.airframe.act.aileron_act.p = aileron
    ac.x.airframe.act.elevator_act.v = 0.0
    ac.x.airframe.act.elevator_act.p = elevator
    ac.x.airframe.act.rudder_act.v = 0.0
    ac.x.airframe.act.rudder_act.p = rudder

    ac.x.airframe.aero.α_filt = α_a #ensures zero state derivative
    ac.x.airframe.aero.β_filt = β_a #ensures zero state derivative
    ac.x.airframe.fuel .= fuel

    f_ode!(ac, env)

    #check assumptions concerning airframe systems states & derivatives
    @assert !any(SVector{3}(leg.strut.wow for leg in ac.airframe.ldg.y))
    @assert ac.x.airframe.pwp.engine.ω > ac.airframe.pwp.engine.params.ω_idle
    @assert ac.ẋ.airframe.pwp.engine.idle[1] .== 0
    @assert ac.ẋ.airframe.pwp.engine.frc[1] .== 0
    @assert abs(ac.ẋ.airframe.aero.α_filt) < 1e-10
    @assert abs(ac.ẋ.airframe.aero.β_filt) < 1e-10

    @assert all(SVector{8,Float64}(ac.ẋ.airframe.act) .== 0)

end

function cost(ac::System{<:C172FBW.Template})

    v_nd_dot = SVector{3}(ac.ẋ.kinematics.vel.v_eOb_b) / norm(ac.y.kinematics.common.v_eOb_b)
    ω_dot = SVector{3}(ac.ẋ.kinematics.vel.ω_eb_b) #ω should already of order 1
    n_eng_dot = ac.ẋ.airframe.pwp.engine.ω / ac.airframe.pwp.engine.params.ω_rated

    sum(v_nd_dot.^2) + sum(ω_dot.^2) + n_eng_dot^2

end

function get_f_target(ac::System{<:C172FBW.Template},
                      env::System{<:AbstractEnvironment},
                      trim_params::TrimParameters)

    let ac = ac, env = env, trim_params = trim_params
        function (x::TrimState)
            assign!(ac, env, trim_params, x)
            return cost(ac)
        end
    end

end


function Aircraft.trim!( ac::System{<:C172FBW.Template};
                env::System{<:AbstractEnvironment} = System(SimpleEnvironment()),
                trim_params::TrimParameters = TrimParameters())

    trim_state = TrimState() #could initial condition as an optional input

    f_target = get_f_target(ac, env, trim_params)

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
        α_a = ac.airframe.aero.params.α_stall[2], #critical AoA is 0.28 < 0.36
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
        println("Warning: Optimization failed with exit_flag $exit_flag")
    end
    trim_state_opt = TrimState(minx)
    assign!(ac, env, trim_params, trim_state_opt)
    return (success = success, result = trim_state_opt)


end

function Aircraft.trim!(
    world::System{<:SimpleWorld{<:C172FBW.Template, <:AbstractEnvironment}};
    trim_params::TrimParameters = TrimParameters())

    trim!(world.ac; env = world.env, trim_params = trim_params)

end

# ################################################################################
# ############################### Linearization ##################################

#labels corresponding to LinearX components within the overall
#Cessna172RBase{NED} state vector
const XLabels = (
        "kinematics.pos.ψ_nb", "kinematics.pos.θ_nb", "kinematics.pos.φ_nb",
        "kinematics.pos.ϕ", "kinematics.pos.λ", "kinematics.pos.h_e",
        "kinematics.vel.ω_eb_b[1]", "kinematics.vel.ω_eb_b[2]", "kinematics.vel.ω_eb_b[3]",
        "kinematics.vel.v_eOb_b[1]", "kinematics.vel.v_eOb_b[2]", "kinematics.vel.v_eOb_b[3]",
        "airframe.aero.α_filt", "airframe.aero.β_filt",
        "airframe.pwp.engine.ω", "airframe.fuel[1]",
        "airframe.act.throttle_act.v", "airframe.act.throttle_act.p",
        "airframe.act.aileron_act.v", "airframe.act.aileron_act.p",
        "airframe.act.elevator_act.v", "airframe.act.elevator_act.p",
        "airframe.act.rudder_act.v", "airframe.act.rudder_act.p",
    )

@kwdef mutable struct LinearX <: FieldVector{24, Float64}
    ψ::Float64 = 0.0 #heading
    θ::Float64 = 0.0 #inclination
    φ::Float64 = 0.0 #bank
    ϕ::Float64 = 0.0 #latitude
    λ::Float64 = 0.0 #longitude
    h::Float64 = 0.0 #ellipsoidal altitude
    p::Float64 = 0.0 #roll rate (ω_eb_b)
    q::Float64 = 0.0 #pitch rate (ω_eb_b)
    r::Float64 = 0.0 #yaw rate (ω_eb_b)
    v_x::Float64 = 0.0 #Ob/ECEF velocity, x-body
    v_y::Float64 = 0.0 #Ob/ECEF velocity, y-body
    v_z::Float64 = 0.0 #Ob/ECEF velocity, z-body
    α_filt::Float64 = 0.0 #filtered AoA
    β_filt::Float64 = 0.0 #filtered AoS
    ω_eng::Float64 = 0.0 #engine speed
    fuel::Float64 = 0.0 #fuel fraction
    thr_v::Float64 = 0.0 #throttle actuator velocity
    thr_p::Float64 = 0.0 #throttle actuator position
    ail_v::Float64 = 0.0 #aileron actuator velocity
    ail_p::Float64 = 0.0 #aileron actuator position
    ele_v::Float64 = 0.0 #elevator actuator velocity
    ele_p::Float64 = 0.0 #elevator actuator position
    rud_v::Float64 = 0.0 #rudder actuator velocity
    rud_p::Float64 = 0.0 #rudder actuator position
end

#flaps and mixture are omitted from the control vector and treated as parameters
@kwdef mutable struct LinearU <: FieldVector{4, Float64}
    throttle_cmd::Float64 = 0.0
    aileron_cmd::Float64 = 0.0
    elevator_cmd::Float64 = 0.0
    rudder_cmd::Float64 = 0.0
end


# @kwdef mutable struct LinearY <: FieldVector
#     ψ = 0.0, θ = 0.0, φ = 0.0, #heading, inclination, bank (body/NED)
#     ϕ = 0.0, λ = 0.0, h = 0.0, #latitude, longitude, ellipsoidal altitude
#     p = 0.0, q = 0.0, r = 0.0, #angular rates (ω_eb_b)
#     TAS = 0.0, α = 0.0, β = 0.0, #airspeed, AoA, AoS
#     f_x = 0.0, f_y = 0.0, f_z = 0.0, #specific force at G (f_iG_b)
#     ω_eng = 0.0, m_fuel = 0.0 #engine speed, fuel mass
#     v_N = 0.0, v_E = 0.0, v_D = 0.0 #Ob/ECEF velocity, NED axes
# end


# function assign!(u::LinearU, ac::System{<:Cessna172RBase})

#     @unpack throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd = ac.airframe.act.u
#     @pack! u = throttle, aileron, elevator, rudder

# end

# function assign!(ac::System{<:Cessna172RBase}, u::LinearU)

#     @unpack throttle, aileron, elevator, rudder = u
#     @pack! ac.airframe.act.u = throttle, aileron, elevator, rudder

# end

# function assign!(y::LinearY, ac::System{<:Cessna172RBase})

#     @unpack q_nb, n_e, h_e, ω_eb_b = ac.y.kinematics
#     @unpack α, β = ac.y.airframe.aero
#     @unpack ψ, θ, φ = REuler(q_nb)
#     @unpack ϕ, λ = LatLon(n_e)

#     h = h_e
#     p, q, r = ω_eb_b
#     f_x, f_y, f_z = ac.y.rigidbody.f_G_b
#     TAS = ac.y.air.TAS
#     ω_eng = ac.y.airframe.pwp.engine.ω
#     m_fuel = ac.y.airframe.fuel.m_avail

#     @pack! y = ψ, θ, φ, ϕ, λ, h, p, q, r, TAS, α, β, f_x, f_y, f_z, ω_eng, m_fuel

# end

# function Aircraft.linearize!(ac::System{<:Cessna172RBase{NED}};
#     env::System{<:AbstractEnvironment} = System(SimpleEnvironment()),
#     trim_params::TrimParameters = TrimParameters())

#     (_, trim_state) = trim!(ac; env, trim_params)

#     #save the trimmed aircraft's ẋ, x, u and y for later
#     ẋ0_full = copy(ac.ẋ)
#     x0_full = copy(ac.x)
#     u0 = similar(LinearUTemplate); assign!(u0, ac) #get reference value from trimmed aircraft
#     y0 = similar(LinearYTemplate); assign!(y0, ac) #idem

#     #function wrapper around f_ode!(), mutates ẋ and y.
#     f_nonlinear! = let ac = ac, env = env,
#                        trim_params = trim_params, trim_state = trim_state,
#                        u_axes = getaxes(u0), y_axes = getaxes(y0)

#         function (ẋ, y, x, u)

#             # cast y and u into ComponentVectors in case we get generic Vectors
#             # from FiniteDiff. these do not allocate, because the underlying
#             #data is already in u and y
#             u_cv = ComponentVector(u, u_axes)
#             y_cv = ComponentVector(y, y_axes)

#             #make sure any input or state not set by x and u is at its reference
#             #trim value. this reverts any potential changes to the aircraft done
#             #by functions sharing the same aircraft instance
#             assign!(ac, env, trim_params, trim_state)

#             assign!(ac, u_cv)
#             ac.x .= x
#             f_ode!(ac, env)

#             ẋ .= ac.ẋ
#             assign!(y_cv, ac) #this also updates y (shares its data with y_cv)

#         end

#     end

#     (A_full, B_full, C_full, D_full) = ss_matrices(f_nonlinear!;
#                                             ẋ0 = ẋ0_full, y0, x0 = x0_full, u0)

#     #once we're done, ensure the aircraft is restored to its trimmed status, so
#     #the response can be compared with that of its linear counterpart
#     assign!(ac, env, trim_params, trim_state)

#     #find the indices for the components in the reduced LinearX state vector
#     #within the overall aircraft state vector
#     x_indices = [ComponentArrays.label2index(x0_full, s)[1] for s in x_labels]

#     x_axis = getaxes(LinearXTemplate)[1]

#     #extract the required elements from ẋ0_full and x0_full and rebuild them
#     #with the LinearX axis
#     ẋ0 = ComponentVector(ẋ0_full[x_indices], x_axis)
#     x0 = ComponentVector(x0_full[x_indices], x_axis)

#     #extract the required rows and columns from A, the required rows from B, and
#     #the required columns from C. then rebuild them with the new x_axis and the
#     #previous u_axis and y_axis
#     A = ComponentMatrix(A_full[x_indices, x_indices], x_axis, x_axis)
#     B = ComponentMatrix(B_full[x_indices, :], x_axis, getaxes(u0)[1])
#     C = ComponentMatrix(C_full[:, x_indices], getaxes(y0)[1], x_axis)
#     D = D_full

#     return LinearStateSpace(ẋ0, x0, u0, y0, A, B, C, D)

# end


# function ss_matrices(f_nonlinear!::Function; ẋ0, y0, x0, u0)

#     f_A! = let u = u0, y = similar(y0) #y is discarded
#         (ẋ, x) -> f_nonlinear!(ẋ, y, x, u)
#     end

#     f_B! = let x = x0, y = similar(y0) #y is discarded
#         (ẋ, u) -> f_nonlinear!(ẋ, y, x, u)
#     end

#     f_C! = let u = u0, ẋ = similar(ẋ0) #ẋ is discarded
#         (y, x) -> f_nonlinear!(ẋ, y, x, u)
#     end

#     f_D! = let x = x0, ẋ = similar(ẋ0) #ẋ is discarded
#         (y, u) -> f_nonlinear!(ẋ, y, x, u)
#     end

#     #preallocate
#     A = x0 * x0'
#     B = x0 * u0'
#     C = y0 * x0'
#     D = y0 * u0'

#     jacobian!(A, f_A!, x0)
#     jacobian!(B, f_B!, u0)
#     jacobian!(C, f_C!, x0)
#     jacobian!(D, f_D!, u0)

#     return (A, B, C, D)

# end

################################################################################
################################## Variants ####################################

include(normpath("variants/base.jl")); @reexport using .C172FBWBase
# include(normpath("variants/cas.jl")); @reexport using .C172FBWCAS

end