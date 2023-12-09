module C172FBW

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, Reexport
using ControlSystems, RobustAndOptimalControl

using Flight.FlightCore
using Flight.FlightCore.Utils

using Flight.FlightPhysics
using Flight.FlightComponents

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

#with an underdamped actuator, the position state can still transiently exceed
#the intended range due to overshoot. the true actuator position should
#therefore be clamped. in the real world, this behaviour could correspond to a
#clutched output actuator, where the output position saturates beyond a given
#opposing torque (for example, if the surface's mechanical limits are hit)

#saturate on command, not on position, which only tends asymptotically to cmd!

struct Actuator{R} <: SystemDefinition #second order linear actuator model
    ω_n::Float64 #natural frequency (default: 10 Hz)
    ζ::Float64 #damping ratio (default: underdamped with minimal resonance)
    function Actuator(; ω_n::Real = 5*2π, ζ::Real = 0.6,
                        range::Tuple{Real, Real} = (-1.0, 1.0))
        new{Ranged{Float64, range[1], range[2]}}(ω_n, ζ)
    end
end

@kwdef struct ActuatorY{R}
    cmd::R = R(0.0)
    pos::R = R(0.0)
    vel::Float64 = 0.0
    sat::Int64 = 0
end

Systems.init(::SystemX, ::Actuator) = ComponentVector(v = 0.0, p = 0.0)
Systems.init(::SystemU, ::Actuator{R}) where {R} = Ref(R(0.0))
Systems.init(::SystemY, ::Actuator{R}) where {R} = ActuatorY{R}()

function Systems.f_ode!(sys::System{Actuator{R}}) where {R}

    @unpack ẋ, x, u, constants = sys
    @unpack ω_n, ζ = constants

    cmd = u[]
    pos = R(x.p)
    vel = x.v
    sat = saturation(cmd)

    ẋ.v = ω_n^2 * (Float64(cmd) - x.p) - 2ζ*ω_n*x.v
    ẋ.p = x.v

    sys.y = ActuatorY(; cmd, pos, vel, sat)

end

function GUI.draw(sys::System{<:Actuator})
    CImGui.Text("Actuator Command"); CImGui.SameLine(200); display_bar("", sys.y.cmd)
    CImGui.Text("Actuator Position"); CImGui.SameLine(200); display_bar("", sys.y.pos)
    @running_plot("Actuator Position", Float64(sys.y.pos), typemin(sys.y.pos), typemax(sys.y.pos), 0.0, 120)
end

function GUI.draw!(sys::System{<:Actuator})
    sys.u[] = safe_slider("Actuator Command", sys.u[], "%.6f")
    CImGui.Text("Actuator Command"); CImGui.SameLine(150); display_bar("", sys.y.cmd)
    CImGui.Text("Actuator Position"); CImGui.SameLine(150); display_bar("", sys.y.pos)
    @running_plot("Actuator Position", Float64(sys.y.pos), typemin(sys.y.pos), typemax(sys.y.pos), 0.0, 120)
end

################################################################################
#################################### Actuation #################################

#Fly-by-wire actuation system. Throttle, steering and aerodynamic surfaces are
#controlled via actuators, the rest of are direct feedthrough

@kwdef struct Actuation <: C172.Actuation
    throttle::Actuator = Actuator(range = (0.0, 1.0))
    aileron::Actuator = Actuator(range = (-1.0, 1.0))
    elevator::Actuator = Actuator(range = (-1.0, 1.0))
    rudder::Actuator = Actuator(range = (-1.0, 1.0))
    flaps::Actuator = Actuator(range = (0.0, 1.0))
    steering::Actuator = Actuator(range = (-1.0, 1.0))
end

function C172.assign!(aero::System{<:C172.Aero},
                    ldg::System{<:C172.Ldg},
                    pwp::System{<:Piston.Thruster},
                    act::System{<:Actuation})

    @unpack throttle, aileron, elevator, rudder, flaps, steering = act.y

    aero.u.e = -elevator.pos
    aero.u.a = aileron.pos
    aero.u.r = -rudder.pos
    aero.u.f = flaps.pos
    pwp.engine.u.throttle = throttle.pos
    ldg.nose.steering.u[] = steering.pos

end

function GUI.draw(sys::System{Actuation}, label::String = "Cessna 172 Fly-By-Wire Actuation")

    CImGui.Begin(label)
    CImGui.PushItemWidth(-60)

    foreach(keys(sys.subsystems), values(sys.subsystems)) do k, ss
        if CImGui.CollapsingHeader(uppercasefirst(string(k)))
            GUI.draw(ss)
            CImGui.Dummy(10.0, 10.0);
        end
    end

    CImGui.PopItemWidth()
    CImGui.End()

end

function GUI.draw!(sys::System{Actuation}, label::String = "Cessna 172 Fly-By-Wire Actuation")

    CImGui.Begin(label)
    CImGui.PushItemWidth(-60)

    foreach(keys(sys.subsystems), values(sys.subsystems)) do k, ss
        if CImGui.CollapsingHeader(uppercasefirst(string(k)))
            GUI.draw!(ss)
            CImGui.Dummy(10.0, 10.0);
        end
    end

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

    u.flaps += 0.3333 * was_released(joystick, :button_3)
    u.flaps -= 0.3333 * was_released(joystick, :button_2)

end


################################################################################
################################# Templates ####################################

const Airframe = C172.Airframe{typeof(PowerPlant()), C172FBW.Actuation}
const Physics{K, T} = Aircraft.Physics{C172FBW.Airframe, K, T} where {K <: AbstractKinematicDescriptor, T <: AbstractTerrain}
const Template{K, T, A} = Aircraft.Template{C172FBW.Physics{K, T}, A} where {K <: AbstractKinematicDescriptor, T <: AbstractTerrain, A <: AbstractAvionics}

function Physics(kinematics = LTF(), terrain = HorizontalTerrain())
    Aircraft.Physics(C172.Airframe(PowerPlant(), Actuation()), kinematics, terrain, LocalAtmosphere())
end

function Template(kinematics = LTF(), terrain = HorizontalTerrain(), avionics = NoAvionics())
    Aircraft.Template(Physics(kinematics, terrain), avionics)
end

############################### Trimming #######################################
################################################################################

#assigns trim state and parameters to aircraft physics, then updates aircraft physics
function Aircraft.assign!(physics::System{<:C172FBW.Physics},
                        trim_params::C172.TrimParameters,
                        trim_state::C172.TrimState)

    @unpack EAS, β_a, x_fuel, flaps, mixture, payload = trim_params
    @unpack n_eng, α_a, throttle, aileron, elevator, rudder = trim_state
    @unpack act, pwp, aero, fuel, ldg, pld = physics.airframe

    atm_data = LocalAtmosphericData(physics.atmosphere)
    Systems.init!(physics.kinematics, Kinematics.Initializer(trim_state, trim_params, atm_data))

    act.throttle.u[] = throttle
    act.aileron.u[] = aileron
    act.elevator.u[] = elevator
    act.rudder.u[] = rudder
    act.flaps.u[] = flaps

    #assign payload
    @unpack m_pilot, m_copilot, m_lpass, m_rpass, m_baggage = payload
    @pack! pld.u = m_pilot, m_copilot, m_lpass, m_rpass, m_baggage

    #engine must be running
    pwp.engine.s.state = Piston.eng_running

    pwp.engine.u.mixture = mixture

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
    #zero, and its position state must be equal to the actuator command.
    act.x.throttle.v = 0.0
    act.x.aileron.v = 0.0
    act.x.elevator.v = 0.0
    act.x.rudder.v = 0.0
    act.x.flaps.v = 0.0

    act.x.throttle.p = throttle
    act.x.aileron.p = aileron
    act.x.elevator.p = elevator
    act.x.rudder.p = rudder
    act.x.flaps.p = flaps

    aero.x.α_filt = α_a #ensures zero state derivative
    aero.x.β_filt = β_a #ensures zero state derivative
    fuel.x .= Float64(x_fuel)

    f_ode!(physics)

    #check essential assumptions about airframe systems states & derivatives
    @assert !any(SVector{3}(leg.strut.wow for leg in ldg.y))
    @assert pwp.x.engine.ω > pwp.engine.constants.ω_idle
    @assert pwp.x.engine.idle[1] .== 0
    @assert pwp.x.engine.frc[1] .== 0
    @assert abs(aero.ẋ.α_filt) < 1e-10
    @assert abs(aero.ẋ.β_filt) < 1e-10

    @assert all(SVector{12,Float64}(act.ẋ) .== 0)

end



################################################################################
############################### Linearization ##################################

@kwdef struct XLinear <: FieldVector{24, Float64}
    p::Float64 = 0.0; q::Float64 = 0.0; r::Float64 = 0.0; #angular rates (ω_eb_b)
    ψ::Float64 = 0.0; θ::Float64 = 0.0; φ::Float64 = 0.0; #heading, inclination, bank (body/NED)
    v_x::Float64 = 0.0; v_y::Float64 = 0.0; v_z::Float64 = 0.0; #aerodynamic velocity, body axes
    ϕ::Float64 = 0.0; λ::Float64 = 0.0; h::Float64 = 0.0; #latitude, longitude, ellipsoidal altitude
    α_filt::Float64 = 0.0; β_filt::Float64 = 0.0; #filtered airflow angles
    ω_eng::Float64 = 0.0; fuel::Float64 = 0.0; #engine speed, fuel fraction
    thr_v::Float64 = 0.0; thr_p::Float64 = 0.0; #throttle actuator states
    ail_v::Float64 = 0.0; ail_p::Float64 = 0.0; #aileron actuator states
    ele_v::Float64 = 0.0; ele_p::Float64 = 0.0; #elevator actuator states
    rud_v::Float64 = 0.0; rud_p::Float64 = 0.0 #rudder actuator states
end

#flaps and mixture are trim parameters and thus omitted from the control vector
@kwdef struct ULinear <: FieldVector{4, Float64}
    throttle_cmd::Float64 = 0.0
    aileron_cmd::Float64 = 0.0
    elevator_cmd::Float64 = 0.0
    rudder_cmd::Float64 = 0.0
end

#all states (for full-state feedback), plus other useful stuff, plus control inputs
@kwdef struct YLinear <: FieldVector{41, Float64}
    p::Float64 = 0.0; q::Float64 = 0.0; r::Float64 = 0.0; #angular rates (ω_eb_b)
    ψ::Float64 = 0.0; θ::Float64 = 0.0; φ::Float64 = 0.0; #heading, inclination, bank (body/NED)
    v_x::Float64 = 0.0; v_y::Float64 = 0.0; v_z::Float64 = 0.0; #aerodynamic velocity, body axes
    ϕ::Float64 = 0.0; λ::Float64 = 0.0; h::Float64 = 0.0; #latitude, longitude, ellipsoidal altitude
    α_filt::Float64 = 0.0; β_filt::Float64 = 0.0; #filtered airflow angles
    ω_eng::Float64 = 0.0; fuel::Float64 = 0.0; #engine speed, available fuel fraction
    thr_v::Float64 = 0.0; thr_p::Float64 = 0.0; #throttle actuator states
    ail_v::Float64 = 0.0; ail_p::Float64 = 0.0; #aileron actuator states
    ele_v::Float64 = 0.0; ele_p::Float64 = 0.0; #elevator actuator states
    rud_v::Float64 = 0.0; rud_p::Float64 = 0.0; #rudder actuator states
    f_x::Float64 = 0.0; f_y::Float64 = 0.0; f_z::Float64 = 0.0; #specific force at G (f_iG_b)
    α::Float64 = 0.0; β::Float64 = 0.0; #unfiltered airflow angles
    EAS::Float64 = 0.0; TAS::Float64 = 0.0; #airspeed
    v_N::Float64 = 0.0; v_E::Float64 = 0.0; v_D::Float64 = 0.0; #Ob/ECEF velocity, NED axes
    χ::Float64 = 0.0; γ::Float64 = 0.0; climb_rate::Float64 = 0.0; #track and flight path angles, climb rate
    throttle_cmd::Float64 = 0.0; aileron_cmd::Float64 = 0.0; #actuator commands
    elevator_cmd::Float64 = 0.0; rudder_cmd::Float64 = 0.0; #actuator commands
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
    thr_v = x_airframe.act.throttle.v
    thr_p = x_airframe.act.throttle.p
    ail_v = x_airframe.act.aileron.v
    ail_p = x_airframe.act.aileron.p
    ele_v = x_airframe.act.elevator.v
    ele_p = x_airframe.act.elevator.p
    rud_v = x_airframe.act.rudder.v
    rud_p = x_airframe.act.rudder.p

    ψ, θ, φ, h = ψ_nb, θ_nb, φ_nb, h_e

    XLinear(;  p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt,
        ω_eng, fuel, thr_v, thr_p, ail_v, ail_p, ele_v, ele_p, rud_v, rud_p)

end

function ULinear(physics::System{<:C172FBW.Physics{NED}})

    @unpack throttle, aileron, elevator, rudder = physics.airframe.act.subsystems
    throttle_cmd = throttle.u[]
    aileron_cmd = aileron.u[]
    elevator_cmd = elevator.u[]
    rudder_cmd = rudder.u[]
    ULinear(; throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd)

end

function YLinear(physics::System{<:C172FBW.Physics{NED}})

    @unpack airframe, air, rigidbody, kinematics = physics.y
    @unpack pwp, fuel, aero,act = airframe

    @unpack e_nb, ϕ_λ, h_e, ω_eb_b, v_eOb_b, v_eOb_n, χ_gnd, γ_gnd = kinematics
    @unpack ψ, θ, φ = e_nb
    @unpack ϕ, λ = ϕ_λ

    h = h_e
    p, q, r = ω_eb_b
    v_x, v_y, v_z = v_eOb_b
    v_N, v_E, v_D = v_eOb_n
    ω_eng = pwp.engine.ω
    fuel = fuel.x_avail
    α_filt = aero.α_filt
    β_filt = aero.β_filt

    thr_v = act.throttle.vel
    thr_p = act.throttle.pos
    ail_v = act.aileron.vel
    ail_p = act.aileron.pos
    ele_v = act.elevator.vel
    ele_p = act.elevator.pos
    rud_v = act.rudder.vel
    rud_p = act.rudder.pos

    f_x, f_y, f_z = physics.y.rigidbody.f_G_b
    EAS = physics.y.air.EAS
    TAS = physics.y.air.TAS
    α = physics.y.air.α_b
    β = physics.y.air.β_b
    χ = χ_gnd
    γ = γ_gnd
    climb_rate = -v_D

    @unpack throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd = ULinear(physics)

    YLinear(; p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt,
            ω_eng, fuel, thr_v, thr_p, ail_v, ail_p, ele_v, ele_p, rud_v, rud_p,
            f_x, f_y, f_z, EAS, TAS, α, β, v_N, v_E, v_D, χ, γ, climb_rate,
            throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd)

end

Aircraft.ẋ_linear(physics::System{<:C172FBW.Physics{NED}}) = XLinear(physics.ẋ)
Aircraft.x_linear(physics::System{<:C172FBW.Physics{NED}}) = XLinear(physics.x)
Aircraft.u_linear(physics::System{<:C172FBW.Physics{NED}}) = ULinear(physics)
Aircraft.y_linear(physics::System{<:C172FBW.Physics{NED}}) = YLinear(physics)

function Aircraft.assign_u!(physics::System{<:C172FBW.Physics{NED}}, u::AbstractVector{Float64})

    #The velocity states in the linearized model are meant to be aerodynamic so
    #they can be readily used for flight control design. Since the velocity
    #states in the nonlinear model are Earth-relative, we need to ensure wind
    #velocity is set to zero for linearization.
    @unpack throttle, aileron, elevator, rudder = physics.airframe.act.subsystems
    @unpack throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd = ULinear(u)
    throttle.u[] = throttle_cmd
    aileron.u[] = aileron_cmd
    elevator.u[] = elevator_cmd
    rudder.u[] = rudder_cmd

    physics.atmosphere.u.v_ew_n .= 0

end

function Aircraft.assign_x!(physics::System{<:C172FBW.Physics{NED}}, x::AbstractVector{Float64})

    @unpack p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng,
            fuel, thr_v, thr_p, ail_v, ail_p, ele_v, ele_p, rud_v, rud_p = XLinear(x)

    x_kinematics = physics.x.kinematics
    x_airframe = physics.x.airframe

    ψ_nb, θ_nb, φ_nb, h_e = ψ, θ, φ, h

    @pack! x_kinematics.pos = ψ_nb, θ_nb, φ_nb, ϕ, λ, h_e
    x_kinematics.vel.ω_eb_b .= p, q, r
    x_kinematics.vel.v_eOb_b .= v_x, v_y, v_z
    x_airframe.aero .= α_filt, β_filt
    x_airframe.pwp.engine.ω = ω_eng
    x_airframe.fuel .= fuel
    x_airframe.act.throttle.v = thr_v
    x_airframe.act.throttle.p = thr_p
    x_airframe.act.aileron.v = ail_v
    x_airframe.act.aileron.p = ail_p
    x_airframe.act.elevator.v = ele_v
    x_airframe.act.elevator.p = ele_p
    x_airframe.act.rudder.v = rud_v
    x_airframe.act.rudder.p = rud_p

end

function Control.Continuous.LinearizedSS(
            physics::System{<:C172FBW.Physics{NED}},
            trim_params::C172.TrimParameters = C172.TrimParameters();
            model::Symbol = :full)

    lm = linearize!(physics, trim_params)

    if model === :full
        return lm

    #preserve the ordering of the complete linearized state and output vectors
    elseif model === :lon
        x_labels = [:q, :θ, :v_x, :v_z, :h, :α_filt, :ω_eng, :thr_v, :thr_p, :ele_v, :ele_p]
        u_labels = [:throttle_cmd, :elevator_cmd]
        y_labels = vcat(x_labels, [:f_x, :f_z, :α, :EAS, :TAS, :γ, :climb_rate, :throttle_cmd, :elevator_cmd])
        return Control.Continuous.submodel(lm; x = x_labels, u = u_labels, y = y_labels)

    elseif model === :lat
        x_labels = [:p, :r, :ψ, :φ, :v_x, :v_y, :β_filt, :ail_v, :ail_p, :rud_v, :rud_p]
        u_labels = [:aileron_cmd, :rudder_cmd]
        y_labels = vcat(x_labels, [:f_y, :β, :χ, :aileron_cmd, :rudder_cmd])
        return Control.Continuous.submodel(lm; x = x_labels, u = u_labels, y = y_labels)

    else
        error("Valid model keyword values: :full, :lon, :lat")

    end

end


################################################################################
################################## Variants ####################################

include(normpath("variants/base.jl")); @reexport using .C172FBWBase
include(normpath("variants/cas/cas.jl")); @reexport using .C172CAS
include(normpath("variants/mcs/mcs.jl")); @reexport using .C172MCS

end