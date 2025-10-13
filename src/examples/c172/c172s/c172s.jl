module C172S

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, Reexport
using ControlSystems, RobustAndOptimalControl
using NLopt
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore
using Flight.FlightLib

using ..C172

export Cessna172S, Cessna172Sv0

################################################################################
################################ Powerplant ####################################

function PowerPlant()

    engine = PistonEngine(
        P_rated= Piston.hp2W(200), #W
        ω_rated = Piston.RPM2radpersec(2700),
        ω_stall = Piston.RPM2radpersec(300),
        ω_max = Piston.RPM2radpersec(3100),
        ω_idle = Piston.RPM2radpersec(600),
        τ_start = 40, #N·m
        J = 0.05, #kg·m²
    )

    propeller = Propeller(FixedPitch(), 2;
        sense = Propellers.CW, d = 2.0, J_xx = 0.3,
        t_bp = FrameTransform(r = [2.055, 0, 0.833]))

    PistonThruster(; engine, propeller)

end

################################################################################
################### Mechanical Actuation System for C172S ######################

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

struct MechanicalActuation <: C172.AbstractActuation end

@kwdef mutable struct MechanicalActuationU
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

@kwdef struct MechanicalActuationY
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

Modeling.U(::MechanicalActuation) = MechanicalActuationU()
Modeling.Y(::MechanicalActuation) = MechanicalActuationY()

function Modeling.f_ode!(act::Model{<:MechanicalActuation})

    @unpack aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset,
            flaps, brake_left, brake_right= act.u

    act.y = MechanicalActuationY(; aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset,
            flaps, brake_left, brake_right)

end

function C172.assign!(aero::Model{<:C172.Aero},
                ldg::Model{<:C172.Ldg},
                pwp::Model{<:PistonThruster},
                act::Model{<:MechanicalActuation})

    @unpack aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset,
            brake_left, brake_right, flaps = act.y

    ldg.nose.steering.u.input = (rudder_offset + rudder)
    ldg.left.braking.u[] = brake_left
    ldg.right.braking.u[] = brake_right
    aero.u.e = -(elevator_offset + elevator)
    aero.u.a = (aileron_offset + aileron)
    aero.u.r = -(rudder_offset + rudder)
    aero.u.f = flaps

    return nothing
end


function GUI.draw!(mdl::Model{<:MechanicalActuation}, p_open::Ref{Bool} = Ref(true),
                label::String = "Cessna 172S Actuation")

    u = mdl.u

    CImGui.Begin(label, p_open)

    CImGui.PushItemWidth(-150)

    u.aileron = safe_slider("Aileron", u.aileron, "%.6f")
    u.elevator = safe_slider("Elevator", u.elevator, "%.6f")
    u.rudder = safe_slider("Rudder", u.rudder, "%.6f")
    u.aileron_offset = safe_input("Aileron Offset", u.aileron_offset, 0.001, 0.1, "%.6f")
    u.elevator_offset = safe_input("Elevator Offset", u.elevator_offset, 0.001, 0.1, "%.6f")
    u.rudder_offset = safe_input("Rudder Offset", u.rudder_offset, 0.001, 0.1, "%.6f")
    u.flaps = safe_slider("Flaps", u.flaps, "%.6f")
    u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
    u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")

    CImGui.PopItemWidth()

    CImGui.End()

end


################################################################################
################################# Templates ####################################

#needed for dispatching
const Systems = C172.Systems{typeof(PowerPlant()), typeof(MechanicalActuation())}
const Vehicle{K} = AircraftBase.Vehicle{Systems, K} where {K <: AbstractKinematicDescriptor}
const Aircraft{K, A} = AircraftBase.Aircraft{Vehicle{K}, A} where {K <: AbstractKinematicDescriptor, A <: AbstractAvionics}
const Cessna172S{K, A} = Aircraft{K, A}

function Vehicle(kinematics = WA())
    AircraftBase.Vehicle(
        C172.Systems(PowerPlant(), MechanicalActuation()),
        kinematics, VehicleDynamics())
end


############################## Initialization ##################################
################################################################################

function Modeling.init!(sys::Model{<:Systems}, init::C172.SystemsInitializer)

    @unpack act, pwp, aero, fuel, ldg, pld = sys

    @unpack engine_state, mixture_ctl, n_eng, throttle, mixture, elevator,
    aileron, rudder, flaps, brake_left, brake_right, fuel_load, payload,
    stall, α_a_filt, β_a_filt = init

    @unpack m_pilot, m_copilot, m_lpass, m_rpass, m_baggage = payload

    #assign payload
    @pack! pld.u = m_pilot, m_copilot, m_lpass, m_rpass, m_baggage

    pwp.engine.u.throttle = throttle
    pwp.engine.u.mixture = mixture
    pwp.engine.u.mixture_ctl = mixture_ctl
    act.u.elevator = elevator
    act.u.aileron = aileron
    act.u.rudder = rudder
    act.u.flaps = flaps
    act.u.aileron_offset = 0
    act.u.elevator_offset = 0
    act.u.rudder_offset = 0
    act.u.brake_left = brake_left
    act.u.brake_right = brake_right

    pwp.engine.s[] = engine_state
    pwp.x.engine.ω = n_eng * pwp.engine.ω_rated

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

    aero.s.stall = stall
    aero.x.α_filt = α_a_filt
    aero.x.β_filt = β_a_filt

    fuel.x .= Float64(fuel_load)

end


############################### Trimming #######################################
################################################################################

#assemble initializer from trim state and parameters, then initialize vehicle
function AircraftBase.assign!(vehicle::Model{<:C172S.Vehicle},
                        trim_params::C172.TrimParameters,
                        trim_state::C172.TrimState,
                        atmosphere::Model{<:AbstractAtmosphere},
                        terrain::Model{<:AbstractTerrain})

    @unpack β_a, fuel_load, flaps, mixture, payload = trim_params
    @unpack n_eng, α_a, throttle, aileron, elevator, rudder = trim_state
    @unpack pwp, aero, ldg = vehicle.systems

    kin_init = KinInit(trim_state, trim_params, atmosphere)

    sys_init = C172.SystemsInitializer(;
        engine_state = Piston.EngineState.running,
        mixture_ctl = Piston.MixtureControl.auto,
        n_eng, throttle, mixture, elevator, aileron, rudder,
        flaps, brake_left = 0, brake_right = 0, fuel_load, payload,
        stall = false,
        α_a_filt = α_a, #ensure zero α_a_filt state derivative
        β_a_filt = β_a #ensure zero β_a_filt state derivative
        )

    vehicle_init = C172.Init(kin_init, sys_init)

    #initialize the vehicle with the setup above. this will call f_ode!
    #internally, no need to do it here
    Modeling.init!(vehicle, vehicle_init, atmosphere, terrain)

    #sanity checks for component states & derivatives
    @assert !any(SVector{3}(leg.strut.wow for leg in ldg.y))
    @assert pwp.x.engine.ω > pwp.engine.ω_idle
    @assert pwp.x.engine.idle[1] .== 0
    @assert pwp.x.engine.frc[1] .== 0
    @assert abs(aero.ẋ.α_filt) < 1e-10
    @assert abs(aero.ẋ.β_filt) < 1e-10

end

################################################################################
############################### Linearization ##################################

@kwdef struct XStateSpace <: FieldVector{16, Float64}
    p::Float64 = 0.0; q::Float64 = 0.0; r::Float64 = 0.0; #angular rates (ω_eb_b)
    ψ::Float64 = 0.0; θ::Float64 = 0.0; φ::Float64 = 0.0; #heading, inclination, bank (body/NED)
    v_x::Float64 = 0.0; v_y::Float64 = 0.0; v_z::Float64 = 0.0; #aerodynamic velocity, body axes
    ϕ::Float64 = 0.0; λ::Float64 = 0.0; h::Float64 = 0.0; #latitude, longitude, ellipsoidal altitude
    α_filt::Float64 = 0.0; β_filt::Float64 = 0.0; #filtered airflow angles
    ω_eng::Float64 = 0.0; fuel::Float64 = 0.0; #engine speed, fuel fraction
end

#flaps and mixture are trim parameters and thus omitted from the control vector
@kwdef struct UStateSpace <: FieldVector{4, Float64}
    throttle::Float64 = 0.0
    aileron::Float64 = 0.0
    elevator::Float64 = 0.0
    rudder::Float64 = 0.0
end

#all states (for full-state feedback), plus other useful stuff, plus control inputs
@kwdef struct YStateSpace <: FieldVector{33, Float64}
    p::Float64 = 0.0; q::Float64 = 0.0; r::Float64 = 0.0; #angular rates (ω_eb_b)
    ψ::Float64 = 0.0; θ::Float64 = 0.0; φ::Float64 = 0.0; #heading, inclination, bank (body/NED)
    v_x::Float64 = 0.0; v_y::Float64 = 0.0; v_z::Float64 = 0.0; #aerodynamic velocity, body axes
    ϕ::Float64 = 0.0; λ::Float64 = 0.0; h::Float64 = 0.0; #latitude, longitude, ellipsoidal altitude
    α_filt::Float64 = 0.0; β_filt::Float64 = 0.0; #filtered airflow angles
    ω_eng::Float64 = 0.0; fuel::Float64 = 0.0; #engine speed, available fuel fraction
    f_x::Float64 = 0.0; f_y::Float64 = 0.0; f_z::Float64 = 0.0; #specific force at G (f_iG_b)
    α::Float64 = 0.0; β::Float64 = 0.0; #unfiltered airflow angles
    EAS::Float64 = 0.0; TAS::Float64 = 0.0; #airspeed
    v_N::Float64 = 0.0; v_E::Float64 = 0.0; v_D::Float64 = 0.0; #b/ECEF velocity, NED axes
    χ::Float64 = 0.0; γ::Float64 = 0.0; c::Float64 = 0.0; #track and flight path angles, climb rate
    throttle_out::Float64 = 0.0; aileron_out::Float64 = 0.0; #control inputs
    elevator_out::Float64 = 0.0; rudder_out::Float64 = 0.0; #control inputs
end


function XStateSpace(x_vehicle::ComponentVector)

    x_kinematics = x_vehicle.kinematics
    x_dynamics = x_vehicle.dynamics
    x_systems = x_vehicle.systems

    @unpack ψ_nb, θ_nb, φ_nb, ϕ, λ, h_e = x_kinematics
    p, q, r = x_dynamics.ω_eb_b
    v_x, v_y, v_z = x_dynamics.v_eb_b
    α_filt, β_filt = x_systems.aero
    ω_eng = x_systems.pwp.engine.ω
    fuel = x_systems.fuel[1]
    ψ, θ, φ, h = ψ_nb, θ_nb, φ_nb, h_e

    XStateSpace(; p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng, fuel)

end

function UStateSpace(vehicle::Model{<:C172S.Vehicle{NED}})

    @unpack aileron, elevator, rudder = vehicle.systems.act.u
    @unpack throttle = vehicle.systems.pwp.engine.u
    UStateSpace(; throttle, aileron, elevator, rudder)

end

function YStateSpace(vehicle::Model{<:C172S.Vehicle{NED}})

    @unpack aileron, elevator, rudder = vehicle.systems.act.u
    @unpack throttle = vehicle.systems.pwp.engine.u
    @unpack systems, airflow, dynamics, kinematics = vehicle.y
    @unpack pwp, fuel, aero, act = systems

    @unpack e_nb, ϕ_λ, h_e, ω_eb_b, v_eb_b, v_eb_n, χ_gnd, γ_gnd = kinematics
    @unpack ψ, θ, φ = e_nb
    @unpack ϕ, λ = ϕ_λ

    h = h_e
    p, q, r = ω_eb_b
    v_x, v_y, v_z = v_eb_b
    v_N, v_E, v_D = v_eb_n
    ω_eng = pwp.engine.ω
    fuel = fuel.x_avail
    α = aero.α
    β = aero.β
    α_filt = aero.α_filt
    β_filt = aero.β_filt

    f_x, f_y, f_z = dynamics.f_c_c
    EAS = airflow.EAS
    TAS = airflow.TAS
    χ = χ_gnd
    γ = γ_gnd
    c = -v_D

    throttle_out = throttle
    aileron_out = aileron
    elevator_out = elevator
    rudder_out = rudder

    YStateSpace(; p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng, fuel,
            f_x, f_y, f_z, EAS, TAS, α, β, v_N, v_E, v_D, χ, γ, c,
            throttle_out, aileron_out, elevator_out, rudder_out)


end

AircraftBase.ẋ_linear(vehicle::Model{<:C172S.Vehicle{NED}}) = XStateSpace(vehicle.ẋ)
AircraftBase.x_linear(vehicle::Model{<:C172S.Vehicle{NED}}) = XStateSpace(vehicle.x)
AircraftBase.u_linear(vehicle::Model{<:C172S.Vehicle{NED}}) = UStateSpace(vehicle)
AircraftBase.y_linear(vehicle::Model{<:C172S.Vehicle{NED}}) = YStateSpace(vehicle)

function AircraftBase.assign_u!(vehicle::Model{<:C172S.Vehicle{NED}}, u::AbstractVector{Float64})

    @unpack throttle, aileron, elevator, rudder = UStateSpace(u)
    @pack! vehicle.systems.act.u = aileron, elevator, rudder
    @pack! vehicle.systems.pwp.engine.u = throttle

end

function AircraftBase.assign_x!(vehicle::Model{<:C172S.Vehicle{NED}}, x::AbstractVector{Float64})

    @unpack p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng, fuel = XStateSpace(x)

    x_kinematics = vehicle.x.kinematics
    x_dynamics = vehicle.x.dynamics
    x_systems = vehicle.x.systems

    ψ_nb, θ_nb, φ_nb, h_e = ψ, θ, φ, h

    @pack! x_kinematics = ψ_nb, θ_nb, φ_nb, ϕ, λ, h_e
    x_dynamics.ω_eb_b .= p, q, r
    x_dynamics.v_eb_b .= v_x, v_y, v_z
    x_systems.aero .= α_filt, β_filt
    x_systems.pwp.engine.ω = ω_eng
    x_systems.fuel .= fuel

end


################################################################################
############################### Versions #######################################

include(normpath("c172s0.jl")); @reexport using .C172Sv0

end #module