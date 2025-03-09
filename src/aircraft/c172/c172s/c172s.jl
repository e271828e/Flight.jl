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
        ω_cutoff = Piston.RPM2radpersec(3100),
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
################################# Templates ####################################

const Components = C172.Components{typeof(C172S.PowerPlant()), C172.MechanicalActuation}
const Vehicle{K} = AircraftBase.Vehicle{C172S.Components, K} where {K <: AbstractKinematicDescriptor}
const Aircraft{K, A} = AircraftBase.Aircraft{C172S.Vehicle{K}, A} where {K <: AbstractKinematicDescriptor, A <: AbstractAvionics}
const Cessna172S{K, A} = C172S.Aircraft{K, A}

function Vehicle(kinematics = WA())
    AircraftBase.Vehicle(
        C172.Components(C172S.PowerPlant(), C172.MechanicalActuation()),
        kinematics, VehicleDynamics())
end

################################################################################
############################## Cessna172Sv0 #####################################

const Cessna172Sv0{K} = Cessna172S{K, NoAvionics} where { K <: AbstractKinematicDescriptor}

function Cessna172Sv0(kinematics = WA())
    AircraftBase.Aircraft(C172S.Vehicle(kinematics), NoAvionics())
end

############################ Joystick Mappings #################################

#with no Avionics, input assignments must go directly to the actuation system
function Systems.assign_input!(sys::System{<:Cessna172Sv0},
                                data::JoystickData,
                                mapping::IOMapping)
    Systems.assign_input!(sys.vehicle.components.act, data, mapping)
end

############################### Trimming #######################################
################################################################################

#assigns trim state and parameters to vehicle, then updates vehicle
function AircraftBase.assign!(vehicle::System{<:C172S.Vehicle},
                        trim_params::C172.TrimParameters,
                        trim_state::C172.TrimState,
                        atmosphere::System{<:AbstractAtmosphere},
                        terrain::System{<:AbstractTerrain})

    @unpack EAS, β_a, x_fuel, flaps, mixture, payload = trim_params
    @unpack n_eng, α_a, throttle, aileron, elevator, rudder = trim_state
    @unpack act, pwp, aero, fuel, ldg, pld = vehicle.components

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

    aero.x.α_filt = α_a #ensures zero state derivative
    aero.x.β_filt = β_a #ensures zero state derivative
    fuel.x .= Float64(x_fuel)

    kin_init = KinInit(trim_state, trim_params, atmosphere)

    #initialize the vehicle with the setup above. this will call f_ode!
    #internally, no need to do it here
    Systems.init!(vehicle, kin_init, atmosphere, terrain)

    #check essential assumptions about components systems states & derivatives
    @assert !any(SVector{3}(leg.strut.wow for leg in ldg.y))
    @assert pwp.x.engine.ω > pwp.engine.constants.ω_idle
    @assert pwp.x.engine.idle[1] .== 0
    @assert pwp.x.engine.frc[1] .== 0
    @assert abs(aero.ẋ.α_filt) < 1e-10
    @assert abs(aero.ẋ.β_filt) < 1e-10

end

################################################################################
############################### Linearization ##################################

@kwdef struct XLinear <: FieldVector{16, Float64}
    p::Float64 = 0.0; q::Float64 = 0.0; r::Float64 = 0.0; #angular rates (ω_eb_b)
    ψ::Float64 = 0.0; θ::Float64 = 0.0; φ::Float64 = 0.0; #heading, inclination, bank (body/NED)
    v_x::Float64 = 0.0; v_y::Float64 = 0.0; v_z::Float64 = 0.0; #aerodynamic velocity, body axes
    ϕ::Float64 = 0.0; λ::Float64 = 0.0; h::Float64 = 0.0; #latitude, longitude, ellipsoidal altitude
    α_filt::Float64 = 0.0; β_filt::Float64 = 0.0; #filtered airflow angles
    ω_eng::Float64 = 0.0; fuel::Float64 = 0.0; #engine speed, fuel fraction
end

#flaps and mixture are trim parameters and thus omitted from the control vector
@kwdef struct ULinear <: FieldVector{4, Float64}
    throttle::Float64 = 0.0
    aileron::Float64 = 0.0
    elevator::Float64 = 0.0
    rudder::Float64 = 0.0
end

#all states (for full-state feedback), plus other useful stuff, plus control inputs
@kwdef struct YLinear <: FieldVector{33, Float64}
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


function XLinear(x_vehicle::ComponentVector)

    x_kinematics = x_vehicle.kinematics
    x_dynamics = x_vehicle.dynamics
    x_components = x_vehicle.components

    @unpack ψ_nb, θ_nb, φ_nb, ϕ, λ, h_e = x_kinematics
    p, q, r = x_dynamics.ω_eb_b
    v_x, v_y, v_z = x_dynamics.v_eb_b
    α_filt, β_filt = x_components.aero
    ω_eng = x_components.pwp.engine.ω
    fuel = x_components.fuel[1]
    ψ, θ, φ, h = ψ_nb, θ_nb, φ_nb, h_e

    XLinear(; p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng, fuel)

end

function ULinear(vehicle::System{<:C172S.Vehicle{NED}})

    @unpack throttle, aileron, elevator, rudder = vehicle.components.act.u
    ULinear(; throttle, aileron, elevator, rudder)

end

function YLinear(vehicle::System{<:C172S.Vehicle{NED}})

    @unpack throttle, aileron, elevator, rudder = vehicle.components.act.u
    @unpack components, airflow, dynamics, kinematics = vehicle.y
    @unpack pwp, fuel, aero, act = components

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

    YLinear(; p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng, fuel,
            f_x, f_y, f_z, EAS, TAS, α, β, v_N, v_E, v_D, χ, γ, c,
            throttle_out, aileron_out, elevator_out, rudder_out)


end

AircraftBase.ẋ_linear(vehicle::System{<:C172S.Vehicle{NED}}) = XLinear(vehicle.ẋ)
AircraftBase.x_linear(vehicle::System{<:C172S.Vehicle{NED}}) = XLinear(vehicle.x)
AircraftBase.u_linear(vehicle::System{<:C172S.Vehicle{NED}}) = ULinear(vehicle)
AircraftBase.y_linear(vehicle::System{<:C172S.Vehicle{NED}}) = YLinear(vehicle)

function AircraftBase.assign_u!(vehicle::System{<:C172S.Vehicle{NED}}, u::AbstractVector{Float64})

    #The velocity states in the linearized model are meant to be aerodynamic so
    #they can be readily used for flight control design. Since the velocity
    #states in the nonlinear model are Earth-relative, we need to ensure wind
    #velocity is set to zero for linearization.
    # vehicle.atmosphere.u.v_ew_n .= 0
    @unpack throttle, aileron, elevator, rudder = ULinear(u)
    @pack! vehicle.components.act.u = throttle, aileron, elevator, rudder

end

function AircraftBase.assign_x!(vehicle::System{<:C172S.Vehicle{NED}}, x::AbstractVector{Float64})

    @unpack p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng, fuel = XLinear(x)

    x_kinematics = vehicle.x.kinematics
    x_dynamics = vehicle.x.dynamics
    x_components = vehicle.x.components

    ψ_nb, θ_nb, φ_nb, h_e = ψ, θ, φ, h

    @pack! x_kinematics = ψ_nb, θ_nb, φ_nb, ϕ, λ, h_e
    x_dynamics.ω_eb_b .= p, q, r
    x_dynamics.v_eb_b .= v_x, v_y, v_z
    x_components.aero .= α_filt, β_filt
    x_components.pwp.engine.ω = ω_eng
    x_components.fuel .= fuel

end

function Control.Continuous.LinearizedSS(
            vehicle::System{<:C172S.Vehicle{NED}},
            trim_params::C172.TrimParameters = C172.TrimParameters();
            model::Symbol = :full)

    lm = linearize!(vehicle, trim_params)

    if model === :full
        return lm

    elseif model === :lon
        x_labels = [:q, :θ, :v_x, :v_z, :h, :α_filt, :ω_eng]
        u_labels = [:throttle, :elevator]
        y_labels = vcat(x_labels, [:f_x, :f_z, :α, :EAS, :TAS, :γ, :c, :throttle_out, :elevator_out])
        return Control.Continuous.submodel(lm; x = x_labels, u = u_labels, y = y_labels)

    elseif model === :lat
        x_labels = [:p, :r, :ψ, :φ, :v_x, :v_y, :β_filt]
        u_labels = [:aileron, :rudder]
        y_labels = vcat(x_labels, [:f_y, :β, :χ, :aileron_out, :rudder_out])
        return Control.Continuous.submodel(lm; x = x_labels, u = u_labels, y = y_labels)

    else
        error("Valid model keyword values: :full, :lon, :lat")

    end

end



end