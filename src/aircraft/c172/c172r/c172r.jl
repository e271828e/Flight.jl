module C172R

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, Reexport
using ControlSystems, RobustAndOptimalControl
using NLopt
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore
using Flight.FlightLib

using ..AircraftBase
using ..C172


################################################################################
################################ Powerplant ####################################

function PowerPlant()

    prop_data = Propellers.Lookup(Propellers.Blade(), 2)

    propeller = Propeller(prop_data;
        sense = Propellers.CW, d = 2.0, J_xx = 0.3,
        t_bp = FrameTransform(r = [2.055, 0, 0.833]))

    PistonThruster(; propeller)

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

Systems.U(::Actuation) = ActuationU()
Systems.Y(::Actuation) = ActuationY()

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
                pwp::System{<:PistonThruster},
                act::System{<:Actuation})

    @unpack eng_start, eng_stop,
            throttle, mixture, aileron, elevator, rudder,
            aileron_offset, elevator_offset, rudder_offset,
            brake_left, brake_right, flaps = act.y

    pwp.engine.u.start = eng_start
    pwp.engine.u.stop = eng_stop
    pwp.engine.u.throttle = throttle
    pwp.engine.u.mixture = mixture
    ldg.nose.steering.u.input = (rudder_offset + rudder)
    ldg.left.braking.u[] = brake_left
    ldg.right.braking.u[] = brake_right
    aero.u.e = -(elevator_offset + elevator)
    aero.u.a = (aileron_offset + aileron)
    aero.u.r = -(rudder_offset + rudder)
    aero.u.f = flaps

    return nothing
end


function GUI.draw(sys::System{Actuation}, p_open::Ref{Bool} = Ref(true),
                label::String = "Cessna 172R Actuation")

    y = sys.y

    CImGui.Begin(label, p_open)

    CImGui.PushItemWidth(-60)

    CImGui.Text("Engine Start: $(y.eng_start)")
    CImGui.Text("Engine Stop: $(y.eng_stop)")

    # @running_plot("Throttle", y.throttle, 0, 1, 0.0, 120)
    display_bar("Throttle", y.throttle, 0, 1)
    # @running_plot("Aileron", y.aileron, -1, 1, 0.0, 120)
    display_bar("Aileron", y.aileron, -1, 1)
    # @running_plot("Elevator", y.elevator, -1, 1, 0.0, 120)
    display_bar("Elevator", y.elevator, -1, 1)
    # @running_plot("Rudder", y.rudder, -1, 1, 0.0, 120)
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

function GUI.draw!(sys::System{Actuation}, p_open::Ref{Bool} = Ref(true),
                label::String = "Cessna 172R Actuation")

    u = sys.u

    CImGui.Begin(label, p_open)

    CImGui.PushItemWidth(-150)

    mode_button("Engine Start", true, u.eng_start, false; HSV_requested = HSV_green); CImGui.SameLine()
    u.eng_start = CImGui.IsItemActive()
    mode_button("Engine Stop", true, u.eng_stop, false; HSV_requested = HSV_red); CImGui.SameLine()
    u.eng_stop = CImGui.IsItemActive()

    u.throttle = safe_slider("Throttle", u.throttle, "%.6f", true)
    # @running_plot("Throttle", u.throttle, 0, 1, 0.0, 120)
    u.aileron = safe_slider("Aileron", u.aileron, "%.6f", true)
    # @running_plot("Aileron", u.aileron, -1, 1, 0.0, 120)
    u.elevator = safe_slider("Elevator", u.elevator, "%.6f", true)
    # @running_plot("Elevator", u.elevator, -1, 1, 0.0, 120, true)
    u.rudder = safe_slider("Rudder", u.rudder, "%.6f", true)
    # @running_plot("Rudder", u.rudder, -1, 1, 0.0, 120)

    u.aileron_offset = safe_input("Aileron Offset", u.aileron_offset, 0.001, 0.1, "%.6f", true)
    u.elevator_offset = safe_input("Elevator Offset", u.elevator_offset, 0.001, 0.1, "%.6f", true)
    u.rudder_offset = safe_input("Rudder Offset", u.rudder_offset, 0.001, 0.1, "%.6f", true)
    u.flaps = safe_slider("Flaps", u.flaps, "%.6f", true)
    u.mixture = safe_slider("Mixture", u.mixture, "%.6f", true)
    u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f", true)
    u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f", true)

    CImGui.PopItemWidth()

    CImGui.End()

end

################################## Joysticks ###################################

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function Systems.assign_input!(sys::System{<:Actuation},
                           data::XBoxControllerData, ::IOMapping)

    u = sys.u

    u.aileron = get_axis_value(data, :right_stick_x) |> aileron_curve
    u.elevator = get_axis_value(data, :right_stick_y) |> elevator_curve
    u.rudder = get_axis_value(data, :left_stick_x) |> rudder_curve
    u.brake_left = get_axis_value(data, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(data, :right_trigger) |> brake_curve

    u.aileron_offset -= 0.01 * was_released(data, :dpad_left)
    u.aileron_offset += 0.01 * was_released(data, :dpad_right)
    u.elevator_offset += 0.01 * was_released(data, :dpad_down)
    u.elevator_offset -= 0.01 * was_released(data, :dpad_up)

    u.flaps += 0.3333 * was_released(data, :right_bumper)
    u.flaps -= 0.3333 * was_released(data, :left_bumper)

    u.throttle += 0.1 * was_released(data, :button_Y)
    u.throttle -= 0.1 * was_released(data, :button_A)
end

function Systems.assign_input!(sys::System{<:Actuation},
                           data::T16000MData, ::IOMapping)

    u = sys.u

    u.throttle = get_axis_value(data, :throttle)
    u.aileron = get_axis_value(data, :stick_x) |> aileron_curve
    u.elevator = get_axis_value(data, :stick_y) |> elevator_curve
    u.rudder = get_axis_value(data, :stick_z) |> rudder_curve

    u.brake_left = is_pressed(data, :button_1)
    u.brake_right = is_pressed(data, :button_1)

    u.aileron_offset -= 2e-4 * is_pressed(data, :hat_left)
    u.aileron_offset += 2e-4 * is_pressed(data, :hat_right)
    u.elevator_offset += 2e-4 * is_pressed(data, :hat_down)
    u.elevator_offset -= 2e-4 * is_pressed(data, :hat_up)

    u.flaps += 0.3333 * was_released(data, :button_3)
    u.flaps -= 0.3333 * was_released(data, :button_2)

end


################################################################################
################################# Template #####################################

const Components = C172.Components{typeof(PowerPlant()), Actuation}
const Vehicle{K, T} = AircraftBase.Vehicle{Components, K, T} where {K <: AbstractKinematicDescriptor, T <: AbstractTerrain}
const Aircraft{K, T, A} = AircraftBase.Aircraft{Vehicle{K, T}, A} where {K <: AbstractKinematicDescriptor, T <: AbstractTerrain, A <: AbstractAvionics}

function Vehicle(kinematics = LTF(), terrain = HorizontalTerrain())
    AircraftBase.Vehicle(C172.Components(PowerPlant(), Actuation()), kinematics, RigidBodyDynamics(), terrain, LocalAtmosphere())
end

function Aircraft(kinematics = LTF(), terrain = HorizontalTerrain(), avionics = NoAvionics())
    AircraftBase.Aircraft(Vehicle(kinematics, terrain), avionics)
end


############################### Trimming #######################################
################################################################################

#assigns trim state and parameters to vehicle, then updates vehicle
function AircraftBase.assign!(vehicle::System{<:C172R.Vehicle},
                        trim_params::C172.TrimParameters,
                        trim_state::C172.TrimState)

    @unpack EAS, β_a, x_fuel, flaps, mixture, payload = trim_params
    @unpack n_eng, α_a, throttle, aileron, elevator, rudder = trim_state
    @unpack act, pwp, aero, fuel, ldg, pld = vehicle.components

    atm_data = AtmData(vehicle.atmosphere)
    Systems.init!(vehicle, KinInit(trim_state, trim_params, atm_data))

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

    f_ode!(vehicle)

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
    v_N::Float64 = 0.0; v_E::Float64 = 0.0; v_D::Float64 = 0.0; #Ob/ECEF velocity, NED axes
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
    v_x, v_y, v_z = x_dynamics.v_eOb_b
    α_filt, β_filt = x_components.aero
    ω_eng = x_components.pwp.engine.ω
    fuel = x_components.fuel[1]
    ψ, θ, φ, h = ψ_nb, θ_nb, φ_nb, h_e

    XLinear(; p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng, fuel)

end

function ULinear(vehicle::System{<:C172R.Vehicle{NED}})

    @unpack throttle, aileron, elevator, rudder = vehicle.components.act.u
    ULinear(; throttle, aileron, elevator, rudder)

end

function YLinear(vehicle::System{<:C172R.Vehicle{NED}})

    @unpack throttle, aileron, elevator, rudder = vehicle.components.act.u
    @unpack components, air, accelerations, kinematics = vehicle.y
    @unpack pwp, fuel, aero,act = components

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

    f_x, f_y, f_z = accelerations.f_Gb_b
    EAS = air.EAS
    TAS = air.TAS
    α = air.α_b
    β = air.β_b
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

AircraftBase.ẋ_linear(vehicle::System{<:C172R.Vehicle{NED}}) = XLinear(vehicle.ẋ)
AircraftBase.x_linear(vehicle::System{<:C172R.Vehicle{NED}}) = XLinear(vehicle.x)
AircraftBase.u_linear(vehicle::System{<:C172R.Vehicle{NED}}) = ULinear(vehicle)
AircraftBase.y_linear(vehicle::System{<:C172R.Vehicle{NED}}) = YLinear(vehicle)

function AircraftBase.assign_u!(vehicle::System{<:C172R.Vehicle{NED}}, u::AbstractVector{Float64})

    #The velocity states in the linearized model are meant to be aerodynamic so
    #they can be readily used for flight control design. Since the velocity
    #states in the nonlinear model are Earth-relative, we need to ensure wind
    #velocity is set to zero for linearization.
    vehicle.atmosphere.u.v_ew_n .= 0
    @unpack throttle, aileron, elevator, rudder = ULinear(u)
    @pack! vehicle.components.act.u = throttle, aileron, elevator, rudder

end

function AircraftBase.assign_x!(vehicle::System{<:C172R.Vehicle{NED}}, x::AbstractVector{Float64})

    @unpack p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng, fuel = XLinear(x)

    x_kinematics = vehicle.x.kinematics
    x_dynamics = vehicle.x.dynamics
    x_components = vehicle.x.components

    ψ_nb, θ_nb, φ_nb, h_e = ψ, θ, φ, h

    @pack! x_kinematics = ψ_nb, θ_nb, φ_nb, ϕ, λ, h_e
    x_dynamics.ω_eb_b .= p, q, r
    x_dynamics.v_eOb_b .= v_x, v_y, v_z
    x_components.aero .= α_filt, β_filt
    x_components.pwp.engine.ω = ω_eng
    x_components.fuel .= fuel

end

function Control.Continuous.LinearizedSS(
            vehicle::System{<:C172R.Vehicle{NED}},
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


################################################################################
############################### Cessna172R #################################

export Cessna172R

#Cessna172R baseline, C172R.Vehicle with NoAvionics
const Cessna172R{K, T} = C172R.Aircraft{K, T, NoAvionics} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain}

function Cessna172R(kinematics = LTF(), terrain = HorizontalTerrain())
    C172R.Aircraft(kinematics, terrain, NoAvionics())
end

############################ Joystick Mappings #################################

#redirect input assignments directly to the actuation system
function Systems.assign_input!(sys::System{<:Cessna172R},
                                data::JoystickData,
                                mapping::IOMapping)
    Systems.assign_input!(sys.vehicle.components.act, data, mapping)
end


end