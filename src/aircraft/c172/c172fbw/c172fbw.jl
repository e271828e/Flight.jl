module C172FBW

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, Reexport
using ControlSystems, RobustAndOptimalControl

using Flight.FlightCore
using Flight.FlightLib

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

    PistonThruster(; propeller)

end

################################################################################
################################## Actuator2 ###################################

#second order linear actuator model

#with an underdamped actuator, the position state can still transiently exceed
#the intended range due to overshoot. the true actuator position should
#therefore be clamped. in the real world, this behaviour could correspond to a
#clutched output actuator, where the output position saturates beyond a given
#opposing torque (for example, if the surface's mechanical limits are hit)

#saturate on command, not on position, which only tends asymptotically to cmd!

struct Actuator2{R} <: SystemDefinition #second order linear actuator model
    ω_n::Float64 #natural frequency (default: 10 Hz)
    ζ::Float64 #damping ratio (default: underdamped with minimal resonance)
    function Actuator2(; ω_n::Real = 5*2π, ζ::Real = 0.6,
                        range::Tuple{Real, Real} = (-1.0, 1.0))
        new{Ranged{Float64, range[1], range[2]}}(ω_n, ζ)
    end
end

@kwdef struct Actuator2Y{R}
    cmd::R = R(0.0)
    pos::R = R(0.0)
    vel::Float64 = 0.0
    sat::Int64 = 0
end

Systems.X(::Actuator2) = ComponentVector(v = 0.0, p = 0.0)
Systems.U(::Actuator2{R}) where {R} = Ref(R(0.0))
Systems.Y(::Actuator2{R}) where {R} = Actuator2Y{R}()

function Systems.f_ode!(sys::System{Actuator2{R}}) where {R}

    @unpack ẋ, x, u, constants = sys
    @unpack ω_n, ζ = constants

    cmd = u[]
    pos = R(x.p)
    vel = x.v
    sat = saturation(cmd)

    ẋ.v = ω_n^2 * (Float64(cmd) - x.p) - 2ζ*ω_n*x.v
    ẋ.p = x.v

    sys.y = Actuator2Y(; cmd, pos, vel, sat)

end


################################################################################
#################################### Actuation #################################

#Fly-by-wire actuation system. Throttle, steering and aerodynamic surfaces are
#controlled via actuators, the rest of are direct feedthrough

@kwdef struct Actuation <: C172.Actuation
    throttle::Actuator2 = Actuator2(range = (0.0, 1.0))
    aileron::Actuator2 = Actuator2(range = (-1.0, 1.0))
    elevator::Actuator2 = Actuator2(range = (-1.0, 1.0))
    rudder::Actuator2 = Actuator2(range = (-1.0, 1.0))
    flaps::Actuator2 = Actuator2(range = (0.0, 1.0))
    steering::Actuator2 = Actuator2(range = (-1.0, 1.0))
end

function C172.assign!(aero::System{<:C172.Aero},
                    ldg::System{<:C172.Ldg},
                    pwp::System{<:PistonThruster},
                    act::System{<:Actuation})

    @unpack throttle, aileron, elevator, rudder, flaps, steering = act.y

    aero.u.e = -elevator.pos
    aero.u.a = aileron.pos
    aero.u.r = -rudder.pos
    aero.u.f = flaps.pos
    pwp.engine.u.throttle = throttle.pos
    ldg.nose.steering.u.input = steering.pos

end

################################### GUI ########################################

function GUI.draw(sys::System{Actuation}, p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172 Fly-By-Wire Actuation")

    CImGui.Begin(label, p_open)
    CImGui.PushItemWidth(-60)

    #because @cstatic allocates storage using global variables, functions using
    #it are not reentrant; calls to such functions will share the same storage.

    labels = uppercasefirst.(string.(keys(sys.subsystems)))
    commands = map(ss->ss.y.cmd, values(sys.subsystems))
    positions = map(ss->ss.y.pos, values(sys.subsystems))

    @cstatic(
        thr_buffer = fill(Cfloat(0),90), thr_offset = Ref(Cint(0)),
        ail_buffer = fill(Cfloat(0),90), ail_offset = Ref(Cint(0)),
        ele_buffer = fill(Cfloat(0),90), ele_offset = Ref(Cint(0)),
        rud_buffer = fill(Cfloat(0),90), rud_offset = Ref(Cint(0)),
        flp_buffer = fill(Cfloat(0),90), flp_offset = Ref(Cint(0)),
        str_buffer = fill(Cfloat(0),90), str_offset = Ref(Cint(0)),

        begin

            buffers = (thr_buffer, ail_buffer, ele_buffer, rud_buffer, flp_buffer, str_buffer)
            offsets = (thr_offset, ail_offset, ele_offset, rud_offset, flp_offset, str_offset)

            foreach(labels, commands, positions, buffers, offsets) do label, command, position, buffer, offset

                if CImGui.CollapsingHeader(label)
                    CImGui.Text("$label Command"); CImGui.SameLine(200); display_bar("", command)
                    CImGui.Text("$label Position"); CImGui.SameLine(200); display_bar("", position)
                    buffer[offset[]+1] = Cfloat(position)
                    offset[] = (offset[]+1) % length(buffer)
                    CImGui.PlotLines("$label Position", buffer, length(buffer), offset[], "$label Position",
                                    Cfloat(typemin(position)), Cfloat(typemax(position)),
                                    (Cint(0), Cint(120)))
                end
            end

        end)

    CImGui.PopItemWidth()
    CImGui.End()

end

function GUI.draw!(sys::System{Actuation}, p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172 Fly-By-Wire Actuation")

    CImGui.Begin(label, p_open)
    CImGui.PushItemWidth(-60)

    #because @cstatic allocates storage using global variables, functions using
    #it are not reentrant; calls to such functions will share the same storage.

    labels = uppercasefirst.(string.(keys(sys.subsystems)))
    inputs = map(ss->ss.u, values(sys.subsystems))
    commands = map(ss->ss.y.cmd, values(sys.subsystems))
    positions = map(ss->ss.y.pos, values(sys.subsystems))

    @cstatic(
        thr_buffer = fill(Cfloat(0),90), thr_offset = Ref(Cint(0)),
        ail_buffer = fill(Cfloat(0),90), ail_offset = Ref(Cint(0)),
        ele_buffer = fill(Cfloat(0),90), ele_offset = Ref(Cint(0)),
        rud_buffer = fill(Cfloat(0),90), rud_offset = Ref(Cint(0)),
        flp_buffer = fill(Cfloat(0),90), flp_offset = Ref(Cint(0)),
        str_buffer = fill(Cfloat(0),90), str_offset = Ref(Cint(0)),

        begin

            buffers = (thr_buffer, ail_buffer, ele_buffer, rud_buffer, flp_buffer, str_buffer)
            offsets = (thr_offset, ail_offset, ele_offset, rud_offset, flp_offset, str_offset)

            foreach(labels, inputs, commands, positions, buffers, offsets) do label, input, command, position, buffer, offset

                if CImGui.CollapsingHeader(label)
                    input[] = safe_slider("$label Command", input[], "%.6f")
                    CImGui.Text("$label Command"); CImGui.SameLine(200); display_bar("", command)
                    CImGui.Text("$label Position"); CImGui.SameLine(200); display_bar("", position)
                    buffer[offset[]+1] = Cfloat(position)
                    offset[] = (offset[]+1) % length(buffer)
                    CImGui.PlotLines("$label Position", buffer, length(buffer), offset[], "$label Position",
                                    Cfloat(typemin(position)), Cfloat(typemax(position)),
                                    (Cint(0), Cint(120)))
                end
            end

        end)

    CImGui.PopItemWidth()
    CImGui.End()

end


################################## Joysticks ###################################

pitch_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
roll_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
yaw_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
# brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function Systems.assign_input!(sys::System{<:Actuation},
                                data::XBoxControllerData, ::IOMapping)

    @unpack throttle, aileron, elevator, rudder, steering, flaps = sys.subsystems

    aileron.u[] = get_axis_value(data, :right_stick_x) |> roll_curve
    elevator.u[] = get_axis_value(data, :right_stick_y) |> pitch_curve
    rudder.u[] = get_axis_value(data, :left_stick_x) |> yaw_curve
    # u.brake_left = get_axis_value(data, :left_trigger) |> brake_curve
    # u.brake_right = get_axis_value(data, :right_trigger) |> brake_curve
    steering.u[] = get_axis_value(data, :left_stick_x) |> yaw_curve

    flaps.u[] += 0.3333 * was_released(data, :right_bumper)
    flaps.u[] -= 0.3333 * was_released(data, :left_bumper)

    throttle.u[] += 0.1 * was_released(data, :button_Y)
    throttle.u[] -= 0.1 * was_released(data, :button_A)
end

function Systems.assign_input!(sys::System{<:Actuation},
                                data::T16000MData, ::IOMapping)

    @unpack throttle, aileron, elevator, rudder, steering, flaps = sys.subsystems

    throttle.u[] = get_axis_value(data, :throttle)
    aileron.u[] = get_axis_value(data, :stick_x) |> roll_curve
    elevator.u[] = get_axis_value(data, :stick_y) |> pitch_curve
    rudder.u[] = get_axis_value(data, :stick_z) |> yaw_curve
    steering.u[] = get_axis_value(data, :stick_z) |> yaw_curve

    # u.brake_left = is_pressed(data, :button_1)
    # u.brake_right = is_pressed(data, :button_1)

    flaps.u[] += 0.3333 * was_released(data, :button_3)
    flaps.u[] -= 0.3333 * was_released(data, :button_2)

end


################################################################################
################################# Templates ####################################

const Components = C172.Components{typeof(PowerPlant()), C172FBW.Actuation}
const Vehicle{K, T} = AircraftBase.Vehicle{C172FBW.Components, K, T} where {K <: AbstractKinematicDescriptor, T <: AbstractTerrain}
const Aircraft{K, T, A} = AircraftBase.Aircraft{C172FBW.Vehicle{K, T}, A} where {K <: AbstractKinematicDescriptor, T <: AbstractTerrain, A <: AbstractAvionics}

function Vehicle(kinematics = WA(), terrain = HorizontalTerrain())
    AircraftBase.Vehicle(
        C172.Components(PowerPlant(), Actuation()),
        kinematics,
        VehicleDynamics(),
        terrain,
        LocalAtmosphere())
end

function Aircraft(kinematics = WA(), terrain = HorizontalTerrain(), avionics = NoAvionics())
    AircraftBase.Aircraft(Vehicle(kinematics, terrain), avionics)
end

############################### Trimming #######################################
################################################################################

#assigns trim state and parameters to vehicle, then updates vehicle
function AircraftBase.assign!(vehicle::System{<:C172FBW.Vehicle},
                        trim_params::C172.TrimParameters,
                        trim_state::C172.TrimState)

    @unpack EAS, β_a, x_fuel, flaps, mixture, payload = trim_params
    @unpack n_eng, α_a, throttle, aileron, elevator, rudder = trim_state
    @unpack act, pwp, aero, fuel, ldg, pld = vehicle.components

    atm_data = AtmData(vehicle.atmosphere)
    Systems.init!(vehicle, KinInit(trim_state, trim_params, atm_data))

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

    f_ode!(vehicle)

    #check essential assumptions about components systems states & derivatives
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
    v_N::Float64 = 0.0; v_E::Float64 = 0.0; v_D::Float64 = 0.0; #b/ECEF velocity, NED axes
    χ::Float64 = 0.0; γ::Float64 = 0.0; climb_rate::Float64 = 0.0; #track and flight path angles, climb rate
    throttle_cmd::Float64 = 0.0; aileron_cmd::Float64 = 0.0; #actuator commands
    elevator_cmd::Float64 = 0.0; rudder_cmd::Float64 = 0.0; #actuator commands
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
    thr_v = x_components.act.throttle.v
    thr_p = x_components.act.throttle.p
    ail_v = x_components.act.aileron.v
    ail_p = x_components.act.aileron.p
    ele_v = x_components.act.elevator.v
    ele_p = x_components.act.elevator.p
    rud_v = x_components.act.rudder.v
    rud_p = x_components.act.rudder.p

    ψ, θ, φ, h = ψ_nb, θ_nb, φ_nb, h_e

    XLinear(;  p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt,
        ω_eng, fuel, thr_v, thr_p, ail_v, ail_p, ele_v, ele_p, rud_v, rud_p)

end

function ULinear(vehicle::System{<:C172FBW.Vehicle{NED}})

    @unpack throttle, aileron, elevator, rudder = vehicle.components.act.subsystems
    throttle_cmd = throttle.u[]
    aileron_cmd = aileron.u[]
    elevator_cmd = elevator.u[]
    rudder_cmd = rudder.u[]
    ULinear(; throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd)

end

function YLinear(vehicle::System{<:C172FBW.Vehicle{NED}})

    @unpack components, air, dynamics, kinematics = vehicle.y
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

    thr_v = act.throttle.vel
    thr_p = act.throttle.pos
    ail_v = act.aileron.vel
    ail_p = act.aileron.pos
    ele_v = act.elevator.vel
    ele_p = act.elevator.pos
    rud_v = act.rudder.vel
    rud_p = act.rudder.pos

    f_x, f_y, f_z = dynamics.f_c_c
    EAS = air.EAS
    TAS = air.TAS
    χ = χ_gnd
    γ = γ_gnd
    climb_rate = -v_D

    @unpack throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd = ULinear(vehicle)

    YLinear(; p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt,
            ω_eng, fuel, thr_v, thr_p, ail_v, ail_p, ele_v, ele_p, rud_v, rud_p,
            f_x, f_y, f_z, EAS, TAS, α, β, v_N, v_E, v_D, χ, γ, climb_rate,
            throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd)

end

AircraftBase.ẋ_linear(vehicle::System{<:C172FBW.Vehicle{NED}}) = XLinear(vehicle.ẋ)
AircraftBase.x_linear(vehicle::System{<:C172FBW.Vehicle{NED}}) = XLinear(vehicle.x)
AircraftBase.u_linear(vehicle::System{<:C172FBW.Vehicle{NED}}) = ULinear(vehicle)
AircraftBase.y_linear(vehicle::System{<:C172FBW.Vehicle{NED}}) = YLinear(vehicle)

function AircraftBase.assign_u!(vehicle::System{<:C172FBW.Vehicle{NED}}, u::AbstractVector{Float64})

    #The velocity states in the linearized model are meant to be aerodynamic so
    #they can be readily used for flight control design. Since the velocity
    #states in the nonlinear model are Earth-relative, we need to ensure wind
    #velocity is set to zero for linearization.
    @unpack throttle, aileron, elevator, rudder = vehicle.components.act.subsystems
    @unpack throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd = ULinear(u)
    throttle.u[] = throttle_cmd
    aileron.u[] = aileron_cmd
    elevator.u[] = elevator_cmd
    rudder.u[] = rudder_cmd

    vehicle.atmosphere.u.v_ew_n .= 0

end

function AircraftBase.assign_x!(vehicle::System{<:C172FBW.Vehicle{NED}}, x::AbstractVector{Float64})

    @unpack p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng,
            fuel, thr_v, thr_p, ail_v, ail_p, ele_v, ele_p, rud_v, rud_p = XLinear(x)

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
    x_components.act.throttle.v = thr_v
    x_components.act.throttle.p = thr_p
    x_components.act.aileron.v = ail_v
    x_components.act.aileron.p = ail_p
    x_components.act.elevator.v = ele_v
    x_components.act.elevator.p = ele_p
    x_components.act.rudder.v = rud_v
    x_components.act.rudder.p = rud_p

end

function Control.Continuous.LinearizedSS(
            vehicle::System{<:C172FBW.Vehicle{NED}},
            trim_params::C172.TrimParameters = C172.TrimParameters();
            model::Symbol = :full)

    lm = linearize!(vehicle, trim_params)

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
############################### Cessna172FBW ###############################

export Cessna172FBW

#C172FBW.Vehicle with NoAvionics
const Cessna172FBW{K, T} = C172FBW.Aircraft{K, T, NoAvionics} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain}

function Cessna172FBW(kinematics = WA(), terrain = HorizontalTerrain())
    C172FBW.Aircraft(kinematics, terrain, NoAvionics())
end

############################ Joystick Mappings #################################

#map input assignments directly to the actuation system
function Systems.assign_input!(sys::System{<:Cessna172FBW},
                                data::JoystickData,
                                mapping::IOMapping)

    Systems.assign_input!(sys.vehicle.components.act, data, mapping)
end

################################################################################
################################## Variants ####################################

include(normpath("control/c172fbw_ctl.jl"))

include(normpath("c172fbw_v1.jl")); @reexport using .C172FBWv1

end