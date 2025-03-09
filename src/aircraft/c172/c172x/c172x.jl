module C172X

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, Reexport
using ControlSystems, RobustAndOptimalControl

using Flight.FlightCore
using Flight.FlightLib

using ..C172

export Cessna172X, Cessna172Xv0

################################################################################
################################## Actuator1 ###################################

#first order linear actuator model

#τ has units of rad/s, not Hz: F(s) = 1/(1+sτ); F(jω) = 1/(1+jωτ).

struct Actuator1{R <: Ranged} <: SystemDefinition #second order linear actuator model
    τ::Float64 #time constant (default: 0.05s)
    function Actuator1(; τ::Real = 1/20, range::Tuple{Real, Real} = (-1.0, 1.0))
        new{Ranged{Float64, range[1], range[2]}}(τ)
    end
end

@kwdef struct Actuator1Y{R}
    cmd::R = R(0.0)
    pos::R = R(0.0)
    sat::Int64 = 0
end

Systems.X(::Actuator1) = ComponentVector(p = 0.0)
Systems.U(::Actuator1{R}) where {R} = Ref(R(0.0))
Systems.Y(::Actuator1{R}) where {R} = Actuator1Y{R}()

function Systems.f_ode!(sys::System{Actuator1{R}}) where {R}

    @unpack ẋ, x, u, constants = sys
    @unpack τ = constants

    cmd = u[]
    pos = R(x.p)
    sat = saturation(cmd)

    ẋ.p =  1/τ * (Float64(cmd) - x.p)

    sys.y = Actuator1Y(; cmd, pos, sat)

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
############################# FlyByWireActuation ###############################

@kwdef struct FlyByWireActuation <: C172.AbstractActuation
    throttle::Actuator1 = Actuator1(range = (0.0, 1.0))
    mixture::Actuator1 = Actuator1(range = (0.0, 1.0))
    aileron::Actuator1 = Actuator1(range = (-1.0, 1.0))
    elevator::Actuator1 = Actuator1(range = (-1.0, 1.0))
    rudder::Actuator1 = Actuator1(range = (-1.0, 1.0))
    flaps::Actuator1 = Actuator1(range = (0.0, 1.0))
    steering::Actuator1 = Actuator1(range = (-1.0, 1.0))
    brake_left::Actuator1 = Actuator1(range = (0.0, 1.0))
    brake_right::Actuator1 = Actuator1(range = (0.0, 1.0))
end

function C172.assign!(aero::System{<:C172.Aero},
                    ldg::System{<:C172.Ldg},
                    pwp::System{<:PistonThruster},
                    act::System{<:FlyByWireActuation})

    @unpack throttle, mixture, aileron, elevator, rudder, flaps, steering,
            brake_left, brake_right = act.y

    aero.u.e = -elevator.pos
    aero.u.a = aileron.pos
    aero.u.r = -rudder.pos
    aero.u.f = flaps.pos
    pwp.engine.u.throttle = throttle.pos
    pwp.engine.u.mixture = mixture.pos
    ldg.nose.steering.u.input = steering.pos
    ldg.left.braking.u[] = brake_left.pos
    ldg.right.braking.u[] = brake_right.pos

end


################################### GUI ########################################

function GUI.draw(sys::System{FlyByWireActuation}, p_open::Ref{Bool} = Ref(true),
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
        mix_buffer = fill(Cfloat(0),90), mix_offset = Ref(Cint(0)),
        ail_buffer = fill(Cfloat(0),90), ail_offset = Ref(Cint(0)),
        ele_buffer = fill(Cfloat(0),90), ele_offset = Ref(Cint(0)),
        rud_buffer = fill(Cfloat(0),90), rud_offset = Ref(Cint(0)),
        flp_buffer = fill(Cfloat(0),90), flp_offset = Ref(Cint(0)),
        str_buffer = fill(Cfloat(0),90), str_offset = Ref(Cint(0)),
        bkl_buffer = fill(Cfloat(0),90), bkl_offset = Ref(Cint(0)),
        bkr_buffer = fill(Cfloat(0),90), bkr_offset = Ref(Cint(0)),

        begin

            buffers = (thr_buffer, mix_buffer, ail_buffer, ele_buffer, rud_buffer, flp_buffer, str_buffer, bkl_buffer, bkr_buffer)
            offsets = (thr_offset, mix_offset, ail_offset, ele_offset, rud_offset, flp_offset, str_offset, bkl_offset, bkr_offset)

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

function GUI.draw!(sys::System{FlyByWireActuation}, p_open::Ref{Bool} = Ref(true),
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
        mix_buffer = fill(Cfloat(0),90), mix_offset = Ref(Cint(0)),
        ail_buffer = fill(Cfloat(0),90), ail_offset = Ref(Cint(0)),
        ele_buffer = fill(Cfloat(0),90), ele_offset = Ref(Cint(0)),
        rud_buffer = fill(Cfloat(0),90), rud_offset = Ref(Cint(0)),
        flp_buffer = fill(Cfloat(0),90), flp_offset = Ref(Cint(0)),
        str_buffer = fill(Cfloat(0),90), str_offset = Ref(Cint(0)),
        bkl_buffer = fill(Cfloat(0),90), bkl_offset = Ref(Cint(0)),
        bkr_buffer = fill(Cfloat(0),90), bkr_offset = Ref(Cint(0)),

        begin

            buffers = (thr_buffer, mix_buffer, ail_buffer, ele_buffer, rud_buffer, flp_buffer, str_buffer, bkl_buffer, bkr_buffer)
            offsets = (thr_offset, mix_offset, ail_offset, ele_offset, rud_offset, flp_offset, str_offset, bkl_offset, bkr_offset)

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
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function Systems.assign_input!(sys::System{<:FlyByWireActuation},
                           joystick::XBoxController,
                           ::IOMapping)

    #mixture and brakes will not be assigned
    @unpack throttle, mixture, aileron, elevator, rudder, steering, flaps,
            brake_left, brake_right = sys.subsystems

    aileron.u[] = get_axis_value(joystick, :right_stick_x) |> roll_curve
    elevator.u[] = get_axis_value(joystick, :right_stick_y) |> pitch_curve
    rudder.u[] = get_axis_value(joystick, :left_stick_x) |> yaw_curve
    steering.u[] = get_axis_value(joystick, :left_stick_x) |> yaw_curve
    # brake_left.u[] = get_axis_value(joystick, :left_trigger) |> brake_curve
    # brake_right.u[] = get_axis_value(joystick, :right_trigger) |> brake_curve

    flaps.u[] += 0.3333 * was_released(joystick, :right_bumper)
    flaps.u[] -= 0.3333 * was_released(joystick, :left_bumper)

    throttle.u[] += 0.1 * was_released(joystick, :button_Y)
    throttle.u[] -= 0.1 * was_released(joystick, :button_A)
end

function Systems.assign_input!(sys::System{<:FlyByWireActuation},
                           joystick::T16000M,
                           ::IOMapping)

    #mixture and brakes will not be assigned
    @unpack throttle, mixture, aileron, elevator, rudder, steering, flaps,
            brake_left, brake_right = sys.subsystems

    throttle.u[] = get_axis_value(joystick, :throttle)
    aileron.u[] = get_axis_value(joystick, :stick_x) |> roll_curve
    elevator.u[] = get_axis_value(joystick, :stick_y) |> pitch_curve
    rudder.u[] = get_axis_value(joystick, :stick_z) |> yaw_curve
    steering.u[] = get_axis_value(joystick, :stick_z) |> yaw_curve
    # brake_left.u[] = is_pressed(joystick, :button_1)
    # brake_right.u[] = is_pressed(joystick, :button_1)

    flaps.u[] += 0.3333 * was_released(joystick, :button_3)
    flaps.u[] -= 0.3333 * was_released(joystick, :button_2)

end


################################################################################
################################# Templates ####################################

#reuse C172S power plant, replace actuation system
const Components = C172.Components{typeof(C172S.PowerPlant()), FlyByWireActuation}
const Vehicle{K} = AircraftBase.Vehicle{C172X.Components, K} where {K <: AbstractKinematicDescriptor}
const Aircraft{K, A} = AircraftBase.Aircraft{C172X.Vehicle{K}, A} where {K <: AbstractKinematicDescriptor, A <: AbstractAvionics}
const Cessna172X{K, A} = C172X.Aircraft{K, A}

function Vehicle(kinematics = WA())
    AircraftBase.Vehicle(
        C172.Components(C172S.PowerPlant(), FlyByWireActuation()),
        kinematics, VehicleDynamics())
end

################################################################################
################################# Cessna172Xv0 ################################

const Cessna172Xv0{K} = Cessna172X{K, NoAvionics} where { K <: AbstractKinematicDescriptor}

function Cessna172Xv0(kinematics = WA())
    AircraftBase.Aircraft(Vehicle(kinematics), NoAvionics())
end

############################ Joystick Mappings #################################

#map input assignments directly to the actuation system
function Systems.assign_input!(sys::System{<:Cessna172Xv0},
                                joystick::JoystickData,
                                mapping::IOMapping)
    Systems.assign_input!(sys.vehicle.components.act, joystick, mapping)
end


############################### Trimming #######################################
################################################################################

#assigns trim state and parameters to vehicle, then updates vehicle
function AircraftBase.assign!(vehicle::System{<:C172X.Vehicle},
                        trim_params::C172.TrimParameters,
                        trim_state::C172.TrimState)

    @unpack EAS, β_a, x_fuel, flaps, mixture, payload = trim_params
    @unpack n_eng, α_a, throttle, aileron, elevator, rudder = trim_state
    @unpack act, pwp, aero, fuel, ldg, pld = vehicle.components

    atm_data = AtmosphericData(vehicle.atmosphere)
    Systems.init!(vehicle, KinInit(trim_state, trim_params, atm_data))

    act.throttle.u[] = throttle
    act.mixture.u[] = mixture
    act.aileron.u[] = aileron
    act.elevator.u[] = elevator
    act.rudder.u[] = rudder
    act.flaps.u[] = flaps

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

    #actuator states: in steady state every actuator's position state must be
    #equal to the actuator command.
    act.x.throttle.p = throttle
    act.x.mixture.p = mixture
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

    @assert all(SVector{9,Float64}(act.ẋ) .== 0)

end



################################################################################
############################### Linearization ##################################

@kwdef struct XLinear <: FieldVector{20, Float64}
    p::Float64 = 0.0; q::Float64 = 0.0; r::Float64 = 0.0; #angular rates (ω_eb_b)
    ψ::Float64 = 0.0; θ::Float64 = 0.0; φ::Float64 = 0.0; #heading, inclination, bank (body/NED)
    v_x::Float64 = 0.0; v_y::Float64 = 0.0; v_z::Float64 = 0.0; #aerodynamic velocity, body axes
    ϕ::Float64 = 0.0; λ::Float64 = 0.0; h::Float64 = 0.0; #latitude, longitude, ellipsoidal altitude
    α_filt::Float64 = 0.0; β_filt::Float64 = 0.0; #filtered airflow angles
    ω_eng::Float64 = 0.0; fuel::Float64 = 0.0; #engine speed, fuel fraction
    thr_p::Float64 = 0.0; #throttle actuator position
    ail_p::Float64 = 0.0; #aileron actuator position
    ele_p::Float64 = 0.0; #elevator actuator position
    rud_p::Float64 = 0.0 #rudder actuator position
end

#flaps and mixture are trim parameters and thus omitted from the control vector
@kwdef struct ULinear <: FieldVector{4, Float64}
    throttle_cmd::Float64 = 0.0
    aileron_cmd::Float64 = 0.0
    elevator_cmd::Float64 = 0.0
    rudder_cmd::Float64 = 0.0
end

#all states (for full-state feedback), plus other useful stuff, plus control inputs
@kwdef struct YLinear <: FieldVector{37, Float64}
    p::Float64 = 0.0; q::Float64 = 0.0; r::Float64 = 0.0; #angular rates (ω_eb_b)
    ψ::Float64 = 0.0; θ::Float64 = 0.0; φ::Float64 = 0.0; #heading, inclination, bank (body/NED)
    v_x::Float64 = 0.0; v_y::Float64 = 0.0; v_z::Float64 = 0.0; #aerodynamic velocity, body axes
    ϕ::Float64 = 0.0; λ::Float64 = 0.0; h::Float64 = 0.0; #latitude, longitude, ellipsoidal altitude
    α_filt::Float64 = 0.0; β_filt::Float64 = 0.0; #filtered airflow angles
    ω_eng::Float64 = 0.0; fuel::Float64 = 0.0; #engine speed, available fuel fraction
    thr_p::Float64 = 0.0; #throttle actuator position
    ail_p::Float64 = 0.0; #aileron actuator position
    ele_p::Float64 = 0.0; #elevator actuator position
    rud_p::Float64 = 0.0; #rudder actuator position
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
    thr_p = x_components.act.throttle.p
    ail_p = x_components.act.aileron.p
    ele_p = x_components.act.elevator.p
    rud_p = x_components.act.rudder.p

    ψ, θ, φ, h = ψ_nb, θ_nb, φ_nb, h_e

    XLinear(;  p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt,
        ω_eng, fuel, thr_p, ail_p, ele_p, rud_p)

end

function ULinear(vehicle::System{<:C172X.Vehicle{NED}})

    @unpack throttle, aileron, elevator, rudder = vehicle.components.act.subsystems
    throttle_cmd = throttle.u[]
    aileron_cmd = aileron.u[]
    elevator_cmd = elevator.u[]
    rudder_cmd = rudder.u[]
    ULinear(; throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd)

end

function YLinear(vehicle::System{<:C172X.Vehicle{NED}})

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

    thr_p = act.throttle.pos
    ail_p = act.aileron.pos
    ele_p = act.elevator.pos
    rud_p = act.rudder.pos

    f_x, f_y, f_z = dynamics.f_c_c
    EAS = air.EAS
    TAS = air.TAS
    χ = χ_gnd
    γ = γ_gnd
    climb_rate = -v_D

    @unpack throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd = ULinear(vehicle)

    YLinear(; p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt,
            ω_eng, fuel, thr_p, ail_p, ele_p, rud_p,
            f_x, f_y, f_z, EAS, TAS, α, β, v_N, v_E, v_D, χ, γ, climb_rate,
            throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd)

end

AircraftBase.ẋ_linear(vehicle::System{<:C172X.Vehicle{NED}}) = XLinear(vehicle.ẋ)
AircraftBase.x_linear(vehicle::System{<:C172X.Vehicle{NED}}) = XLinear(vehicle.x)
AircraftBase.u_linear(vehicle::System{<:C172X.Vehicle{NED}}) = ULinear(vehicle)
AircraftBase.y_linear(vehicle::System{<:C172X.Vehicle{NED}}) = YLinear(vehicle)

function AircraftBase.assign_u!(vehicle::System{<:C172X.Vehicle{NED}}, u::AbstractVector{Float64})

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

function AircraftBase.assign_x!(vehicle::System{<:C172X.Vehicle{NED}}, x::AbstractVector{Float64})

    @unpack p, q, r, ψ, θ, φ, v_x, v_y, v_z, ϕ, λ, h, α_filt, β_filt, ω_eng,
            fuel, thr_p, ail_p, ele_p, rud_p = XLinear(x)

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
    x_components.act.throttle.p = thr_p
    x_components.act.aileron.p = ail_p
    x_components.act.elevator.p = ele_p
    x_components.act.rudder.p = rud_p

end

function Control.Continuous.LinearizedSS(
            vehicle::System{<:C172X.Vehicle{NED}},
            trim_params::C172.TrimParameters = C172.TrimParameters();
            model::Symbol = :full)

    lm = linearize!(vehicle, trim_params)

    if model === :full
        return lm

    #preserve the ordering of the complete linearized state and output vectors
    elseif model === :lon
        x_labels = [:q, :θ, :v_x, :v_z, :h, :α_filt, :ω_eng, :thr_p, :ele_p]
        u_labels = [:throttle_cmd, :elevator_cmd]
        y_labels = vcat(x_labels, [:f_x, :f_z, :α, :EAS, :TAS, :γ, :climb_rate, :throttle_cmd, :elevator_cmd])
        return Control.Continuous.submodel(lm; x = x_labels, u = u_labels, y = y_labels)

    elseif model === :lat
        x_labels = [:p, :r, :ψ, :φ, :v_x, :v_y, :β_filt, :ail_p, :rud_p]
        u_labels = [:aileron_cmd, :rudder_cmd]
        y_labels = vcat(x_labels, [:f_y, :β, :χ, :aileron_cmd, :rudder_cmd])
        return Control.Continuous.submodel(lm; x = x_labels, u = u_labels, y = y_labels)

    else
        error("Valid model keyword values: :full, :lon, :lat")

    end

end


################################################################################
################################## Variants ####################################

include(normpath("control/c172x_ctl.jl"))
include(normpath("navigation/c172x_nav.jl"))

include(normpath("c172x1.jl")); @reexport using .C172Xv1
include(normpath("c172x2.jl")); @reexport using .C172Xv2

end