module C172Rv2

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI
using Flight.FlightCore.IODevices
using Flight.FlightCore.Joysticks
using Flight.FlightCore.Utils: Ranged

using Flight.FlightPhysics.Kinematics
using Flight.FlightPhysics.RigidBody
using Flight.FlightPhysics.Environment

using Flight.FlightAircraft.Control
using Flight.FlightAircraft.Aircraft
using Flight.FlightAircraft.World
using Flight.FlightAircraft.Control: PIDDiscreteY

using ..Airframe

export Cessna172Rv2

################################################################################
############################### Avionics #################################

Base.@kwdef struct PitchRateControl <: Component
    c1::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 0, k_i = 1, k_d = 0) #pure integrator
    c2::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 10, k_i = 20, k_d = 0.5, τ_d = 0.05, β_p = 1, β_d = 1) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
Base.@kwdef mutable struct PitchRateControlU
    q_cmd::Float64 = 0.0
    q_fbk::Float64 = 0.0
    reset::Bool = false
end

Base.@kwdef struct PitchRateControlY
    q_cmd::Float64 = 0.0
    q_fbk::Float64 = 0.0
    e_cmd::Float64 = 0.0 #elevator command
    reset::Bool = false
    c1::PIDDiscreteY{1} = PIDDiscreteY{1}()
    c2::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::PitchRateControl) = PitchRateControlU()
Systems.init(::SystemY, ::PitchRateControl) = PitchRateControlY()

function Systems.init!(sys::System{PitchRateControl})
    @unpack c1, c2 = sys.subsystems
    c1.u.bound_lo .= -Inf
    c1.u.bound_hi .= Inf
    c1.u.anti_windup .= true

    c2.u.bound_lo .= -1 #lower bound for MechanicalActuation's normalized elevator input
    c2.u.bound_hi .= 1 #upper bound for MechanicalActuation's normalized elevator input
    c2.u.anti_windup .= true

    println("Try setting β_d = 0 in c2")
end

function Systems.f_disc!(sys::System{PitchRateControl}, Δt::Real)
    @unpack q_cmd, q_fbk, reset = sys.u
    @unpack c1, c2 = sys.subsystems

    c1.u.setpoint .= q_cmd
    c1.u.feedback .= q_fbk
    c1.u.reset .= reset
    c1.u.sat_ext .= c2.y.sat_out
    f_disc!(c1, Δt)

    c2.u.setpoint .= c1.y.out #connected to c1's output
    c2.u.feedback .= 0.0 #no feedback, just feedforward path
    c2.u.reset .= reset
    c2.u.sat_ext .= 0 #only output saturation required
    f_disc!(c2, Δt)

    e_cmd = c2.y.out[1]

    sys.y = PitchRateControlY(; q_cmd, q_fbk, e_cmd, reset, c1 = c1.y, c2 = c2.y)

end

struct RollRateControl <: Component end
struct RollRateControlY end

#rationale for beta control in the inner CAS is that the user is likely to
#desire automatic turn combination in conjunction with roll rate and pitch rate
#augmentation, while yaw rate augmentation is not useful by itself
struct SideslipControl <: Component end
struct SideslipControlY end

################################ Avionics ######################################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@enum CASState begin
    CAS_disabled = 0
    CAS_standby = 1
    CAS_active = 2
end

Base.@kwdef mutable struct AvionicsInterfaceU
    eng_start::Bool = false
    eng_stop::Bool = false
    CAS_enable::Bool = false
    throttle::Ranged{Float64, 0, 1} = 0.0
    mixture::Ranged{Float64, 0, 1} = 0.5
    roll_input::Ranged{Float64, -1, 1} = 0.0
    pitch_input::Ranged{Float64, -1, 1} = 0.0
    yaw_input::Ranged{Float64, -1, 1} = 0.0
    aileron_trim::Ranged{Float64, -1, 1} = 0.0 #only relevant with CAS disabled
    elevator_trim::Ranged{Float64, -1, 1} = 0.0 #only relevant with CAS disabled
    rudder_trim::Ranged{Float64, -1, 1} = 0.0 #only relevant with CAS disabled
    flaps::Ranged{Float64, 0, 1} = 0.0
    brake_left::Ranged{Float64, 0, 1} = 0.0
    brake_right::Ranged{Float64, 0, 1} = 0.0
end

Base.@kwdef struct AvionicsInterfaceY
    eng_start::Bool = false
    eng_stop::Bool = false
    CAS_enable::Bool = false
    throttle::Float64 = 0.0
    mixture::Float64 = 0.5
    roll_input::Float64 = 0.0
    pitch_input::Float64 = 0.0
    yaw_input::Float64 = 0.0
    aileron_trim::Float64 = 0.0
    elevator_trim::Float64 = 0.0
    rudder_trim::Float64 = 0.0
    flaps::Float64 = 0.0
    brake_left::Float64 = 0.0
    brake_right::Float64 = 0.0
end

Base.@kwdef struct AvionicsLogicY
    flight_phase::FlightPhase = phase_gnd
    CAS_state::CASState = CAS_disabled
end

Base.@kwdef struct Avionics <: AbstractAvionics
    p_cmd_sf::Float64 = 0.2 #roll_input to p_cmd scale factor (roll_input ∈ [-1, 1])
    q_cmd_sf::Float64 = 0.2 #pitch_input to q_cmd scale factor (pitch_input ∈ [-1, 1])
    β_cmd_sf::Float64 = 0.1 #yaw_input to β_cmd scale factor (yaw_input ∈ [-1, 1])
    p_control::RollRateControl = RollRateControl()
    q_control::PitchRateControl = PitchRateControl()
    β_control::SideslipControl = SideslipControl()
end

const AvionicsU = AvionicsInterfaceU

Base.@kwdef struct AvionicsY
    interface::AvionicsInterfaceY = AvionicsInterfaceY()
    logic::AvionicsLogicY = AvionicsLogicY()
    p_control::RollRateControlY = RollRateControlY()
    q_control::PitchRateControlY = PitchRateControlY()
    β_control::SideslipControlY = SideslipControlY()
end

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()
Systems.init(::SystemS, ::Avionics) = nothing #keep subsystems local


########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:Avionics}, Δt::Real,
                        airframe::System{<:C172RAirframe}, kinematics::KinematicData,
                        ::RigidBodyData, air::AirData, ::TerrainData)

    @unpack p_cmd_sf, q_cmd_sf, β_cmd_sf = avionics.params
    @unpack p_control, q_control, β_control = avionics.subsystems
    @unpack eng_start, eng_stop, CAS_enable,
            throttle, mixture, roll_input, pitch_input, yaw_input,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right = avionics.u

    nlg_wow = airframe.y.ldg.nose.strut.wow
    lmain_wow = airframe.y.ldg.left.strut.wow
    rmain_wow = airframe.y.ldg.right.strut.wow

    flight_phase = (nlg_wow && lmain_wow && rmain_wow) ? phase_air : phase_gnd

    if !CAS_enable
        CAS_state = CAS_disabled
    else #CAS enabled
        CAS_state = (flight_phase == phase_gnd ? CAS_standby : CAS_active)
    end

    p, q, _ = kinematics.ω_lb_b
    β = air.β_b

    # p_control.u.reset = (CAS_state === CAS_active ? false : true)
    # p_control.u.p_cmd = p_cmd_sf * Float64(roll_input)
    # p_control.u.p_fbk = p
    f_disc!(p_control, Δt)

    q_control.u.reset = (CAS_state === CAS_active ? false : true)
    q_control.u.q_cmd = q_cmd_sf * Float64(pitch_input)
    q_control.u.q_fbk = q
    f_disc!(q_control, Δt)

    # β_control.u.reset = (CAS_state === CAS_active ? false : true)
    # β_control.u.β_cmd = β_cmd * Float64(yaw_input)
    # β_control.u.β_fbk = β
    f_disc!(β_control, Δt)

    interface_y = AvionicsInterfaceY(;
            eng_start, eng_stop, CAS_enable, throttle, mixture,
            roll_input, pitch_input, yaw_input,
            aileron_trim, elevator_trim, rudder_trim,
            flaps, brake_left, brake_right)

    logic_y = AvionicsLogicY(; flight_phase, CAS_state)

    avionics.y = AvionicsY( interface = interface_y,
                            logic = logic_y,
                            p_control = p_control.y,
                            q_control = q_control.y,
                            β_control = β_control.y)

    return false

end

function Aircraft.map_controls!(airframe::System{<:C172RAirframe},
                                avionics::System{Avionics})

    @unpack eng_start, eng_stop, throttle, mixture,
            roll_input, pitch_input, yaw_input,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right = avionics.y.interface

    u_act = airframe.act.u

    if avionics.y.logic.CAS_state === CAS_active
        # u_act.aileron = avionics.y.p_control.a_cmd
        u_act.aileron = roll_input #roll_input maps directly to elevator
        u_act.elevator = avionics.y.q_control.e_cmd
        # u_act.rudder = avionics.y.β_control.r_cmd
        u_act.rudder = yaw_input #yaw_input maps directly to elevator
    else #disabled or standby
        u_act.aileron = roll_input #roll_input maps directly to elevator
        u_act.elevator = pitch_input #pitch_input maps directly to elevator
        u_act.rudder = yaw_input #yaw_input maps directly to elevator
    end

    #assign the rest of stuff
    @pack!  u_act = eng_start, eng_stop, throttle, mixture,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right

end


################################## GUI #########################################

function GUI.draw!(avionics::System{<:Avionics}, airframe::System{<:C172RAirframe},
                    label::String = "Cessna 172R Direct Controls")

    u = avionics.u

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    u.eng_start = dynamic_button("Engine Start", 0.4); CImGui.SameLine()
    u.eng_stop = dynamic_button("Engine Stop", 0.0)
    u.throttle = safe_slider("Throttle", u.throttle, "%.6f")
    u.aileron = safe_slider("Aileron", u.aileron, "%.6f")
    u.elevator = safe_slider("Elevator", u.elevator, "%.6f")
    u.rudder = safe_slider("Rudder", u.rudder, "%.6f")
    u.aileron_trim = safe_input("Aileron Trim", u.aileron_trim, 0.001, 0.1, "%.6f")
    u.elevator_trim = safe_input("Elevator Trim", u.elevator_trim, 0.001, 0.1, "%.6f")
    u.rudder_trim = safe_input("Rudder Trim", u.rudder_trim, 0.001, 0.1, "%.6f")
    u.flaps = safe_slider("Flaps", u.flaps, "%.6f")
    u.mixture = safe_slider("Mixture", u.mixture, "%.6f")
    u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
    u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")

    CImGui.Text(@sprintf("Engine Speed: %.3f RPM",
                        Piston.radpersec2RPM(airframe.y.pwp.engine.ω)))

    CImGui.PopItemWidth()

    CImGui.End()

end

################################################################################
############################# Cessna172Rv2 #####################################

#Cessna172R variant with Avionics avionics
const Cessna172Rv2{K, F} = AircraftTemplate{K, F, Avionics} where {K, F <: C172RAirframe}
Cessna172Rv2(kinematics = LTF()) = AircraftTemplate(kinematics, C172RAirframe(), Avionics())


############################ Joystick Mappings #################################

function IODevices.assign!(sys::System{<:Cessna172Rv2}, joystick::Joystick,
                           mapping::InputMapping)
    IODevices.assign!(sys.avionics, joystick, mapping)
end

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign!(sys::System{Avionics},
                           joystick::XBoxController,
                           ::DefaultMapping)

    u = sys.u

    u.aileron = get_axis_value(joystick, :right_analog_x) |> aileron_curve
    u.elevator = get_axis_value(joystick, :right_analog_y) |> elevator_curve
    u.rudder = get_axis_value(joystick, :left_analog_x) |> rudder_curve
    u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

    u.aileron_trim -= 0.01 * was_released(joystick, :dpad_left)
    u.aileron_trim += 0.01 * was_released(joystick, :dpad_right)
    u.elevator_trim += 0.01 * was_released(joystick, :dpad_down)
    u.elevator_trim -= 0.01 * was_released(joystick, :dpad_up)

    u.throttle += 0.1 * was_released(joystick, :button_Y)
    u.throttle -= 0.1 * was_released(joystick, :button_A)

    u.flaps += 0.3333 * was_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

end

function IODevices.assign!(sys::System{Avionics},
                           joystick::T16000M,
                           ::DefaultMapping)

    u = sys.u

    u.throttle = get_axis_value(joystick, :throttle)
    u.aileron = get_axis_value(joystick, :stick_x) |> aileron_curve
    u.elevator = get_axis_value(joystick, :stick_y) |> elevator_curve
    u.rudder = get_axis_value(joystick, :stick_z) |> rudder_curve

    u.brake_left = is_pressed(joystick, :button_1)
    u.brake_right = is_pressed(joystick, :button_1)

    u.aileron_trim -= 2e-4 * is_pressed(joystick, :hat_left)
    u.aileron_trim += 2e-4 * is_pressed(joystick, :hat_right)
    u.elevator_trim += 2e-4 * is_pressed(joystick, :hat_down)
    u.elevator_trim -= 2e-4 * is_pressed(joystick, :hat_up)

    u.flaps += 0.3333 * was_released(joystick, :button_3)
    u.flaps -= 0.3333 * was_released(joystick, :button_2)

end


end #module