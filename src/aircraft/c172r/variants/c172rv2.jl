module C172Rv2

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI
using Flight.FlightCore.IODevices
using Flight.FlightCore.Joysticks
using Flight.FlightCore.Utils: Ranged

using Flight.FlightPhysics.Attitude
using Flight.FlightPhysics.Kinematics
using Flight.FlightPhysics.RigidBody
using Flight.FlightPhysics.Environment

using Flight.FlightComponents.Control
using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.World
using Flight.FlightComponents.Control: PIDDiscreteY

using ..Airframe

export Cessna172Rv2

################################################################################
################################# Avionics #####################################

################################################################################
################################ PitchRateComp #################################

@enum PitchControlMode begin
    elevator_mode = 0
    pitch_rate_mode = 1
    inclination_mode = 2
end

@kwdef struct PitchRateComp <: SystemDefinition
    c1::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 0, k_i = 1, k_d = 0) #pure integrator
    c2::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 10, k_i = 20, k_d = 0.5, τ_d = 0.05, β_p = 1, β_d = 1) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct PitchRateCompU
    setpoint::Float64 = 0.0
    feedback::Float64 = 0.0
    reset::Bool = false
    sat_ext::Int64 = 0.0
end

@kwdef struct PitchRateCompY
    setpoint::Float64 = 0.0
    feedback::Float64 = 0.0
    reset::Bool = false
    sat_ext::Int64 = 0.0
    out::Float64 = 0.0 #elevator command
    c1::PIDDiscreteY{1} = PIDDiscreteY{1}()
    c2::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::PitchRateComp) = PitchRateCompU()
Systems.init(::SystemY, ::PitchRateComp) = PitchRateCompY()

#here, we are leaving the compensators' outputs unbounded, relying only on the
#external saturation input coming from the elevator to stop the integrators
function Systems.init!(sys::System{PitchRateComp})
    sys.c1.u.reset .= true
    sys.c2.u.reset .= true
    sys.c1.u.anti_windup .= true
    sys.c2.u.anti_windup .= true
end

#if c2's output saturates positively, then both c2 and c1's integrators must
#stop accumulating positively, and viceversa.
function Systems.f_disc!(sys::System{PitchRateComp}, Δt::Real)
    @unpack setpoint, feedback, reset, sat_ext = sys.u
    @unpack c1, c2 = sys.subsystems

    c1.u.setpoint .= setpoint
    c1.u.feedback .= feedback
    c1.u.reset .= reset
    c1.u.sat_ext .= sat_ext
    f_disc!(c1, Δt)

    c2.u.setpoint .= c1.y.out #connected to c1's output
    c2.u.feedback .= 0.0 #no feedback, just feedforward path
    c2.u.reset .= reset
    c2.u.sat_ext .= sat_ext
    f_disc!(c2, Δt)

    out = c2.y.out[1]

    sys.y = PitchRateCompY(; setpoint, feedback, reset, sat_ext, out, c1 = c1.y, c2 = c2.y)

end

################################################################################
############################### PitchAxisControl ###############################

@kwdef struct PitchControl <: SystemDefinition
    q_comp::PitchRateComp = PitchRateComp()
    θ_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 4, k_i = 0.3, k_d = 0.3, τ_d = 0.05, β_p = 1, β_d = 1) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct PitchControlU
    mode::PitchControlMode = elevator_mode
    e_cmd::Float64 = 0.0
    q_cmd::Float64 = 0.0
    θ_cmd::Float64 = 0.0
    q::Float64 = 0.0
    θ::Float64 = 0.0
end

@kwdef mutable struct PitchControlS
    mode_prev::PitchControlMode = elevator_mode
end

@kwdef struct PitchControlY
    mode::PitchControlMode = elevator_mode
    mode_prev::PitchControlMode = elevator_mode
    e_out::Ranged{Float64, -1., 1.} = 0.0 #elevator output
    e_sat::Int64 = 0 #elevator saturation state
    q_comp::PitchRateCompY = PitchRateCompY()
    θ_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::PitchControl) = PitchControlU()
Systems.init(::SystemY, ::PitchControl) = PitchControlY()
Systems.init(::SystemS, ::PitchControl) = PitchControlS()

function Systems.init!(sys::System{PitchControl})
    Systems.init!(sys.q_comp)
    sys.θ_comp.u.reset .= true
    sys.θ_comp.u.anti_windup .= true
end

function Systems.f_disc!(sys::System{PitchControl}, Δt::Real)

    @unpack mode, e_cmd, q_cmd, θ_cmd, q, θ = sys.u
    @unpack mode_prev = sys.s
    @unpack q_comp, θ_comp = sys.subsystems

    if mode != mode_prev #reset compensators on mode change
        q_comp.u.reset = true; f_disc!(q_comp, Δt)
        θ_comp.u.reset .= true; f_disc!(θ_comp, Δt)
    end

    if mode === elevator_mode
        e_out = Ranged(e_cmd, -1., 1.)
    else #pitch rate compensator active
        q_comp.u.reset = false
        q_comp.u.feedback = q
        if mode === pitch_rate_mode
            q_comp.u.setpoint = q_cmd
        else #inclination_mode
            θ_comp.u.reset .= false
            θ_comp.u.setpoint .= θ_cmd
            θ_comp.u.feedback .= θ
            f_disc!(θ_comp, Δt)
            q_comp.u.setpoint = θ_comp.y.out[1]
        end
        f_disc!(q_comp, Δt)
        e_out = Ranged(q_comp.y.out, -1., 1.)
    end

    e_sat = (e_out == typemax(e_out)) - (e_out == typemin(e_out))
    q_comp.u.sat_ext = e_sat #will take effect on the next call
    θ_comp.u.sat_ext .= e_sat #will take effect on the next call

    mode_prev = mode
    sys.s.mode_prev = mode_prev
    sys.y = PitchControlY(; mode, mode_prev, e_out, e_sat, q_comp = q_comp.y, θ_comp = θ_comp.y)

end

################################################################################
############################### RollAxisControl ################################

@enum RollControlMode begin
    aileron_mode = 0
    roll_rate_mode = 1
    bank_mode = 2
end

@kwdef struct RollControl <: SystemDefinition
    p_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 0.5, k_i = 10, k_d = 0.05, τ_d = 0.05) #roll rate compensator, see notebook
    k_φ::Float64 = 2.0 #bank gain (no compensation, only proportional control needed)
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct RollControlU
    mode::RollControlMode = aileron_mode
    a_cmd::Float64 = 0.0
    p_cmd::Float64 = 0.0
    φ_cmd::Float64 = 0.0
    p::Float64 = 0.0
    φ::Float64 = 0.0
end

@kwdef struct RollControlY
    mode::RollControlMode = aileron_mode
    a_out::Ranged{Float64, -1., 1.} = 0.0 #direct aileron command
    a_sat::Int64 = 0 #aileron saturation state
    p_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::RollControl) = RollControlU()
Systems.init(::SystemY, ::RollControl) = RollControlY()

function Systems.init!(sys::System{RollControl})
    sys.p_comp.u.reset .= true
    sys.p_comp.u.anti_windup .= true
    println("Bank gain adjustment pending")
end

function Systems.f_disc!(sys::System{RollControl}, Δt::Real)

    @unpack mode, a_cmd, p_cmd, φ_cmd, p, φ = sys.u
    @unpack p_comp = sys.subsystems
    @unpack k_φ = sys.params

    if mode == aileron_mode
        p_comp.u.reset .= true
        f_disc!(p_comp, Δt)
        a_out = Ranged(a_cmd, -1., 1.)
    else
        p_comp.u.reset .= false
        if mode == roll_rate_mode
            p_comp.u.setpoint .= p_cmd
        else #mode == inclination_mode
            p_comp.u.setpoint .= k_φ * (φ_cmd - φ)
        end
        p_comp.u.feedback .= p
        f_disc!(p_comp, Δt)
        a_out = Ranged(p_comp.y.out[1], -1., 1.)
    end

    a_sat = (a_out == typemax(a_out)) - (a_out == typemin(a_out))
    p_comp.u.sat_ext .= a_sat #will take effect on the next call

    sys.y = RollControlY(; mode, a_out, a_sat, p_comp = p_comp.y)

end

################################################################################
############################## YawAxisControl ##################################

@enum YawControlMode begin
    rudder_mode = 0
    sideslip_mode = 1
end

@kwdef struct YawControl <: SystemDefinition
    β_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 10, k_i = 25, k_d = 5, τ_d = 0.05) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct YawControlU
    mode::YawControlMode = rudder_mode
    r_cmd::Float64 = 0.0
    β_cmd::Float64 = 0.0
    β::Float64 = 0.0
end

@kwdef struct YawControlY
    mode::YawControlMode = rudder_mode
    r_out::Ranged{Float64, -1., 1.} = 0.0 #rudder output
    r_sat::Int64 = 0 #rudder saturation state
    β_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::YawControl) = YawControlU()
Systems.init(::SystemY, ::YawControl) = YawControlY()

function Systems.init!(sys::System{YawControl})
    sys.β_comp.u.reset .= true
    sys.β_comp.u.anti_windup .= true
end

function Systems.f_disc!(sys::System{YawControl}, Δt::Real)
    @unpack mode, r_cmd, β_cmd, β = sys.u
    @unpack β_comp = sys.subsystems

    if mode == rudder_mode
        β_comp.u.reset .= true
        f_disc!(β_comp, Δt)
        r_out = Ranged(r_cmd, -1., 1.)
    else #mode == sideslip_mode
        β_comp.u.reset .= false
        β_comp.u.setpoint .= β_cmd
        β_comp.u.feedback .= β
        f_disc!(β_comp, Δt)
        r_out = Ranged(-β_comp.y.out[1], -1., 1.) #note sign inversion, see design notebook
    end

    #rudder output is inverted from β_comp's output, so we need to invert the
    #rudder saturation signal as well before assigning it back to β_comp
    r_sat = (r_out == typemax(r_out)) - (r_out == typemin(r_out))
    β_comp.u.sat_ext .= -r_sat #will take effect on the next call

    sys.y = YawControlY(; mode, r_out, r_sat, β_comp = β_comp.y)

end

##################################################################################
################################## Avionics ######################################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@enum CASState begin
    CAS_disabled = 0
    CAS_standby = 1
    CAS_active = 2
end

@kwdef mutable struct AvionicsInterfaceU
    eng_start::Bool = false
    eng_stop::Bool = false
    CAS_enable::Bool = false
    roll_control_mode_select::RollControlMode = aileron_mode #selected roll axis mode
    pitch_control_mode_select::PitchControlMode = elevator_mode #selected pitch axis mode
    yaw_control_mode_select::YawControlMode = rudder_mode #selected yaw axis mode
    throttle::Ranged{Float64, 0., 1.} = 0.0
    mixture::Ranged{Float64, 0., 1.} = 0.5
    roll_input::Ranged{Float64, -1., 1.} = 0.0
    pitch_input::Ranged{Float64, -1., 1.} = 0.0
    yaw_input::Ranged{Float64, -1., 1.} = 0.0
    aileron_trim::Ranged{Float64, -1., 1.} = 0.0 #only relevant with CAS disabled
    elevator_trim::Ranged{Float64, -1., 1.} = 0.0 #only relevant with CAS disabled
    rudder_trim::Ranged{Float64, -1., 1.} = 0.0 #only relevant with CAS disabled
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end

@kwdef struct AvionicsInterfaceY
    eng_start::Bool = false
    eng_stop::Bool = false
    CAS_enable::Bool = false
    roll_control_mode::RollControlMode = aileron_mode #actual roll axis mode
    pitch_control_mode::PitchControlMode = elevator_mode #actual pitch axis mode
    yaw_control_mode::YawControlMode = rudder_mode #actual yaw axis mode
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

@kwdef struct AvionicsLogicY
    flight_phase::FlightPhase = phase_gnd
    CAS_state::CASState = CAS_disabled
end

@kwdef struct Avionics <: AbstractAvionics
    p_input_sf::Float64 = 0.2 #external roll axis input to p_cmd scale factor
    q_input_sf::Float64 = 0.2 #external pitch axis input to q_cmd scale factor
    φ_input_sf::Float64 = deg2rad(60) #external roll axis input to φ_cmd scale factor
    θ_input_sf::Float64 = deg2rad(20) #external pitch axis input to θ_cmd scale factor
    β_input_sf::Float64 = deg2rad(10) #external yaw axis input to β_cmd scale factor
    roll_control::RollControl = RollControl()
    pitch_control::PitchControl = PitchControl()
    yaw_control::YawControl = YawControl()
end

const AvionicsU = AvionicsInterfaceU

@kwdef struct AvionicsY
    interface::AvionicsInterfaceY = AvionicsInterfaceY()
    logic::AvionicsLogicY = AvionicsLogicY()
    roll_control::RollControlY = RollControlY()
    pitch_control::PitchControlY = PitchControlY()
    yaw_control::YawControlY = YawControlY()
end

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()
Systems.init(::SystemS, ::Avionics) = nothing #keep subsystems local


########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:Avionics}, Δt::Real,
                        airframe::System{<:C172RAirframe}, kinematics::KinematicData,
                        ::RigidBodyData, air::AirData, ::TerrainData)

    @unpack roll_control, pitch_control, yaw_control = avionics.subsystems
    @unpack eng_start, eng_stop, CAS_enable,
            roll_control_mode_select, pitch_control_mode_select, yaw_control_mode_select,
            throttle, mixture, roll_input, pitch_input, yaw_input,
            aileron_trim, elevator_trim, rudder_trim, flaps, brake_left, brake_right = avionics.u
    @unpack p_input_sf, q_input_sf, φ_input_sf, θ_input_sf, β_input_sf = avionics.params

    nlg_wow = airframe.y.ldg.nose.strut.wow
    lmain_wow = airframe.y.ldg.left.strut.wow
    rmain_wow = airframe.y.ldg.right.strut.wow

    flight_phase = (!nlg_wow && !lmain_wow && !rmain_wow) ? phase_air : phase_gnd

    if !CAS_enable
        CAS_state = CAS_disabled
    else #CAS enabled
        CAS_state = (flight_phase == phase_gnd ? CAS_standby : CAS_active)
    end

    if CAS_state == CAS_active
        roll_control_mode = roll_control_mode_select
        pitch_control_mode = pitch_control_mode_select
        yaw_control_mode = yaw_control_mode_select
    else
        roll_control_mode = aileron_mode
        pitch_control_mode = elevator_mode
        yaw_control_mode = rudder_mode
    end

    e_nb = REuler(kinematics.q_nb)
    @unpack θ, φ = e_nb
    p, q, _ = kinematics.ω_lb_b
    β = air.β_b

    a_cmd = Float64(roll_input) #already ∈ [-1, 1]
    p_cmd = p_input_sf * Float64(roll_input)
    φ_cmd = φ_input_sf * Float64(roll_input)
    roll_control.u.mode = roll_control_mode
    @pack! roll_control.u = a_cmd, p_cmd, φ_cmd, p, φ
    f_disc!(roll_control, Δt)

    e_cmd = Float64(pitch_input) #already ∈ [-1, 1]
    q_cmd = q_input_sf * Float64(pitch_input)
    θ_cmd = θ_input_sf * Float64(pitch_input)
    pitch_control.u.mode = pitch_control_mode
    @pack! pitch_control.u = e_cmd, q_cmd, θ_cmd, q, θ
    f_disc!(pitch_control, Δt)

    r_cmd = Float64(yaw_input) #already ∈ [-1, 1]
    β_cmd = β_input_sf * Float64(yaw_input)
    yaw_control.u.mode = yaw_control_mode
    @pack! yaw_control.u = r_cmd, β_cmd, β
    f_disc!(yaw_control, Δt)

    # @show β
    # @show θ

    interface_y = AvionicsInterfaceY(;
            eng_start, eng_stop, CAS_enable,
            roll_control_mode, pitch_control_mode, yaw_control_mode,
            throttle, mixture, roll_input, pitch_input, yaw_input,
            aileron_trim, elevator_trim, rudder_trim,
            flaps, brake_left, brake_right)

    logic_y = AvionicsLogicY(; flight_phase, CAS_state)

    avionics.y = AvionicsY( interface = interface_y, logic = logic_y,
                            roll_control = roll_control.y,
                            pitch_control = pitch_control.y,
                            yaw_control = yaw_control.y)

    return false

end

function Aircraft.map_controls!(airframe::System{<:C172RAirframe},
                                avionics::System{Avionics})

    @unpack eng_start, eng_stop, throttle, mixture,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right = avionics.y.interface

    u_act = airframe.act.u

    u_act.aileron = avionics.y.roll_control.a_out
    u_act.elevator = avionics.y.pitch_control.e_out
    u_act.rudder = avionics.y.yaw_control.r_out

    @pack!  u_act = eng_start, eng_stop, throttle, mixture,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right

end


################################## GUI #########################################

function GUI.draw!(avionics::System{<:Avionics}, airframe::System{<:C172RAirframe},
                    label::String = "Cessna 172R CAS Avionics")

    u = avionics.u
    y = avionics.y

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    u.eng_start = dynamic_button("Engine Start", 0.4); CImGui.SameLine()
    u.eng_stop = dynamic_button("Engine Stop", 0.0); CImGui.SameLine()
    CImGui.Text(@sprintf("Engine Speed: %.3f RPM", Piston.radpersec2RPM(airframe.y.pwp.engine.ω)))

    u.CAS_enable = toggle_switch("CAS", 0.4, u.CAS_enable)
    CImGui.Text("Flight Phase: $(y.logic.flight_phase)")
    CImGui.Text("CAS State: $(y.logic.CAS_state)")

    #maybe make the displayed variables depend on CAS state and mode
    #(aileron input vs roll rate demand vs bank angle demand)

    u.throttle = safe_slider("Throttle", u.throttle, "%.6f")
    u.roll_input = safe_slider("Roll Input", u.roll_input, "%.6f")
    u.pitch_input = safe_slider("Pitch Input", u.pitch_input, "%.6f")
    u.yaw_input = safe_slider("Yaw Input", u.yaw_input, "%.6f")
    u.aileron_trim = safe_input("Aileron Trim", u.aileron_trim, 0.001, 0.1, "%.6f")
    u.elevator_trim = safe_input("Elevator Trim", u.elevator_trim, 0.001, 0.1, "%.6f")
    u.rudder_trim = safe_input("Rudder Trim", u.rudder_trim, 0.001, 0.1, "%.6f")
    u.flaps = safe_slider("Flaps", u.flaps, "%.6f")
    u.mixture = safe_slider("Mixture", u.mixture, "%.6f")
    u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
    u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")


    CImGui.PopItemWidth()

    CImGui.End()

end

################################################################################
############################# Cessna172Rv2 #####################################

#Cessna172R variant with Avionics avionics
const Cessna172Rv2{K, F} = AircraftTemplate{K, F, Avionics} where {K, F <: C172RAirframe}
Cessna172Rv2(kinematics = LTF()) = AircraftTemplate(kinematics, C172RAirframe(), Avionics())


# ############################ Joystick Mappings #################################

# function IODevices.assign!(sys::System{<:Cessna172Rv2}, joystick::Joystick,
#                            mapping::InputMapping)
#     IODevices.assign!(sys.avionics, joystick, mapping)
# end

# elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
# aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
# rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
# brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

# function IODevices.assign!(sys::System{Avionics},
#                            joystick::XBoxController,
#                            ::DefaultMapping)

#     u = sys.u

#     u.roll_input = get_axis_value(joystick, :right_analog_x) |> aileron_curve
#     u.pitch_input = get_axis_value(joystick, :right_analog_y) |> elevator_curve
#     u.yaw_input = get_axis_value(joystick, :left_analog_x) |> rudder_curve
#     u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
#     u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

#     u.aileron_trim -= 0.01 * was_released(joystick, :dpad_left)
#     u.aileron_trim += 0.01 * was_released(joystick, :dpad_right)
#     u.elevator_trim += 0.01 * was_released(joystick, :dpad_down)
#     u.elevator_trim -= 0.01 * was_released(joystick, :dpad_up)

#     u.throttle += 0.1 * was_released(joystick, :button_Y)
#     u.throttle -= 0.1 * was_released(joystick, :button_A)

#     u.flaps += 0.3333 * was_released(joystick, :right_bumper)
#     u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

# end

# function IODevices.assign!(sys::System{Avionics},
#                            joystick::T16000M,
#                            ::DefaultMapping)

#     u = sys.u

#     u.throttle = get_axis_value(joystick, :throttle)
#     u.roll_input = get_axis_value(joystick, :stick_x) |> aileron_curve
#     u.pitch_input = get_axis_value(joystick, :stick_y) |> elevator_curve
#     u.yaw_input = get_axis_value(joystick, :stick_z) |> rudder_curve

#     u.brake_left = is_pressed(joystick, :button_1)
#     u.brake_right = is_pressed(joystick, :button_1)

#     u.aileron_trim -= 2e-4 * is_pressed(joystick, :hat_left)
#     u.aileron_trim += 2e-4 * is_pressed(joystick, :hat_right)
#     u.elevator_trim += 2e-4 * is_pressed(joystick, :hat_down)
#     u.elevator_trim -= 2e-4 * is_pressed(joystick, :hat_up)

#     u.flaps += 0.3333 * was_released(joystick, :button_3)
#     u.flaps -= 0.3333 * was_released(joystick, :button_2)

# end


end #module