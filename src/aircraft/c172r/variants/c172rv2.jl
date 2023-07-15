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

using Flight.FlightComponents.Control
using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.World
using Flight.FlightComponents.Control: PIDDiscreteY

using ..Airframe

export Cessna172Rv2

################################################################################
################################# Avionics #####################################

################################################################################
################################# PitchAxisControl #############################

@enum PitchAxisMode begin
    elevator_mode = 0
    pitch_rate_mode = 1
    inclination_mode = 2
end

############################# PitchRateComp #################################

@kwdef struct PitchRateComp <: SystemDefinition
    c1::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 0, k_i = 1, k_d = 0) #pure integrator
    c2::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 10, k_i = 20, k_d = 0.5, τ_d = 0.05, β_p = 1, β_d = 1) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct PitchRateCompU
    setpoint::Float64 = 0.0
    feedback::Float64 = 0.0
    reset::Bool = false
end

@kwdef struct PitchRateCompY
    setpoint::Float64 = 0.0
    feedback::Float64 = 0.0
    reset::Bool = false
    out::Float64 = 0.0 #elevator command
    sat_out::Int64 = 0 #elevator saturation
    c1::PIDDiscreteY{1} = PIDDiscreteY{1}()
    c2::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::PitchRateComp) = PitchRateCompU()
Systems.init(::SystemY, ::PitchRateComp) = PitchRateCompY()

function Systems.init!(sys::System{PitchRateComp})
    @unpack c1, c2 = sys.subsystems
    c1.u.bound_lo .= -Inf
    c1.u.bound_hi .= Inf
    c1.u.anti_windup .= true

    c2.u.bound_lo .= -1 #lower bound for MechanicalActuation's normalized elevator input
    c2.u.bound_hi .= 1 #upper bound for MechanicalActuation's normalized elevator input
    c2.u.anti_windup .= true
end

function Systems.f_disc!(sys::System{PitchRateComp}, Δt::Real)
    @unpack setpoint, feedback, reset = sys.u
    @unpack c1, c2 = sys.subsystems

    c1.u.setpoint .= setpoint
    c1.u.feedback .= feedback
    c1.u.reset .= reset
    c1.u.sat_ext .= c2.y.sat_out
    f_disc!(c1, Δt)

    c2.u.setpoint .= c1.y.out #connected to c1's output
    c2.u.feedback .= 0.0 #no feedback, just feedforward path
    c2.u.reset .= reset
    c2.u.sat_ext .= 0 #no external saturation signal required
    f_disc!(c2, Δt)

    out = c2.y.out[1]
    sat_out = c2.y.sat_out[1]

    sys.y = PitchRateCompY(; setpoint, feedback, reset, out, sat_out, c1 = c1.y, c2 = c2.y)

end


############################## PitchAxisControl ################################

@kwdef struct PitchAxisControl <: SystemDefinition
    q_comp::PitchRateCompensator = PitchRateCompensator()
    k_θ::Float64 = 2.8 #θ loop gain, see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct PitchAxisControlU
    mode::PitchAxisMode = elevator_mode
    e_cmd::Float64 = 0.0 #elevator command
    q_cmd::Float64 = 0.0
    θ_cmd::Float64 = 0.0
    q::Float64 = 0.0
    θ::Float64 = 0.0
end

@kwdef struct PitchAxisControlY
    mode::PitchAxisMode = elevator_mode
    e_out::Float64 = 0.0 #elevator output
    e_sat::Int64 = 0 #elevator saturation state
    q_cmp::PitchRateCompY = PitchRateCompY()
end

Systems.init(::SystemU, ::PitchAxisControl) = PitchAxisControlU()
Systems.init(::SystemY, ::PitchAxisControl) = PitchAxisControlY()

function Systems.f_disc!(sys::System{PitchAxisControl}, Δt::Real)

    @unpack mode, e_cmd, q_cmd, θ_cmd, q, θ = sys.u
    @unpack q_comp = sys.subsystems
    @unpack k_θ = sys.params

    if mode == elevator_mode
        q_comp.u.reset = true
        f_disc!(q_comp, Δt)
        e_out = u.e_cmd #elevator output is directly elevator command
        e_sat = 0 #no saturation for direct elevator command
    else
        q_comp.u.reset = false
        if mode == pitch_rate_mode
            q_comp.u.setpoint = q_cmd
        else #mode == inclination_mode
            q_comp.u.setpoint = k_θ * (θ_cmd - θ)
        end
        q_comp.u.feedback = q
        f_disc!(q_comp, Δt)
        e_out = q_comp.y.out
        e_sat = q_comp.y.sat_out
    end

    sys.y = PitchAxisControlY(; mode, e_cmd, e_sat, q_cmp = q_cmp.y)

end


################################################################################
###############################

@enum RollAxisMode begin
    aileron_mode = 0
    roll_rate_mode = 1
    bank_mode = 2
end

############################## RollRateComp #################################

@kwdef struct RollRateCmp <: SystemDefinition
    c::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 0.5, k_i = 10, k_d = 0.05, τ_d = 0.05) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct RollRateCmpU
    p_cmd::Float64 = 0.0
    p_fbk::Float64 = 0.0
    reset::Bool = false
end

@kwdef struct RollRateCmpY
    p_cmd::Float64 = 0.0
    p_fbk::Float64 = 0.0
    a_cmd::Float64 = 0.0 #aileron command
    a_sat::Int64 = 0
    reset::Bool = false
    c::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::RollRateCmp) = RollRateCmpU()
Systems.init(::SystemY, ::RollRateCmp) = RollRateCmpY()

function Systems.init!(sys::System{RollRateCmp})
    @unpack c = sys.subsystems
    c.u.bound_lo .= -1 #lower bound for MechanicalActuation's normalized aileron input
    c.u.bound_hi .= 1 #upper bound for MechanicalActuation's normalized aileron input
    c.u.sat_ext .= 0 #only output saturation required
    c.u.anti_windup .= true
end

function Systems.f_disc!(sys::System{RollRateCmp}, Δt::Real)
    @unpack p_cmd, p_fbk, reset = sys.u
    @unpack c = sys.subsystems

    c.u.setpoint .= p_cmd
    c.u.feedback .= p_fbk
    c.u.reset .= reset
    f_disc!(c, Δt)

    a_cmd = c.y.out[1]

    sys.y = RollRateCmpY(; p_cmd, p_fbk, a_cmd, reset, c = c.y)

end

############################## BankCmp ##################################

#this is just a proportional compensator, we don't need all the bells and
#whistles of a PIDDiscrete System
@kwdef struct BankCmp <: SystemDefinition
    k_φ::Float64 = 1 #see notebook
end

@kwdef mutable struct BankCmpU
    φ_cmd::Float64 = 0.0
    φ_fbk::Float64 = 0.0
end

@kwdef struct BankCmpY
    φ_cmd::Float64 = 0.0
    φ_fbk::Float64 = 0.0
    p_cmd::Float64 = 0.0
end

Systems.init(::SystemU, ::BankCmp) = BankCmpU()
Systems.init(::SystemY, ::BankCmp) = BankCmpY()

function Systems.init!(::System{BankCmp})
    println("Bank gain adjustment pending")
end

function Systems.f_disc!(sys::System{BankCmp}, ::Real)
    @unpack φ_cmd, φ_fbk = sys.u
    p_cmd = sys.params.k_φ * (φ_cmd - φ_fbk)
    sys.y = BankCmpY(; φ_cmd, φ_fbk, p_cmd)
end

############################## RollAxisControl ################################

@kwdef struct RollAxisControl <: SystemDefinition
    p_input_sf::Float64 = 0.2 #external roll axis input to p_cmd scale factor
    φ_input_sf::Float64 = 0.4 #external roll axis input to φ_cmd scale factor
    p_cmp::RollRateCmp = RollRateCmp()
    φ_cmp::BankCmp = BankCmp()
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct RollAxisControlU
    mode::RollAxisMode = aileron_mode
    input::Ranged{Float64, -1, 1} = 0.0
    p_fbk::Float64 = 0.0
    φ_fbk::Float64 = 0.0
    reset::Bool = false
end

@kwdef struct RollAxisControlY
    a_cmd::Float64 = 0.0 #aileron command
    a_sat::Int64 = 0 #aileron saturation
    reset::Bool = false
    p_cmp::RollRateCmpY = RollRateCmpY()
    φ_cmp::BankCmpY = BankCmpY()
end

Systems.init(::SystemU, ::RollAxisControl) = RollAxisControlU()
Systems.init(::SystemY, ::RollAxisControl) = RollAxisControlY()

function Systems.f_disc!(sys::System{RollAxisControl}, Δt::Real)

    @unpack mode, input, p_fbk, φ_fbk, reset = sys.u
    @unpack p_cmp, φ_cmp = sys.subsystems
    @unpack p_input_sf, φ_input_sf = sys.params

    if mode == aileron_mode
        p_cmp.u.reset = true
        f_disc!(p_cmp, Δt)
        a_cmd = Float64(input) #input ∈ [-1,1]
        a_sat = 0 #not applicable with direct aileron command

    elseif mode == roll_rate_mode
        p_cmp.u.reset = reset #only reset on external reset input
        p_cmp.u.p_cmd = p_input_sf * Float64(input)
        p_cmp.u.p_fbk = p_fbk
        f_disc!(p_cmp, Δt)
        a_cmd = p_cmp.y.a_cmd
        a_sat = p_cmp.y.a_sat

    else #bank_mode
        φ_cmp.u.φ_cmd = φ_input_sf * Float64(input)
        φ_cmp.u.φ_fbk = φ_fbk
        f_disc!(φ_cmp, Δt)
        p_cmp.u.reset = reset #only reset on external reset input
        p_cmp.u.p_cmd = φ_cmp.y.p_cmd
        p_cmp.u.p_fbk = p_fbk
        f_disc!(p_cmp, Δt)
        a_cmd = p_cmp.y.a_cmd
        a_sat = p_cmp.y.a_sat
    end

    sys.y = RollAxisControlY(; a_cmd, a_sat, reset,
                                p_cmp = p_cmp.y,
                                φ_cmp = φ_cmp.y)

end


################################################################################
############################## YawAxisControl ##################################

@enum YawAxisMode begin
    rudder_mode = 0
    sideslip_mode = 1
end

############################## SideslipCmp #################################

#rationale for beta control in the inner CAS is that the user is likely to
#desire automatic turn combination in conjunction with roll rate and pitch rate
#augmentation, while yaw rate augmentation is not useful by itself
@kwdef struct SideslipCmp <: SystemDefinition
    c::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 10, k_i = 25, k_d = 5, τ_d = 0.05) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct SideslipCmpU
    β_cmd::Float64 = 0.0
    β_fbk::Float64 = 0.0
    reset::Bool = false
end

@kwdef struct SideslipCmpY
    β_cmd::Float64 = 0.0
    β_fbk::Float64 = 0.0
    r_cmd::Float64 = 0.0 #rudder command
    reset::Bool = false
    c::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::SideslipCmp) = SideslipCmpU()
Systems.init(::SystemY, ::SideslipCmp) = SideslipCmpY()

function Systems.init!(sys::System{SideslipCmp})
    @unpack c = sys.subsystems
    c.u.bound_lo .= -1 #lower bound for MechanicalActuation's normalized aileron input
    c.u.bound_hi .= 1 #upper bound for MechanicalActuation's normalized aileron input
    c.u.sat_ext .= 0 #only output saturation required
    c.u.anti_windup .= true
end

function Systems.f_disc!(sys::System{SideslipCmp}, Δt::Real)
    @unpack β_cmd, β_fbk, reset = sys.u
    @unpack c = sys.subsystems

    c.u.setpoint .= β_cmd
    c.u.feedback .= β_fbk
    c.u.reset .= reset
    f_disc!(c, Δt)

    #note the sign inversion, see design notebook!
    r_cmd = -c.y.out[1]

    sys.y = SideslipCmpY(; β_cmd, β_fbk, r_cmd, reset, c = c.y)

end

############################ YawAxisControl ####################################

@kwdef struct YawAxisControl <: SystemDefinition
    β_input_sf::Float64 = 0.2 #external yaw axis input to β_cmd scale factor
    β_cmp::SideslipCmp = SideslipCmp()
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct YawAxisControlU
    mode::YawAxisMode = rudder_mode
    input::Ranged{Float64, -1, 1} = 0.0
    β_fbk::Float64 = 0.0
    reset::Bool = false
end

@kwdef struct YawAxisControlY
    r_cmd::Float64 = 0.0 #rudder command
    r_sat::Int64 = 0 #rudder saturation
    reset::Bool = false
    β_cmp::SideslipCmpY = SideslipCmpY()
end

Systems.init(::SystemU, ::YawAxisControl) = YawAxisControlU()
Systems.init(::SystemY, ::YawAxisControl) = YawAxisControlY()

function Systems.init!(sys::System{YawAxisControl})
    Systems.init!(sys.β_cmp)
end

function Systems.f_disc!(sys::System{YawAxisControl}, Δt::Real)

    @unpack mode, input, β_fbk, reset = sys.u
    @unpack β_cmp = sys.subsystems
    @unpack β_input_sf = sys.params

    if mode == rudder_mode
        β_cmp.u.reset = true
        f_disc!(β_cmp, Δt)
        r_cmd = Float64(input) #input ∈ [-1,1]
        r_sat = 0 #not applicable with direct aileron command

    else # mode == sideslip_mode
        β_cmp.u.reset = reset #only reset on external reset input
        β_cmp.u.p_cmd = β_input_sf * Float64(input)
        β_cmp.u.p_fbk = β_fbk
        f_disc!(β_cmp, Δt)
        r_cmd = β_cmp.y.r_cmd
        r_sat = β_cmp.y.r_sat
    end

    sys.y = YawAxisControlY(; r_cmd, r_sat, reset, β_cmp = β_cmp.y)

end

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


@kwdef mutable struct AvionicsInterfaceU
    eng_start::Bool = false
    eng_stop::Bool = false
    throttle::Ranged{Float64, 0, 1} = 0.0
    mixture::Ranged{Float64, 0, 1} = 0.5
    CAS_enable::Bool = false
    roll_axis_mode_select::RollAxisMode = aileron_mode #selected roll axis mode
    pitch_axis_mode_select::PitchAxisMode = elevator_mode #selected pitch axis mode
    yaw_axis_mode_select::YawAxisMode = rudder_mode #selected yaw axis mode
    roll_axis_input::Ranged{Float64, -1, 1} = 0.0
    pitch_axis_input::Ranged{Float64, -1, 1} = 0.0
    yaw_axis_input::Ranged{Float64, -1, 1} = 0.0
    q_input_sf::Float64 = 0.2 #external pitch axis input to q_cmd scale factor
    θ_input_sf::Float64 = 0.4 #external pitch axis input to θ_cmd scale factor
    aileron_trim::Ranged{Float64, -1, 1} = 0.0 #only relevant with CAS disabled
    elevator_trim::Ranged{Float64, -1, 1} = 0.0 #only relevant with CAS disabled
    rudder_trim::Ranged{Float64, -1, 1} = 0.0 #only relevant with CAS disabled
    flaps::Ranged{Float64, 0, 1} = 0.0
    brake_left::Ranged{Float64, 0, 1} = 0.0
    brake_right::Ranged{Float64, 0, 1} = 0.0
end

@kwdef struct AvionicsInterfaceY
    eng_start::Bool = false
    eng_stop::Bool = false
    throttle::Float64 = 0.0
    mixture::Float64 = 0.5
    CAS_enable::Bool = false
    roll_axis_mode::RollAxisMode = aileron_mode #actual roll axis mode
    pitch_axis_mode::PitchAxisMode = elevator_mode #actual pitch axis mode
    yaw_axis_mode::PitchAxisMode = rudder_mode #actual yaw axis mode
    roll_axis_input::Float64 = 0.0
    pitch_axis_input::Float64 = 0.0
    yaw_axis_input::Float64 = 0.0
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
    roll_axis_control::RollAxisControl = RollAxisControl()
    pitch_axis_control::PitchAxisControl = PitchAxisControl()
    yaw_axis_control::YawAxisControl = YawAxisControl()
end

const AvionicsU = AvionicsInterfaceU

@kwdef struct AvionicsY
    interface::AvionicsInterfaceY = AvionicsInterfaceY()
    logic::AvionicsLogicY = AvionicsLogicY()
    roll_axis_control::RollAxisControlY = RollAxisControlY()
    pitch_axis_control::PitchAxisControlY = PitchAxisControlY()
    yaw_axis_control::YawAxisControlY = YawAxisControlY()
end

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()
Systems.init(::SystemS, ::Avionics) = nothing #keep subsystems local


########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:Avionics}, Δt::Real,
                        airframe::System{<:C172RAirframe}, kinematics::KinematicData,
                        ::RigidBodyData, air::AirData, ::TerrainData)

    @unpack roll_axis_control, pitch_axis_control, yaw_axis_control = avionics.subsystems
    @unpack eng_start, eng_stop, throttle, mixture,
            CAS_enable, roll_axis_mode, pitch_axis_mode, yaw_axis_mode,
            roll_axis_input, pitch_axis_input, yaw_axis_input,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right = avionics.u

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
        roll_axis_mode = roll_axis_mode_select
        pitch_axis_mode = pitch_axis_mode_select
        yaw_axis_mode = yaw_axis_mode_select
    else
        roll_axis_mode = aileron_mode
        pitch_axis_mode = elevator_mode
        yaw_axis_mode = rudder_mode
    end

    e_nb = REuler(kinematics.q_nb)
    @unpack θ, φ = e_nb
    p, q, _ = kinematics.ω_lb_b
    β = air.β_b

    # roll_axis_control.u.mode = roll_axis_mode
    # roll_axis_control.u.input = roll_axis_input
    # roll_axis_control.u.p_fbk = p
    # roll_axis_control.u.φ_fbk = φ
    # f_disc!(roll_axis_control, Δt)

    pitch_axis.u.mode = pitch_axis_mode
    pitch_axis.u.e_cmd = Float64(pitch_axis_input) #already ∈ [-1, 1]
    pitch_axis.u.q_cmd = q_input_sf * Float64(pitch_axis_input)
    pitch_axis.u.θ_cmd = θ_input_sf * Float64(pitch_axis_input)
    pitch_axis.u.q = q
    pitch_axis.u.θ = θ
    f_disc!(pitch_axis, Δt)

    # yaw_axis_control.u.mode = yaw_axis_mode
    # yaw_axis_control.u.input = yaw_axis_input
    # yaw_axis_control.u.β_fbk = β
    # f_disc!(yaw_axis_control, Δt)

    # @show CAS_state
    # @show q_control.y.reset
    # @show q_control.y.e_cmd
    # @show p_control.y.a_cmd
    # @show β_control.y.β_cmd
    # @show β_control.y.β_fbk
    # @show β_control.y.r_cmd
    # @show q_control.y.q_cmd
    # @show q_control.y.q_fbk
    # @show q_control.y.e_cmd

    interface_y = AvionicsInterfaceY(;
            eng_start, eng_stop, throttle, mixture,
            CAS_enable, roll_axis_mode, pitch_axis_mode, yaw_axis_mode,
            roll_axis_input, pitch_axis_input, yaw_axis_input,
            aileron_trim, elevator_trim, rudder_trim,
            flaps, brake_left, brake_right)

    logic_y = AvionicsLogicY(; flight_phase, CAS_state)

    avionics.y = AvionicsY( interface = interface_y,
                            logic = logic_y,
                            roll_axis_control = roll_axis_control.y,
                            pitch_axis_control = pitch_axis_control.y,
                            yaw_axis_control = yaw_axis_control.y)

    return false

end

function Aircraft.map_controls!(airframe::System{<:C172RAirframe},
                                avionics::System{Avionics})

    @unpack eng_start, eng_stop, throttle, mixture,
            roll_input, pitch_input, yaw_input,
            aileron_trim, elevator_trim, rudder_trim, flaps,
            brake_left, brake_right = avionics.y.interface

    u_act = airframe.act.u

    u_act.aileron = avionics.y.roll_axis_control.a_cmd
    u_act.elevator = avionics.y.pitch_axis_control.e_cmd
    u_act.rudder = avionics.y.yaw_axis_control.r_cmd

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

    u.roll_input = get_axis_value(joystick, :right_analog_x) |> aileron_curve
    u.pitch_input = get_axis_value(joystick, :right_analog_y) |> elevator_curve
    u.yaw_input = get_axis_value(joystick, :left_analog_x) |> rudder_curve
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
    u.roll_input = get_axis_value(joystick, :stick_x) |> aileron_curve
    u.pitch_input = get_axis_value(joystick, :stick_y) |> elevator_curve
    u.yaw_input = get_axis_value(joystick, :stick_z) |> rudder_curve

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