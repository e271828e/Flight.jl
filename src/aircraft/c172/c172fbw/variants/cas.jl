module C172FBWCAS

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
using Flight.FlightComponents.Piston
using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.World
using Flight.FlightComponents.Control: PIDDiscreteY

using ...C172
using ..C172FBW

export Cessna172FBWCAS

################################################################################
############################### PitchRateCmp ##################################

@kwdef struct PitchRateCmp <: SystemDefinition
    c1::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 0, k_i = 1, k_d = 0) #pure integrator
    c2::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 5.2, k_i = 25, k_d = 0.45, τ_d = 0.04, β_p = 1, β_d = 1) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct PitchRateCmpU
    setpoint::Float64 = 0.0
    feedback::Float64 = 0.0
    reset::Bool = false
    sat_ext::Int64 = 0.0
end

@kwdef struct PitchRateCmpY
    setpoint::Float64 = 0.0
    feedback::Float64 = 0.0
    reset::Bool = false
    sat_ext::Int64 = 0.0
    out::Float64 = 0.0 #elevator command
    c1::PIDDiscreteY{1} = PIDDiscreteY{1}()
    c2::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::PitchRateCmp) = PitchRateCmpU()
Systems.init(::SystemY, ::PitchRateCmp) = PitchRateCmpY()

#here, we are leaving the compensators' outputs unbounded, relying only on the
#external saturation input coming from the elevator to stop the integrators
function Systems.init!(sys::System{PitchRateCmp})
    sys.c1.u.reset .= true
    sys.c1.u.anti_windup .= true
    sys.c2.u.reset .= true
    sys.c2.u.anti_windup .= true
end

#if c2's output saturates positively, then both c2 and c1's integrators must
#stop accumulating positively, and viceversa.
function Systems.f_disc!(sys::System{PitchRateCmp}, Δt::Real)
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

    sys.y = PitchRateCmpY(; setpoint, feedback, reset, sat_ext, out, c1 = c1.y, c2 = c2.y)

end

function GUI.draw(pitch_rate_comp::System{<:PitchRateCmp})
    if CImGui.TreeNode("Integrator")
        GUI.draw(pitch_rate_comp.c1)
        CImGui.TreePop()
    end
    if CImGui.TreeNode("PID")
        GUI.draw(pitch_rate_comp.c2)
        CImGui.TreePop()
    end
end

################################################################################
################################### CAS ########################################

@enum RollMode begin
    aileron_mode = 0
    roll_rate_mode = 1
end

@enum PitchMode begin
    elevator_mode = 0
    pitch_rate_mode = 1
end

@enum YawMode begin
    rudder_mode = 0
    sideslip_mode = 1
end

################################################################################

@kwdef struct CASInner <: SystemDefinition
    p_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 0.8, k_i = 10, k_d = 0.05, τ_d = 0.04) #roll rate compensator, see notebook
    q_comp::PitchRateCmp = PitchRateCmp(
        PIDDiscrete{1}(k_p = 0, k_i = 1, k_d = 0),
        PIDDiscrete{1}(k_p = 4.8, k_i = 35, k_d = 0.6, τ_d = 0.04)) #see notebook
    β_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 8, k_i = 20, k_d = 5, τ_d = 0.04) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct CASInnerU
    roll_mode::RollMode = aileron_mode
    a_cmd::Float64 = 0.0
    p_cmd::Float64 = 0.0
    p::Float64 = 0.0
    pitch_mode::PitchMode = elevator_mode
    e_cmd::Float64 = 0.0
    q_cmd::Float64 = 0.0
    q::Float64 = 0.0
    yaw_mode::YawMode = sideslip_mode
    r_cmd::Float64 = 0.0
    β_cmd::Float64 = 0.0
    β::Float64 = 0.0
end

@kwdef mutable struct CASInnerS
    roll_mode_prev::RollMode = aileron_mode
    pitch_mode_prev::PitchMode = elevator_mode
    yaw_mode_prev::PitchMode = rudder_mode
end

@kwdef struct CASInnerY
    roll_mode::RollMode = aileron_mode
    roll_mode_prev::RollMode = aileron_mode
    a_cmd::Float64 = 0.0
    p_cmd::Float64 = 0.0
    a_out::Ranged{Float64, -1., 1.} = 0.0
    a_sat::Int64 = 0 #aileron saturation state
    p_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
    pitch_mode::PitchMode = elevator_mode
    pitch_mode_prev::PitchMode = elevator_mode
    e_cmd::Float64 = 0.0
    q_cmd::Float64 = 0.0
    e_out::Ranged{Float64, -1., 1.} = 0.0
    e_sat::Int64 = 0 #elevator saturation state
    q_comp::PitchRateCmpY = PitchRateCmpY()
    yaw_mode::YawMode = rudder_mode
    yaw_mode_prev::YawMode = rudder_mode
    r_cmd::Float64 = 0.0
    β_cmd::Float64 = 0.0
    r_out::Ranged{Float64, -1., 1.} = 0.0
    r_sat::Int64 = 0 #rudder saturation state
    β_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::CASInner) = CASInnerU()
Systems.init(::SystemY, ::CASInner) = CASInnerY()
Systems.init(::SystemS, ::CASInner) = CASInnerS()

function Systems.init!(sys::System{CASInner})
    Systems.init!(sys.q_comp)
    sys.p_comp.u.reset .= true
    sys.p_comp.u.anti_windup .= true
    sys.β_comp.u.reset .= true
    sys.β_comp.u.anti_windup .= true
end

function Systems.f_disc!(sys::System{PitchControl}, Δt::Real)

    @unpack roll_mode, a_cmd, p_cmd, p,
            pitch_mode, e_cmd, q_cmd, q,
            yaw_mode, r_cmd, β_cmd, β = sys.u
    @unpack roll_mode_prev, pitch_mode_prev, yaw_mode_prev = sys.s
    @unpack p_comp, q_comp, β_comp = sys.subsystems

    #manage mode transitions; since we only have two modes per axis, any mode
    #transition means the compensator should reset
    if roll_mode != roll_mode_prev
        p_comp.u.reset = true; f_disc!(q_comp, Δt)
    end

    if pitch_mode != pitch_mode_prev
        q_comp.u.reset = true; f_disc!(q_comp, Δt)
    end

    if yaw_mode != yaw_mode_prev
        β_comp.u.reset .= true; f_disc!(β_comp, Δt)
    end

    if roll_mode == aileron_mode
        a_out = Ranged(a_cmd, -1., 1.)
    else #roll rate compensator active
        p_comp.u.reset .= false
        p_comp.u.feedback .= p
        p_comp.u.setpoint .= p_cmd
        f_disc!(p_comp, Δt)
        a_out = Ranged(p_comp.y.out[1], -1., 1.)
    end

    if mode === elevator_mode
        e_out = Ranged(e_cmd, -1., 1.)
    else #pitch rate compensator active
        q_comp.u.reset = false
        q_comp.u.feedback = q
        q_comp.u.setpoint = q_cmd
        f_disc!(q_comp, Δt)
        e_out = Ranged(q_comp.y.out, -1., 1.)
    end

    if mode == rudder_mode
        r_out = Ranged(r_cmd, -1., 1.)
    else #sideslip compensator active
        β_comp.u.reset .= false
        β_comp.u.setpoint .= β_cmd
        β_comp.u.feedback .= β
        f_disc!(β_comp, Δt)
        r_out = Ranged(-β_comp.y.out[1], -1., 1.) #note sign inversion, see design notebook
    end

    #determine saturation states and assign them back to the compensators, will
    #take effect on the next call. since rudder output is inverted from
    #β_comp's output, we need to invert the saturation signal as well
    a_sat = (a_out == typemax(a_out)) - (a_out == typemin(a_out))
    e_sat = (e_out == typemax(e_out)) - (e_out == typemin(e_out))
    r_sat = (r_out == typemax(r_out)) - (r_out == typemin(r_out))

    p_comp.u.sat_ext .= a_sat
    q_comp.u.sat_ext = e_sat
    β_comp.u.sat_ext .= -r_sat

    #update FSM states
    roll_mode_prev = roll_mode
    pitch_mode_prev = pitch_mode
    yaw_mode_prev = yaw_mode
    @pack! sys.s = roll_mode_prev, pitch_mode_prev, yaw_mode_prev

    sys.y = PitchControlY(;
        roll_mode, roll_mode_prev, a_cmd, p_cmd, a_out, a_sat, p_comp = p_comp.y,
        pitch_mode, pitch_mode_prev, e_cmd, q_cmd, e_out, e_sat, q_comp = q_comp.y,
        yaw_mode, yaw_mode_prev, r_cmd, β_cmd, r_out, r_sat, β_comp = β_comp.y)


function GUI.draw(sys::System{<:CASInner})
    @unpack p_comp, q_comp, β_comp = sys.subsystems
    CImGui.Begin("CAS")
    if CImGui.TreeNode("Roll Rate Compensator")
        GUI.draw(p_comp)
        CImGui.TreePop()
    end
    if CImGui.TreeNode("Pitch Rate Compensator")
        GUI.draw(q_comp)
        CImGui.TreePop()
    end
    if CImGui.TreeNode("Sideslip Compensator")
        GUI.draw(β_comp)
        CImGui.TreePop()
    end
    CImGui.End()
end


##################################################################################
################################## Avionics ######################################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

#surface command offsets should be applied downstream of the CAS-computed
#surface commands to enable smooth transitions from a manually trimmed flight
#condition to a CAS control mode. if these offsets are modified with the CAS
#modes enabled, they will be handled as disturbances by the CAS and modify the
#computed surface commands to track the required demands

@kwdef mutable struct AvionicsInterfaceU
    eng_start::Bool = false
    eng_stop::Bool = false
    roll_mode_select::RollMode = aileron_mode #selected roll axis mode
    pitch_mode_select::PitchMode = elevator_mode #selected pitch axis mode
    yaw_mode_select::YawMode = rudder_mode #selected yaw axis mode
    throttle::Ranged{Float64, 0., 1.} = 0.0
    mixture::Ranged{Float64, 0., 1.} = 0.5
    aileron_cmd::Ranged{Float64, -1., 1.} = 0.0
    elevator_cmd::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd::Ranged{Float64, -1., 1.} = 0.0
    aileron_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    elevator_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    roll_rate_cmd::Ranged{Float64, -1., 1.} = 0.0 #normalized
    pitch_rate_cmd::Ranged{Float64, -1., 1.} = 0.0 #normalized
    sideslip_cmd::Ranged{Float64, -1., 1.} = 0.0 #normalized
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end

@kwdef struct AvionicsInterfaceY
    eng_start::Bool = false
    eng_stop::Bool = false
    roll_mode::RollMode = aileron_mode #actual roll axis mode
    pitch_mode::PitchMode = elevator_mode #actual pitch axis mode
    yaw_mode::YawMode = rudder_mode #actual yaw axis mode
    throttle::Float64 = 0.0
    mixture::Float64 = 0.5
    aileron_cmd::Float64 = 0.0
    elevator_cmd::Float64 = 0.0
    rudder_cmd::Float64 = 0.0
    aileron_cmd_offset::Float64 = 0.0
    elevator_cmd_offset::Float64 = 0.0
    rudder_cmd_offset::Float64 = 0.0
    roll_rate_cmd::Float64 = 0.0 #normalized
    pitch_rate_cmd::Float64 = 0.0 #normalized
    sideslip_cmd::Float64 = 0.0 #normalized
    flaps::Float64 = 0.0
    brake_left::Float64 = 0.0
    brake_right::Float64 = 0.0
end

@kwdef struct AvionicsLogicY
    flight_phase::FlightPhase = phase_gnd
end

@kwdef struct Avionics <: AbstractAvionics
    p_cmd_sf::Float64 = 0.2 #p_cmd scale factor
    q_cmd_sf::Float64 = 0.2 #q_cmd scale factor
    β_cmd_sf::Float64 = deg2rad(10) #β_cmd scale factor
    cas::CASInner = CASInner()
end

const AvionicsU = AvionicsInterfaceU

@kwdef struct AvionicsY
    interface::AvionicsInterfaceY = AvionicsInterfaceY()
    logic::AvionicsLogicY = AvionicsLogicY()
    cas::CASInnerY = CASInnerY()
end

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()
Systems.init(::SystemS, ::Avionics) = nothing #keep subsystems local


########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:Avionics}, Δt::Real,
                        airframe::System{<:C172.Airframe}, kinematics::KinematicData,
                        ::RigidBodyData, air::AirData, ::TerrainData)

    @unpack cas = avionics.subsystems
    # @unpack eng_start, eng_stop,
    #         roll_mode_select, pitch_mode_select, yaw_mode_select,
    #         throttle, mixture, aileron_cmd, elevator_cmd, rudder_cmd,
    #         aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
    #         flaps, brake_left, brake_right = avionics.u
    @unpack p_cmd_sf, q_cmd_sf, β_cmd_sf = avionics.params

    any_wow = any(SVector{3}(leg.strut.wow for leg in airframe.ldg.y))
    flight_phase = any_wow ? phase_gnd : phase_air

    if flight_phase == phase_gnd
        roll_mode = aileron_mode
        pitch_mode = elevator_mode
        yaw_mode = rudder_mode
    else
        # if ap.y.state === ap_enabled
        #     roll_mode = roll_rate_mode
        #     pitch_mode = pitch_rate_mode
        #     yaw_mode = sideslip_mode
        # else
            roll_mode =  roll_mode_select
            pitch_mode = pitch_mode_select
            yaw_mode = yaw_mode_select
        # end
    end

    # e_nb = REuler(kinematics.q_nb)
    # @unpack θ, φ = e_nb
    p, q, _ = kinematics.ω_lb_b
    β = air.β_b

    #the sideslip controller provides a positive β from a positive β_cmd input.
    #when in direct rudder command mode, a positive yaw input causes a positive
    #yaw rate. however, a positive β_cmd increment initially requires a negative
    #yaw rate. the sign inversion from the external sideslip_cmd to the CAS
    #β_cmd input keeps the consistency in the perceived behaviour between both
    #yaw control modes.
    a_cmd = Float64(aileron_cmd)
    e_cmd = Float64(elevator_cmd)
    r_cmd = Float64(rudder_cmd)

    p_cmd = p_cmd_sf * Float64(roll_rate_cmd)
    q_cmd = q_cmd_sf * Float64(pitch_rate_cmd)
    β_cmd = -β_cmd_sf * Float64(sideslip_cmd)

    @pack! cas.u = roll_mode, pitch_mode, yaw_mode, a_cmd, e_cmd, r_cmd, p_cmd, q_cmd, β_cmd, p, q, β
    f_disc!(cas, Δt)

    # interface_y = AvionicsInterfaceY(;
    #         eng_start, eng_stop, CAS_state,
    #         roll_mode, pitch_mode, yaw_mode,
    #         throttle, mixture, roll_input, pitch_input, yaw_input,
    #         aileron_offset, elevator_offset, rudder_offset,
    #         flaps, brake_left, brake_right)

    # logic_y = AvionicsLogicY(; flight_phase)

    # avionics.y = AvionicsY( interface = interface_y, logic = logic_y,
    #                         roll_control = roll_control.y,
    #                         pitch_control = pitch_control.y,
    #                         yaw_control = yaw_control.y)

    return false

end

# function Aircraft.map_controls!(airframe::System{<:C172.Airframe},
#                                 avionics::System{Avionics})

#     @unpack eng_start, eng_stop, throttle, mixture,
#             aileron_offset, elevator_offset, rudder_offset, flaps,
#             brake_left, brake_right = avionics.y.interface

#     u_act = airframe.act.u

#     u_act.aileron = avionics.y.roll_control.a_out
#     u_act.elevator = avionics.y.pitch_control.e_out
#     u_act.rudder = avionics.y.yaw_control.r_out

#     @pack!  u_act = eng_start, eng_stop, throttle, mixture,
#             aileron_offset, elevator_offset, rudder_offset, flaps,
#             brake_left, brake_right

# end


# ################################## GUI #########################################

# function control_mode_HSV(mode, selected_mode, active_mode)
#     if active_mode === mode
#         return HSV_green
#     elseif selected_mode === mode
#         return HSV_amber
#     else
#         return HSV_gray
#     end
# end

# function GUI.draw!(avionics::System{<:Avionics}, airframe::System{<:C172.Airframe},
#                     label::String = "Cessna 172R CAS Avionics")

#     u = avionics.u

#     CImGui.Begin(label)

#     CImGui.PushItemWidth(-60)

#     if airframe.y.pwp.engine.state === Piston.eng_off
#         eng_start_HSV = HSV_gray
#     elseif airframe.y.pwp.engine.state === Piston.eng_starting
#         eng_start_HSV = HSV_amber
#     else
#         eng_start_HSV = HSV_green
#     end
#     dynamic_button("Engine Start", eng_start_HSV, 0.1, 0.2)
#     u.eng_start = CImGui.IsItemActive()
#     CImGui.SameLine()
#     dynamic_button("Engine Stop", HSV_gray, (HSV_gray[1], HSV_gray[2], HSV_gray[3] + 0.1), (0.0, 0.8, 0.8))
#     u.eng_stop = CImGui.IsItemActive()
#     CImGui.SameLine()
#     CImGui.Text(@sprintf("%.3f RPM", Piston.radpersec2RPM(airframe.y.pwp.engine.ω)))
#     CImGui.Separator()

#     if avionics.y.interface.CAS_state === CAS_disabled
#         CAS_HSV = HSV_gray
#     elseif avionics.y.interface.CAS_state === CAS_standby
#         CAS_HSV = HSV_amber
#     else
#         CAS_HSV = HSV_green
#     end
#     dynamic_button("CAS", CAS_HSV, 0.1, 0.1)
#     CImGui.IsItemActivated() ? u.CAS_enable = !u.CAS_enable : nothing
#     CImGui.SameLine()
#     CImGui.Text("Flight Phase: $(avionics.y.logic.flight_phase)")

#     @unpack roll_mode, pitch_mode, yaw_mode = avionics.y.interface

#     CImGui.Text("Roll Control Mode: "); CImGui.SameLine()
#     dynamic_button("Aileron", control_mode_HSV(aileron_mode, u.roll_mode_select, roll_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.roll_mode_select = aileron_mode : nothing
#     dynamic_button("Roll Rate", control_mode_HSV(roll_rate_mode, u.roll_mode_select, roll_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.roll_mode_select = roll_rate_mode : nothing
#     dynamic_button("Roll Angle", control_mode_HSV(roll_angle_mode, u.roll_mode_select, roll_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.roll_mode_select = roll_angle_mode : nothing

#     CImGui.Separator()
#     CImGui.Text("Pitch Control Mode: "); CImGui.SameLine()
#     dynamic_button("Elevator", control_mode_HSV(elevator_mode, u.pitch_mode_select, pitch_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.pitch_mode_select = elevator_mode : nothing
#     dynamic_button("Pitch Rate", control_mode_HSV(pitch_rate_mode, u.pitch_mode_select, pitch_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.pitch_mode_select = pitch_rate_mode : nothing
#     dynamic_button("Pitch Angle", control_mode_HSV(pitch_angle_mode, u.pitch_mode_select, pitch_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.pitch_mode_select = pitch_angle_mode : nothing

#     CImGui.Separator()
#     CImGui.Text("Yaw Control Mode: "); CImGui.SameLine()
#     dynamic_button("Rudder", control_mode_HSV(rudder_mode, u.yaw_mode_select, yaw_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.yaw_mode_select = rudder_mode : nothing
#     dynamic_button("Sideslip", control_mode_HSV(sideslip_mode, u.yaw_mode_select, yaw_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.yaw_mode_select = sideslip_mode : nothing

#     CImGui.Separator()
#     u.throttle = safe_slider("Throttle", u.throttle, "%.6f")
#     u.roll_input = safe_slider("Roll Input", u.roll_input, "%.6f")
#     u.pitch_input = safe_slider("Pitch Input", u.pitch_input, "%.6f")
#     u.yaw_input = safe_slider("Yaw Input", u.yaw_input, "%.6f")
#     u.aileron_offset = safe_input("Aileron Offset", u.aileron_offset, 0.001, 0.1, "%.6f")
#     u.elevator_offset = safe_input("Elevator Offset", u.elevator_offset, 0.001, 0.1, "%.6f")
#     u.rudder_offset = safe_input("Rudder Offset", u.rudder_offset, 0.001, 0.1, "%.6f")
#     u.flaps = safe_slider("Flaps", u.flaps, "%.6f")
#     u.mixture = safe_slider("Mixture", u.mixture, "%.6f")
#     u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
#     u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")

#     #Internals
#     CImGui.Separator()
#     @unpack roll_control, pitch_control, yaw_control = avionics.subsystems

#     if CImGui.TreeNode("Internals")
#         show_roll_control = @cstatic check=false @c CImGui.Checkbox("Roll Control", &check); CImGui.SameLine()
#         show_roll_control && GUI.draw(roll_control)
#         show_pitch_control = @cstatic check=false @c CImGui.Checkbox("Pitch Control", &check); CImGui.SameLine()
#         show_pitch_control && GUI.draw(pitch_control)
#         show_yaw_control = @cstatic check=false @c CImGui.Checkbox("Yaw Control", &check); CImGui.SameLine()
#         show_yaw_control && GUI.draw(yaw_control)
#         CImGui.TreePop()
#     end


#     CImGui.PopItemWidth()

#     CImGui.End()

# end

# ################################################################################
# ############################# Cessna172RCAS #####################################

# #Cessna172R with control augmenting Avionics
# const Cessna172RCAS{K} = C172R.Template{K, Avionics} where {K}
# Cessna172RCAS(kinematics = LTF()) = C172R.Template(kinematics, Avionics())


# # ############################ Joystick Mappings #################################

# function IODevices.assign!(sys::System{<:Cessna172RCAS}, joystick::Joystick,
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

#     u.aileron_offset -= 0.01 * was_released(joystick, :dpad_left)
#     u.aileron_offset += 0.01 * was_released(joystick, :dpad_right)
#     u.elevator_offset += 0.01 * was_released(joystick, :dpad_down)
#     u.elevator_offset -= 0.01 * was_released(joystick, :dpad_up)

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

#     u.aileron_offset -= 2e-4 * is_pressed(joystick, :hat_left)
#     u.aileron_offset += 2e-4 * is_pressed(joystick, :hat_right)
#     u.elevator_offset += 2e-4 * is_pressed(joystick, :hat_down)
#     u.elevator_offset -= 2e-4 * is_pressed(joystick, :hat_up)

#     u.flaps += 0.3333 * was_released(joystick, :button_3)
#     u.flaps -= 0.3333 * was_released(joystick, :button_2)

# end

# function IODevices.assign!(sys::System{Avionics},
#                            joystick::GladiatorNXTEvo,
#                            ::DefaultMapping)

#     u = sys.u

#     u.throttle = get_axis_value(joystick, :throttle)
#     u.roll_input = get_axis_value(joystick, :stick_x) |> aileron_curve
#     u.pitch_input = get_axis_value(joystick, :stick_y) |> elevator_curve
#     u.yaw_input = get_axis_value(joystick, :stick_z) |> rudder_curve

#     u.brake_left = is_pressed(joystick, :red_trigger_half)
#     u.brake_right = is_pressed(joystick, :red_trigger_half)

#     u.aileron_offset -= 2e-4 * is_pressed(joystick, :A3_left)
#     u.aileron_offset += 2e-4 * is_pressed(joystick, :A3_right)
#     u.elevator_offset += 2e-4 * is_pressed(joystick, :A3_down)
#     u.elevator_offset -= 2e-4 * is_pressed(joystick, :A3_up)

#     if is_pressed(joystick, :A3_press)
#         u.aileron_offset = 0
#         u.elevator_offset = 0
#     end

#     u.flaps += 0.3333 * was_released(joystick, :switch_down)
#     u.flaps -= 0.3333 * was_released(joystick, :switch_up)

# end



end #module