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
    c2::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 5.2, k_i = 25, k_d = 0.45, τ_d = 0.04) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct PitchRateCmpU
    setpoint::MVector{1,Float64} = zeros(MVector{1})
    feedback::MVector{1,Float64} = zeros(MVector{1})
    reset::MVector{1,Bool} = zeros(MVector{1, Bool})
    sat_ext::MVector{1,Int64} = zeros(MVector{1, Int64})
end

@kwdef struct PitchRateCmpY
    setpoint::SVector{1,Float64} = zeros(SVector{1})
    feedback::SVector{1,Float64} = zeros(SVector{1})
    reset::SVector{1,Bool} = zeros(SVector{1, Bool})
    sat_ext::SVector{1,Int64} = zeros(SVector{1, Int64})
    out::SVector{1,Float64} = zeros(SVector{1})
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

    out = c2.y.out

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
############################### PitchAxisControl ###############################

@enum PitchMode begin
    elevator_mode = 0
    pitch_rate_mode = 1
    pitch_angle_mode = 2
    climb_rate_mode = 3
    altitude_hold_mode = 4
    airspeed_mode = 5
end

@kwdef struct PitchControl <: SystemDefinition
    q_comp::PitchRateCmp = PitchRateCmp()
    θ_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 2.5, k_i = 1.7, k_d = 0.18, τ_d = 0.04) #replace design with pure pitch rate feedback
    c_comp::PIDDiscrete{1} = PIDDiscrete{1}() #NOT IMPLEMENTED
    h_comp::PIDDiscrete{1} = PIDDiscrete{1}()#NOT IMPLEMENTED
    v_comp::PIDDiscrete{1} = PIDDiscrete{1}()#NOT IMPLEMENTED
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct PitchControlU
    mode::PitchMode = elevator_mode
    e_ext::Float64 = 0.0
    q_ext::Float64 = 0.0
    θ_ext::Float64 = 0.0
    c_ext::Float64 = 0.0 #external climb rate command
    h_ext::Float64 = 0.0 #not available as external mode
    v_ext::Float64 = 0.0 #not available as external mode
    q::Float64 = 0.0 #pitch rate
    θ::Float64 = 0.0 #pitch angle
    c::Float64 = 0.0 #climb rate
    h::Float64 = 0.0 #altitude
    v::Float64 = 0.0 #airspeed
end

@kwdef struct PitchControlY
    mode::PitchMode = elevator_mode
    e_out::Ranged{Float64, -1., 1.} = 0.0 #elevator output
    e_sat::Int64 = 0 #elevator saturation state
    q_comp::PitchRateCmpY = PitchRateCmpY()
    θ_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
    c_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
    h_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
    v_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::PitchControl) = PitchControlU()
Systems.init(::SystemY, ::PitchControl) = PitchControlY()
Systems.init(::SystemS, ::PitchControl) = PitchControlS()

function Systems.init!(sys::System{PitchControl})
    Systems.init!(sys.q_comp)
    sys.θ_comp.u.reset .= true
    sys.θ_comp.u.anti_windup .= true
    sys.c_comp.u.reset .= true
    sys.c_comp.u.anti_windup .= true
    sys.h_comp.u.reset .= true
    sys.h_comp.u.anti_windup .= true
    sys.v_comp.u.reset .= true
    sys.v_comp.u.anti_windup .= true
end

function Systems.f_disc!(sys::System{PitchControl}, Δt::Real)

    @unpack mode, e_ext, q_ext, θ_ext, c_ext, h_ext, v_ext, q, θ, c, h, v = sys.u
    @unpack q_comp, θ_comp, c_comp, h_comp, v_comp = sys.subsystems

    v_comp.u.feedback .= v
    h_comp.u.feedback .= h
    c_comp.u.feedback .= c
    θ_comp.u.feedback .= θ
    q_comp.u.feedback .= q

    v_comp.u.reset .= (mode != airspeed_mode)
    h_comp.u.reset .= (mode != altitude_hold_mode)
    c_comp.u.reset .= (mode != altitude_hold_mode && mode != climb_rate_mode)
    θ_comp.u.reset .= (mode === pitch_rate_mode || mode === elevator_mode)
    q_comp.u.reset .= (mode === elevator_mode)

    v_comp.u.setpoint .= v_ext
    f_disc!(v_comp, Δt)
    θ_from_v = v_comp.y.out[1]

    h_comp.u.setpoint .= h_ext
    f_disc!(h_comp, Δt)
    c_from_h = h_comp.y.out[1]

    if mode === altitude_hold_mode
        c_comp.u.setpoint .= c_from_h
    elseif mode === climb_rate_mode
        c_comp.u.setpoint .= c_ext
    end
    f_disc!(c_comp, Δt)
    θ_from_c = c_comp.y.out[1]

    if (mode === altitude_hold_mode || mode === climb_rate_mode)
        θ_comp.u.setpoint .= θ_from_c
    elseif mode === airspeed_mode
        θ_comp.u.setpoint .= θ_from_v
    elseif mode === pitch_angle_mode
        θ_comp.u.setpoint .= θ_ext
    end
    f_disc!(θ_comp, Δt)
    q_from_θ = θ_comp.y.out[1]

    if (mode != pitch_rate_mode)
        q_comp.u.setpoint .= q_from_θ
    else
        q_comp.u.setpoint .= q_ext
    end
    f_disc!(q_comp, Δt)
    e_from_q = q_comp.y.out[1]

    if mode === elevator_mode
        e_out = Ranged(e_ext, -1., 1.)
    else
        e_out = Ranged(e_from_q, -1., 1.)
    end

    #determine elevator saturation state
    e_sat = (e_out == typemax(e_out)) - (e_out == typemin(e_out))

    #assign to compensators (will take effect on the next call)
    q_comp.u.sat_ext .= e_sat
    θ_comp.u.sat_ext .= e_sat
    c_comp.u.sat_ext .= e_sat
    h_comp.u.sat_ext .= e_sat
    v_comp.u.sat_ext .= e_sat

    sys.y = PitchControlY(; mode, e_out, e_sat, q_comp = q_comp.y,
                            θ_comp = θ_comp.y, c_comp = c_comp.y,
                            h_comp = h_comp.y, v_comp = v_comp.y)

end

function GUI.draw(pitch_control::System{<:PitchControl})
    CImGui.Begin("Pitch Control")
    if CImGui.TreeNode("Pitch Rate Compensator")
        GUI.draw(pitch_control.q_comp)
        CImGui.TreePop()
    end
    if CImGui.TreeNode("Pitch Angle Compensator")
        GUI.draw(pitch_control.θ_comp)
        CImGui.TreePop()
    end
    CImGui.End()
end


################################################################################
################################## ThrottleControl #############################

@enum ThrottleMode begin
    direct_throttle_mode = 0
    airspeed_throttle_mode = 1
end

##################################################################################
################################## Avionics ######################################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@kwdef mutable struct PhysicalControlsU
    eng_start::Bool = false
    eng_stop::Bool = false
    mixture::Ranged{Float64, 0., 1.} = 0.5
    throttle_input::Ranged{Float64, 0., 1.} = 0.0 #used in direct_throttle_mode
    roll_input::Ranged{Float64, -1., 1.} = 0.0 #used in aileron_mode and roll_rate_mode
    pitch_input::Ranged{Float64, -1., 1.} = 0.0 #used in elevator_mode and pitch_rate_mode
    yaw_input::Ranged{Float64, -1., 1.} = 0.0 #used in rudder_mode and sideslip_mode
    aileron_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    elevator_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end

#these can all be found as outputs in PhysicalControls, no need to duplicate
#them here
# @kwdef struct PhysicalControlsY
#     eng_start::Bool = false
#     eng_stop::Bool = false
#     mixture::Float64 = 0.5
#     throttle_input::Float64 = 0.0
#     roll_input::Float64 = 0.0
#     pitch_input::Float64 = 0.0
#     yaw_input::Float64 = 0.0
#     aileron_cmd_offset::Float64 = 0.0
#     elevator_cmd_offset::Float64 = 0.0
#     rudder_cmd_offset::Float64 = 0.0
#     flaps::Float64 = 0.0
#     brake_left::Float64 = 0.0
#     brake_right::Float64 = 0.0
# end

@kwdef mutable struct DigitalControlsU
    # throttle_mode_sel::ThrottleMode = direct_throttle_mode #selected throttle mode
    # roll_mode_sel::RollMode = aileron_mode #selected roll axis mode
    pitch_mode_sel::PitchMode = elevator_mode #selected pitch axis mode
    # yaw_mode_sel::YawMode = rudder_mode #selected yaw axis mode
    # vertical_mode_sel::VerticalMode = no_vertical_mode #altitude_control_mode
    # horizontal_mode_sel::LateralMode = no_horizontal_mode #track_angle_mode
    p_cmd_sf::Float64 = 0.2 #roll_input to p_cmd scale factor
    q_cmd_sf::Float64 = 0.2 #pitch_input to q_cmd scale factor
    β_cmd_sf::Float64 = -deg2rad(10) #yaw_input β_cmd scale factor, sign inverted
end

#the β control loop tracks the β_cmd input. a positive β_cmd increment initially
#produces a negative yaw rate. the sign inversion β_cmd_sf keeps consistency in
#the perceived behaviour between direct rudder and β control modes.

#these are directly reflected in AvionicsLogicY
# @kwdef struct DigitalControlsY
#     throttle_mode_sel::ThrottleMode = direct_throttle_mode #selected throttle mode
#     roll_mode_sel::RollMode = aileron_mode #selected roll axis mode
#     pitch_mode_sel::PitchMode = elevator_mode #selected pitch axis mode
#     yaw_mode_sel::YawMode = rudder_mode #selected yaw axis mode
#     altitude_mode_sel::Bool = false
#     track_mode_sel::Bool = false
# end

@kwdef struct AvionicsInternalsY
    flight_phase::FlightPhase = phase_gnd
    # throttle_mode::ThrottleMode = direct_throttle_mode
    # roll_mode::RollMode = aileron_mode #actual roll axis mode
    pitch_mode::PitchMode = elevator_mode #actual pitch axis mode
    # yaw_mode::YawMode = rudder_mode #actual yaw axis mode
    # vertical_mode::VerticalMode = no_vertical_mode
    # horizontal_mode::LateralMode = no_horizontal_mode
end

@kwdef struct Avionics <: AbstractAvionics
    # throttle_control::ThrottleControl = ThrottleControl()
    # roll_control::RollControl = RollControl()
    pitch_control::PitchControl = PitchControl()
    # yaw_control::YawControl = YawControl()
end

@kwdef struct AvionicsU
    physical::PhysicalControlsU = PhysicalControlsU()
    digital::DigitalControlsU = DigitalControlsU()
end

@kwdef struct AvionicsY
    internals::AvionicsInternalsY = AvionicsInternalsY()
    # throttle_control::ThrottleControlY = ThrottleControlY()
    # roll_control::RollControlY = RollControlY()
    pitch_control::PitchControlY = PitchControlY()
    # yaw_control::YawControlY = YawControlY()
end

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()
Systems.init(::SystemS, ::Avionics) = nothing #keep subsystems local


# ########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:Avionics}, Δt::Real,
                        airframe::System{<:C172.Airframe}, kinematics::KinematicData,
                        ::RigidBodyData, air::AirData, ::TerrainData)

    @unpack eng_start, eng_stop, throttle, mixture,
            roll_input, pitch_input, rudder_input,
            aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
            flaps, brake_left, brake_right = avionics.u.controls
    # @unpack roll_mode_select, pitch_mode_select, yaw_mode_select = avionics.u.logic
    @unpack p_cmd_sf, q_cmd_sf, β_cmd_sf = avionics.params

    @unpack θ, φ = REuler(kinematics.q_nb)
    p, q, _ = kinematics.ω_lb_b
    β = air.β_b

    any_wow = any(SVector{3}(leg.strut.wow for leg in airframe.ldg.y))
    flight_phase = any_wow ? phase_gnd : phase_air

    if flight_phase == phase_gnd
        # throttle_control.u.mode = direct_throttle_mode
        # roll_control.u.mode = aileron_mode
        pitch_control.u.mode = elevator_mode
        # yaw_control.u.mode = rudder_mode
        # vertical_control.u.mode = no_vertical_control
        # horizontal_control.u.mode = no_horizontal_control
    else
        # throttle_control.u.mode = throttle_mode_sel
        # roll_control.u.mode = roll_mode_sel
        pitch_control.u.mode = pitch_mode_sel
        # yaw_control.u.mode = yaw_mode_sel
        # vertical_control.u.mode = vertical_mode_sel
        # horizontal_control.u.mode = horizontal_control
    end

    #assign external commands to control channels
    pitch_control.u.e_ext = Float64(pitch_input)
    pitch_control.u.q_ext = q_cmd_sf * Float64(pitch_input)
    pitch_control.u.θ_ext = u.digital.θ_cmd
    pitch_control.u.c_ext = u.digital.c_cmd
    #h and v not externally available

    #high level vertical_control is called after the assignment of control
    #modes and commands from the external interfaces, it can override them, but
    #it should only do so if its active mode requires it

    #assign inputs to high level control
    # @pack! vertical_control.u = h
    # f_disc!(vertical_control, throttle_control.u, pitch_control.u, Δt)
    # f_disc!(horizontal_control, roll_control.u, yaw_control.u, Δt)

    #assign feedback to control channels
    # @pack! throttle_control.u = v
    # @pack! roll_control.u = p, φ
    @pack! pitch_control.u = q, θ, c, h, v
    # @pack! yaw_control.u = β

    # f_disc!(throttle_control, Δt)
    # f_disc!(roll_control, Δt)
    f_disc!(pitch_control, Δt)
    # f_disc!(yaw_control, Δt)

    avionics.y = AvionicsY( controls = controls_y, logic = logic_y, cas = cas.y)

    return false

end

# function Aircraft.map_controls!(airframe::System{<:C172.Airframe},
#                                 avionics::System{Avionics})

#     @unpack eng_start, eng_stop, throttle, mixture,
#             aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
#             flaps, brake_left, brake_right = avionics.y.interface

#     @unpack a_out, e_out, r_out = avionics.y.cas

#     airframe.act.u.aileron_cmd = a_out
#     airframe.act.u.elevator_cmd = e_out
#     airframe.act.u.rudder_cmd = r_out

#     @pack!  u_act = eng_start, eng_stop, throttle, mixture,
#             aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
#             flaps, brake_left, brake_right

#     #surface command offsets will be applied downstream of the CAS-computed
#     #surface commands. this enables smooth transitions from a manually trimmed
#     #flight condition to a CAS control mode. if these offsets are modified with
#     #the CAS modes enabled, they will be handled as disturbances by the CAS and
#     #modify the computed surface commands to track the required demands

# end


# # ################################## GUI #########################################

# # function control_mode_HSV(mode, selected_mode, active_mode)
# #     if active_mode === mode
# #         return HSV_green
# #     elseif selected_mode === mode
# #         return HSV_amber
# #     else
# #         return HSV_gray
# #     end
# # end

# function GUI.draw!(avionics::System{<:Avionics}, airframe::System{<:C172.Airframe},
#                     label::String = "Cessna 172 FBW CAS Avionics")

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