module C172FBWMCS

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays, HDF5, Interpolations

using Flight.FlightCore
using Flight.FlightCore.Utils

using Flight.FlightPhysics
using Flight.FlightComponents
using Flight.FlightComponents.Control.Discrete: Integrator, IntegratorOutput, PID, PIDOutput, LQRTracker, LQRTrackerOutput

using ...C172
using ..C172FBW

export Cessna172FBWMCS


################################################################################
########################## AbstractControlMode #################################

#a discrete system implementing a specific longitudinal or lateral control mode
abstract type AbstractControlMode <: SystemDefinition end

################################################################################
############################# Longitudinal Control #############################

#state vector for all longitudinal controllers
@kwdef struct XLon <: FieldVector{10, Float64}
    q::Float64 = 0.0; θ::Float64 = 0.0; #pitch rate, pitch angle
    v_x::Float64 = 0.0; v_z::Float64 = 0.0; #aerodynamic velocity, body axes
    α_filt::Float64 = 0.0; ω_eng::Float64 = 0.0; #filtered AoS, engine speed
    thr_v::Float64 = 0.0; thr_p::Float64 = 0.0; #throttle actuator states
    ele_v::Float64 = 0.0; ele_p::Float64 = 0.0; #elevator actuator states
end

#control input vector for longitudinal controllers
@kwdef struct ULon{T} <: FieldVector{2, T}
    throttle_cmd::T = 0.0
    elevator_cmd::T = 0.0
end

#assemble state vector from aircraft physics
function XLon(physics::System{<:C172FBW.Physics})

    @unpack throttle_cmd, elevator_cmd = physics.airframe.act.u
    @unpack airframe, air, rigidbody, kinematics = physics.y
    @unpack pwp, aero, act = airframe
    @unpack e_nb, ω_eb_b = kinematics

    q = ω_eb_b[2]
    θ = e_nb.θ
    v_x, _, v_z = air.v_wOb_b
    ω_eng = pwp.engine.ω
    α_filt = aero.α_filt
    thr_v = act.throttle_act.vel
    thr_p = act.throttle_act.pos
    ele_v = act.elevator_act.vel
    ele_p = act.elevator_act.pos

    XLon(; q, θ, v_x, v_z, α_filt, ω_eng, thr_v, thr_p, ele_v, ele_p)

end

#command vector for throttle + elevator SAS mode
@kwdef struct ZLonThrEle <: FieldVector{2, Float64}
    throttle_cmd::Float64 = 0.0
    elevator_cmd::Float64 = 50.0
end

function ZLonThrEle(physics::System{<:C172FBW.Physics})
    @unpack act = physics.y.airframe
    throttle_cmd = act.throttle_act.cmd
    elevator_cmd = act.elevator_act.cmd
    ZLonThrEle(; throttle_cmd, elevator_cmd)
end

#command vector for EAS + climb rate mode
@kwdef struct ZLonEASClm <: FieldVector{2, Float64}
    EAS::Float64 = 50.0
    climb_rate::Float64 = 0.0
end

function ZLonEASClm(physics::System{<:C172FBW.Physics})
    EAS = physics.y.air.EAS
    climb_rate = -physics.y.kinematics.common.v_eOb_n[3]
    ZLonEASClm(; EAS, climb_rate)
end

#command vector for EAS + throttle mode
@kwdef struct ZLonEASThr <: FieldVector{2, Float64}
    EAS::Float64 = 50.0
    throttle_cmd::Float64 = 0.0
end

function ZLonEASThr(physics::System{<:C172FBW.Physics})
    EAS = physics.y.air.EAS
    throttle_cmd = physics.y.airframe.act.throttle_act.cmd
    ZLonEASThr(; EAS, throttle_cmd)
end

################################################################################
############################### Lateral Control ################################

#state vector for all lateral controllers
@kwdef struct XLat <: FieldVector{10, Float64}
    p::Float64 = 0.0; r::Float64 = 0.0; φ::Float64 = 0.0; #roll rate, yaw rate, bank angle
    v_x::Float64 = 0.0; v_y::Float64 = 0.0; β_filt::Float64 = 0.0; #aerodynamic velocity, body axes, filtered AoS
    ail_v::Float64 = 0.0; ail_p::Float64 = 0.0; #aileron actuator states
    rud_v::Float64 = 0.0; rud_p::Float64 = 0.0; #rudder actuator states
end

#control input vector for lateral controllers
@kwdef struct ULat{T} <: FieldVector{2, T}
    aileron_cmd::T = 0.0
    rudder_cmd::T = 0.0
end

#assemble state vector from aircraft physics
function XLat(physics::System{<:C172FBW.Physics})

    @unpack aileron_cmd, rudder_cmd = physics.airframe.act.u
    @unpack airframe, air, rigidbody, kinematics = physics.y
    @unpack aero, act = airframe
    @unpack e_nb, ω_eb_b = kinematics

    p, _, r = ω_eb_b
    φ = e_nb.φ
    v_x, v_y, _ = air.v_wOb_b
    β_filt = aero.β_filt
    ail_v = act.aileron_act.vel
    ail_p = act.aileron_act.pos
    rud_v = act.rudder_act.vel
    rud_p = act.rudder_act.pos

    XLat(; p, r, φ, v_x, v_y, β_filt, ail_v, ail_p, rud_v, rud_p)

end

#command vector for throttle + elevator SAS mode
@kwdef struct ZLatPhiBeta <: FieldVector{2, Float64}
    φ::Float64 = 0.0
    β::Float64 = 0.0
end

function ZLatPhiBeta(physics::System{<:C172FBW.Physics})
    φ = physics.y.kinematics.common.e_nb.φ
    β = physics.y.air.β_b
    ZLatPhiBeta(; φ, β)
end


################################################################################
################################## Avionics ####################################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@enum LonMode begin
    lon_SAS_off = 0
    lon_thr_ele = 1
    lon_thr_q = 2
    lon_EAS_q = 3
    lon_EAS_clm = 4
    lon_EAS_thr = 5
    lon_EAS_alt = 6
end

@enum LatMode begin
    lat_SAS_off = 0
    lat_p_β = 1 #SISO roll rate over φ + β
    lat_φ_β = 2 #MIMO
    lat_AP_χ_β = 4
end

@kwdef struct Avionics <: AbstractAvionics
    thr_ele::LQRTracker{10, 2, 2, 20, 4} = LQRTracker{10, 2, 2}()
    EAS_clm::LQRTracker{10, 2, 2, 20, 4} = LQRTracker{10, 2, 2}()
    EAS_thr::LQRTracker{10, 2, 2, 20, 4} = LQRTracker{10, 2, 2}()
    q2e_int::Integrator = Integrator()
    q2e_pid::PID = PID()
    v2t_pid::PID = PID()
end

@kwdef mutable struct Inceptors
    eng_start::Bool = false
    eng_stop::Bool = false
    mixture::Ranged{Float64, 0., 1.} = 0.5
    throttle_input::Ranged{Float64, 0., 1.} = 0.0 #sets throttle_sp
    roll_input::Ranged{Float64, -1., 1.} = 0.0 #sets aileron_sp or p_sp
    pitch_input::Ranged{Float64, -1., 1.} = 0.0 #sets elevator_sp or q_sp
    yaw_input::Ranged{Float64, -1., 1.} = 0.0 #sets rudder_sp or β_sp
    throttle_sp_offset::Ranged{Float64, 0., 1.} = 0.0 #for direct throttle only
    aileron_sp_offset::Ranged{Float64, -1., 1.} = 0.0 #for direct aileron only
    elevator_sp_offset::Ranged{Float64, -1., 1.} = 0.0 #for direct elevator only
    rudder_sp_offset::Ranged{Float64, -1., 1.} = 0.0 #for direct rudder only
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end


@kwdef mutable struct DigitalInputs
    lon_mode_sel::LonMode = lon_SAS_off #selected longitudinal control mode
    lat_mode_sel::LatMode = lat_SAS_off #selected lateral control mode
    p_sf::Float64 = 1.0 #roll input to roll rate scale factor (0.1)
    q_sf::Float64 = 1.0 #pitch input to pitch rate scale factor (0.1)
    β_sf::Float64 = 1.0 #yaw input to β_sp scale factor (0.1)
    EAS_sp::Float64 = 50.0 #equivalent airspeed setpoint
    clm_sp::Float64 = 0.0 #climb rate setpoint
#     φ_dmd::Float64 = 0.0 #bank angle demand
#     χ_dmd::Float64 = 0.0 #course angle demand
#     h_dmd::Float64 = 0.0 #altitude demand
#     h_ref::AltitudeRef = ellipsoidal #altitude reference
end

@kwdef struct AvionicsU
    inceptors::Inceptors = Inceptors()
    digital::DigitalInputs = DigitalInputs()
end

@kwdef mutable struct AvionicsS
    lon_mode_prev::LonMode = lon_SAS_off
    lat_mode_prev::LatMode = lat_SAS_off
end

@kwdef struct Moding
    lon_mode_sel::LonMode = lon_SAS_off
    lat_mode_sel::LatMode = lat_SAS_off
    lon_mode::LonMode = lon_SAS_off
    lat_mode::LatMode = lat_SAS_off
    EAS_sp::Float64 = 0.0
    clm_sp::Float64 = 0.0
end

@kwdef struct ActuationCommands
    eng_start::Bool = false
    eng_stop::Bool = false
    mixture::Ranged{Float64, 0., 1.} = 0.5
    throttle_cmd::Ranged{Float64, 0., 1.} = 0.0
    aileron_cmd::Ranged{Float64, -1., 1.} = 0.0
    elevator_cmd::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd::Ranged{Float64, -1., 1.} = 0.0
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end

@kwdef struct AvionicsY
    moding::Moding = Moding()
    actuation::ActuationCommands = ActuationCommands()
    qθ_int::IntegratorOutput = IntegratorOutput()
    θEAS_ctl::LQRTrackerOutput{10, 2, 2, 20, 4} = LQRTrackerOutput{10, 2, 2}()
end

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()
Systems.init(::SystemS, ::Avionics) = AvionicsS()

function Systems.init!(sys::System{<:Avionics})

    @unpack θEAS_ctl = sys.subsystems
    θEAS_ctl.u.bound_lo .= ULon(; throttle_cmd = 0, elevator_cmd = -1)
    θEAS_ctl.u.bound_hi .=  ULon(; throttle_cmd = 1, elevator_cmd = 1)

    #these will be moved to f_disc! for gain scheduling
    θEAS_ctl.u.x_trim .= XLon(; q = -8.283340369112567e-6, θ = 0.023553480271892083,
                    v_x = 52.472488387730245, v_z = 1.236138810758686,
                    α_filt = 0.023553489660139648, ω_eng = 248.097186418114,
                    thr_v = 0.0, thr_p = 0.6489929542635975,
                    ele_v = 0.0, ele_p = -0.2424771622783136)

    θEAS_ctl.u.u_trim .= ULon(; throttle_cmd = 0.6489929542635975, elevator_cmd = -0.2424771622783136)
    θEAS_ctl.u.z_trim .= ZθEAS(; θ = 0.023553480271892083, EAS = 50)

    θEAS_ctl.u.C_fbk .= @SMatrix[0.0324562  -0.56216   0.0767795   0.00818343  -0.0297093   0.000135349   0.000166313   0.00628894  -3.31639e-5  -0.00286731;
                                2.59501     9.94617  -0.269466   -0.117464     7.67692    -1.95015e-5   -0.000331639  -0.00602097   0.0294353    1.53726]
    θEAS_ctl.u.C_fwd .= @SMatrix[4.04034   0.108129;
                                9.50847  -0.327379]
    θEAS_ctl.u.C_int .= 0
    error("Assign all gains!!!")

end

function Control.Discrete.reset!(sys::System{<:C172FBWMCS.Avionics})
    Control.Discrete.reset!(sys.θEAS_ctl)
end

########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:C172FBWMCS.Avionics},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack subsystems, u, s = avionics

    @unpack qθ_int, θEAS_ctl = avionics.subsystems

    @unpack eng_start, eng_stop, mixture,
            throttle_input, roll_input, pitch_input, yaw_input,
            throttle_sp_offset, aileron_sp_offset, elevator_sp_offset, rudder_sp_offset,
            flaps, brake_left, brake_right = u.inceptors

    @unpack lon_mode_sel, lat_mode_sel, p_sf, q_sf, β_sf, EAS_sp, θ_sp = u.digital

    @unpack lon_mode_prev, lat_mode_prev = s

    @unpack airframe, air = physics.y

    kinematics = physics.y.kinematics.common
    p, q, r = kinematics.ω_lb_b
    @unpack θ, φ = kinematics.e_nb
    @unpack EAS = air

    φ_bnd = clamp(φ, -π/3, π/3)

    throttle_sp = throttle_input + throttle_sp_offset
    elevator_sp = pitch_input + elevator_sp_offset
    aileron_sp = roll_input + aileron_sp_offset
    rudder_sp = yaw_input + rudder_sp_offset

    p_sp = p_sf * Float64(roll_input)
    q_sp = q_sf * Float64(pitch_input)
    β_sp = β_sf * Float64(yaw_input)

    any_wow = any(SVector{3}(leg.strut.wow for leg in airframe.ldg))
    flight_phase = any_wow ? phase_gnd : phase_air

    if flight_phase === phase_gnd
        lon_mode = lon_SAS_off
        lat_mode = lat_SAS_off
    else #air
        lon_mode = lon_mode_sel
        lat_mode = lat_mode_sel
    end


    ########################## Longitudinal Modes ##############################

    #actuation commands from pilot inceptor set points
    if lon_mode === lon_SAS_off

        throttle_cmd = throttle_sp
        elevator_cmd = elevator_sp

    #actuation commands computed by SAS
    elseif lon_mode === lon_thr_ele || lon_mode === lon_thr_q || lon_mode === lon_EAS_q

        #elevator_sp overridden by q tracker
        if lon_mode === lon_thr_q || lon_mode === lon_EAS_q

            elevator_cmd_sat = ULon(thr_ele.y.out_sat).elevator_cmd

            q2e_int.u.input = q_sp - q
            q2e_int.u.sat_ext = elevator_cmd_sat
            f_disc!(q2e_int, Δt)

            q2e.u.input = q2e_int.y.output
            q2e.u.sat_ext = elevator_cmd_sat
            # assign!(q2e, q2e_lookup)
            f_disc!(q2e, Δt)
            elevator_sp = q2e.y.output

            #throttle_sp overridden by EAS tracker
            if lon_mode === lon_EAS_q

                throttle_cmd_sat = ULon(thr_ele.y.out_sat).throttle_cmd
                v2t.u.input = EAS_sp
                v2t.u.sat_ext = throttle_cmd_sat
                # assign!(v2t, v2t_lookup)
                f_disc!(v2t, Δt)
                throttle_sp = v2t.y.output

            end

        end

        #SAS has no integral action, so it doesn't need to be reset
        thr_ele.u.x .= XLon(physics) #state feedback
        thr_ele.u.z .= ZLonThrEle(physics) #command variable feedback
        thr_ele.u.z_sp .= ZLonThrEle(; throttle_sp, elevator_sp) #command variable setpoint
        # assign!(thr_ele, thr_ele_lookup)
        f_disc!(thr_ele, Δt)
        @unpack throttle_cmd, elevator_cmd = ULon(thr_ele.y.output)

    elseif lon_mode === lon_EAS_clm || lon_mode === lon_EAS_thr || lon_mode == lon_EAS_alt

        #altitude mode overrides the other two
        if lon_mode === lon_EAS_alt
            f_disc!(alt_ctl, kinematics, air, Δt)

            if bla
                lon_mode = lon_EAS_clm
            else
                lon_mode = lon_EAS_thr
                # throttle_sp =1...
            end
            #otherwise, lon_mode_sel is respected
        end

        if lon_mode == lon_EAS_clm

            #actuation commands computed by EAS/climb_rate LQRTracker
            if lon_mode != lon_mode_prev
                reset!(EAS_clm)
                # z_sp = ZLonEASClm(physics)
                # EAS_sp, throttle_sp = z_sp.EAS, z_sp.throttle_cmd
            end

            EAS_clm.u.x .= XLon(physics) #state feedback
            EAS_clm.u.z .= ZLonEASClm(physics) #command variable feedback
            EAS_clm.u.z_sp .= ZLonEASClm(; EAS = EAS_sp, climb_rate = clm_sp) #command variable setpoint
            # assign!(EAS_clm, EAS_clm_lookup)
            f_disc!(EAS_clm, Δt)
            @unpack throttle_cmd, elevator_cmd = ULon(EAS_clm.y.output)

        else #lon_mode == lon_EAS_thr

            if lon_mode != lon_mode_prev
                reset!(EAS_thr)
                # z_sp = ZLonEASThr(physics)
                # EAS_sp, throttle_sp = z_sp.EAS, z_sp.throttle_cmd
            end

            EAS_thr.u.x .= XLon(physics) #state feedback
            EAS_thr.u.z .= ZLonEASThr(physics) #command variable feedback
            EAS_thr.u.z_sp .= ZLonEASThr(; EAS = EAS_sp, throttle_cmd = throttle_sp) #command variable setpoint
            # assign!(EAS_thr, EAS_thr_lookup)
            f_disc!(EAS_thr, Δt)
            @unpack throttle_cmd, elevator_cmd = ULon(EAS_thr.y.output)

        end

    end

    ############################# Lateral Modes ################################

    if lat_mode === lat_direct
        aileron_cmd = aileron_sp
        rudder_cmd = rudder_sp
    else
        aileron_cmd = aileron_sp
        rudder_cmd = rudder_sp
    end

    ########################### State & Outputs ################################

    s.lon_mode_prev = lon_mode
    s.lat_mode_prev = lat_mode

    error("Assign all outputs!!!")
    moding = Moding(; lon_mode_sel, lat_mode_sel, lon_mode, lat_mode,
                    EAS_sp, clm_sp )

    actuation = ActuationCommands(; eng_start, eng_stop, mixture,
                throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd,
                flaps, brake_left, brake_right)

    avionics.y = AvionicsY(; moding, actuation,
                            qθ_int = qθ_int.y,
                            θEAS_ctl = θEAS_ctl.y)

    return false

end

function Aircraft.assign!(airframe::System{<:C172FBW.Airframe},
                          avionics::System{Avionics})

    @unpack eng_start, eng_stop, mixture, throttle_cmd, aileron_cmd,
            elevator_cmd, rudder_cmd, flaps, brake_left, brake_right = avionics.y.actuation

    @pack! airframe.act.u = eng_start, eng_stop, mixture, throttle_cmd, aileron_cmd,
           elevator_cmd, rudder_cmd, flaps, brake_left, brake_right

end



################################## GUI #########################################

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
        SameLine, NewLine, IsItemActive, Separator, Text, Checkbox, RadioButton

function mode_button_HSV(button_mode, selected_mode, active_mode)
    if active_mode === button_mode
        return HSV_green
    elseif selected_mode === button_mode
        return HSV_amber
    else
        return HSV_gray
    end
end

function GUI.draw!(avionics::System{<:C172FBWMCS.Avionics},
                    physics::System{<:C172FBW.Physics},
                    label::String = "Cessna 172 SS CAS Avionics")
end

function GUI.draw(moding::Moding)


end

################################################################################
############################# Cessna172FBWMCS ##################################

const Cessna172FBWMCS{K, T} = C172FBW.Template{K, T, C172FBWMCS.Avionics} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain}

function Cessna172FBWMCS(kinematics = LTF(), terrain = HorizontalTerrain())
    C172FBW.Template(kinematics, terrain, C172FBWMCS.Avionics())
end


##################################### Tools ####################################

function Aircraft.trim!(ac::System{<:Cessna172FBWMCS},
                        trim_params::C172.TrimParameters = C172.TrimParameters())

    result = trim!(ac.physics, trim_params)
    trim_state = result[2]

    #makes Avionics inputs consistent with the trim solution obtained for the
    #aircraft physics so the trim condition is preserved upon simulation start
    #when direct control modes are selected
    @unpack mixture, flaps = trim_params
    @unpack throttle, aileron, elevator, rudder = trim_state

    u = ac.avionics.u
    u.inceptors.throttle_input = 0
    u.inceptors.roll_input = 0
    u.inceptors.pitch_input = 0
    u.inceptors.yaw_input = 0
    u.inceptors.throttle_sp_offset = throttle
    u.inceptors.aileron_sp_offset = aileron
    u.inceptors.elevator_sp_offset = elevator
    u.inceptors.rudder_sp_offset = rudder
    u.inceptors.mixture = mixture
    u.inceptors.flaps = flaps

    u.digital.lon_mode_sel = lon_direct
    u.digital.lat_mode_sel = lat_direct

    #update avionics outputs
    f_disc!(ac.avionics, ac.physics, 1)

    return result

end

function Aircraft.linearize!(ac::System{<:Cessna172FBWMCS}, args...; kwargs...)
    linearize!(ac.physics, args...; kwargs...)
end


# ############################ Joystick Mappings #################################

function IODevices.assign!(sys::System{<:Cessna172FBWMCS}, joystick::Joystick,
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

    u = sys.u.inceptors

    u.roll_input = get_axis_value(joystick, :right_analog_x) |> aileron_curve
    u.pitch_input = get_axis_value(joystick, :right_analog_y) |> elevator_curve
    u.yaw_input = get_axis_value(joystick, :left_analog_x) |> rudder_curve
    u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

    u.aileron_sp_offset -= 0.01 * was_released(joystick, :dpad_left)
    u.aileron_sp_offset += 0.01 * was_released(joystick, :dpad_right)
    u.elevator_sp_offset += 0.01 * was_released(joystick, :dpad_down)
    u.elevator_sp_offset -= 0.01 * was_released(joystick, :dpad_up)

    u.throttle_input += 0.1 * was_released(joystick, :button_Y)
    u.throttle_input -= 0.1 * was_released(joystick, :button_A)

    u.flaps += 0.3333 * was_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

end

function IODevices.assign!(sys::System{Avionics},
                           joystick::T16000M,
                           ::DefaultMapping)

    u = sys.u.inceptors

    u.throttle_input = get_axis_value(joystick, :throttle)
    u.roll_input = get_axis_value(joystick, :stick_x) |> aileron_curve
    u.pitch_input = get_axis_value(joystick, :stick_y) |> elevator_curve
    u.yaw_input = get_axis_value(joystick, :stick_z) |> rudder_curve

    u.brake_left = is_pressed(joystick, :button_1)
    u.brake_right = is_pressed(joystick, :button_1)

    u.aileron_sp_offset -= 2e-4 * is_pressed(joystick, :hat_left)
    u.aileron_sp_offset += 2e-4 * is_pressed(joystick, :hat_right)
    u.elevator_sp_offset += 2e-4 * is_pressed(joystick, :hat_down)
    u.elevator_sp_offset -= 2e-4 * is_pressed(joystick, :hat_up)

    u.flaps += 0.3333 * was_released(joystick, :button_3)
    u.flaps -= 0.3333 * was_released(joystick, :button_2)

end

function IODevices.assign!(sys::System{Avionics},
                           joystick::GladiatorNXTEvo,
                           ::DefaultMapping)

    u = sys.u.inceptors

    u.throttle_input = get_axis_value(joystick, :throttle)
    u.roll_input = get_axis_value(joystick, :stick_x) |> aileron_curve
    u.pitch_input = get_axis_value(joystick, :stick_y) |> elevator_curve
    u.yaw_input = get_axis_value(joystick, :stick_z) |> rudder_curve

    u.brake_left = is_pressed(joystick, :red_trigger_half)
    u.brake_right = is_pressed(joystick, :red_trigger_half)

    u.aileron_sp_offset -= 2e-4 * is_pressed(joystick, :A3_left)
    u.aileron_sp_offset += 2e-4 * is_pressed(joystick, :A3_right)
    u.elevator_sp_offset += 2e-4 * is_pressed(joystick, :A3_down)
    u.elevator_sp_offset -= 2e-4 * is_pressed(joystick, :A3_up)

    if is_pressed(joystick, :A3_press)
        u.aileron_offset = 0
        u.elevator_offset = 0
    end

    u.flaps += 0.3333 * was_released(joystick, :switch_down)
    u.flaps -= 0.3333 * was_released(joystick, :switch_up)

end



end #module