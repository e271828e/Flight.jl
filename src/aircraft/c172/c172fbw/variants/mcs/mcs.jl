module C172FBWMCS

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays, HDF5, Interpolations

using Flight.FlightCore
using Flight.FlightCore.Utils

using Flight.FlightPhysics
using Flight.FlightComponents
using Flight.FlightComponents.Control.Discrete: Integrator, IntegratorOutput, LQRTracker, LQRTrackerOutput

using ...C172
using ..C172FBW

export Cessna172FBWMCS


################################################################################
########################## AbstractControlChannel ##############################

abstract type AbstractControlChannel <: SystemDefinition end

################################################################################
############################# Longitudinal Control #############################

@kwdef struct LonControl <: AbstractControlChannel end

@kwdef mutable struct LonControlU
#     mode::PitchMode = direct_elevator_mode
    e_dmd::Ranged{Float64, -1., 1.} = 0.0 #aileron actuation demand
    thr_dmd::Ranged{Float64, 0., 1.} = 0.0 #throttle actuation demand
end

@kwdef struct LonControlY
#     mode::PitchMode = direct_elevator_mode
#     θ_dmd::Float64 = 0.0
#     EAS_dmd::Float64 = 0.0
    e_cmd::Ranged{Float64, -1., 1.} = 0.0 #elevator actuation command
    thr_cmd::Ranged{Float64, 0., 1.} = 0.0
    e_sat::Int64 = 0 #elevator saturation state
    thr_sat::Int64 = 0 #throttle saturation state
end

Systems.init(::SystemU, ::LonControl) = LonControlU()
Systems.init(::SystemY, ::LonControl) = LonControlY()

function Systems.init!(sys::System{<:LonControl})
end

function Control.reset!(sys::System{<:LonControl})
end

function Systems.f_disc!(sys::System{<:LonControl},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    #these are the control inputs and states for the longitudinal control
    #design, and therefore the ones that we need to store in the trim condition.
    #for the inner control modes we remove h from the state vector,
    #because it worsens F's conditioning and we have no need for it until we
    #design the high-level altitude mode

    # u_labels = [:elevator_cmd, :throttle_cmd]
    # x_labels = [:q, :θ, :v_x, :v_z, :α_filt, :ω_eng, :ele_v, :ele_p, :thr_v, :thr_p]
    @unpack airframe, air = physics.y
    kinematics = physics.y.kinematics.common

    x_trim = ComponentVector()
    u_trim = ComponentVector()

    # x =
    Δx = x - x_trim
    Δu = -K * Δx
    u = Δu + u_trim + u_ext
    e_cmd = Ranged(u.elevator_cmd, -1., 1.)
    thr_cmd = Ranged(u.elevator_cmd, -1., 1.)

    e_sat = saturation(e_cmd)
    thr_sat = saturation(thr_cmd)

    sys.y = LonControlY(; e_cmd, thr_cmd, e_sat, thr_sat)

end

function GUI.draw(sys::System{<:LonControl})
end


################################################################################
############################# Longitudinal Control #############################

@kwdef struct LatControl <: AbstractControlChannel end

@kwdef mutable struct LatControlU
#     mode::PitchMode = direct_elevator_mode
    a_dmd::Ranged{Float64, -1., 1.} = 0.0 #aileron actuation demand
    r_dmd::Ranged{Float64, -1., 1.} = 0.0 #throttle actuation demand
end

@kwdef struct LatControlY
#     mode::PitchMode = direct_elevator_mode
#     θ_dmd::Float64 = 0.0
#     EAS_dmd::Float64 = 0.0
    a_cmd::Ranged{Float64, -1., 1.} = 0.0 #aileron actuation command
    r_cmd::Ranged{Float64, -1., 1.} = 0.0 #rudder actuation command
    a_sat::Int64 = 0 #aileron saturation state
    r_sat::Int64 = 0 #rudder saturation state
end

Systems.init(::SystemU, ::LatControl) = LatControlU()
Systems.init(::SystemY, ::LatControl) = LatControlY()

function Systems.init!(sys::System{<:LatControl})
end

function Control.reset!(sys::System{<:LatControl})
end

function Systems.f_disc!(sys::System{<:LatControl},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    # u_labels = [:elevator_cmd, :throttle_cmd]
    # x_labels = [:q, :θ, :v_x, :v_z, :α_filt, :ω_eng, :ele_v, :ele_p, :thr_v, :thr_p]

    # x_trim = ComponentVector()
    # u_trim = ComponentVector()

    # Δx = x - x_trim
    # Δu = -K * Δx
    # u = Δu + u_trim + u_ext
    # a_cmd = Ranged(u.aileron_cmd, -1., 1.)
    # r_cmd = Ranged(u.rudder_cmd, -1., 1.)

    # a_sat = saturation(a_cmd)
    # r_sat = saturation(r_cmd)

    # sys.y = LatControlY(; a_cmd, r_cmd, a_sat, r_sat)

end

function GUI.draw(sys::System{<:LatControl})
end
####################### TO SIMPLIFY FROM HERE ON DOWN ##########################

################################################################################
################################## Avionics ####################################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@kwdef struct Avionics <: AbstractAvionics
    lon_ctl::LonControl = LonControl()
    lat_ctl::LatControl = LatControl()
end

@kwdef mutable struct Inceptors
    eng_start::Bool = false
    eng_stop::Bool = false
    mixture::Ranged{Float64, 0., 1.} = 0.5
    throttle_input::Ranged{Float64, 0., 1.} = 0.0 #used in direct_throttle_mode
    roll_input::Ranged{Float64, -1., 1.} = 0.0 #used in aileron_mode and roll_rate_mode
    pitch_input::Ranged{Float64, -1., 1.} = 0.0 #used in direct_elevator_mode and pitch_rate_mode
    yaw_input::Ranged{Float64, -1., 1.} = 0.0 #used in rudder_mode and sideslip_mode
    aileron_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    elevator_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end

@kwdef mutable struct DigitalInputs
#     throttle_mode_sel::ThrottleMode = direct_throttle_mode #selected throttle channel mode
#     roll_mode_sel::RollMode = direct_aileron_mode #selected roll channel mode
#     pitch_mode_sel::PitchMode = direct_elevator_mode #selected pitch channel mode
#     yaw_mode_sel::YawMode = direct_rudder_mode #selected yaw channel mode
#     lon_mode_sel::LonMode = lon_mode_semi #selected longitudinal control mode
#     lat_mode_sel::LatMode = lat_mode_semi #selected lateral control mode
#     EAS_dmd::Float64 = 40.0 #equivalent airspeed demand
#     θ_dmd::Float64 = 0.0 #pitch angle demand
#     c_dmd::Float64 = 0.0 #climb rate demand
#     φ_dmd::Float64 = 0.0 #bank angle demand
#     χ_dmd::Float64 = 0.0 #course angle demand
#     h_dmd::Float64 = 0.0 #altitude demand
#     h_ref::AltitudeRef = ellipsoidal #altitude reference
#     p_dmd_sf::Float64 = 1.0 #roll_input to p_dmd scale factor (0.2)
#     q_dmd_sf::Float64 = 1.0 #pitch_input to q_dmd scale factor (0.2)
end

@kwdef struct AvionicsU
    inceptors::Inceptors = Inceptors()
    digital::DigitalInputs = DigitalInputs()
end

@kwdef struct AvionicsModing
    # flight_phase::FlightPhase = phase_gnd
    # throttle_mode::ThrottleMode = direct_throttle_mode
    # roll_mode::RollMode = direct_aileron_mode
    # pitch_mode::PitchMode = direct_elevator_mode
    # yaw_mode::YawMode = direct_rudder_mode
    # lon_mode::LonMode = lon_mode_semi
    # lat_mode::LatMode = lat_mode_semi
end

@kwdef struct ActuationCommands
    # eng_start::Bool = false
    # eng_stop::Bool = false
    # mixture::Ranged{Float64, 0., 1.} = 0.5
    # throttle_cmd::Ranged{Float64, 0., 1.} = 0.0
    # aileron_cmd::Ranged{Float64, -1., 1.} = 0.0
    # elevator_cmd::Ranged{Float64, -1., 1.} = 0.0
    # rudder_cmd::Ranged{Float64, -1., 1.} = 0.0
    # flaps::Ranged{Float64, 0., 1.} = 0.0
    # brake_left::Ranged{Float64, 0., 1.} = 0.0
    # brake_right::Ranged{Float64, 0., 1.} = 0.0
end

@kwdef struct AvionicsY
    moding::AvionicsModing = AvionicsModing()
    actuation::ActuationCommands = ActuationCommands()
end

# Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()
Systems.init(::SystemS, ::Avionics) = nothing #keep subsystems local


# ########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:C172FBWMCS.Avionics},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack eng_start, eng_stop, mixture, throttle_input,
            roll_input, pitch_input, yaw_input,
            aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
            flaps, brake_left, brake_right = avionics.u.inceptors

    # @unpack lon_mode_sel, lat_mode_sel, EAS_dmd, θ_dmd, c_dmd, φ_dmd, χ_dmd, h_dmd,
    #         h_ref, p_dmd_sf, q_dmd_sf = avionics.u.digital

    @unpack lon_ctl = avionics.subsystems

    @unpack airframe, air = physics.y
    kinematics = physics.y.kinematics.common

    #direct surface and inner loop demands always come from inceptors
    lon_ctl.u.thr_dmd = throttle_input
    lon_ctl.u.e_dmd = pitch_input + elevator_cmd_offset
    lat_ctl.u.a_dmd = roll_input + aileron_cmd_offset
    lat_ctl.u.r_dmd = yaw_input + rudder_cmd_offset

    # roll_ctl.u.p_dmd = p_dmd_sf * Float64(roll_input)
    # pitch_ctl.u.q_dmd = q_dmd_sf * Float64(pitch_input)

    #digital inputs, may be overridden by high level modes (like AltControl)
    # throttle_ctl.u.EAS_dmd = EAS_dmd
    # pitch_ctl.u.θ_dmd = θ_dmd
    # pitch_ctl.u.c_dmd = c_dmd
    # pitch_ctl.u.EAS_dmd = EAS_dmd
    # roll_ctl.u.φ_dmd = φ_dmd
    # roll_ctl.u.χ_dmd = χ_dmd
    # alt_ctl.u.h_dmd = h_dmd
    # alt_ctl.u.h_ref = h_ref

    any_wow = any(SVector{3}(leg.strut.wow for leg in airframe.ldg))
    flight_phase = any_wow ? phase_gnd : phase_air

    if flight_phase === phase_gnd

        # #these are irrelevant on ground, but must be defined in all paths
        # lon_mode = lon_mode_semi
        # lat_mode = lat_mode_semi

        # throttle_ctl.u.mode = direct_throttle_mode
        # roll_ctl.u.mode = direct_aileron_mode
        # pitch_ctl.u.mode = direct_elevator_mode
        # yaw_ctl.u.mode = direct_rudder_mode

    else #air

        # lon_mode = lon_mode_sel
        # lat_mode = lat_mode_sel

        # if lon_mode === lon_mode_semi

            # #prioritize airspeed control via throttle
            # if (throttle_mode_sel === EAS_throttle_mode) && (pitch_mode_sel === EAS_pitch_mode)
            #     throttle_ctl.u.mode = EAS_throttle_mode
            #     pitch_ctl.u.mode = pitch_rate_mode
            # else
            #     throttle_ctl.u.mode = throttle_mode_sel
            #     pitch_ctl.u.mode = pitch_mode_sel
            # end

        # else #lon_mode === lon_mode_alt

            #we may need to reset altcontrol on mode change if it uses an integrator
        #     f_disc!(alt_ctl, kinematics, air, Δt)

        #     throttle_ctl.u.mode = alt_ctl.y.throttle_mode
        #     throttle_ctl.u.thr_dmd = alt_ctl.y.thr_dmd

        #     pitch_ctl.u.mode = alt_ctl.y.pitch_mode
        #     pitch_ctl.u.c_dmd = alt_ctl.y.c_dmd

        # end

        # if lat_mode === lat_mode_semi

        #     roll_ctl.u.mode = roll_mode_sel
        #     yaw_ctl.u.mode = yaw_mode_sel

        # end

    end

    f_disc!(lon_ctl, physics, Δt)
    f_disc!(lat_ctl, physics, Δt)

    elevator_cmd = lon_ctl.y.e_cmd
    throttle_cmd = lon_ctl.y.thr_cmd
    aileron_cmd = lat_ctl.y.a_cmd
    rudder_cmd = lat_ctl.y.r_cmd

    moding = AvionicsModing(;
      )

    actuation = ActuationCommands(; eng_start, eng_stop, mixture,
                throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd,
                flaps, brake_left, brake_right)

    avionics.y = AvionicsY(; moding, actuation,
                            lon_ctl = throttle_ctl.y,
                            )

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

    @unpack airframe = physics
    @unpack throttle_ctl, roll_ctl, pitch_ctl, yaw_ctl = avionics.subsystems

    u_inc = avionics.u.inceptors
    u_dig = avionics.u.digital
    y_mod = avionics.y.moding
    y_act = avionics.y.actuation

    Begin(label)

    PushItemWidth(-60)

    show_inceptors = @cstatic check=false @c Checkbox("Inceptors", &check); SameLine()
    show_digital = @cstatic check=false @c Checkbox("Digital", &check); SameLine()
    show_internals = @cstatic check=false @c Checkbox("Internals", &check)

    if show_inceptors
        Separator()
        if airframe.y.pwp.engine.state === Piston.eng_off
            eng_start_HSV = HSV_gray
        elseif airframe.y.pwp.engine.state === Piston.eng_starting
            eng_start_HSV = HSV_amber
        else
            eng_start_HSV = HSV_green
        end
        dynamic_button("Engine Start", eng_start_HSV, 0.1, 0.2)
        u_inc.eng_start = IsItemActive()
        SameLine()
        dynamic_button("Engine Stop", HSV_gray, (HSV_gray[1], HSV_gray[2], HSV_gray[3] + 0.1), (0.0, 0.8, 0.8))
        u_inc.eng_stop = IsItemActive()
        SameLine()
        u_inc.mixture = safe_slider("Mixture", u_inc.mixture, "%.6f")
        # Text(@sprintf("%.3f RPM", Piston.radpersec2RPM(airframe.y.pwp.engine.ω)))
        Separator()
        u_inc.throttle = safe_slider("Throttle", u_inc.throttle, "%.6f")
        u_inc.roll_input = safe_slider("Roll Input", u_inc.roll_input, "%.6f")
        u_inc.pitch_input = safe_slider("Pitch Input", u_inc.pitch_input, "%.6f")
        u_inc.yaw_input = safe_slider("Yaw Input", u_inc.yaw_input, "%.6f")
        Separator()
        u_inc.aileron_cmd_offset = safe_input("Aileron Offset", u_inc.aileron_cmd_offset, 0.001, 0.1, "%.6f")
        u_inc.elevator_cmd_offset = safe_input("Elevator Offset", u_inc.elevator_cmd_offset, 0.001, 0.1, "%.6f")
        u_inc.rudder_cmd_offset = safe_input("Rudder Offset", u_inc.rudder_cmd_offset, 0.001, 0.1, "%.6f")
        u_inc.flaps = safe_slider("Flaps", u_inc.flaps, "%.6f")
        Separator()
        u_inc.brake_left = safe_slider("Left Brake", u_inc.brake_left, "%.6f")
        u_inc.brake_right = safe_slider("Right Brake", u_inc.brake_right, "%.6f")
    end

    if show_digital
        Separator()
        AlignTextToFramePadding()
        Text("Longitudinal Control Mode")
        SameLine()
        dynamic_button("Semi-Automatic", mode_button_HSV(lon_mode_semi, u_dig.lon_mode_sel, y_mod.lon_mode), 0.1, 0.1)
        IsItemActive() ? u_dig.lon_mode_sel = lon_mode_semi : nothing
        SameLine()
        dynamic_button("Automatic", mode_button_HSV(lon_mode_alt, u_dig.lon_mode_sel, y_mod.lon_mode), 0.1, 0.1)
        IsItemActive() ? u_dig.lon_mode_sel = lon_mode_alt : nothing

        AlignTextToFramePadding()
        Text("Throttle Control Mode")
        SameLine()
        dynamic_button("Direct", mode_button_HSV(direct_throttle_mode, u_dig.throttle_mode_sel, y_mod.throttle_mode), 0.1, 0.1)
        IsItemActive() ? u_dig.throttle_mode_sel = direct_throttle_mode : nothing
        SameLine()
        dynamic_button("EAS##Throttle", mode_button_HSV(EAS_throttle_mode, u_dig.throttle_mode_sel, y_mod.throttle_mode), 0.1, 0.1)
        IsItemActive() ? u_dig.throttle_mode_sel = EAS_throttle_mode : nothing

        AlignTextToFramePadding()
        Text("Pitch Control Mode")
        SameLine()
        foreach(("Elevator", "Pitch Rate", "Pitch Angle", "Climb Rate", "EAS##Pitch"),
                (direct_elevator_mode, pitch_rate_mode, pitch_angle_mode, climb_rate_mode, EAS_pitch_mode)) do label, mode
            dynamic_button(label, mode_button_HSV(mode, u_dig.pitch_mode_sel, y_mod.pitch_mode), 0.1, 0.1)
            IsItemActive() ? u_dig.pitch_mode_sel = mode : nothing
            SameLine()
        end
        NewLine()

        u_dig.q_dmd_sf = safe_input("Pitch Rate Sensitivity (s/deg)", rad2deg(u_dig.q_dmd_sf), 0.01, 1.0, "%.3f") |> deg2rad
        u_dig.θ_dmd = safe_input("Pitch Angle Demand (deg)", rad2deg(u_dig.θ_dmd), 0.01, 1.0, "%.3f") |> deg2rad
        u_dig.c_dmd = safe_input("Climb Rate Demand (m/s)", u_dig.c_dmd, 0.01, 1.0, "%.3f")
        u_dig.EAS_dmd = safe_input("EAS Demand (m/s)", u_dig.EAS_dmd, 0.1, 1.0, "%.3f")
        u_dig.h_dmd = safe_input("Altitude Demand (m)", u_dig.h_dmd, 0.1, 1.0, "%.3f")
        AlignTextToFramePadding()
        Text("Altitude Reference")
        SameLine()
        RadioButton("Ellipsoidal", u_dig.h_ref === ellipsoidal) ? u_dig.h_ref = ellipsoidal : nothing
        SameLine()
        RadioButton("Orthometric", u_dig.h_ref === orthometric) ? u_dig.h_ref = orthometric : nothing

        Separator()
        AlignTextToFramePadding()
        Text("Lateral Control")
        SameLine()
        dynamic_button("Semi-Automatic", mode_button_HSV(lat_mode_semi, u_dig.lat_mode_sel, y_mod.lat_mode), 0.1, 0.1)
        IsItemActive() ? u_dig.lat_mode_sel = lat_mode_semi : nothing

        AlignTextToFramePadding()
        Text("Roll Control Mode")
        SameLine()
        foreach(("Aileron", "Roll Rate", "Bank Angle", "Course Angle"),
                (direct_aileron_mode, roll_rate_mode, bank_angle_mode, course_angle_mode)) do label, mode
            dynamic_button(label, mode_button_HSV(mode, u_dig.roll_mode_sel, y_mod.roll_mode), 0.1, 0.1)
            IsItemActive() ? u_dig.roll_mode_sel = mode : nothing
            SameLine()
        end
        NewLine()

        AlignTextToFramePadding()
        Text("Yaw Control Mode")
        SameLine()
        dynamic_button("Rudder", mode_button_HSV(direct_rudder_mode, u_dig.yaw_mode_sel, y_mod.yaw_mode), 0.1, 0.1)
        IsItemActive() ? u_dig.yaw_mode_sel = direct_rudder_mode : nothing

        u_dig.p_dmd_sf = safe_input("Roll Rate Sensitivity (s/deg)", rad2deg(u_dig.p_dmd_sf), 0.1, 1.0, "%.3f") |> deg2rad
        u_dig.φ_dmd = safe_input("Bank Angle Demand (deg)", rad2deg(u_dig.φ_dmd), 0.1, 1.0, "%.3f") |> deg2rad
        u_dig.χ_dmd = safe_input("Course Angle Demand (deg)", rad2deg(u_dig.χ_dmd), 0.1, 1.0, "%.3f") |> deg2rad
    end

    if show_internals
        Begin("Internals")
        Separator()
        show_throttle_ctl = @cstatic check=false @c Checkbox("Throttle Control", &check); SameLine()
        show_roll_ctl = @cstatic check=false @c Checkbox("Roll Control", &check); SameLine()
        show_pitch_ctl = @cstatic check=false @c Checkbox("Pitch Control", &check); SameLine()
        show_yaw_ctl = @cstatic check=false @c Checkbox("Yaw Control", &check); SameLine()
        show_moding = @cstatic check=false @c Checkbox("Moding", &check); SameLine()
        # show_actuation = @cstatic check=false @c Checkbox("Actuation", &check); SameLine()
        show_throttle_ctl && GUI.draw(throttle_ctl)
        show_roll_ctl && GUI.draw(roll_ctl)
        show_pitch_ctl && GUI.draw(pitch_ctl)
        show_yaw_ctl && GUI.draw(yaw_ctl)
        show_moding && GUI.draw(y_mod)
        # show_actuation && GUI.draw(y_act)
        End()
    end


    PopItemWidth()

    End()

end

function GUI.draw(moding::AvionicsModing)

    @unpack flight_phase, throttle_mode, roll_mode, pitch_mode, yaw_mode, lon_mode, lat_mode = moding

    Begin("Moding")
    Text("Flight Phase: $flight_phase")
    Text("Throttle Mode: $throttle_mode")
    Text("Roll Mode: $roll_mode")
    Text("Pitch Mode: $pitch_mode")
    Text("Yaw Mode: $yaw_mode")
    Text("Longitudinal Mode: $lon_mode")
    Text("Lateral Mode: $lat_mode")

    CImGui.End()

end

################################################################################
############################# Cessna172FBWMCS ##################################

#Cessna172R with control augmenting Avionics
const Cessna172FBWMCS{K} = C172FBW.Template{K, Avionics} where {K <: AbstractKinematicDescriptor}
Cessna172FBWMCS(kinematics = LTF()) = C172FBW.Template(kinematics, Avionics())


##################################### Tools ####################################

function Aircraft.trim!(ac::System{<:Cessna172FBWMCS},
                        trim_params::C172.TrimParameters = C172.TrimParameters())

    result = trim!(ac.physics, trim_params)
    trim_state = result[2]

    #makes Avionics inputs consistent with the trim solution obtained for the
    #aircraft physics so the trim condition is preserved during simulation
    @unpack mixture, flaps = trim_params
    @unpack throttle, aileron, elevator, rudder = trim_state

    u = ac.avionics.u
    u.inceptors.throttle = throttle
    u.inceptors.roll_input = aileron
    u.inceptors.pitch_input = elevator
    u.inceptors.yaw_input = rudder
    u.inceptors.mixture = mixture
    u.inceptors.flaps = flaps

    u.digital.lon_mode_sel = C172FBWMCS.lon_mode_semi
    u.digital.lat_mode_sel = C172FBWMCS.lat_mode_semi
    u.digital.throttle_mode_sel = C172FBWMCS.direct_throttle_mode
    u.digital.roll_mode_sel = C172FBWMCS.direct_aileron_mode
    u.digital.pitch_mode_sel = C172FBWMCS.direct_elevator_mode
    u.digital.yaw_mode_sel = C172FBWMCS.direct_rudder_mode

    #update avionics outputs
    f_disc!(ac.avionics, 1, ac.physics)

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

    u.aileron_cmd_offset -= 0.01 * was_released(joystick, :dpad_left)
    u.aileron_cmd_offset += 0.01 * was_released(joystick, :dpad_right)
    u.elevator_cmd_offset += 0.01 * was_released(joystick, :dpad_down)
    u.elevator_cmd_offset -= 0.01 * was_released(joystick, :dpad_up)

    u.throttle += 0.1 * was_released(joystick, :button_Y)
    u.throttle -= 0.1 * was_released(joystick, :button_A)

    u.flaps += 0.3333 * was_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

end

function IODevices.assign!(sys::System{Avionics},
                           joystick::T16000M,
                           ::DefaultMapping)

    u = sys.u.inceptors

    u.throttle = get_axis_value(joystick, :throttle)
    u.roll_input = get_axis_value(joystick, :stick_x) |> aileron_curve
    u.pitch_input = get_axis_value(joystick, :stick_y) |> elevator_curve
    u.yaw_input = get_axis_value(joystick, :stick_z) |> rudder_curve

    u.brake_left = is_pressed(joystick, :button_1)
    u.brake_right = is_pressed(joystick, :button_1)

    u.aileron_cmd_offset -= 2e-4 * is_pressed(joystick, :hat_left)
    u.aileron_cmd_offset += 2e-4 * is_pressed(joystick, :hat_right)
    u.elevator_cmd_offset += 2e-4 * is_pressed(joystick, :hat_down)
    u.elevator_cmd_offset -= 2e-4 * is_pressed(joystick, :hat_up)

    u.flaps += 0.3333 * was_released(joystick, :button_3)
    u.flaps -= 0.3333 * was_released(joystick, :button_2)

end

function IODevices.assign!(sys::System{Avionics},
                           joystick::GladiatorNXTEvo,
                           ::DefaultMapping)

    u = sys.u.inceptors

    u.throttle = get_axis_value(joystick, :throttle)
    u.roll_input = get_axis_value(joystick, :stick_x) |> aileron_curve
    u.pitch_input = get_axis_value(joystick, :stick_y) |> elevator_curve
    u.yaw_input = get_axis_value(joystick, :stick_z) |> rudder_curve

    u.brake_left = is_pressed(joystick, :red_trigger_half)
    u.brake_right = is_pressed(joystick, :red_trigger_half)

    u.aileron_cmd_offset -= 2e-4 * is_pressed(joystick, :A3_left)
    u.aileron_cmd_offset += 2e-4 * is_pressed(joystick, :A3_right)
    u.elevator_cmd_offset += 2e-4 * is_pressed(joystick, :A3_down)
    u.elevator_cmd_offset -= 2e-4 * is_pressed(joystick, :A3_up)

    if is_pressed(joystick, :A3_press)
        u.aileron_offset = 0
        u.elevator_offset = 0
    end

    u.flaps += 0.3333 * was_released(joystick, :switch_down)
    u.flaps -= 0.3333 * was_released(joystick, :switch_up)

end



end #module