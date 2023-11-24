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
########################## AbstractControlMode #################################

#a discrete system implementing a specific longitudinal or lateral control mode
abstract type AbstractControlMode <: SystemDefinition end

################################################################################
############################# Longitudinal Control #############################

#state vector for all longitudinal controllers
@kwdef struct XLon <: FieldVector{10, Float64}
    q::Float64 = 0.0; θ::Float64 = 0.0; #pitch rate, pitch angle
    v_x::Float64 = 0.0; v_z::Float64 = 0.0; #aerodynamic velocity, body axes
    α_filt::Float64 = 0.0; ω_eng::Float64 = 0.0; #filtered airflow angles
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

    @unpack throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd = physics.airframe.act.u
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

################################################################################
################################## θEASCtl #####################################

#command vector for θ/EAS controller
@kwdef struct ZθEAS <: FieldVector{2, Float64}
    θ::Float64 = 0.0
    EAS::Float64 = 50.0
end

#assemble command vector from aircraft physics
function ZθEAS(physics::System{<:C172FBW.Physics})
    @unpack air, kinematics = physics.y
    θ = kinematics.e_nb.θ
    EAS = air.EAS
    ZθEAS(; θ, EAS)
end

###############################

@kwdef struct θEASCtl <: AbstractControlMode
    lqr::LQRTracker{10, 2, 2, 20, 4} = LQRTracker{10, 2, 2}()
end

@kwdef struct θEASCtlOutputs
    u::ULon{Float64} = ULon(0.0, 0.0)
    u_sat::ULon{Int64} = ULon(0, 0)
    lqr::LQRTrackerOutput{10, 2, 2, 20, 4} = LQRTrackerOutput{10, 2, 2}()
end

Systems.init(::SystemU, ::θEASCtl) = Ref(ZθEAS())
Systems.init(::SystemY, ::θEASCtl) = θEASCtlOutputs()

function Systems.init!(sys::System{<:θEASCtl})
    #let the LQRTracker handle saturation
    @unpack bound_lo, bound_hi = sys.lqr.u
    bound_lo .= ULon(; throttle_cmd = 0, elevator_cmd = -1)
    bound_hi .=  ULon(; throttle_cmd = 1, elevator_cmd = 1)
end

Control.Discrete.reset!(sys::System{<:θEASCtl}) = Control.Discrete.reset!(sys.lqr)

#smooth reset, assigns the current values of the command variables as inputs
function Control.Discrete.reset!(sys::System{<:θEASCtl}, physics::System{<:C172FBW.Physics})
    reset!(sys)
    sys.u[] = ZθEAS(physics)
end

function Systems.f_disc!(sys::System{<:θEASCtl},
                         physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack lqr = sys
    # @unpack x_trim, u_trim, z_trim, C_fbk, C_fwd, C_int, z_sp, z, x = lqr.u

    #Lookup will have to be parametric with parameters XTrim, CFBK, CFWD, etc.
    #We will do:
    # x_trim .= lookup.x_trim(EAS, h_o)
    # u_trim .= lookup.u_trim(EAS, h_o)

    #create a function assign!(lqr, lookup, args...) that internally does
    #lqr.x_trim .= lookup.x_trim(args...), etc. try different options, see how
    #far it can be streamlined without allocating

    #these will eventually come from lookup tables
    lqr.u.x_trim .= XLon(; q = -8.283340369112567e-6, θ = 0.023553480271892083,
                    v_x = 52.472488387730245, v_z = 1.236138810758686,
                    α_filt = 0.023553489660139648, ω_eng = 248.097186418114,
                    thr_v = 0.0, thr_p = 0.6489929542635975,
                    ele_v = 0.0, ele_p = -0.2424771622783136)

    lqr.u.u_trim .= ULon(; throttle_cmd = 0.6489929542635975, elevator_cmd = -0.2424771622783136)
    lqr.u.z_trim .= ZθEAS(; θ = 0.023553480271892083, EAS = 50)

    lqr.u.C_fbk .= @SMatrix[0.0324562  -0.56216   0.0767795   0.00818343  -0.0297093   0.000135349   0.000166313   0.00628894  -3.31639e-5  -0.00286731;
                  2.59501     9.94617  -0.269466   -0.117464     7.67692    -1.95015e-5   -0.000331639  -0.00602097   0.0294353    1.53726]
    lqr.u.C_fwd .= @SMatrix[4.04034   0.108129;
                  9.50847  -0.327379]
    lqr.u.C_int .= 0

    # lqr.u.z_sp .= sys.u[]
    lqr.u.z_sp .= lqr.u.z_trim
    lqr.u.z .= ZθEAS(physics)
    lqr.u.x .= XLon(physics)

    f_disc!(lqr, Δt)

    u = ULon(lqr.y.output)
    u_sat = ULon(lqr.y.out_sat)

    sys.y = θEASCtlOutputs(; u, u_sat, lqr = lqr.y)

end

function GUI.draw(sys::System{<:θEASCtl})
end


################################################################################
################################## Avionics ####################################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@enum LonMode begin
    lon_direct = 0 #throttle + elevator
    lon_q_EAS = 1 #pitch rate + EAS
    lon_θ_EAS = 2 #pitch angle + EAS
    lon_c_EAS = 3 #climb rate + EAS
    lon_t_EAS = 4 #throttle + EAS
    lon_h_EAS = 5 #altitude + EAS
end

@enum LatMode begin
    lat_direct = 0 #aileron + rudder
    lat_p_β = 1 #roll rate + sideslip
    lat_φ_β = 2 #bank angle + sideslip
    lat_χ_β = 3 #course angle + sideslip
end

@kwdef struct Avionics <: AbstractAvionics
    θEAS_ctl::θEASCtl = θEASCtl()
end

@kwdef mutable struct Inceptors
    eng_start::Bool = false
    eng_stop::Bool = false
    mixture::Ranged{Float64, 0., 1.} = 0.5
    throttle_input::Ranged{Float64, 0., 1.} = 0.0 #used in direct_throttle_mode
    roll_input::Ranged{Float64, -1., 1.} = 0.0 #used in aileron_mode and roll_rate_mode
    pitch_input::Ranged{Float64, -1., 1.} = 0.0 #used in direct_elevator_mode and pitch_rate_mode
    yaw_input::Ranged{Float64, -1., 1.} = 0.0 #used in rudder_mode and sideslip_mode
    aileron_sp_offset::Ranged{Float64, -1., 1.} = 0.0
    elevator_sp_offset::Ranged{Float64, -1., 1.} = 0.0
    rudder_sp_offset::Ranged{Float64, -1., 1.} = 0.0
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end


@kwdef mutable struct DigitalInputs
    lon_mode_sel::LonMode = lon_mode_semi #selected longitudinal control mode
    lat_mode_sel::LatMode = lat_mode_semi #selected lateral control mode
    p_sf::Float64 = 1.0 #roll_input to p_sp scale factor (0.2)
    q_sf::Float64 = 1.0 #pitch_input to q_sp scale factor (0.2)
    β_sf::Float64 = 1.0 #yaw input to β_sp scale factor (0.1)
    EAS_sp::Float64 = 40.0 #equivalent airspeed setpoint
    θ_sp::Float64 = 0.0 #pitch angle setpoint
#     c_dmd::Float64 = 0.0 #climb rate demand
#     φ_dmd::Float64 = 0.0 #bank angle demand
#     χ_dmd::Float64 = 0.0 #course angle demand
#     h_dmd::Float64 = 0.0 #altitude demand
#     h_ref::AltitudeRef = ellipsoidal #altitude reference
end

@kwdef struct AvionicsU
    inceptors::Inceptors = Inceptors()
    digital::DigitalInputs = DigitalInputs()
end

@kwdef struct AvionicsModing
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
    moding::AvionicsModing = AvionicsModing()
    actuation::ActuationCommands = ActuationCommands()
end

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()
Systems.init(::SystemS, ::Avionics) = nothing #keep subsystems local


########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:C172FBWMCS.Avionics},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack θEAS_ctl = avionics.subsystems

    @unpack eng_start, eng_stop, mixture, throttle_input,
            roll_input, pitch_input, yaw_input,
            aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
            flaps, brake_left, brake_right = avionics.u.inceptors

    @unpack lon_mode_sel, lat_mode_sel, p_sf, q_sf, β_sf, EAS_sp, θ_sp = avionics.u.digital

    @unpack airframe, air = physics.y
    kinematics = physics.y.kinematics.common

    #direct surface and inner loop demands always come from inceptors
    throttle_sp = throttle_input
    elevator_sp = pitch_input + elevator_sp_offset
    aileron_sp = roll_input + aileron_sp_offset
    rudder_sp = yaw_input + rudder_sp_offset

    p_sp = p_sf * Float64(roll_input)
    q_sp = q_sf * Float64(pitch_input)
    β_sp = β_sf * Float64(yaw_input)

    any_wow = any(SVector{3}(leg.strut.wow for leg in airframe.ldg))
    flight_phase = any_wow ? phase_gnd : phase_air

    if flight_phase === phase_gnd

        # #these are irrelevant on ground, but must be defined in all paths
        lon_mode = lon_direct
        lat_mode = lat_direct

    else #air

        lon_mode = lon_mode_sel
        lat_mode = lat_mode_sel

    end

    #we don't write to our own u. this should be a rule. we write to our y, and
    #the GUI can read it and decide whether it respects it and reassigns it to u
    #on the next call, or overwrites with a different value

    if lon_mode === lon_direct
        throttle_cmd = throttle_sp
        elevator_cmd = elevator_sp
    else #lon_mode === lon_θ_EAS
        if lon_mode != lon_mode_prev
            reset!(θEAS_ctl, physics)
            @unpack θ_sp, EAS_sp = ZθEAS(physics)
        end
        θEAS_ctl.u[] = ZθEAS(θ = θ_sp, EAS = EAS_sp)
        f_disc!(θEAS_ctl, Δt)
        @unpack throttle_cmd, elevator_cmd = θEAS_ctl.u
    end

    # if lat_mode === lat_direct
        aileron_cmd = aileron_sp
        rudder_cmd = rudder_sp
    # else
    # end

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
end

function GUI.draw(moding::AvionicsModing)


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