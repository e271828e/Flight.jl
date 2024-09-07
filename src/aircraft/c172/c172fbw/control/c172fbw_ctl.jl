module FlightControl

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays

using Flight.FlightCore
using Flight.FlightPhysics
using Flight.FlightComponents
using Flight.FlightComponents.Control.Discrete: Integrator, IntegratorOutput,
    PID, PIDOutput, PIDParams, LQRTracker, LQRTrackerOutput, LQRTrackerParams,
    PIDLookup, LQRTrackerLookup, load_pid_lookup, load_lqr_tracker_lookup

using ...AircraftBase
using ...C172
using ..C172FBW


################################################################################
####################### AbstractControlChannel #################################

#a discrete system implementing a specific longitudinal or lateral control mode
abstract type AbstractControlChannel <: SystemDefinition end


################################################################################
################################## LonControl ##################################

@enum LonControlMode begin
    lon_direct = 0
    lon_thr_ele = 1
    lon_thr_q = 2
    lon_thr_θ = 3
    lon_thr_EAS = 4
    lon_EAS_q = 5
    lon_EAS_θ = 6
    lon_EAS_clm = 7
end

function te2te_enabled(mode::LonControlMode)
    mode === lon_thr_ele || mode === lon_thr_q || mode === lon_thr_θ ||
    mode === lon_thr_EAS || mode === lon_EAS_q || mode === lon_EAS_θ
end

function q2e_enabled(mode::LonControlMode)
    mode === lon_thr_q || mode === lon_thr_θ || mode === lon_thr_EAS ||
    mode === lon_EAS_q || mode === lon_EAS_θ
end

function θ2q_enabled(mode::LonControlMode)
    mode === lon_thr_θ || mode === lon_thr_EAS || mode === lon_EAS_θ
end

v2θ_enabled(mode::LonControlMode) = (mode === lon_thr_EAS)
v2t_enabled(mode::LonControlMode) = (mode === lon_EAS_q || mode === lon_EAS_θ)

############################## FieldVectors ####################################

#state vector for all longitudinal controllers
@kwdef struct XLon <: FieldVector{10, Float64}
    q::Float64 = 0.0; θ::Float64 = 0.0; #pitch rate, pitch angle
    v_x::Float64 = 0.0; v_z::Float64 = 0.0; #aerodynamic velocity, body axes
    α_filt::Float64 = 0.0; ω_eng::Float64 = 0.0; #filtered AoS, engine speed
    thr_v::Float64 = 0.0; thr_p::Float64 = 0.0; #throttle actuator states
    ele_v::Float64 = 0.0; ele_p::Float64 = 0.0; #elevator actuator states
end

#assemble state vector from vehicle
function XLon(vehicle::System{<:C172FBW.Vehicle})

    @unpack components, air, kinematics = vehicle.y
    @unpack pwp, aero, act = components
    @unpack e_nb, ω_eb_b = kinematics

    q = ω_eb_b[2]
    θ = e_nb.θ
    v_x, _, v_z = air.v_wOb_b
    ω_eng = pwp.engine.ω
    α_filt = aero.α_filt
    thr_v = act.throttle.vel
    thr_p = act.throttle.pos
    ele_v = act.elevator.vel
    ele_p = act.elevator.pos

    XLon(; q, θ, v_x, v_z, α_filt, ω_eng, thr_v, thr_p, ele_v, ele_p)

end

#control input vector for longitudinal controllers
@kwdef struct ULon{T} <: FieldVector{2, T}
    throttle_cmd::T = 0.0
    elevator_cmd::T = 0.0
end

function ULon{Float64}(vehicle::System{<:C172FBW.Vehicle})
    @unpack act = vehicle.y.components
    throttle_cmd = Float64(act.throttle.cmd)
    elevator_cmd = Float64(act.elevator.cmd)
    ULon(; throttle_cmd, elevator_cmd)
end

#command vector for throttle + elevator SAS mode
const ZLonThrEle = ULon{Float64}

#command vector for EAS + climb rate mode
@kwdef struct ZLonEASClm <: FieldVector{2, Float64}
    EAS::Float64 = C172.TrimParameters().EAS
    climb_rate::Float64 = 0.0
end

function ZLonEASClm(vehicle::System{<:C172FBW.Vehicle})
    EAS = vehicle.y.air.EAS
    climb_rate = -vehicle.y.kinematics.data.v_eOb_n[3]
    ZLonEASClm(; EAS, climb_rate)
end


################################## System ######################################

#since all LQRTrackers have the same dimensions, all LQRTrackerLookup parametric
#types should be the same. to be confirmed. same with PIDLookup types.
@kwdef struct LonControl{LQ <: LQRTrackerLookup, LP <: PIDLookup} <: AbstractControlChannel
    te2te_lookup::LQ = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "te2te_lookup.h5"))
    vc2te_lookup::LQ = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "vc2te_lookup.h5"))
    q2e_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "q2e_lookup.h5"))
    v2θ_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "v2θ_lookup.h5"))
    v2t_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "v2t_lookup.h5"))
    te2te_lqr::LQRTracker{10, 2, 2, 20, 4} = LQRTracker{10, 2, 2}()
    vc2te_lqr::LQRTracker{10, 2, 2, 20, 4} = LQRTracker{10, 2, 2}()
    q2e_int::Integrator = Integrator()
    q2e_pid::PID = PID()
    v2θ_pid::PID = PID()
    v2t_pid::PID = PID()
end

@kwdef mutable struct LonControlU
    mode::LonControlMode = lon_direct #selected control mode
    throttle_sp::Float64 = 0.0
    elevator_sp::Float64 = 0.0
    q_sp::Float64 = 0.0
    θ_sp::Float64 = 0.0
    EAS_sp::Float64 = C172.TrimParameters().EAS #equivalent airspeed setpoint
    clm_sp::Float64 = 0.0 #climb rate setpoint
end

@kwdef struct LonControlY
    mode::LonControlMode = lon_direct
    throttle_sp::Float64 = 0.0
    elevator_sp::Float64 = 0.0
    q_sp::Float64 = 0.0
    θ_sp::Float64 = 0.0
    EAS_sp::Float64 = C172.TrimParameters().EAS #equivalent airspeed setpoint
    clm_sp::Float64 = 0.0 #climb rate setpoint
    throttle_cmd::Ranged{Float64, 0., 1.} = 0.0
    elevator_cmd::Ranged{Float64, -1., 1.} = 0.0
    te2te_lqr::LQRTrackerOutput{10, 2, 2, 20, 4} = LQRTrackerOutput{10, 2, 2}()
    vc2te_lqr::LQRTrackerOutput{10, 2, 2, 20, 4} = LQRTrackerOutput{10, 2, 2}()
    q2e_int::IntegratorOutput = IntegratorOutput()
    q2e_pid::PIDOutput = PIDOutput()
    v2θ_pid::PIDOutput = PIDOutput()
    v2t_pid::PIDOutput = PIDOutput()
end

Systems.U(::LonControl) = LonControlU()
Systems.Y(::LonControl) = LonControlY()

function Systems.init!(sys::System{<:LonControl})
    #set output bounds in all LQR Trackers so they can handle saturation
    foreach((sys.te2te_lqr, sys.vc2te_lqr)) do lqr
        lqr.u.bound_lo .= ULon(; throttle_cmd = 0, elevator_cmd = -1)
        lqr.u.bound_hi .= ULon(; throttle_cmd = 1, elevator_cmd = 1)
    end

end


function Systems.f_disc!(::NoScheduling, sys::System{<:LonControl},
                            vehicle::System{<:C172FBW.Vehicle})

    @unpack mode, throttle_sp, elevator_sp, q_sp, θ_sp, EAS_sp, clm_sp = sys.u
    @unpack te2te_lqr, vc2te_lqr, q2e_int, q2e_pid, v2θ_pid, v2t_pid = sys.subsystems
    @unpack te2te_lookup, vc2te_lookup, q2e_lookup, v2θ_lookup, v2t_lookup = sys.constants

    EAS = vehicle.y.air.EAS
    h_e = Float64(vehicle.y.kinematics.h_e)
    _, q, r = vehicle.y.kinematics.ω_lb_b
    @unpack θ, φ = vehicle.y.kinematics.e_nb
    mode_prev = sys.y.mode


    if mode === lon_direct #actuation commands set from pilot setpoints

        throttle_cmd = throttle_sp
        elevator_cmd = elevator_sp

    elseif te2te_enabled(mode) #actuation commands computed by te2te SAS

        u_lon_sat = ULon(te2te_lqr.y.out_sat)

        if v2t_enabled(mode) #throttle_sp overridden by v2t

            Control.Discrete.assign!(v2t_pid, v2t_lookup(EAS, h_e))

            if mode != mode_prev
                Systems.reset!(v2t_pid)
                k_i = v2t_pid.u.k_i
                (k_i != 0) && (v2t_pid.s.x_i0 = Float64(sys.y.throttle_cmd))
            end

            v2t_pid.u.input = EAS_sp - EAS
            v2t_pid.u.sat_ext = u_lon_sat.throttle_cmd
            f_disc!(v2t_pid)
            throttle_sp = v2t_pid.y.output

        end

        if q2e_enabled(mode) #elevator_sp overridden by q2e

            Control.Discrete.assign!(q2e_pid, q2e_lookup(EAS, h_e))

            if mode != mode_prev
                Systems.reset!(q2e_int)
                Systems.reset!(q2e_pid)
                k_i = q2e_pid.u.k_i
                (k_i != 0) && (q2e_pid.s.x_i0 = ZLonThrEle(te2te_lqr.u.z_sp).elevator_cmd)
            end

            if θ2q_enabled(mode) #q_sp overridden by θ2q

                if v2θ_enabled(mode) #θ_sp overridden by v2θ

                    Control.Discrete.assign!(v2θ_pid, v2θ_lookup(EAS, h_e))

                    if mode != mode_prev
                        Systems.reset!(v2θ_pid)
                        k_i = v2θ_pid.u.k_i
                        (k_i != 0) && (v2θ_pid.s.x_i0 = -θ) #sign inversion!
                    end

                    v2θ_pid.u.input = EAS_sp - EAS
                    v2θ_pid.u.sat_ext = -u_lon_sat.elevator_cmd #sign inversion!
                    f_disc!(v2θ_pid)
                    θ_sp = -v2θ_pid.y.output #sign inversion!

                end

                k_p_θ = 1.0
                θ_dot_sp = k_p_θ * (θ_sp - θ)
                φ_bnd = clamp(φ, -π/3, π/3)
                q_sp = 1/cos(φ_bnd) * θ_dot_sp + r * tan(φ_bnd)

            end

            q2e_int.u.input = q_sp - q
            q2e_int.u.sat_ext = u_lon_sat.elevator_cmd
            f_disc!(q2e_int)

            q2e_pid.u.input = q2e_int.y.output
            q2e_pid.u.sat_ext = u_lon_sat.elevator_cmd
            f_disc!(q2e_pid)
            elevator_sp = q2e_pid.y.output

        end

        Control.Discrete.assign!(te2te_lqr, te2te_lookup(EAS, h_e))

        #te2te is purely proportional, so it doesn't need resetting

        te2te_lqr.u.x .= XLon(vehicle) #state feedback
        te2te_lqr.u.z .= ZLonThrEle(vehicle) #command variable feedback
        te2te_lqr.u.z_sp .= ZLonThrEle(; throttle_cmd = throttle_sp, elevator_cmd = elevator_sp) #command variable setpoint
        f_disc!(te2te_lqr)
        @unpack throttle_cmd, elevator_cmd = ULon(te2te_lqr.y.output)

    else #mode === lon_EAS_clm #actuation commands computed by vc2te

        Control.Discrete.assign!(vc2te_lqr, vc2te_lookup(EAS, h_e))

        #we don't need anything else for a smooth mode transition; as long as we
        #are in steady state and EAS_sp and clm_sp are set to the current EAS
        #and clm values, the throttle and elevator commands from the LQR should
        #match its current values
        (mode != mode_prev) && Systems.reset!(vc2te_lqr)

        vc2te_lqr.u.x .= XLon(vehicle) #state feedback
        vc2te_lqr.u.z .= ZLonEASClm(vehicle) #command variable feedback
        vc2te_lqr.u.z_sp .= ZLonEASClm(; EAS = EAS_sp, climb_rate = clm_sp) #command variable setpoint
        f_disc!(vc2te_lqr)
        @unpack throttle_cmd, elevator_cmd = ULon(vc2te_lqr.y.output)

    end

    sys.y = LonControlY(; mode, throttle_sp, elevator_sp, q_sp, θ_sp, EAS_sp, clm_sp,
        throttle_cmd, elevator_cmd, te2te_lqr = te2te_lqr.y, vc2te_lqr = vc2te_lqr.y,
        q2e_int = q2e_int.y, q2e_pid = q2e_pid.y, v2θ_pid = v2θ_pid.y, v2t_pid = v2t_pid.y)

end


################################################################################
################################# LatControl ###################################

@enum LatControlMode begin
    lat_direct = 0 #direct aileron_cmd + rudder_cmd
    lat_p_β = 1 #roll rate + sideslip
    lat_φ_β = 2 #bank angle + sideslip
    lat_χ_β = 4 #couse angle + sideslip
end

################################# FieldVectors #################################

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

#command vector for φ + β SAS
@kwdef struct ZLatPhiBeta <: FieldVector{2, Float64}
    φ::Float64 = 0.0; β::Float64 = 0.0
end


#assemble state vector from vehicle
function XLat(vehicle::System{<:C172FBW.Vehicle})

    @unpack components, air, kinematics = vehicle.y
    @unpack aero, act = components
    @unpack e_nb, ω_eb_b = kinematics

    p, _, r = ω_eb_b
    φ = e_nb.φ
    v_x, v_y, _ = air.v_wOb_b
    β_filt = aero.β_filt
    ail_v = act.aileron.vel
    ail_p = act.aileron.pos
    rud_v = act.rudder.vel
    rud_p = act.rudder.pos

    XLat(; p, r, φ, v_x, v_y, β_filt, ail_v, ail_p, rud_v, rud_p)

end

function ZLatPhiBeta(vehicle::System{<:C172FBW.Vehicle})
    φ = vehicle.y.kinematics.data.e_nb.φ
    β = vehicle.y.air.β_b
    ZLatPhiBeta(; φ, β)
end

################################## System ######################################

@kwdef struct LatControl{LQ <: LQRTrackerLookup, LP <: PIDLookup} <: AbstractControlChannel
    φβ2ar_lookup::LQ = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "φβ2ar_lookup.h5"))
    p2φ_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "p2φ_lookup.h5"))
    χ2φ_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "χ2φ_lookup.h5"))
    φβ2ar_lqr::LQRTracker{10, 2, 2, 20, 4} = LQRTracker{10, 2, 2}()
    p2φ_int::Integrator = Integrator()
    p2φ_pid::PID = PID()
    χ2φ_pid::PID = PID()
end

@kwdef mutable struct LatControlU
    mode::LatControlMode = lat_direct #lateral control mode
    aileron_sp::Float64 = 0.0 #aileron command setpoint
    rudder_sp::Float64 = 0.0 #rudder command setpoint
    p_sp::Float64 = 0.0 #roll rate setpoint
    β_sp::Float64 = 0.0 #sideslip angle setpoint
    φ_sp::Float64 = 0.0 #bank angle setpoint
    χ_sp::Float64 = 0.0 #course angle setpoint
end

@kwdef struct LatControlY
    mode::LatControlMode = lat_direct
    aileron_sp::Float64 = 0.0 #aileron command setpoint
    rudder_sp::Float64 = 0.0 #rudder command setpoint
    p_sp::Float64 = 0.0 #roll rate setpoint
    β_sp::Float64 = 0.0 #sideslip angle setpoint
    φ_sp::Float64 = 0.0 #bank angle setpoint
    χ_sp::Float64 = 0.0 #course angle setpoint
    aileron_cmd::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd::Ranged{Float64, -1., 1.} = 0.0
    φβ2ar_lqr::LQRTrackerOutput{10, 2, 2, 20, 4} = LQRTrackerOutput{10, 2, 2}()
    p2φ_int::IntegratorOutput = IntegratorOutput()
    p2φ_pid::PIDOutput = PIDOutput()
    χ2φ_pid::PIDOutput = PIDOutput()
end

Systems.U(::LatControl) = LatControlU()
Systems.Y(::LatControl) = LatControlY()

function Systems.init!(sys::System{<:LatControl})

    foreach((sys.φβ2ar_lqr, )) do lqr
        lqr.u.bound_lo .= ULat(; aileron_cmd = -1, rudder_cmd = -1)
        lqr.u.bound_hi .= ULat(; aileron_cmd = 1, rudder_cmd = 1)
    end

    #set φ setpoint limits for the course angle compensator output
    sys.χ2φ_pid.u.bound_lo = -π/4
    sys.χ2φ_pid.u.bound_hi = π/4

end

function Systems.f_disc!(::NoScheduling, sys::System{<:LatControl},
                        vehicle::System{<:C172FBW.Vehicle})

    @unpack mode, aileron_sp, rudder_sp, p_sp, β_sp, φ_sp, χ_sp = sys.u
    @unpack φβ2ar_lqr, p2φ_int, p2φ_pid, χ2φ_pid = sys.subsystems
    @unpack φβ2ar_lookup, p2φ_lookup, χ2φ_lookup = sys.constants
    @unpack air, kinematics = vehicle.y

    EAS = air.EAS
    h_e = Float64(kinematics.h_e)
    φ = kinematics.e_nb.φ
    mode_prev = sys.y.mode

    if mode === lat_direct
        aileron_cmd = aileron_sp
        rudder_cmd = rudder_sp

    #compute φ_sp depending on active roll path
    else #lat_p_β || lat_φ_β || lat_χ_β

        u_lat_sat = ULat(φβ2ar_lqr.y.out_sat)

        #φ_sp computed by roll rate tracker
        if mode === lat_p_β

            Control.Discrete.assign!(p2φ_pid, p2φ_lookup(EAS, Float64(h_e)))

            if mode != mode_prev
                Systems.reset!(p2φ_int)
                Systems.reset!(p2φ_pid)
                k_i = p2φ_pid.u.k_i
                (k_i != 0) && (p2φ_pid.s.x_i0 = ZLatPhiBeta(φβ2ar_lqr.u.z_sp).φ)
            end

            p = kinematics.ω_lb_b[1]
            p2φ_int.u.input = p_sp - p
            p2φ_int.u.sat_ext = u_lat_sat.aileron_cmd
            f_disc!(p2φ_int)

            p2φ_pid.u.input = p2φ_int.y.output
            p2φ_pid.u.sat_ext = u_lat_sat.aileron_cmd
            f_disc!(p2φ_pid)
            φ_sp = p2φ_pid.y.output

        elseif mode === lat_χ_β

            Control.Discrete.assign!(χ2φ_pid, χ2φ_lookup(EAS, Float64(h_e)))

            if mode != mode_prev
                Systems.reset!(χ2φ_pid)
                k_i = χ2φ_pid.u.k_i
                (k_i != 0) && (χ2φ_pid.s.x_i0 = ZLatPhiBeta(φβ2ar_lqr.u.z_sp).φ)
            end

            χ = kinematics.χ_gnd
            χ2φ_pid.u.input = wrap_to_π(χ_sp - χ)
            χ2φ_pid.u.sat_ext = u_lat_sat.aileron_cmd
            f_disc!(χ2φ_pid)
            φ_sp = χ2φ_pid.y.output

        else #mode === lat_φ_β

            #φ_sp and β_sp directly set by input values, nothing to do here

        end

        Control.Discrete.assign!(φβ2ar_lqr, φβ2ar_lookup(EAS, Float64(h_e)))

        (mode != mode_prev) && Systems.reset!(φβ2ar_lqr)

        φβ2ar_lqr.u.x .= XLat(vehicle)
        φβ2ar_lqr.u.z .= ZLatPhiBeta(vehicle)
        φβ2ar_lqr.u.z_sp .= ZLatPhiBeta(; φ = φ_sp, β = β_sp)
        f_disc!(φβ2ar_lqr)
        @unpack aileron_cmd, rudder_cmd = ULat(φβ2ar_lqr.y.output)

    end

    sys.y = LatControlY(; mode, aileron_sp, rudder_sp, p_sp, β_sp, φ_sp, χ_sp,
        aileron_cmd, rudder_cmd, φβ2ar_lqr = φβ2ar_lqr.y, p2φ_int = p2φ_int.y,
        p2φ_pid = p2φ_pid.y, χ2φ_pid = χ2φ_pid.y)

end



################################################################################
############################### Guidance Modes #################################

@enum AltGuidanceState begin
    alt_acquire = 0
    alt_hold = 1
end

@kwdef struct AltitudeGuidance <: AbstractControlChannel
    k_h2c::Float64 = 0.2 #good margins for the whole envelope
end

@kwdef mutable struct AltitudeGuidanceU
    h_sp::Union{HEllip, HOrth} = HEllip(0.0) #altitude setpoint
end

@kwdef mutable struct AltitudeGuidanceS
    state::AltGuidanceState = alt_hold
    h_thr::Float64 = 10.0 #current a
end

@kwdef struct AltitudeGuidanceY
    state::AltGuidanceState = alt_hold
    h_thr::Float64 = 0.0 #current altitude switching threshold
    lon_ctl_mode::LonControlMode = lon_EAS_clm
    throttle_sp::Float64 = 0.0
    clm_sp::Float64 = 0.0
end

Systems.U(::AltitudeGuidance) = AltitudeGuidanceU()
Systems.S(::AltitudeGuidance) = AltitudeGuidanceS()
Systems.Y(::AltitudeGuidance) = AltitudeGuidanceY()

get_Δh(h_sp::HEllip, vehicle::System{<:C172FBW.Vehicle}) = h_sp - vehicle.y.kinematics.h_e
get_Δh(h_sp::HOrth, vehicle::System{<:C172FBW.Vehicle}) = h_sp - vehicle.y.kinematics.h_o

function Systems.f_disc!(::NoScheduling, sys::System{<:AltitudeGuidance},
                        vehicle::System{<:C172FBW.Vehicle})

    Δh = get_Δh(sys.u.h_sp, vehicle)
    clm_sp = sys.constants.k_h2c * Δh
    @unpack state, h_thr = sys.s
    # @show Δh
    # @show state

    if state === alt_acquire

        lon_ctl_mode = lon_thr_EAS
        throttle_sp = Δh > 0 ? 1.0 : 0.0 #full throttle to climb, idle to descend
        # (abs(Δh) + 1 < h_thr) && (sys.s.state = alt_hold)
        sys.s.h_thr = abs(Δh)
        clm = -vehicle.y.kinematics.v_eOb_n[3]
        (abs(clm_sp) < abs(clm)) && (sys.s.state = alt_hold)

    else #alt_hold

        lon_ctl_mode = lon_EAS_clm
        throttle_sp = 0.0 #no effect, controlled by EAS_clm
        (abs(Δh) - 1 > h_thr) && (sys.s.state = alt_acquire)

    end

    sys.y = AltitudeGuidanceY(; state, h_thr, lon_ctl_mode, throttle_sp, clm_sp)

end

@kwdef struct SegmentGuidance <: SystemDefinition end
@kwdef struct SegmentGuidanceY end


################################################################################
#################################### Controller ##################################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@enum VerticalGuidanceMode begin
    vrt_gdc_off = 0
    vrt_gdc_alt = 1
end

@enum HorizontalGuidanceMode begin
    hor_gdc_off = 0
    hor_gdc_line = 1
end

################################################################################

@kwdef struct Controller{C1 <: LonControl, C2 <: LatControl} <: SystemDefinition
    lon_ctl::C1 = LonControl()
    lat_ctl::C2 = LatControl()
    alt_gdc::AltitudeGuidance = AltitudeGuidance()
    seg_gdc::SegmentGuidance = SegmentGuidance()
end

#CockpitInputs
@kwdef mutable struct ControllerU
    eng_start::Bool = false #passthrough
    eng_stop::Bool = false #passthrough
    mixture::Ranged{Float64, 0., 1.} = 0.5 #passthrough
    flaps::Ranged{Float64, 0., 1.} = 0.0 #passthrough
    steering::Ranged{Float64, -1., 1.} = 0.0 #passthrough
    brake_left::Ranged{Float64, 0., 1.} = 0.0 #passthrough
    brake_right::Ranged{Float64, 0., 1.} = 0.0 #passthrough
    throttle_sp_input::Ranged{Float64, 0., 1.} = 0.0 #sets throttle_sp
    aileron_sp_input::Ranged{Float64, -1., 1.} = 0.0 #sets aileron_sp or p_sp
    elevator_sp_input::Ranged{Float64, -1., 1.} = 0.0 #sets elevator_sp or q_sp
    rudder_sp_input::Ranged{Float64, -1., 1.} = 0.0 #sets rudder_sp or β_sp
    throttle_sp_offset::Ranged{Float64, 0., 1.} = 0.0 #for direct throttle control only
    aileron_sp_offset::Ranged{Float64, -1., 1.} = 0.0 #for direct aileron control only
    elevator_sp_offset::Ranged{Float64, -1., 1.} = 0.0 #for direct elevator control only
    rudder_sp_offset::Ranged{Float64, -1., 1.} = 0.0 #for direct rudder control only
    vrt_gdc_mode_req::VerticalGuidanceMode = vrt_gdc_off #requested vertical guidance mode
    hor_gdc_mode_req::HorizontalGuidanceMode = hor_gdc_off #requested horizontal guidance mode
    lon_ctl_mode_req::LonControlMode = lon_direct #requested longitudinal control mode
    lat_ctl_mode_req::LatControlMode = lat_direct #requested lateral control mode
    EAS_sp::Float64 = C172.TrimParameters().EAS #equivalent airspeed setpoint
    q_sp::Float64 = 0.0 #pitch rate setpoint
    θ_sp::Float64 = 0.0 #pitch angle setpoint
    clm_sp::Float64 = 0.0 #climb rate setpoint
    p_sp::Float64 = 0.0 #roll rate setpoint
    φ_sp::Float64 = 0.0 #bank angle setpoint
    χ_sp::Float64 = 0.0 #course angle setpoint
    β_sp::Float64 = 0.0 #sideslip angle setpoint
    h_sp::Union{HEllip, HOrth} = HEllip(0.0) #altitude setpoint
end

@kwdef struct ControllerY
    flight_phase::FlightPhase = phase_gnd
    vrt_gdc_mode::VerticalGuidanceMode = vrt_gdc_off #active vertical guidance mode
    hor_gdc_mode::HorizontalGuidanceMode = hor_gdc_off #active horizontal guidance mode
    lon_ctl_mode::LonControlMode = lon_direct #active longitudinal control mode
    lat_ctl_mode::LatControlMode = lat_direct #active lateral control mode
    alt_gdc::AltitudeGuidanceY = AltitudeGuidanceY()
    seg_gdc::SegmentGuidanceY = SegmentGuidanceY()
    lon_ctl::LonControlY = LonControlY()
    lat_ctl::LatControlY = LatControlY()
end

Systems.U(::Controller) = ControllerU()
Systems.Y(::Controller) = ControllerY()


########################### Update Methods #####################################


function Systems.f_disc!(::NoScheduling, sys::System{<:Controller},
                        vehicle::System{<:C172FBW.Vehicle})

    @unpack eng_start, eng_stop, mixture, flaps, steering, brake_left, brake_right,
            throttle_sp_input, aileron_sp_input, elevator_sp_input, rudder_sp_input,
            throttle_sp_offset, aileron_sp_offset, elevator_sp_offset, rudder_sp_offset,
            vrt_gdc_mode_req, hor_gdc_mode_req, lon_ctl_mode_req, lat_ctl_mode_req,
            q_sp, EAS_sp, θ_sp, clm_sp, p_sp, φ_sp, χ_sp, β_sp, h_sp = sys.u

    @unpack lon_ctl, lat_ctl, alt_gdc, seg_gdc = sys.subsystems

    throttle_sp = throttle_sp_input + throttle_sp_offset
    elevator_sp = elevator_sp_input + elevator_sp_offset
    aileron_sp = aileron_sp_input + aileron_sp_offset
    rudder_sp = rudder_sp_input + rudder_sp_offset

    any_wow = any(SVector{3}(leg.strut.wow for leg in vehicle.y.components.ldg))
    flight_phase = any_wow ? phase_gnd : phase_air

    if flight_phase === phase_gnd

        vrt_gdc_mode = vrt_gdc_off
        hor_gdc_mode = hor_gdc_off
        lon_ctl_mode = lon_direct
        lat_ctl_mode = lat_direct

    elseif flight_phase === phase_air

        vrt_gdc_mode = vrt_gdc_mode_req
        hor_gdc_mode = hor_gdc_mode_req

        if vrt_gdc_mode === vrt_gdc_off

            lon_ctl_mode = lon_ctl_mode_req

        else #vrt_gdc_mode === vrt_gdc_alt

            alt_gdc.u.h_sp = h_sp
            f_disc!(alt_gdc, vehicle)

            lon_ctl_mode = alt_gdc.y.lon_ctl_mode
            throttle_sp = alt_gdc.y.throttle_sp
            clm_sp = alt_gdc.y.clm_sp

        end

        if hor_gdc_mode === hor_gdc_off

            lat_ctl_mode = lat_ctl_mode_req

            #below a v_gnd threshold, override χ mode and revert to φ
            if (lat_ctl_mode === lat_χ_β) && (vehicle.y.kinematics.v_gnd < 10.0)
                lat_ctl_mode = lat_φ_β
            end

        else #hor_gdc_mode === hor_gdc_line

            # seg_gdc.u.line_sp = line_sp
            # f_disc!(seg_gdc, vehicle)

            # lat_ctl_mode = seg_gdc.y.lat_ctl_mode
            # χ_sp = seg_gdc.y.χ_sp
            # β_sp = 0.0

        end

    end

    lon_ctl.u.mode = lon_ctl_mode
    @pack! lon_ctl.u = throttle_sp, elevator_sp, q_sp, θ_sp, EAS_sp, clm_sp
    f_disc!(lon_ctl, vehicle)

    lat_ctl.u.mode = lat_ctl_mode
    @pack! lat_ctl.u = aileron_sp, rudder_sp, p_sp, φ_sp, β_sp, χ_sp
    f_disc!(lat_ctl, vehicle)

    sys.y = ControllerY(; flight_phase,
        vrt_gdc_mode, hor_gdc_mode, lon_ctl_mode, lat_ctl_mode,
        lon_ctl = lon_ctl.y, lat_ctl = lat_ctl.y,
        alt_gdc = alt_gdc.y, seg_gdc = seg_gdc.y)

end

function AircraftBase.assign!(components::System{<:C172FBW.Components},
                          sys::System{<:Controller})

    @unpack act, pwp, ldg = components.subsystems
    @unpack eng_start, eng_stop, mixture, flaps, steering, brake_left, brake_right = sys.u
    @unpack throttle_cmd, elevator_cmd = sys.lon_ctl.y
    @unpack aileron_cmd, rudder_cmd = sys.lat_ctl.y

    act.throttle.u[] = throttle_cmd
    act.aileron.u[] = aileron_cmd
    act.elevator.u[] = elevator_cmd
    act.rudder.u[] = rudder_cmd
    act.flaps.u[] = flaps
    act.steering.u[] = steering
    pwp.engine.u.start = eng_start
    pwp.engine.u.stop = eng_stop
    pwp.engine.u.mixture = mixture
    ldg.left.braking.u[] = brake_left
    ldg.right.braking.u[] = brake_right

end





##################################### Tools ####################################

function AircraftBase.trim!(sys::System{<:Controller},
                            vehicle::System{<:C172FBW.Vehicle})

    #here we assume that the vehicle's y has already been updated to its trim
    #value by trim!(vehicle, params)
    y_act = vehicle.y.components.act
    @unpack ω_lb_b, v_eOb_n, e_nb, χ_gnd, h_e = vehicle.y.kinematics
    @unpack EAS, β_b = vehicle.y.air

    #makes Controller inputs consistent with the trim solution obtained for the
    #vehicle, so the trim condition is preserved upon simulation start
    #when the corresponding control modes are selected
    Systems.reset!(sys)

    #in a fly-by-wire implementation, it makes more sense to assign the trim
    #values to the inputs rather to the offsets
    u = sys.u
    u.throttle_sp_input = y_act.throttle.pos
    u.aileron_sp_input = y_act.aileron.pos
    u.elevator_sp_input = y_act.elevator.pos
    u.rudder_sp_input = y_act.rudder.pos
    u.throttle_sp_offset = 0
    u.aileron_sp_offset = 0
    u.elevator_sp_offset = 0
    u.rudder_sp_offset = 0
    u.flaps = y_act.flaps.pos
    u.mixture = vehicle.y.components.pwp.engine.mixture

    u.q_sp = ω_lb_b[2]
    u.θ_sp = e_nb.θ
    u.EAS_sp = EAS
    u.clm_sp = -v_eOb_n[3]
    u.p_sp = ω_lb_b[1]
    u.φ_sp = e_nb.φ
    u.β_sp = β_b
    u.χ_sp = χ_gnd
    u.h_sp = h_e

    u.vrt_gdc_mode_req = vrt_gdc_off
    u.hor_gdc_mode_req = hor_gdc_off

    #do an update with the inner SAS loops enabled so that their internal
    #setpoints are made consistent with the trim values. this should propagate
    #to the actuator commands produced by the control laws, making them
    #consistent with the trim state ones
    u.lon_ctl_mode_req = lon_thr_ele
    u.lat_ctl_mode_req = lat_φ_β
    f_disc!(sys, vehicle)

    #restore direct modes
    u.lon_ctl_mode_req = lon_direct
    u.lat_ctl_mode_req = lat_direct
    f_disc!(sys, vehicle)

end


################################################################################
############################ Joystick Mappings #################################

pitch_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
roll_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
yaw_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)


function Systems.assign_input!(sys::System{<:Controller},
                           data::XBoxController, ::IOMapping)

    u = sys.u

    p_sf = 0.5 #roll rate sensitivity
    q_sf = 0.5 #pitch rate sensitivity

    roll_input = get_axis_value(data, :right_stick_x) |> roll_curve
    pitch_input = get_axis_value(data, :right_stick_y) |> pitch_curve
    yaw_input = get_axis_value(data, :left_stick_x) |> yaw_curve

    u.throttle_sp_input += 0.1 * was_released(data, :button_Y)
    u.throttle_sp_input -= 0.1 * was_released(data, :button_A)

    u.aileron_sp_input = roll_input
    u.elevator_sp_input = pitch_input
    u.rudder_sp_input = yaw_input

    u.p_sp = p_sf * roll_input
    u.q_sp = q_sf * pitch_input

    u.steering = get_axis_value(data, :left_stick_x) |> yaw_curve
    u.brake_left = get_axis_value(data, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(data, :right_trigger) |> brake_curve

    u.aileron_sp_offset -= 1e-3 * was_released(data, :dpad_left)
    u.aileron_sp_offset += 1e-3 * was_released(data, :dpad_right)
    u.elevator_sp_offset += 1e-3 * was_released(data, :dpad_down)
    u.elevator_sp_offset -= 1e-3 * was_released(data, :dpad_up)

    u.flaps += 0.3333 * was_released(data, :right_bumper)
    u.flaps -= 0.3333 * was_released(data, :left_bumper)

end

function Systems.assign_input!(sys::System{<:Controller},
                                data::T16000MData, ::IOMapping)

    u = sys.u

    p_sf = 0.5 #roll rate sensitivity
    q_sf = 0.5 #pitch rate sensitivity

    throttle_input = get_axis_value(data, :throttle)
    roll_input = get_axis_value(data, :stick_x) |> roll_curve
    pitch_input = get_axis_value(data, :stick_y) |> pitch_curve
    yaw_input = get_axis_value(data, :stick_z) |> yaw_curve

    u.throttle_sp_input = throttle_input
    u.aileron_sp_input = roll_input
    u.elevator_sp_input = pitch_input
    u.rudder_sp_input = yaw_input

    u.p_sp = p_sf * roll_input
    u.q_sp = q_sf * pitch_input

    u.steering = get_axis_value(data, :stick_z) |> yaw_curve
    u.brake_left = is_pressed(data, :button_1)
    u.brake_right = is_pressed(data, :button_1)

    u.aileron_sp_offset -= 1e-3 * is_pressed(data, :hat_left)
    u.aileron_sp_offset += 1e-3 * is_pressed(data, :hat_right)
    u.elevator_sp_offset += 1e-3 * is_pressed(data, :hat_down)
    u.elevator_sp_offset -= 1e-3 * is_pressed(data, :hat_up)

    u.flaps += 0.3333 * was_released(data, :button_3)
    u.flaps -= 0.3333 * was_released(data, :button_2)

end

function Systems.assign_input!(sys::System{<:Controller},
                           data::GladiatorNXTEvoData, ::IOMapping)

    u = sys.u

    p_sf = 0.5 #roll rate sensitivity
    q_sf = 0.5 #pitch rate sensitivity

    throttle_input = get_axis_value(data, :throttle)
    roll_input = get_axis_value(data, :stick_x) |> roll_curve
    pitch_input = get_axis_value(data, :stick_y) |> pitch_curve
    yaw_input = get_axis_value(data, :stick_z) |> yaw_curve

    u.throttle_sp_input = throttle_input
    u.aileron_sp_input = roll_input
    u.elevator_sp_input = pitch_input
    u.rudder_sp_input = yaw_input

    u.p_sp = p_sf * roll_input
    u.q_sp = q_sf * pitch_input

    u.steering = get_axis_value(data, :stick_z) |> yaw_curve
    u.brake_left = is_pressed(data, :red_trigger_half)
    u.brake_right = is_pressed(data, :red_trigger_half)

    u.aileron_sp_offset -= 1e-3 * is_pressed(data, :A3_left)
    u.aileron_sp_offset += 1e-3 * is_pressed(data, :A3_right)
    u.elevator_sp_offset += 1e-3 * is_pressed(data, :A3_down)
    u.elevator_sp_offset -= 1e-3 * is_pressed(data, :A3_up)

    if is_pressed(data, :A3_press)
        u.aileron_sp_offset = 0
        u.elevator_sp_offset = 0
    end

    u.flaps += 0.3333 * was_released(data, :switch_down)
    u.flaps -= 0.3333 * was_released(data, :switch_up)

end

################################################################################
################################## GUI #########################################

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
        Dummy, SameLine, NewLine, IsItemActive, Separator, Text, Checkbox, RadioButton


function GUI.draw(sys::System{<:LonControl}, p_open::Ref{Bool} = Ref(true))
    Begin("Longitudinal Control", p_open)
        Text("Mode: $(sys.y.mode)")
        foreach(keys(sys.subsystems), values(sys.subsystems)) do label, ss
            if CImGui.TreeNode(string(label))
                GUI.draw(ss)
                CImGui.TreePop()
            end
        end
    End()
end

function GUI.draw(sys::System{<:LatControl}, p_open::Ref{Bool} = Ref(true))
    Begin("Lateral Control", p_open)
        Text("Mode: $(sys.y.mode)")
        foreach(keys(sys.subsystems), values(sys.subsystems)) do label, ss
            if CImGui.TreeNode(string(label))
                GUI.draw(ss)
                CImGui.TreePop()
            end
        end
    End()
end

function GUI.draw(sys::System{<:AltitudeGuidance}, p_open::Ref{Bool} = Ref(true))
    Begin("Altitude Guidance", p_open)
        Text("State: $(sys.y.state)")
        Text("Control Mode: $(sys.y.lon_ctl_mode)")
        Text("Altitude Threshold: $(sys.y.h_thr)")
        Text("Throttle Setpoint: $(sys.y.throttle_sp)")
        Text("Climb Rate Setpoint: $(sys.y.clm_sp)")
    End()
end


function GUI.draw!(sys::System{<:Controller},
                    vehicle::System{<:C172FBW.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172 FBW Controller")

    u = sys.u
    y = sys.y

    @unpack components, kinematics, dynamics, air = vehicle.y
    @unpack act, pwp, fuel, ldg = components

    @unpack e_nb, ω_lb_b, n_e, ϕ_λ, h_e, h_o, v_gnd, χ_gnd, γ_gnd, v_eOb_n = kinematics
    @unpack CAS, EAS, TAS, α_b, β_b, T, p, pt = air
    @unpack ψ, θ, φ = e_nb
    @unpack ϕ, λ = ϕ_λ

    p, q, r = ω_lb_b
    clm = -v_eOb_n[3]
    Δh = h_o - TerrainData(vehicle.constants.terrain, n_e).altitude

    Begin(label, p_open)

    ############################# Engine Control ###############################

    if CImGui.CollapsingHeader("Engine")
        PushItemWidth(-100)
        mode_button("Engine Start", true, pwp.engine.state === Piston.eng_starting, pwp.engine.state === Piston.eng_running)
        u.eng_start = IsItemActive()
        SameLine()
        mode_button("Engine Stop", true, u.eng_stop, false; HSV_requested = HSV_red)
        u.eng_stop = IsItemActive()
        SameLine()
        u.mixture = safe_slider("Mixture", u.mixture, "%.6f", true)
        @unpack state, throttle, ω, MAP, M_shaft, P_shaft, ṁ = pwp.engine

        if CImGui.BeginTable("Engine Data", 4)
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("State: $(state)")
                CImGui.TableNextColumn(); Text(@sprintf("Throttle: %.3f %%", 100*throttle))
                CImGui.TableNextColumn(); Text(@sprintf("Speed: %.3f RPM", Piston.radpersec2RPM(ω)))
                CImGui.TableNextColumn(); Text(@sprintf("Manifold Pressure: %.0f Pa", MAP))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text(@sprintf("Shaft Torque: %.3f Nm", M_shaft))
                CImGui.TableNextColumn(); Text(@sprintf("Shaft Power: %.3f kW", P_shaft/1e3))
                CImGui.TableNextColumn(); Text(@sprintf("Fuel Flow: %.3f g/s", ṁ*1e3)); SameLine(250)
                CImGui.TableNextColumn(); Text(@sprintf("Remaining Fuel: %.3f kg", fuel.m_avail))
            CImGui.EndTable()
        end
        Separator()
        PopItemWidth()
    end

    ############################### Guidance ###################################

    @cstatic h_datum=Int32(0) begin
    if CImGui.CollapsingHeader("Vertical Guidance")

        if CImGui.BeginTable("VerticalGuidance", 3, CImGui.ImGuiTableFlags_SizingStretchProp )#| CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    Text("Mode")
                CImGui.TableNextColumn();
                    mode_button("Off", vrt_gdc_off, u.vrt_gdc_mode_req, y.vrt_gdc_mode)
                    IsItemActive() && (u.vrt_gdc_mode_req = vrt_gdc_off); SameLine()
                    mode_button("Altitude", vrt_gdc_alt, u.vrt_gdc_mode_req, y.vrt_gdc_mode)
                    if IsItemActive()
                        u.vrt_gdc_mode_req = vrt_gdc_alt
                        u.h_sp = h_datum == 0 ? h_e : h_o
                    end
                CImGui.TableNextColumn();
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Altitude (m)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    h_val = safe_input("Altitude Setpoint", Float64(u.h_sp), 1, 1.0, "%.3f")
                    PopItemWidth()
                CImGui.TableNextColumn();
                    isa(u.h_sp, HEllip) && Text(@sprintf("%.3f", Float64(h_e)))
                    isa(u.h_sp, HOrth) && Text(@sprintf("%.3f", Float64(h_o)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                CImGui.TableNextColumn();
                    @c RadioButton("Ellipsoidal", &h_datum, 0); SameLine()
                    @c RadioButton("Orthometric", &h_datum, 1)
                    u.h_sp = h_datum == 0 ? HEllip(h_val) : HOrth(h_val)
            CImGui.EndTable()
        end #table

        Separator()

    end #header
    end #cstatic


    ########################## Longitudinal Control ########################

    if CImGui.CollapsingHeader("Longitudinal Control")

        if CImGui.BeginTable("LonCtlModes", 3, CImGui.ImGuiTableFlags_SizingStretchProp )#| CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    Text("Mode")
                CImGui.TableNextColumn();
                    mode_button("Direct##Lon", lon_direct, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_direct); SameLine()
                    mode_button("Throttle + Pitch SAS", lon_thr_ele, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_thr_ele); SameLine()
                    mode_button("Throttle + Pitch Rate", lon_thr_q, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_thr_q; u.q_sp = 0); SameLine()
                    mode_button("Throttle + Pitch Angle", lon_thr_θ, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_thr_θ; u.θ_sp = θ)
                    mode_button("EAS + Throttle", lon_thr_EAS, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_thr_EAS; u.EAS_sp = EAS); SameLine()
                    mode_button("EAS + Pitch Rate", lon_EAS_q, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_EAS_q; u.q_sp = 0; u.EAS_sp = EAS); SameLine()
                    mode_button("EAS + Pitch Angle", lon_EAS_θ, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_EAS_θ; u.EAS_sp = EAS; u.θ_sp = θ); SameLine()
                    mode_button("EAS + Climb Rate", lon_EAS_clm, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_EAS_clm; u.EAS_sp = EAS; u.clm_sp = clm)
                CImGui.TableNextColumn();
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    AlignTextToFramePadding(); Text("Throttle Input")
                    AlignTextToFramePadding(); Text("Throttle Offset")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.throttle_sp_input = safe_slider("Throttle Input", u.throttle_sp_input, "%.6f")
                    u.throttle_sp_offset = safe_input("Throttle_Offset", u.throttle_sp_offset, 0.1, 0.1, "%.3f")
                    PopItemWidth()
                CImGui.TableNextColumn();
                    Text(@sprintf("%.3f", Float64(y.lon_ctl.throttle_sp)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    AlignTextToFramePadding(); Text("Elevator Input")
                    AlignTextToFramePadding(); Text("Elevator Offset")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.elevator_sp_input = safe_slider("Elevator Input", u.elevator_sp_input, "%.6f")
                    u.elevator_sp_offset = safe_input("Elevator Offset", u.elevator_sp_offset, 0.1, 0.1, "%.3f")
                    PopItemWidth()
                CImGui.TableNextColumn();
                    Text(@sprintf("%.3f", Float64(y.lon_ctl.elevator_sp)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Pitch Rate (deg/s)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.q_sp = safe_input("Pitch Rate", rad2deg(u.q_sp), 0.1, 1.0, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(q)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Pitch Angle (deg)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.θ_sp = safe_input("Pitch Angle", rad2deg(u.θ_sp), 0.1, 1.0, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(θ)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("EAS (m/s)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.EAS_sp = safe_input("EAS", u.EAS_sp, 0.1, 1.0, "%.3f")
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", EAS))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Climb Rate (m/s)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.clm_sp = safe_input("Climb Rate", u.clm_sp, 0.1, 1.0, "%.3f")
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", clm))
            CImGui.EndTable()
        end

        Separator()

    end

    ############################### Lateral Control ############################

    if CImGui.CollapsingHeader("Lateral Control")

        if CImGui.BeginTable("LatCtlModes", 3, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    Text("Mode")
                CImGui.TableNextColumn();
                    mode_button("Direct##Lat", lat_direct, u.lat_ctl_mode_req, y.lat_ctl_mode); SameLine()
                    IsItemActive() && (u.lat_ctl_mode_req = lat_direct)
                    mode_button("Roll Rate + AoS", lat_p_β, u.lat_ctl_mode_req, y.lat_ctl_mode); SameLine()
                    IsItemActive() && (u.lat_ctl_mode_req = lat_p_β; u.p_sp = 0; u.β_sp = 0)
                    mode_button("Bank Angle + AoS", lat_φ_β, u.lat_ctl_mode_req, y.lat_ctl_mode); SameLine()
                    IsItemActive() && (u.lat_ctl_mode_req = lat_φ_β; u.φ_sp = φ; u.β_sp = 0)
                    mode_button("Course Angle + AoS", lat_χ_β, u.lat_ctl_mode_req, y.lat_ctl_mode); SameLine()
                    IsItemActive() && (u.lat_ctl_mode_req = lat_χ_β; u.χ_sp = χ_gnd; u.β_sp = 0)
                CImGui.TableNextColumn();
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    AlignTextToFramePadding(); Text("Aileron Input")
                    AlignTextToFramePadding(); Text("Aileron Offset")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.aileron_sp_input = safe_slider("Aileron Input", u.aileron_sp_input, "%.6f")
                    u.aileron_sp_offset = safe_input("Aileron Offset", u.aileron_sp_offset, 0.001, 0.1, "%.3f")
                    PopItemWidth()
                CImGui.TableNextColumn();
                    Text(@sprintf("%.3f", Float64(y.lat_ctl.aileron_sp)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    AlignTextToFramePadding(); Text("Rudder Input")
                    AlignTextToFramePadding(); Text("Rudder Offset")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.rudder_sp_input = safe_slider("Rudder Input", u.rudder_sp_input, "%.6f")
                    u.rudder_sp_offset = safe_input("Rudder Offset", u.rudder_sp_offset, 0.001, 0.1, "%.3f")
                    PopItemWidth()
                CImGui.TableNextColumn();
                    Text(@sprintf("%.3f", Float64(y.lat_ctl.rudder_sp)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Roll Rate (deg/s)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.p_sp = safe_input("Roll Rate", rad2deg(u.p_sp), 0.1, 1.0, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(p)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Bank Angle (deg)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.φ_sp = safe_input("Bank Angle", rad2deg(u.φ_sp), 0.1, 1.0, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(φ)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Course Angle (deg)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.χ_sp = safe_input("Course Angle", rad2deg(u.χ_sp), 0.1, 1.0, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(χ_gnd)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Sideslip Angle (deg)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.β_sp = safe_input("Sideslip Angle", rad2deg(u.β_sp), 0.1, 1.0, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(β_b)))
            CImGui.EndTable()
        end

        Separator()

    end


    ############################################################################

    if CImGui.CollapsingHeader("Secondary Actuation")

        if CImGui.BeginTable("SecondaryActuation", 2, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Flaps")
                CImGui.TableNextColumn();
                PushItemWidth(-10)
                u.flaps = safe_slider("Flaps Input", u.flaps, "%.6f")
                PopItemWidth()
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Steering")
                CImGui.TableNextColumn();
                PushItemWidth(-10)
                u.steering = safe_slider("Steering", u.steering, "%.6f")
                PopItemWidth()
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Left Brake")
                CImGui.TableNextColumn();
                PushItemWidth(-10)
                u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
                PopItemWidth()
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Right Brake")
                CImGui.TableNextColumn();
                PushItemWidth(-10)
                u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")
                PopItemWidth()
            CImGui.EndTable()
        end

        Separator()

    end


    ############################################################################

    if CImGui.CollapsingHeader("Primary Actuator Data")
        if CImGui.BeginTable("Primary Actuator Data", 5, CImGui.ImGuiTableFlags_SizingStretchSame)# | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                CImGui.TableNextColumn(); Text("Throttle")
                CImGui.TableNextColumn(); Text("Elevator")
                CImGui.TableNextColumn(); Text("Aileron")
                CImGui.TableNextColumn(); Text("Rudder")
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Setpoint")
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(y.lon_ctl.throttle_sp)))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(y.lon_ctl.elevator_sp)))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(y.lat_ctl.aileron_sp)))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(y.lat_ctl.rudder_sp)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Command")
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(y.lon_ctl.throttle_cmd)))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(y.lon_ctl.elevator_cmd)))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(y.lat_ctl.aileron_cmd)))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(y.lat_ctl.rudder_cmd)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Position")
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(act.throttle.pos)))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(act.elevator.pos)))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(act.aileron.pos)))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", Float64(act.rudder.pos)))
            CImGui.EndTable()
        end
    end


    if CImGui.CollapsingHeader("Flight Data")

        if CImGui.BeginTable("Flight Data", 2, CImGui.ImGuiTableFlags_SizingStretchSame | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();

                # Text("Flight Phase"); SameLine(240)
                # Text("$(y.flight_phase)")
                Text("Airspeed (Equivalent)"); SameLine(240)
                Text(@sprintf("%.3f m/s | %.3f kts", EAS, Atmosphere.SI2kts(EAS)))
                Text("Airspeed (True)"); SameLine(240)
                Text(@sprintf("%.3f m/s | %.3f kts", TAS, Atmosphere.SI2kts(TAS)))
                Text("Angle of Attack"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(α_b)))
                Text("Angle of Sideslip"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(β_b)))
                Separator()

                Text("Ground Speed"); SameLine(240)
                Text(@sprintf("%.3f m/s | %.3f kts", v_gnd, Atmosphere.SI2kts(v_gnd)))
                Text("Course Angle"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(χ_gnd)))
                Text("Flight Path Angle"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(γ_gnd)))
                Text("Climb Rate"); SameLine(240)
                Text(@sprintf("%.3f m/s", clm))

            CImGui.TableNextColumn();

                Text("Latitude"); SameLine(240)
                Text(@sprintf("%.6f deg", rad2deg(ϕ)))
                Text("Longitude"); SameLine(240)
                Text(@sprintf("%.6f deg", rad2deg(λ)))
                Text("Altitude (Ellipsoidal)"); SameLine(240)
                Text(@sprintf("%.3f m | %.3f ft", Float64(h_e), Float64(h_e)/0.3048))
                Text("Altitude (Orthometric)"); SameLine(240)
                Text(@sprintf("%.3f m | %.3f ft", Float64(h_o), Float64(h_o)/0.3048))
                Text("Height Over Ground"); SameLine(240)
                Text(@sprintf("%.3f m | %.3f ft", Δh, Δh/0.3048))
                Separator()

                Text("Heading"); SameLine(240);
                Text(@sprintf("%.3f deg", rad2deg(ψ)))
                Text("Inclination"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(θ)))
                Text("Bank"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(φ)))

        CImGui.TableNextRow()
            CImGui.TableNextColumn();

                Text("Roll Rate"); SameLine(240)
                Text(@sprintf("%.3f deg/s", rad2deg(p)))
                Text("Pitch Rate"); SameLine(240)
                Text(@sprintf("%.3f deg/s", rad2deg(q)))
                Text("Yaw Rate"); SameLine(240)
                Text(@sprintf("%.3f deg/s", rad2deg(r)))

            CImGui.TableNextColumn();
                Text("Specific Force (x)"); SameLine(240)
                Text(@sprintf("%.3f g", dynamics.f_Gb_b[1]/Dynamics.g₀))
                Text("Specific Force (y)"); SameLine(240)
                Text(@sprintf("%.3f g", dynamics.f_Gb_b[2]/Dynamics.g₀))
                Text("Specific Force (z)"); SameLine(240)
                Text(@sprintf("%.3f g", dynamics.f_Gb_b[3]/Dynamics.g₀))

            CImGui.EndTable()
        end

        Separator()
    end

    if CImGui.CollapsingHeader("Internals")
        @cstatic c_lon=false c_lat=false c_alt=false begin
            @c Checkbox("Longitudinal Control##Internals", &c_lon)
            SameLine()
            @c Checkbox("Lateral Control##Internals", &c_lat)
            SameLine()
            @c Checkbox("Altitude Guidance##Internals", &c_alt)
            c_lon && @c GUI.draw(sys.lon_ctl, &c_lon)
            c_lat && @c GUI.draw(sys.lat_ctl, &c_lat)
            c_alt && @c GUI.draw(sys.alt_gdc, &c_alt)
        end
    end

end

end #module