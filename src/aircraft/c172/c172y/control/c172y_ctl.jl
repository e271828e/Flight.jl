module C172YControl

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes
using EnumX

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightLib.Control.Discrete: Integrator, IntegratorOutput,
    PID, PIDOutput, PIDParams, LQRTracker, LQRTrackerOutput, LQRTrackerParams,
    PIDLookup, LQRTrackerLookup, load_pid_lookup, load_lqr_tracker_lookup

using ...AircraftBase
using ...C172
using ..C172Y

wrap_to_π(x) = x + 2π*floor((π-x)/(2π))

function flaps_schedule(EAS::Real)
    EAS_lo = 30; EAS_hi = 35
    EAS < EAS_lo && return 1.0
    EAS > EAS_hi && return 0.0
    return 1 - (EAS - EAS_lo)/(EAS_hi - EAS_lo)
end


################################################################################
############################### Control ########################################

#a discrete model implementing a specific longitudinal or lateral control mode
abstract type AbstractControlChannel <: ModelDefinition end


################################################################################
################################## ControllerLon ##################################

@enumx T=ModeControlLonEnum ModeControlLon begin
    direct = 0 #direct throttle & elevator
    sas = 1 #direct throttle + pitch SAS
    thr_q = 2 #direct throttle + pitch rate
    thr_θ = 3 #direct throttle + pitch angle
    thr_EAS = 4 #direct throttle + EAS
    EAS_q = 5 #EAS + pitch rate
    EAS_θ = 6 #EAS + pitch angle
    EAS_clm = 7 #EAS + climb rate
end

using .ModeControlLon: ModeControlLonEnum

#elevator pitch SAS always enabled except in direct mode
e2e_enabled(mode::ModeControlLonEnum) = (mode != ModeControlLon.direct)

function q2e_enabled(mode::ModeControlLonEnum)
    mode === ModeControlLon.thr_q || mode === ModeControlLon.thr_θ ||
    mode === ModeControlLon.thr_EAS || mode === ModeControlLon.EAS_q ||
    mode === ModeControlLon.EAS_θ || mode === ModeControlLon.EAS_clm
end

function θ2q_enabled(mode::ModeControlLonEnum)
    mode === ModeControlLon.thr_θ || mode === ModeControlLon.thr_EAS ||
    mode === ModeControlLon.EAS_θ || mode === ModeControlLon.EAS_clm
end

function v2t_enabled(mode::ModeControlLonEnum)
    mode === ModeControlLon.EAS_q || mode === ModeControlLon.EAS_θ ||
    mode === ModeControlLon.EAS_clm
end

c2θ_enabled(mode::ModeControlLonEnum) = (mode === ModeControlLon.EAS_clm)
v2θ_enabled(mode::ModeControlLonEnum) = (mode === ModeControlLon.thr_EAS)

############################## FieldVectors ####################################

#state vector for pitch LQR SAS
@kwdef struct XPitch <: FieldVector{6, Float64}
    q::Float64 = 0.0 #pitch rate
    θ::Float64 = 0.0 #pitch angle
    EAS::Float64 = 0.0 #equivalent airspeed
    α::Float64 = 0.0 #AoA
    α_filt::Float64 = 0.0 #filtered AoA (from aerodynamics model)
    ele_p::Float64 = 0.0 #elevator actuator state
end

#assemble state vector from vehicle
function XPitch(vehicle::Model{<:C172Y.Vehicle})

    @unpack systems, airflow, kinematics = vehicle.y
    @unpack pwp, aero, act = systems
    @unpack e_nb, ω_eb_b = kinematics

    q = ω_eb_b[2]
    θ = e_nb.θ
    EAS = airflow.EAS
    α = aero.α
    α_filt = aero.α_filt
    ele_p = act.elevator.pos

    XPitch(; q, θ, EAS, α, α_filt, ele_p)

end


################################## Model ######################################

@kwdef struct ControllerLon{LQ <: LQRTrackerLookup, LP <: PIDLookup} <: AbstractControlChannel
    e2e_lookup::LQ = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "e2e_lookup.h5"))
    q2e_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "q2e_lookup.h5"))
    v2θ_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "v2θ_lookup.h5"))
    c2θ_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "c2θ_lookup.h5"))
    v2t_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "v2t_lookup.h5"))
    e2e_lqr::LQRTracker{6, 1, 1, 6, 1} = LQRTracker{6, 1, 1}()
    q2e_int::Integrator = Integrator()
    q2e_pid::PID = PID()
    v2θ_pid::PID = PID()
    c2θ_pid::PID = PID()
    v2t_pid::PID = PID()
end

@kwdef mutable struct ControllerLonU
    mode::ModeControlLonEnum = ModeControlLon.direct
    throttle_ref::Float64 = 0.0
    elevator_ref::Float64 = 0.0
    q_ref::Float64 = 0.0
    θ_ref::Float64 = 0.0
    EAS_ref::Float64 = C172.TrimParameters().EAS
    clm_ref::Float64 = 0.0 #climb rate reference
end

@kwdef struct ControllerLonY
    mode::ModeControlLonEnum = ModeControlLon.direct
    throttle_ref::Float64 = 0.0
    elevator_ref::Float64 = 0.0
    q_ref::Float64 = 0.0
    θ_ref::Float64 = 0.0
    EAS_ref::Float64 = C172.TrimParameters().EAS
    clm_ref::Float64 = 0.0 #climb rate reference
    throttle_cmd::Ranged{Float64, 0., 1.} = 0.0
    elevator_cmd::Ranged{Float64, -1., 1.} = 0.0
    e2e_lqr::LQRTrackerOutput{6, 1, 1, 6, 1} = LQRTrackerOutput{6, 1, 1}()
    q2e_int::IntegratorOutput = IntegratorOutput()
    q2e_pid::PIDOutput = PIDOutput()
    v2θ_pid::PIDOutput = PIDOutput()
    c2θ_pid::PIDOutput = PIDOutput()
    v2t_pid::PIDOutput = PIDOutput()
end

Modeling.U(::ControllerLon) = ControllerLonU()
Modeling.Y(::ControllerLon) = ControllerLonY()

function Modeling.init!(mdl::Model{<:ControllerLon})
    #e2e determines elevator saturation for all pitch control loops (q2e, c2θ,
    #v2θ), so we don't need to set bounds for those
    mdl.e2e_lqr.u.bound_lo .= -1
    mdl.e2e_lqr.u.bound_hi .= 1
    #we do need to set bounds for v2t, so it can handle throttle saturation
    mdl.v2t_pid.u.bound_lo = 0
    mdl.v2t_pid.u.bound_hi = 1
end


function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:ControllerLon},
                        vehicle::Model{<:C172Y.Vehicle})

    @unpack mode, throttle_ref, elevator_ref, q_ref, θ_ref, EAS_ref, clm_ref = mdl.u
    @unpack e2e_lqr, q2e_int, q2e_pid, v2θ_pid, c2θ_pid, v2t_pid = mdl.submodels
    @unpack e2e_lookup, q2e_lookup, v2θ_lookup, c2θ_lookup, v2t_lookup = mdl.constants

    EAS = vehicle.y.airflow.EAS
    h_e = Float64(vehicle.y.kinematics.h_e)
    _, q, r = vehicle.y.kinematics.ω_wb_b
    @unpack θ, φ = vehicle.y.kinematics.e_nb
    clm = -vehicle.y.kinematics.v_eb_n[3]
    mode_prev = mdl.y.mode

    #if not overridden by the control modes, actuation commands are simply
    #their respective reference values
    throttle_cmd = throttle_ref
    elevator_cmd = elevator_ref

    if v2t_enabled(mode) #throttle_cmd overridden by v2t

        Control.Discrete.assign!(v2t_pid, v2t_lookup(EAS, h_e))

        if mode != mode_prev
            Control.reset!(v2t_pid)
            k_i = v2t_pid.u.k_i
            (k_i != 0) && (v2t_pid.s.x_i0 = Float64(mdl.y.throttle_cmd))
        end

        v2t_pid.u.input = EAS_ref - EAS
        f_periodic!(v2t_pid)
        throttle_cmd = v2t_pid.y.output

    end

    if e2e_enabled(mode) #elevator_cmd overridden by e2e SAS

        elevator_cmd_sat = e2e_lqr.y.out_sat[1]

        if q2e_enabled(mode) #elevator_ref overridden by q2e

            Control.Discrete.assign!(q2e_pid, q2e_lookup(EAS, h_e))

            if mode != mode_prev

                Control.reset!(q2e_int)
                Control.reset!(q2e_pid)
                k_i = q2e_pid.u.k_i
                (k_i != 0) && (q2e_pid.s.x_i0 = e2e_lqr.u.z_ref[1])

                # (k_i != 0) && (q2e_pid.s.x_i0 = Float64(mdl.y.elevator_cmd))
            end

            if θ2q_enabled(mode) #q_ref overridden by θ2q

                if v2θ_enabled(mode) #θ_ref overridden by v2θ

                    Control.Discrete.assign!(v2θ_pid, v2θ_lookup(EAS, h_e))

                    if mode != mode_prev
                        Control.reset!(v2θ_pid)
                        k_i = v2θ_pid.u.k_i
                        (k_i != 0) && (v2θ_pid.s.x_i0 = -θ) #sign inversion!
                    end

                    v2θ_pid.u.input = EAS_ref - EAS
                    v2θ_pid.u.sat_ext = -elevator_cmd_sat #sign inversion!
                    f_periodic!(v2θ_pid)
                    θ_ref = -v2θ_pid.y.output #sign inversion!

                elseif c2θ_enabled(mode) #θ_ref overridden by c2θ

                    Control.Discrete.assign!(c2θ_pid, c2θ_lookup(EAS, h_e))

                    if mode != mode_prev
                        Control.reset!(c2θ_pid)
                        k_i = c2θ_pid.u.k_i
                        (k_i != 0) && (c2θ_pid.s.x_i0 = θ)
                    end

                    c2θ_pid.u.input = clm_ref - clm
                    c2θ_pid.u.sat_ext = elevator_cmd_sat
                    f_periodic!(c2θ_pid)
                    θ_ref = c2θ_pid.y.output

                else #ModeControlLon.EAS_θ || ModeControlLon.thr_θ

                    #θ_ref unmodified, input value is kept

                end

                k_p_θ = 1.0
                θ_dot_ref = k_p_θ * (θ_ref - θ)
                φ_bnd = clamp(φ, -π/3, π/3)
                q_ref = 1/cos(φ_bnd) * θ_dot_ref + r * tan(φ_bnd)

            end

            q2e_int.u.input = q_ref - q
            q2e_int.u.sat_ext = elevator_cmd_sat
            f_periodic!(q2e_int)

            q2e_pid.u.input = q2e_int.y.output
            q2e_pid.u.sat_ext = elevator_cmd_sat
            f_periodic!(q2e_pid)
            elevator_ref = q2e_pid.y.output

        end

        Control.Discrete.assign!(e2e_lqr, e2e_lookup(EAS, h_e))

        #e2e is purely proportional, so it doesn't need resetting

        e2e_lqr.u.x .= XPitch(vehicle) #state feedback
        e2e_lqr.u.z .= Float64(vehicle.y.systems.act.elevator.cmd) #command variable feedback
        e2e_lqr.u.z_ref .= elevator_ref #command variable reference
        f_periodic!(e2e_lqr)
        elevator_cmd = e2e_lqr.y.output[1]

    end

    mdl.y = ControllerLonY(; mode, throttle_ref, elevator_ref, q_ref, θ_ref, EAS_ref, clm_ref,
        throttle_cmd, elevator_cmd, e2e_lqr = e2e_lqr.y,
        q2e_int = q2e_int.y, q2e_pid = q2e_pid.y, v2θ_pid = v2θ_pid.y,
        c2θ_pid = c2θ_pid.y, v2t_pid = v2t_pid.y)

end


################################################################################
################################# ControllerLat ###################################

@enumx T=ModeControlLatEnum ModeControlLat begin
    direct = 0 #direct aileron & rudder
    sas = 1 #roll & yaw SAS
    p_β = 2 #roll rate + sideslip
    φ_β = 3 #bank angle + sideslip
    χ_β = 4 #course angle + sideslip
end

using .ModeControlLat: ModeControlLatEnum

ar2ar_enabled(mode::ModeControlLatEnum) = (mode === ModeControlLat.sas)

function φβ2ar_enabled(mode::ModeControlLatEnum)
    mode === ModeControlLat.p_β || mode === ModeControlLat.ModeControlLat.φ_β ||
    mode === ModeControlLat.χ_β
end

p2φ_enabled(mode::ModeControlLatEnum) = (mode === ModeControlLat.p_β)
χ2φ_enabled(mode::ModeControlLatEnum) = (mode === ModeControlLat.χ_β)

################################# FieldVectors #################################

#state vector for φβ LQR tracker
@kwdef struct XLat <: FieldVector{8, Float64}
    p::Float64 = 0.0 #roll rate
    r::Float64 = 0.0 #yaw rate
    φ::Float64 = 0.0; #bank angle
    EAS::Float64 = 0.0 #equivalent airspeed
    β::Float64 = 0.0 #AoS
    β_filt::Float64 = 0.0; #filtered AoS (from aerodynamics model)
    ail_p::Float64 = 0.0; #aileron actuator states
    rud_p::Float64 = 0.0; #rudder actuator states
end

function XLat(vehicle::Model{<:C172Y.Vehicle})

    @unpack systems, airflow, kinematics = vehicle.y
    @unpack aero, act = systems
    @unpack e_nb, ω_eb_b = kinematics

    p, _, r = ω_eb_b
    φ = e_nb.φ
    EAS = airflow.EAS
    β = aero.β
    β_filt = aero.β_filt
    ail_p = act.aileron.pos
    rud_p = act.rudder.pos

    XLat(; p, r, φ, EAS, β, β_filt, ail_p, rud_p)

end

@kwdef struct ULat{T} <: FieldVector{2, T}
    aileron_cmd::T = 0.0
    rudder_cmd::T = 0.0
end

@kwdef struct ZLat <: FieldVector{2, Float64}
    φ::Float64 = 0.0
    β::Float64 = 0.0
end


################################## Model ######################################

@kwdef struct ControllerLat{LQ <: LQRTrackerLookup, LP <: PIDLookup} <: AbstractControlChannel
    ar2ar_lookup::LQ = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "ar2ar_lookup.h5"))
    φβ2ar_lookup::LQ = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "φβ2ar_lookup.h5"))
    p2φ_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "p2φ_lookup.h5"))
    χ2φ_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "χ2φ_lookup.h5"))
    ar2ar_lqr::LQRTracker{8, 2, 2, 16, 4} = LQRTracker{8, 2, 2}()
    φβ2ar_lqr::LQRTracker{8, 2, 2, 16, 4} = LQRTracker{8, 2, 2}()
    p2φ_int::Integrator = Integrator()
    p2φ_pid::PID = PID()
    χ2φ_pid::PID = PID()
end

@kwdef mutable struct ControllerLatU
    mode::ModeControlLatEnum = ModeControlLat.direct #lateral control mode
    aileron_ref::Float64 = 0.0 #aileron command reference
    rudder_ref::Float64 = 0.0 #rudder command reference
    p_ref::Float64 = 0.0 #roll rate reference
    β_ref::Float64 = 0.0 #sideslip angle reference
    φ_ref::Float64 = 0.0 #bank angle reference
    χ_ref::Float64 = 0.0 #course angle reference
end

@kwdef struct ControllerLatY
    mode::ModeControlLatEnum = ModeControlLat.direct
    aileron_ref::Float64 = 0.0 #aileron command reference
    rudder_ref::Float64 = 0.0 #rudder command reference
    p_ref::Float64 = 0.0 #roll rate reference
    β_ref::Float64 = 0.0 #sideslip angle reference
    φ_ref::Float64 = 0.0 #bank angle reference
    χ_ref::Float64 = 0.0 #course angle reference
    aileron_cmd::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd::Ranged{Float64, -1., 1.} = 0.0
    ar2ar_lqr::LQRTrackerOutput{8, 2, 2, 16, 4} = LQRTrackerOutput{8, 2, 2}()
    φβ2ar_lqr::LQRTrackerOutput{8, 2, 2, 16, 4} = LQRTrackerOutput{8, 2, 2}()
    p2φ_int::IntegratorOutput = IntegratorOutput()
    p2φ_pid::PIDOutput = PIDOutput()
    χ2φ_pid::PIDOutput = PIDOutput()
end

Modeling.U(::ControllerLat) = ControllerLatU()
Modeling.Y(::ControllerLat) = ControllerLatY()

function Modeling.init!(mdl::Model{<:ControllerLat})

    foreach((mdl.φβ2ar_lqr, mdl.ar2ar_lqr)) do lqr
        lqr.u.bound_lo .= ULat(; aileron_cmd = -1, rudder_cmd = -1)
        lqr.u.bound_hi .= ULat(; aileron_cmd = 1, rudder_cmd = 1)
    end

    #set φ reference limits for the course angle compensator output
    mdl.χ2φ_pid.u.bound_lo = -π/4
    mdl.χ2φ_pid.u.bound_hi = π/4

end

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:ControllerLat},
                        vehicle::Model{<:C172Y.Vehicle})

    @unpack mode, aileron_ref, rudder_ref, p_ref, β_ref, φ_ref, χ_ref = mdl.u
    @unpack ar2ar_lqr, φβ2ar_lqr, p2φ_int, p2φ_pid, χ2φ_pid = mdl.submodels
    @unpack ar2ar_lookup, φβ2ar_lookup, p2φ_lookup, χ2φ_lookup = mdl.constants
    @unpack airflow, kinematics, systems = vehicle.y

    EAS = airflow.EAS
    h_e = Float64(kinematics.h_e)

    @unpack θ, φ = kinematics.e_nb
    p, _, _ = kinematics.ω_wb_b
    β = systems.aero.β
    mode_prev = mdl.y.mode

    aileron_cmd = aileron_ref
    rudder_cmd = rudder_ref

    if ar2ar_enabled(mode) #aileron_cmd and #rudder_cmd overridden by ar2ar

        #no integral control, so no need for reset on mode change
        Control.Discrete.assign!(ar2ar_lqr, ar2ar_lookup(EAS, Float64(h_e)))

        ar2ar_lqr.u.x .= XLat(vehicle)
        ar2ar_lqr.u.z .= ULat(; aileron_cmd = systems.act.aileron.cmd,
                                rudder_cmd = systems.act.rudder.cmd)
        ar2ar_lqr.u.z_ref .= ULat(; aileron_cmd = aileron_ref, rudder_cmd = rudder_ref)
        f_periodic!(ar2ar_lqr)
        @unpack aileron_cmd, rudder_cmd = ULat(ar2ar_lqr.y.output)

    end

    if φβ2ar_enabled(mode) #aileron_cmd and #rudder_cmd overridden by φβ2ar

        u_lat_sat = ULat(φβ2ar_lqr.y.out_sat)

        if p2φ_enabled(mode) #φ_ref overridden by roll rate tracker

            Control.Discrete.assign!(p2φ_pid, p2φ_lookup(EAS, Float64(h_e)))

            if mode != mode_prev
                #our next φ output must match φ reference at φβ2ar input
                Control.reset!(p2φ_int)
                Control.reset!(p2φ_pid)
                k_i = p2φ_pid.u.k_i
                (k_i != 0) && (p2φ_pid.s.x_i0 = ZLat(φβ2ar_lqr.u.z_ref).φ)
            end

            p2φ_int.u.input = p_ref - p
            p2φ_int.u.sat_ext = u_lat_sat.aileron_cmd
            f_periodic!(p2φ_int)

            p2φ_pid.u.input = p2φ_int.y.output
            p2φ_pid.u.sat_ext = u_lat_sat.aileron_cmd
            f_periodic!(p2φ_pid)
            φ_ref = p2φ_pid.y.output

        elseif χ2φ_enabled(mode) #φ_ref overridden by course angle tracker

            Control.Discrete.assign!(χ2φ_pid, χ2φ_lookup(EAS, Float64(h_e)))

            if mode != mode_prev
                #our next φ output must match φ reference at φβ2ar input
                Control.reset!(χ2φ_pid)
                k_i = χ2φ_pid.u.k_i
                (k_i != 0) && (χ2φ_pid.s.x_i0 = ZLat(φβ2ar_lqr.u.z_ref).φ)
            end

            χ = kinematics.χ_gnd
            χ2φ_pid.u.input = wrap_to_π(χ_ref - χ)
            χ2φ_pid.u.sat_ext = u_lat_sat.aileron_cmd
            f_periodic!(χ2φ_pid)
            φ_ref = χ2φ_pid.y.output

        else #ModeControlLat.ModeControlLat.φ_β

            #φ_ref and β_ref directly set by input values, nothing to do here

        end

        Control.Discrete.assign!(φβ2ar_lqr, φβ2ar_lookup(EAS, Float64(h_e)))

        if mode != mode_prev
            Control.reset!(φβ2ar_lqr)
        end

        φβ2ar_lqr.u.x .= XLat(vehicle)
        φβ2ar_lqr.u.z .= ZLat(; φ, β)
        φβ2ar_lqr.u.z_ref .= ZLat(; φ = φ_ref, β = β_ref)
        f_periodic!(φβ2ar_lqr)
        @unpack aileron_cmd, rudder_cmd = ULat(φβ2ar_lqr.y.output)

    end

    mdl.y = ControllerLatY(; mode, aileron_ref, rudder_ref, p_ref, β_ref, φ_ref, χ_ref,
        aileron_cmd, rudder_cmd, ar2ar_lqr = ar2ar_lqr.y, φβ2ar_lqr = φβ2ar_lqr.y,
        p2φ_int = p2φ_int.y, p2φ_pid = p2φ_pid.y, χ2φ_pid = χ2φ_pid.y)

end


################################################################################
############################### Guidance #######################################

#a discrete model implementing a specific longitudinal or lateral guidance mode
abstract type AbstractGuidanceChannel <: ModelDefinition end


################################################################################
############################### Vertical Guidance ##############################

@enumx T=AltTrackingStateEnum AltTrackingState begin
    acquire = 0
    hold = 1
end

using .AltTrackingState: AltTrackingStateEnum

@kwdef struct AltitudeTracking <: AbstractGuidanceChannel
    k_h2c::Float64 = 0.2 #h_err to climb rate gain
    h_thr::Float64 = 10.0 #altitude threshold for state switching
end

@kwdef mutable struct AltitudeTrackingU
    h_ref::Union{HEllip, HOrth} = HEllip(0.0) #altitude reference
end

@kwdef mutable struct AltitudeTrackingS
    state::AltTrackingStateEnum = AltTrackingState.hold
end

@kwdef struct AltitudeTrackingY
    state::AltTrackingStateEnum = AltTrackingState.hold
    mode_ctl_lon::ModeControlLonEnum = ModeControlLon.EAS_clm
    h_err::Float64 = 0.0 #current altitude error
    throttle_ref::Float64 = 0.0
    clm_ref::Float64 = 0.0
end

Modeling.U(::AltitudeTracking) = AltitudeTrackingU()
Modeling.S(::AltitudeTracking) = AltitudeTrackingS()
Modeling.Y(::AltitudeTracking) = AltitudeTrackingY()

altitude_error(h_ref::HEllip, kin_data::KinData) = h_ref - kin_data.h_e
altitude_error(h_ref::HOrth, kin_data::KinData) = h_ref - kin_data.h_o

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:AltitudeTracking},
                        vehicle::Model{<:C172Y.Vehicle})

    @unpack state = mdl.s
    @unpack h_ref = mdl.u
    @unpack k_h2c, h_thr = mdl.constants

    h_err = altitude_error(h_ref, vehicle.y.kinematics)

    if state === AltTrackingState.acquire

        mode_ctl_lon = ModeControlLon.thr_EAS
        throttle_ref = h_err > 0 ? 1.0 : 0.0 #full throttle to climb, idle to descend
        clm_ref = 0.0 #no effect
        (abs(h_err) < h_thr - 1) && (mdl.s.state = AltTrackingState.hold)

    else #AltTrackingState.hold

        mode_ctl_lon = ModeControlLon.EAS_clm
        throttle_ref = 0.0 #no effect
        clm_ref = k_h2c * h_err
        (abs(h_err) > h_thr + 1) && (mdl.s.state = AltTrackingState.acquire)

    end

    mdl.y = AltitudeTrackingY(; state, mode_ctl_lon, h_err, throttle_ref, clm_ref)

end


################################################################################
############################# HorizontalGuidance ###############################

struct Segment
    p1::Geographic{LatLon, Ellipsoidal}
    p2::Geographic{LatLon, Ellipsoidal}
    u_12_h::SVector{3,Float64}
    l_12_h::Float64
    χ_12::Float64
    γ_12::Float64
    q_en::RQuat

    function Segment(ap1::Abstract3DPosition, ap2::Abstract3DPosition)

        p1 = Geographic{LatLon,Ellipsoidal}(ap1)
        p2 = Geographic{LatLon,Ellipsoidal}(ap2)

        r_e1_e = Cartesian(ap1)
        r_e2_e = Cartesian(ap2)
        r_12_e = r_e2_e - r_e1_e

        n_e1 = NVector(p1)
        n_e2 = NVector(p2)
        n_e = NVector(0.5*(n_e1[:] + n_e2[:]))
        q_en = ltf(n_e) #LTF at midpoint

        r_12_n = q_en'(r_12_e)
        r_12_h = SVector{3,Float64}(r_12_n[1], r_12_n[2], 0.0)
        l_12_h = norm(r_12_h)
        @assert l_12_h > 0 "Invalid guidance segment" #must not be vertical or degenerate

        Δh_12 = HEllip(p2) - HEllip(p1)
        u_12_h = r_12_h / l_12_h
        χ_12 = azimuth(u_12_h) #azimuth
        γ_12 = atan(Δh_12, l_12_h) #flight path angle (0 for equal altitude)

        new(p1, p2, u_12_h, l_12_h, χ_12, γ_12, q_en)

    end

end

#construct a Segment from:
#1) initial point p1
#2) an azimuth χ measured on NED(p1)
#3) a distance l along χ measured on the horizontal plane of NED(p1)
#4) an altitude increment Δh

#note: l is the straight-line distance from p1 to p2t, NOT from p1 to p2.
#however, as long as p1 and p2 are relatively close, both their straight-line
#and geodesic distances will be close to l

function Segment(p1::Abstract3DPosition; χ::Real, l::Real, Δh::Real = 0.0)
    geo1 = Geographic{LatLon, Ellipsoidal}(p1)
    q_en1 = ltf(geo1)
    r_12_n1 = SVector{3,Float64}(l*cos(χ), l*sin(χ), 0)
    r_12_e = q_en1(r_12_n1)
    r_e1_e = Cartesian(p1)
    r_e2_e = r_e1_e + r_12_e
    geo2 = Geographic(LatLon(r_e2_e), HEllip(geo1) + Δh)
    Segment(geo1, geo2)
end


Segment() = Segment(Geographic(LatLon(), HEllip()),
                    Geographic(LatLon(0.001, 0), HEllip()))

################################################################################


@kwdef struct SegmentGuidance <: AbstractGuidanceChannel
    Δχ_inf::Ranged{Float64, 0., π/2.} = π/2 #intercept angle for x2d_1p → ∞ (rad)
    e_sf::Float64 = 1000 #cross-track distance scale factor (m)
end

@kwdef mutable struct SegmentGuidanceU
    seg_ref::Segment = Segment()
end

@kwdef struct SegmentGuidanceY
    seg_ref::Segment = Segment()
    l_1b_h::Float64 = 0.0 #along-track horizontal distance
    e_1b_h::Float64 = 0.0 #cross-track horizontal distance, positive right
    χ_ref::Float64 = 0.0 #course angle reference
end

Modeling.U(::SegmentGuidance) = SegmentGuidanceU()
Modeling.Y(::SegmentGuidance) = SegmentGuidanceY()

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:SegmentGuidance},
                        vehicle::Model{<:C172Y.Vehicle})

    @unpack seg_ref = mdl.u
    @unpack Δχ_inf, e_sf = mdl.constants
    @unpack p1, q_en, χ_12, u_12_h = seg_ref

    r_eb_e = Cartesian(vehicle.y.kinematics.r_eb_e)
    r_e1_e = Cartesian(p1)
    r_1b_e = r_eb_e - r_e1_e
    r_1b_n = q_en'(r_1b_e)
    r_1b_h = SVector{3,Float64}(r_1b_n[1], r_1b_n[2], 0)
    l_1b_h = u_12_h ⋅ r_1b_h
    e_1b_h = (u_12_h × r_1b_h)[3]
    Δχ = -Float64(Δχ_inf)/(π/2) * atan(e_1b_h / e_sf)
    χ_ref = wrap_to_π(χ_12 + Δχ)

    mdl.y = SegmentGuidanceY(; seg_ref, l_1b_h, e_1b_h, χ_ref)

end
    # #these go in the longitudinal guidance law
    # r_eb_e = Cartesian(vehicle.y.kinematics.r_eb_e)
    # r_1b_e = r_eb_e - r_e1_e #position vector from p1 to p
    # r_1b_n1 = q_n1e(r_1b_e)
    # r_1b_h = SVector{3,Float64}(r_1b_n1[1], r_1b_n1[2], 0)
    # l_1b_h = u_12_h ⋅ r_1b_h
    # h_ref = h1 + l_1b_h * tan(γ_12) #nominal altitude at our location
    # clm_ff = V_gnd * sin(γ_12) #nominal climb rate at our current V_gnd
    # clm_ref = clm_ff + mdl.k_h2c * (h_ref - h)

    #we need to compute a normalized parameter p_1b.
    #If p_1b<0, we're behind p1, and we set clm_ff = 0, h_nom =
    #if p_1b > 1, we're beyond p2 and we set clm = 0, h_nom = h2.


################################################################################
################################# Controller ##################################

@enumx T=FlightPhaseEnum FlightPhase begin
    gnd = 0
    air = 1
end

@enumx T=ModeGuidanceLonEnum ModeGuidanceLon begin
    off = 0
    alt = 1
end

@enumx T=ModeGuidanceLatEnum ModeGuidanceLat begin
    off = 0
    seg = 1
end

using .FlightPhase: FlightPhaseEnum
using .ModeGuidanceLon: ModeGuidanceLonEnum
using .ModeGuidanceLat: ModeGuidanceLatEnum

################################################################################

@kwdef struct Controller{C1 <: ControllerLon, C2 <: ControllerLat} <: ModelDefinition
    ctl_lon::C1 = ControllerLon()
    ctl_lat::C2 = ControllerLat()
    gdc_lon_alt::AltitudeTracking = AltitudeTracking()
    gdc_lat_seg::SegmentGuidance = SegmentGuidance()
end

#CockpitInputs
@kwdef mutable struct ControllerU
    mode_gdc_lon_req::ModeGuidanceLonEnum = ModeGuidanceLon.off #requested longitudinal guidance mode
    mode_gdc_lat_req::ModeGuidanceLatEnum = ModeGuidanceLat.off #requested lateral guidance mode
    mode_ctl_lon_req::ModeControlLonEnum = ModeControlLon.direct #requested longitudinal control mode
    mode_ctl_lat_req::ModeControlLatEnum = ModeControlLat.direct #requested lateral control mode
    eng_start::Bool = false #passthrough
    eng_stop::Bool = false #passthrough
    mixture::Ranged{Float64, 0., 1.} = 0.5 #passthrough
    flaps::Ranged{Float64, 0., 1.} = 0.0 #passthrough
    brake_left::Ranged{Float64, 0., 1.} = 0.0 #passthrough
    brake_right::Ranged{Float64, 0., 1.} = 0.0 #passthrough
    throttle_axis::Ranged{Float64, 0., 1.} = 0.0
    aileron_axis::Ranged{Float64, -1., 1.} = 0.0
    elevator_axis::Ranged{Float64, -1., 1.} = 0.0
    rudder_axis::Ranged{Float64, -1., 1.} = 0.0
    throttle_offset::Ranged{Float64, 0., 1.} = 0.0
    aileron_offset::Ranged{Float64, -1., 1.} = 0.0
    elevator_offset::Ranged{Float64, -1., 1.} = 0.0
    rudder_offset::Ranged{Float64, -1., 1.} = 0.0
    EAS_ref::Float64 = C172.TrimParameters().EAS #equivalent airspeed reference
    q_ref::Float64 = 0.0 #pitch rate reference
    θ_ref::Float64 = 0.0 #pitch angle reference
    clm_ref::Float64 = 0.0 #climb rate reference
    p_ref::Float64 = 0.0 #roll rate reference
    φ_ref::Float64 = 0.0 #bank angle reference
    χ_ref::Float64 = 0.0 #course angle reference
    β_ref::Float64 = 0.0 #sideslip angle reference
    h_ref::Union{HEllip, HOrth} = HEllip(0.0) #altitude reference
    seg_ref::Segment = Segment() #reference segment
end

@kwdef struct ControllerY
    flight_phase::FlightPhaseEnum = FlightPhase.gnd
    mode_gdc_lon::ModeGuidanceLonEnum = ModeGuidanceLon.off #active longitudinal guidance mode
    mode_gdc_lat::ModeGuidanceLatEnum = ModeGuidanceLat.off #active lateral guidance mode
    mode_ctl_lon::ModeControlLonEnum = ModeControlLon.direct #active longitudinal control mode
    mode_ctl_lat::ModeControlLatEnum = ModeControlLat.direct #active lateral control mode
    ctl_lon::ControllerLonY = ControllerLonY()
    ctl_lat::ControllerLatY = ControllerLatY()
    gdc_lon_alt::AltitudeTrackingY = AltitudeTrackingY()
    gdc_lat_seg::SegmentGuidanceY = SegmentGuidanceY()
end

Modeling.U(::Controller) = ControllerU()
Modeling.Y(::Controller) = ControllerY()


########################### Update Methods #####################################


function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:Controller},
                        vehicle::Model{<:C172Y.Vehicle})
    @unpack mode_gdc_lon_req, mode_gdc_lat_req, mode_ctl_lon_req, mode_ctl_lat_req,
            eng_start, eng_stop, mixture, flaps, brake_left, brake_right,
            throttle_axis, aileron_axis, elevator_axis, rudder_axis,
            throttle_offset, aileron_offset, elevator_offset, rudder_offset,
            q_ref, EAS_ref, θ_ref, clm_ref, p_ref, φ_ref, χ_ref, β_ref, h_ref, seg_ref = mdl.u

    @unpack ctl_lon, ctl_lat, gdc_lon_alt, gdc_lat_seg = mdl.submodels

    throttle_ref = throttle_axis + throttle_offset
    elevator_ref = elevator_axis + elevator_offset
    aileron_ref = aileron_axis + aileron_offset
    rudder_ref = rudder_axis + rudder_offset

    any_wow = any(SVector{3}(leg.strut.wow for leg in vehicle.y.systems.ldg))
    flight_phase = any_wow ? FlightPhase.gnd : FlightPhase.air

    if flight_phase === FlightPhase.gnd

        mode_gdc_lon = ModeGuidanceLon.off
        mode_gdc_lat = ModeGuidanceLat.off
        mode_ctl_lon = ModeControlLon.direct
        mode_ctl_lat = ModeControlLat.direct

    elseif flight_phase === FlightPhase.air

        mode_gdc_lon = mode_gdc_lon_req
        mode_gdc_lat = mode_gdc_lat_req

        if mode_gdc_lon === ModeGuidanceLon.off

            mode_ctl_lon = mode_ctl_lon_req

        else #mode_gdc_lon === ModeGuidanceLon.alt

            gdc_lon_alt.u.h_ref = h_ref
            f_periodic!(gdc_lon_alt, vehicle)

            mode_ctl_lon = gdc_lon_alt.y.mode_ctl_lon
            throttle_ref = gdc_lon_alt.y.throttle_ref
            clm_ref = gdc_lon_alt.y.clm_ref

        end

        if mode_gdc_lat === ModeGuidanceLat.off

            mode_ctl_lat = mode_ctl_lat_req

            #below a v_gnd threshold, override χ mode and revert to φ
            if (mode_ctl_lat === ModeControlLat.χ_β) && (vehicle.y.kinematics.v_gnd < 10.0)
                mode_ctl_lat = ModeControlLat.ModeControlLat.φ_β
            end

        else #mode_gdc_lat === ModeGuidanceLat.seg

            gdc_lat_seg.u.seg_ref = seg_ref
            f_periodic!(gdc_lat_seg, vehicle)
            mode_ctl_lat = ModeControlLat.χ_β
            χ_ref = gdc_lat_seg.y.χ_ref

        end

    end

    ctl_lon.u.mode = mode_ctl_lon
    @pack! ctl_lon.u = throttle_ref, elevator_ref, q_ref, θ_ref, EAS_ref, clm_ref
    f_periodic!(ctl_lon, vehicle)

    ctl_lat.u.mode = mode_ctl_lat
    @pack! ctl_lat.u = aileron_ref, rudder_ref, p_ref, φ_ref, β_ref, χ_ref
    f_periodic!(ctl_lat, vehicle)

    mdl.y = ControllerY(; flight_phase,
        mode_gdc_lon, mode_gdc_lat, mode_ctl_lon, mode_ctl_lat,
        ctl_lon = ctl_lon.y, ctl_lat = ctl_lat.y,
        gdc_lon_alt = gdc_lon_alt.y, gdc_lat_seg = gdc_lat_seg.y)

end

function AircraftBase.assign!(systems::Model{<:C172Y.Systems},
                          mdl::Model{<:Controller})

    @unpack act, pwp, ldg = systems.submodels
    @unpack eng_start, eng_stop, mixture, flaps, brake_left, brake_right = mdl.u
    @unpack throttle_cmd, elevator_cmd = mdl.ctl_lon.y
    @unpack aileron_cmd, rudder_cmd = mdl.ctl_lat.y

    act.throttle.u[] = throttle_cmd
    act.aileron.u[] = aileron_cmd
    act.elevator.u[] = elevator_cmd
    act.rudder.u[] = rudder_cmd
    act.flaps.u[] = flaps
    act.mixture.u[] = mixture
    act.brake_left.u[] = brake_left
    act.brake_right.u[] = brake_right
    pwp.engine.u.start = eng_start
    pwp.engine.u.stop = eng_stop

end


##################################### Tools ####################################

function Modeling.init!(mdl::Model{<:Controller},
                            vehicle::Model{<:C172Y.Vehicle})

    #here we assume that the vehicle's y has already been updated to its trim
    #value by init!(vehicle, params)
    y_act = vehicle.y.systems.act
    @unpack ω_wb_b, v_eb_n, e_nb, χ_gnd, ϕ_λ, h_e = vehicle.y.kinematics
    @unpack EAS = vehicle.y.airflow
    @unpack β = vehicle.y.systems.aero

    #we need to make Controller inputs consistent with the vehicle status, so
    #that trim conditions are preserved upon simulation start when different
    #control modes are selected
    Control.reset!(mdl)

    u = mdl.u
    u.throttle_axis = y_act.throttle.pos
    u.aileron_axis = y_act.aileron.pos
    u.elevator_axis = y_act.elevator.pos
    u.rudder_axis = y_act.rudder.pos
    u.throttle_offset = 0
    u.aileron_offset = 0
    u.elevator_offset = 0
    u.rudder_offset = 0
    u.flaps = y_act.flaps.pos
    u.mixture = y_act.mixture.pos

    u.q_ref = ω_wb_b[2]
    u.θ_ref = e_nb.θ
    u.EAS_ref = EAS
    u.clm_ref = -v_eb_n[3]
    u.p_ref = ω_wb_b[1]
    u.φ_ref = e_nb.φ
    u.β_ref = β
    u.χ_ref = χ_gnd
    u.h_ref = h_e
    u.seg_ref = Segment(Geographic(ϕ_λ, h_e); χ = χ_gnd, l = 1000)

    u.mode_gdc_lon_req = ModeGuidanceLon.off
    u.mode_gdc_lat_req = ModeGuidanceLat.off

    #for the trim condition to be preserved when the simulation is started with
    #sas (rather than direct) modes enabled, we need a post-trim update of the
    #controller with the inner LQR SAS loops enabled. this loads the z_trim,
    #u_trim and x_trim values corresponding to the trim state into the LQR
    #trackers, and runs them once. after this, their outputs will match the
    #actuator commands required by the trim condition. this match is only
    #approximate, because in general, the trim values loaded from the lookup
    #tables will be interpolated, rather than exactly computed at specific
    #controller design points, but it is good enough.
    u.mode_ctl_lon_req = ModeControlLon.sas
    u.mode_ctl_lat_req = ModeControlLat.sas
    f_periodic!(mdl, vehicle)

    u.mode_ctl_lon_req = ModeControlLon.sas
    u.mode_ctl_lat_req = ModeControlLat.ModeControlLat.φ_β
    f_periodic!(mdl, vehicle)

    #restore direct modes
    u.mode_ctl_lon_req = ModeControlLon.direct
    u.mode_ctl_lat_req = ModeControlLat.direct
    f_periodic!(mdl, vehicle)

end


################################################################################
################################## GUI #########################################

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
        Dummy, SameLine, NewLine, IsItemActive, IsItemActivated,
        Separator, Text, Checkbox, RadioButton


function GUI.draw(mdl::Model{<:ControllerLon}, p_open::Ref{Bool} = Ref(true))
    Begin("Longitudinal Control", p_open)
        Text("Mode: $(mdl.y.mode)")
        foreach(keys(mdl.submodels), values(mdl.submodels)) do label, ss
            if CImGui.TreeNode(string(label))
                GUI.draw(ss)
                CImGui.TreePop()
            end
        end
    End()
end

function GUI.draw(mdl::Model{<:ControllerLat}, p_open::Ref{Bool} = Ref(true))
    Begin("Lateral Control", p_open)
        Text("Mode: $(mdl.y.mode)")
        foreach(keys(mdl.submodels), values(mdl.submodels)) do label, ss
            if CImGui.TreeNode(string(label))
                GUI.draw(ss)
                CImGui.TreePop()
            end
        end
    End()
end

function GUI.draw(mdl::Model{<:AltitudeTracking}, p_open::Ref{Bool} = Ref(true))

    @unpack state, mode_ctl_lon, h_err, throttle_ref, clm_ref = mdl.y
    @unpack h_ref = mdl.u
    @unpack h_thr = mdl.constants

    Begin("Altitude Guidance", p_open)

        Text("State: $state")
        Text("Control Mode: $mode_ctl_lon")
        Text("Altitude Ref: $(Float64(h_ref))")
        Text("Altitude Datum: $(typeof(h_ref))")
        Text("Altitude Error: $h_err")
        Text("Altitude Threshold: $h_thr")
        Text("Throttle Setpoint: $throttle_ref")
        Text("Climb Rate Setpoint: $clm_ref")
    End()
end


function GUI.draw!(ctl::Model{<:Controller},
                    vehicle::Model{<:C172Y.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna172Yv1 Flight Control")

    @unpack u, y, Δt, submodels = ctl
    @unpack ctl_lon, ctl_lat, gdc_lon_alt = submodels

    @unpack systems, kinematics, dynamics, airflow = vehicle.y
    @unpack act, pwp, fuel, ldg, aero = systems

    @unpack e_nb, ω_wb_b, n_e, ϕ_λ, h_e, h_o, v_gnd, χ_gnd, γ_gnd, v_eb_n = kinematics
    @unpack α, β = aero
    @unpack CAS, EAS, TAS, T, p, pt = airflow
    @unpack ψ, θ, φ = e_nb
    @unpack ϕ, λ = ϕ_λ

    p, q, r = ω_wb_b
    clm = -v_eb_n[3]
    hog = (ldg.left.strut.Δh + ldg.right.strut.Δh + ldg.nose.strut.Δh) / 3

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
        u.mixture = safe_slider("Mixture", u.mixture, "%.6f"; show_label = true)
        @unpack state, throttle, ω, MAP, τ_shaft, P_shaft, ṁ = pwp.engine

        if CImGui.BeginTable("Engine Data", 4)
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("State: $(state)")
                CImGui.TableNextColumn(); Text(@sprintf("Throttle: %.3f %%", 100*throttle))
                CImGui.TableNextColumn(); Text(@sprintf("Speed: %.3f RPM", Piston.radpersec2RPM(ω)))
                CImGui.TableNextColumn(); Text(@sprintf("Manifold Pressure: %.0f Pa", MAP))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text(@sprintf("Shaft Torque: %.3f Nm", τ_shaft))
                CImGui.TableNextColumn(); Text(@sprintf("Shaft Power: %.3f kW", P_shaft/1e3))
                CImGui.TableNextColumn(); Text(@sprintf("Fuel Flow: %.3f g/s", ṁ*1e3)); SameLine(250)
                CImGui.TableNextColumn(); Text(@sprintf("Remaining Fuel: %.3f kg", fuel.m_avail))
            CImGui.EndTable()
        end
        Separator()
        PopItemWidth()
    end

    ########################## Longitudinal Control ########################

    if CImGui.CollapsingHeader("Longitudinal Control Channel")

        if CImGui.BeginTable("LonCtlModes", 3, CImGui.ImGuiTableFlags_SizingStretchProp )#| CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    Text("Mode")
                CImGui.TableNextColumn();
                    mode_button("Direct##Lon", ModeControlLon.direct, u.mode_ctl_lon_req, y.mode_ctl_lon)
                    IsItemActive() && (u.mode_ctl_lon_req = ModeControlLon.direct); SameLine()
                    mode_button("Pitch SAS", ModeControlLon.sas, u.mode_ctl_lon_req, y.mode_ctl_lon)
                    IsItemActive() && (u.mode_ctl_lon_req = ModeControlLon.sas); SameLine()
                    mode_button("Throttle + Pitch Rate", ModeControlLon.thr_q, u.mode_ctl_lon_req, y.mode_ctl_lon)
                    IsItemActive() && (u.mode_ctl_lon_req = ModeControlLon.thr_q; u.q_ref = 0); SameLine()
                    mode_button("Throttle + Pitch Angle", ModeControlLon.thr_θ, u.mode_ctl_lon_req, y.mode_ctl_lon)
                    IsItemActive() && (u.mode_ctl_lon_req = ModeControlLon.thr_θ; u.θ_ref = θ); SameLine()
                    mode_button("Throttle + EAS", ModeControlLon.thr_EAS, u.mode_ctl_lon_req, y.mode_ctl_lon)
                    IsItemActive() && (u.mode_ctl_lon_req = ModeControlLon.thr_EAS; u.EAS_ref = EAS)
                    mode_button("EAS + Pitch Rate", ModeControlLon.EAS_q, u.mode_ctl_lon_req, y.mode_ctl_lon)
                    IsItemActive() && (u.mode_ctl_lon_req = ModeControlLon.EAS_q; u.q_ref = 0; u.EAS_ref = EAS); SameLine()
                    mode_button("EAS + Pitch Angle", ModeControlLon.EAS_θ, u.mode_ctl_lon_req, y.mode_ctl_lon)
                    IsItemActive() && (u.mode_ctl_lon_req = ModeControlLon.EAS_θ; u.EAS_ref = EAS; u.θ_ref = θ); SameLine()
                    mode_button("EAS + Climb Rate", ModeControlLon.EAS_clm, u.mode_ctl_lon_req, y.mode_ctl_lon)
                    IsItemActive() && (u.mode_ctl_lon_req = ModeControlLon.EAS_clm; u.EAS_ref = EAS; u.clm_ref = clm)
                CImGui.TableNextColumn();
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    AlignTextToFramePadding(); Text("Throttle Axis")
                    AlignTextToFramePadding(); Text("Throttle Offset")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.throttle_axis = safe_slider("Throttle Axis", u.throttle_axis, "%.6f")
                    u.throttle_offset = safe_input("Throttle_Offset", u.throttle_offset, 0.01, 0.1, "%.3f")
                    PopItemWidth()
                # CImGui.TableNextColumn();
                #     Text(@sprintf("%.3f", Float64(u.throttle_axis)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    AlignTextToFramePadding(); Text("Elevator Axis")
                    AlignTextToFramePadding(); Text("Elevator Offset")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.elevator_axis = safe_slider("Elevator Axis", u.elevator_axis, "%.6f")
                    u.elevator_offset = safe_input("Elevator Offset", u.elevator_offset, 0.01, 0.1, "%.3f")
                    PopItemWidth()
                # CImGui.TableNextColumn();
                #     Text(@sprintf("%.3f", Float64(u.elevator_axis)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Pitch Rate (deg/s)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.q_ref = safe_slider("Pitch Rate", rad2deg(u.q_ref), -10, 10, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(q)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Pitch Angle (deg)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.θ_ref = safe_slider("Pitch Angle", rad2deg(u.θ_ref), -15, 15, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(θ)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("EAS (m/s)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.EAS_ref = safe_input("EAS", u.EAS_ref, 0.1, 1.0, "%.3f")
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", EAS))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Climb Rate (m/s)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.clm_ref = safe_input("Climb Rate", u.clm_ref, 0.1, 1.0, "%.3f")
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", clm))
            CImGui.EndTable()
        end

        Separator()

        if CImGui.BeginTable("Actuator Data", 5, CImGui.ImGuiTableFlags_SizingStretchSame)# | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                CImGui.TableNextColumn(); Text("Input")
                CImGui.TableNextColumn(); Text("Command")
                CImGui.TableNextColumn(); Text("Position")
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Throttle")
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.ctl_lon.throttle_ref)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.ctl_lon.throttle_cmd)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(act.throttle.pos)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Elevator")
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.ctl_lon.elevator_ref)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.ctl_lon.elevator_cmd)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(act.elevator.pos)))
            CImGui.EndTable()
        end

        Separator()

    end

    ############################### Lateral Control ############################

    if CImGui.CollapsingHeader("Lateral Control Channel")

        if CImGui.BeginTable("LatCtlModes", 3, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    Text("Mode")
                CImGui.TableNextColumn();
                    mode_button("Direct##Lat", ModeControlLat.direct, u.mode_ctl_lat_req, y.mode_ctl_lat); SameLine()
                    IsItemActive() && (u.mode_ctl_lat_req = ModeControlLat.direct)
                    mode_button("Roll/Yaw SAS", ModeControlLat.sas, u.mode_ctl_lat_req, y.mode_ctl_lat); SameLine()
                    IsItemActive() && (u.mode_ctl_lat_req = ModeControlLat.sas)
                    mode_button("Roll Rate + AoS", ModeControlLat.p_β, u.mode_ctl_lat_req, y.mode_ctl_lat); SameLine()
                    IsItemActive() && (u.mode_ctl_lat_req = ModeControlLat.p_β; u.p_ref = 0; u.β_ref = β)
                    mode_button("Bank Angle + AoS", ModeControlLat.ModeControlLat.φ_β, u.mode_ctl_lat_req, y.mode_ctl_lat); SameLine()
                    IsItemActive() && (u.mode_ctl_lat_req = ModeControlLat.ModeControlLat.φ_β; u.φ_ref = φ; u.β_ref = β)
                    mode_button("Course Angle + AoS", ModeControlLat.χ_β, u.mode_ctl_lat_req, y.mode_ctl_lat); SameLine()
                    IsItemActive() && (u.mode_ctl_lat_req = ModeControlLat.χ_β; u.χ_ref = χ_gnd; u.β_ref = β)
                CImGui.TableNextColumn();
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    AlignTextToFramePadding(); Text("Aileron Axis")
                    AlignTextToFramePadding(); Text("Aileron Offset")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.aileron_axis = safe_slider("Aileron Axis", u.aileron_axis, "%.6f")
                    u.aileron_offset = safe_input("Aileron Offset", u.aileron_offset, 0.01, 0.1, "%.3f")
                    PopItemWidth()
                # CImGui.TableNextColumn();
                #     Text(@sprintf("%.3f", Float64(u.aileron_axis)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    AlignTextToFramePadding(); Text("Rudder Axis")
                    AlignTextToFramePadding(); Text("Rudder Offset")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.rudder_axis = safe_slider("Rudder Axis", u.rudder_axis, "%.6f")
                    u.rudder_offset = safe_input("Rudder Offset", u.rudder_offset, 0.01, 0.1, "%.3f")
                    PopItemWidth()
                # CImGui.TableNextColumn();
                #     Text(@sprintf("%.3f", Float64(u.rudder_axis)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Roll Rate (deg/s)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.p_ref = safe_slider("Roll Rate", rad2deg(u.p_ref), -30, 30, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(p)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Bank Angle (deg)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.φ_ref = safe_slider("Bank Angle", rad2deg(u.φ_ref), -60, 60, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(φ)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Course Angle (deg)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.χ_ref = safe_slider("Course Angle", rad2deg(u.χ_ref), -180, 180, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(χ_gnd)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Sideslip Angle (deg)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.β_ref = safe_slider("Sideslip Angle", rad2deg(u.β_ref), -10, 10, "%.3f") |> deg2rad
                    PopItemWidth()
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(β)))
            CImGui.EndTable()
        end

        Separator()

        if CImGui.BeginTable("Actuator Data", 5, CImGui.ImGuiTableFlags_SizingStretchSame)# | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                CImGui.TableNextColumn(); Text("Input")
                CImGui.TableNextColumn(); Text("Command")
                CImGui.TableNextColumn(); Text("Position")
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Aileron")
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.ctl_lat.aileron_ref)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.ctl_lat.aileron_cmd)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(act.aileron.pos)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Rudder")
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.ctl_lat.rudder_cmd)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.ctl_lat.rudder_ref)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(act.rudder.pos)))
            CImGui.EndTable()
        end

        Separator()

    end

    ############################### Guidance ###################################

    @cstatic h_datum=Int32(0) begin #ellipsoidal: 0, orthometric: 1
    if CImGui.CollapsingHeader("Vertical Guidance")

        if CImGui.BeginTable("Longitudinal Guidance", 3, CImGui.ImGuiTableFlags_SizingStretchProp )#| CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                    Text("Mode")
                CImGui.TableNextColumn();
                    mode_button("Off", ModeGuidanceLon.off, u.mode_gdc_lon_req, y.mode_gdc_lon)
                    IsItemActive() && (u.mode_gdc_lon_req = ModeGuidanceLon.off); SameLine()
                    mode_button("Altitude", ModeGuidanceLon.alt, u.mode_gdc_lon_req, y.mode_gdc_lon)
                    if IsItemActive()
                        u.mode_gdc_lon_req = ModeGuidanceLon.alt
                        u.h_ref = (h_datum == 0 ? h_e : h_o)
                        u.EAS_ref = EAS
                    end
                CImGui.TableNextColumn();
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Altitude (m)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    h_ref_f64 = safe_input("Altitude Setpoint", Float64(u.h_ref), 1, 1.0, "%.3f")
                    u.h_ref = (h_datum == 0 ? HEllip(h_ref_f64) : HOrth(h_ref_f64))
                    PopItemWidth()
                CImGui.TableNextColumn();
                    h_datum == 0 && Text(@sprintf("%.3f", Float64(h_e)))
                    h_datum == 1 && Text(@sprintf("%.3f", Float64(h_o)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                CImGui.TableNextColumn();
                    @c RadioButton("Ellipsoidal", &h_datum, 0); SameLine()
                    IsItemActive() && (u.h_ref = h_e)
                    @c RadioButton("Orthometric", &h_datum, 1)
                    IsItemActive() && (u.h_ref = h_o)
            CImGui.EndTable()
        end #table

        Separator()

    end #header
    end #cstatic


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
                Text(@sprintf("%.3f deg", rad2deg(α)))
                Text("Sideslip Angle"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(β)))
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
                Text(@sprintf("%.3f m | %.3f ft", hog, hog/0.3048))
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
                Text(@sprintf("%.3f g", dynamics.f_c_c[1]/Dynamics.g₀))
                Text("Specific Force (y)"); SameLine(240)
                Text(@sprintf("%.3f g", dynamics.f_c_c[2]/Dynamics.g₀))
                Text("Specific Force (z)"); SameLine(240)
                Text(@sprintf("%.3f g", dynamics.f_c_c[3]/Dynamics.g₀))

            CImGui.EndTable()
        end

        Separator()
    end

    if CImGui.CollapsingHeader("Internals")
        @cstatic c_lon=false c_lat=false c_alt=false begin
            Text("Sampling Period: $Δt")
            @c Checkbox("Longitudinal Control##Internals", &c_lon)
            SameLine()
            @c Checkbox("Lateral Control##Internals", &c_lat)
            SameLine()
            @c Checkbox("Altitude Guidance##Internals", &c_alt)
            c_lon && @c GUI.draw(ctl_lon, &c_lon)
            c_lat && @c GUI.draw(ctl_lat, &c_lat)
            c_alt && @c GUI.draw(gdc_lon_alt, &c_alt)
        end
    end

    End()

end

################################################################################
################################## JSON3  ######################################

#declare ControllerU as mutable
StructTypes.StructType(::Type{ControllerU}) = StructTypes.Mutable()
#replace Greek characters from ControllerU fields in the JSON string
StructTypes.names(::Type{ControllerU}) = ((:χ_ref, :chi_ref), (:θ_ref, :theta_ref),
    (:φ_ref, :phi_ref), (:β_ref, :beta_ref))

#enable JSON parsing of integers as ModeControlLonEnum
StructTypes.StructType(::Type{ModeControlLonEnum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{ModeControlLonEnum}) = Int32 #default enum type
StructTypes.lower(x::ModeControlLonEnum) = Int32(x)

#enable JSON parsing of integers as ModeControlLatEnum
StructTypes.StructType(::Type{ModeControlLatEnum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{ModeControlLatEnum}) = Int32 #default enum type
StructTypes.lower(x::ModeControlLatEnum) = Int32(x)

#enable JSON parsing of integers as ModeGuidanceLonEnum
StructTypes.StructType(::Type{ModeGuidanceLonEnum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{ModeGuidanceLonEnum}) = Int32 #default enum type
StructTypes.lower(x::ModeGuidanceLonEnum) = Int32(x)

#enable JSON parsing of integers as ModeGuidanceLatEnum
StructTypes.StructType(::Type{ModeGuidanceLatEnum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{ModeGuidanceLatEnum}) = Int32 #default enum type
StructTypes.lower(x::ModeGuidanceLatEnum) = Int32(x)

StructTypes.StructType(::Type{Segment}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{Segment}) = NTuple{2, Geographic{LatLon, Ellipsoidal}}
StructTypes.lower(seg::Segment) = (seg.p1, seg.p2)
function StructTypes.construct(::Type{Segment}, ps::NTuple{2, Geographic{LatLon, Ellipsoidal}})
    Segment(ps[1], ps[2])
end

#now we can do:
# JSON3.read(JSON3.write(ControllerU()), ControllerU)
# JSON3.read!(JSON3.write(ControllerU()), ControllerU())

end #module