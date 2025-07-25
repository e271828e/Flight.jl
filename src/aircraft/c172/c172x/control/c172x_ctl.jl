module C172XControl

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightLib.Control.Discrete: Integrator, IntegratorOutput,
    PID, PIDOutput, PIDParams, LQRTracker, LQRTrackerOutput, LQRTrackerParams,
    PIDLookup, LQRTrackerLookup, load_pid_lookup, load_lqr_tracker_lookup

using ...AircraftBase
using ...C172
using ..C172X


wrap_to_π(x) = x + 2π*floor((π-x)/(2π))

################################################################################
####################### AbstractControlChannel #################################

#a discrete system implementing a specific longitudinal or lateral control mode
abstract type AbstractControlChannel <: ModelDefinition end


################################################################################
################################## LonControl ##################################

@enum LonControlMode begin
    lon_direct = 0 #direct throttle & elevator
    lon_sas = 1 #direct throttle + pitch SAS
    lon_thr_q = 2 #direct throttle + pitch rate
    lon_thr_θ = 3 #direct throttle + pitch angle
    lon_thr_EAS = 4 #direct throttle + EAS
    lon_EAS_q = 5 #EAS + pitch rate
    lon_EAS_θ = 6 #EAS + pitch angle
    lon_EAS_clm = 7 #EAS + climb rate
end

#elevator pitch SAS always enabled except in direct mode
e2e_enabled(mode::LonControlMode) = (mode != lon_direct)

function q2e_enabled(mode::LonControlMode)
    mode === lon_thr_q || mode === lon_thr_θ || mode === lon_thr_EAS ||
    mode === lon_EAS_q || mode === lon_EAS_θ || mode === lon_EAS_clm
end

function θ2q_enabled(mode::LonControlMode)
    mode === lon_thr_θ || mode === lon_thr_EAS || mode === lon_EAS_θ ||
    mode === lon_EAS_clm
end

function v2t_enabled(mode::LonControlMode)
    mode === lon_EAS_q || mode === lon_EAS_θ || mode === lon_EAS_clm
end

c2θ_enabled(mode::LonControlMode) = (mode === lon_EAS_clm)
v2θ_enabled(mode::LonControlMode) = (mode === lon_thr_EAS)

############################## FieldVectors ####################################

#state vector for pitch LQR SAS
@kwdef struct XPitch <: FieldVector{6, Float64}
    q::Float64 = 0.0 #pitch rate
    θ::Float64 = 0.0 #pitch angle
    v_x::Float64 = 0.0 #aerodynamic velocity, x body
    v_z::Float64 = 0.0 #aerodynamic velocity, z body
    α_filt::Float64 = 0.0 #filtered AoA
    ele_p::Float64 = 0.0 #elevator actuator state
end

#assemble state vector from vehicle
function XPitch(vehicle::Model{<:C172X.Vehicle})

    @unpack systems, airflow, kinematics = vehicle.y
    @unpack pwp, aero, act = systems
    @unpack e_nb, ω_eb_b = kinematics

    q = ω_eb_b[2]
    θ = e_nb.θ
    v_x, _, v_z = airflow.v_wb_b
    α_filt = aero.α_filt
    ele_p = act.elevator.pos

    XPitch(; q, θ, v_x, v_z, α_filt, ele_p)

end


################################## Model ######################################

@kwdef struct LonControl{LQ <: LQRTrackerLookup, LP <: PIDLookup} <: AbstractControlChannel
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

@kwdef mutable struct LonControlU
    mode::LonControlMode = lon_direct
    throttle_ref::Float64 = 0.0
    elevator_ref::Float64 = 0.0
    q_ref::Float64 = 0.0
    θ_ref::Float64 = 0.0
    EAS_ref::Float64 = C172.TrimParameters().EAS
    clm_ref::Float64 = 0.0 #climb rate reference
end

@kwdef struct LonControlY
    mode::LonControlMode = lon_direct
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

Modeling.U(::LonControl) = LonControlU()
Modeling.Y(::LonControl) = LonControlY()

function Modeling.init!(mdl::Model{<:LonControl})
    #e2e determines elevator saturation for all pitch control loops (q2e, c2θ,
    #v2θ), so we don't need to set bounds for those
    mdl.e2e_lqr.u.bound_lo .= -1
    mdl.e2e_lqr.u.bound_hi .= 1
    #we do need to set bounds for v2t, so it can handle throttle saturation
    mdl.v2t_pid.u.bound_lo = 0
    mdl.v2t_pid.u.bound_hi = 1
end


function Modeling.f_disc!(::NoScheduling, mdl::Model{<:LonControl},
                        vehicle::Model{<:C172X.Vehicle})

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
        f_disc!(v2t_pid)
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
                    f_disc!(v2θ_pid)
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
                    f_disc!(c2θ_pid)
                    θ_ref = c2θ_pid.y.output

                else #lon_EAS_θ || lon_thr_θ

                    #θ_ref unmodified, input value is kept

                end

                k_p_θ = 1.0
                θ_dot_ref = k_p_θ * (θ_ref - θ)
                φ_bnd = clamp(φ, -π/3, π/3)
                q_ref = 1/cos(φ_bnd) * θ_dot_ref + r * tan(φ_bnd)

            end

            q2e_int.u.input = q_ref - q
            q2e_int.u.sat_ext = elevator_cmd_sat
            f_disc!(q2e_int)

            q2e_pid.u.input = q2e_int.y.output
            q2e_pid.u.sat_ext = elevator_cmd_sat
            f_disc!(q2e_pid)
            elevator_ref = q2e_pid.y.output

        end

        Control.Discrete.assign!(e2e_lqr, e2e_lookup(EAS, h_e))

        #e2e is purely proportional, so it doesn't need resetting

        e2e_lqr.u.x .= XPitch(vehicle) #state feedback
        e2e_lqr.u.z .= Float64(vehicle.y.systems.act.elevator.cmd) #command variable feedback
        e2e_lqr.u.z_ref .= elevator_ref #command variable reference
        f_disc!(e2e_lqr)
        elevator_cmd = e2e_lqr.y.output[1]

    end

    mdl.y = LonControlY(; mode, throttle_ref, elevator_ref, q_ref, θ_ref, EAS_ref, clm_ref,
        throttle_cmd, elevator_cmd, e2e_lqr = e2e_lqr.y,
        q2e_int = q2e_int.y, q2e_pid = q2e_pid.y, v2θ_pid = v2θ_pid.y,
        c2θ_pid = c2θ_pid.y, v2t_pid = v2t_pid.y)

end


################################################################################
################################# LatControl ###################################

@enum LatControlMode begin
    lat_direct = 0 #direct aileron & rudder
    lat_sas = 1 #roll & yaw SAS
    lat_p_β = 2 #roll rate + sideslip
    lat_φ_β = 3 #bank angle + sideslip
    lat_χ_β = 4 #course angle + sideslip
end

ar2ar_enabled(mode::LatControlMode) = (mode === lat_sas)
φβ2ar_enabled(mode::LatControlMode) = (mode === lat_p_β || mode === lat_φ_β || mode === lat_χ_β)
p2φ_enabled(mode::LatControlMode) = (mode === lat_p_β)
χ2φ_enabled(mode::LatControlMode) = (mode === lat_χ_β)

################################# FieldVectors #################################

#state vector for φβ LQR tracker
@kwdef struct XLat <: FieldVector{8, Float64}
    p::Float64 = 0.0 #roll rate
    r::Float64 = 0.0 #yaw rate
    φ::Float64 = 0.0; #bank angle
    v_x::Float64 = 0.0 #aerodynamic velocity, x body
    v_y::Float64 = 0.0 #aerodynamic velocity, y body
    β_filt::Float64 = 0.0; #filtered AoS
    ail_p::Float64 = 0.0; #aileron actuator states
    rud_p::Float64 = 0.0; #rudder actuator states
end

function XLat(vehicle::Model{<:C172X.Vehicle})

    @unpack systems, airflow, kinematics = vehicle.y
    @unpack aero, act = systems
    @unpack e_nb, ω_eb_b = kinematics

    p, _, r = ω_eb_b
    φ = e_nb.φ
    v_x, v_y, _ = airflow.v_wb_b
    β_filt = aero.β_filt
    ail_p = act.aileron.pos
    rud_p = act.rudder.pos

    XLat(; p, r, φ, v_x, v_y, β_filt, ail_p, rud_p)

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

@kwdef struct LatControl{LQ <: LQRTrackerLookup, LP <: PIDLookup} <: AbstractControlChannel
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

@kwdef mutable struct LatControlU
    mode::LatControlMode = lat_direct #lateral control mode
    aileron_ref::Float64 = 0.0 #aileron command reference
    rudder_ref::Float64 = 0.0 #rudder command reference
    p_ref::Float64 = 0.0 #roll rate reference
    β_ref::Float64 = 0.0 #sideslip angle reference
    φ_ref::Float64 = 0.0 #bank angle reference
    χ_ref::Float64 = 0.0 #course angle reference
end

@kwdef struct LatControlY
    mode::LatControlMode = lat_direct
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

Modeling.U(::LatControl) = LatControlU()
Modeling.Y(::LatControl) = LatControlY()

function Modeling.init!(mdl::Model{<:LatControl})

    foreach((mdl.φβ2ar_lqr, mdl.ar2ar_lqr)) do lqr
        lqr.u.bound_lo .= ULat(; aileron_cmd = -1, rudder_cmd = -1)
        lqr.u.bound_hi .= ULat(; aileron_cmd = 1, rudder_cmd = 1)
    end

    #set φ reference limits for the course angle compensator output
    mdl.χ2φ_pid.u.bound_lo = -π/4
    mdl.χ2φ_pid.u.bound_hi = π/4

end

function Modeling.f_disc!(::NoScheduling, mdl::Model{<:LatControl},
                        vehicle::Model{<:C172X.Vehicle})

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
        f_disc!(ar2ar_lqr)
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
            f_disc!(p2φ_int)

            p2φ_pid.u.input = p2φ_int.y.output
            p2φ_pid.u.sat_ext = u_lat_sat.aileron_cmd
            f_disc!(p2φ_pid)
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
            f_disc!(χ2φ_pid)
            φ_ref = χ2φ_pid.y.output

        else #lat_φ_β

            #φ_ref and β_ref directly set by input values, nothing to do here

        end

        Control.Discrete.assign!(φβ2ar_lqr, φβ2ar_lookup(EAS, Float64(h_e)))

        if mode != mode_prev
            Control.reset!(φβ2ar_lqr)
        end

        φβ2ar_lqr.u.x .= XLat(vehicle)
        φβ2ar_lqr.u.z .= ZLat(; φ, β)
        φβ2ar_lqr.u.z_ref .= ZLat(; φ = φ_ref, β = β_ref)
        f_disc!(φβ2ar_lqr)
        @unpack aileron_cmd, rudder_cmd = ULat(φβ2ar_lqr.y.output)

    end

    mdl.y = LatControlY(; mode, aileron_ref, rudder_ref, p_ref, β_ref, φ_ref, χ_ref,
        aileron_cmd, rudder_cmd, ar2ar_lqr = ar2ar_lqr.y, φβ2ar_lqr = φβ2ar_lqr.y,
        p2φ_int = p2φ_int.y, p2φ_pid = p2φ_pid.y, χ2φ_pid = χ2φ_pid.y)

end



################################################################################
############################### Guidance Modes #################################

@enum AltGuidanceState begin
    alt_acquire = 0
    alt_hold = 1
end

@enum AltDatum begin
    ellipsoidal = 0
    orthometric = 1
end

################################################################################

@kwdef struct AltitudeGuidance <: AbstractControlChannel
    k_h2c::Float64 = 0.2 #good margins for the whole envelope
end

@kwdef mutable struct AltitudeGuidanceU
    h_ref::Float64 = 0.0 #altitude reference
    h_datum::AltDatum = ellipsoidal
end

@kwdef mutable struct AltitudeGuidanceS
    state::AltGuidanceState = alt_hold
    h_thr::Float64 = 10.0 #current a
end

@kwdef struct AltitudeGuidanceY
    state::AltGuidanceState = alt_hold
    lon_ctl_mode::LonControlMode = lon_EAS_clm
    Δh::Float64 = 0.0 #current altitude error
    h_thr::Float64 = 0.0 #current altitude switching threshold
    throttle_ref::Float64 = 0.0
    clm_ref::Float64 = 0.0
end

Modeling.U(::AltitudeGuidance) = AltitudeGuidanceU()
Modeling.S(::AltitudeGuidance) = AltitudeGuidanceS()
Modeling.Y(::AltitudeGuidance) = AltitudeGuidanceY()


function Modeling.f_disc!(::NoScheduling, mdl::Model{<:AltitudeGuidance},
                        vehicle::Model{<:C172X.Vehicle})

    @unpack state, h_thr = mdl.s
    @unpack h_ref, h_datum = mdl.u
    @unpack h_e, h_o, v_eb_n = vehicle.y.kinematics

    h = h_datum === ellipsoidal ? Float64(h_e) : Float64(h_o)
    Δh = h_ref - h
    clm_ref = mdl.k_h2c * Δh

    if state === alt_acquire

        lon_ctl_mode = lon_thr_EAS
        throttle_ref = Δh > 0 ? 1.0 : 0.0 #full throttle to climb, idle to descend
        (abs(Δh) + 1 < h_thr) && (mdl.s.state = alt_hold)
        # mdl.s.h_thr = abs(Δh)
        # clm = -vehicle.y.kinematics.v_eb_n[3]
        # (abs(clm_ref) < abs(clm)) && (mdl.s.state = alt_hold)

    else #alt_hold

        lon_ctl_mode = lon_EAS_clm
        throttle_ref = 0.0 #no effect, controlled by EAS_clm
        (abs(Δh) - 1 > h_thr) && (mdl.s.state = alt_acquire)

    end

    mdl.y = AltitudeGuidanceY(; state, lon_ctl_mode, Δh, h_thr, throttle_ref, clm_ref)

end

@kwdef struct SegmentGuidance <: AbstractControlChannel end
@kwdef struct SegmentGuidanceY end


################################################################################
################################# Controller ##################################

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

@kwdef struct Controller{C1 <: LonControl, C2 <: LatControl} <: ModelDefinition
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
    vrt_gdc_mode_req::VerticalGuidanceMode = vrt_gdc_off #requested vertical guidance mode
    hor_gdc_mode_req::HorizontalGuidanceMode = hor_gdc_off #requested horizontal guidance mode
    lon_ctl_mode_req::LonControlMode = lon_direct #requested longitudinal control mode
    lat_ctl_mode_req::LatControlMode = lat_direct #requested lateral control mode
    EAS_ref::Float64 = C172.TrimParameters().EAS #equivalent airspeed reference
    q_ref::Float64 = 0.0 #pitch rate reference
    θ_ref::Float64 = 0.0 #pitch angle reference
    clm_ref::Float64 = 0.0 #climb rate reference
    p_ref::Float64 = 0.0 #roll rate reference
    φ_ref::Float64 = 0.0 #bank angle reference
    χ_ref::Float64 = 0.0 #course angle reference
    β_ref::Float64 = 0.0 #sideslip angle reference
    h_ref::Float64 = 0.0 #altitude reference
    h_datum::AltDatum = ellipsoidal #altitude datum
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

Modeling.U(::Controller) = ControllerU()
Modeling.Y(::Controller) = ControllerY()


########################### Update Methods #####################################


function Modeling.f_disc!(::NoScheduling, mdl::Model{<:Controller},
                        vehicle::Model{<:C172X.Vehicle})

    # println("Hi")

    @unpack eng_start, eng_stop, mixture, flaps, brake_left, brake_right,
            throttle_axis, aileron_axis, elevator_axis, rudder_axis,
            throttle_offset, aileron_offset, elevator_offset, rudder_offset,
            vrt_gdc_mode_req, hor_gdc_mode_req, lon_ctl_mode_req, lat_ctl_mode_req,
            q_ref, EAS_ref, θ_ref, clm_ref, p_ref, φ_ref, χ_ref, β_ref, h_ref, h_datum = mdl.u

    @unpack lon_ctl, lat_ctl, alt_gdc, seg_gdc = mdl.submodels

    throttle_ref = throttle_axis + throttle_offset
    elevator_ref = elevator_axis + elevator_offset
    aileron_ref = aileron_axis + aileron_offset
    rudder_ref = rudder_axis + rudder_offset

    any_wow = any(SVector{3}(leg.strut.wow for leg in vehicle.y.systems.ldg))
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

            alt_gdc.u.h_ref = h_ref
            alt_gdc.u.h_datum = h_datum
            f_disc!(alt_gdc, vehicle)

            lon_ctl_mode = alt_gdc.y.lon_ctl_mode
            throttle_ref = alt_gdc.y.throttle_ref
            clm_ref = alt_gdc.y.clm_ref

        end

        if hor_gdc_mode === hor_gdc_off

            lat_ctl_mode = lat_ctl_mode_req

            #below a v_gnd threshold, override χ mode and revert to φ
            if (lat_ctl_mode === lat_χ_β) && (vehicle.y.kinematics.v_gnd < 10.0)
                lat_ctl_mode = lat_φ_β
            end

        else #hor_gdc_mode === hor_gdc_line

            #remove this when implemented
            lat_ctl_mode = lat_ctl_mode_req
            # seg_gdc.u.line_ref = line_ref
            # f_disc!(seg_gdc, vehicle)

            # lat_ctl_mode = seg_gdc.y.lat_ctl_mode
            # χ_ref = seg_gdc.y.χ_ref
            # β_ref = 0.0

        end

    end

    lon_ctl.u.mode = lon_ctl_mode
    @pack! lon_ctl.u = throttle_ref, elevator_ref, q_ref, θ_ref, EAS_ref, clm_ref
    f_disc!(lon_ctl, vehicle)

    lat_ctl.u.mode = lat_ctl_mode
    @pack! lat_ctl.u = aileron_ref, rudder_ref, p_ref, φ_ref, β_ref, χ_ref
    f_disc!(lat_ctl, vehicle)

    mdl.y = ControllerY(; flight_phase,
        vrt_gdc_mode, hor_gdc_mode, lon_ctl_mode, lat_ctl_mode,
        lon_ctl = lon_ctl.y, lat_ctl = lat_ctl.y,
        alt_gdc = alt_gdc.y, seg_gdc = seg_gdc.y)

end

function AircraftBase.assign!(systems::Model{<:C172X.Systems},
                          mdl::Model{<:Controller})

    @unpack act, pwp, ldg = systems.submodels
    @unpack eng_start, eng_stop, mixture, flaps, brake_left, brake_right = mdl.u
    @unpack throttle_cmd, elevator_cmd = mdl.lon_ctl.y
    @unpack aileron_cmd, rudder_cmd = mdl.lat_ctl.y

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
                            vehicle::Model{<:C172X.Vehicle})

    #here we assume that the vehicle's y has already been updated to its trim
    #value by init!(vehicle, params)
    y_act = vehicle.y.systems.act
    @unpack ω_wb_b, v_eb_n, e_nb, χ_gnd, h_e = vehicle.y.kinematics
    @unpack EAS = vehicle.y.airflow
    @unpack β = vehicle.y.systems.aero

    #we need to make Controller inputs consistent with the vehicle status, so
    #trim conditions are preserved upon simulation start when the corresponding
    #control modes are selected
    Control.reset!(mdl)

    #in a fly-by-wire implementation, it makes more sense to assign the trim
    #values to the inputs rather than the offsets
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
    u.h_ref = Float64(h_e)
    u.h_datum = ellipsoidal

    u.vrt_gdc_mode_req = vrt_gdc_off
    u.hor_gdc_mode_req = hor_gdc_off

    #do an update with the inner SAS loops enabled so that their internal
    #reference values are made consistent with the trim values. this will then make the
    #actuator commands output by the Controller consistent with the trim state
    #values
    u.lon_ctl_mode_req = lon_sas
    u.lat_ctl_mode_req = lat_φ_β
    f_disc!(mdl, vehicle)

    #restore direct modes
    u.lon_ctl_mode_req = lon_direct
    u.lat_ctl_mode_req = lat_direct
    f_disc!(mdl, vehicle)

end


################################################################################
################################## GUI #########################################

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
        Dummy, SameLine, NewLine, IsItemActive, Separator, Text, Checkbox, RadioButton


function GUI.draw(mdl::Model{<:LonControl}, p_open::Ref{Bool} = Ref(true))
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

function GUI.draw(mdl::Model{<:LatControl}, p_open::Ref{Bool} = Ref(true))
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

function GUI.draw(mdl::Model{<:AltitudeGuidance}, p_open::Ref{Bool} = Ref(true))

    @unpack state, lon_ctl_mode, Δh, h_thr, throttle_ref, clm_ref = mdl.y

    Begin("Altitude Guidance", p_open)

        Text("State: $state")
        Text("Control Mode: $lon_ctl_mode")
        Text("Altitude Error: $Δh")
        Text("Altitude Threshold: $h_thr")
        Text("Throttle Setpoint: $throttle_ref")
        Text("Climb Rate Setpoint: $clm_ref")
    End()
end


function GUI.draw!(ctl::Model{<:Controller},
                    vehicle::Model{<:C172X.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna172Xv1 Flight Control")

    @unpack u, y, Δt, submodels = ctl
    @unpack lon_ctl, lat_ctl, alt_gdc = submodels

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
                    mode_button("Direct##Lon", lon_direct, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_direct); SameLine()
                    mode_button("Pitch SAS", lon_sas, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_sas); SameLine()
                    mode_button("Throttle + Pitch Rate", lon_thr_q, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_thr_q; u.q_ref = 0); SameLine()
                    mode_button("Throttle + Pitch Angle", lon_thr_θ, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_thr_θ; u.θ_ref = θ)
                    mode_button("EAS + Throttle", lon_thr_EAS, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_thr_EAS; u.EAS_ref = EAS); SameLine()
                    mode_button("EAS + Pitch Rate", lon_EAS_q, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_EAS_q; u.q_ref = 0; u.EAS_ref = EAS); SameLine()
                    mode_button("EAS + Pitch Angle", lon_EAS_θ, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_EAS_θ; u.EAS_ref = EAS; u.θ_ref = θ); SameLine()
                    mode_button("EAS + Climb Rate", lon_EAS_clm, u.lon_ctl_mode_req, y.lon_ctl_mode)
                    IsItemActive() && (u.lon_ctl_mode_req = lon_EAS_clm; u.EAS_ref = EAS; u.clm_ref = clm)
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
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.lon_ctl.throttle_ref)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.lon_ctl.throttle_cmd)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(act.throttle.pos)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Elevator")
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.lon_ctl.elevator_ref)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.lon_ctl.elevator_cmd)))
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
                    mode_button("Direct##Lat", lat_direct, u.lat_ctl_mode_req, y.lat_ctl_mode); SameLine()
                    IsItemActive() && (u.lat_ctl_mode_req = lat_direct)
                    mode_button("Roll/Yaw SAS", lat_sas, u.lat_ctl_mode_req, y.lat_ctl_mode); SameLine()
                    IsItemActive() && (u.lat_ctl_mode_req = lat_sas)
                    mode_button("Roll Rate + AoS", lat_p_β, u.lat_ctl_mode_req, y.lat_ctl_mode); SameLine()
                    IsItemActive() && (u.lat_ctl_mode_req = lat_p_β; u.p_ref = 0; u.β_ref = β)
                    mode_button("Bank Angle + AoS", lat_φ_β, u.lat_ctl_mode_req, y.lat_ctl_mode); SameLine()
                    IsItemActive() && (u.lat_ctl_mode_req = lat_φ_β; u.φ_ref = φ; u.β_ref = β)
                    mode_button("Course Angle + AoS", lat_χ_β, u.lat_ctl_mode_req, y.lat_ctl_mode); SameLine()
                    IsItemActive() && (u.lat_ctl_mode_req = lat_χ_β; u.χ_ref = χ_gnd; u.β_ref = β)
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
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.lat_ctl.aileron_ref)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.lat_ctl.aileron_cmd)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(act.aileron.pos)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Rudder")
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.lat_ctl.rudder_cmd)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(y.lat_ctl.rudder_ref)))
                CImGui.TableNextColumn(); Text(@sprintf("%.6f", Float64(act.rudder.pos)))
            CImGui.EndTable()
        end

        Separator()

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
                        u.h_ref = (u.h_datum === ellipsoidal ? Float64(h_e) : Float64(h_o))
                        u.EAS_ref = EAS
                    end
                CImGui.TableNextColumn();
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); AlignTextToFramePadding(); Text("Altitude (m)")
                CImGui.TableNextColumn();
                    PushItemWidth(-10)
                    u.h_ref = safe_input("Altitude Setpoint", Float64(u.h_ref), 1, 1.0, "%.3f")
                    PopItemWidth()
                CImGui.TableNextColumn();
                    u.h_datum === ellipsoidal && Text(@sprintf("%.3f", Float64(h_e)))
                    u.h_datum === orthometric && Text(@sprintf("%.3f", Float64(h_o)))
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                CImGui.TableNextColumn();
                    @c RadioButton("Ellipsoidal", &h_datum, 0); SameLine()
                    @c RadioButton("Orthometric", &h_datum, 1)
                    u.h_datum = (h_datum == 0 ? ellipsoidal : orthometric)
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
            c_lon && @c GUI.draw(lon_ctl, &c_lon)
            c_lat && @c GUI.draw(lat_ctl, &c_lat)
            c_alt && @c GUI.draw(alt_gdc, &c_alt)
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

#enable JSON parsing of integers as LonControlMode
StructTypes.StructType(::Type{LonControlMode}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{LonControlMode}) = Int32 #default enum type
StructTypes.lower(x::LonControlMode) = Int32(x)

#enable JSON parsing of integers as LatControlMode
StructTypes.StructType(::Type{LatControlMode}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{LatControlMode}) = Int32 #default enum type
StructTypes.lower(x::LatControlMode) = Int32(x)

#enable JSON parsing of integers as VerticalGuidanceMode
StructTypes.StructType(::Type{VerticalGuidanceMode}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{VerticalGuidanceMode}) = Int32 #default enum type
StructTypes.lower(x::VerticalGuidanceMode) = Int32(x)

#enable JSON parsing of integers as HorizontalGuidanceMode
StructTypes.StructType(::Type{HorizontalGuidanceMode}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{HorizontalGuidanceMode}) = Int32 #default enum type
StructTypes.lower(x::HorizontalGuidanceMode) = Int32(x)

#enable JSON parsing of integers as HorizontalGuidanceMode
StructTypes.StructType(::Type{AltDatum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{AltDatum}) = Int32 #default enum type
StructTypes.lower(x::AltDatum) = Int32(x)

#now we can do:
# JSON3.read(JSON3.write(ControllerU()), ControllerU)
# JSON3.read!(JSON3.write(ControllerU()), ControllerU())

end #module