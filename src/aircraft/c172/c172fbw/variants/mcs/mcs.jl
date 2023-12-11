module C172MCS

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays, HDF5,
    Interpolations, StructArrays

using Flight.FlightCore
using Flight.FlightCore.Utils

using Flight.FlightPhysics
using Flight.FlightComponents
using Flight.FlightComponents.Control.Discrete: Integrator, IntegratorOutput,
    PID, PIDOutput, PIDParams, LQRTracker, LQRTrackerOutput, LQRTrackerParams

using ...C172
using ..C172FBW

export Cessna172MCS

################################################################################
############################# Lookups ##########################################

const PIDLookup = PIDParams{T} where {T <: AbstractInterpolation}

const LQRTrackerLookup = LQRTrackerParams{CB, CF, CI, X, U, Z} where {
    CB <: AbstractInterpolation,
    CF <: AbstractInterpolation,
    CI <: AbstractInterpolation,
    X <: AbstractInterpolation,
    U <: AbstractInterpolation,
    Z <: AbstractInterpolation}

function (lookup::PIDLookup)(EAS::Real, h::Real)
    @unpack k_p, k_i, k_d, τ_f = lookup
    PIDParams(; k_p = k_p(EAS, h),
                k_i = k_i(EAS, h),
                k_d = k_d(EAS, h),
                τ_f = τ_f(EAS, h)
                )
end

function (lookup::LQRTrackerLookup)(EAS::Real, h::Real)
    @unpack C_fbk, C_fwd, C_int, x_trim, u_trim, z_trim = lookup
    LQRTrackerParams(;
        C_fbk = C_fbk(EAS, h),
        C_fwd = C_fwd(EAS, h),
        C_int = C_int(EAS, h),
        x_trim = x_trim(EAS, h),
        u_trim = u_trim(EAS, h),
        z_trim = z_trim(EAS, h),
        )
end

#save and load functions are agnostic about the number and lengths of
#interpolation dimensions
function save_lookup(params::Union{Array{<:LQRTrackerParams, N}, Array{<:PIDParams, N}},
                    bounds::NTuple{N, Tuple{Real, Real}},
                    fname::String = joinpath(@__DIR__, "test.h5")) where {N}

    params_nt = StructArrays.components(StructArray(params))

    fid = h5open(fname, "w")

    create_group(fid, "params")
    foreach(keys(params_nt), values(params_nt)) do k, v
        fid["params"][string(k)] = stack(v)
    end

    fid["bounds"] = stack(bounds) #2xN matrix

    close(fid)
end


function load_pid_lookup(fname::String = joinpath(@__DIR__, "data", "p2φ_lookup.h5"))

    fid = h5open(fname, "r")

    #read fieldnames as ordered in LQRTrackerParams and into an instance
    params_stacked = map(name -> read(fid["params"][string(name)]), fieldnames(PIDParams))
    bounds_stacked = read(fid["bounds"])

    close(fid)

    #arrange bounds back into a Tuple of 2-Tuples
    bounds = mapslices(x->tuple(x...), bounds_stacked, dims = 1) |> vec |> Tuple
    N = length(bounds) #number of interpolation dimensions

    #PID parameters are scalars, so these are already N-dimensional arrays
    params = params_stacked
    @assert allequal(size.(params))
    itp_lengths = size(params[1])

    #define interpolation mode and ranges, handling singleton dimensions
    itp_args = map(bounds, itp_lengths) do b, l
        r = range(b..., length = l)
        (mode, scaling) = length(r) > 1 ? (BSpline(Linear()), r) : (NoInterp(), 1:1)
        return (mode = mode, scaling = scaling)
    end |> collect |> StructArray

    @unpack mode, scaling = itp_args
    interps = [extrapolate(scale(interpolate(p, tuple(mode...)), scaling...), Flat())
                for p in params]

    return PIDParams(interps...)

end


function load_lqr_tracker_lookup(fname::String = joinpath(@__DIR__, "data", "φβ_lookup.h5"))

    fid = h5open(fname, "r")

    #read fieldnames as ordered in LQRTrackerParams and into an instance
    params_stacked = map(name -> read(fid["params"][string(name)]), fieldnames(LQRTrackerParams))
    bounds_stacked = read(fid["bounds"])

    close(fid)

    #arrange bounds back into a Tuple of 2-Tuples
    bounds = mapslices(x->tuple(x...), bounds_stacked, dims = 1) |> vec |> Tuple
    N = length(bounds) #number of interpolation dimensions

    #generate N-dimensional arrays of either SVectors (for x_trim, u_trim and
    #z_trim) or SMatrices (for C_fbk, C_fwd and C_int)
    params_static = map(params_stacked) do p_stacked
        if ndims(p_stacked) == N+1 #vector parameter
            return map(SVector{size(p_stacked)[1]}, eachslice(p_stacked; dims = Tuple(2:N+1)))
        elseif ndims(p_stacked) == N+2 #matrix parameter
            return map(SMatrix{size(p_stacked)[1:2]...}, eachslice(p_stacked; dims = Tuple(3:N+2)))
        else
            error("Number of interpolation dimensions was determined as $N, "*
                "so stacked arrays must be either either $(N+1)-dimensional "*
                "for vector parameters or $(N+2)-dimensional for matrix parameters. "*
                "Stacked array is $(ndims(p_stacked))-dimensional for $p_name")
        end
    end

    #lengths of N interpolation dimensions must be consistent among params
    @assert allequal(size.(params_static))
    itp_lengths = size(params_static[1])

    #define interpolation mode and ranges, handling singleton dimensions
    itp_args = map(bounds, itp_lengths) do b, l
        r = range(b..., length = l)
        (mode, scaling) = length(r) > 1 ? (BSpline(Linear()), r) : (NoInterp(), 1:1)
        return (mode = mode, scaling = scaling)
    end |> collect |> StructArray

    @unpack mode, scaling = itp_args
    interps = [extrapolate(scale(interpolate(p, tuple(mode...)), scaling...), Flat())
                for p in params_static]

    return LQRTrackerParams(interps...)

end

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

#assemble state vector from aircraft physics
function XLon(physics::System{<:C172FBW.Physics})

    @unpack airframe, air, kinematics = physics.y
    @unpack pwp, aero, act = airframe
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

function ULon{Float64}(physics::System{<:C172FBW.Physics})
    @unpack act = physics.y.airframe
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

function ZLonEASClm(physics::System{<:C172FBW.Physics})
    EAS = physics.y.air.EAS
    climb_rate = -physics.y.kinematics.common.v_eOb_n[3]
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

Systems.init(::SystemU, ::LonControl) = LonControlU()
Systems.init(::SystemY, ::LonControl) = LonControlY()

function Systems.init!(sys::System{<:LonControl})
    #set output bounds in all LQR Trackers so they can handle saturation
    foreach((sys.te2te_lqr, sys.vc2te_lqr)) do lqr
        lqr.u.bound_lo .= ULon(; throttle_cmd = 0, elevator_cmd = -1)
        lqr.u.bound_hi .= ULon(; throttle_cmd = 1, elevator_cmd = 1)
    end

end


function Systems.f_disc!(sys::System{<:LonControl},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack mode, throttle_sp, elevator_sp, q_sp, θ_sp, EAS_sp, clm_sp = sys.u
    @unpack te2te_lqr, vc2te_lqr, q2e_int, q2e_pid, v2θ_pid, v2t_pid = sys.subsystems
    @unpack te2te_lookup, vc2te_lookup, q2e_lookup, v2θ_lookup, v2t_lookup = sys.constants

    EAS = physics.y.air.EAS
    h_e = Float64(physics.y.kinematics.h_e)
    _, q, r = physics.y.kinematics.ω_lb_b
    @unpack θ, φ = physics.y.kinematics.e_nb
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
            f_disc!(v2t_pid, Δt)
            throttle_sp = v2t_pid.y.output

        end

        if q2e_enabled(mode) #elevator_sp overridden by q2e

            Control.Discrete.assign!(q2e_pid, q2e_lookup(EAS, h_e))

            if mode != mode_prev
                Systems.reset!(q2e_int)
                Systems.reset!(q2e_pid)
                k_i = q2e_pid.u.k_i
                (k_i != 0) && (q2e_pid.s.x_i0 = Float64(sys.y.elevator_cmd))
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
                    f_disc!(v2θ_pid, Δt)
                    θ_sp = -v2θ_pid.y.output #sign inversion!

                end

                k_p_θ = 1.0
                θ_dot_sp = k_p_θ * (θ_sp - θ)
                φ_bnd = clamp(φ, -π/3, π/3)
                q_sp = 1/cos(φ_bnd) * θ_dot_sp + r * tan(φ_bnd)

            end

            q2e_int.u.input = q_sp - q
            q2e_int.u.sat_ext = u_lon_sat.elevator_cmd
            f_disc!(q2e_int, Δt)

            q2e_pid.u.input = q2e_int.y.output
            q2e_pid.u.sat_ext = u_lon_sat.elevator_cmd
            f_disc!(q2e_pid, Δt)
            elevator_sp = q2e_pid.y.output

        end

        Control.Discrete.assign!(te2te_lqr, te2te_lookup(EAS, h_e))

        #te2te is purely proportional, so it doesn't need resetting

        te2te_lqr.u.x .= XLon(physics) #state feedback
        te2te_lqr.u.z .= ZLonThrEle(physics) #command variable feedback
        te2te_lqr.u.z_sp .= ZLonThrEle(; throttle_cmd = throttle_sp, elevator_cmd = elevator_sp) #command variable setpoint
        f_disc!(te2te_lqr, Δt)
        @unpack throttle_cmd, elevator_cmd = ULon(te2te_lqr.y.output)

    else #mode === lon_EAS_clm #actuation commands computed by vc2te

        Control.Discrete.assign!(vc2te_lqr, vc2te_lookup(EAS, h_e))

        (mode != mode_prev) && Systems.reset!(vc2te_lqr)

        vc2te_lqr.u.x .= XLon(physics) #state feedback
        vc2te_lqr.u.z .= ZLonEASClm(physics) #command variable feedback
        vc2te_lqr.u.z_sp .= ZLonEASClm(; EAS = EAS_sp, climb_rate = clm_sp) #command variable setpoint
        f_disc!(vc2te_lqr, Δt)
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


#assemble state vector from aircraft physics
function XLat(physics::System{<:C172FBW.Physics})

    @unpack airframe, air, kinematics = physics.y
    @unpack aero, act = airframe
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

function ZLatPhiBeta(physics::System{<:C172FBW.Physics})
    φ = physics.y.kinematics.common.e_nb.φ
    β = physics.y.air.β_b
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

Systems.init(::SystemU, ::LatControl) = LatControlU()
Systems.init(::SystemY, ::LatControl) = LatControlY()

function Systems.init!(sys::System{<:LatControl})

    foreach((sys.φβ2ar_lqr, )) do lqr
        lqr.u.bound_lo .= ULat(; aileron_cmd = -1, rudder_cmd = -1)
        lqr.u.bound_hi .= ULat(; aileron_cmd = 1, rudder_cmd = 1)
    end

    #set φ setpoint limits for the course angle compensator output
    sys.χ2φ_pid.u.bound_lo = -π/4
    sys.χ2φ_pid.u.bound_hi = π/4

end

function Systems.f_disc!(sys::System{<:LatControl},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack mode, aileron_sp, rudder_sp, p_sp, β_sp, φ_sp, χ_sp = sys.u
    @unpack φβ2ar_lqr, p2φ_int, p2φ_pid, χ2φ_pid = sys.subsystems
    @unpack φβ2ar_lookup, p2φ_lookup, χ2φ_lookup = sys.constants
    @unpack air, kinematics = physics.y

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
                (k_i != 0) && (p2φ_pid.s.x_i0 = φ)
            end

            p = kinematics.ω_lb_b[1]
            p2φ_int.u.input = p_sp - p
            p2φ_int.u.sat_ext = u_lat_sat.aileron_cmd
            f_disc!(p2φ_int, Δt)

            p2φ_pid.u.input = p2φ_int.y.output
            p2φ_pid.u.sat_ext = u_lat_sat.aileron_cmd
            f_disc!(p2φ_pid, Δt)
            φ_sp = p2φ_pid.y.output

        elseif mode === lat_χ_β

            Control.Discrete.assign!(χ2φ_pid, χ2φ_lookup(EAS, Float64(h_e)))

            if mode != mode_prev
                Systems.reset!(χ2φ_pid)
                k_i = χ2φ_pid.u.k_i
                (k_i != 0) && (χ2φ_pid.s.x_i0 = φ)
            end

            χ = kinematics.χ_gnd
            χ2φ_pid.u.input = wrap_to_π(χ_sp - χ)
            χ2φ_pid.u.sat_ext = u_lat_sat.aileron_cmd
            f_disc!(χ2φ_pid, Δt)
            φ_sp = χ2φ_pid.y.output

        else #mode === lat_φ_β

            #φ_sp and β_sp directly set by input values, nothing to do here

        end

        Control.Discrete.assign!(φβ2ar_lqr, φβ2ar_lookup(EAS, Float64(h_e)))

        (mode != mode_prev) && Systems.reset!(φβ2ar_lqr)

        φβ2ar_lqr.u.x .= XLat(physics)
        φβ2ar_lqr.u.z .= ZLatPhiBeta(physics)
        φβ2ar_lqr.u.z_sp .= ZLatPhiBeta(; φ = φ_sp, β = β_sp)
        f_disc!(φβ2ar_lqr, Δt)
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

Systems.init(::SystemU, ::AltitudeGuidance) = AltitudeGuidanceU()
Systems.init(::SystemS, ::AltitudeGuidance) = AltitudeGuidanceS()
Systems.init(::SystemY, ::AltitudeGuidance) = AltitudeGuidanceY()

get_Δh(h_sp::HEllip, physics::System{<:C172FBW.Physics}) = h_sp - physics.y.kinematics.h_e
get_Δh(h_sp::HOrth, physics::System{<:C172FBW.Physics}) = h_sp - physics.y.kinematics.h_o

function Systems.f_disc!(sys::System{<:AltitudeGuidance},
                        physics::System{<:C172FBW.Physics}, ::Real)

    Δh = get_Δh(sys.u.h_sp, physics)
    clm_sp = sys.constants.k_h2c * Δh
    clm = -physics.y.kinematics.v_eOb_n[3]
    @unpack state, h_thr = sys.s
    # @show Δh
    # @show state

    if state === alt_acquire

        lon_ctl_mode = lon_thr_EAS
        throttle_sp = Δh > 0 ? 1.0 : 0.0 #full throttle to climb, idle to descend
        sys.s.h_thr = abs(Δh)
        (abs(clm_sp) < abs(clm)) && (sys.s.state = alt_hold)

    else #alt_hold

        lon_ctl_mode = lon_EAS_clm
        throttle_sp = 0.0 #no effect, controlled by EAS_clm
        (abs(Δh) > h_thr) && (sys.s.state = alt_acquire)

    end

    sys.y = AltitudeGuidanceY(; state, h_thr, lon_ctl_mode, throttle_sp, clm_sp)

end

@kwdef struct SegmentGuidance <: SystemDefinition end
@kwdef struct SegmentGuidanceY end


################################################################################
#################################### Avionics ##################################

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

@kwdef struct Avionics <: AbstractAvionics
    lon_ctl::LonControl = LonControl()
    lat_ctl::LatControl = LatControl()
    alt_gdc::AltitudeGuidance = AltitudeGuidance()
    seg_gdc::SegmentGuidance = SegmentGuidance()
end

#CockpitInputs
@kwdef mutable struct AvionicsU
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

@kwdef struct AvionicsY
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

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()


########################### Update Methods #####################################


function Systems.f_disc!(avionics::System{<:C172MCS.Avionics},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack eng_start, eng_stop, mixture, flaps, steering, brake_left, brake_right,
            throttle_sp_input, aileron_sp_input, elevator_sp_input, rudder_sp_input,
            throttle_sp_offset, aileron_sp_offset, elevator_sp_offset, rudder_sp_offset,
            vrt_gdc_mode_req, hor_gdc_mode_req, lon_ctl_mode_req, lat_ctl_mode_req,
            q_sp, EAS_sp, θ_sp, clm_sp, p_sp, φ_sp, χ_sp, β_sp, h_sp = avionics.u

    @unpack lon_ctl, lat_ctl, alt_gdc, seg_gdc = avionics.subsystems

    throttle_sp = throttle_sp_input + throttle_sp_offset
    elevator_sp = elevator_sp_input + elevator_sp_offset
    aileron_sp = aileron_sp_input + aileron_sp_offset
    rudder_sp = rudder_sp_input + rudder_sp_offset

    any_wow = any(SVector{3}(leg.strut.wow for leg in physics.y.airframe.ldg))
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
            f_disc!(alt_gdc, physics, Δt)

            lon_ctl_mode = alt_gdc.y.lon_ctl_mode
            throttle_sp = alt_gdc.y.throttle_sp
            clm_sp = alt_gdc.y.clm_sp

        end

        if hor_gdc_mode === hor_gdc_off

            lat_ctl_mode = lat_ctl_mode_req

        else #hor_gdc_mode === hor_gdc_line

            # seg_gdc.u.line_sp = line_sp
            # f_disc!(seg_gdc, physics, Δt)

            # lat_ctl_mode = seg_gdc.y.lat_ctl_mode
            # χ_sp = seg_gdc.y.χ_sp
            # β_sp = 0.0

        end

    end

    lon_ctl.u.mode = lon_ctl_mode
    @pack! lon_ctl.u = throttle_sp, elevator_sp, q_sp, θ_sp, EAS_sp, clm_sp
    f_disc!(lon_ctl, physics, Δt)

    lat_ctl.u.mode = lat_ctl_mode
    @pack! lat_ctl.u = aileron_sp, rudder_sp, p_sp, φ_sp, β_sp, χ_sp
    f_disc!(lat_ctl, physics, Δt)

    avionics.y = AvionicsY(; flight_phase,
        vrt_gdc_mode, hor_gdc_mode, lon_ctl_mode, lat_ctl_mode,
        lon_ctl = lon_ctl.y, lat_ctl = lat_ctl.y,
        alt_gdc = alt_gdc.y, seg_gdc = seg_gdc.y)

    return false

end

function Aircraft.assign!(airframe::System{<:C172FBW.Airframe},
                          avionics::System{Avionics})

    @unpack act, pwp, ldg = airframe.subsystems
    @unpack eng_start, eng_stop, mixture, flaps, steering, brake_left, brake_right = avionics.u
    @unpack throttle_cmd, elevator_cmd = avionics.lon_ctl.y
    @unpack aileron_cmd, rudder_cmd = avionics.lat_ctl.y

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




################################################################################
############################# Cessna172MCS ##################################

const Cessna172MCS{K, T} = C172FBW.Template{K, T, C172MCS.Avionics} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain}

function Cessna172MCS(kinematics = LTF(), terrain = HorizontalTerrain())
    C172FBW.Template(kinematics, terrain, C172MCS.Avionics())
end


##################################### Tools ####################################

function Aircraft.trim!(ac::System{<:Cessna172MCS},
                        trim_params::C172.TrimParameters = C172.TrimParameters())


    result = trim!(ac.physics, trim_params)

    #only ac.physics.y has been updated by the previous call, we need to update
    #ac.y as well
    f_ode!(ac)

    trim_state = result[2]
    @unpack mixture, flaps, EAS = trim_params
    @unpack throttle, aileron, elevator, rudder = trim_state
    @unpack ω_lb_b, v_eOb_n, e_nb, χ_gnd, h_e = ac.y.physics.kinematics
    @unpack β_b = ac.y.physics.air

    #makes Avionics inputs consistent with the trim solution obtained for the
    #aircraft physics, so the trim condition is preserved upon simulation start
    #when the corresponding control modes are selected
    Systems.reset!(ac.avionics)

    u = ac.avionics.u
    u.throttle_sp_input = 0
    u.aileron_sp_input = 0
    u.elevator_sp_input = 0
    u.rudder_sp_input = 0
    u.throttle_sp_offset = throttle
    u.aileron_sp_offset = aileron
    u.elevator_sp_offset = elevator
    u.rudder_sp_offset = rudder
    u.flaps = flaps
    u.mixture = mixture

    u.vrt_gdc_mode_req = vrt_gdc_off
    u.hor_gdc_mode_req = hor_gdc_off
    u.lon_ctl_mode_req = lon_direct
    u.lat_ctl_mode_req = lat_direct

    u.q_sp = ω_lb_b[2]
    u.θ_sp = e_nb.θ
    u.EAS_sp = EAS
    u.clm_sp = -v_eOb_n[3]
    u.p_sp = ω_lb_b[1]
    u.φ_sp = e_nb.φ
    u.β_sp = β_b
    u.χ_sp = χ_gnd
    u.h_sp = h_e

    f_disc!(ac.avionics, ac.physics, 1) #IMPORTANT: update avionics outputs

    return result

end

function Aircraft.linearize!(ac::System{<:Cessna172MCS}, args...; kwargs...)
    linearize!(ac.physics, args...; kwargs...)
end


################################################################################
############################ Joystick Mappings #################################

function IODevices.assign!(sys::System{<:Cessna172MCS}, joystick::Joystick,
                           mapping::InputMapping)
    IODevices.assign!(sys.avionics, joystick, mapping)
end

pitch_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
roll_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
yaw_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)


function IODevices.assign!(sys::System{Avionics},
                           joystick::XBoxController,
                           ::DefaultMapping)

    u = sys.u

    p_sf = 0.5 #roll rate sensitivity
    q_sf = 0.5 #pitch rate sensitivity

    roll_input = get_axis_value(joystick, :right_analog_x) |> roll_curve
    pitch_input = get_axis_value(joystick, :right_analog_y) |> pitch_curve
    yaw_input = get_axis_value(joystick, :left_analog_x) |> yaw_curve

    u.throttle_sp_input += 0.1 * was_released(joystick, :button_Y)
    u.throttle_sp_input -= 0.1 * was_released(joystick, :button_A)

    u.aileron_sp_input = roll_input
    u.elevator_sp_input = pitch_input
    u.rudder_sp_input = yaw_input

    u.p_sp = p_sf * roll_input
    u.q_sp = q_sf * pitch_input

    u.steering = get_axis_value(joystick, :left_analog_x) |> yaw_curve
    u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

    u.aileron_sp_offset -= 0.01 * was_released(joystick, :dpad_left)
    u.aileron_sp_offset += 0.01 * was_released(joystick, :dpad_right)
    u.elevator_sp_offset += 0.01 * was_released(joystick, :dpad_down)
    u.elevator_sp_offset -= 0.01 * was_released(joystick, :dpad_up)

    u.flaps += 0.3333 * was_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

end

function IODevices.assign!(sys::System{Avionics},
                           joystick::T16000M,
                           ::DefaultMapping)

    u = sys.u

    p_sf = 0.5 #roll rate sensitivity
    q_sf = 0.5 #pitch rate sensitivity

    throttle_input = get_axis_value(joystick, :throttle)
    roll_input = get_axis_value(joystick, :stick_x) |> roll_curve
    pitch_input = get_axis_value(joystick, :stick_y) |> pitch_curve
    yaw_input = get_axis_value(joystick, :stick_z) |> yaw_curve

    u.throttle_sp_input = throttle_input
    u.aileron_sp_input = roll_input
    u.elevator_sp_input = pitch_input
    u.rudder_sp_input = yaw_input

    u.p_sp = p_sf * roll_input
    u.q_sp = q_sf * pitch_input

    u.steering = get_axis_value(joystick, :stick_z) |> yaw_curve
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

    u = sys.u

    throttle_input = get_axis_value(joystick, :throttle)
    roll_input = get_axis_value(joystick, :stick_x) |> roll_curve
    pitch_input = get_axis_value(joystick, :stick_y) |> pitch_curve
    yaw_input = get_axis_value(joystick, :stick_z) |> yaw_curve

    u.throttle_sp_input = throttle_input
    u.aileron_sp_input = roll_input
    u.elevator_sp_input = pitch_input
    u.rudder_sp_input = yaw_input

    u.p_sp = p_sf * roll_input
    u.q_sp = q_sf * pitch_input

    u.steering = get_axis_value(joystick, :stick_z) |> yaw_curve
    u.brake_left = is_pressed(joystick, :red_trigger_half)
    u.brake_right = is_pressed(joystick, :red_trigger_half)

    u.aileron_sp_offset -= 2e-4 * is_pressed(joystick, :A3_left)
    u.aileron_sp_offset += 2e-4 * is_pressed(joystick, :A3_right)
    u.elevator_sp_offset += 2e-4 * is_pressed(joystick, :A3_down)
    u.elevator_sp_offset -= 2e-4 * is_pressed(joystick, :A3_up)

    if is_pressed(joystick, :A3_press)
        u.aileron_sp_offset = 0
        u.elevator_sp_offset = 0
    end

    u.flaps += 0.3333 * was_released(joystick, :switch_down)
    u.flaps -= 0.3333 * was_released(joystick, :switch_up)

end

################################################################################
################################## GUI #########################################

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
        Dummy, SameLine, NewLine, IsItemActive, Separator, Text, Checkbox, RadioButton

function mode_button_HSV(button_mode, selected_mode, active_mode)
    if active_mode === button_mode
        return HSV_green
    elseif selected_mode === button_mode
        return HSV_amber
    else
        return HSV_gray
    end
end


function GUI.draw(sys::System{<:LonControl}, p_open::Ref{Bool})
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

function GUI.draw(sys::System{<:LatControl}, p_open::Ref{Bool})
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

function GUI.draw(sys::System{<:AltitudeGuidance}, p_open::Ref{Bool})
    Begin("Altitude Guidance", p_open)
        Text("State: $(sys.y.state)")
        Text("Control Mode: $(sys.y.lon_ctl_mode)")
        Text("Altitude Threshold: $(sys.y.h_thr)")
        Text("Throttle Setpoint: $(sys.y.throttle_sp)")
        Text("Climb Rate Setpoint: $(sys.y.clm_sp)")
    End()
end


function GUI.draw!(avionics::System{<:C172MCS.Avionics},
                    physics::System{<:C172FBW.Physics},
                    label::String = "Cessna 172 MCS Avionics")

    u = avionics.u
    y = avionics.y

    @unpack airframe, kinematics, rigidbody, air = physics.y
    @unpack act, pwp, fuel, ldg = airframe

    @unpack e_nb, ω_lb_b, n_e, ϕ_λ, h_e, h_o, v_gnd, χ_gnd, γ_gnd, v_eOb_n = kinematics
    @unpack CAS, EAS, TAS, α_b, β_b, T, p, pt = air
    @unpack ψ, θ, φ = e_nb
    @unpack ϕ, λ = ϕ_λ

    p, q, r = ω_lb_b
    clm = -v_eOb_n[3]
    Δh = h_o - TerrainData(physics.constants.terrain, n_e).altitude

    Begin(label)

    ############################# Engine Control ###############################

    if pwp.engine.state === Piston.eng_off
        eng_start_HSV = HSV_gray
    elseif pwp.engine.state === Piston.eng_starting
        eng_start_HSV = HSV_amber
    else
        eng_start_HSV = HSV_green
    end

    if CImGui.CollapsingHeader("Engine")
        PushItemWidth(-100)
        dynamic_button("Engine Start", eng_start_HSV, 0.1, 0.1)
        u.eng_start = IsItemActive()
        SameLine()
        dynamic_button("Engine Stop", HSV_gray, (HSV_gray[1], HSV_gray[2], HSV_gray[3] + 0.1), (0.0, 0.8, 0.8))
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

    if CImGui.CollapsingHeader("Vertical Guidance")

        AlignTextToFramePadding(); Text("Mode"); SameLine(160)

        dynamic_button("Off", mode_button_HSV(vrt_gdc_off, u.vrt_gdc_mode_req, y.vrt_gdc_mode), 0.1, 0.1)
        IsItemActive() && (u.vrt_gdc_mode_req = vrt_gdc_off); SameLine()

        @cstatic h_datum=Int32(0) begin
            dynamic_button("Altitude", mode_button_HSV(vrt_gdc_alt, u.vrt_gdc_mode_req, y.vrt_gdc_mode), 0.1, 0.1)
            if IsItemActive()
                u.vrt_gdc_mode_req = vrt_gdc_alt
                u.h_sp = h_datum == 0 ? h_e : h_o
            end
            AlignTextToFramePadding(); Text("Altitude (m)"); SameLine(160)
            CImGui.BeginGroup()
                h_val = safe_input("Altitude Setpoint", Float64(u.h_sp), 1, 1.0, "%.3f")
                SameLine()
                isa(u.h_sp, HEllip) && Text(@sprintf("%.3f", Float64(h_e)))
                isa(u.h_sp, HOrth) && Text(@sprintf("%.3f", Float64(h_o)))
                @c RadioButton("Ellipsoidal", &h_datum, 0); SameLine()
                @c RadioButton("Orthometric", &h_datum, 1)
                u.h_sp = h_datum == 0 ? HEllip(h_val) : HOrth(h_val)
            CImGui.EndGroup()
        end

        Separator()

    end

    # AlignTextToFramePadding()
    # Text("Horizontal Guidance")
    # SameLine(160)
    # dynamic_button("Off", mode_button_HSV(hor_gdc_off, u.hor_gdc_mode_req, y.hor_gdc_mode), 0.1, 0.1)
    # IsItemActive() && (u.hor_gdc_mode_req = hor_gdc_off)

    ########################## Longitudinal Control ########################

    if CImGui.CollapsingHeader("Longitudinal Control")
        AlignTextToFramePadding(); Text("Mode"); SameLine(160)

        CImGui.BeginGroup()

            dynamic_button("Direct##Lon", mode_button_HSV(lon_direct, u.lon_ctl_mode_req, y.lon_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lon_ctl_mode_req = lon_direct); SameLine()

            dynamic_button("Throttle + Pitch SAS", mode_button_HSV(lon_thr_ele, u.lon_ctl_mode_req, y.lon_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lon_ctl_mode_req = lon_thr_ele); SameLine()

            dynamic_button("Throttle + Pitch Rate", mode_button_HSV(lon_thr_q, u.lon_ctl_mode_req, y.lon_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lon_ctl_mode_req = lon_thr_q; u.q_sp = 0); SameLine()

            dynamic_button("Throttle + Pitch Angle", mode_button_HSV(lon_thr_θ, u.lon_ctl_mode_req, y.lon_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lon_ctl_mode_req = lon_thr_θ; u.θ_sp = θ)

            dynamic_button("EAS + Throttle", mode_button_HSV(lon_thr_EAS, u.lon_ctl_mode_req, y.lon_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lon_ctl_mode_req = lon_thr_EAS; u.EAS_sp = EAS); SameLine()

            dynamic_button("EAS + Pitch Rate", mode_button_HSV(lon_EAS_q, u.lon_ctl_mode_req, y.lon_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lon_ctl_mode_req = lon_EAS_q; u.q_sp = 0; u.EAS_sp = EAS); SameLine()

            dynamic_button("EAS + Pitch Angle", mode_button_HSV(lon_EAS_θ, u.lon_ctl_mode_req, y.lon_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lon_ctl_mode_req = lon_EAS_θ; u.EAS_sp = EAS; u.θ_sp = θ); SameLine()

            dynamic_button("EAS + Climb Rate", mode_button_HSV(lon_EAS_clm, u.lon_ctl_mode_req, y.lon_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lon_ctl_mode_req = lon_EAS_clm; u.EAS_sp = EAS; u.clm_sp = clm)

        CImGui.EndGroup()

        AlignTextToFramePadding(); Text("Throttle Input"); SameLine(160)
        u.throttle_sp_input = safe_slider("Throttle Input", u.throttle_sp_input, "%.6f")

        AlignTextToFramePadding(); Text("Throttle Offset"); SameLine(160)
        u.throttle_sp_offset = safe_input("Throttle_Offset", u.throttle_sp_offset, 0.001, 0.1, "%.3f")

        AlignTextToFramePadding(); Text("Elevator Input"); SameLine(160)
        u.elevator_sp_input = safe_slider("Elevator Input", u.elevator_sp_input, "%.6f")

        AlignTextToFramePadding(); Text("Elevator Offset"); SameLine(160)
        u.elevator_sp_offset = safe_input("Elevator Offset", u.elevator_sp_offset, 0.001, 0.1, "%.3f")

        AlignTextToFramePadding(); Text("Pitch Rate (deg/s)"); SameLine(160)
        u.q_sp = safe_input("Pitch Rate", rad2deg(u.q_sp), 0.1, 1.0, "%.3f") |> deg2rad
        SameLine(); Text(@sprintf("%.3f", rad2deg(q)))

        AlignTextToFramePadding(); Text("Pitch Angle (deg)"); SameLine(160)
        u.θ_sp = safe_input("Pitch Angle", rad2deg(u.θ_sp), 0.1, 1.0, "%.3f") |> deg2rad
        SameLine(); Text(@sprintf("%.3f", rad2deg(θ)))

        AlignTextToFramePadding(); Text("EAS (m/s)"); SameLine(160)
        u.EAS_sp = safe_input("EAS", u.EAS_sp, 0.1, 1.0, "%.3f")
        SameLine(); Text(@sprintf("%.3f", EAS))

        AlignTextToFramePadding(); Text("Climb Rate (m/s)"); SameLine(160)
        u.clm_sp = safe_input("Climb Rate", u.clm_sp, 0.1, 1.0, "%.3f")
        SameLine(); Text(@sprintf("%.3f", clm))

        Separator()

    end

    ############################### Lateral Control ########################

    if CImGui.CollapsingHeader("Lateral Control")

        AlignTextToFramePadding(); Text("Mode"); SameLine(160)

        CImGui.BeginGroup()

            dynamic_button("Direct##Lat", mode_button_HSV(lat_direct, u.lat_ctl_mode_req, y.lat_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lat_ctl_mode_req = lat_direct); SameLine()

            dynamic_button("Roll Rate + AoS", mode_button_HSV(lat_p_β, u.lat_ctl_mode_req, y.lat_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lat_ctl_mode_req = lat_p_β; u.p_sp = 0; u.β_sp = 0); SameLine()

            dynamic_button("Bank Angle + AoS", mode_button_HSV(lat_φ_β, u.lat_ctl_mode_req, y.lat_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lat_ctl_mode_req = lat_φ_β; u.φ_sp = φ; u.β_sp = 0); SameLine()

            dynamic_button("Course Angle + AoS", mode_button_HSV(lat_χ_β, u.lat_ctl_mode_req, y.lat_ctl_mode), 0.1, 0.1)
            IsItemActive() && (u.lat_ctl_mode_req = lat_χ_β; u.χ_sp = χ_gnd; u.β_sp = 0)

        CImGui.EndGroup()


        AlignTextToFramePadding(); Text("Aileron Input"); SameLine(160)
        u.aileron_sp_input = safe_slider("Aileron Input", u.aileron_sp_input, "%.6f")

        AlignTextToFramePadding(); Text("Aileron Offset"); SameLine(160)
        u.aileron_sp_offset = safe_input("Aileron Offset", u.aileron_sp_offset, 0.001, 0.1, "%.3f")

        AlignTextToFramePadding(); Text("Rudder Input"); SameLine(160)
        u.rudder_sp_input = safe_slider("Rudder Input", u.rudder_sp_input, "%.6f")

        AlignTextToFramePadding(); Text("Rudder Offset"); SameLine(160)
        u.rudder_sp_offset = safe_input("Rudder Offset", u.rudder_sp_offset, 0.001, 0.1, "%.3f")

        AlignTextToFramePadding(); Text("Roll Rate (deg/s)"); SameLine(160)
        u.p_sp = safe_input("Roll Rate", rad2deg(u.p_sp), 0.1, 1.0, "%.3f") |> deg2rad
        SameLine(); Text(@sprintf("%.3f", rad2deg(p)))

        AlignTextToFramePadding(); Text("Bank Angle (deg)"); SameLine(160)
        u.φ_sp = safe_input("Bank Angle", rad2deg(u.φ_sp), 0.1, 1.0, "%.3f") |> deg2rad
        SameLine(); Text(@sprintf("%.3f", rad2deg(φ)))

        AlignTextToFramePadding(); Text("Course Angle (deg)"); SameLine(160)
        u.χ_sp = safe_input("Course Angle", rad2deg(u.χ_sp), 0.1, 1.0, "%.3f") |> deg2rad
        SameLine(); Text(@sprintf("%.3f", rad2deg(χ_gnd)))

        AlignTextToFramePadding(); Text("Sideslip Angle (deg)"); SameLine(160)
        u.β_sp = safe_input("Sideslip Angle", rad2deg(u.β_sp), 0.1, 1.0, "%.3f") |> deg2rad
        SameLine(); Text(@sprintf("%.3f", rad2deg(β_b)))

        Separator()

    end

    ############################################################################

    if CImGui.CollapsingHeader("Secondary Actuation")

        AlignTextToFramePadding(); Text("Flaps"); SameLine(160)
        u.flaps = safe_slider("Flaps Input", u.flaps, "%.6f")

        AlignTextToFramePadding(); Text("Steering"); SameLine(160)
        u.steering = safe_slider("Steering", u.steering, "%.6f")

        AlignTextToFramePadding(); Text("Left Brake"); SameLine(160)
        u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")

        AlignTextToFramePadding(); Text("Right Brake"); SameLine(160)
        u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")
        Separator()

    end


    ############################################################################

    if CImGui.CollapsingHeader("Primary Actuator Data")
        if CImGui.BeginTable("Primary Actuator Data", 5, CImGui.ImGuiTableFlags_SizingStretchSame | CImGui.ImGuiTableFlags_BordersInner)
            CImGui.TableNextRow()
                CImGui.TableNextColumn();
                CImGui.TableNextColumn(); Text("Throttle")
                CImGui.TableNextColumn(); Text("Elevator")
                CImGui.TableNextColumn(); Text("Aileron")
                CImGui.TableNextColumn(); Text("Rudder")
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Setpoint")
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(y.lon_ctl.throttle_sp))))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(y.lon_ctl.elevator_sp))))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(y.lat_ctl.aileron_sp))))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(y.lat_ctl.rudder_sp))))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Command")
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(y.lon_ctl.throttle_cmd))))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(y.lon_ctl.elevator_cmd))))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(y.lat_ctl.aileron_cmd))))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(y.lat_ctl.rudder_cmd))))
            CImGui.TableNextRow()
                CImGui.TableNextColumn(); Text("Position")
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(act.throttle.pos))))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(act.elevator.pos))))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(act.aileron.pos))))
                CImGui.TableNextColumn(); Text(@sprintf("%.3f", rad2deg(Float64(act.rudder.pos))))
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

                Text("Ground Speed"); SameLine(240)
                Text(@sprintf("%.3f m/s | %.3f kts", v_gnd, Atmosphere.SI2kts(v_gnd)))
                Text("Course Angle"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(χ_gnd)))
                Text("Flight Path Angle"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(γ_gnd)))

            CImGui.TableNextColumn();
                Text("Heading"); SameLine(240);
                Text(@sprintf("%.3f deg", rad2deg(ψ)))
                Text("Inclination"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(θ)))
                Text("Bank"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(φ)))

        Separator()
                Text("Roll Rate"); SameLine(240)
                Text(@sprintf("%.3f deg/s", rad2deg(p)))
                Text("Pitch Rate"); SameLine(240)
                Text(@sprintf("%.3f deg/s", rad2deg(q)))
                Text("Yaw Rate"); SameLine(240)
                Text(@sprintf("%.3f deg/s", rad2deg(r)))
                Separator()

                Text("Specific Force (x)"); SameLine(240)
                Text(@sprintf("%.3f g", rigidbody.f_G_b[1]/RigidBody.g₀))
                Text("Specific Force (y)"); SameLine(240)
                Text(@sprintf("%.3f g", rigidbody.f_G_b[2]/RigidBody.g₀))
                Text("Specific Force (z)"); SameLine(240)
                Text(@sprintf("%.3f g", rigidbody.f_G_b[3]/RigidBody.g₀))
                Separator()

            CImGui.EndTable()
        end

        Separator()
    end

    if CImGui.CollapsingHeader("Internals")
        @cstatic check=false begin
            @c Checkbox("Longitudinal Control##Internals", &check)
            check && @c GUI.draw(avionics.lon_ctl, &check)
        end
        SameLine()
        @cstatic check=false begin
            @c Checkbox("Lateral Control##Internals", &check)
            check && @c GUI.draw(avionics.lat_ctl, &check)
        end
        SameLine()
        @cstatic check=false begin
            @c Checkbox("Altitude Guidance##Internals", &check)
            check && @c GUI.draw(avionics.alt_gdc, &check)
        end
    end

end

# include(joinpath(@__DIR__, "design", "mcs_design.jl")); using .MCSDesign

end #module