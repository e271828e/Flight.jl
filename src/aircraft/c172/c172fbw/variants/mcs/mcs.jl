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
    lon_SAS_off = 0 #direct throttle_cmd + elevator_cmd
    lon_thr_ele = 1 #throttle_cmd + elevator_cmd via SAS
    lon_thr_q = 2 #throttle_cmd + pitch rate
    lon_EAS_q = 3 #EAS + pitch rate
    lon_EAS_clm = 4 #EAS + climb rate
    lon_EAS_thr = 5 #EAS + throttle
end

############################## FieldVectors ####################################

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
    throttle_cmd::T = 0.0; elevator_cmd::T = 0.0
end

#command vector for throttle + elevator SAS mode
@kwdef struct ZLonThrEle <: FieldVector{2, Float64}
    throttle_cmd::Float64 = 0.0; elevator_cmd::Float64 = 50.0
end

#command vector for EAS + climb rate mode
@kwdef struct ZLonEASClm <: FieldVector{2, Float64}
    EAS::Float64 = 50.0; climb_rate::Float64 = 0.0
end

#command vector for EAS + throttle mode
@kwdef struct ZLonEASThr <: FieldVector{2, Float64}
    EAS::Float64 = 50.0; throttle_cmd::Float64 = 0.0
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
    thr_v = act.throttle_act.vel
    thr_p = act.throttle_act.pos
    ele_v = act.elevator_act.vel
    ele_p = act.elevator_act.pos

    XLon(; q, θ, v_x, v_z, α_filt, ω_eng, thr_v, thr_p, ele_v, ele_p)

end

function ULon(physics::System{<:C172FBW.Physics})
    @unpack throttle_cmd, elevator_cmd = physics.y.airframe.act
    ULon(; throttle_cmd, elevator_cmd)
end

function ZLonThrEle(physics::System{<:C172FBW.Physics})
    @unpack act = physics.y.airframe
    throttle_cmd = act.throttle_act.cmd
    elevator_cmd = act.elevator_act.cmd
    ZLonThrEle(; throttle_cmd, elevator_cmd)
end

function ZLonEASClm(physics::System{<:C172FBW.Physics})
    EAS = physics.y.air.EAS
    climb_rate = -physics.y.kinematics.common.v_eOb_n[3]
    ZLonEASClm(; EAS, climb_rate)
end

function ZLonEASThr(physics::System{<:C172FBW.Physics})
    EAS = physics.y.air.EAS
    throttle_cmd = physics.y.airframe.act.throttle_act.cmd
    ZLonEASThr(; EAS, throttle_cmd)
end

################################## System ######################################

#since all LQRTrackers have the same dimensions, all LQRTrackerLookup parametric
#types should be the same. to be confirmed. same with PIDLookup types.
@kwdef struct LonControl{LQ <: LQRTrackerLookup, LP <: PIDLookup} <: AbstractControlChannel
    te2te_lookup::LQ = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "te2te_lookup.h5"))
    vc2te_lookup::LQ = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "vc2te_lookup.h5"))
    vt2te_lookup::LQ = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "vt2te_lookup.h5"))
    q2e_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "q2e_lookup.h5"))
    v2t_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "v2t_lookup.h5"))
    te2te_lqr::LQRTracker{10, 2, 2, 20, 4} = LQRTracker{10, 2, 2}()
    vc2te_lqr::LQRTracker{10, 2, 2, 20, 4} = LQRTracker{10, 2, 2}()
    vt2te_lqr::LQRTracker{10, 2, 2, 20, 4} = LQRTracker{10, 2, 2}()
    q2e_int::Integrator = Integrator()
    q2e_pid::PID = PID()
    v2t_pid::PID = PID()
end

@kwdef mutable struct LonControlU
    mode::LonControlMode = lon_SAS_off #selected control mode
    throttle_sp::Float64 = 0.0
    elevator_sp::Float64 = 0.0
    q_sp::Float64 = 0.0
    EAS_sp::Float64 = 50.0 #equivalent airspeed setpoint
    clm_sp::Float64 = 0.0 #climb rate setpoint
end

@kwdef struct LonControlY
    mode::LonControlMode = lon_SAS_off
    throttle_cmd::Ranged{Float64, 0., 1.} = 0.0
    elevator_cmd::Ranged{Float64, -1., 1.} = 0.0
    te2te_lqr::LQRTrackerOutput{10, 2, 2, 20, 4} = LQRTrackerOutput{10, 2, 2}()
    vc2te_lqr::LQRTrackerOutput{10, 2, 2, 20, 4} = LQRTrackerOutput{10, 2, 2}()
    vt2te_lqr::LQRTrackerOutput{10, 2, 2, 20, 4} = LQRTrackerOutput{10, 2, 2}()
    q2e_int::IntegratorOutput = IntegratorOutput()
    q2e_pid::PIDOutput = PIDOutput()
    v2t_pid::PIDOutput = PIDOutput()
end

Systems.init(::SystemU, ::LonControl) = LonControlU()
Systems.init(::SystemY, ::LonControl) = LonControlY()

function Systems.init!(sys::System{<:LonControl})
    #set output bounds in all LQR Trackers so they can handle saturation
    foreach((sys.te2te_lqr, sys.vc2te_lqr, sys.vt2te_lqr)) do lqr
        lqr.u.bound_lo .= ULon(; throttle_cmd = 0, elevator_cmd = -1)
        lqr.u.bound_hi .= ULon(; throttle_cmd = 1, elevator_cmd = 1)
    end

end

function Systems.f_disc!(sys::System{<:LonControl},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack mode, throttle_sp, elevator_sp, q_sp, EAS_sp, clm_sp = sys.u
    @unpack te2te_lqr, vc2te_lqr, vt2te_lqr, q2e_int, q2e_pid, v2t_pid= sys.subsystems
    @unpack te2te_lookup, vc2te_lookup, vt2te_lookup, q2e_lookup, v2t_lookup = sys.constants

    EAS = physics.y.air.EAS
    h_e = Float64(physics.y.kinematics.h_e)
    _, q, _ = physics.y.kinematics.ω_lb_b
    mode_prev = sys.y.mode

    #actuation commands from pilot inceptor set points
    if mode === lon_SAS_off

        throttle_cmd = throttle_sp
        elevator_cmd = elevator_sp

    #actuation commands computed by thr+ele SAS
    elseif mode === lon_thr_ele || mode === lon_thr_q || mode === lon_EAS_q

        u_lon_sat = ULon(te2te_lqr.y.out_sat)

        #elevator_sp computed by q tracker
        if mode === lon_thr_q || mode === lon_EAS_q

            Control.Discrete.assign!(q2e_pid, q2e_lookup(EAS, h_e))

            if mode != mode_prev
                Systems.reset!(q2e_int)
                Systems.reset!(q2e_pid)
                #set the PID integrator's state to the value that yields a
                #steady-state output equal to the previous elevator_cmd
                q2e_pid.s.x_i0 = Float64(sys.y.elevator_cmd) / q2e_pid.u.k_i
            end

            q2e_int.u.input = q_sp - q
            q2e_int.u.sat_ext = u_lon_sat.elevator_cmd
            f_disc!(q2e_int, Δt)

            q2e_pid.u.input = q2e_int.y.output
            q2e_pid.u.sat_ext = u_lon_sat.elevator_cmd
            f_disc!(q2e_pid, Δt)
            elevator_sp = q2e_pid.y.output

            #throttle_sp computed by EAS tracker
            if mode === lon_EAS_q

                Control.Discrete.assign!(v2t_pid, v2t_lookup(EAS, h_e))

                if mode != mode_prev
                    Systems.reset!(v2t_pid)
                    #set the PID integrator's state to the value that yields a
                    #steady-state output equal to the previous throttle_cmd
                    v2t_pid.s.x_i0 = Float64(sys.y.throttle_cmd) / v2t_pid.u.k_i
                end

                v2t_pid.u.input = EAS_sp - EAS
                v2t_pid.u.sat_ext = u_lon_sat.throttle_cmd
                f_disc!(v2t_pid, Δt)
                throttle_sp = v2t_pid.y.output

            end

        end

        #thr+ele SAS is purely proportional, so it doesn't need resetting
        Control.Discrete.assign!(te2te_lqr, te2te_lookup(EAS, h_e))
        te2te_lqr.u.x .= XLon(physics) #state feedback
        te2te_lqr.u.z .= ZLonThrEle(physics) #command variable feedback
        te2te_lqr.u.z_sp .= ZLonThrEle(; throttle_cmd = throttle_sp, elevator_cmd = elevator_sp) #command variable setpoint
        f_disc!(te2te_lqr, Δt)
        @unpack throttle_cmd, elevator_cmd = ULon(te2te_lqr.y.output)

    #actuation commands computed by EAS+climb_rate to thr+ele LQRTracker
    elseif mode === lon_EAS_clm

        if mode != mode_prev
            Systems.reset!(vc2te_lqr)
        end

        Control.Discrete.assign!(vc2te_lqr, vc2te_lookup(EAS, h_e))
        vc2te_lqr.u.x .= XLon(physics) #state feedback
        vc2te_lqr.u.z .= ZLonEASClm(physics) #command variable feedback
        vc2te_lqr.u.z_sp .= ZLonEASClm(; EAS = EAS_sp, climb_rate = clm_sp) #command variable setpoint
        f_disc!(vc2te_lqr, Δt)
        @unpack throttle_cmd, elevator_cmd = ULon(vc2te_lqr.y.output)

    #actuation commands computed by EAS+throttle to thr+ele LQRTracker
    else #mode == lon_EAS_thr

        if mode != mode_prev
            Systems.reset!(vt2te_lookup)
        end

        Control.Discrete.assign!(vt2te_lqr, vt2te_lookup(EAS, h_e))
        vt2te_lqr.u.x .= XLon(physics) #state feedback
        vt2te_lqr.u.z .= ZLonEASThr(physics) #command variable feedback
        vt2te_lqr.u.z_sp .= ZLonEASThr(; EAS = EAS_sp, throttle_cmd = throttle_sp) #command variable setpoint
        f_disc!(vt2te_lqr, Δt)
        @unpack throttle_cmd, elevator_cmd = ULon(vt2te_lqr.y.output)

    end

    sys.y = LonControlY(; mode, throttle_cmd, elevator_cmd,
                        te2te_lqr = te2te_lqr.y,
                        vc2te_lqr = vc2te_lqr.y,
                        vt2te_lqr = vt2te_lqr.y,
                        q2e_int = q2e_int.y,
                        q2e_pid = q2e_pid.y,
                        v2t_pid = v2t_pid.y)

end


################################################################################
################################# LatControl ###################################

@enum LatControlMode begin
    lat_SAS_off = 0 #direct aileron_cmd + rudder_cmd
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
    ail_v = act.aileron_act.vel
    ail_p = act.aileron_act.pos
    rud_v = act.rudder_act.vel
    rud_p = act.rudder_act.pos

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
    mode::LatControlMode = lat_SAS_off #lateral control mode
    aileron_sp::Float64 = 0.0 #aileron command setpoint
    rudder_sp::Float64 = 0.0 #rudder command setpoint
    p_sp::Float64 = 0.0 #roll rate setpoint
    β_sp::Float64 = 0.0 #sideslip angle setpoint
    φ_sp::Float64 = 0.0 #bank angle setpoint
    χ_sp::Float64 = 0.0 #course angle setpoint
end

@kwdef struct LatControlY
    mode::LatControlMode = lat_SAS_off
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

    #set φ demand limits for the course angle compensator output
    sys.χ2φ_pid.u.bound_lo = -π/4
    sys.χ2φ_pid.u.bound_hi = π/4

end

function Systems.f_disc!(sys::System{<:LatControl},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack mode, aileron_sp, rudder_sp, p_sp, φ_sp, χ_sp = sys.u
    @unpack φβ2ar_lqr, p2φ_int, p2φ_pid, χ2φ_pid = sys.subsystems
    @unpack φβ2ar_lookup, p2φ_lookup, χ2φ_lookup = sys.constants

    EAS = physics.y.air.EAS
    h_e = Float64(physics.y.kinematics.h_e)
    @unpack φ = physics.y.kinematics.e_nb
    mode_prev = sys.y.mode

    if mode === lat_SAS_off
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
                #set the PID integrator's state to the value that yields a
                #steady-state output equal to the current φ
                p2φ_pid.s.x_i0 = φ / p2φ_pid.u.k_i
            end

            p2φ_int.u.input = p_sp - p
            p2φ_int.u.sat_ext = u_lat_sat.aileron_cmd
            f_disc!(p2φ_int, Δt)

            p2φ_pid.u.input = p2φ_int.y.output
            p2φ_pid.u.sat_ext = u_lat_sat.aileron_cmd
            f_disc!(p2φ_pid, Δt)
            φ_sp = p2φ_pid.y.output

        elseif mode === lat_χ_β

            if mode != mode_prev
                Systems.reset!(χ2φ_pid)
                # χ2φ_pid.s.x_i0 = φ / χ2φ_pid.u.k_i #not really essential
            end

            Control.Discrete.assign!(χ2φ_pid, χ2φ_lookup(EAS, Float64(h_e)))
            χ2φ_pid.u.input = wrap_to_π(χ_sp - χ)
            χ2φ_pid.u.sat_ext = u_lat_sat.aileron_cmd
            f_disc!(χ2φ_pid, Δt)
            φ_sp = χ2φ_pid.y.output

        else #mode === lat_φ_β

            #φ_sp and β_sp directly set by input values, nothing to do here

        end

        if mode != mode_prev
            Systems.reset!(φβ2ar_lqr)
        end

        Control.Discrete.assign!(φβ2ar_lqr, φβ2ar_lookup(EAS, Float64(h_e)))
        φβ2ar_lqr.u.x .= XLat(physics)
        φβ2ar_lqr.u.z .= ZLatPhiBeta(physics)
        φβ2ar_lqr.u.z_sp .= ZLatPhiBeta(; φ = φ_sp, β = β_sp)
        f_disc!(φβ2ar_lqr, Δt)
        @unpack aileron_cmd, rudder_cmd = ULat(φβ2ar_lqr.y.output)

    end

    sys.y = LatControlY(; mode, aileron_cmd, rudder_cmd,
                        φβ2ar_lqr = φβ2ar_lqr.y,
                        p2φ_int = p2φ_int.y,
                        p2φ_pid = p2φ_pid.y,
                        χ2φ_pid = χ2φ_pid.y)

end



################################################################################
############################### Guidance Modes #################################

@kwdef struct AltitudeGuidance <: SystemDefinition end

@kwdef struct AltitudeGuidanceU
    h_sp::Float64 = 0.0 #climb rate setpoint
end

@kwdef struct AltitudeGuidanceY end

@kwdef struct SegmentGuidance <: SystemDefinition end
@kwdef struct SegmentGuidanceY end


################################################################################
################################# FlightGuidance ###############################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@enum VerticalGuidanceMode begin
    ver_gdc_off = 0
    ver_gdc_alt = 1
end

@enum HorizontalGuidanceMode begin
    hor_gdc_off = 0
    hor_gdc_line = 1
end

@kwdef struct FlightGuidance <: SystemDefinition
    lon_ctl::LonControl = LonControl()
    lat_ctl::LatControl = LatControl()
    alt_gdc::AltitudeGuidance = AltitudeGuidance()
    seg_gdc::SegmentGuidance = SegmentGuidance()
end

@kwdef mutable struct FlightGuidanceInputs
    flight_phase::FlightPhase = phase_gnd #slave to sensor might go here
    ver_gdc_mode_req::VerticalGuidanceMode = ver_gdc_off #requested vertical guidance mode
    hor_gdc_mode_req::HorizontalGuidanceMode = hor_gdc_off #requested horizontal guidance mode
    lon_ctl_mode_req::LonControlMode = lon_SAS_off #requested longitudinal control mode
    lat_ctl_mode_req::LatControlMode = lat_SAS_off #requested lateral control mode
    throttle_sp::Float64 = 0.0 #throttle command setpoint
    elevator_sp::Float64 = 0.0 #elevator command setpoint
    aileron_sp::Float64 = 0.0 #aileron command setpoint
    rudder_sp::Float64 = 0.0 #rudder command setpoint
    p_sp::Float64 = 0.0 #roll rate setpoint
    q_sp::Float64 = 0.0 #pitch rate setpoint
    β_sp::Float64 = 0.0 #sideslip angle setpoint
    EAS_sp::Float64 = 50.0 #equivalent airspeed setpoint
    clm_sp::Float64 = 0.0 #climb rate setpoint
    φ_sp::Float64 = 0.0 #bank angle demand
    χ_sp::Float64 = 0.0 #course angle demand
    h_sp::Union{HEllip, HOrth} = HEllip(0.0) #altitude setpoint
    # seg_sp::Segment = 0.0 #line segment setpoint
end

@kwdef struct FlightGuidanceOutputs
    ver_gdc_mode::VerticalGuidanceMode = ver_gdc_off #active vertical guidance mode
    hor_gdc_mode::HorizontalGuidanceMode = hor_gdc_off #active horizontal guidance mode
    lon_ctl_mode::LonControlMode = lon_SAS_off #active longitudinal control mode
    lat_ctl_mode::LatControlMode = lat_SAS_off #active lateral control mode
    throttle_cmd::Ranged{Float64, 0., 1.} = 0.0 #throttle command
    elevator_cmd::Ranged{Float64, -1., 1.} = 0.0 #elevator command
    aileron_cmd::Ranged{Float64, -1., 1.} = 0.0 #aileron command
    rudder_cmd::Ranged{Float64, -1., 1.} = 0.0 #rudder command
    lon_ctl::LonControlY = LonControlY()
    lat_ctl::LatControlY = LatControlY()
    alt_gdc::AltitudeGuidanceY = AltitudeGuidanceY()
    seg_gdc::SegmentGuidanceY = SegmentGuidanceY()
end

Systems.init(::SystemU, ::FlightGuidance) = FlightGuidanceInputs()
Systems.init(::SystemY, ::FlightGuidance) = FlightGuidanceOutputs()

function Systems.f_disc!(sys::System{<:FlightGuidance},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack flight_phase, ver_gdc_mode_req, hor_gdc_mode_req,
            lon_ctl_mode_req, lat_ctl_mode_req,
            throttle_sp, elevator_sp, aileron_sp, rudder_sp,
            p_sp, q_sp, β_sp, EAS_sp, clm_sp, φ_sp, χ_sp, h_sp = sys.u
    @unpack lon_ctl, lat_ctl, alt_gdc, seg_gdc = sys.subsystems

    if flight_phase === phase_gnd

        ver_gdc_mode = ver_gdc_off
        hor_gdc_mode = hor_gdc_off
        lon_ctl_mode = lon_SAS_off
        lat_ctl_mode = lat_SAS_off

    elseif flight_phase === phase_air

        ver_gdc_mode = ver_gdc_mode_req
        hor_gdc_mode = hor_gdc_mode_req

        if ver_gdc_mode === ver_gdc_off

            lon_ctl_mode = lon_ctl_mode_req

        else #ver_gdc_mode === ver_gdc_alt

            # alt_gdc.u.h_sp = h_sp
            # f_disc!(alt_gdc, physics, Δt)

            # lon_ctl_mode = alt_gdc.y.lon_ctl_mode
            # EAS_sp = alt_gdc.y.EAS_sp
            # clm_sp = alt_gdc.y.clm_sp

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
    @pack! lon_ctl.u = throttle_sp, elevator_sp, q_sp, EAS_sp, clm_sp
    f_disc!(lon_ctl, physics, Δt)
    @unpack throttle_cmd, elevator_cmd = lon_ctl.y

    lat_ctl.u.mode = lat_ctl_mode
    @pack! lat_ctl.u = aileron_sp, rudder_sp, p_sp, φ_sp, β_sp, χ_sp
    f_disc!(lat_ctl, physics, Δt)
    @unpack aileron_cmd, rudder_cmd = lat_ctl.y

    sys.y = FlightGuidanceOutputs(;
        ver_gdc_mode, hor_gdc_mode, lon_ctl_mode, lat_ctl_mode,
        throttle_cmd, elevator_cmd, aileron_cmd, rudder_cmd,
        lon_ctl = lon_ctl.y, lat_ctl = lat_ctl.y,
        alt_gdc = alt_gdc.y, seg_gdc = seg_gdc.y)

end


################################################################################

@kwdef struct Avionics <: AbstractAvionics
    flight::FlightGuidance = FlightGuidance()
end

#CockpitInputs
@kwdef mutable struct AvionicsU
    eng_start::Bool = false
    eng_stop::Bool = false
    mixture::Ranged{Float64, 0., 1.} = 0.5
    throttle_input::Ranged{Float64, 0., 1.} = 0.0 #sets throttle_sp
    roll_input::Ranged{Float64, -1., 1.} = 0.0 #sets aileron_sp or p_sp
    pitch_input::Ranged{Float64, -1., 1.} = 0.0 #sets elevator_sp or q_sp
    yaw_input::Ranged{Float64, -1., 1.} = 0.0 #sets rudder_sp or β_sp
    throttle_sp_offset::Ranged{Float64, 0., 1.} = 0.0 #for direct throttle control only
    aileron_sp_offset::Ranged{Float64, -1., 1.} = 0.0 #for direct aileron control only
    elevator_sp_offset::Ranged{Float64, -1., 1.} = 0.0 #for direct elevator control only
    rudder_sp_offset::Ranged{Float64, -1., 1.} = 0.0 #for direct rudder control only
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
    p_sf::Float64 = 1.0 #roll input to roll rate scale factor (0.1)
    q_sf::Float64 = 1.0 #pitch input to pitch rate scale factor (0.1)
    β_sf::Float64 = 1.0 #yaw input to β_sp scale factor (0.1)
    ver_gdc_mode_req::VerticalGuidanceMode = ver_gdc_off #requested vertical guidance mode
    hor_gdc_mode_req::HorizontalGuidanceMode = hor_gdc_off #requested horizontal guidance mode
    lon_ctl_mode_req::LonControlMode = lon_SAS_off #requested longitudinal control mode
    lat_ctl_mode_req::LatControlMode = lat_SAS_off #requested lateral control mode
    EAS_sp::Float64 = 50.0 #equivalent airspeed setpoint
    clm_sp::Float64 = 0.0 #climb rate setpoint
    β_sp::Float64 = 0.0 #sideslip angle setpoint
    φ_sp::Float64 = 0.0 #bank angle demand
    χ_sp::Float64 = 0.0 #course angle demand
    h_sp::Union{HEllip, HOrth} = HEllip(0.0) #altitude setpoint
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
    flight_phase::FlightPhase = phase_gnd
    actuation::ActuationCommands = ActuationCommands()
    flight::FlightGuidanceOutputs = FlightGuidanceOutputs()
end

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()


########################### Update Methods #####################################


function Systems.f_disc!(avionics::System{<:C172MCS.Avionics},
                        physics::System{<:C172FBW.Physics}, Δt::Real)

    @unpack subsystems, u = avionics
    flight = avionics.flight

    @unpack airframe, air = physics.y

    @unpack eng_start, eng_stop, mixture, throttle_input, roll_input, pitch_input, yaw_input,
            throttle_sp_offset, aileron_sp_offset, elevator_sp_offset, rudder_sp_offset,
            flaps, brake_left, brake_right, ver_gdc_mode_req, hor_gdc_mode_req,
            lon_ctl_mode_req, lat_ctl_mode_req, p_sf, q_sf, β_sf, EAS_sp, clm_sp,
            β_sp, φ_sp, χ_sp, h_sp = avionics.u

    throttle_sp = throttle_input + throttle_sp_offset
    elevator_sp = pitch_input + elevator_sp_offset
    aileron_sp = roll_input + aileron_sp_offset
    rudder_sp = yaw_input + rudder_sp_offset

    p_sp = p_sf * Float64(roll_input)
    q_sp = q_sf * Float64(pitch_input)
    β_sp = β_sf * Float64(yaw_input)

    any_wow = any(SVector{3}(leg.strut.wow for leg in airframe.ldg))
    flight_phase = any_wow ? phase_gnd : phase_air

    @pack! flight.u = flight_phase, throttle_sp, elevator_sp, #computed
        aileron_sp, rudder_sp, p_sp, q_sp, β_sp, #computed
        ver_gdc_mode_req, hor_gdc_mode_req, #forwarded
        lon_ctl_mode_req, lat_ctl_mode_req, EAS_sp, clm_sp, φ_sp, χ_sp, h_sp #forwarded

    #for now, brakes, nws, flaps and mixture do not go through FlightGuidance
    f_disc!(flight, physics, Δt)
    @unpack throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd = flight.y

    actuation = ActuationCommands(; eng_start, eng_stop, mixture,
                throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd,
                flaps, brake_left, brake_right)

    avionics.y = AvionicsY(; flight_phase, actuation, flight = flight.y)

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

function GUI.draw!(avionics::System{<:C172MCS.Avionics},
                    physics::System{<:C172FBW.Physics},
                    label::String = "Cessna 172 SS CAS Avionics")
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
    trim_state = result[2]

    #makes Avionics inputs consistent with the trim solution obtained for the
    #aircraft physics so the trim condition is preserved upon simulation start
    #when direct control modes are selected
    @unpack mixture, flaps = trim_params
    @unpack throttle, aileron, elevator, rudder = trim_state

    Systems.reset!(ac.avionics)

    u = ac.avionics.u
    u.throttle_input = 0
    u.roll_input = 0
    u.pitch_input = 0
    u.yaw_input = 0
    u.throttle_sp_offset = throttle
    u.aileron_sp_offset = aileron
    u.elevator_sp_offset = elevator
    u.rudder_sp_offset = rudder
    u.flaps = flaps
    u.mixture = mixture

    u.ver_gdc_mode_req = ver_gdc_off
    u.hor_gdc_mode_req = hor_gdc_off
    u.lon_ctl_mode_req = lon_SAS_off
    u.lat_ctl_mode_req = lat_SAS_off

    #IMPORTANT: update avionics outputs
    f_disc!(ac.avionics, ac.physics, 1)

    return result

end

function Aircraft.linearize!(ac::System{<:Cessna172MCS}, args...; kwargs...)
    linearize!(ac.physics, args...; kwargs...)
end


# ############################ Joystick Mappings #################################

function IODevices.assign!(sys::System{<:Cessna172MCS}, joystick::Joystick,
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

# include(joinpath(@__DIR__, "design", "mcs_design.jl")); using .MCSDesign

end #module