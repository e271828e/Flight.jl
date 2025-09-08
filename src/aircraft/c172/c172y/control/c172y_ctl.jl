module C172YControl

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes
using EnumX

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
    Dummy, SameLine, NewLine, IsItemActive, IsItemActivated, Separator, Text,
    Checkbox, RadioButton, TableNextColumn, TableNextRow, BeginTable, EndTable

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

function is_on_gnd(mdl::Model{<:C172.Vehicle})
    any(SVector{3}(leg.strut.wow for leg in mdl.y.systems.ldg))
end


################################################################################
############################## ControlLawsLon ##################################

@enumx T=ModeControlLonEnum ModeControlLon begin
    direct = 0 #direct throttle & elevator
    sas = 1 #direct throttle + pitch SAS
    thr_q = 2 #direct throttle + pitch rate
    thr_θ = 3 #direct throttle + pitch angle
    thr_EAS = 4 #direct throttle + EAS
    EAS_q = 5 #EAS + pitch rate
    EAS_θ = 6 #EAS + pitch angle
    EAS_clm = 7 #EAS + climb rate
    EAS_alt = 8 #EAS + altitude hold
end

@enumx T=AltTrackingStateEnum AltTrackingState begin
    acquire = 0
    hold = 1
end

using .ModeControlLon: ModeControlLonEnum
using .AltTrackingState: AltTrackingStateEnum

function te2te_enabled(mode::ModeControlLonEnum)
    mode === ModeControlLon.sas || mode === ModeControlLon.thr_q ||
    mode === ModeControlLon.thr_θ || mode === ModeControlLon.EAS_q ||
    mode === ModeControlLon.EAS_θ || mode === ModeControlLon.EAS_clm
end

function q2e_enabled(mode::ModeControlLonEnum)
    mode === ModeControlLon.thr_q || mode === ModeControlLon.thr_θ ||
    mode === ModeControlLon.EAS_q || mode === ModeControlLon.EAS_θ ||
    mode === ModeControlLon.EAS_clm
end

function θ2q_enabled(mode::ModeControlLonEnum)
    mode === ModeControlLon.thr_θ || mode === ModeControlLon.EAS_θ ||
    mode === ModeControlLon.EAS_clm
end

function v2t_enabled(mode::ModeControlLonEnum)
    mode === ModeControlLon.EAS_q || mode === ModeControlLon.EAS_θ ||
    mode === ModeControlLon.EAS_clm
end

function c2θ_enabled(mode::ModeControlLonEnum)
    mode === ModeControlLon.EAS_clm
end

tv2te_enabled(mode::ModeControlLonEnum) = (mode === ModeControlLon.thr_EAS)

vh2te_enabled(mode::ModeControlLonEnum) = (mode === ModeControlLon.EAS_alt)


################################################################################
############################### LQR Vectors ####################################

#state vector for complete longitudinal dynamics
@kwdef struct XLonFull <: FieldVector{9, Float64}
    q::Float64 = 0.0 #pitch rate
    θ::Float64 = 0.0 #pitch angle
    EAS::Float64 = 0.0 #equivalent airspeed
    α::Float64 = 0.0 #AoA
    h::Float64 = 0.0 #altitude
    α_filt::Float64 = 0.0 #filtered AoA (from aerodynamics model)
    n_eng::Float64 = 0.0 #engine speed (rad/s)
    thr_p::Float64 = 0.0 #throttle actuator state
    ele_p::Float64 = 0.0 #elevator actuator state
end

#assemble state vector from vehicle
function XLonFull(vehicle::Model{<:C172Y.Vehicle})

    @unpack systems, airflow, kinematics = vehicle.y
    @unpack pwp, aero, act = systems
    @unpack e_nb, ω_eb_b, h_e = kinematics

    q = ω_eb_b[2]
    θ = e_nb.θ
    EAS = airflow.EAS
    α = aero.α
    h = h_e
    α_filt = aero.α_filt
    n_eng = pwp.engine.n
    thr_p = act.throttle.pos
    ele_p = act.elevator.pos

    XLonFull(; q, θ, EAS, α, h, α_filt, n_eng, thr_p, ele_p)

end

################################################################################

#state vector for reduced longitudinal dynamics
@kwdef struct XLonRed <: FieldVector{8, Float64}
    q::Float64 = 0.0 #pitch rate
    θ::Float64 = 0.0 #pitch angle
    EAS::Float64 = 0.0 #equivalent airspeed
    α::Float64 = 0.0 #AoA
    α_filt::Float64 = 0.0 #filtered AoA (from aerodynamics model)
    n_eng::Float64 = 0.0 #engine speed (rad/s)
    thr_p::Float64 = 0.0 #throttle actuator state
    ele_p::Float64 = 0.0 #elevator actuator state
end

#assemble state vector from vehicle
function XLonRed(vehicle::Model{<:C172Y.Vehicle})

    @unpack systems, airflow, kinematics = vehicle.y
    @unpack pwp, aero, act = systems
    @unpack e_nb, ω_eb_b, h_e = kinematics

    q = ω_eb_b[2]
    θ = e_nb.θ
    EAS = airflow.EAS
    α = aero.α
    α_filt = aero.α_filt
    n_eng = pwp.engine.n
    thr_p = act.throttle.pos
    ele_p = act.elevator.pos

    XLonRed(; q, θ, EAS, α, α_filt, n_eng, thr_p, ele_p)

end

################################################################################

#control vector for reduced and full longitudinal dynamics
@kwdef struct ULon{T} <: FieldVector{2, T}
    throttle_cmd::T = 0.0
    elevator_cmd::T = 0.0
end

################################################################################

#command vector for throttle + elevator SAS LQR tracker
@kwdef struct Zte <: FieldVector{2, Float64}
    throttle_cmd::Float64 = 0.0
    elevator_cmd::Float64 = 0.0
end

function Zte(vehicle::Model{<:C172Y.Vehicle})
    throttle_cmd = vehicle.y.systems.act.throttle.cmd
    elevator_cmd = vehicle.y.systems.act.elevator.cmd
    Zte(; throttle_cmd, elevator_cmd)
end

################################################################################

#command vector for EAS + altitude LQR tracker
@kwdef struct Zvh <: FieldVector{2, Float64}
    EAS::Float64 = 0.0
    h::Float64 = 0.0
end

function Zvh(vehicle::Model{<:C172Y.Vehicle})
    EAS = vehicle.y.airflow.EAS
    h = vehicle.y.kinematics.h_e
    Zvh(; EAS, h)
end

################################################################################

#command vector for throttle + EAS LQR tracker
@kwdef struct Ztv <: FieldVector{2, Float64}
    throttle_cmd::Float64 = 0.0
    EAS::Float64 = 0.0
end

function Ztv(vehicle::Model{<:C172Y.Vehicle})
    throttle_cmd = vehicle.y.systems.act.throttle.cmd
    EAS = vehicle.y.airflow.EAS
    Ztv(; throttle_cmd, EAS)
end

################################## Model ######################################

@kwdef struct ControlLawsLon{LQ1 <: LQRTrackerLookup, LQ2 <: LQRTrackerLookup, LP <: PIDLookup} <: ModelDefinition
    te2te_lookup::LQ1 = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "te2te_lookup.h5"))
    tv2te_lookup::LQ1 = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "tv2te_lookup.h5"))
    vh2te_lookup::LQ2 = load_lqr_tracker_lookup(joinpath(@__DIR__, "data", "vh2te_lookup.h5"))
    q2e_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "q2e_lookup.h5"))
    c2θ_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "c2θ_lookup.h5"))
    v2t_lookup::LP = load_pid_lookup(joinpath(@__DIR__, "data", "v2t_lookup.h5"))
    te2te_lqr::LQRTracker{8, 2, 2, 16, 4} = LQRTracker{8, 2, 2}()
    tv2te_lqr::LQRTracker{8, 2, 2, 16, 4} = LQRTracker{8, 2, 2}()
    vh2te_lqr::LQRTracker{9, 2, 2, 18, 4} = LQRTracker{9, 2, 2}()
    q2e_int::Integrator = Integrator()
    q2e_pid::PID = PID()
    c2θ_pid::PID = PID()
    v2t_pid::PID = PID()
    k_p_θ::Float64 = 1.0
    h_thr::Float64 = 10.0 #error threshold for altitude tracking mode switch
    h_hys::Float64 = 1.0 #hysteresis for altitude tracking mode switch (h_hys < h_thr)
end

@kwdef mutable struct ControlLawsLonS
    h_state::AltTrackingStateEnum = AltTrackingState.hold
end

@kwdef mutable struct ControlLawsLonU
    mode_req::ModeControlLonEnum = ModeControlLon.direct
    throttle_axis::Ranged{Float64, 0., 1.} = 0.0
    throttle_offset::Ranged{Float64, 0., 1.} = 0.0
    elevator_axis::Ranged{Float64, -1., 1.} = 0.0
    elevator_offset::Ranged{Float64, -1., 1.} = 0.0
    q_ref::Float64 = 0.0 #pitch rate reference
    θ_ref::Float64 = 0.0 #pitch angle reference
    EAS_ref::Float64 = C172.TrimParameters().EAS #equivalent airspeed reference
    clm_ref::Float64 = 0.0 #climb rate reference
    h_ref::HEllip = HEllip(0.0)
end

@kwdef struct ControlLawsLonY
    mode::ModeControlLonEnum = ModeControlLon.direct
    throttle_ref::Ranged{Float64, 0., 1.} = 0.0
    elevator_ref::Ranged{Float64, -1., 1.} = 0.0
    q_ref::Float64 = 0.0
    θ_ref::Float64 = 0.0
    EAS_ref::Float64 = C172.TrimParameters().EAS
    clm_ref::Float64 = 0.0
    h_ref::HEllip = HEllip(0.0) #avoid h_ref in output to keep output type isbits
    h_state::AltTrackingStateEnum = AltTrackingState.hold
    throttle_cmd::Ranged{Float64, 0., 1.} = 0.0
    elevator_cmd::Ranged{Float64, -1., 1.} = 0.0
    te2te_lqr::LQRTrackerOutput{8, 2, 2, 16, 4} = LQRTrackerOutput{8, 2, 2}()
    tv2te_lqr::LQRTrackerOutput{8, 2, 2, 16, 4} = LQRTrackerOutput{8, 2, 2}()
    vh2te_lqr::LQRTrackerOutput{9, 2, 2, 18, 4} = LQRTrackerOutput{9, 2, 2}()
    q2e_int::IntegratorOutput = IntegratorOutput()
    q2e_pid::PIDOutput = PIDOutput()
    c2θ_pid::PIDOutput = PIDOutput()
    v2t_pid::PIDOutput = PIDOutput()
end

Modeling.S(::ControlLawsLon) = ControlLawsLonS()
Modeling.U(::ControlLawsLon) = ControlLawsLonU()
Modeling.Y(::ControlLawsLon) = ControlLawsLonY()

function Modeling.init!(mdl::Model{<:ControlLawsLon})
    #when te2te_lqr is active, the actual throttle and elevator saturation
    #states are observed at its output. every non-LQR compensator built upon
    #te2te_lqr (q2e_int, q2e_pid, c2θ_pid) will receive these saturation states
    #as external saturation inputs to halt integration. therefore, we don't need
    #to set explicit output bounds for these compensators.
    foreach((mdl.te2te_lqr, mdl.tv2te_lqr, mdl.vh2te_lqr)) do lqr
        lqr.u.bound_lo .= ULon(; throttle_cmd = 0, elevator_cmd = -1)
        lqr.u.bound_hi .= ULon(; throttle_cmd = 1, elevator_cmd = 1)
    end
end


function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:ControlLawsLon},
                        vehicle::Model{<:C172Y.Vehicle})

    @unpack mode_req, throttle_axis, throttle_offset, elevator_axis,
            elevator_offset, q_ref, θ_ref, EAS_ref, clm_ref, h_ref = mdl.u
    @unpack te2te_lqr, tv2te_lqr, vh2te_lqr, q2e_int, q2e_pid,
            c2θ_pid, v2t_pid = mdl.submodels
    @unpack te2te_lookup, tv2te_lookup, vh2te_lookup, q2e_lookup,
            c2θ_lookup, v2t_lookup, k_p_θ, h_thr, h_hys = mdl.constants

    EAS = vehicle.y.airflow.EAS
    h_e = vehicle.y.kinematics.h_e
    _, q, r = vehicle.y.kinematics.ω_wb_b
    @unpack θ, φ = vehicle.y.kinematics.e_nb
    clm = -vehicle.y.kinematics.v_eb_n[3]
    h_err = h_ref - h_e
    h_state = mdl.s.h_state
    mode_prev = mdl.y.mode

    throttle_ref = Float64(throttle_axis + throttle_offset)
    elevator_ref = Float64(elevator_axis + elevator_offset)

    #if not overridden by te2te, vh2te or tv2te, actuation commands are simply
    #their respective reference values
    throttle_cmd = throttle_ref
    elevator_cmd = elevator_ref

    if is_on_gnd(vehicle)
        mode = ModeControlLon.direct
    else #air
        if mode_req === ModeControlLon.EAS_alt
            if h_state === AltTrackingState.acquire
                mode = ModeControlLon.thr_EAS
                throttle_ref = h_err > 0 ? 1.0 : 0.0 #full throttle to climb, idle to descend
                (abs(h_err) < h_thr - h_hys) && (mdl.s.h_state = AltTrackingState.hold)
            else #AltTrackingState.hold
                mode = ModeControlLon.EAS_alt
                (abs(h_err) > h_thr + h_hys) && (mdl.s.h_state = AltTrackingState.acquire)
            end
        else #mode_req != EAS_alt
            mode = mode_req
        end
    end

    if te2te_enabled(mode) #throttle_cmd and elevator_cmd overridden by te2te SAS

        u_lon_sat = ULon(te2te_lqr.y.out_sat)

        if v2t_enabled(mode) #throttle_ref overridden by v2t

            Control.Discrete.assign!(v2t_pid, v2t_lookup(EAS, Float64(h_e)))

            if mode != mode_prev
                Control.reset!(v2t_pid)
                k_i = v2t_pid.u.k_i
                (k_i != 0) && (v2t_pid.s.x_i0 = Float64(mdl.y.throttle_cmd))
            end

            v2t_pid.u.input = EAS_ref - EAS
            v2t_pid.u.sat_ext = u_lon_sat.throttle_cmd
            f_periodic!(v2t_pid)
            throttle_ref = v2t_pid.y.output

        end

        if q2e_enabled(mode) #elevator_ref overridden by q2e

            Control.Discrete.assign!(q2e_pid, q2e_lookup(EAS, Float64(h_e)))

            if mode != mode_prev

                Control.reset!(q2e_int)
                Control.reset!(q2e_pid)
                k_i = q2e_pid.u.k_i
                (k_i != 0) && (q2e_pid.s.x_i0 = Zte(te2te_lqr.u.z_ref).elevator_cmd)

            end

            if θ2q_enabled(mode) #q_ref overridden by θ2q

                if c2θ_enabled(mode) #θ_ref overridden by c2θ

                    Control.Discrete.assign!(c2θ_pid, c2θ_lookup(EAS, Float64(h_e)))

                    if mode != mode_prev
                        Control.reset!(c2θ_pid)
                        k_i = c2θ_pid.u.k_i
                        (k_i != 0) && (c2θ_pid.s.x_i0 = θ)
                    end

                    c2θ_pid.u.input = clm_ref - clm
                    c2θ_pid.u.sat_ext = u_lon_sat.elevator_cmd
                    f_periodic!(c2θ_pid)
                    θ_ref = c2θ_pid.y.output

                end

                θ_dot_ref = k_p_θ * (θ_ref - θ)
                φ_bnd = clamp(φ, -π/3, π/3)
                q_ref = 1/cos(φ_bnd) * θ_dot_ref + r * tan(φ_bnd)

            end

            q2e_int.u.input = q_ref - q
            q2e_int.u.sat_ext = u_lon_sat.elevator_cmd
            f_periodic!(q2e_int)

            q2e_pid.u.input = q2e_int.y.output
            q2e_pid.u.sat_ext = u_lon_sat.elevator_cmd
            f_periodic!(q2e_pid)
            elevator_ref = q2e_pid.y.output

        end

        Control.Discrete.assign!(te2te_lqr, te2te_lookup(EAS, Float64(h_e)))

        #te2te is purely proportional, so it doesn't need resetting
        te2te_lqr.u.x .= XLonRed(vehicle) #state feedback
        te2te_lqr.u.z .= Zte(vehicle) #command variable feedback
        te2te_lqr.u.z_ref .= Zte(; throttle_cmd = throttle_ref,
            elevator_cmd = elevator_ref) #command variable reference
        f_periodic!(te2te_lqr)
        @unpack throttle_cmd, elevator_cmd = ULon(te2te_lqr.y.output)

    end

    if tv2te_enabled(mode) #throttle_cmd and elevator_cmd overridden by tv2te

        Control.Discrete.assign!(tv2te_lqr, tv2te_lookup(EAS, Float64(h_e)))

        (mode != mode_prev) && Control.reset!(tv2te_lqr)

        tv2te_lqr.u.x .= XLonRed(vehicle)
        tv2te_lqr.u.z .= Ztv(vehicle)
        tv2te_lqr.u.z_ref .= Ztv(; throttle_cmd = throttle_ref, EAS = EAS_ref)
        f_periodic!(tv2te_lqr)
        @unpack throttle_cmd, elevator_cmd = ULon(tv2te_lqr.y.output)

    end

    if vh2te_enabled(mode) #throttle_cmd and elevator_cmd overridden by vh2te

        Control.Discrete.assign!(vh2te_lqr, vh2te_lookup(EAS, Float64(h_e)))

        (mode != mode_prev) && Control.reset!(vh2te_lqr)

        vh2te_lqr.u.x .= XLonFull(vehicle)
        vh2te_lqr.u.z .= Zvh(vehicle)
        vh2te_lqr.u.z_ref .= Zvh(; EAS = EAS_ref, h = h_ref)
        f_periodic!(vh2te_lqr)
        @unpack throttle_cmd, elevator_cmd = ULon(vh2te_lqr.y.output)

    end

    mdl.y = ControlLawsLonY(; mode, throttle_ref, elevator_ref, q_ref, θ_ref,
        EAS_ref, clm_ref, h_ref, h_state, throttle_cmd, elevator_cmd,
        te2te_lqr = te2te_lqr.y, tv2te_lqr = tv2te_lqr.y, vh2te_lqr = vh2te_lqr.y,
        q2e_int = q2e_int.y, q2e_pid = q2e_pid.y,
        c2θ_pid = c2θ_pid.y, v2t_pid = v2t_pid.y)

end


function AircraftBase.assign!(systems::Model{<:C172Y.Systems},
                          mdl::Model{<:ControlLawsLon})

    @unpack act = systems.submodels
    @unpack throttle_cmd, elevator_cmd = mdl.y

    act.throttle.u[] = throttle_cmd
    act.elevator.u[] = elevator_cmd

end


############################# Initialization ###################################

function Modeling.init!(lon::Model{<:ControlLawsLon},
                        vehicle::Model{<:C172Y.Vehicle})

    #we assume that the vehicle's y has already been updated to its trim
    #value by init!(vehicle, params)
    @unpack u = lon
    @unpack throttle, elevator = vehicle.y.systems.act
    @unpack ω_wb_b, v_eb_n, e_nb, h_e = vehicle.y.kinematics
    @unpack EAS = vehicle.y.airflow

    #reset all controller submodels
    Control.reset!(lon)

    #make inputs consistent with the vehicle status, so
    u.throttle_axis = throttle.pos
    u.elevator_axis = elevator.pos
    u.throttle_offset = 0
    u.elevator_offset = 0
    u.q_ref = ω_wb_b[2]
    u.θ_ref = e_nb.θ
    u.EAS_ref = EAS
    u.clm_ref = -v_eb_n[3]
    u.h_ref = h_e

    #for the trim condition to be preserved when the simulation is started with
    #SAS-based modes enabled, we need a post-trim update of the control laws
    #with the inner LQR SAS loops enabled. this loads the z_trim, u_trim and
    #x_trim values corresponding to the trim state into the LQR trackers, and
    #runs them once.

    #after this, their outputs will match the actuator commands required by the
    #trim condition. this match is only approximate, because in general, the
    #trim values loaded from the lookup tables will be interpolated, rather than
    #exactly computed at specific controller design points. however, it is close
    #enough.

    #without this step, on the first call to f_ode!(::Model{<:Aircraft},...),
    #the outputs from the uninitialized SAS would be assigned to the vehicle,
    #overwriting the actuator commands set by the trim function with the
    #incorrect default values

    #initialize te2te outputs
    u.mode_req = ModeControlLon.sas
    f_periodic!(lon, vehicle)

    #initialize tv2te outputs
    u.mode_req = ModeControlLon.thr_EAS
    f_periodic!(lon, vehicle)

    #initialize vh2te outputs
    u.mode_req = ModeControlLon.EAS_alt
    f_periodic!(lon, vehicle)

    #restore direct mode
    u.mode_req = ModeControlLon.direct
    f_periodic!(lon, vehicle)

end


################################# JSON3 ########################################

StructTypes.StructType(::Type{ControlLawsLonU}) = StructTypes.Mutable()

#replace Unicode characters from struct fields in the JSON string
StructTypes.names(::Type{ControlLawsLonU}) = ((:θ_ref, :theta_ref),)

#enable JSON parsing of integers as ModeControlLonEnum
StructTypes.StructType(::Type{ModeControlLonEnum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{ModeControlLonEnum}) = Int32 #default enum type
StructTypes.lower(x::ModeControlLonEnum) = Int32(x)


################################### GUI ########################################

function GUI.draw!(lon::Model{<:ControlLawsLon},
                    vehicle::Model{<:C172Y.Vehicle})

    @unpack u, y, Δt = lon

    @unpack systems, kinematics, dynamics, airflow = vehicle.y
    @unpack e_nb, ω_wb_b, v_eb_n, h_e, h_o, n_e = kinematics
    @unpack EAS = airflow
    @unpack θ = e_nb

    q = ω_wb_b[2]
    clm = -v_eb_n[3]

    if CImGui.CollapsingHeader("Longitudinal Control")

        if BeginTable("LonCtlModes", 3, CImGui.ImGuiTableFlags_SizingStretchProp )#| CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
            TableNextRow()
                TableNextColumn();
                    Text("Mode")
                TableNextColumn();
                    mode_button("Direct##Lon", ModeControlLon.direct, lon.u.mode_req, y.mode)
                    IsItemActive() && (lon.u.mode_req = ModeControlLon.direct); SameLine()
                    mode_button("Throttle/Pitch SAS", ModeControlLon.sas, lon.u.mode_req, y.mode)
                    IsItemActive() && (lon.u.mode_req = ModeControlLon.sas); SameLine()
                    mode_button("Throttle + Pitch Rate", ModeControlLon.thr_q, lon.u.mode_req, y.mode)
                    IsItemActive() && (lon.u.mode_req = ModeControlLon.thr_q; lon.u.q_ref = 0); SameLine()
                    mode_button("Throttle + Pitch Angle", ModeControlLon.thr_θ, lon.u.mode_req, y.mode)
                    IsItemActive() && (lon.u.mode_req = ModeControlLon.thr_θ; lon.u.θ_ref = θ); SameLine()
                    mode_button("Throttle + EAS", ModeControlLon.thr_EAS, lon.u.mode_req, y.mode)
                    IsItemActive() && (lon.u.mode_req = ModeControlLon.thr_EAS; lon.u.EAS_ref = EAS)
                    mode_button("EAS + Pitch Rate", ModeControlLon.EAS_q, lon.u.mode_req, y.mode)
                    IsItemActive() && (lon.u.mode_req = ModeControlLon.EAS_q; lon.u.q_ref = 0; lon.u.EAS_ref = EAS); SameLine()
                    mode_button("EAS + Pitch Angle", ModeControlLon.EAS_θ, lon.u.mode_req, y.mode)
                    IsItemActive() && (lon.u.mode_req = ModeControlLon.EAS_θ; lon.u.EAS_ref = EAS; lon.u.θ_ref = θ); SameLine()
                    mode_button("EAS + Climb Rate", ModeControlLon.EAS_clm, lon.u.mode_req, y.mode)
                    IsItemActive() && (lon.u.mode_req = ModeControlLon.EAS_clm; lon.u.EAS_ref = EAS; lon.u.clm_ref = clm); SameLine()
                    mode_button("EAS + Altitude Hold", ModeControlLon.EAS_alt, lon.u.mode_req, y.mode)
                    IsItemActive() && (lon.u.mode_req = ModeControlLon.EAS_alt; lon.u.EAS_ref = EAS; lon.u.h_ref = h_e); SameLine()
                TableNextColumn();
            TableNextRow()
                TableNextColumn();
                    AlignTextToFramePadding(); Text("Throttle Axis")
                    AlignTextToFramePadding(); Text("Throttle Offset")
                TableNextColumn();
                    PushItemWidth(-10)
                    lon.u.throttle_axis = safe_slider("Throttle Axis", lon.u.throttle_axis, "%.6f")
                    lon.u.throttle_offset = safe_input("Throttle_Offset", lon.u.throttle_offset, 0.01, 0.1, "%.3f")
                    PopItemWidth()
            TableNextRow()
                TableNextColumn();
                    AlignTextToFramePadding(); Text("Elevator Axis")
                    AlignTextToFramePadding(); Text("Elevator Offset")
                TableNextColumn();
                    PushItemWidth(-10)
                    lon.u.elevator_axis = safe_slider("Elevator Axis", lon.u.elevator_axis, "%.6f")
                    lon.u.elevator_offset = safe_input("Elevator Offset", lon.u.elevator_offset, 0.01, 0.1, "%.3f")
                    PopItemWidth()
            TableNextRow()
                TableNextColumn(); AlignTextToFramePadding(); Text("Pitch Rate (deg/s)")
                TableNextColumn();
                    PushItemWidth(-10)
                    lon.u.q_ref = safe_slider("Pitch Rate", rad2deg(lon.u.q_ref), -10, 10, "%.3f") |> deg2rad
                    PopItemWidth()
                TableNextColumn(); Text(@sprintf("%.3f", rad2deg(q)))
            TableNextRow()
                TableNextColumn(); AlignTextToFramePadding(); Text("Pitch Angle (deg)")
                TableNextColumn();
                    PushItemWidth(-10)
                    lon.u.θ_ref = safe_slider("Pitch Angle", rad2deg(lon.u.θ_ref), -15, 15, "%.3f") |> deg2rad
                    PopItemWidth()
                TableNextColumn(); Text(@sprintf("%.3f", rad2deg(θ)))
            TableNextRow()
                TableNextColumn(); AlignTextToFramePadding(); Text("EAS (m/s)")
                TableNextColumn();
                    PushItemWidth(-10)
                    lon.u.EAS_ref = safe_input("EAS", lon.u.EAS_ref, 0.1, 1.0, "%.3f")
                    PopItemWidth()
                TableNextColumn(); Text(@sprintf("%.3f", EAS))
            TableNextRow()
                TableNextColumn(); AlignTextToFramePadding(); Text("Climb Rate (m/s)")
                TableNextColumn();
                    PushItemWidth(-10)
                    lon.u.clm_ref = safe_input("Climb Rate", lon.u.clm_ref, 0.1, 1.0, "%.3f")
                    PopItemWidth()
                TableNextColumn(); Text(@sprintf("%.3f", clm))
            TableNextRow()
                TableNextColumn(); AlignTextToFramePadding(); Text("Altitude (m)")
                TableNextColumn();
                @cstatic h_datum = :ellip begin
                    RadioButton("Ellipsoidal", h_datum === :ellip); SameLine()
                    IsItemActive() && (h_datum = :ellip)
                    RadioButton("Orthometric", h_datum === :orth)
                    IsItemActive() && (h_datum = :orth)
                    SameLine()
                    PushItemWidth(-10)
                    h_ref_ellip = lon.u.h_ref
                    h_ref_orth = HOrth(h_ref_ellip, n_e)
                    h_ref_f64 = (h_datum === :ellip ? Float64(h_ref_ellip) : Float64(h_ref_orth))
                    h_ref_f64 = safe_input("Altitude Setpoint", h_ref_f64, 1, 1.0, "%.3f")
                    lon.u.h_ref = (h_datum === :ellip ? HEllip(h_ref_f64) : HEllip(HOrth(h_ref_f64), n_e))
                    PopItemWidth()
                TableNextColumn();
                    (h_datum === :ellip) && Text(@sprintf("%.3f", Float64(h_e)))
                    (h_datum === :orth) && Text(@sprintf("%.3f", Float64(h_o)))
                end
            TableNextRow()
                TableNextColumn();
                TableNextColumn();
            EndTable()
        end

        Separator()

        if BeginTable("Actuator Data", 5, CImGui.ImGuiTableFlags_SizingStretchSame)# | CImGui.ImGuiTableFlags_BordersInner)
            TableNextRow()
                TableNextColumn();
                TableNextColumn(); Text("Reference")
                TableNextColumn(); Text("Command")
                TableNextColumn(); Text("Position")
            TableNextRow()
                TableNextColumn(); Text("Throttle")
                TableNextColumn(); Text(@sprintf("%.6f", Float64(y.throttle_ref)))
                TableNextColumn(); Text(@sprintf("%.6f", Float64(y.throttle_cmd)))
                TableNextColumn(); Text(@sprintf("%.6f", Float64(systems.act.throttle.pos)))
            TableNextRow()
                TableNextColumn(); Text("Elevator")
                TableNextColumn(); Text(@sprintf("%.6f", Float64(y.elevator_ref)))
                TableNextColumn(); Text(@sprintf("%.6f", Float64(y.elevator_cmd)))
                TableNextColumn(); Text(@sprintf("%.6f", Float64(systems.act.elevator.pos)))
            EndTable()
        end

        Separator()

    end

end

function draw_internals(lon::Model{<:ControlLawsLon}, p_open::Ref{Bool} = Ref(true))
    Begin("Lateral Control", p_open)
        Text("Sampling Period: $(lon.Δt)")
        Text("Mode: $(lon.y.mode)")
        foreach(keys(lon.submodels), values(lon.submodels)) do label, ss
            if CImGui.TreeNode(string(label))
                GUI.draw(ss)
                CImGui.TreePop()
            end
        end
    End()
end



################################################################################
############################# ControlLawsLat ###################################

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
@kwdef struct XLatRed <: FieldVector{8, Float64}
    p::Float64 = 0.0 #roll rate
    r::Float64 = 0.0 #yaw rate
    φ::Float64 = 0.0; #bank angle
    EAS::Float64 = 0.0 #equivalent airspeed
    β::Float64 = 0.0 #AoS
    β_filt::Float64 = 0.0; #filtered AoS (from aerodynamics model)
    ail_p::Float64 = 0.0; #aileron actuator states
    rud_p::Float64 = 0.0; #rudder actuator states
end

function XLatRed(vehicle::Model{<:C172Y.Vehicle})

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

    XLatRed(; p, r, φ, EAS, β, β_filt, ail_p, rud_p)

end

@kwdef struct ULatRed{T} <: FieldVector{2, T}
    aileron_cmd::T = 0.0
    rudder_cmd::T = 0.0
end

@kwdef struct Zφβ <: FieldVector{2, Float64}
    φ::Float64 = 0.0
    β::Float64 = 0.0
end

function Zφβ(vehicle::Model{<:C172Y.Vehicle})
    φ = vehicle.y.kinematics.e_nb.φ
    β = vehicle.y.systems.aero.β
    Zφβ(; φ, β)
end

@kwdef struct Zar <: FieldVector{2, Float64}
    aileron_cmd::Float64 = 0.0
    rudder_cmd::Float64 = 0.0
end

function Zar(vehicle::Model{<:C172Y.Vehicle})
    aileron_cmd = vehicle.y.systems.act.aileron.cmd
    rudder_cmd = vehicle.y.systems.act.rudder.cmd
    Zar(; aileron_cmd, rudder_cmd)
end

################################## Model ######################################

@kwdef struct ControlLawsLat{LQ <: LQRTrackerLookup, LP <: PIDLookup} <: ModelDefinition
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

@kwdef mutable struct ControlLawsLatU
    mode_req::ModeControlLatEnum = ModeControlLat.direct #lateral control mode
    aileron_axis::Ranged{Float64, -1., 1.} = 0.0
    aileron_offset::Ranged{Float64, -1., 1.} = 0.0
    rudder_axis::Ranged{Float64, -1., 1.} = 0.0
    rudder_offset::Ranged{Float64, -1., 1.} = 0.0
    p_ref::Float64 = 0.0 #roll rate reference
    β_ref::Float64 = 0.0 #sideslip angle reference
    φ_ref::Float64 = 0.0 #bank angle reference
    χ_ref::Float64 = 0.0 #course angle reference
end

@kwdef struct ControlLawsLatY
    mode::ModeControlLatEnum = ModeControlLat.direct
    aileron_ref::Ranged{Float64, -1., 1.} = 0.0
    rudder_ref::Ranged{Float64, -1., 1.} = 0.0
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

Modeling.U(::ControlLawsLat) = ControlLawsLatU()
Modeling.Y(::ControlLawsLat) = ControlLawsLatY()

function Modeling.init!(mdl::Model{<:ControlLawsLat})

    foreach((mdl.φβ2ar_lqr, mdl.ar2ar_lqr)) do lqr
        lqr.u.bound_lo .= ULatRed(; aileron_cmd = -1, rudder_cmd = -1)
        lqr.u.bound_hi .= ULatRed(; aileron_cmd = 1, rudder_cmd = 1)
    end

    #set φ reference limits for the course angle compensator output
    mdl.χ2φ_pid.u.bound_lo = -π/4
    mdl.χ2φ_pid.u.bound_hi = π/4

end

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:ControlLawsLat},
                        vehicle::Model{<:C172Y.Vehicle})

    @unpack mode_req, aileron_axis, aileron_offset, rudder_axis, rudder_offset,
            p_ref, β_ref, φ_ref, χ_ref = mdl.u
    @unpack ar2ar_lqr, φβ2ar_lqr, p2φ_int, p2φ_pid, χ2φ_pid = mdl.submodels
    @unpack ar2ar_lookup, φβ2ar_lookup, p2φ_lookup, χ2φ_lookup = mdl.constants
    @unpack airflow, kinematics, systems = vehicle.y

    EAS = airflow.EAS
    h_e = Float64(kinematics.h_e)

    @unpack θ, φ = kinematics.e_nb
    p, _, _ = kinematics.ω_wb_b
    β = systems.aero.β
    mode_prev = mdl.y.mode

    mode = (is_on_gnd(vehicle) ? ModeControlLat.direct : mode_req)

    aileron_ref = Float64(aileron_axis + aileron_offset)
    rudder_ref = Float64(rudder_axis + rudder_offset)

    aileron_cmd = aileron_ref
    rudder_cmd = rudder_ref

    if ar2ar_enabled(mode) #aileron_cmd and #rudder_cmd overridden by ar2ar

        #no integral control, so no need for reset on mode change
        Control.Discrete.assign!(ar2ar_lqr, ar2ar_lookup(EAS, Float64(h_e)))

        ar2ar_lqr.u.x .= XLatRed(vehicle)
        ar2ar_lqr.u.z .= Zar(vehicle)
        ar2ar_lqr.u.z_ref .= Zar(; aileron_cmd = aileron_ref, rudder_cmd = rudder_ref)
        f_periodic!(ar2ar_lqr)
        @unpack aileron_cmd, rudder_cmd = ULatRed(ar2ar_lqr.y.output)

    end

    if φβ2ar_enabled(mode) #aileron_cmd and #rudder_cmd overridden by φβ2ar

        u_lat_sat = ULatRed(φβ2ar_lqr.y.out_sat)

        if p2φ_enabled(mode) #φ_ref overridden by roll rate tracker

            Control.Discrete.assign!(p2φ_pid, p2φ_lookup(EAS, Float64(h_e)))

            if mode != mode_prev
                #our next φ output must match φ reference at φβ2ar input
                Control.reset!(p2φ_int)
                Control.reset!(p2φ_pid)
                k_i = p2φ_pid.u.k_i
                (k_i != 0) && (p2φ_pid.s.x_i0 = Zφβ(φβ2ar_lqr.u.z_ref).φ)
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
                (k_i != 0) && (χ2φ_pid.s.x_i0 = Zφβ(φβ2ar_lqr.u.z_ref).φ)
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

        φβ2ar_lqr.u.x .= XLatRed(vehicle)
        φβ2ar_lqr.u.z .= Zφβ(vehicle)
        φβ2ar_lqr.u.z_ref .= Zφβ(; φ = φ_ref, β = β_ref)
        f_periodic!(φβ2ar_lqr)
        @unpack aileron_cmd, rudder_cmd = ULatRed(φβ2ar_lqr.y.output)

    end

    mdl.y = ControlLawsLatY(; mode, aileron_ref, rudder_ref,
        p_ref, β_ref, φ_ref, χ_ref, aileron_cmd, rudder_cmd,
        ar2ar_lqr = ar2ar_lqr.y, φβ2ar_lqr = φβ2ar_lqr.y,
        p2φ_int = p2φ_int.y, p2φ_pid = p2φ_pid.y, χ2φ_pid = χ2φ_pid.y)

end


function AircraftBase.assign!(systems::Model{<:C172Y.Systems},
                          mdl::Model{<:ControlLawsLat})

    @unpack act = systems.submodels
    @unpack aileron_cmd, rudder_cmd = mdl.y

    act.aileron.u[] = aileron_cmd
    act.rudder.u[] = rudder_cmd

end


function Modeling.init!(lat::Model{<:ControlLawsLat},
                        vehicle::Model{<:C172Y.Vehicle})

    #we assume that the vehicle's y has already been updated to its trim value
    #by init!(vehicle, params)
    @unpack u = lat
    @unpack ω_wb_b, v_eb_n, e_nb, χ_gnd, ϕ_λ, h_e = vehicle.y.kinematics
    @unpack aileron, rudder = vehicle.y.systems.act
    @unpack EAS = vehicle.y.airflow
    @unpack β = vehicle.y.systems.aero

    #reset all controller submodels
    Control.reset!(lat)

    #make ControlLaws inputs consistent with vehicle status
    u.aileron_axis = aileron.pos
    u.rudder_axis = rudder.pos
    u.aileron_offset = 0
    u.rudder_offset = 0
    u.p_ref = ω_wb_b[1]
    u.φ_ref = e_nb.φ
    u.β_ref = β
    u.χ_ref = χ_gnd

    #initialize ar2ar outputs
    u.mode_req = ModeControlLat.sas
    f_periodic!(lat, vehicle)

    #initialize φβ2ar outputs
    u.mode_req = ModeControlLat.ModeControlLat.φ_β
    f_periodic!(lat, vehicle)

    #restore direct mode
    u.mode_req = ModeControlLat.direct
    f_periodic!(lat, vehicle)

end


################################# JSON3 ########################################

StructTypes.StructType(::Type{ControlLawsLatU}) = StructTypes.Mutable()

#replace Unicode characters from struct fields in the JSON string
StructTypes.names(::Type{ControlLawsLatU}) = ((:χ_ref, :chi_ref),
    (:φ_ref, :phi_ref), (:β_ref, :beta_ref))

#enable JSON parsing of integers as ModeControlLatEnum
StructTypes.StructType(::Type{ModeControlLatEnum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{ModeControlLatEnum}) = Int32 #default enum type
StructTypes.lower(x::ModeControlLatEnum) = Int32(x)

#now we can do:
# JSON3.read(JSON3.write(ControlLawsLatU()), ControlLawsLatU)
# JSON3.read!(JSON3.write(ControlLawsLatU()), ControlLawsLatU())

################################### GUI ########################################

function GUI.draw!(lat::Model{<:ControlLawsLat},
                    vehicle::Model{<:C172Y.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna172Y Control Laws")

    @unpack u, y, Δt = lat

    @unpack systems, kinematics = vehicle.y
    @unpack e_nb, ω_wb_b, χ_gnd, v_eb_n = kinematics
    @unpack β = systems.aero
    @unpack ψ, φ = e_nb

    p = ω_wb_b[1]

    Begin(label, p_open)

    if CImGui.CollapsingHeader("Lateral Control")

        if BeginTable("LatCtlModes", 3, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
            TableNextRow()
                TableNextColumn();
                    Text("Mode")
                TableNextColumn();
                    mode_button("Direct##Lat", ModeControlLat.direct, lat.u.mode_req, y.mode); SameLine()
                    IsItemActive() && (lat.u.mode_req = ModeControlLat.direct)
                    mode_button("Roll/Yaw SAS", ModeControlLat.sas, lat.u.mode_req, y.mode); SameLine()
                    IsItemActive() && (lat.u.mode_req = ModeControlLat.sas)
                    mode_button("Roll Rate + AoS", ModeControlLat.p_β, lat.u.mode_req, y.mode); SameLine()
                    IsItemActive() && (lat.u.mode_req = ModeControlLat.p_β; lat.u.p_ref = 0; lat.u.β_ref = β)
                    mode_button("Bank Angle + AoS", ModeControlLat.ModeControlLat.φ_β, lat.u.mode_req, y.mode); SameLine()
                    IsItemActive() && (lat.u.mode_req = ModeControlLat.ModeControlLat.φ_β; lat.u.φ_ref = φ; lat.u.β_ref = β)
                    mode_button("Course Angle + AoS", ModeControlLat.χ_β, lat.u.mode_req, y.mode); SameLine()
                    IsItemActive() && (lat.u.mode_req = ModeControlLat.χ_β; lat.u.χ_ref = χ_gnd; lat.u.β_ref = β)
                TableNextColumn();
            TableNextRow()
                TableNextColumn();
                    AlignTextToFramePadding(); Text("Aileron Axis")
                    AlignTextToFramePadding(); Text("Aileron Offset")
                TableNextColumn();
                    PushItemWidth(-10)
                    lat.u.aileron_axis = safe_slider("Aileron Axis", lat.u.aileron_axis, "%.6f")
                    lat.u.aileron_offset = safe_input("Aileron Offset", lat.u.aileron_offset, 0.01, 0.1, "%.3f")
                    PopItemWidth()
            TableNextRow()
                TableNextColumn();
                    AlignTextToFramePadding(); Text("Rudder Axis")
                    AlignTextToFramePadding(); Text("Rudder Offset")
                TableNextColumn();
                    PushItemWidth(-10)
                    lat.u.rudder_axis = safe_slider("Rudder Axis", lat.u.rudder_axis, "%.6f")
                    lat.u.rudder_offset = safe_input("Rudder Offset", lat.u.rudder_offset, 0.01, 0.1, "%.3f")
                    PopItemWidth()
            TableNextRow()
                TableNextColumn(); AlignTextToFramePadding(); Text("Roll Rate (deg/s)")
                TableNextColumn();
                    PushItemWidth(-10)
                    lat.u.p_ref = safe_slider("Roll Rate", rad2deg(lat.u.p_ref), -30, 30, "%.3f") |> deg2rad
                    PopItemWidth()
                TableNextColumn(); Text(@sprintf("%.3f", rad2deg(p)))
            TableNextRow()
                TableNextColumn(); AlignTextToFramePadding(); Text("Bank Angle (deg)")
                TableNextColumn();
                    PushItemWidth(-10)
                    lat.u.φ_ref = safe_slider("Bank Angle", rad2deg(lat.u.φ_ref), -60, 60, "%.3f") |> deg2rad
                    PopItemWidth()
                TableNextColumn(); Text(@sprintf("%.3f", rad2deg(φ)))
            TableNextRow()
                TableNextColumn(); AlignTextToFramePadding(); Text("Course Angle (deg)")
                TableNextColumn();
                    PushItemWidth(-10)
                    lat.u.χ_ref = safe_slider("Course Angle", rad2deg(lat.u.χ_ref), -180, 180, "%.3f") |> deg2rad
                    PopItemWidth()
                TableNextColumn(); Text(@sprintf("%.3f", rad2deg(χ_gnd)))
            TableNextRow()
                TableNextColumn(); AlignTextToFramePadding(); Text("Sideslip Angle (deg)")
                TableNextColumn();
                    PushItemWidth(-10)
                    lat.u.β_ref = safe_slider("Sideslip Angle", rad2deg(lat.u.β_ref), -10, 10, "%.3f") |> deg2rad
                    PopItemWidth()
                TableNextColumn(); Text(@sprintf("%.3f", rad2deg(β)))
            EndTable()
        end

        Separator()

        if BeginTable("Actuator Data", 5, CImGui.ImGuiTableFlags_SizingStretchSame)# | CImGui.ImGuiTableFlags_BordersInner)
            TableNextRow()
                TableNextColumn();
                TableNextColumn(); Text("Reference")
                TableNextColumn(); Text("Command")
                TableNextColumn(); Text("Position")
            TableNextRow()
                TableNextColumn(); Text("Aileron")
                TableNextColumn(); Text(@sprintf("%.6f", Float64(y.aileron_ref)))
                TableNextColumn(); Text(@sprintf("%.6f", Float64(y.aileron_cmd)))
                TableNextColumn(); Text(@sprintf("%.6f", Float64(systems.act.aileron.pos)))
            TableNextRow()
                TableNextColumn(); Text("Rudder")
                TableNextColumn(); Text(@sprintf("%.6f", Float64(y.rudder_cmd)))
                TableNextColumn(); Text(@sprintf("%.6f", Float64(y.rudder_ref)))
                TableNextColumn(); Text(@sprintf("%.6f", Float64(systems.act.rudder.pos)))
            EndTable()
        end

        Separator()

    end

    End()

end

function draw_internals(lat::Model{<:ControlLawsLat}, p_open::Ref{Bool} = Ref(true))
    Begin("Lateral Control", p_open)
        Text("Sampling Period: $(lat.Δt)")
        Text("Mode: $(lat.y.mode)")
        foreach(keys(lat.submodels), values(lat.submodels)) do label, ss
            if CImGui.TreeNode(string(label))
                GUI.draw(ss)
                CImGui.TreePop()
            end
        end
    End()
end


################################################################################
################################# ControlLaws ##################################

@kwdef struct ControlLaws{C1 <: ControlLawsLon, C2 <: ControlLawsLat} <: AbstractAvionics
    lon::C1 = ControlLawsLon()
    lat::C2 = ControlLawsLat()
end

#define fields with submodel's names to create parent-child input linkage
@kwdef struct ControlLawsU
    lon::ControlLawsLonU = ControlLawsLonU()
    lat::ControlLawsLatU = ControlLawsLatU()
end

@kwdef struct ControlLawsY
    lon::ControlLawsLonY = ControlLawsLonY()
    lat::ControlLawsLatY = ControlLawsLatY()
end

Modeling.U(::ControlLaws) = ControlLawsU()
Modeling.Y(::ControlLaws) = ControlLawsY()

function Modeling.f_output!(mdl::Model{<:ControlLaws})
    mdl.y = ControlLawsY(mdl.lon.y, mdl.lat.y)
end

@no_ode ControlLaws
@no_step ControlLaws
@sm_periodic ControlLaws


########################### Update Methods #####################################

function AircraftBase.assign!(systems::Model{<:C172Y.Systems},
                                avionics::Model{<:ControlLaws})
    AircraftBase.assign!(systems, avionics.lon)
    AircraftBase.assign!(systems, avionics.lat)
end

############################# Initialization ###################################

function Modeling.init!(avionics::Model{<:ControlLaws},
                        vehicle::Model{<:C172Y.Vehicle})
    Modeling.init!(avionics.lon, vehicle)
    Modeling.init!(avionics.lat, vehicle)
end

################################## GUI #########################################

function GUI.draw!(ctl::Model{<:ControlLaws},
                    vehicle::Model{<:C172Y.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna172Y Control Laws")

    @unpack lon, lat = ctl.submodels

    Begin(label, p_open)

    GUI.draw!(lon, vehicle)
    GUI.draw!(lat, vehicle)
    GUI.draw(vehicle.y)

    if CImGui.CollapsingHeader("Internals")
        @cstatic c_lon=false c_lat=false c_alt=false begin
            @c Checkbox("Longitudinal Control##Internals", &c_lon)
            SameLine()
            @c Checkbox("Lateral Control##Internals", &c_lat)
            c_lon && @c draw_internals(lon, &c_lon)
            c_lat && @c draw_internals(lon, &c_lat)
        end
    end

    End()

end

end #module