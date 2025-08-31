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
    mode_req::ModeControlLonEnum = ModeControlLon.direct
    throttle_axis::Ranged{Float64, 0., 1.} = 0.0
    throttle_offset::Ranged{Float64, 0., 1.} = 0.0
    elevator_axis::Ranged{Float64, -1., 1.} = 0.0
    elevator_offset::Ranged{Float64, -1., 1.} = 0.0
    q_ref::Float64 = 0.0 #pitch rate reference
    θ_ref::Float64 = 0.0 #pitch angle reference
    EAS_ref::Float64 = C172.TrimParameters().EAS #equivalent airspeed reference
    clm_ref::Float64 = 0.0 #climb rate reference
end

@kwdef struct ControllerLonY
    mode::ModeControlLonEnum = ModeControlLon.direct
    throttle_ref::Ranged{Float64, 0., 1.} = 0.0
    elevator_ref::Ranged{Float64, -1., 1.} = 0.0
    q_ref::Float64 = 0.0
    θ_ref::Float64 = 0.0
    EAS_ref::Float64 = C172.TrimParameters().EAS
    clm_ref::Float64 = 0.0
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

    @unpack mode_req, throttle_axis, throttle_offset, elevator_axis,
            elevator_offset, q_ref, θ_ref, EAS_ref, clm_ref = mdl.u
    @unpack e2e_lqr, q2e_int, q2e_pid, v2θ_pid, c2θ_pid, v2t_pid = mdl.submodels
    @unpack e2e_lookup, q2e_lookup, v2θ_lookup, c2θ_lookup, v2t_lookup = mdl.constants

    EAS = vehicle.y.airflow.EAS
    h_e = Float64(vehicle.y.kinematics.h_e)
    _, q, r = vehicle.y.kinematics.ω_wb_b
    @unpack θ, φ = vehicle.y.kinematics.e_nb
    clm = -vehicle.y.kinematics.v_eb_n[3]
    mode_prev = mdl.y.mode

    mode = (is_on_gnd(vehicle) ? ModeControlLon.direct : mode_req)

    throttle_ref = Float64(throttle_axis + throttle_offset)
    elevator_ref = Float64(elevator_axis + elevator_offset)

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

    mdl.y = ControllerLonY(; mode, throttle_ref, elevator_ref,
        q_ref, θ_ref, EAS_ref, clm_ref, throttle_cmd, elevator_cmd,
        e2e_lqr = e2e_lqr.y, q2e_int = q2e_int.y, q2e_pid = q2e_pid.y,
        v2θ_pid = v2θ_pid.y, c2θ_pid = c2θ_pid.y, v2t_pid = v2t_pid.y)

end


function AircraftBase.assign!(systems::Model{<:C172Y.Systems},
                          mdl::Model{<:ControllerLon})

    @unpack act = systems.submodels
    @unpack throttle_cmd, elevator_cmd = mdl.y

    act.throttle.u[] = throttle_cmd
    act.elevator.u[] = elevator_cmd

end


############################# Initialization ###################################

function Modeling.init!(lon::Model{<:ControllerLon},
                        vehicle::Model{<:C172Y.Vehicle})

    #we assume that the vehicle's y has already been updated to its trim
    #value by init!(vehicle, params)
    y_act = vehicle.y.systems.act
    @unpack ω_wb_b, v_eb_n, e_nb = vehicle.y.kinematics
    @unpack EAS = vehicle.y.airflow

    #reset all controller submodels
    Control.reset!(lon)

    #make inputs consistent with the vehicle status, so
    lon.u.throttle_axis = y_act.throttle.pos
    lon.u.elevator_axis = y_act.elevator.pos
    lon.u.throttle_offset = 0
    lon.u.elevator_offset = 0
    lon.u.q_ref = ω_wb_b[2]
    lon.u.θ_ref = e_nb.θ
    lon.u.EAS_ref = EAS
    lon.u.clm_ref = -v_eb_n[3]

    #for the trim condition to be preserved when the simulation is started with
    #sas (rather than direct) modes enabled, we need a post-trim update of the
    #controller with the inner LQR SAS loops enabled. this loads the z_trim,
    #u_trim and x_trim values corresponding to the trim state into the LQR
    #trackers, and runs them once. after this, their outputs will match the
    #actuator commands required by the trim condition. this match is only
    #approximate, because in general, the trim values loaded from the lookup
    #tables will be interpolated, rather than exactly computed at specific
    #controller design points, but it is good enough. without this step, all
    #modes that build upon the SAS will break the trim equilibrium when selected
    lon.u.mode_req = ModeControlLon.sas
    f_periodic!(lon, vehicle)

    #restore direct modes
    lon.u.mode_req = ModeControlLon.direct
    f_periodic!(lon, vehicle)

end

function GUI.draw!(lon::Model{<:ControllerLon},
                    vehicle::Model{<:C172Y.Vehicle})

    @unpack u, y, Δt = lon

    @unpack systems, kinematics, dynamics, airflow = vehicle.y
    @unpack e_nb, ω_wb_b, v_eb_n = kinematics
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
                    mode_button("Pitch SAS", ModeControlLon.sas, lon.u.mode_req, y.mode)
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
                    IsItemActive() && (lon.u.mode_req = ModeControlLon.EAS_clm; lon.u.EAS_ref = EAS; lon.u.clm_ref = clm)
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

function draw_internals(lon::Model{<:ControllerLon}, p_open::Ref{Bool} = Ref(true))
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

@kwdef struct ControllerLatY
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

        ar2ar_lqr.u.x .= XLat(vehicle)
        ar2ar_lqr.u.z .= ULat(; aileron_cmd = systems.act.aileron.cmd,
                                rudder_cmd = systems.act.rudder.cmd)
        ar2ar_lqr.u.z_ref .= ULat(; aileron_cmd, rudder_cmd)
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

    mdl.y = ControllerLatY(; mode, aileron_ref, rudder_ref,
        p_ref, β_ref, φ_ref, χ_ref, aileron_cmd, rudder_cmd,
        ar2ar_lqr = ar2ar_lqr.y, φβ2ar_lqr = φβ2ar_lqr.y,
        p2φ_int = p2φ_int.y, p2φ_pid = p2φ_pid.y, χ2φ_pid = χ2φ_pid.y)

end


function AircraftBase.assign!(systems::Model{<:C172Y.Systems},
                          mdl::Model{<:ControllerLat})

    @unpack act = systems.submodels
    @unpack aileron_cmd, rudder_cmd = mdl.y

    act.aileron.u[] = aileron_cmd
    act.rudder.u[] = rudder_cmd

end


function Modeling.init!(lat::Model{<:ControllerLat},
                        vehicle::Model{<:C172Y.Vehicle})

    #we assume that the vehicle's y has already been updated to its trim value
    #by init!(vehicle, params)
    y_act = vehicle.y.systems.act
    @unpack ω_wb_b, v_eb_n, e_nb, χ_gnd, ϕ_λ, h_e = vehicle.y.kinematics
    @unpack EAS = vehicle.y.airflow
    @unpack β = vehicle.y.systems.aero

    #reset all controller submodels
    Control.reset!(lat)

    #make Controller inputs consistent with vehicle status
    lat.u.aileron_axis = y_act.aileron.pos
    lat.u.rudder_axis = y_act.rudder.pos
    lat.u.aileron_offset = 0
    lat.u.rudder_offset = 0
    lat.u.p_ref = ω_wb_b[1]
    lat.u.φ_ref = e_nb.φ
    lat.u.β_ref = β
    lat.u.χ_ref = χ_gnd

    #for the trim condition to be preserved when the simulation is started with
    #sas (rather than direct) modes enabled, we need a post-trim update of the
    #controller with the inner LQR SAS loops enabled. this loads the z_trim,
    #u_trim and x_trim values corresponding to the trim state into the LQR
    #trackers, and runs them once. after this, their outputs will match the
    #actuator commands required by the trim condition. this match is only
    #approximate, because in general, the trim values loaded from the lookup
    #tables will be interpolated, rather than exactly computed at specific
    #controller design points, but it is good enough.
    lat.u.mode_req = ModeControlLat.sas
    f_periodic!(lat, vehicle)

    lat.u.mode_req = ModeControlLat.ModeControlLat.φ_β
    f_periodic!(lat, vehicle)

    #restore direct mode
    lat.u.mode_req = ModeControlLat.direct
    f_periodic!(lat, vehicle)

end


function GUI.draw!(lat::Model{<:ControllerLat},
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


function draw_internals(lat::Model{<:ControllerLat}, p_open::Ref{Bool} = Ref(true))
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
################################# Controller ##################################

@kwdef struct Controller{C1 <: ControllerLon, C2 <: ControllerLat} <: AbstractAvionics
    lon::C1 = ControllerLon()
    lat::C2 = ControllerLat()
end

#define fields with submodel's names to create parent-child input linkage
@kwdef mutable struct ControllerU
    lon = ControllerLonU()
    lat = ControllerLatU()
end

Modeling.U(::Controller) = ControllerU()

@no_ode Controller
@no_step Controller
@sm_periodic Controller


########################### Update Methods #####################################

function AircraftBase.assign!(systems::Model{<:C172Y.Systems}, mdl::Model{<:Controller})
    AircraftBase.assign!(systems, mdl.lon)
    AircraftBase.assign!(systems, mdl.lat)
end

############################# Initialization ###################################

function Modeling.init!(mdl::Model{<:Controller}, vehicle::Model{<:C172Y.Vehicle})
    Modeling.init!(mdl.lon, vehicle)
    Modeling.init!(mdl.lat, vehicle)
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

#now we can do:
# JSON3.read(JSON3.write(ControllerU()), ControllerU)
# JSON3.read!(JSON3.write(ControllerU()), ControllerU())


################################################################################
################################## GUI #########################################



function GUI.draw!(ctl::Model{<:Controller},
                    vehicle::Model{<:C172Y.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna172Y Control Laws")

    @unpack u, y, Δt, submodels = ctl
    @unpack lon, lat = submodels

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

    GUI.draw!(lon, vehicle)
    GUI.draw!(lat, vehicle)


    if CImGui.CollapsingHeader("Flight Data")

        if BeginTable("Flight Data", 2, CImGui.ImGuiTableFlags_SizingStretchSame | CImGui.ImGuiTableFlags_BordersInner)
            TableNextRow()
                TableNextColumn();

                Text("Airspeed (Calibrated)"); SameLine(240)
                Text(@sprintf("%.3f m/s | %.3f kts", CAS, Atmosphere.SI2kts(CAS)))
                Text("Airspeed (Equivalent)"); SameLine(240)
                Text(@sprintf("%.3f m/s | %.3f kts", EAS, Atmosphere.SI2kts(EAS)))
                Text("Airspeed (True)"); SameLine(240)
                Text(@sprintf("%.3f m/s | %.3f kts", TAS, Atmosphere.SI2kts(TAS)))
                Text("Angle of Attack"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(α)))
                Text("Sideslip Angle"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(β)))

                Separator()

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


            TableNextColumn();

                Text("Latitude"); SameLine(240)
                Text(@sprintf("%.6f deg", rad2deg(ϕ)))
                Text("Longitude"); SameLine(240)
                Text(@sprintf("%.6f deg", rad2deg(λ)))
                Text("Altitude (Ellipsoidal)"); SameLine(240)
                Text(@sprintf("%.3f m | %.3f ft", Float64(h_e), Float64(h_e)/0.3048))
                Text("Altitude (Orthometric)"); SameLine(240)
                Text(@sprintf("%.3f m | %.3f ft", Float64(h_o), Float64(h_o)/0.3048))
                # Text("Height Over Ground"); SameLine(240)
                # Text(@sprintf("%.3f m | %.3f ft", hog, hog/0.3048))

                Separator()

                Text("Ground Speed"); SameLine(240)
                Text(@sprintf("%.3f m/s | %.3f kts", v_gnd, Atmosphere.SI2kts(v_gnd)))
                Text("Course Angle"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(χ_gnd)))
                Text("Flight Path Angle"); SameLine(240)
                Text(@sprintf("%.3f deg", rad2deg(γ_gnd)))
                Text("Climb Rate"); SameLine(240)
                Text(@sprintf("%.3f m/s", clm))

                Separator()

                Text("Specific Force (x)"); SameLine(240)
                Text(@sprintf("%.3f g", dynamics.f_c_c[1]/Dynamics.g₀))
                Text("Specific Force (y)"); SameLine(240)
                Text(@sprintf("%.3f g", dynamics.f_c_c[2]/Dynamics.g₀))
                Text("Specific Force (z)"); SameLine(240)
                Text(@sprintf("%.3f g", dynamics.f_c_c[3]/Dynamics.g₀))

            EndTable()
        end

        Separator()
    end

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