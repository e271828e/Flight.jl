module C172FBWCAS

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays, HDF5, Interpolations

using Flight.FlightCore
using Flight.FlightCore.Utils

using Flight.FlightPhysics

using Flight.FlightComponents
using Flight.FlightComponents.Control.Discrete: Integrator, IntegratorOutput, PID, PIDOutput, PIDParams
import Flight.FlightComponents.Control.PIDOpt

using ...C172
using ..C172FBW

export Cessna172FBWCAS


################################################################################
################################ Lookup ########################################

struct Lookup{T <: Interpolations.Extrapolation}
    interps::PIDParams{T}
    data::PIDParams{Array{Float64,2}}
    EAS_bounds::NTuple{2, Float64}
    h_bounds::NTuple{2, Float64}
end

function Lookup(data::PIDParams{Array{Float64, 2}},
                EAS_bounds::NTuple{2, Float64},
                h_bounds::NTuple{2, Float64})

    EAS_length, h_length = size(data.k_p)
    EAS_range = range(EAS_bounds..., length = EAS_length)
    h_range = range(h_bounds..., length = h_length)

    (EAS, h) =  map((EAS_range, h_range)) do range
        (mode, scaling) = length(range) > 1 ? (BSpline(Linear()), range) : (NoInterp(), 1:1)
        return (mode = mode, scaling = scaling)
    end

    interps = [extrapolate(scale(interpolate(coef, (EAS.mode, h.mode)),
        EAS.scaling, h.scaling), (Flat(), Flat())) for coef in NamedTuple(data)]

    Lookup(PIDParams(interps...), data, EAS_bounds, h_bounds)

end

Base.getproperty(lookup::Lookup, s::Symbol) = getproperty(lookup, Val(s))
@generated function Base.getproperty(lookup::Lookup, ::Val{S}) where {S}
    if S ∈ fieldnames(Lookup)
        return :(getfield(lookup, $(QuoteNode(S))))
    elseif S ∈ fieldnames(PIDParams)
        return :(getfield(getfield(lookup, :interps), $(QuoteNode(S))))
    else
        error("Lookup has no property $S")
    end
end

function Control.Discrete.PIDParams(lookup::Lookup, EAS::Real, h::Real)

    @unpack k_p, k_i, k_d, τ_f = lookup.interps
    PIDParams(; k_p = k_p(EAS, h),
                k_i = k_i(EAS, h),
                k_d = k_d(EAS, h),
                τ_f = τ_f(EAS, h))
end


(lookup::Lookup)(EAS::Real, h::Real) = PIDParams(lookup, EAS, h)

function save_lookup(lookup::Lookup, fname::String)

    fid = h5open(fname, "w")

    foreach(pairs(NamedTuple(lookup.data))) do (name, data)
        fid[string(name)] = data
    end

    fid["EAS_start"], fid["EAS_end"] = lookup.EAS_bounds
    fid["h_start"], fid["h_end"] = lookup.h_bounds

    close(fid)
end

function load_lookup(fname::String)

    fid = h5open(fname, "r")

    params_data = map(fieldnames(PIDParams)) do name
        read(fid, string(name))
    end
    data = PIDParams(params_data...)

    EAS_bounds = (read(fid["EAS_start"]), read(fid["EAS_end"]))
    h_bounds = (read(fid["h_start"]), read(fid["h_end"]))

    close(fid)

    return Lookup(data, EAS_bounds, h_bounds)

end

################################################################################
########################## AbstractControlChannel ##############################

abstract type AbstractControlChannel <: SystemDefinition end

################################################################################
################################## ThrottleControl #############################

@enum ThrottleMode begin
    direct_throttle_mode = 0
    EAS_throttle_mode = 1
end

@kwdef struct ThrottleControl{L <: Lookup} <: AbstractControlChannel
    v2t_lookup::L = load_lookup(joinpath(@__DIR__, "data", "v2t_lookup.h5"))
    v2t::PID = PID()
end

@kwdef mutable struct ThrottleControlU
    mode::ThrottleMode = direct_throttle_mode #throttle control mode
    thr_dmd::Ranged{Float64, 0., 1.} = 0.0 #throttle actuation demand
    EAS_dmd::Float64 = 0.0 #equivalent airspeed demand
end

@kwdef struct ThrottleControlY
    mode::ThrottleMode = direct_throttle_mode
    thr_dmd::Ranged{Float64, 0., 1.} = 0.0
    EAS_dmd::Float64 = 0.0
    thr_cmd::Ranged{Float64, 0., 1.} = 0.0 #throttle actuation command
    thr_sat::Int64 = 0 #throttle saturation state
    v2t::PIDOutput = PIDOutput()
end

Systems.init(::SystemU, ::ThrottleControl) = ThrottleControlU()
Systems.init(::SystemY, ::ThrottleControl) = ThrottleControlY()

function Systems.init!(sys::System{<:ThrottleControl})
    #set throttle command limits (not strictly necessary, since the
    #compensator's output is converted to a Ranged type downstream anyway)
    v2t = sys.v2t
    v2t.u.bound_lo = 0
    v2t.u.bound_hi = 1
end

function Control.Discrete.reset!(sys::System{<:ThrottleControl})
    #set default inputs and states
    sys.u.mode = direct_throttle_mode
    sys.u.thr_dmd = 0
    sys.u.EAS_dmd = 0
    #reset compensators
    Control.Discrete.reset!.(values(sys.subsystems))
    #set default outputs
    sys.y = ThrottleControlY()
end

function Systems.f_disc!(sys::System{<:ThrottleControl}, kin::KinematicData, air::AirData, Δt::Real)

    @unpack mode, thr_dmd, EAS_dmd = sys.u
    @unpack v2t = sys.subsystems
    @unpack v2t_lookup = sys.constants

    Control.Discrete.assign!(v2t, v2t_lookup(air.EAS, Float64(kin.h_o)))

    if mode === direct_throttle_mode
        thr_cmd = thr_dmd
    else #equivalent airspeed mode
        v2t.u.input = EAS_dmd - air.EAS
        f_disc!(v2t, Δt)
        thr_cmd = Ranged(v2t.y.output, 0., 1.) #compensator sign inversion
    end

    thr_sat = saturation(thr_cmd)
    v2t.u.sat_ext = thr_sat #saturation must also be inverted!

    sys.y = ThrottleControlY(; mode, thr_dmd, EAS_dmd,
                            thr_cmd, thr_sat, v2t = v2t.y)

end

function GUI.draw(sys::System{<:ThrottleControl})

    @unpack v2t = sys.subsystems
    @unpack mode, thr_dmd, EAS_dmd, thr_cmd, thr_sat = sys.y

    CImGui.Begin("Throttle Control")

    CImGui.Text("Mode: $mode")
    CImGui.Text(@sprintf("Throttle demand: %.3f", Float64(thr_dmd)))
    CImGui.Text(@sprintf("EAS demand: %.3f m/s", rad2deg(EAS_dmd)))

    CImGui.Text(@sprintf("Throttle command: %.3f", Float64(thr_cmd)))
    CImGui.Text("Throttle saturation: $thr_sat")

    if @cstatic check=false @c CImGui.Checkbox("EAS Compensator", &check)
        CImGui.Begin("EAS Compensator"); GUI.draw(v2t); CImGui.End()
    end

    CImGui.End()

end

################################################################################
############################### PitchControl ################################

@enum PitchMode begin
    direct_elevator_mode = 0
    pitch_rate_mode = 1
    pitch_angle_mode = 2
    climb_rate_mode = 3
    EAS_pitch_mode = 4
end

################################################################################

@kwdef struct PitchControl{L <: Lookup} <: AbstractControlChannel
    q2e_lookup::L = load_lookup(joinpath(@__DIR__, "data", "q2e_lookup.h5"))
    θ2q_lookup::L = load_lookup(joinpath(@__DIR__, "data", "θ2q_lookup.h5"))
    c2θ_lookup::L = load_lookup(joinpath(@__DIR__, "data", "c2θ_lookup.h5"))
    v2θ_lookup::L = load_lookup(joinpath(@__DIR__, "data", "v2θ_lookup.h5"))
    q2e_int::Integrator = Integrator() #pitch rate integrator
    q2e::PID = PID() #pitch rate to elevator compensator
    θ2q::PID = PID() #pitch angle to pitch rate compensator
    c2θ::PID = PID() #climb rate to pitch angle compensator
    v2θ::PID = PID() #EAS to pitch angle compensator
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct PitchControlU
    mode::PitchMode = direct_elevator_mode #pitch control mode
    e_dmd::Ranged{Float64, -1., 1.} = 0.0 #elevator actuation demand
    q_dmd::Float64 = 0.0 #pitch rate demand
    θ_dmd::Float64 = 0.0 #pitch angle demand
    c_dmd::Float64 = 0.0 #climb rate demand
    EAS_dmd::Float64 = 0.0 #EAS demand
end

@kwdef struct PitchControlY
    mode::PitchMode = direct_elevator_mode
    e_dmd::Ranged{Float64, -1., 1.} = 0.0
    q_dmd::Float64 = 0.0
    θ_dmd::Float64 = 0.0
    c_dmd::Float64 = 0.0
    EAS_dmd::Float64 = 0.0
    e_cmd::Ranged{Float64, -1., 1.} = 0.0 #elevator actuation command
    e_sat::Int64 = 0 #elevator saturation state
    q2e_int::IntegratorOutput = IntegratorOutput()
    q2e::PIDOutput = PIDOutput()
    θ2q::PIDOutput = PIDOutput()
    c2θ::PIDOutput = PIDOutput()
    v2θ::PIDOutput = PIDOutput()
end

Systems.init(::SystemU, ::PitchControl) = PitchControlU()
Systems.init(::SystemY, ::PitchControl) = PitchControlY()

function Systems.init!(sys::System{<:PitchControl})
    @unpack q2e, θ2q, c2θ, v2θ = sys.subsystems
    #gains not set here, they will be assigned by gain scheduling
    q2e.u.bound_lo = -1
    q2e.u.bound_hi = 1
end

function Control.Discrete.reset!(sys::System{<:PitchControl})
    #set default inputs and states
    sys.u.mode = direct_elevator_mode
    sys.u.e_dmd = 0
    sys.u.q_dmd = 0
    sys.u.θ_dmd = 0
    sys.u.c_dmd = 0
    #reset compensators
    Control.Discrete.reset!.(values(sys.subsystems))
    #set default outputs
    sys.y = PitchControlY()
end

function Systems.f_disc!(sys::System{<:PitchControl}, kin::KinematicData, air::AirData, Δt::Real)

    @unpack mode, e_dmd, q_dmd, θ_dmd, c_dmd, EAS_dmd = sys.u
    @unpack q2e_int, q2e, θ2q, c2θ, v2θ = sys.subsystems
    @unpack q2e_lookup, θ2q_lookup, c2θ_lookup, v2θ_lookup = sys.constants

    Control.Discrete.assign!(q2e, q2e_lookup(air.EAS, Float64(kin.h_o)))
    Control.Discrete.assign!(θ2q, θ2q_lookup(air.EAS, Float64(kin.h_o)))
    Control.Discrete.assign!(c2θ, c2θ_lookup(air.EAS, Float64(kin.h_o)))
    Control.Discrete.assign!(v2θ, v2θ_lookup(air.EAS, Float64(kin.h_o)))

    _, q, r = kin.ω_lb_b
    @unpack θ, φ = kin.e_nb
    c = -kin.v_eOb_n[3] #climb rate
    EAS = air.EAS

    if mode === direct_elevator_mode
        e_cmd = e_dmd
    else
        if mode === pitch_rate_mode
            q2e_int.u.input = q_dmd - q
        else #pitch_angle, climb_rate, EAS
            if mode === pitch_angle_mode
                θ2q.u.input = θ_dmd - θ
            else #climb rate to θ
                if mode === climb_rate_mode
                    c2θ.u.input = c_dmd - c
                    f_disc!(c2θ, Δt)
                    θ2q.u.input = c2θ.y.output - θ
                else #EAS to θ
                    v2θ.u.input = EAS_dmd - EAS
                    f_disc!(v2θ, Δt)
                    θ2q.u.input = -v2θ.y.output - θ #sign inversion!
                end
            end
            f_disc!(θ2q, Δt)
            θ_dot_dmd = θ2q.y.output
            φ_bnd = clamp(φ, -π/3, π/3)
            q_θ_dmd = 1/cos(φ_bnd) * θ_dot_dmd + r * tan(φ_bnd)
            q2e_int.u.input = q_θ_dmd - q
        end
        f_disc!(q2e_int, Δt)
        q2e.u.input = q2e_int.y.output
        f_disc!(q2e, Δt)
        e_cmd = Ranged(q2e.y.output, -1., 1.)
    end

    e_sat = saturation(e_cmd)

    #will take effect on the next call
    q2e_int.u.sat_ext = e_sat
    q2e.u.sat_ext = e_sat
    θ2q.u.sat_ext = e_sat
    c2θ.u.sat_ext = e_sat
    v2θ.u.sat_ext = -e_sat #sign inversion!

    sys.y = PitchControlY(; mode, e_dmd, q_dmd, θ_dmd, c_dmd, EAS_dmd,
                            e_cmd, e_sat, q2e_int = q2e_int.y, q2e = q2e.y,
                            θ2q = θ2q.y, c2θ = c2θ.y, v2θ = v2θ.y)

end


function GUI.draw(sys::System{<:PitchControl})

    @unpack q2e_int, q2e, θ2q, c2θ, v2θ = sys.subsystems
    @unpack mode, e_dmd, q_dmd, θ_dmd, c_dmd, EAS_dmd, e_cmd, e_sat = sys.y

    CImGui.Begin("Pitch Control")

    CImGui.Text("Mode: $mode")
    CImGui.Text(@sprintf("Elevator demand: %.3f", Float64(e_dmd)))
    CImGui.Text(@sprintf("Pitch rate demand: %.3f deg/s", rad2deg(q_dmd)))
    CImGui.Text(@sprintf("Pitch angle demand: %.3f deg", rad2deg(θ_dmd)))
    CImGui.Text(@sprintf("Climb rate demand: %.3f m/s", c_dmd))
    CImGui.Text(@sprintf("EAS demand: %.3f m/s", EAS_dmd))

    CImGui.Text(@sprintf("Elevator command: %.3f", Float64(e_cmd)))
    CImGui.Text("Elevator saturation: $e_sat")

    if @cstatic check=false @c CImGui.Checkbox("Pitch Rate Integrator", &check)
        CImGui.Begin("Pitch Rate Integrator"); GUI.draw(q2e_int); CImGui.End()
    end

    if @cstatic check=false @c CImGui.Checkbox("Pitch Rate Compensator", &check)
        CImGui.Begin("Pitch Rate Compensator"); GUI.draw(q2e); CImGui.End()
    end

    if @cstatic check=false @c CImGui.Checkbox("Pitch Angle Compensator", &check)
        CImGui.Begin("Pitch Angle Compensator"); GUI.draw(θ2q); CImGui.End()
    end

    if @cstatic check=false @c CImGui.Checkbox("Climb Rate Compensator", &check)
        CImGui.Begin("Climb Rate Compensator"); GUI.draw(c2θ); CImGui.End()
    end

    if @cstatic check=false @c CImGui.Checkbox("EAS Compensator", &check)
        CImGui.Begin("EAS Compensator"); GUI.draw(v2θ); CImGui.End()
    end

    CImGui.End()

end

################################################################################
################################## RollControl ###############################

@enum RollMode begin
    direct_aileron_mode = 0
    roll_rate_mode = 1
    bank_angle_mode = 2
    course_angle_mode = 3
end

@kwdef struct RollControl{L <: Lookup} <: AbstractControlChannel
    p2a_lookup::L = load_lookup(joinpath(@__DIR__, "data", "p2a_lookup.h5"))
    φ2p_lookup::L = load_lookup(joinpath(@__DIR__, "data", "φ2p_lookup.h5"))
    χ2φ_lookup::L = load_lookup(joinpath(@__DIR__, "data", "χ2φ_lookup.h5"))
    p2a::PID = PID()
    φ2p::PID = PID()
    χ2φ::PID = PID()
end

@kwdef mutable struct RollControlU
    mode::RollMode = direct_aileron_mode
    a_dmd::Ranged{Float64, -1., 1.} = 0.0 #aileron actuation demand
    p_dmd::Float64 = 0.0 #roll rate demand
    φ_dmd::Float64 = 0.0 #bank angle demand
    χ_dmd::Float64 = 0.0 #course angle demand
end

@kwdef struct RollControlY
    mode::RollMode = direct_aileron_mode
    a_dmd::Ranged{Float64, -1., 1.} = Ranged(0.0, -1.0, 1.0) #aileron actuation demand
    p_dmd::Float64 = 0.0
    φ_dmd::Float64 = 0.0
    χ_dmd::Float64 = 0.0
    a_cmd::Ranged{Float64, -1., 1.} = Ranged(0.0, -1.0, 1.0) #aileron actuation command
    a_sat::Int64 = 0 #aileron saturation state
    p2a::PIDOutput = PIDOutput()
    φ2p::PIDOutput = PIDOutput()
    χ2φ::PIDOutput = PIDOutput()
end

Systems.init(::SystemU, ::RollControl) = RollControlU()
Systems.init(::SystemY, ::RollControl) = RollControlY()

function Systems.init!(sys::System{<:RollControl})
    @unpack p2a, φ2p, χ2φ = sys.subsystems

    #gains not set here, they will be assigned by gain scheduling
    p2a.u.bound_lo = -1
    p2a.u.bound_hi = 1

    #set φ demand limits for the course angle compensator output
    χ2φ.u.bound_lo = -π/4
    χ2φ.u.bound_hi = π/4
end

function Control.Discrete.reset!(sys::System{<:RollControl})
    #set default inputs and states
    sys.u.mode = direct_aileron_mode
    sys.u.a_dmd = 0
    sys.u.p_dmd = 0
    sys.u.φ_dmd = 0
    sys.u.χ_dmd = 0
    #reset compensators
    Control.Discrete.reset!.(values(sys.subsystems))
    #set default outputs
    sys.y = RollControlY()
end

function Systems.f_disc!(sys::System{<:RollControl}, kin::KinematicData, air::AirData, Δt::Real)

    @unpack mode, a_dmd, p_dmd, φ_dmd, χ_dmd = sys.u
    @unpack p2a, φ2p, χ2φ = sys.subsystems
    @unpack p2a_lookup, φ2p_lookup, χ2φ_lookup = sys.constants

    Control.Discrete.assign!(p2a, p2a_lookup(air.EAS, Float64(kin.h_o)))
    Control.Discrete.assign!(φ2p, φ2p_lookup(air.EAS, Float64(kin.h_o)))
    Control.Discrete.assign!(χ2φ, χ2φ_lookup(air.EAS, Float64(kin.h_o)))

    p = kin.ω_lb_b[1]
    φ = kin.e_nb.φ
    χ = kin.χ_gnd

    if mode === direct_aileron_mode
        a_cmd = a_dmd
    else
        if mode === roll_rate_mode
            p2a.u.input = p_dmd - p
        else #bank angle, course angle
            if mode === bank_angle_mode
                φ2p.u.input = φ_dmd - φ
            else #course angle
                χ2φ.u.input = wrap_to_π(χ_dmd - χ)
                f_disc!(χ2φ, Δt)
                φ2p.u.input = χ2φ.y.output - φ
            end
            f_disc!(φ2p, Δt)
            p2a.u.input = φ2p.y.output - p
        end
        f_disc!(p2a, Δt)
        a_cmd = Ranged(p2a.y.output, -1., 1.)
    end

    a_sat = saturation(a_cmd)

    #will take effect on the next call
    p2a.u.sat_ext = a_sat
    φ2p.u.sat_ext = a_sat
    χ2φ.u.sat_ext = a_sat

    sys.y = RollControlY(; mode, a_dmd, p_dmd, φ_dmd, χ_dmd, a_cmd, a_sat,
                             p2a = p2a.y, φ2p = φ2p.y, χ2φ = χ2φ.y)

end


function GUI.draw(sys::System{<:RollControl})

    @unpack p2a, φ2p, χ2φ = sys.subsystems
    @unpack mode, a_dmd, p_dmd, φ_dmd, χ_dmd, a_cmd, a_sat = sys.y

    CImGui.Begin("Roll Control")

    CImGui.Text("Mode: $mode")
    CImGui.Text(@sprintf("Aileron demand: %.3f", Float64(a_dmd)))
    CImGui.Text(@sprintf("Roll rate demand: %.3f deg/s", rad2deg(p_dmd)))
    CImGui.Text(@sprintf("Bank angle demand: %.3f deg", rad2deg(φ_dmd)))
    CImGui.Text(@sprintf("Track angle demand: %.3f deg", rad2deg(χ_dmd)))
    CImGui.Text(@sprintf("Aileron command: %.3f", Float64(a_cmd)))
    CImGui.Text("Aileron saturation: $a_sat")

    if @cstatic check=false @c CImGui.Checkbox("Roll Rate Compensator", &check)
        CImGui.Begin("Roll Rate Compensator"); GUI.draw(p2a); CImGui.End()
    end

    if @cstatic check=false @c CImGui.Checkbox("Bank Angle Compensator", &check)
        CImGui.Begin("Bank Angle Compensator"); GUI.draw(φ2p); CImGui.End()
    end

    if @cstatic check=false @c CImGui.Checkbox("Course Angle Compensator", &check)
        CImGui.Begin("Course Angle Compensator"); GUI.draw(χ2φ); CImGui.End()
    end

    CImGui.End()

end


################################################################################
############################ Altitude Control ##################################

@enum AltControlState begin
    altitude_acquire = 0
    altitude_hold = 1
end

@enum AltitudeRef begin
    ellipsoidal = 0
    orthometric = 1
end

#we don't need a dynamic compensator for altitude hold, a simple gain is sufficient
@kwdef struct AltControl <: AbstractControlChannel
    k_h2c::Float64 = 0.2
end

@kwdef mutable struct AltControlU
    h_dmd::Float64 = 0.0 #altitude demand
    h_ref::AltitudeRef = ellipsoidal #altitude reference
end

@kwdef struct AltControlY
    state::AltControlState = altitude_acquire
    throttle_mode::ThrottleMode = direct_throttle_mode
    pitch_mode::PitchMode = pitch_angle_mode
    thr_dmd::Float64 = 0.0
    c_dmd::Float64 = 0.0
end

Systems.init(::SystemU, ::AltControl) = AltControlU()
Systems.init(::SystemY, ::AltControl) = AltControlY()

function Systems.f_disc!(sys::System{<:AltControl}, kin::KinematicData, ::AirData, Δt::Real)

    @unpack h_dmd, h_ref = sys.u
    @unpack k_h2c = sys.constants

    h_threshold = 20 #within 20 m we switch to altitude_hold
    h = (h_ref === ellipsoidal) ? Float64(kin.h_e) : Float64(kin.h_o)
    state = abs(h_dmd - h) > h_threshold ? altitude_acquire : altitude_hold

    if state === altitude_acquire

        throttle_mode = direct_throttle_mode
        pitch_mode = EAS_pitch_mode #EAS_dmd set independently

        thr_dmd = h_dmd > h ? 1.0 : 0.0 #full throttle to climb, idle to descend
        c_dmd = 0.0 #no effect, v2θ active

    else #altitude_hold

        throttle_mode = EAS_throttle_mode #EAS_dmd set independently
        pitch_mode = climb_rate_mode

        thr_dmd = 0.0 #no effect, v2t active
        c_dmd = k_h2c * (h_dmd - h)

    end

    sys.y = AltControlY(; state, throttle_mode, pitch_mode, thr_dmd, c_dmd)

end

################################################################################
#################################### YawControl ################################

@enum YawMode begin
    direct_rudder_mode = 0
end

@kwdef struct YawControl <: AbstractControlChannel end

@kwdef mutable struct YawControlU
    mode::YawMode = direct_rudder_mode #yaw control mode
    r_dmd::Ranged{Float64, -1., 1.} = 0.0 #aileron actuation demand
end

@kwdef struct YawControlY
    mode::YawMode = direct_rudder_mode
    r_dmd::Ranged{Float64, -1., 1.} = 0.0
    r_cmd::Ranged{Float64, -1., 1.} = 0.0 #rudder actuation command
    r_sat::Int64 = 0 #rudder saturation state
end

Systems.init(::SystemU, ::YawControl) = YawControlU()
Systems.init(::SystemY, ::YawControl) = YawControlY()

function Systems.init!(sys::System{<:YawControl})
end

function Control.Discrete.reset!(sys::System{<:YawControl})
    #set default inputs and states
    sys.u.mode = direct_rudder_mode
    sys.u.r_dmd = 0
    #reset pid
    Control.Discrete.reset!.(values(sys.subsystems))
    #set default outputs
    sys.y = YawControlY()
end

function Systems.f_disc!(sys::System{<:YawControl}, kin::KinematicData, air::AirData, Δt::Real)

    @unpack mode, r_dmd = sys.u

    if mode === direct_rudder_mode
        r_cmd = r_dmd
    else
        r_cmd = r_dmd
    end

    r_sat = saturation(r_cmd)
    sys.y = YawControlY(; mode, r_dmd, r_cmd, r_sat)

end


function GUI.draw(sys::System{<:YawControl})

    @unpack mode, r_dmd, r_cmd, r_sat = sys.y

    CImGui.Begin("Yaw Control")

    CImGui.Text("Mode: $mode")
    CImGui.Text(@sprintf("Rudder demand: %.3f", Float64(r_dmd)))
    CImGui.Text(@sprintf("Rudder command: %.3f", Float64(r_cmd)))
    CImGui.Text("Rudder saturation: $r_sat")

    CImGui.End()

end

##################################################################################
################################## Avionics ######################################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@enum LonMode begin
    lon_mode_semi = 0
    lon_mode_alt = 1
end

@enum LatMode begin
    lat_mode_semi = 0
end

@kwdef struct Avionics <: AbstractAvionics
    throttle_ctl::ThrottleControl = ThrottleControl()
    roll_ctl::RollControl = RollControl()
    pitch_ctl::PitchControl = PitchControl()
    yaw_ctl::YawControl = YawControl()
    alt_ctl::AltControl = AltControl()
end

@kwdef mutable struct Inceptors
    eng_start::Bool = false
    eng_stop::Bool = false
    mixture::Ranged{Float64, 0., 1.} = 0.5
    throttle::Ranged{Float64, 0., 1.} = 0.0 #used in direct_throttle_mode
    roll_input::Ranged{Float64, -1., 1.} = 0.0 #used in aileron_mode and roll_rate_mode
    pitch_input::Ranged{Float64, -1., 1.} = 0.0 #used in direct_elevator_mode and pitch_rate_mode
    yaw_input::Ranged{Float64, -1., 1.} = 0.0 #used in rudder_mode and sideslip_mode
    aileron_cmd_offset::Ranged{Float64, -1., 1.} = 0.0 #only for direct mode
    elevator_cmd_offset::Ranged{Float64, -1., 1.} = 0.0 #only for direct mode
    rudder_cmd_offset::Ranged{Float64, -1., 1.} = 0.0 #only for direct mode
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end

@kwdef mutable struct DigitalInputs
    throttle_mode_sel::ThrottleMode = direct_throttle_mode #selected throttle channel mode
    roll_mode_sel::RollMode = direct_aileron_mode #selected roll channel mode
    pitch_mode_sel::PitchMode = direct_elevator_mode #selected pitch channel mode
    yaw_mode_sel::YawMode = direct_rudder_mode #selected yaw channel mode
    lon_mode_sel::LonMode = lon_mode_semi #selected longitudinal control mode
    lat_mode_sel::LatMode = lat_mode_semi #selected lateral control mode
    EAS_dmd::Float64 = 40.0 #equivalent airspeed demand
    θ_dmd::Float64 = 0.0 #pitch angle demand
    c_dmd::Float64 = 0.0 #climb rate demand
    φ_dmd::Float64 = 0.0 #bank angle demand
    χ_dmd::Float64 = 0.0 #course angle demand
    h_dmd::Float64 = 0.0 #altitude demand
    h_ref::AltitudeRef = ellipsoidal #altitude reference
    p_dmd_sf::Float64 = 1.0 #roll_input to p_dmd scale factor (0.2)
    q_dmd_sf::Float64 = 1.0 #pitch_input to q_dmd scale factor (0.2)
end

@kwdef struct AvionicsU
    inceptors::Inceptors = Inceptors()
    digital::DigitalInputs = DigitalInputs()
end

@kwdef struct AvionicsModing
    flight_phase::FlightPhase = phase_gnd
    throttle_mode::ThrottleMode = direct_throttle_mode
    roll_mode::RollMode = direct_aileron_mode
    pitch_mode::PitchMode = direct_elevator_mode
    yaw_mode::YawMode = direct_rudder_mode
    lon_mode::LonMode = lon_mode_semi
    lat_mode::LatMode = lat_mode_semi
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
    throttle_ctl::ThrottleControlY = ThrottleControlY()
    roll_ctl::RollControlY = RollControlY()
    pitch_ctl::PitchControlY = PitchControlY()
    yaw_ctl::YawControlY = YawControlY()
    alt_ctl::AltControlY = AltControlY()
end

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()
Systems.init(::SystemS, ::Avionics) = nothing #keep subsystems local


# ########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:C172FBWCAS.Avionics},
                        physics::System{<:C172FBW.Physics},
                        Δt::Real)

    @unpack eng_start, eng_stop, mixture, throttle,
            roll_input, pitch_input, yaw_input,
            aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
            flaps, brake_left, brake_right = avionics.u.inceptors

    @unpack throttle_mode_sel, roll_mode_sel, pitch_mode_sel, yaw_mode_sel,
            lon_mode_sel, lat_mode_sel, EAS_dmd, θ_dmd, c_dmd, φ_dmd, χ_dmd, h_dmd,
            h_ref, p_dmd_sf, q_dmd_sf = avionics.u.digital

    @unpack throttle_ctl, roll_ctl, pitch_ctl, yaw_ctl, alt_ctl = avionics.subsystems

    @unpack airframe, air = physics.y
    kinematics = physics.y.kinematics.common

    #direct surface and inner loop demands always come from inceptors
    roll_ctl.u.a_dmd = roll_input + aileron_cmd_offset
    pitch_ctl.u.e_dmd = pitch_input + elevator_cmd_offset
    yaw_ctl.u.r_dmd = yaw_input + rudder_cmd_offset
    throttle_ctl.u.thr_dmd = throttle

    roll_ctl.u.p_dmd = p_dmd_sf * Float64(roll_input)
    pitch_ctl.u.q_dmd = q_dmd_sf * Float64(pitch_input)

    #digital inputs, may be overridden by high level modes (like AltControl)
    throttle_ctl.u.EAS_dmd = EAS_dmd
    pitch_ctl.u.θ_dmd = θ_dmd
    pitch_ctl.u.c_dmd = c_dmd
    pitch_ctl.u.EAS_dmd = EAS_dmd
    roll_ctl.u.φ_dmd = φ_dmd
    roll_ctl.u.χ_dmd = χ_dmd
    alt_ctl.u.h_dmd = h_dmd
    alt_ctl.u.h_ref = h_ref

    any_wow = any(SVector{3}(leg.strut.wow for leg in airframe.ldg))
    flight_phase = any_wow ? phase_gnd : phase_air

    if flight_phase === phase_gnd

        #these are irrelevant on ground, but must be defined in all paths
        lon_mode = lon_mode_semi
        lat_mode = lat_mode_semi

        throttle_ctl.u.mode = direct_throttle_mode
        roll_ctl.u.mode = direct_aileron_mode
        pitch_ctl.u.mode = direct_elevator_mode
        yaw_ctl.u.mode = direct_rudder_mode

    else #air

        lon_mode = lon_mode_sel
        lat_mode = lat_mode_sel

        if lon_mode === lon_mode_semi

            #prioritize airspeed control via throttle
            if (throttle_mode_sel === EAS_throttle_mode) && (pitch_mode_sel === EAS_pitch_mode)
                throttle_ctl.u.mode = EAS_throttle_mode
                pitch_ctl.u.mode = pitch_rate_mode
            else
                throttle_ctl.u.mode = throttle_mode_sel
                pitch_ctl.u.mode = pitch_mode_sel
            end

        else #lon_mode === lon_mode_alt

            #we may need to reset altcontrol on mode change if it uses an integrator
            f_disc!(alt_ctl, kinematics, air, Δt)

            throttle_ctl.u.mode = alt_ctl.y.throttle_mode
            throttle_ctl.u.thr_dmd = alt_ctl.y.thr_dmd

            pitch_ctl.u.mode = alt_ctl.y.pitch_mode
            pitch_ctl.u.c_dmd = alt_ctl.y.c_dmd

        end

        if lat_mode === lat_mode_semi

            roll_ctl.u.mode = roll_mode_sel
            yaw_ctl.u.mode = yaw_mode_sel

        end

    end

    f_disc!(throttle_ctl, kinematics, air, Δt)
    f_disc!(roll_ctl, kinematics, air, Δt)
    f_disc!(pitch_ctl, kinematics, air, Δt)
    f_disc!(yaw_ctl, kinematics, air, Δt)

    throttle_cmd = throttle_ctl.y.thr_cmd
    aileron_cmd = roll_ctl.y.a_cmd
    elevator_cmd = pitch_ctl.y.e_cmd
    rudder_cmd = yaw_ctl.y.r_cmd

    moding = AvionicsModing(;
        flight_phase,
        throttle_mode = throttle_ctl.y.mode,
        roll_mode = roll_ctl.y.mode,
        pitch_mode = pitch_ctl.y.mode,
        yaw_mode = yaw_ctl.y.mode,
        lon_mode,
        lat_mode
      )

    #all signal  except for throttle, roll_input, pitch_input and yaw_input pass through
    actuation = ActuationCommands(; eng_start, eng_stop, mixture,
                throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd,
                flaps, brake_left, brake_right)

    avionics.y = AvionicsY(; moding, actuation,
                            throttle_ctl = throttle_ctl.y,
                            roll_ctl = roll_ctl.y,
                            pitch_ctl = pitch_ctl.y,
                            yaw_ctl = yaw_ctl.y,
                            alt_ctl = alt_ctl.y
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

function GUI.draw!(avionics::System{<:C172FBWCAS.Avionics},
                    physics::System{<:C172FBW.Physics},
                    label::String = "Cessna 172 FBW CAS Avionics")

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
############################# Cessna172FBWCAS ##################################

const Cessna172FBWCAS{K, T} = C172FBW.Template{K, T, C172FBWCAS.Avionics} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain}

function Cessna172FBWCAS(kinematics = LTF(), terrain = HorizontalTerrain())
    C172FBW.Template(kinematics, terrain, C172FBWCAS.Avionics())
end


##################################### Tools ####################################

function Aircraft.trim!(ac::System{<:Cessna172FBWCAS},
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

    u.digital.lon_mode_sel = C172FBWCAS.lon_mode_semi
    u.digital.lat_mode_sel = C172FBWCAS.lat_mode_semi
    u.digital.throttle_mode_sel = C172FBWCAS.direct_throttle_mode
    u.digital.roll_mode_sel = C172FBWCAS.direct_aileron_mode
    u.digital.pitch_mode_sel = C172FBWCAS.direct_elevator_mode
    u.digital.yaw_mode_sel = C172FBWCAS.direct_rudder_mode

    #update avionics outputs
    f_disc!(ac.avionics, 1, ac.physics)

    return result

end

function Aircraft.linearize!(ac::System{<:Cessna172FBWCAS}, args...; kwargs...)
    linearize!(ac.physics, args...; kwargs...)
end


# ############################ Joystick Mappings #################################

function IODevices.assign!(sys::System{<:Cessna172FBWCAS}, joystick::Joystick,
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
        u.aileron_cmd_offset = 0
        u.elevator_cmd_offset = 0
    end

    u.flaps += 0.3333 * was_released(joystick, :switch_down)
    u.flaps -= 0.3333 * was_released(joystick, :switch_up)

end



end #module