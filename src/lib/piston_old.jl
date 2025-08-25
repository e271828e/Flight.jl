module Piston

using Interpolations, StaticArrays, StructArrays, ComponentArrays, UnPack

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightLib.Atmosphere: R

using ..Propellers: AbstractPropeller, Propeller
using ..Control.Continuous: PIVector, PIVectorU, PIVectorY

export PistonEngine, PistonThruster


################################################################################
########################### AbstractPistonEngine ###############################

abstract type AbstractPistonEngine <: ModelDefinition end

const β = ISA_layers[1].β

inHg2Pa(p) = 3386.389p
ft2m(h) = 0.3048h
hp2W(P) = 735.49875P
RPM2radpersec(ω) = ω*π/30
radpersec2RPM(ω) = ω/(π/30)

T_ISA(p) = T_std * (p / p_std) ^ (-β * R / g_std)

#inlet parameter from static pressure assuming ISA conditions
p2δ(p) = (p/p_std) * (T_ISA(p)/T_std)^(-0.5)

#inlet parameter from altitude, assuming ISA conditions
function h2δ(h)
    @unpack p, T = ISAData(HGeop(h))
    p / p_std / √(T / T_std)
end

#by default, engine mass is assumed to be accounted for in the airframe
Dynamics.get_mp_b(::Model{<:AbstractPistonEngine}) = MassProperties()
#its internal angular momentum is neglected
Dynamics.get_wr_b(::Model{<:AbstractPistonEngine}) = Wrench()
#and it receives no direct external wrench
Dynamics.get_hr_b(::Model{<:AbstractPistonEngine}) = zeros(SVector{3})


################################################################################
############################## PistonEngine ####################################

#represents a family of naturally aspirated, fuel-injected aviation engines.
#based on performance data available for the Lycoming IO-360A engine, normalized
#with rated power and rated speed to allow for arbitrary engine sizing
struct PistonEngine{L} <: AbstractPistonEngine
    P_rated::Float64
    ω_rated::Float64
    ω_stall::Float64 #speed below which engine shuts down
    ω_max::Float64 #speed above ω_rated for which power output drops back to zero
    ω_idle::Float64 #target idle speed
    τ_start::Float64 #starter torque
    J::Float64 #equivalent axial moment of inertia of the engine shaft
    idle::PIVector{1} #idle MAP control compensator
    frc::PIVector{1} #friction constraint compensator
    lookup::L #performance lookup table
end

function PistonEngine(;
    P_rated= hp2W(200),
    ω_rated = RPM2radpersec(2700), #IO 360
    ω_stall = RPM2radpersec(300),
    ω_max = RPM2radpersec(3100),
    ω_idle = RPM2radpersec(600),
    τ_start = 40,
    J = 0.05,
    idle = PIVector{1}(),
    frc =  PIVector{1}()
    )

    n_stall = ω_stall / ω_rated
    n_max = ω_max / ω_rated
    lookup = generate_lookup(; n_stall, n_max)

    PistonEngine{typeof(lookup)}(
        P_rated, ω_rated, ω_stall, ω_max, ω_idle, τ_start, J, idle, frc, lookup)
end

@enum EngineState begin
    eng_off = 0
    eng_starting = 1
    eng_running = 2
end

@kwdef mutable struct PistonEngineS
    state::EngineState = eng_off
end

#best economy: mixture around 0.1
#best power: mixture around 0.5
@kwdef mutable struct PistonEngineU
    start::Bool = false
    stop::Bool = false
    throttle::Ranged{Float64, 0., 1.} = 0.0 #throttle setting
    mixture::Ranged{Float64, 0., 1.} = 0.5 #mixture setting
    τ_load::Float64 = 0.0
    J_load::Float64 = 0.0
end

@kwdef struct PistonEngineY
    start::Bool = false #start control
    stop::Bool = false #stop control
    throttle::Float64 = 0.0 #throttle setting
    mixture::Float64 = 0.0 #mixture setting
    state::EngineState = eng_off #engine discrete state
    MAP::Float64 = 0.0 #manifold air pressure
    ω::Float64 = 0.0 #angular velocity (crankshaft)
    τ_shaft::Float64 = 0.0 #shaft output torque
    P_shaft::Float64 = 0.0 #shaft power
    SFC::Float64 = 0.0 #specific fuel consumption
    ṁ::Float64 = 0.0 #fuel consumption
    idle::PIVectorY{1} = PIVectorY{1}()
    frc::PIVectorY{1} = PIVectorY{1}()
end

Modeling.X(eng::PistonEngine) = ComponentVector( ω = 0.0,
                                          idle = Modeling.X(eng.idle),
                                          frc = Modeling.X(eng.frc))
Modeling.U(::PistonEngine) = PistonEngineU()
Modeling.Y(::PistonEngine) = PistonEngineY()
Modeling.S(::PistonEngine) = PistonEngineS()

function Modeling.init!(mdl::Model{<:PistonEngine})
    #set up friction constraint compensator
    @unpack idle, frc = mdl.submodels

    idle.u.k_p .= 4.0
    idle.u.k_i .= 2.0
    idle.u.bound_lo .= -0.5
    idle.u.bound_hi .= 0.5

    frc.u.k_p .= 5.0
    frc.u.k_i .= 200.0
    frc.u.bound_lo .= -1
    frc.u.bound_hi .= 1
end

function Modeling.f_ode!(eng::Model{<:PistonEngine}, air_data::AirData)

    @unpack ω_rated, ω_idle, P_rated, J, τ_start, lookup = eng.constants
    @unpack idle, frc = eng.submodels
    @unpack start, stop, τ_load, J_load = eng.u

    throttle = Float64(eng.u.throttle)
    mixture = Float64(eng.u.mixture)
    state = eng.s.state
    ω = eng.x.ω

    #update friction constraint compensator
    frc.u.input .= -ω
    f_ode!(frc)

    #update idle compensator
    idle.u.input .= 1 - ω / ω_idle #normalized ω error
    f_ode!(idle)

    #normalized MAP at idle throttle from idle compensator output
    μ_ratio_idle = 0.5 + idle.y.output[1]

    #normalized engine speed
    n = ω / ω_rated

    #inlet air parameter for ISA conditions
    δ = p2δ(air_data.p)

    #normalized MAP at wide open throttle
    μ_wot = lookup.μ_wot(n, δ)

    #part throttle normalized MAP
    μ = μ_wot * (μ_ratio_idle + throttle * (1 - μ_ratio_idle))

    @show lookup.π_ratio(mixture)
    @show lookup.sfc_ratio(mixture)

    if state === eng_off

        #with the engine off, introduce a friction constraint to make the
        #propeller actually stop instead of slowing down asymptotically. with
        #the engine running, all friction is assumed to be already accounted for
        #in the performance tables
        τ_fr_max = 0.01 * P_rated / ω_rated #1% of rated torque
        τ_fr = frc.y.output[1] .* τ_fr_max #scale τ_fr_max with compensator output

        MAP = air_data.p
        τ_shaft = τ_fr
        P_shaft = 0.0
        SFC = 0.0
        ṁ = 0.0

    elseif state === eng_starting

        MAP = μ * p_std
        τ_shaft = τ_start
        P_shaft = τ_shaft * ω
        SFC = 0.0
        ṁ = 0.0

    else #state === eng_running

        #normalized power at part throttle, altitude, ISA conditions and maximum
        #power mixture
        π_ISA_pow = compute_π_ISA_pow(lookup, n, μ, δ)

        #correction for non-ISA conditions
        π_pow = π_ISA_pow * √(T_ISA(air_data.p)/air_data.T)

        #correction for arbitrary mixture setting
        π_actual = π_pow * lookup.π_ratio(mixture)

        MAP = μ * p_std
        P_shaft = P_rated * π_actual
        τ_shaft = (ω > 0 ? P_shaft / ω : 0.0) #for ω < ω_stall we should not even be here
        SFC = lookup.sfc_pow(n, π_actual) * lookup.sfc_ratio(mixture)
        ṁ = SFC * P_shaft

    end

    Στ = τ_shaft + τ_load
    ΣJ = J + J_load
    ω_dot = Στ / ΣJ

    eng.ẋ.ω = ω_dot

    eng.y = PistonEngineY(; start, stop, throttle, mixture, state,
                            MAP, ω, τ_shaft, P_shaft, SFC, ṁ,
                            idle = eng.idle.y, frc = eng.frc.y)

end

function Modeling.f_step!(eng::Model{<:PistonEngine}, fuel_available::Bool = true)

    @unpack idle, frc = eng.submodels
    @unpack ω_stall, ω_idle = eng.constants

    ω = eng.x.ω

    if eng.s.state === eng_off

        eng.u.start && (eng.s.state = eng_starting)

    elseif eng.s.state === eng_starting

        !eng.u.start && (eng.s.state = eng_off)
        (ω > ω_idle && fuel_available) && (eng.s.state = eng_running)

    else #eng_running

        (eng.u.stop || ω < ω_stall || !fuel_available) && (eng.s.state = eng_off)

    end

    f_step!(idle)
    f_step!(frc)

end


function generate_lookup(; n_stall, n_max)

    @assert n_stall < 0.667
    @assert n_max > 1.074

    #δ for which a given μ is the wide open throttle μ
    δ_wot = let

        n_range = range(0.667, 1, length = 2)
        μ_range = range(0.401, 0.936, length = 9)

        δ_data = [0.455 0.523 0.587 0.652 0.718 0.781 0.844 0.906 0.965;
                0.464 0.530 0.596 0.662 0.727 0.792 0.855 0.921 0.981]

        extrapolate(scale(interpolate(δ_data, BSpline(Linear())), n_range, μ_range), Line())

    end

    #wide open throttle normalized manifold pressure for a given δ
    μ_wot = let

        n_range = range(0.667, 1, length = 2)
        δ_range = range(0.441, 1, length = 9)

        μ_knots = range(0.401, 0.936, length = 9)
        μ_data = Array{Float64}(undef, length(n_range), length(δ_range))
        for (i,n) in enumerate(n_range)
            #inverse interpolation μ(δ)
            μ_1D = linear_interpolation(δ_wot(n, μ_knots), μ_knots, extrapolation_bc = Line())
            μ_data[i, :] = μ_1D.(δ_range)
        end

        extrapolate(scale(interpolate(μ_data, BSpline(Linear())), n_range, δ_range), Line())

    end

    #part throttle normalized power at sea level for ISA conditions (δ = 1) and maximum power mixture
    π_std = let

        n_data = [n_stall, 0.667, 0.704, 0.741, 0.778, 0.815, 0.852, 0.889, 0.926, 0.963, 1.000, 1.074, n_max]
        μ_data = [0, 0.568, 1.0]

        μ_knots = [zeros(1, length(n_data))
                   0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568
                   1.000 0.836 0.854 0.874 0.898 0.912 0.939 0.961 0.959 0.958 0.956 0.953 1.000]

        #power is zero at MAP = 0 regardless of n (first row)
        #power is zero at n_stall and n_max regardless of MAP (first and last columns)
        π_knots = [zeros(1, length(n_data))
                   0 0.270 0.305 0.335 0.360 0.380 0.405 0.428 0.450 0.476 0.498 0.498 0
                   0 0.489 0.548 0.609 0.680 0.729 0.810 0.880 0.920 0.965 1.000 0.950 0]

        π_data = Array{Float64,2}(undef, (length(n_data), length(μ_data)))

        for i in 1:length(n_data)
            π_std_1D = linear_interpolation(
                μ_knots[:,i], π_knots[:,i], extrapolation_bc = Line())
            π_data[i,:] = π_std_1D.(μ_data)
        end

        linear_interpolation(
            (n_data, μ_data), π_data, extrapolation_bc = ((Flat(), Flat()), (Flat(), Flat())))

    end

    #wide-open throttle normalized power at altitude for ISA conditions and maximum power mixture
    π_wot = let

        n_data = [n_stall, 0.667, 1.000, 1.074, n_max]
        δ_data = [0, 0.441, 1]

        π_data = Array{Float64,2}(undef, length(n_data), length(δ_data))

        π_data[:, 1] .= 0 #power should vanish for δ → 0 ∀n
        π_data[:, 2] .= [0, 0.23, 0.409, 0.409, 0] #it also vanishes at n_stall and n_max ∀δ
        π_data[:, 3] .= [π_std(n, μ_wot(n, 1)) for n in n_data] #at δ=1 by definition it's π_std(μ_wot)

        linear_interpolation((n_data,δ_data), π_data, extrapolation_bc = ((Flat(), Flat()), (Flat(), Line())))

    end

    #actual normalized power over normalized power at maximum power mixture
    π_ratio = let

        m_data = [0.0, 0.071, 0.137, 0.281, 0.402, 0.472, 0.542, 0.637, 0.816, 1.0]
        π_ratio_data = [0.860, 0.931, 0.961, 0.989, 0.999, 1.000, 0.999, 0.994, 0.976, 0.950]

        linear_interpolation(m_data, π_ratio_data, extrapolation_bc = Flat())

    end

    #actual sfc over sfc at maximum power mixture
    sfc_ratio = let

        m_data = [0.0, 0.071, 0.137, 0.281, 0.402, 0.472, 0.542, 0.637, 0.816, 1.0]
        sfc_ratio_data = [0.870, 0.850, 0.854, 0.901, 0.959, 1.0, 1.042, 1.105, 1.243, 1.428]

        linear_interpolation(m_data, sfc_ratio_data, extrapolation_bc = Flat())

    end

    #sfc at maximum power mixture
    sfc_pow = let

        n_data = [2000, 2200, 2400, 2600, 2700] / 2700
        π_data = 10 .^ range(-1, 0, length = 8)

        sfc_data = 1e-7 * [
            1.7671   1.43728  1.19992  1.02909  0.906153  0.817674  0.753997  0.708169
            1.83791  1.49664  1.25103  1.07427  0.947056  0.855503  0.789613  0.742193
            1.98614  1.60588  1.3322   1.13524  0.993496  0.891482  0.818064  0.765226
            2.11663  1.70062  1.40123  1.18576  1.03069   0.919083  0.838765  0.780961
            2.33484  1.85418  1.50825  1.2593   1.08012   0.951177  0.858376  0.791588]

        linear_interpolation((n_data, π_data), sfc_data, extrapolation_bc = Line())

    end

    return (δ_wot = δ_wot, μ_wot = μ_wot, π_std = π_std, π_wot = π_wot,
            π_ratio = π_ratio, sfc_ratio = sfc_ratio, sfc_pow = sfc_pow )

end

function compute_π_ISA_pow(lookup, n, μ, δ)

        δ_wot = lookup.δ_wot(n, μ) #δ at which our μ would be μ_wot

        #normalized power at part throttle, sea level, ISA conditions (δ = 1)
        #and maximum power mixture
        π_std = lookup.π_std(n, μ)

        #normalized power at full throttle, altitude, ISA conditions and maximum
        #power mixture
        π_wot = lookup.π_wot(n, δ_wot)

        if abs(δ_wot - 1) < 5e-3
            π_ISA_pow = π_std
        else
            π_ISA_pow = π_std + (π_wot - π_std) / (δ_wot - 1) * (δ - 1)
        end

        return max(π_ISA_pow, 0)

end


function GUI.draw(mdl::Model{<:PistonEngine}, p_open::Ref{Bool} = Ref(true),
                window_label::String = "Piston Engine")

    @unpack u, y, constants = mdl
    @unpack idle, frc = mdl
    @unpack start, stop, state, throttle, mixture, MAP, ω, τ_shaft, P_shaft, ṁ, SFC = y

    CImGui.Begin(window_label, p_open)

        CImGui.Text("Start: $start")
        CImGui.Text("Stop: $stop")
        CImGui.Text("State: $state")
        CImGui.Text(@sprintf("Throttle: %.3f", throttle))
        CImGui.Text(@sprintf("Mixture: %.3f", mixture))
        CImGui.Text(@sprintf("Manifold Pressure: %.3f Pa", MAP))
        CImGui.Text(@sprintf("Speed: %.3f RPM", radpersec2RPM(ω)))
        CImGui.Text(@sprintf("Shaft Torque: %.3f N*m", τ_shaft))
        CImGui.Text(@sprintf("Shaft Power: %.3f kW", P_shaft/1e3))
        CImGui.Text(@sprintf("Fuel Consumption: %.3f g/s", ṁ*1e3))
        CImGui.Text(@sprintf("Specific Fuel Consumption: %.3f g/(s*kW)", SFC*1e6))

        if CImGui.TreeNode("Idle RPM Controller")
            GUI.draw(idle, window_label)
            CImGui.TreePop()
        end

        if CImGui.TreeNode("Friction Regulator")
            GUI.draw(frc, window_label)
            CImGui.TreePop()
        end

    CImGui.End()

end


# ################################################################################
# ########################## PistonThruster ######################################

#τ_shaft is always positive. for a CW thruster, the gear ratio should be
#positive as well and, under normal operating conditions, τ_prop will be
#negative. for a CCW thruster, the gear ratio should be negative and, under
#normal operating conditions, τ_prop will be positive

#the sign of τ_prop may be inverted under negative propeller thrust conditions,
#with the propeller driving the engine instead of the other way around

@kwdef struct PistonThruster{E <: AbstractPistonEngine,
                        P <: AbstractPropeller} <: ModelDefinition
    engine::E = PistonEngine()
    propeller::P = Propeller()
    gear_ratio::Float64 = 1.0

    function PistonThruster(eng::E, prop::P, gear_ratio) where {E, P}
        @assert(sign(gear_ratio) * Int(prop.sense) > 0,
            "PistonThruster gear ratio sign does not match propeller turn sign")
        new{E,P}(eng, prop, gear_ratio)
    end

end


function Modeling.f_ode!(mdl::Model{<:PistonThruster}, air_data::AirData, kin_data::KinData)

    @unpack engine, propeller = mdl
    @unpack gear_ratio = mdl.constants

    ω_eng = engine.x.ω
    ω_prop = gear_ratio * ω_eng
    f_ode!(propeller, kin_data, air_data, ω_prop)

    τ_prop = propeller.y.wr_p.τ[1]
    τ_eq = gear_ratio * τ_prop #load torque seen from the engine shaft

    J_prop = propeller.J_xx
    J_eq = gear_ratio^2 * J_prop #load moment of inertia seen from the engine side

    engine.u.τ_load = τ_eq
    engine.u.J_load = J_eq
    f_ode!(engine, air_data)

    update_output!(mdl)

end

function Modeling.f_step!(mdl::Model{<:PistonThruster}, fuel_available::Bool = true)

    @unpack engine, propeller = mdl

    f_step!(engine, fuel_available)
    f_step!(propeller)

end

Dynamics.get_mp_b(::Model{<:PistonThruster}) = MassProperties()
Dynamics.get_wr_b(mdl::Model{<:PistonThruster}) = get_wr_b(mdl.propeller) #only external
Dynamics.get_hr_b(mdl::Model{<:PistonThruster}) = get_hr_b(mdl.propeller)


################################################################################
################################# GUI ##########################################

function GUI.draw(mdl::Model{<:PistonThruster}, p_open::Ref{Bool} = Ref(true),
                window_label::String = "Piston Thruster")

    CImGui.Begin(window_label, p_open)
        @cstatic c_eng=false c_prop=false begin
            @c CImGui.Checkbox("Engine", &c_eng)
            @c CImGui.Checkbox("Propeller", &c_prop)
            c_eng && @c GUI.draw(mdl.engine, &c_eng)
            c_prop && @c GUI.draw(mdl.propeller, &c_prop)
        end
    CImGui.End()

end


end #module