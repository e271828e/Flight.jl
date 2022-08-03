module Piston

using Interpolations, Unitful, Plots, StaticArrays, StructArrays, ComponentArrays, UnPack

using Flight.Systems, Flight.Utils
using Flight.Kinematics, Flight.RigidBody, Flight.Atmosphere
using Flight.Atmosphere: ISA_layers, ISAData, p_std, T_std, g_std, R
using Flight.Geodesy: HGeop
using Flight.Propellers: AbstractPropeller, Propeller
using Flight.Friction

import Flight.Systems: init, f_cont!, f_disc!
import Flight.RigidBody: MassTrait, WrenchTrait, AngularMomentumTrait, get_hr_b, get_wr_b
import Flight.RigidBody: get_mp_b

########################### AbstractFuelSupply #################################

abstract type AbstractFuelSupply <: SystemDescriptor end

MassTrait(::System{<:AbstractFuelSupply}) = HasMass()
WrenchTrait(::System{<:AbstractFuelSupply}) = GetsNoExternalWrench()
AngularMomentumTrait(::System{<:AbstractFuelSupply}) = HasNoAngularMomentum()

#to be extended by concrete subtypes
fuel_available(f::System{<:AbstractFuelSupply}) = throw(
    MethodError(fuel_available, (f,)))

#a massless, infinite fuel supply, for testing purposes
struct MagicFuelSupply <: AbstractFuelSupply end

init(::SystemU, ::MagicFuelSupply) = Ref(true)

f_cont!(::System{MagicFuelSupply}) = nothing
f_disc!(::System{MagicFuelSupply}) = false

get_mp_b(::System{MagicFuelSupply}) = MassProperties()
fuel_available(f::System{MagicFuelSupply}) = f.u[]

########################### AbstractPistonEngine ###############################

abstract type AbstractPistonEngine <: SystemDescriptor end

# can't figure out how to register these so that they are seen from module
# methods
# @unit inHg "inHg" InchOfMercury 3386.389*Unitful.Pa false
# @unit hp "hp" Horsepower 735.49875*Unitful.W false

const β = ISA_layers[1].β

inHg2Pa(p) = 3386.389p
ft2m(h) = 0.3048h
hp2W(P) = 735.49875P

T_ISA(p) = T_std * (p / p_std) ^ (-β * R / g_std)

#inlet parameter from static pressure assuming ISA conditions
p2δ(p) = (p/p_std) * (T_ISA(p)/T_std)^(-0.5)

function h2δ(h)
    @unpack p, T = ISAData(HGeop(h))
    p / p_std / √(T / T_std)
end

#by default, engine mass is assumed to be accounted for in the vehicle's
#structure
MassTrait(::System{<:AbstractPistonEngine}) = HasNoMass()
#the propeller gets it instead
WrenchTrait(::System{<:AbstractPistonEngine}) = GetsNoExternalWrench()
#assumed to be negligible by default
AngularMomentumTrait(::System{<:AbstractPistonEngine}) = HasNoAngularMomentum()

########################### IdleController #####################################

#a simple PI controller to maintain the desired target idle RPM
Base.@kwdef struct IdleController <: SystemDescriptor
    k_p::Float64 = 4
    k_i::Float64 = 2
    ω_target::Float64 = 60
end

Base.@kwdef struct IdleControllerY
    ϵ::Float64 = 0.0
    ϵ_int::Float64 = 0.0
    output_raw::Float64 = 0.0
    output::Float64 = 0.0
    sat::Bool = false
end

init(::SystemX, ::IdleController) = [0.0]
init(::SystemY, ::IdleController) = IdleControllerY()

function f_cont!(sys::System{IdleController}, ω::Real)

    @unpack k_p, k_i, ω_target = sys.params

    ϵ = (ω - ω_target)/ω_target
    ϵ_int = sys.x[]
    output_raw = -(k_p * ϵ + k_i * ϵ_int)
    output = min(max(0.0, output_raw), 1.0)
    sat = ((output_raw > 0) && (output_raw < 1) ? false : true)
    sys.ẋ .= !sat * ϵ

    sys.y = IdleControllerY(ϵ, ϵ_int, output_raw, output, sat)

end

f_disc!(::System{IdleController}, args...) = false

############################ Engine ###############################

#represents a family of naturally aspirated, fuel-injected aviation engines.
#based on performance data available for the Lycoming IO-360A engine. data is
#normalized with rated power and rated speed to allow for arbitrary engine
#sizing
struct Engine{L} <: AbstractPistonEngine
    P_rated::Float64
    ω_rated::Float64
    ω_stall::Float64 #speed below which engine shuts down
    ω_cutoff::Float64 #speed above ω_rated for which power output drops back to zero
    M_start::Float64 #starter torque
    J::Float64 #equivalent axial moment of inertia of the engine shaft
    idle::IdleController
    lookup::L
end

function Engine(;
    P_rated= 200 |> hp2W,
    ω_rated = ustrip(u"rad/s", 2700u"rpm"), #IO 360
    ω_stall = ustrip(u"rad/s", 300u"rpm"),
    ω_cutoff = ustrip(u"rad/s", 3100u"rpm"),
    M_start = 40,
    J = 0.05,
    idle = IdleController(ω_target = 2ω_stall),
    )

    n_stall = ω_stall / ω_rated
    n_cutoff = ω_cutoff / ω_rated
    lookup = generate_lookup(; n_stall, n_cutoff)

    Engine{typeof(lookup)}(
        P_rated, ω_rated, ω_stall, ω_cutoff, M_start, J, idle, lookup)
end

@enum EngineState begin
    eng_off = 0
    eng_starting = 1
    eng_running = 2
end

Base.@kwdef mutable struct PistonEngineS
    state::EngineState = eng_off
end

#best economy: mixture around 0.1
#best power: mixture around 0.5
Base.@kwdef mutable struct PistonEngineU
    start::Bool = false
    stop::Bool = false
    thr::Ranged{Float64, 0, 1} = 0.0 #throttle setting
    mix::Ranged{Float64, 0, 1} = 0.5 #mixture setting
end

Base.@kwdef struct PistonEngineY
    start::Bool = false #start control
    stop::Bool = false #stop control
    throttle::Float64 = 0.0 #throttle setting
    mixture::Float64 = 0.0 #mixture setting
    state::EngineState = eng_off #engine discrete state
    MAP::Float64 = 0.0 #manifold air pressure
    ω::Float64 = 0.0 #angular velocity (crankshaft)
    M_shaft::Float64 = 0.0 #shaft output torque
    P_shaft::Float64 = 0.0 #shaft power
    SFC::Float64 = 0.0 #specific fuel consumption
    ṁ::Float64 = 0.0 #fuel consumption
    idle::IdleControllerY = IdleControllerY()
end

init(::SystemX, eng::Engine) = ComponentVector(ω = 0.0, idle = init_x(eng.idle))
init(::SystemU, ::Engine) = PistonEngineU()
init(::SystemY, ::Engine) = PistonEngineY()
init(::SystemS, ::Engine) = PistonEngineS()

function f_cont!(eng::System{<:Engine}, air::AirflowData; M_load::Real, J_load::Real)

    @unpack ω_rated, P_rated, J, M_start, lookup = eng.params
    @unpack thr, mix, start, stop = eng.u
    state = eng.s.state
    ω = eng.x.ω

    throttle = Float64(thr)
    mixture = Float64(mix)

    f_cont!(eng.idle, ω) #update idle controller

    μ_ratio_idle = eng.idle.y.output

    #normalized engine speed
    n = ω / ω_rated

    #inlet air parameter for ISA conditions
    δ = p2δ(air.p)

    #normalized MAP at wide open throttle
    μ_wot = lookup.μ_wot(n, δ)

    #actual, part throttle normalized MAP
    μ = μ_wot * (μ_ratio_idle + throttle * (1 - μ_ratio_idle))

    if state === eng_off

        MAP = air.p
        M_shaft = 0.0
        P_shaft = 0.0
        SFC = 0.0
        ṁ = 0.0

    elseif state === eng_starting

        MAP = μ * p_std
        M_shaft = M_start
        P_shaft = M_shaft * ω
        SFC = 0.0
        ṁ = 0.0

    else #state === eng_running

        #normalized power at part throttle, altitude, ISA conditions and maximum
        #power mixture
        π_ISA_pow = compute_π_ISA_pow(lookup, n, μ, δ)

        #correction for non-ISA conditions
        π_pow = π_ISA_pow * √(T_ISA(air.p)/air.T)

        #correction for arbitrary mixture setting
        π_actual = π_pow * lookup.π_ratio(mixture)

        MAP = μ * p_std
        P_shaft = P_rated * π_actual
        M_shaft = (ω > 0 ? P_shaft / ω : 0.0) #for ω < ω_stall we should not even be here
        SFC = lookup.sfc_pow(n, π_actual) * lookup.sfc_ratio(mixture)
        ṁ = SFC * P_shaft

    end

    ΣM = M_shaft + M_load
    ΣJ = J + J_load
    ω_dot = ΣM / ΣJ

    eng.ẋ.ω = ω_dot

    eng.y = PistonEngineY(; start, stop, throttle, mixture, state,
                            MAP, ω, M_shaft, P_shaft, SFC, ṁ, idle = eng.idle.y)

end

function f_disc!(eng::System{<:Engine}, fuel::System{<:AbstractFuelSupply})

    ω = eng.x.ω
    ω_stall = eng.params.ω_stall
    ω_target = eng.idle.params.ω_target

    if eng.s.state === eng_off

        eng.u.start ? eng.s.state = eng_starting : nothing

    elseif eng.s.state === eng_starting

        !eng.u.start ? eng.s.state = eng_off : nothing

        (ω > ω_target && fuel_available(fuel) ) ? eng.s.state = eng_running : nothing

    else #eng_running

        (ω < ω_stall || eng.u.stop || !fuel_available(fuel)) ? eng.s.state = eng_off : nothing

    end

    x_mod = false
    x_mod = x_mod || f_disc!(eng.idle)
    return x_mod

end


function generate_lookup(; n_stall, n_cutoff)

    @assert n_stall < 1
    @assert n_cutoff > 1

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
            μ_1D = LinearInterpolation(δ_wot(n, μ_knots), μ_knots, extrapolation_bc = Line())
            μ_data[i, :] = μ_1D.(δ_range)
        end

        extrapolate(scale(interpolate(μ_data, BSpline(Linear())), n_range, δ_range), Line())

    end

    #part throttle normalized power at sea level for ISA conditions (δ = 1) and maximum power mixture
    π_std = let

        n_data = [n_stall, 0.667, 0.704, 0.741, 0.778, 0.815, 0.852, 0.889, 0.926, 0.963, 1.000, 1.074, n_cutoff]
        μ_data = [0, 0.568, 1.0]

        μ_knots = [zeros(1, length(n_data))
                   0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568
                   1.000 0.836 0.854 0.874 0.898 0.912 0.939 0.961 0.959 0.958 0.956 0.953 1.000]

        #power is zero at MAP = 0 regardless of n (first row)
        #power is zero at n_stall and n_cutoff regardless of MAP (first and last columns)
        π_knots = [zeros(1, length(n_data))
                   0 0.270 0.305 0.335 0.360 0.380 0.405 0.428 0.450 0.476 0.498 0.498 0
                   0 0.489 0.548 0.609 0.680 0.729 0.810 0.880 0.920 0.965 1.000 0.950 0]

        π_data = Array{Float64,2}(undef, (length(n_data), length(μ_data)))

        for i in 1:length(n_data)
            π_std_1D = LinearInterpolation(
                μ_knots[:,i], π_knots[:,i], extrapolation_bc = Line())
            π_data[i,:] = π_std_1D.(μ_data)
        end

        LinearInterpolation(
            (n_data, μ_data), π_data, extrapolation_bc = ((Flat(), Flat()), (Flat(), Flat())))

    end

    #wide-open throttle normalized power at altitude for ISA conditions and maximum power mixture
    π_wot = let

        n_data = [n_stall, 0.667, 1.000, 1.074, n_cutoff]
        δ_data = [0, 0.441, 1]

        π_data = Array{Float64,2}(undef, length(n_data), length(δ_data))

        π_data[:, 1] .= 0 #power should vanish for δ → 0 ∀n
        π_data[:, 2] .= [0, 0.23, 0.409, 0.409, 0] #it also vanishes at n_stall and n_cutoff ∀δ
        π_data[:, 3] .= [π_std(n, μ_wot(n, 1)) for n in n_data] #at δ=1 by definition it's π_std(μ_wot)

        LinearInterpolation((n_data,δ_data), π_data, extrapolation_bc = ((Flat(), Flat()), (Flat(), Line())))

    end

    #actual normalized power over normalized power at maximum power mixture
    π_ratio = let

        m_data = [0.0, 0.071, 0.137, 0.281, 0.402, 0.472, 0.542, 0.637, 0.816, 1.0]
        π_ratio_data = [0.860, 0.931, 0.961, 0.989, 0.999, 1.000, 0.999, 0.994, 0.976, 0.950]

        LinearInterpolation(m_data, π_ratio_data, extrapolation_bc = Flat())

    end

    #actual sfc over sfc at maximum power mixture
    sfc_ratio = let

        m_data = [0.0, 0.071, 0.137, 0.281, 0.402, 0.472, 0.542, 0.637, 0.816, 1.0]
        sfc_ratio_data = [0.870, 0.850, 0.854, 0.901, 0.959, 1.0, 1.042, 1.105, 1.243, 1.428]

        LinearInterpolation(m_data, sfc_ratio_data, extrapolation_bc = Flat())

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

        LinearInterpolation((n_data, π_data), sfc_data, extrapolation_bc = Line())

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

############################# Thruster ###################################

#M_eng is always positive. therefore, for a CW thruster, n should be positive as
#well and, under normal operating, conditions M_prop will be negative. for a CCW
#thruster, n should be negative and M_prop will be positive under normal
#operating conditions.

#however, the sign of M_prop may be inverted under negative propeller thrust
#conditions (generative), with the propeller driving the engine instead of
#the other way around

Base.@kwdef struct Thruster{E <: AbstractPistonEngine,
                            P <: AbstractPropeller} <: SystemDescriptor
    engine::E = Engine()
    propeller::P = Propeller()
    gear_ratio::Float64 = 1.0 #gear ratio
    friction::Friction.Regulator{1} = Friction.Regulator{1}()
    M_fr_max::Float64 = 5.0 #maximum friction torque
    function Thruster(eng::E, prop::P, gear_ratio, friction, M_fr_max) where {E, P}
        @assert sign(gear_ratio) * Int(prop.sense) > 0 "Thruster gear ratio sign "*
        "does not match propeller turn sign"
        new{E,P}(eng, prop, gear_ratio, friction, M_fr_max)
    end
end


function f_cont!(thr::System{<:Thruster}, air::AirflowData, kin::KinematicData)

    @unpack engine, propeller, friction = thr
    @unpack gear_ratio, M_fr_max = thr.params

    ω_eng = engine.x.ω
    ω_prop = gear_ratio * ω_eng

    f_cont!(friction, SVector{1, Float64}(ω_eng))
    f_cont!(propeller, kin, air, ω_prop)

    M_prop = propeller.y.wr_p.M[1]
    M_eq = gear_ratio * M_prop #load torque seen from the engine shaft
    M_fr = friction.y.α[1] .* M_fr_max

    #when the engine is stopped, introduce a friction constraint to make the
    #propeller actually stop instead of slowing down asymptotically, and once
    #stopped to keep it so (unless there is a lot of headwind...). with the
    #engine running, all friction is assumed to be already accounted for by the
    #performance tables, which give shaft torque and power
    if engine.s.state === eng_off
        M_eq += M_fr
    end

    J_prop = propeller.params.J_xx
    J_eq = gear_ratio^2 * J_prop #load moment of inertia seen from the engine side

    f_cont!(engine, air; M_load = M_eq, J_load = J_eq)

    Systems.assemble_y!(thr)

end

function f_disc!(thr::System{<:Thruster}, fuel::System{<:AbstractFuelSupply})

    @unpack engine, propeller, friction = thr

    x_mod = false
    x_mod = x_mod || f_disc!(engine, fuel)
    x_mod = x_mod || f_disc!(propeller)
    x_mod = x_mod || f_disc!(friction)
    return x_mod

end

MassTrait(::System{<:Thruster}) = HasNoMass()
WrenchTrait(::System{<:Thruster}) = GetsExternalWrench()
AngularMomentumTrait(::System{<:Thruster}) = HasAngularMomentum()

get_wr_b(thr::System{<:Thruster}) = get_wr_b(thr.propeller) #only external
get_hr_b(thr::System{<:Thruster}) = get_hr_b(thr.propeller)
end #module