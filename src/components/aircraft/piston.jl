module Piston

using Interpolations, StaticArrays, StructArrays, ComponentArrays, UnPack

using Flight.Engine.Systems
using Flight.Engine.Utils: Ranged

using Flight.Physics.Kinematics
using Flight.Physics.RigidBody
using Flight.Physics.Geodesy: HGeop

using ..Atmosphere
using ..Atmosphere: ISA_layers, ISAData, p_std, T_std, g_std, R
using ..Propellers: AbstractPropeller, Propeller
using ..Control: PICompensator, PICompensatorU, PICompensatorY


################################################################################
########################### AbstractFuelSupply #################################

abstract type AbstractFuelSupply <: Component end

RigidBody.MassTrait(::System{<:AbstractFuelSupply}) = HasMass()
RigidBody.AngMomTrait(::System{<:AbstractFuelSupply}) = HasNoAngularMomentum()
RigidBody.WrenchTrait(::System{<:AbstractFuelSupply}) = GetsNoExternalWrench()

#to be extended by concrete subtypes
fuel_available(f::System{<:AbstractFuelSupply}) = throw(
    MethodError(fuel_available, (f,)))

#a massless, infinite fuel supply, for testing purposes
struct MagicFuelSupply <: AbstractFuelSupply end

Systems.init(::SystemU, ::MagicFuelSupply) = Ref(true)
Systems.f_ode!(::System{MagicFuelSupply}) = nothing

RigidBody.get_mp_Ob(::System{MagicFuelSupply}) = MassProperties()
fuel_available(f::System{MagicFuelSupply}) = f.u[]


################################################################################
########################### AbstractPistonEngine ###############################

abstract type AbstractPistonEngine <: Component end

const β = ISA_layers[1].β

inHg2Pa(p) = 3386.389p
ft2m(h) = 0.3048h
hp2W(P) = 735.49875P
RPM2radpersec(ω) = ω*π/30

T_ISA(p) = T_std * (p / p_std) ^ (-β * R / g_std)

#inlet parameter from static pressure assuming ISA conditions
p2δ(p) = (p/p_std) * (T_ISA(p)/T_std)^(-0.5)

function h2δ(h)
    @unpack p, T = ISAData(HGeop(h))
    p / p_std / √(T / T_std)
end

#by default, engine mass is assumed to be accounted for in the vehicle's
#structure
RigidBody.MassTrait(::System{<:AbstractPistonEngine}) = HasNoMass()
#assumed to be negligible by default
RigidBody.AngMomTrait(::System{<:AbstractPistonEngine}) = HasNoAngularMomentum()
#no direct external wrench on the engine
RigidBody.WrenchTrait(::System{<:AbstractPistonEngine}) = GetsNoExternalWrench()


################################################################################
################################# Engine #######################################

#represents a family of naturally aspirated, fuel-injected aviation engines.
#based on performance data available for the Lycoming IO-360A engine, normalized
#with rated power and rated speed to allow for arbitrary engine sizing
struct Engine{L} <: AbstractPistonEngine
    P_rated::Float64
    ω_rated::Float64
    ω_stall::Float64 #speed below which engine shuts down
    ω_cutoff::Float64 #speed above ω_rated for which power output drops back to zero
    ω_idle::Float64 #target idle speed
    M_start::Float64 #starter torque
    J::Float64 #equivalent axial moment of inertia of the engine shaft
    idle::PICompensator{1} #idle MAP control compensator
    frc::PICompensator{1} #friction constraint compensator
    lookup::L #performance lookup table
end

function Engine(;
    P_rated= hp2W(200),
    ω_rated = RPM2radpersec(2700), #IO 360
    ω_stall = RPM2radpersec(300),
    ω_cutoff = RPM2radpersec(3100),
    ω_idle = RPM2radpersec(600),
    M_start = 40,
    J = 0.05,
    idle = PICompensator{1}(k_p = 4.0, k_i = 2.0, bounds = (-0.5, 0.5)),
    frc =  PICompensator{1}(k_p = 5.0, k_i = 400.0, k_l = 0.2, bounds = (-1.0, 1.0))
    )

    n_stall = ω_stall / ω_rated
    n_cutoff = ω_cutoff / ω_rated
    lookup = generate_lookup(; n_stall, n_cutoff)

    Engine{typeof(lookup)}(
        P_rated, ω_rated, ω_stall, ω_cutoff, ω_idle, M_start, J, idle, frc, lookup)
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
    idle::PICompensatorU{1} = PICompensatorU{1}()
    frc::PICompensatorU{1} = PICompensatorU{1}()
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
    idle::PICompensatorY{1} = PICompensatorY{1}()
    frc::PICompensatorY{1} = PICompensatorY{1}()
end

Systems.init(::SystemX, eng::Engine) = ComponentVector(
    ω = 0.0, idle = init_x(eng.idle), frc = init_x(eng.frc))
Systems.init(::SystemU, ::Engine) = PistonEngineU()
Systems.init(::SystemY, ::Engine) = PistonEngineY()
Systems.init(::SystemS, ::Engine) = PistonEngineS()

function Systems.f_ode!(eng::System{<:Engine}, air::AirData;
                        M_load::Real, J_load::Real)

    @unpack ω_rated, ω_idle, P_rated, J, M_start, lookup = eng.params
    @unpack frc, idle = eng.subsystems
    @unpack thr, mix, start, stop = eng.u

    throttle = Float64(thr)
    mixture = Float64(mix)
    state = eng.s.state
    ω = eng.x.ω

    #update friction constraint compensator
    frc.u.input .= ω
    frc.u.sat_enable .= true
    f_ode!(frc)

    idle.u.input .= (ω_idle - ω) / ω_idle #normalized ω idle error
    idle.u.sat_enable .= true
    f_ode!(idle) #update idle compensator

    μ_ratio_idle = 0.5 + idle.y.out[1]

    #normalized engine speed
    n = ω / ω_rated

    #inlet air parameter for ISA conditions
    δ = p2δ(air.p)

    #normalized MAP at wide open throttle
    μ_wot = lookup.μ_wot(n, δ)

    #actual, part throttle normalized MAP
    μ = μ_wot * (μ_ratio_idle + throttle * (1 - μ_ratio_idle))

    if state === eng_off

        #with the engine off, introduce a friction constraint to make the
        #propeller actually stop instead of slowing down asymptotically. with
        #the engine running, all friction is assumed to be already accounted for
        #in the performance tables
        M_fr_max = 0.01 * P_rated / ω_rated #1% of rated torque
        M_fr = -frc.y.out[1] .* M_fr_max #scale M_fr_max with compensator feedback

        MAP = air.p
        M_shaft = M_fr
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

function Systems.f_step!(eng::System{<:Engine}, fuel::System{<:AbstractFuelSupply})

    @unpack idle, frc = eng.subsystems
    @unpack ω_stall, ω_idle = eng.params

    ω = eng.x.ω

    if eng.s.state === eng_off

        eng.u.start && (eng.s.state = eng_starting)

    elseif eng.s.state === eng_starting

        !eng.u.start && (eng.s.state = eng_off)
        (ω > ω_idle && fuel_available(fuel)) && (eng.s.state = eng_running)

    else #eng_running

        (eng.u.stop || ω < ω_stall || !fuel_available(fuel)) && (eng.s.state = eng_off)

    end

    x_mod = false
    x_mod |= f_step!(idle)
    x_mod |= f_step!(frc)
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
            μ_1D = linear_interpolation(δ_wot(n, μ_knots), μ_knots, extrapolation_bc = Line())
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
            π_std_1D = linear_interpolation(
                μ_knots[:,i], π_knots[:,i], extrapolation_bc = Line())
            π_data[i,:] = π_std_1D.(μ_data)
        end

        linear_interpolation(
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


################################################################################
############################### Thruster #######################################

#M_shaft is always positive. for a CW thruster, the gear ratio should be
#positive as well and, under normal operating, conditions M_prop will be
#negative. for a CCW thruster, the gear ratio should be negative and M_prop will
#be positive under normal operating conditions.

#the sign of M_prop may be inverted under negative propeller thrust conditions,
#with the propeller driving the engine instead of the other way around

Base.@kwdef struct Thruster{E <: AbstractPistonEngine,
                            P <: AbstractPropeller} <: Component
    engine::E = Engine()
    propeller::P = Propeller()
    gear_ratio::Float64 = 1.0 #gear ratio
    function Thruster(eng::E, prop::P, gear_ratio) where {E, P}
        @assert sign(gear_ratio) * Int(prop.sense) > 0 "Thruster gear ratio sign "*
        "does not match propeller turn sign"
        new{E,P}(eng, prop, gear_ratio)
    end
end


function Systems.f_ode!(thr::System{<:Thruster}, air::AirData, kin::KinematicData)

    @unpack engine, propeller = thr
    @unpack gear_ratio = thr.params

    ω_eng = engine.x.ω
    ω_prop = gear_ratio * ω_eng

    f_ode!(propeller, kin, air, ω_prop)

    M_prop = propeller.y.wr_p.M[1]
    M_eq = gear_ratio * M_prop #load torque seen from the engine shaft

    J_prop = propeller.params.J_xx
    J_eq = gear_ratio^2 * J_prop #load moment of inertia seen from the engine side

    f_ode!(engine, air; M_load = M_eq, J_load = J_eq)

    update_y!(thr)

end

function Systems.f_step!(thr::System{<:Thruster}, fuel::System{<:AbstractFuelSupply})

    @unpack engine, propeller = thr

    x_mod = false
    x_mod |= f_step!(engine, fuel)
    x_mod |= f_step!(propeller)
    return x_mod

end

RigidBody.MassTrait(::System{<:Thruster}) = HasNoMass()
RigidBody.AngMomTrait(::System{<:Thruster}) = HasAngularMomentum()
RigidBody.WrenchTrait(::System{<:Thruster}) = GetsExternalWrench()

RigidBody.get_wr_b(thr::System{<:Thruster}) = get_wr_b(thr.propeller) #only external
RigidBody.get_hr_b(thr::System{<:Thruster}) = get_hr_b(thr.propeller)

end #module