module Piston

using Interpolations, Unitful, Plots, StructArrays, ComponentArrays, UnPack

using Flight.Modeling, Flight.Misc
using Flight.Dynamics, Flight.Airdata
using Flight.Atmosphere: ISA_layers, ISAData, p_std, T_std, g_std, R
using Flight.Geodesy: AltGeop

import Flight.Modeling: init, f_cont!, f_disc!
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_hr_b, get_wr_b
import Flight.Plotting: plots

export PistonEngine

# can't figure out how to register these so that they are seen from module
# methods
# @unit inHg "inHg" InchOfMercury 3386.389*Unitful.Pa false
# @unit hp "hp" Horsepower 735.49875*Unitful.W false

println("Reminder: Define AngularMomentum trait for PistonPowerplant!")

const β = ISA_layers[1].β

inHg2Pa(p) = 3386.389p
ft2m(h) = 0.3048h
hp2W(P) = 735.49875P

T_ISA(p) = T_std * (p / p_std) ^ (-β * R / g_std)

#inlet parameter from static pressure assuming ISA conditions
p2δ(p) = (p/p_std) * (T_ISA(p)/T_std)^(-0.5)

function h2δ(h)
    @unpack p, T = ISAData(AltGeop(h))
    p / p_std / √(T / T_std)
end


########################### AbstractPistonEngine ###############################

abstract type AbstractPistonEngine <: SystemDescriptor end

MassTrait(::System{<:AbstractPistonEngine}) = HasNoMass()
WrenchTrait(::System{<:AbstractPistonEngine}) = GetsNoExternalWrench()
AngularMomentumTrait(::System{<:AbstractPistonEngine}) = HasNoAngularMomentum()

############################ PistonEngine ###############################

#represents a family of naturally aspirated, fuel-injected aviation engines.
#based on performance data available for the Lycoming IO-360A engine. data is
#normalized with rated power and rated speed to allow for arbitrary engine
#sizing
struct PistonEngine{D} <: AbstractPistonEngine
    P_rated::Float64
    ω_rated::Float64
    ω_shutdown::Float64 #speed below which engine shuts down
    ω_cutoff::Float64 #speed above ω_rated for which power output drops back to zero
    μ_ratio_idle::Float64 #μ_idle/μ_wot, may be used to adjust idle RPM for a given load
    dataset::D
end

function PistonEngine(;
    P_rated= 200 |> hp2W, ω_rated = ustrip(u"rad/s", 2700u"rpm"), #IO 360
    ω_shutdown = ustrip(u"rad/s", 400u"rpm"), ω_cutoff = ustrip(u"rad/s", 2701*1.4u"rpm"),
    μ_ratio_idle = 0.25)

    n_shutdown = ω_shutdown / ω_rated
    n_cutoff = ω_cutoff / ω_rated
    dataset = generate_dataset(; n_shutdown, n_cutoff)

    PistonEngine{typeof(dataset)}(P_rated, ω_rated, ω_shutdown, ω_cutoff, μ_ratio_idle, dataset)
end

Base.@kwdef mutable struct PistonEngineD
    running::Bool = true
end

#best economy: mixture around 0.1
#best power: mixture around 0.5
Base.@kwdef mutable struct PistonEngineU
    start::Bool = false
    shutdown::Bool = false
    thr::Bounded{Float64, 0, 1} = 0.0 #throttle setting
    mix::Bounded{Float64, 0, 1} = 0.5 #mixture setting
end

Base.@kwdef struct PistonEngineY
    start::Bool = false #start control
    shutdown::Bool = false #shutdown control
    throttle::Float64 = 0.0 #throttle setting
    mixture::Float64 = 0.0 #mixture setting
    MAP::Float64 = 0.0 #manifold air pressure
    ω::Float64 = 0.0 #angular velocity (crankshaft)
    M::Float64 = 0.0 #output torque
    P::Float64 = 0.0 #output power
    SFC::Float64 = 0.0 #specific fuel consumption
    ṁ::Float64 = 0.0 #fuel consumption
end

init(c::PistonEngine, ::SystemX) = ComponentVector(ω = 1.5 * c.ω_shutdown)
init(::PistonEngine, ::SystemY) = PistonEngineY()
init(::PistonEngine, ::SystemU) = PistonEngineU()
init(::PistonEngine, ::SystemD) = PistonEngineD()

function generate_dataset(; n_shutdown, n_cutoff)

    @assert n_shutdown < 1
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

        n_data = [n_shutdown, 0.667, 0.704, 0.741, 0.778, 0.815, 0.852, 0.889, 0.926, 0.963, 1.000, 1.074, n_cutoff]
        μ_data = [0, 0.568, 1.0]

        μ_knots = [zeros(1, length(n_data))
                   0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568
                   1.000 0.836 0.854 0.874 0.898 0.912 0.939 0.961 0.959 0.958 0.956 0.953 1.000]

        #power is zero at MAP = 0 regardless of n (first row)
        #power is zero at n_shutdown and n_cutoff regardless of MAP (first and last columns)
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

        n_data = [n_shutdown, 0.667, 1.000, 1.074, n_cutoff]
        δ_data = [0, 0.441, 1]

        π_data = Array{Float64,2}(undef, length(n_data), length(δ_data))

        π_data[:, 1] .= 0 #power should vanish for δ → 0 ∀n
        π_data[:, 2] .= [0, 0.23, 0.409, 0.409, 0] #it also vanishes at n_shutdown and n_cutoff ∀δ
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

function compute_π_ISA_pow(dataset, n, μ, δ)

        δ_wot = dataset.δ_wot(n, μ) #δ at which our μ would be μ_wot

        #normalized power at part throttle, sea level, ISA conditions (δ = 1)
        #and maximum power mixture
        π_std = dataset.π_std(n, μ)

        #normalized power at full throttle, altitude, ISA conditions and maximum
        #power mixture
        π_wot = dataset.π_wot(n, δ_wot)

        if abs(δ_wot - 1) < 5e-3
            π_ISA_pow = π_std
        else
            π_ISA_pow = π_std + (π_wot - π_std) / (δ_wot - 1) * (δ - 1)
        end

end


function f_cont!(sys::System{<:PistonEngine}, air::AirData; M_load::Real, J_load::Real)

    @unpack ω_rated, P_rated, dataset, μ_ratio_idle = sys.params

    @assert J_load > 0 "Equivalent moment of inertia at the engine shaft must be positive"
    #M_load should be negative under normal operation (for a CCW propeller, this
    #is achieved by means of a gearbox). however, in negative propeller thrust
    #conditions (generative), it may become positive, and drive the engine
    #instead of being driven by it

    @unpack thr, mix, start, shutdown = sys.u
    throttle = Float64(thr)
    mixture = Float64(mix)
    ω_crankshaft = sys.x.ω

    if sys.d.running

        #normalized engine speed
        n = ω / ω_rated

        #inlet air parameter for ISA conditions
        δ = p2δ(air.p)

        #normalized MAP at wide open throttle
        μ_wot = dataset.μ_wot(n, δ)

        #actual, part throttle normalized MAP
        μ = μ_wot * (μ_ratio_idle + throttle * (1 - μ_ratio_idle))

        #normalized power at part throttle, altitude, ISA conditions and maximum
        #power mixture
        π_ISA_pow = compute_π_ISA_pow(dataset, n, μ, δ)

        #correction for non-ISA conditions
        π_pow = π_ISA_pow * √(T_ISA(air.p)/air.T)

        #correction for arbitrary mixture setting
        π_actual = π_pow * dataset.π_ratio(mixture)

        MAP = μ * p_std
        P = P_rated * π_actual
        M = (ω > 0 ? P / ω : 0.0) #for ω < ω_shutdown we should not even be here
        SFC = dataset.sfc_pow(n, π_actual) * dataset.sfc_ratio(mixture)

    else

        MAP = air.p
        P = 0.0
        M = 0.0
        SFC = 0.0

    end

    ṁ =  SFC * P

    sys.ẋ.ω = (M + M_load) / J_load #M_load must have the appropriate sign!

    sys.y = PistonEngineY(; start, shutdown, throttle, mixture, MAP, ω, M, P, SFC, ṁ)

end

function f_disc!(sys::System{<:PistonEngine}, fuel::Bool)

    ω_shutdown = sys.params.ω_shutdown

    x_mod = false

    if !sys.d.running && sys.u.start && fuel
        sys.x.ω = 1.5 * ω_shutdown
        x_mod = true
        sys.d.running = true
    end

    (sys.x.ω < ω_shutdown || sys.u.shutdown || !fuel) ? sys.d.running = false : nothing

    return x_mod

end

end #module