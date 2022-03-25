module Piston

using Interpolations, Unitful, Plots, StructArrays, ComponentArrays, UnPack

using Flight.Modeling, Flight.Misc
using Flight.Dynamics, Flight.Airdata
using Flight.Atmosphere: ISA_layers, ISAData, p_std, T_std, g_std, R
using Flight.Geodesy: AltGeop

import Flight.Modeling: init_x, init_y, init_u, init_d, f_cont!, f_disc!
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_hr_b, get_wr_b
import Flight.Plotting: plots

export PistonEngine

# can't figure out how to register these so that they are seen from module
# methods
# @unit inHg "inHg" InchOfMercury 3386.389*Unitful.Pa false
# @unit hp "hp" Horsepower 735.49875*Unitful.W false

println("Remember to define AngularMomentum trait for PistonPowerplant!")

const β = ISA_layers[1].β

inHg2Pa(p) = 3386.389p
ft2m(h) = 0.3048h
hp2W(P) = 735.49875P

T_ISA(p) = T_std * (p / p_std) ^ (-β * R / g_std)
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
#based on performance data available for the Lycoming IO360 engine, normalized
#with rated power and RPMs, and extended for wider RPM and MAP ranges
struct PistonEngine{D} <: AbstractPistonEngine
    P_rated::Float64
    ω_rated::Float64
    ω_shutdown::Float64 #speed below which engine shuts down
    ω_cutoff::Float64 #speed above ω_rated for which power output drops back to zero
    idle_ratio::Float64 #μ_idle/μ_wot, may be used to adjust idle RPM for a given load
    dataset::D
end

function PistonEngine(;
    P_rated= 200 |> hp2W, ω_rated = ustrip(u"rad/s", 2700u"rpm"), #IO 360
    ω_shutdown = ustrip(u"rad/s", 400u"rpm"), ω_cutoff = ustrip(u"rad/s", 2701*1.4u"rpm"),
    idle_ratio = 0.25)

    n_shutdown = ω_shutdown / ω_rated
    n_cutoff = ω_cutoff / ω_rated
    dataset = generate_dataset(; n_shutdown, n_cutoff)

    PistonEngine{typeof(dataset)}(P_rated, ω_rated, ω_shutdown, ω_cutoff, idle_ratio, dataset)
end

Base.@kwdef mutable struct PistonEngineD
    running::Bool = true
end

Base.@kwdef mutable struct PistonEngineU
    throttle::Bounded{Float64, 0, 1} = 0.0
    start::Bool = false
    stop::Bool = false
end

Base.@kwdef struct PistonEngineY
    throttle::Float64 = 0.0
    n::Float64 = 0.0 #normalized engine speed, n = ω/ω_rated
    δ::Float64 = 0.0 #normalized inlet parameter, δ = p/p_std/sqrt(T_ISA/T_std)
    μ::Float64 = 0.0 #normalized manifold pressure, μ = MAP/p_std
    ω::Float64 = 0.0 #angular rate
    M::Float64 = 0.0 #output torque
    P::Float64 = 0.0 #output power
    # ṁ::Float64 = 0.0 #fuel consumption
end

init_x(c::PistonEngine) = ComponentVector(ω = 1.5 * c.ω_shutdown)
init_d(::PistonEngine) = PistonEngineD()
init_u(::PistonEngine) = PistonEngineU()
init_y(::PistonEngine) = PistonEngineY()

function generate_dataset(; n_shutdown, n_cutoff)

    @assert n_shutdown < 1
    @assert n_cutoff > 1

    δ_wot = let

        n_range = range(0.667, 1, length = 2)
        μ_range = range(0.401, 0.936, length = 9)

        δ_data = [0.455 0.523 0.587 0.652 0.718 0.781 0.844 0.906 0.965;
                0.464 0.530 0.596 0.662 0.727 0.792 0.855 0.921 0.981]

        extrapolate(scale(interpolate(δ_data, BSpline(Linear())), n_range, μ_range), Line())

    end

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

    π_ISA_std = let

        n_range = [n_shutdown, 0.667, 0.704, 0.741, 0.778, 0.815, 0.852, 0.889, 0.926, 0.963, 1.000, 1.074, n_cutoff]
        μ_range = [0, 0.568, 1.0]

        μ_knots = [zeros(1, length(n_range))
                   0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568
                   1.000 0.836 0.854 0.874 0.898 0.912 0.939 0.961 0.959 0.958 0.956 0.953 1.000]

        #power is zero at MAP = 0 regardless of n (first row)
        #power is zero at n_shutdown and n_cutoff regardless of MAP (first and last columns)
        π_knots = [zeros(1, length(n_range))
                   0 0.270 0.305 0.335 0.360 0.380 0.405 0.428 0.450 0.476 0.498 0.498 0
                   0 0.489 0.548 0.609 0.680 0.729 0.810 0.880 0.920 0.965 1.000 0.950 0]

        π_data = Array{Float64,2}(undef, (length(n_range), length(μ_range)))

        for i in 1:length(n_range)
            π_ISA_std_1D = LinearInterpolation(
                μ_knots[:,i], π_knots[:,i], extrapolation_bc = Line())
            π_data[i,:] = π_ISA_std_1D.(μ_range)
        end

        LinearInterpolation(
            (n_range, μ_range), π_data, extrapolation_bc = ((Flat(), Flat()), (Flat(), Flat())))

    end

   π_ISA_wot = let

        n_range = [n_shutdown, 0.667, 1.000, 1.074, n_cutoff]
        δ_range = [0, 0.441, 1]

        π_data = Array{Float64,2}(undef, length(n_range), length(δ_range))

        π_data[:, 1] .= 0 #power should vanish for δ → 0 ∀n
        π_data[:, 2] .= [0, 0.23, 0.409, 0.409, 0] #it also vanishes at n_shutdown and n_cutoff ∀δ
        π_data[:, 3] .= [π_ISA_std(n, μ_wot(n, 1)) for n in n_range] #at δ=1 by definition it's π_ISA_std(μ_wot)

        LinearInterpolation((n_range,δ_range), π_data, extrapolation_bc = ((Flat(), Flat()), (Flat(), Line())))

   end

   return (δ_wot = δ_wot, μ_wot = μ_wot, π_ISA_std = π_ISA_std, π_ISA_wot = π_ISA_wot)

end

function compute_π_ISA(dataset, n, μ, δ)

    δ_wot = dataset.δ_wot(n, μ) #δ at which our μ would be μ_wot
    π_ISA_std = dataset.π_ISA_std(n, μ)
    π_ISA_wot = dataset.π_ISA_wot(n, δ_wot)

    if abs(δ_wot - 1) < 1e-3
        π_ISA = π_ISA_std
    else
        π_ISA = π_ISA_std + (π_ISA_wot - π_ISA_std) / (δ_wot - 1) * (δ - 1)
    end

end

function f_cont!(sys::System{<:PistonEngine}, air::AirData; M_load::Real, J_load::Real)

    @unpack ω_rated, P_rated, dataset, idle_ratio = sys.params

    @assert J_load > 0 "Equivalent moment of inertia at the engine shaft must be positive"
    #M_load should be negative under normal operation (for a CCW propeller, this
    #is achieved by means of a gearbox). however, in negative propeller thrust
    #conditions (generative), it may become positive, and drive the engine
    #instead of being driven by it

    throttle = sys.u.throttle |> Float64
    ω = sys.x.ω
    n = ω / ω_rated
    δ = p2δ(air.p)

    if !sys.d.running
        μ = air.p / p_std #manifold pressure equals ambient pressure
        P = 0.0 #engine power ≈ 0
        M = 0.0 #engine torque ≈ 0
    else
        μ = dataset.μ_wot(n, δ) * (idle_ratio + throttle * (1 - idle_ratio))
        π_ISA = compute_π_ISA(dataset, n, μ, δ)
        π = π_ISA * √(T_ISA(air.p)/air.T)

        #correct for mixture setting

        #for ω < ω_shutdown we should not be here, but just to be safe let's handle
        #division by zero
        P = π * P_rated
        M = (ω > 0 ? P / ω : 0.0)
    end

    sys.ẋ.ω = (M - M_load) / J_load

    sys.y = PistonEngineY(; throttle, n, δ, μ, ω, M, P)

end

function f_disc!(sys::System{<:PistonEngine})

    ω_shutdown = sys.params.ω_shutdown

    x_mod = false

    if !sys.d.running && sys.u.start
        sys.x.ω = 1.5 * ω_shutdown
        x_mod = true
        sys.d.running = true
    end

    (sys.u.stop || sys.x.ω < ω_shutdown) ? sys.d.running = false : nothing

    return x_mod

end

end #module