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

# can't figure out how to register these so they are seen from module methods
# @unit inHg "inHg" InchOfMercury 3386.389*Unitful.Pa false
# @unit hp "hp" Horsepower 735.49875*Unitful.W false

const β = ISA_layers[1].β

inHg2Pa(p) = 3386.389p
ft2m(h) = 0.3048h
hp2W(P) = 735.49875P

function h2δ(h)
    @unpack p, T = ISAData(AltGeop(h))
    p / p_std / √(T / T_std)
end

T_ISA(p) = T_std * (p / p_std) ^ (-β * R / g_std)
p2δ(p) = (p/p_std) * (T_ISA(p)/T_std)^(-0.5)

########################### AbstractPistonEngine ###############################

abstract type AbstractPistonEngine <: SystemDescriptor end

MassTrait(::System{<:AbstractPistonEngine}) = HasNoMass()
WrenchTrait(::System{<:AbstractPistonEngine}) = GetsNoExternalWrench()
AngularMomentumTrait(::System{<:AbstractPistonEngine}) = HasNoAngularMomentum()

############################ PistonEngine ###############################

#represents a family of naturally aspirated, fuel-injected aviation engines.
#based on performance data available for the Lycoming IO360 engine, normalized
#with rated power and RPMs, and extended to a wider RPM range
Base.@kwdef struct PistonEngine <: AbstractPistonEngine
    ω_rated::Float64 = ustrip(u"rad/s", 2700u"rpm") #IO360
    P_rated::Float64 = 200 |> hp2W
end

const GPE_constants = (
    n_idle = 0.25, #normalized idle speed
    n_max = 1.2, #normalized speed (n_max > 1) for which power output drops to zero
    μ_idle_std = 0.3, #normalized idle manifold pressure at idle, standard conditions
    π_idle_std = 0.01, #normalized idle residual power output, standard conditions
)

Base.@kwdef struct PistonEngineY
    throttle::Float64 = 0.0
    n::Float64 = 0.0 #normalized engine speed, n = ω/ω_rated
    δ::Float64 = 0.0 #normalized intake parameter, δ = p/p_std/sqrt(T_ISA/T_std)
    μ::Float64 = 0.0 #normalized manifold pressure, μ = MAP/p_std
    ω::Float64 = 0.0 #angular rate
    M::Float64 = 0.0 #output torque
    P::Float64 = 0.0 #output power
    # ṁ::Float64 = 0.0 #fuel consumption
end

init_x(::Type{PistonEngine}) = ComponentVector(ω = 0.0)
init_u(::Type{PistonEngine}) = Ref(Bounded(0.0, 0, 1))
init_y(::Type{PistonEngine}) = PistonEngineY()

function generate_data(::Type{PistonEngine})

    @unpack n_idle, n_max, μ_idle_std, π_idle_std = GPE_constants

    f_δ_wot = let

        n_range = range(0.667, 1, length = 2)
        μ_range = range(0.401, 0.936, length = 9)

        δ_data = [0.455 0.523 0.587 0.652 0.718 0.781 0.844 0.906 0.965;
                0.464 0.530 0.596 0.662 0.727 0.792 0.855 0.921 0.981]

        extrapolate(scale(interpolate(δ_data, BSpline(Linear())), n_range, μ_range), Line())

    end

    f_μ_wot = let

        n_range = range(0.667, 1, length = 2)
        δ_range = range(0.441, 1, length = 9)

        μ_knots = range(0.401, 0.936, length = 9)
        μ_data = Array{Float64}(undef, length(n_range), length(δ_range))
        for (i,n) in enumerate(n_range)
            #inverse interpolation μ(δ)
            f_μ_1D = LinearInterpolation(f_δ_wot(n, μ_knots), μ_knots, extrapolation_bc = Line())
            μ_data[i, :] = f_μ_1D.(δ_range)
        end

        extrapolate(scale(interpolate(μ_data, BSpline(Linear())), n_range, δ_range), Line())

    end

    f_π_ISA_std = let

        n_range = [n_idle, 0.667, 0.704, 0.741, 0.778, 0.815, 0.852, 0.889, 0.926, 0.963, 1.000, 1.074, n_max]
        μ_range = [μ_idle_std, 0.568, 1.0]

        #except for n_idle or n_max, π values correspond to MAP = 17 inHg and MAP_lim(n)
        μ_knots = [μ_idle_std 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 0.568 μ_idle_std
                  1.000      0.836 0.854 0.874 0.898 0.912 0.939 0.961 0.959 0.958 0.956 0.953 1.000]

        π_knots = [π_idle_std    0.270 0.305 0.335 0.360 0.380 0.405 0.428 0.450 0.476 0.498 0.498 π_idle_std
                  1.8π_idle_std 0.489 0.548 0.609 0.680 0.729 0.810 0.880 0.920 0.965 1.000 1.000 π_idle_std]

        π_data = Array{Float64,2}(undef, (length(n_range), length(μ_range)))

        for i in 1:length(n_range)
            f_π_ISA_std_1D = LinearInterpolation(
                [μ_knots[1,i], μ_knots[2,i]], [π_knots[1,i], π_knots[2,i]], extrapolation_bc = Line())
            π_data[i,:] = f_π_ISA_std_1D.(μ_range)
        end

        LinearInterpolation(
            (n_range, μ_range), π_data, extrapolation_bc = ((Flat(), Flat()), (Flat(), Line())))

    end

   f_π_ISA_wot = let

        n_range = [n_idle, 0.667, 1.000, 1.074, n_max]
        δ_range = [0, 0.441, 1]

        π_data = Array{Float64,2}(undef, length(n_range), length(δ_range))

        π_data[:, 1] .= 0 #power should vanish for δ → 0 ∀n
        π_data[:, 2] .= [π_idle_std, 0.23, 0.409, 0.409, π_idle_std]
        π_data[:, 3] .= [f_π_ISA_std(n, f_μ_wot(n, 1)) for n in n_range]

        LinearInterpolation((n_range,δ_range), π_data, extrapolation_bc = ((Flat(), Flat()), (Flat(), Line())))

   end

   μ_idle_ratio = μ_idle_std / f_μ_wot(n_idle, 1)

   return (f_δ_wot = f_δ_wot, f_μ_wot = f_μ_wot,
           f_π_ISA_std = f_π_ISA_std, f_π_ISA_wot = f_π_ISA_wot,
           μ_idle_ratio = μ_idle_ratio)

end

const GPE_data = generate_data(PistonEngine)


function compute_μ(thr, n, δ)
    @unpack f_μ_wot, μ_idle_ratio = GPE_data
    return f_μ_wot(n, δ) * (μ_idle_ratio + thr * (1 - μ_idle_ratio))
end


function compute_π_ISA(n, μ, δ)

    @unpack f_δ_wot, f_μ_wot, f_π_ISA_std, f_π_ISA_wot, μ_idle_ratio = GPE_data

    #δ at which our μ would be μ_wot
    δ_wot = f_δ_wot(n, μ)
    π_ISA_std = f_π_ISA_std(n, μ)
    π_ISA_wot = f_π_ISA_wot(n, δ_wot)

    if abs(δ_wot - 1) < 1e-3
        π_ISA = π_ISA_std
    else
        π_ISA = π_ISA_std + (π_ISA_wot - π_ISA_std) / (δ_wot - 1) * (δ - 1)
    end

    return π_ISA

end


function f_cont!(sys::System{<:PistonEngine}, air::AirData; T::Real, J::Real)

    @unpack ω_rated, P_rated = sys.params

    @assert J > 0 "Equivalent moment of inertia at the engine shaft must be positive"
    #T_load should be negative under normal operation (for a CCW propeller, this
    #is achieved by means of a gearbox). but in some cases it may become
    #positive, and drive the engine instead of being driven by it

    @show throttle = sys.u[] |> Float64
    @show ω = sys.x.ω

    @show n = ω / ω_rated
    @show δ = p2δ(air.p)
    @show μ = compute_μ(throttle, n, δ)
    @show π_ISA = compute_π_ISA(n, μ, δ)
    @show π_true = π_ISA * √(T_ISA(air.p)/air.T)

    @show P = π_true * P_rated
    @show M = (ω > 0 ? P / ω : 0)

    sys.y = PistonEngineY(; throttle, n, δ, μ, ω, M, P)

end


end #module