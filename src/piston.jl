module Piston

using Interpolations, Plots, Unitful, StructArrays, UnPack

using Flight.Modeling, Flight.Misc
using Flight.Dynamics, Flight.Airdata
using Flight.Atmosphere: ISA_layers, ISAData, p_std, T_std, g_std, R
using Flight.Geodesy: AltGeop

import Flight.Modeling: init_x, init_y, init_u, init_d, f_cont!, f_disc!
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_hr_b, get_wr_b
import Flight.Plotting: plots

export GenericPistonEngine

@unit inHg "inHg" InchOfMercury 3386.389u"Pa" false
@unit hp "hp" Horsepower 735.49875u"W" false
Unitful.register(Piston)

function δ_from_h(h)
    @unpack p, T = ISAData(AltGeop(h))
    p / p_std / √(T / T_std)
end

δ_from_p(p) = (p/p_std) ^ (1+ISA_layers[1].β*R/(2*g_std)) #computes normalized δ

################################################################################
################################ AbstractPistonEngine ########################

abstract type AbstractPistonEngine <: SystemDescriptor end

MassTrait(::System{<:AbstractPistonEngine}) = HasNoMass()
WrenchTrait(::System{<:AbstractPistonEngine}) = GetsNoExternalWrench()
AngularMomentumTrait(::System{<:AbstractPistonEngine}) = HasNoAngularMomentum()

################################################################################
###################### GenericPistonEngine #####################################

Base.@kwdef struct GenericPistonEngine <: AbstractPistonEngine
    ω_rated::Float64 = ustrip(u"rad/s", 2700u"rpm") #IO360
    P_rated::Float64 = ustrip(u"W", 200u"hp")
    n_idle::Float64 = 0.25 #normalized idle speed
    n_max::Float64 = 1.2 #normalized speed above 1 for which power output drops to zero
    μ_idle_std::Float64 = 0.3 #normalized manifold pressure at idle, standard conditions
    π_idle_std::Float64 = 0.01 #normalized residual power output at idle, standard conditions
end


function build_interp_π_ISA_std(eng::GenericPistonEngine)

    @unpack n_idle, n_max, μ_idle_std, π_idle_std = eng

    n_range = [n_idle, 0.667, 0.704, 0.741, 0.778, 0.815, 0.852, 0.889, 0.926, 0.963, 1.0, 1.074, n_max]
    μ_range = [μ_idle_std, 0.568, 1.0]

    μ_std_1 = [μ_idle_std, 0.568, 0.568, 0.568, 0.568, 0.568, 0.568, 0.568, 0.568, 0.568, 0.568, 0.568, μ_idle_std] #each value in n_range_π
    π_ISA_std_1 = [π_idle_std, 0.27, 0.305, 0.335, 0.36, 0.38, 0.405, 0.4275, 0.45, 0.4755, 0.4975, 0.4975, π_idle_std]

    μ_std_2 = [1.0, 0.836, 0.854, 0.874, 0.898, 0.912, 0.939, 0.961, 0.959, 0.958, 0.956, 0.953, 1.0]
    π_ISA_std_2 =  [1.8π_idle_std, 0.489, 0.5475, 0.609, 0.68, 0.72875, 0.81, 0.88, 0.92, 0.965, 1.0, 1.0, π_idle_std]

    #for each normalized engine speed n, build an interpolator π_ISA_std(μ) from
    #the two data points (μ_1(n), π_ISA_std_1(n)), (μ_2(n), π_ISA_std_2(n)).
    #this assumes a linear variation of π_ISA_std(n, μ) with μ, which is
    #apparent from the graph. then evaluate the interpolator at the common set
    #of μ values μ_range.
    π_ISA_std = Array{Float64,2}(undef, (length(n_range), length(μ_range)))
    for i in 1:length(n_range)
        interp_π_ISA_std_1D = LinearInterpolation([μ_std_1[i], μ_std_2[i]], [π_ISA_std_1[i], π_ISA_std_2[i]], extrapolation_bc = Line())
        π_ISA_std[i,:] = interp_π_ISA_std_1D.(μ_range)
    end

    LinearInterpolation((n_range, μ_range), π_ISA_std, extrapolation_bc = ((Flat(), Flat()), (Flat(), Line())))

end

function build_interp_π_ISA_wot(eng::GenericPistonEngine)

    @unpack n_idle, n_max, π_idle_std = eng

    interp_μ_wot = build_interp_μ_wot(eng)
    interp_π_ISA_std = build_interp_π_ISA_std(eng)

    n_range = [n_idle, 0.667, 1.000, 1.074, n_max]
    δ_range = [0, 0.441, 1]

    π_data = Array{Float64,2}(undef, length(n_range), length(δ_range))

    π_data[:, 1] .= 0 #power should vanish for δ → 0 ∀n
    π_data[:, 2] .= [π_idle_std, 0.23, 0.409, 0.409, π_idle_std]
    π_data[:, 3] .= [interp_π_ISA_std(n, interp_μ_wot(n, 1)) for n in n_range]

    LinearInterpolation((n_range,δ_range), π_data, extrapolation_bc = ((Flat(), Flat()), (Flat(), Line())))

end

function build_interp_δ_wot(::GenericPistonEngine)

    n_range = range(0.667, 1, length = 2)
    μ_range = range(0.401, 0.936, length = 9)

    δ_data = [0.455 0.523 0.587 0.652 0.718 0.781 0.844 0.906 0.965;
              0.464 0.530 0.596 0.662 0.727 0.792 0.855 0.921 0.981]

    extrapolate(scale(interpolate(δ_data, BSpline(Linear())), n_range, μ_range), Line())
end

function build_interp_μ_wot(eng::GenericPistonEngine)

    interp_δ_wot = build_interp_δ_wot(eng)

    n_range = range(0.667, 1, length = 2)
    δ_range = range(0.441, 1, length = 9)

    μ_knots = range(0.401, 0.936, length = 9)
    μ_data = Array{Float64}(undef, length(n_range), length(δ_range))
    for (i,n) in enumerate(n_range)
        #inverse interpolation μ(δ)
        interp_μ_1D = LinearInterpolation(interp_δ_wot(n, μ_knots), μ_knots, extrapolation_bc = Line())
        μ_data[i, :] = interp_μ_1D.(δ_range)
    end

    extrapolate(scale(interpolate(μ_data, BSpline(Linear())), n_range, δ_range), Line())

end

# function save_interp_data()
function test_interp()

    @show p_std_inHg = ustrip(u"inHg", p_std*u"Pa")

    eng = GenericPistonEngine()

    interp_π_ISA_std = build_interp_π_ISA_std(eng)
    # @show interp_π_ISA_std(1800/ω_rated, 22/p_std_inHg) * P_rated
    # @show interp_π_ISA_std(2300/ω_rated, 26/p_std_inHg) * P_rated
    # @show interp_π_ISA_std(1300/ω_rated, 28.6/p_std_inHg) * P_rated

    interp_π_ISA_wot = build_interp_π_ISA_wot(eng)
    # @show interp_π_ISA_std(1800/ω_rated, 22/p_std_inHg) * P_rated
    # @show interp_π_ISA_std(2300/ω_rated, 26/p_std_inHg) * P_rated
    # @show interp_π_ISA_std(1300/ω_rated, 28.6/p_std_inHg) * P_rated

    # interp_μ_wot = build_interp_μ_wot()
    # interp_δ_wot = build_interp_δ_wot()

    n_plot = range(0.15, 1.5, length = 100)
    δ_plot = range(1, 0, length = 100)
    μ_plot = range(0.25p_std_inHg, 30, length = 10)/p_std_inHg

    # π_std_plot = [interp_π_ISA_std(n, μ) for (n, μ) in Iterators.product(n_plot, μ_plot)]
    # # plot(μ_plot, π_std_plot')
    # plot(n_plot, π_std_plot)

    π_wot_plot = [interp_π_ISA_wot(n,p) for (n,p) in Iterators.product(n_plot, δ_plot)]
    plot(δ_plot, π_wot_plot')
    # plot(n_plot, π_wot_plot)


end
"""
function f_cont!()

    #ω es propiedad del motor? no, porque no puede calcularla solo con sus
    #datos. para eso necesita el momento de inercia, que no depende solo de el.
    #no, el motor no tiene ω como estado. la cuestion es: el motor deberia ser
    #System o ni siquiera? quiza si, porque lo que si tiene como estado es
    #on/off. a lo mejor debe tener un estado discreto

    #el motor debe ser System

    #engine inputs: throttle
    #engine outputs: throttle, generated_torque, power, fuel consumption, throttle, MAP

    #aunque realmente si podriamos hacer que omega fuera del motor. si le
    #pasamos el torque y el momento de inercia equivalente que esta moviendo. y
    #eso ya si lo calcula su parent System powerplant. llamando a la helice.
    #entonces la ecuacion de momento cinetico axial la puede aplicar el propio
    #motor. tiene esto sentido??

    #en este caso a f_cont! habria que pasarle solo el torque y el momento de
    #inercia en el eje. que se encarga el parent System
    #f_cont!(sys::System{<:GenericPistonEngine}, air::AirData, T_load::Real, Ixx_load::Real)
    #T_load deberia ser negativa en x. O sea, la gearbox se las tiene que apanar
    #para que el motor siempre vea un par resistente CCW, porque un motor es
    #siempre CW.

    #get x, from x omega, omega to n, from n and n_rated, n̄.
    #get air data. from p, delta. from n̄ and delta, M̃_wot

    #compute the wide-open throttle MAP for the given RPMs and altitude
    @show Mbar_wot = interp_Mbar_wot(nbar, δ)
    idle_MAP_ratio = 0.4
    #this can be tuned so that the engine idles at appropriate RPMs with the
    #chosen propeller
    @show Mbar = Mbar_wot * (idle_MAP_ratio + thr * (1 - idle_MAP_ratio))

    #δ at which our Mbar would be Mbar_wot

    Pbar_ISA_std = interp_Pbar_ISA_std(nbar, Mbar)

    @show δ_wot = interp_δ_wot(nbar, Mbar)
    Pbar_ISA_wot = interp_Pbar_ISA_wot(nbar, δ_wot)
    @show P_ISA_std = Pbar_ISA_std * P_rated
    @show P_ISA_wot = Pbar_ISA_wot * P_rated

    #when p_wot is close to p_std with MAP = MAP_wot (thr = 1),
    #p_wot(MAP_wot(p_std)) = p_std and P̃_wot = P̃_std. we need to avoid the
    #division by zero
    @show abs(δ_wot - 1)
    if abs(δ_wot - 1) < 1e-3
        @show "Hi"
        Pbar_ISA = Pbar_ISA_std
    else
        Pbar_ISA = Pbar_ISA_std + (Pbar_ISA_wot - Pbar_ISA_std) / (δ_wot - 1) * (δ - 1)
        @show Pbar_ISA_std
    end

    @show P_ISA = Pbar_ISA * P_rated
    @show MAP = Mbar * MAP_rated
    @show δ = Mbar * MAP_rated

    # P =

    return max(0, P_ISA)

    #we need to return manifold pressure

end

"""
end #module