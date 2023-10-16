module PIDDesign

using StaticArrays
using UnPack
using NLopt
using ControlSystems
using Trapz: trapz

@kwdef struct Params{T}
    k_p::T = 1.0
    k_i::T = 1.0
    k_d::T = 0.1
    τ_f::T = 0.01
end

T_i(pid::Params{Float64}) = pid.k_p / pid.k_i
T_d(pid::Params{Float64}) = pid.k_d / pid.k_p

function Base.Vector(params::Params)
    @unpack k_p, k_i, k_d, τ_f = params
    return [k_p, k_i, k_d, τ_f]
end

@kwdef struct Settings
    t_sim::Float64 = 5.0
    maxeval::Int64 = 5000
    lower_bounds::Params = Params(; k_p = 0.0, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
    upper_bounds::Params = Params(; k_p = 50.0, k_i = 50.0, k_d = 10.0, τ_f = 0.01)
    initial_step::Params = Params(; k_p = 0.1, k_i = 0.1, k_d = 0.1, τ_f = 0.01)
end

@kwdef struct Metrics <: FieldVector{3, Float64}
    Ms::Float64 #maximum sensitivity
    ∫e::Float64 #integrated absolute error
    ef::Float64 #steady-state error
end

function build_PID(params::Params)
    @unpack k_p, k_i, k_d, τ_f = params
    (k_p + k_i * tf(1, [1,0]) + k_d * tf([1, 0], [τ_f, 1])) |> ss
end

function Metrics(P::AbstractStateSpace, C::AbstractStateSpace, settings::Settings)

    S = sensitivity(P, C) #sensitivity function
    T = output_comp_sensitivity(P, C) #complementary sensitivity function (AKA closed loop)

    T_step = step(T, settings.t_sim)
    t = T_step.t
    y = T_step.y |> vec
    e = abs.(y .- 1.0)

    #hinfnorm appears to be quite brittle, so instead we brute force the
    #computation of maximum sensitivity transfer function magnitude
    S_tf = tf(S)
    iω_range = ((10^x)*im for x in range(-3, 3, length=1000))
    S_range = [abs(S_tf.(iω)[1]) for iω in iω_range]
    Ms = maximum(S_range)

    #integrated error should be normalized with respect to the length of the
    #time window
    ∫e = trapz(t, e)/t[end]
    ef = e[end]

    Metrics(; Ms, ∫e, ef)

end

function cost(P::AbstractStateSpace, C::AbstractStateSpace, settings::Settings, weights::Metrics)
    costs = Metrics(P, C, settings)
    return sum(costs .* weights) / sum(weights)
end


function optimize_PID(  P::LTISystem;
                        params::Params = Params(), #initial condition
                        settings::Settings = Settings(),
                        weights::Metrics = Metrics(ones(3)))

    x0 = params |> Vector
    maxeval = settings.maxeval
    lower_bounds = settings.lower_bounds |> Vector
    upper_bounds = settings.upper_bounds |> Vector
    initial_step = settings.initial_step |> Vector

    P = ss(P)
    f_opt = let P = P, settings = settings, weights = weights
        function (x::Vector{Float64}, ::Vector{Float64})
            C = build_PID(Params(x...))
            cost(P, C, settings, weights)
        end
    end

    minx = x0

    #global search
    opt_glb = Opt(:GN_DIRECT_L, length(x0))
    opt_glb.min_objective = f_opt
    opt_glb.maxeval = maxeval
    opt_glb.stopval = 1e-6
    opt_glb.lower_bounds = lower_bounds
    opt_glb.upper_bounds = upper_bounds
    opt_glb.initial_step = initial_step

    (minf, minx, exit_flag) = optimize(opt_glb, x0)

    #second pass with local optimizer
    opt_loc = Opt(:LN_BOBYQA, length(x0))
    opt_loc.min_objective = f_opt
    opt_loc.maxeval = maxeval
    opt_loc.stopval = 1e-6
    opt_loc.lower_bounds = lower_bounds
    opt_loc.upper_bounds = upper_bounds
    opt_loc.initial_step = initial_step

    (minf, minx, exit_flag) = optimize(opt_loc, minx)

    params_opt = Params(minx...)
    pid_opt = build_PID(params_opt)
    metrics_opt = Metrics(P, pid_opt, settings)

    return (exit_flag, minf, metrics_opt, params_opt, pid_opt)


end

function test()


end


end