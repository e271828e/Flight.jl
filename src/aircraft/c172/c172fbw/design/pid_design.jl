module PIDDesign

using StaticArrays
using UnPack
using NLopt
using ControlSystems
using Trapz: trapz

@kwdef struct Params{T} <: FieldVector{4, T}
    k_p::T = 1.0
    k_i::T = 0.0
    k_d::T = 0.1
    τ_f::T = 0.01
end

function Base.getproperty(params::Params, s::Symbol)
    if s ∈ fieldnames(Params)
        getfield(params, s)
    elseif s === :T_i
        params.k_p / params.k_i
    elseif s === :T_d
        params.k_d / params.k_p
    else
        error("$(typeof(params)) has no property $s")
    end
end

@kwdef struct Settings
    t_sim::Float64 = 5.0
    maxeval::Int64 = 5000
    lower_bounds::Params = Params(; k_p = 0.0, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
    upper_bounds::Params = Params(; k_p = 50.0, k_i = 50.0, k_d = 10.0, τ_f = 0.01)
    initial_step::Params = Params(; k_p = 0.1, k_i = 0.1, k_d = 0.1, τ_f = 0.01)
end

@kwdef struct Metrics{T} <: FieldVector{3, T}
    Ms::T #maximum sensitivity
    ∫e::T #integrated absolute error
    ef::T #steady-state error
end

@kwdef struct Results
    exit_flag::Symbol
    cost::Float64
    metrics::Metrics{Float64}
    params::Params{Float64}
end

function build_PID(params::Params{<:Real})
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

function cost(P::AbstractStateSpace, C::AbstractStateSpace, settings::Settings, weights::Metrics{<:Real})
    costs = Metrics(P, C, settings)
    return sum(costs .* weights) / sum(weights)
end


function optimize_PID(  P::LTISystem;
                        params_0::Params = Params(), #initial condition
                        settings::Settings = Settings(),
                        weights::Metrics{<:Real} = Metrics(ones(3)),
                        global_search::Bool = true)

    x0 = params_0 |> Vector
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

    if global_search
        opt_glb = Opt(:GN_DIRECT_L, length(x0))
        opt_glb.min_objective = f_opt
        opt_glb.maxeval = maxeval
        opt_glb.stopval = 1e-5
        opt_glb.lower_bounds = lower_bounds
        opt_glb.upper_bounds = upper_bounds
        opt_glb.initial_step = initial_step

        (minf, minx, exit_flag) = optimize(opt_glb, x0)
    end

    #second pass with local optimizer
    opt_loc = Opt(:LN_BOBYQA, length(x0))
    # opt_loc = Opt(:LN_COBYLA, length(x0))
    opt_loc.min_objective = f_opt
    opt_loc.maxeval = 5000
    opt_loc.stopval = 1e-5
    opt_loc.lower_bounds = lower_bounds
    opt_loc.upper_bounds = upper_bounds
    opt_loc.initial_step = initial_step

    (minf, minx, exit_flag) = optimize(opt_loc, minx)

    params_opt = Params(minx...)
    pid_opt = build_PID(params_opt)
    metrics_opt = Metrics(P, pid_opt, settings)

    return Results(exit_flag, minf, metrics_opt, params_opt)


end

function check_results(results::Results,
            thresholds::Metrics{Float64} = PIDDesign.Metrics(; Ms = 1.3, ∫e = 0.05, ef = 0.01))

    @unpack exit_flag, metrics = results

    @assert metrics.Ms < thresholds.Ms #sensitivity function maximum magnitude
    @assert metrics.∫e < thresholds.∫e #normalized absolute integrated error
    @assert metrics.ef < thresholds.ef #remaining error after t_sim

    @assert (exit_flag === :ROUNDOFF_LIMITED) | (exit_flag === :STOPVAL_REACHED)

end



end