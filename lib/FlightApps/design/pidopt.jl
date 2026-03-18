module PIDOpt

using StaticArrays, NLopt, ControlSystems
using RobustAndOptimalControl: hinfnorm2
using Trapz: trapz
using ..Control: PIDData, PIDDataPoint

@kwdef struct Metrics{T} <: FieldVector{5, T}
    Ms::T #maximum sensitivity
    ∫e::T #integrated absolute error
    ef::T #absolute steady-state error
    ∫u::T #integrated absolute control effort
    up::T #absolute peak control effort
end

@kwdef struct Settings
    t_sim::Float64 = 5.0
    maxeval::Int64 = 5000
    lower_bounds::PIDDataPoint = PIDDataPoint(; k_p = 0.0, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
    upper_bounds::PIDDataPoint = PIDDataPoint(; k_p = 50.0, k_i = 50.0, k_d = 10.0, τ_f = 0.01)
    initial_step::PIDDataPoint = PIDDataPoint(; k_p = 0.01, k_i = 0.01, k_d = 0.01, τ_f = 0.01)
end

@kwdef struct Results
    exit_flag::Symbol
    cost::Float64
    metrics::Metrics{Float64}
    data::PIDDataPoint
end

function build_PID(data::PIDDataPoint)
    (; k_p, k_i, k_d, τ_f) = data
    (k_p + k_i * tf(1, [1,0]) + k_d * tf([1, 0], [τ_f, 1])) |> ss
end

function Metrics(plant::AbstractStateSpace, pid::AbstractStateSpace,
                       settings::Settings)

    S = sensitivity(plant, pid) #sensitivity function

    #robust computation of H-∞ norm
    Ms, _ = hinfnorm2(S)
    Ms = min(Ms, 1e3) #allow optimizing unstable systems

    T = output_comp_sensitivity(plant, pid) #complementary sensitivity function (AKA closed loop)
    T_step = step(T, settings.t_sim)
    t = T_step.t
    y = T_step.y |> vec
    abs_e = abs.(y .- 1.0)
    #integrated error normalized with respect to the length of the time window
    ∫e = trapz(t, abs_e)/t[end]
    ef = abs_e[end]

    CS = G_CS(plant, pid) #reference to control input
    CS_step = step(CS, settings.t_sim)
    t = CS_step.t
    y = CS_step.y |> vec
    abs_u = abs.(y .- 1.0)
    ∫u = trapz(t, abs_u)/t[end]
    up = maximum(abs_u)

    Metrics(; Ms, ∫e, ef, ∫u, up)

end

function cost(plant::AbstractStateSpace, pid::AbstractStateSpace,
              settings::Settings, weights::Metrics{<:Real})
    costs = Metrics(plant, pid, settings)
    return sum(costs .* weights) / sum(weights)
end


function optimize_PID(  plant::LTISystem;
                    data_0::PIDDataPoint = PIDDataPoint(), #initial condition
                    settings::Settings = Settings(),
                    weights::Metrics{<:Real} = Metrics(ones(5)),
                    global_search::Bool = true)

    x0 = data_0 |> Vector
    lower_bounds = settings.lower_bounds |> Vector
    upper_bounds = settings.upper_bounds |> Vector
    initial_step = settings.initial_step |> Vector
    maxeval = settings.maxeval

    x0 = clamp.(x0, lower_bounds, upper_bounds)

    plant = ss(plant)
    f_opt = let plant = plant, settings = settings, weights = weights
        function (x::Vector{Float64}, ::Vector{Float64})
            pid = build_PID(PIDDataPoint(x...))
            cost(plant, pid, settings, weights)
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

    data_opt = PIDDataPoint(minx...)
    pid_opt = build_PID(data_opt)
    metrics_opt = Metrics(plant, pid_opt, settings)

    return Results(exit_flag, minf, metrics_opt, data_opt)


end

function check_results(results::Results, thresholds::Metrics{Float64})

    (; exit_flag, metrics) = results

    success = true
    if !((exit_flag === :ROUNDOFF_LIMITED) | (exit_flag === :STOPVAL_REACHED))
        @warn "Unexpected exit flag: $exit_flag"
        success = false
    end
    if !(metrics.Ms < thresholds.Ms) #sensitivity function maximum magnitude
        @warn "Sensitivity function maximum magnitude exceeded: $(metrics.Ms) > $(thresholds.Ms)"
        success = false
    end
    if !(metrics.∫e < thresholds.∫e) #normalized absolute integrated error
        @warn "Normalized absolute integrated error exceeded: $(metrics.∫e) > $(thresholds.∫e)"
        success = false
    end
    if !(metrics.ef < thresholds.ef) #remaining error after t_sim
        @warn "Absolute final error exceeded: $(metrics.ef) > $(thresholds.ef)"
        success = false
    end
    return success

end

end #submodule
