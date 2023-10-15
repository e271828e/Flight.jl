module PIDDesign

using StaticArrays
using UnPack
using NLopt
using ControlSystems
using Trapz: trapz


@kwdef struct PIDParams <: FieldVector{3, Float64}
    k_p::Float64 = 1.0
    k_i::Float64 = 1.0
    k_d::Float64 = 0.1
end

T_i(pid::PIDParams) = pid.k_p / pid.k_i
T_d(pid::PIDParams) = pid.k_d / pid.k_p

@kwdef struct OptParams
    τ_f::Float64 = 0.01
    t_sim::Float64 = 5.0
    maxeval::Int64 = 5000
    lower_bounds::PIDParams = PIDParams(; k_p = 0, k_i = 0, k_d = 0)
    upper_bounds::PIDParams = PIDParams(; k_p = 50, k_i = 50, k_d = 10)
    initial_step::PIDParams = PIDParams(; k_p = 0.1, k_i = 0.1, k_d = 0.1)
end

@kwdef struct OptMetrics <: FieldVector{3, Float64}
    Ms::Float64 #maximum sensitivity
    ∫e::Float64 #integrated absolute error
    ef::Float64 #steady-state error
end

function build_PID(pid_params::PIDParams, opt_params::OptParams)
    @unpack k_p, k_i, k_d = pid_params
    τ_f = opt_params.τ_f
    (k_p + k_i * tf(1, [1,0]) + k_d * tf([1, 0], [τ_f, 1])) |> ss
end

function OptMetrics(P::AbstractStateSpace, C::AbstractStateSpace, opt_params::OptParams)

    S = sensitivity(P, C) #sensitivity function
    T = output_comp_sensitivity(P, C) #complementary sensitivity function (AKA closed loop)

    T_step = step(T, opt_params.t_sim)
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

    OptMetrics(; Ms, ∫e, ef)

end

function cost(P::AbstractStateSpace, C::AbstractStateSpace, opt_params::OptParams, weights::OptMetrics)
    costs = OptMetrics(P, C, opt_params)
    return sum(costs .* weights) / sum(weights)
end


function optimize_PID(  P::LTISystem;
                        pid_params::PIDParams = PIDParams(), #initial condition
                        opt_params::OptParams = OptParams(),
                        weights::OptMetrics = OptMetrics(ones(3)))

    x0 = pid_params |> Vector{Float64}
    maxeval = opt_params.maxeval
    lower_bounds = opt_params.lower_bounds |> Vector{Float64}
    upper_bounds = opt_params.upper_bounds |> Vector{Float64}
    initial_step = opt_params.initial_step |> Vector{Float64}

    P = ss(P)
    f_opt = let P = P, opt_params = opt_params, weights = weights
        function (x::Vector{Float64}, ::Vector{Float64})
            C = build_PID(PIDParams(x), opt_params)
            cost(P, C, opt_params, weights)
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

    pid_params_opt = PIDParams(minx)
    pid_opt = build_PID(pid_params_opt, opt_params)
    metrics_opt = OptMetrics(P, pid_opt, opt_params)

    return (exit_flag, minf, metrics_opt, pid_params_opt, pid_opt)


end

function test()


end


end