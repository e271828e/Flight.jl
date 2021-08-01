using DifferentialEquations
using Plots

function example1()

    f(u,p,t) = 1.01u
    u0 = 1/2
    tspan = (0.0, 1.0)
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

    plot(sol, linewidth = 5, title = "Solution to the linear ODE with a thick line", xaxis = "Time (t)", yaxis = "u(t) (in μm)", label = "My Thick Line!")
    plot!(sol.t, t->0.5exp(1.01t), lw=3, ls=:dash, label="True solution")

end


function example2()

    u0 = [1.0, 0, 0]
    p = [10.0, 28, 8/3]
    tspan = (0.0, 100.0)

    function lorenz!(du,u,p,t)
        x, y, z = u
        σ, ρ, β = p
        du[1] = dx = σ * (y - x)
        du[2] = dy = x * (ρ - z) - y
        du[3] = dz = x*y - β*z
    end

    problem = ODEProblem(lorenz!, u0, tspan, p)

    sol = solve(problem)

    return sol

end


function example3()

    l = 1.0
    m = 1.0
    g = 9.81

    function pendulum!(du,u,p,t)
        du[1] = u[2] #dθ = ω
        du[2] = -3g/(2l) * sin(u[1]) + 3/(m * l^2) * p(t)
    end

    θ₀ = 0.01
    ω₀ = 0.0
    u₀ = [θ₀, ω₀]
    tspan = (0.0, 10.0)

    M(t) = 0.2sin(t)

    problem = ODEProblem(pendulum!, u₀, tspan, M)

    sol = solve(problem)

    # plot(sol, linewidth=2, xaxis="t", label=["θ (rad)" "ω (rad/s)"], layout=(2,1))

    return prob, sol

end

function example5()

    l = 1.0
    m = 1.0
    g = 9.81

    function pendulum!(du,u,p,t)
        du[1] = u[2] #dθ = ω
        du[2] = -3g/(2l) * sin(u[1]) + 3/(m * l^2) * p(t)
    end

    θ₀ = 0.01
    ω₀ = 0.0
    u₀ = [θ₀, ω₀]
    tspan = (0.0, 10.0)

    M(t) = 0.2sin(t)

    #this adds periodic tstops to which the integrator steps (in addition to the
    #constant, non-adaptive ones), and then executes #the callback function
    function f_periodic(integrator)
        println("Periodic function run at t = $(integrator.t)")
        println("t = $(integrator.t), u = $(integrator.u)")
    end
    pcb = PeriodicCallback(f_periodic, 1, initial_affect = false)


    problem = ODEProblem(pendulum!, u₀, tspan, M, callback = pcb)
    # problem = ODEProblem(pendulum!, u₀, tspan, M)

    sol = solve(problem)

    # plot(sol, linewidth=2, xaxis="t", label=["θ (rad)" "ω (rad/s)"], layout=(2,1))

    return problem, sol

end

mutable struct MyParams
    forcing_function::Function
    twall_last::Union{Float64, Nothing}
end

function example6()

    l = 1.0
    m = 1.0
    g = 9.81

    function pendulum!(du,u,p,t)
        du[1] = u[2] #dθ = ω
        du[2] = -3g/(2l) * sin(u[1]) + 3/(m * l^2) * p.forcing_function(t)
    end

    θ₀ = 0.01
    ω₀ = 0.0
    u₀ = [θ₀, ω₀]
    tspan = (0.0, 10.0)

    params = MyParams(t->0.2sin(t), nothing)

    function f_realtime!(integrator)
        params = integrator.p
        if params.twall_last === nothing
            params.twall_last = time()
        end
        dt_sim_next = get_proposed_dt(integrator)
        twall_next = params.twall_last + dt_sim_next
        twall_now = time()
        sleep_time = twall_next - twall_now
        sleep(sleep_time)
        sleep_actual = time() - twall_now
        excess_sleep = sleep_actual - sleep_time
        println("dt_sim_next = $dt_sim_next, $sleep_time to sleep until next t_stop")
        println("Slept for $sleep_actual, excess: $excess_sleep")
        params.twall_last = twall_next
        u_modified!(integrator, false)
    end

    # pcb = PeriodicCallback(f_realtime!, 1, initial_affect = false)
    cond(u, t, integrator) = true #always call the real time callback
    dcb = DiscreteCallback(cond, f_realtime!)

    problem = ODEProblem(pendulum!, u₀, tspan, params, callback = dcb)
    # integrator = init(problem, Tsit5())
    integrator = init(problem, Heun(), adaptive = false, dt = 0.02)
    # problem = ODEProblem(pendulum!, u₀, tspan, params)

    return integrator



end
