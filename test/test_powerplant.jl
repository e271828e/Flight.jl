# module TestPowerplant

using Test
using DifferentialEquations
using Flight.Powerplant

export test_pwp

function test_pwp(thruster::ElectricThruster = ElectricThruster(), method = Tsit5(); kwargs...)

    x₀ = init_x(thruster)
    u₀ = init_u(thruster)
    params = (u = u₀, d = thruster)

    y₀ = f_output(x₀, u₀, 0, thruster)
    log = SavedValues(Float64, typeof(y₀))
    scb = SavingCallback(f_output, log)

    problem = ODEProblem{true}(f_update!, x₀, (0, Inf), params)
    integrator = init(problem, method; callback = scb, save_everystep = false, kwargs...)
    return (log = log, int = integrator)

end