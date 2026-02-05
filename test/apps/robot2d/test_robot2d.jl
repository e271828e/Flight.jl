module TestRobot2D

using Test, BenchmarkTools
using Revise

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightApps

includet(joinpath(dirname(@__FILE__), "../../../src/apps/robot2d/robot2d.jl")); using ..Robot2D
using ..Robot2D: Vehicle, Robot, InitParameters, mode_m, mode_v, mode_η

function test_robot2d(; alloc::Bool = true)
    @testset verbose = true "Robot2D" begin
        test_vehicle(; alloc)
        test_controller(; alloc)
    end #testset
end

function test_vehicle(; alloc::Bool = true)

    @testset verbose = true "Vehicle" begin

    mdl = Vehicle() |> Model
    sim = Simulation(mdl; t_end = 20, dt = 0.01)
    (; k_m, b_m, R) = mdl.parameters

    #stationary initialization, vehicle must remain vertical and stationary
    ip = InitParameters()
    init!(sim, ip)
    run!(sim)
    @test all(isapprox.(mdl.x, 0.0; atol = 1e-3))

    #initialization with steady-state velocity, vehicle must remain vertical
    u_m = 0.7
    ip = InitParameters(; u_m)
    init!(sim, ip)
    run!(sim)
    @test isapprox(mdl.x.v, k_m * u_m * R / b_m)
    @test isapprox(mdl.x.ω, 0; atol = 1e-3)
    @test isapprox(mdl.x.θ, 0; atol = 1e-3)
    @test mdl.x.η > 0

    #initialization with zero motor input and forward tilt, must converge at θ = π
    ip = InitParameters(; θ = 0.1)
    init!(sim, ip)
    step!(sim, 0.1, true)
    @test mdl.x.ω > 0 #body must be falling forward
    run!(sim)
    @test isapprox(mdl.x.ω, 0; atol = 1e-3)
    @test isapprox(mdl.x.v, 0; atol = 1e-3)
    @test isapprox(mdl.x.θ, π; atol = 1e-3)
    @test mdl.x.η > 0

    #initialization with zero motor input and backward angular velocity, must
    #converge at θ = -π
    ip = InitParameters(; ω = -0.1)
    init!(sim, ip)
    run!(sim)
    @test isapprox(mdl.x.ω, 0; atol = 1e-3)
    @test isapprox(mdl.x.v, 0; atol = 1e-3)
    @test isapprox(mdl.x.θ, -π; atol = 1e-3)
    @test mdl.x.η < 0

    alloc && @test @ballocated(f_ode!($mdl)) == 0
    alloc && @test @ballocated(f_periodic!(NoScheduling(), $mdl)) == 0
    alloc && @test @ballocated(f_step!($mdl)) == 0

    end #testset

end

function test_controller(; alloc::Bool = true)

    @testset verbose = true "Controller" begin

    # mdl = Robot() |> Model
    mdl = Robot(vehicle = Vehicle(; L = 0.1, R = 0.08, m_1 = 0.5)) |> Model
    sim = Simulation(mdl; t_end = 60, dt = 0.01)
    init!(sim)

    (; vehicle, controller) = mdl.submodels
    (; u) = controller

    u.mode = Robot2D.mode_m
    u.m_ref = 0.1
    step!(sim, 0.1, true)
    @test u.m_ref == vehicle.y.u_m
    @test vehicle.y.θ < 0 #vehicle must be tilting backward
    u.mode = Robot2D.mode_v
    u.v_ref = 0.3
    step!(sim, 10, true)
    @test vehicle.y.v ≈ u.v_ref atol = 1e-3
    u.v_ref = -Inf
    step!(sim, 10, true)
    @test vehicle.y.v ≈ -controller.parameters.v_lim[] atol = 1e-3
    u.mode = Robot2D.mode_η
    u.η_ref = 1.0
    step!(sim, 20, true)
    @test vehicle.y.η ≈ u.η_ref atol = 1e-3

    #test in position mode, wherein both controllers are active
    alloc && @test @ballocated(f_ode!($mdl)) == 0
    alloc && @test @ballocated(f_periodic!(NoScheduling(), $mdl)) == 0
    alloc && @test @ballocated(f_step!($mdl)) == 0

    end #testset

end


function sim_robot(gui::Bool = false)

    mdl = Model(Robot())
    sim = Simulation(mdl; t_end = 100, dt = 0.01)
    init!(sim)
    # mdl.controller.u.mode = Robot2D.mode_η
    # mdl.controller.u.η_ref = 5
    run!(sim; gui)
    ts = TimeSeries(sim)
    return ts

end



end