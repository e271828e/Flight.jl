module TestFriction

using Test
using BenchmarkTools
using StaticArrays
using OrdinaryDiffEq

using Flight
using Flight.Friction: Parameters, Regulator

export test_friction

function test_friction()

    @testset verbose = true "Friction" begin

        @testset verbose = true "Parameters" begin

            μ_s = 0.03; μ_d = 0.02; v_s = 0.005; v_d = 0.01
            fp = Parameters(; μ_s, μ_d, v_s, v_d)
            @test get_μ(fp, 0.001) === μ_s
            @test get_μ(fp, 0.1) === μ_d
            @test get_μ(fp, 0.5(v_s + v_d)) === 0.5(μ_s + μ_d)

        end

        @testset verbose = true "Regulator" begin

            @test Regulator{2}(k_p = [2, 4], k_i = [200.0, 400.0], k_l = [0.1, 0.1]).k_p == [2.0, 4.0]
            @test Regulator{3}(k_p = 3.0, k_i = 10.0, k_l = 0.2).k_i == fill(10.0, 3)
            @test_throws DimensionMismatch Regulator{1}(k_p = [2, 4], k_i = [200.0, 400.0], k_l = zeros(2)) |> System

            reg = Regulator{2}(k_p = [2, 4], k_i = [200.0, 400.0], k_l = [0.1, 0.1])
            sys = reg |> System
            @test sys.u.reset isa MVector{2,Bool}
            @test sys.y.v isa SVector{2,Float64}
            @test sys.y.sat isa SVector{2,Bool}
            @test length(sys.x) === 2

            v = [1, -0.1]
            f_cont!(sys, v)
            @test sys.y.sat == [true, false]
            @test sys.y.α == [-1, 0.4]

            sys.x .= [0.1, 0]
            sys.u.reset .= [false, true]
            @test f_disc!(sys) == false
            @test sys.x == [0.1, 0]

            sys.u.reset .= [true, true]
            @test f_disc!(sys) == true
            @test sys.x == zeros(2)

            @test @ballocated(f_cont!($sys, $v)) == 0
            @test @ballocated(f_disc!($sys)) == 0

        end #testset

    end #testset

end

function friction_regulator_plots()

    reg = Regulator{2}(k_p = [1, 1], k_i = [1, 1], k_l = [0., 0.])
    sys = reg |> System
    v = [0.5, -0.5]

    mdl = Model(sys, (v,); adaptive = false, solver = Heun(), dt = 0.002, y_saveat = 0.002)
    Modeling.step!(mdl, 2, true)
    mdl.u.reset[2] = true
    Modeling.step!(mdl)
    mdl.u.reset[2] = false
    Modeling.step!(mdl, 2, true)

    plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12)
    make_plots(mdl; save_path = joinpath("tmp", "test_friction"), plot_settings...)

end


end #module