module TestPiston

using Test
using Plots
using UnPack
using BenchmarkTools
using LinearAlgebra

using Flight
using Flight.Piston: PistonEngine, inHg2Pa, ft2m, h2δ, p2δ, compute_π_ISA
using Flight.Atmosphere: Atmosphere, p_std, T_std
using Flight.Airdata

export test_piston

function test_piston()
    @testset verbose = true "Piston" begin
        test_dataset()
        test_propagation()
        test_allocations()
    end
end

function test_dataset()
    n_shutdown = 0.15
    n_cutoff = 1.4
    dataset = Piston.generate_dataset(; n_shutdown, n_cutoff)
    ω_rated = 2700
    P_rated = 200

    @testset verbose = true "Dataset" begin

        @testset verbose = true "δ_wot" begin

            let δ_wot = dataset.δ_wot
            #test some values here
            end

        end #testset

        @testset verbose = true "π_ISA_std" begin

            let π_ISA_std = dataset.π_ISA_std
            #test some values here
            end

        end #testset

        @testset verbose = true "π_ISA_wot" begin

            let π_ISA_wot = dataset.π_ISA_wot

                @show π_ISA_wot(1800/ω_rated, 3e3 |> ft2m |> h2δ) * P_rated
                @show π_ISA_wot(2300/ω_rated, 2.4e3 |> ft2m |> h2δ) * P_rated
                @show π_ISA_wot(2500/ω_rated, 10e3 |> ft2m |> h2δ) * P_rated

            end

        end #testset

        @testset verbose = true "π_ISA" begin

            π_ISA = let dataset = dataset
                (n, μ, δ) -> compute_π_ISA(dataset, n, μ, δ)
            end

            #at n_shutdown and below, power is zero regardless of MAP value
            @test π_ISA(n_shutdown, 0, 1) ≈ 0
            @test π_ISA(n_shutdown, dataset.μ_wot(n_shutdown, 1), 1) ≈ 0
            @test π_ISA(0.5*n_shutdown, 0.5, 1) ≈ 0

            #as soon as n rises above n_shutdown, power starts increasing with MAP
            @test π_ISA(1.5*n_shutdown, 0.5, 1) > π_ISA(1.5*n_shutdown, 0.3, 1)

            #sanity checks against IO360 performance charts
            @test 71 <  π_ISA(1800/ω_rated, inHg2Pa(20)/p_std, 3e3 |> ft2m |> h2δ) * P_rated < 84
            @test 131 < π_ISA(2310/ω_rated, inHg2Pa(23.6)/p_std, 2.4e3 |> ft2m |> h2δ) * P_rated < 139
            @test 102 < π_ISA(2500/ω_rated, inHg2Pa(18)/p_std, 10e3 |> ft2m |> h2δ) * P_rated < 119

        end #testset

    end #testset

end #function

function test_propagation()

    @testset verbose = true "Propagation Functions" begin

        kin = KinInit(v_eOb_b = [50, 0, 5]) |> KinData #positive α
        atm = AtmosphereDescriptor() |> System
        atm.u.static.T_sl = T_std + 10
        air = AirData(kin, atm)

        sys = PistonEngine(idle_ratio = 0.2) |> System
        sys.x.ω = 100


        M_load = 0
        J_load = 1
        f_cont!(sys, air; M_load, J_load)

    end #testset

end #function

function test_allocations()

    @test_broken false

end

function plot_dataset()

    dataset = PistonEngine().dataset

    n_plot = range(0, 1.5, length = 100)
    δ_plot = range(1, 0, length = 100)
    @show μ_plot = range(0.1p_std, inHg2Pa(30), length = 10)/p_std

    π_std_plot = [dataset.π_ISA_std(n, μ) for (n, μ) in Iterators.product(n_plot, μ_plot)]
    # plot(μ_plot, π_std_plot')
    plot(n_plot, π_std_plot)

    # π_wot_plot = [dataset.π_ISA_wot(n,p) for (n,p) in Iterators.product(n_plot, δ_plot)]
    # plot(δ_plot, π_wot_plot')
    # # plot(n_plot, π_wot_plot)


end



end #module