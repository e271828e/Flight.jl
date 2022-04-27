module TestPiston

using Test
using Plots
using UnPack
using BenchmarkTools
using LinearAlgebra

using Flight
using Flight.Piston: PistonEngine, inHg2Pa, ft2m, h2δ, p2δ, ft2m, compute_π_ISA_pow
using Flight.Piston: eng_off, eng_starting, eng_running
using Flight.Atmosphere: Atmosphere, p_std, T_std
using Flight.Air

export test_piston

function test_piston()
    @testset verbose = true "PistonEngine" begin
        test_engine_dataset()
        test_propagation()
    end
end

function test_engine_dataset()
    n_stall = 0.15
    n_cutoff = 1.4
    dataset = Piston.generate_dataset(; n_stall, n_cutoff)
    ω_rated = 2700
    P_rated = 200

    @testset verbose = true "Dataset" begin

        @testset verbose = true "δ_wot" begin

            let δ_wot = dataset.δ_wot
                #these graphs have been retouched, so allow more leeway here
                @test δ_wot(1800/ω_rated, inHg2Pa(20)/p_std) ≈ (9500 |> ft2m |> h2δ) atol = 0.1
                @test δ_wot(2700/ω_rated, inHg2Pa(22)/p_std) ≈ (7000 |> ft2m |> h2δ) atol = 0.1
                @test δ_wot(2100/ω_rated, inHg2Pa(16)/p_std) ≈ (15250 |> ft2m |> h2δ) atol = 0.1
                @test δ_wot(2300/ω_rated, inHg2Pa(12)/p_std) ≈ (22000 |> ft2m |> h2δ) atol = 0.1
            end

        end #testset

        @testset verbose = true "π_std" begin

            let π_std = dataset.π_std
                @test π_std(1800/ω_rated, inHg2Pa(20)/p_std) * P_rated ≈ 71 atol = 1
                @test π_std(2050/ω_rated, inHg2Pa(24)/p_std) * P_rated ≈ 113 atol = 1
                @test π_std(2400/ω_rated, inHg2Pa(17)/p_std) * P_rated ≈ 85 atol = 1
                @test π_std(2400/ω_rated, inHg2Pa(28.8)/p_std) * P_rated ≈ 176 atol = 1
            #test some values here
            end

        end #testset

        @testset verbose = true "π_wot" begin

            let π_wot = dataset.π_wot
                #these graphs have been retouched, so allow more leeway here
                @test π_wot(1800/ω_rated, 3e3 |> ft2m |> h2δ) * P_rated ≈ 108 atol = 3
                @test π_wot(2300/ω_rated, 2.4e3 |> ft2m |> h2δ) * P_rated ≈ 153 atol = 3
                @test π_wot(2500/ω_rated, 10e3 |> ft2m |> h2δ) * P_rated ≈ 129 atol = 3
                @test π_wot(2000/ω_rated, 20e3 |> ft2m |> h2δ) * P_rated ≈ 65 atol = 3

            end

        end #testset

        @testset verbose = true "π_ISA_pow" begin

            π_ISA_pow = let dataset = dataset
                (n, μ, δ) -> compute_π_ISA_pow(dataset, n, μ, δ)
            end

            #at n_stall and below, power is zero regardless of MAP value
            @test π_ISA_pow(n_stall, 0, 1) ≈ 0
            @test π_ISA_pow(n_stall, dataset.μ_wot(n_stall, 1), 1) ≈ 0
            @test π_ISA_pow(0.5*n_stall, 0.5, 1) ≈ 0

            #as soon as n rises above n_stall, power starts increasing with MAP
            @test π_ISA_pow(1.5*n_stall, 0.5, 1) > π_ISA_pow(1.5*n_stall, 0.3, 1)

            #sanity checks against IO360 performance charts
            @test 71 <  π_ISA_pow(1800/ω_rated, inHg2Pa(20)/p_std, 3e3 |> ft2m |> h2δ) * P_rated < 84
            @test 131 < π_ISA_pow(2310/ω_rated, inHg2Pa(23.6)/p_std, 2.4e3 |> ft2m |> h2δ) * P_rated < 139
            @test 102 < π_ISA_pow(2500/ω_rated, inHg2Pa(18)/p_std, 10e3 |> ft2m |> h2δ) * P_rated < 119

        end #testset

    end #testset

end #function

function test_propagation()

    @testset verbose = true "Propagation" begin

        kin = KinInit(v_eOb_n = [50, 0, 0]) |> KinData
        atm = AtmosphereDescriptor() |> System
        atm.u.static.T_sl = T_std + 10
        air = AirData(kin, atm)
        eng = PistonEngine(μ_ratio_idle = 0.2) |> System
        fuel = System(MagicFuelSupply())

        ω = 100.0
        y_init = eng.y
        f_cont!(eng, air, ω)
        @test eng.y != y_init #y must have been updated

        ω = 0.0
        eng.d.state = eng_off
        f_cont!(eng, air, ω)
        @test eng.y.M == 0

        eng.u.start = true
        f_disc!(eng, fuel, ω)
        @test eng.d.state == eng_starting

        ω = 1.4eng.params.ω_stall
        f_disc!(eng, fuel, ω) #with ω <= 1.5ω_stall, engine won't leave the starting state
        @test eng.d.state == eng_starting
        f_cont!(eng, air, ω)
        @test eng.y.M > 0 #it should output the starter torque

        ω = 1.6eng.params.ω_stall
        f_disc!(eng, fuel, ω) #engine should start now
        @test eng.d.state == eng_running
        f_cont!(eng, air, ω)
        @test eng.y.M > 0 #and generate its own torque

        #commanded shutdown
        eng.d.state = eng_running
        eng.u.shutdown = true
        f_disc!(eng, fuel, ω) #engine should start now
        eng.u.shutdown = false
        @test eng.d.state == eng_off

        #stall shutdown
        eng.d.state = eng_running
        ω = 0.95eng.params.ω_stall
        f_disc!(eng, fuel, ω) #engine should start now
        @test eng.d.state == eng_off
        ω = 1.6eng.params.ω_stall
        eng.d.state = eng_running

        #without fuel, the engine should shut down
        fuel.u[] = false
        f_disc!(eng, fuel, ω)
        @test eng.d.state == eng_off

        #and then fail to start, even above the required speed
        eng.u.start = true
        f_disc!(eng, fuel, ω)
        @test eng.d.state == eng_starting
        f_disc!(eng, fuel, ω)
        @test eng.d.state != eng_running

        #when fuel is available, the engine starts
        fuel.u[] = true
        f_disc!(eng, fuel, ω)
        @test eng.d.state == eng_running

        @test @ballocated(f_cont!($eng, $air, $ω)) == 0
        @test @ballocated(f_disc!($eng, $fuel, $ω)) == 0

    end #testset

end #function


function plot_dataset()

    dataset = PistonEngine().dataset

    n_plot = range(0, 1.5, length = 100)
    δ_plot = range(1, 0, length = 100)
    @show μ_plot = range(0.1p_std, inHg2Pa(30), length = 10)/p_std

    π_std_plot = [dataset.π_std(n, μ) for (n, μ) in Iterators.product(n_plot, μ_plot)]
    # plot(μ_plot, π_std_plot')
    plot(n_plot, π_std_plot)

    # π_wot_plot = [dataset.π_wot(n,p) for (n,p) in Iterators.product(n_plot, δ_plot)]
    # plot(δ_plot, π_wot_plot')
    # # plot(n_plot, π_wot_plot)


end #function

end #module