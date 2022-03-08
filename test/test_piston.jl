module TestPiston

using Test
using Plots
using UnPack
using BenchmarkTools
using LinearAlgebra

using Flight
using Flight.Piston: PistonEngine, GPE_data, GPE_constants, compute_π_ISA, inHg2Pa, ft2m, h2δ, p2δ
using Flight.Atmosphere: p_std
using Flight.Airdata

function test_piston()
    @testset verbose = true "Piston" begin
        test_π_ISA()
    end
end

function test_π_ISA()

    @testset verbose = true "π_ISA" begin

        @unpack f_δ_wot, f_μ_wot, f_π_ISA_std, f_π_ISA_wot = GPE_data
        @unpack n_idle, n_max, μ_idle_std, π_idle_std = GPE_constants

        ω_rated = 2700
        P_rated = 200

        @test compute_π_ISA(n_idle, μ_idle_std, 1) ≈ π_idle_std

        #power output should remain at π_idle_std for n < n_idle
        @test compute_π_ISA(0, μ_idle_std, 1) ≈ compute_π_ISA(n_idle, μ_idle_std, 1)

        #power output should increase with MAP at n_idle (otherwise it would not
        #respond to acceleration)
        @test compute_π_ISA(n_idle, f_μ_wot(n_idle, 1), 0.8) > compute_π_ISA(n_idle, μ_idle_std, 0.8)

        #sanity checks against IO360 performance charts
        @test 71 < compute_π_ISA(1800/ω_rated, inHg2Pa(20)/p_std, 3e3 |> ft2m |> h2δ) * P_rated < 84
        @test 131 < compute_π_ISA(2310/ω_rated, inHg2Pa(23.6)/p_std, 2.4e3 |> ft2m |> h2δ) * P_rated < 139
        @test 102 < compute_π_ISA(2500/ω_rated, inHg2Pa(18)/p_std, 10e3 |> ft2m |> h2δ) * P_rated < 119

    end #testset

end #function

function test_propagation()

    @testset verbose = true "Propagation Functions" begin

        kin = KinInit(v_eOb_b = [50, 0, 5]) |> KinData #positive α
        atm = AtmosphereDescriptor() |> System
        air = AirData(kin, atm)
        ω = 300

        sys = PistonEngine() |> System

        T = 0
        J = 1
        f_cont!(sys, air; T, J)

    end #testset

end #function

function plot_GPE_data(data = GPE_data)

    @unpack f_δ_wot, f_μ_wot, f_π_ISA_std, f_π_ISA_wot = data

    # p_std_inHg = inHg2Pa(p_std)

    n_plot = range(0.15, 1.5, length = 100)
    δ_plot = range(1, 0, length = 100)
    @show μ_plot = range(0.25, inHg2Pa(30), length = 10)/p_std

    π_std_plot = [f_π_ISA_std(n, μ) for (n, μ) in Iterators.product(n_plot, μ_plot)]
    # plot(μ_plot, π_std_plot')
    plot(n_plot, π_std_plot)

    # π_wot_plot = [f_π_ISA_wot(n,p) for (n,p) in Iterators.product(n_plot, δ_plot)]
    # plot(δ_plot, π_wot_plot')
    # # plot(n_plot, π_wot_plot)


end



end #module