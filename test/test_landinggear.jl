module TestLandingGear

using Test, Plots, UnPack, BenchmarkTools, LinearAlgebra

using Flight
using Flight.Terrain
using Flight.LandingGear: μStaticDynamic, Rolling, Skidding, DryTarmac, WetTarmac, IcyTarmac
using Flight.LandingGear: get_μ

function test_landing_gear()
    @testset verbose = true "LandingGearUnit" begin
        test_friction_parameters()
    end
end

function test_default_friction()
    @testset verbose = true "LandingGearFriction" begin
        @test get_μ(DefaultFriction(), Rolling(), DryTarmac()).static == 0.03
        @test get_μ(DefaultFriction(), Rolling(), WetTarmac()).dynamic == 0.02
        @test get_μ(DefaultFriction(), Skidding(), IcyTarmac()) === μStaticDynamic(0.075, 0.025)
    end
end

function test_friction_regulator()
    @testset verbose = true "FrictionRegulator" begin
        μ_sd = get_μ(DefaultFriction(), Skidding(), DryTarmac())
        fr = FrictionRegulator() |> System
        f_cont!(fr, 0.5, μ_sd)
        fr.y |> pwf

    end
end

end