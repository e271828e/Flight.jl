module TestLandingGear

using Test
using LinearAlgebra
using BenchmarkTools
using StaticArrays

using Flight.FlightCore
using Flight.FlightPhysics

using Flight.FlightComponents.LandingGear
using Flight.FlightComponents.LandingGear: Rolling, Skidding, FrictionCoefficients, get_μ

export test_landing_gear

function test_landing_gear()
    @testset verbose = true "LandingGear" begin
        test_steering()
        test_braking()
        test_simple_damper()
        test_strut()
        test_landing_gear_unit()
    end
end

function test_braking()

    @testset verbose = true "Braking" begin

        nb_sys = System(NoBraking())
        @test isnothing(nb_sys.x)
        @test isnothing(nb_sys.u)
        @test isnothing(nb_sys.y)
        @test LandingGear.get_braking_factor(nb_sys) == 0
        @test @ballocated(f_ode!($nb_sys)) == 0
        @test @ballocated(f_step!($nb_sys)) == 0

        db_sys = System(DirectBraking(η_br = 0.8))
        db_sys.u[] = -2
        @test db_sys.u[] == 0
        db_sys.u[] = 1.3
        @test db_sys.u[] == 1
        db_sys.u[] = 0.5
        @test db_sys.u[] == 0.5
        f_ode!(db_sys)
        @test LandingGear.get_braking_factor(db_sys) == 0.4
        @test @ballocated(f_ode!($db_sys)) == 0
        @test @ballocated(f_step!($db_sys)) == 0
    end

end

function test_steering()

    @testset verbose = true "Steering" begin

        ns_sys = System(NoSteering())
        @test isnothing(ns_sys.x)
        @test isnothing(ns_sys.u)
        @test isnothing(ns_sys.y)
        @test LandingGear.get_steering_angle(ns_sys) == 0
        @test @ballocated(f_ode!($ns_sys)) == 0
        @test @ballocated(f_step!($ns_sys)) == 0

        ds_sys = System(DirectSteering(ψ_max = π/8))
        ds_sys.u[] = -2
        @test ds_sys.u[] == -1
        ds_sys.u[] = 1.3
        @test ds_sys.u[] == 1
        ds_sys.u[] = 0.5
        @test ds_sys.u[] == 0.5
        f_ode!(ds_sys)
        @test LandingGear.get_steering_angle(ds_sys) == π/8 * 0.5
        @test @ballocated(f_ode!($ds_sys)) == 0
        @test @ballocated(f_step!($ds_sys)) == 0
    end

end

function test_simple_damper()

    @testset verbose = true "SimpleDamper" begin
        damper = LandingGear.SimpleDamper()
        @test LandingGear.get_force(damper, -0.1, 0) > 0
        @test LandingGear.get_force(damper, 0, -1) > 0
        # @test_throws AssertionError LandingGear.get_force(damper, -5, 0)
    end

end

function test_strut()

    @testset verbose = true "Strut" begin

        #friction parameters
        @test get_μ(FrictionCoefficients(Rolling(), DryTarmac), 0.0075) ≈ 0.025
        @test get_μ(FrictionCoefficients(Skidding(), DryTarmac), 0.0075) ≈ 0.5
        @test get_μ(FrictionCoefficients(Skidding(), WetTarmac), 1e-5) ≈ 0.25
        @test get_μ(FrictionCoefficients(Skidding(), IcyTarmac), 10) ≈ 0.025

        damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)
        strut = Strut(l_0 = 1.0, damper = damper) |> System

        @test length(strut.x) == 2
        @test length(strut.frc.u.reset) == 2

        steering = System(DirectSteering(ψ_max = π/6))
        braking = System(DirectBraking())
        terrain = HorizontalTerrain()
        loc = NVector()

        #set the initial 2D Location
        h_trn = TerrainData(terrain, loc).altitude
        h = h_trn + 0.9

        #wow = false
        kin = KinematicInit(; h = h_trn + 1.1) |> KinematicData
        f_ode!(strut, steering, braking, terrain, kin) #update Strut
        @test strut.y.wow === false
        #when f_step! executes after the simulation step, the friction
        #compensator reset input will be set
        @test f_step!(strut) == false
        @test strut.frc.u.reset == [true, true]
        #but it will not take effect until the next call to f_ode!
        @test strut.y.frc.reset == [false, false]
        f_ode!(strut, steering, braking, terrain, kin)
        @test strut.y.frc.reset == [true, true]
        strut.x .= 1
        @test all(strut.x .== 1)
        strut.x .= 0
        @test all(strut.x .== 0) #state has not been modified yet
        @test f_step!(strut) == false #x was already 0, not modified
        strut.x[1] = 1
        @test f_step!(strut) == true #now it has
        @test all(strut.x .== 0)

        #normal static load
        kin = KinematicInit(; h ) |> KinematicData
        f_ode!(strut, steering, braking, terrain, kin)
        @test strut.y.ξ ≈ -0.1 #strut is compressed
        @test strut.y.F ≈ 2500 #pure elastic force, positive along the strut's z axis
        @test strut.y.μ_eff == [0, 0] #no motion, no effective friction
        @test strut.y.f_c[1:2] == [0, 0] #no effective friction, no tangential force
        @test strut.y.f_c[3] < 0 #ground reaction negative along contact frame z-axis

        #reset input was set in the previous wow==false test, it takes a call to
        #f_step! to set it back to false
        f_step!(strut)
        @test strut.frc.u.reset == [false, false]

        #oblique static load
        kin = KinematicInit(; h, q_nb = REuler(0, 0, π/12)) |> KinematicData
        f_ode!(strut, steering, braking, terrain, kin)
        @test strut.y.ξ > -0.1 #strut is less compressed
        @test strut.y.F < 2500 #force is smaller

        #axial compressing load
        kin = KinematicInit(; h, v_eOb_n = [0,0,1]) |> KinematicData
        f_ode!(strut, steering, braking, terrain, kin)
        @test strut.y.ξ_dot < 0 #strut is compressing
        @test strut.y.F ≈ 3500 #damping force is added to elastic force

        #low positive longitudinal velocity
        kin = KinematicInit(; h, v_eOb_n = [1e-4, 0, 0] ) |> KinematicData
        f_ode!(strut, steering, braking, terrain, kin)
        @test isapprox.(strut.y.v_eOc_c, [1e-4, 0], atol = 1e-6) |> all
        @test strut.y.μ_max[1] <= strut.y.μ_roll
        @test strut.y.μ_eff[2] == 0 #no lateral motion, no lateral effective friction
        #longitudinal effective friction should be small and negative
        @test strut.y.μ_eff[1] < 0 && abs(strut.y.μ_eff[1]) < strut.y.μ_max[1]
        @test strut.ẋ[1] < 0 #longitudinal velocity error integral should be increasing

        #low positive lateral velocity
        kin = KinematicInit(; h, q_nb = REuler(φ = 0), v_eOb_n = [0, -1e-4, 0] ) |> KinematicData
        f_ode!(strut, steering, braking, terrain, kin)
        @test strut.y.wr_b.F[2] > 0
        @test strut.ẋ[2] > 0 #lateral velocity error integral should be increasing

        #large positive lateral velocity
        kin = KinematicInit(; h, q_nb = REuler(φ = 0), v_eOb_n = [0, -1, 0] ) |> KinematicData
        f_ode!(strut, steering, braking, terrain, kin)
        @test strut.frc.y.sat_out[2] == 1 #large velocity saturates
        @test strut.ẋ[2] == 0 #lateral velocity integral should be decreasing

        #advancing motion with compression
        kin = KinematicInit(; h, v_eOb_n = [10,0,1]) |> KinematicData
        f_ode!(strut, steering, braking, terrain, kin)
        @test isapprox.(strut.y.v_eOc_c, [10, 0], atol = 1e-5) |> all

        #lateral motion with compression
        kin = KinematicInit(; h, ω_lb_b = [1,0,0]) |> KinematicData
        f_ode!(strut, steering, braking, terrain, kin)
        @test isapprox.(strut.y.v_eOc_c, [0, -0.9], atol = 1e-5) |> all

        #off-axis load, forward motion
        kin = KinematicInit(; h, q_nb = REuler(φ = π/12), v_eOb_n = [1e-4, 0, 0] ) |> KinematicData
        f_ode!(strut, steering, braking, terrain, kin)
        @test strut.y.wr_b.F[2] < 0

        #steering
        kin = KinematicInit(; h, v_eOb_n = [10,0,1]) |> KinematicData
        steering.u[] = 0.5
        f_ode!(steering)
        f_ode!(strut, steering, braking, terrain, kin)
        @test strut.y.v_eOc_c[1] < 10
        @test strut.y.v_eOc_c[2] < 0

        #braking
        kin = KinematicInit(; h, q_nb = REuler(φ = 0), v_eOb_n = [1e-4, 0, 0] ) |> KinematicData
        braking.u[] = 1
        f_ode!(braking)
        f_ode!(strut, steering, braking, terrain, kin)
        @test strut.y.μ_max[1] > strut.y.μ_roll

        #check for f_ode! allocations with wow = true (wow = false exits
        #prematurely)
        f_ode!(strut, steering, braking, terrain, kin)
        @test strut.y.wow == true
        @test @ballocated(f_ode!($strut, $steering, $braking, $terrain, $kin)) == 0
        @test @ballocated(f_step!($strut)) == 0

    end

end


function test_landing_gear_unit()

    @testset verbose = true "LandingGearUnit" begin

        ldg = System(LandingGearUnit())

        trn = HorizontalTerrain()

        loc = LatLon()
        h_trn = Terrain.TerrainData(trn, loc).altitude

        #wow = true
        kin = KinematicInit(; h = h_trn + 0.9, v_eOb_n = [0,1,0]) |> KinematicData

        @test (@ballocated f_ode!($ldg, $kin, $trn)) == 0
        @test @ballocated(f_step!($ldg)) == 0

    end

end

function test_harness()

    trn = HorizontalTerrain()
    loc = LatLon()
    h_trn = Terrain.TerrainData(trn, loc).altitude

    damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)
    strut = Strut(l_0 = 1.0, damper = damper)
    ldg = LandingGearUnit(; strut) |> System

    #by default LandingGearUnit is initialized with r_ObOs_b = [0,0,0], so
    #h_Os=h_Ob
    h = h_trn + 0.8
    θ = deg2rad(0)
    φ = deg2rad(0)
    q_nb = REuler(; θ, φ)
    v_eOb_n = [0,0,0]
    ω_lb_b = [0,0,0]
    kin = KinematicInit(; h, v_eOb_n, ω_lb_b, q_nb) |> KinematicData
    f_ode!(ldg, kin, trn)
    f_step!(ldg)
    @show ldg.strut.y.wow
    @show ldg.strut.y.wr_b.F
    @show ldg.strut.y.wr_b.M
    return

end


end #module