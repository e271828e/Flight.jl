module TestLandingGear

using Test
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using UnPack

using Flight.FlightCore
using Flight.FlightLib

#non-exported stuff
using Flight.FlightLib.LandingGear: Rolling, Skidding, FrictionCoefficients, get_μ

export test_landing_gear

function test_landing_gear()
    @testset verbose = true "LandingGear" begin
        test_steering()
        test_braking()
        test_simple_damper()
        test_landing_gear_unit()
    end
end

function test_braking()

    @testset verbose = true "Braking" begin

        nb_sys = System(NoBraking())
        @test isnothing(nb_sys.x)
        @test isnothing(nb_sys.u)
        # @test isnothing(nb_sys.y)
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
        # @test isnothing(ns_sys.y)
        @test LandingGear.get_steering_angle(ns_sys) == 0
        @test @ballocated(f_ode!($ns_sys)) == 0
        @test @ballocated(f_step!($ns_sys)) == 0

        ds_sys = System(DirectSteering(ψ_max = π/8))
        ds_sys.u.engaged = true
        ds_sys.u.input = -2
        @test ds_sys.u.input == -1
        ds_sys.u.input = 1.3
        @test ds_sys.u.input == 1
        ds_sys.u.input = 0.5
        @test ds_sys.u.input == 0.5
        f_ode!(ds_sys)
        @test LandingGear.get_steering_angle(ds_sys, 0.54) == π/8 * 0.5
        ds_sys.u.engaged = false
        f_ode!(ds_sys)
        @test LandingGear.get_steering_angle(ds_sys, 1.0) == 1.0
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

function test_landing_gear_unit()

    @testset verbose = true "Landing Gear Unit" begin

        #friction parameters
        @test get_μ(FrictionCoefficients(Rolling(), DryTarmac), 0.0075) ≈ 0.025
        @test get_μ(FrictionCoefficients(Skidding(), DryTarmac), 0.0075) ≈ 0.5
        @test get_μ(FrictionCoefficients(Skidding(), WetTarmac), 1e-5) ≈ 0.25
        @test get_μ(FrictionCoefficients(Skidding(), IcyTarmac), 10) ≈ 0.025


        ldg = LandingGearUnit(;
            steering = DirectSteering(ψ_max = π/6),
            braking = DirectBraking(),
            strut = Strut(l_0 = 1.0, damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)),
            contact = Contact()) |> System

        @unpack steering, braking, strut, contact = ldg

        @test length(contact.x) == 2
        @test length(contact.frc.u.reset) == 2

        terrain = HorizontalTerrain() |> System
        loc = NVector()

        #set the initial 2D Location
        h_trn = TerrainData(terrain, loc).altitude
        h = h_trn + 0.9

        #wow = false
        kin = KinInit(; h = h_trn + 1.1) |> KinData
        f_ode!(ldg, kin, terrain)

        @test ldg.y.strut.wow === false

        #when f_step! executes after the simulation step, the friction
        #compensator reset input will be set
        f_step!(ldg)
        @test contact.frc.u.reset == [true, true]

        #but it will not take effect until the next call to f_ode!
        @test ldg.y.contact.frc.reset == [false, false]
        f_ode!(ldg, kin, terrain)
        @test ldg.y.contact.frc.reset == [true, true]

        #ensure that the friction regulator does reset
        ldg.x.contact.frc .= 1
        @test all(ldg.x.contact.frc .== 1)
        ldg.x.contact.frc .= 0
        @test all(ldg.x.contact.frc .== 0) #state has not been modified yet

        f_step!(ldg) #x was already 0, not modified
        ldg.x.contact.frc[1] = 1
        f_step!(ldg) == true #now it has been
        @test all(ldg.x.contact.frc .== 0)

        #normal static load
        kin = KinInit(; h ) |> KinData
        f_ode!(ldg, kin, terrain)
        @test ldg.y.strut.ξ ≈ -0.1 #strut is compressed
        @test ldg.y.strut.F_dmp_zs ≈ 2500 #pure elastic force, positive along the strut's z axis
        @test ldg.y.contact.μ_eff == [0, 0] #no motion, no effective friction
        @test ldg.y.contact.f_c[1:2] == [0, 0] #no effective friction, no tangential force
        @test ldg.y.contact.f_c[3] < 0 #ground reaction negative along contact frame z-axis

        #reset input was set in the previous wow==false test, it takes a call to
        #f_step! to set it back to false
        f_step!(ldg)
        @test ldg.contact.frc.u.reset == [false, false]

        #oblique static load
        kin = KinInit(; h, q_nb = REuler(0, 0, π/12)) |> KinData
        f_ode!(ldg, kin, terrain)
        @test ldg.y.strut.ξ > -0.1 #strut is less compressed
        @test ldg.y.strut.F_dmp_zs < 2500 #force is smaller

        #axial compressing load
        kin = KinInit(; h, v_eb_n = [0,0,1]) |> KinData
        f_ode!(ldg, kin, terrain)
        @test ldg.y.strut.ξ_dot < 0 #strut is compressing
        @test ldg.y.strut.F_dmp_zs ≈ 3500 #damping force is added to elastic force

        #low positive longitudinal velocity
        kin = KinInit(; h, v_eb_n = [1e-4, 0, 0] ) |> KinData
        f_ode!(ldg, kin, terrain)
        @test isapprox.(ldg.y.strut.v_ec_xy, [1e-4, 0], atol = 1e-6) |> all
        @test ldg.y.contact.μ_max[1] <= ldg.y.contact.μ_roll
        @test ldg.y.contact.μ_eff[2] == 0 #no lateral motion, no lateral effective friction
        #longitudinal effective friction should be small and negative
        @test ldg.y.contact.μ_eff[1] < 0 && abs(ldg.y.contact.μ_eff[1]) < ldg.y.contact.μ_max[1]
        @test contact.ẋ[1] < 0 #longitudinal velocity error integral should be increasing

        #low positive lateral velocity
        kin = KinInit(; h, q_nb = REuler(φ = 0), v_eb_n = [0, -1e-4, 0] ) |> KinData
        f_ode!(ldg, kin, terrain)
        @test ldg.y.contact.wr_b.F[2] > 0
        @test contact.ẋ[2] > 0 #lateral velocity error integral should be increasing

        #large positive lateral velocity
        kin = KinInit(; h, q_nb = REuler(φ = 0), v_eb_n = [0, -1, 0] ) |> KinData
        f_ode!(ldg, kin, terrain)
        @test contact.frc.y.sat_out[2] == 1 #large velocity saturates
        @test contact.ẋ[2] == 0 #lateral velocity integral should be decreasing

        #advancing motion with compression
        kin = KinInit(; h, v_eb_n = [10,0,1]) |> KinData
        f_ode!(ldg, kin, terrain)
        @test isapprox.(ldg.y.strut.v_ec_xy, [10, 0], atol = 1e-5) |> all

        #lateral motion with compression
        kin = KinInit(; h, ω_wb_b = [1,0,0]) |> KinData
        f_ode!(ldg, kin, terrain)
        @test isapprox.(ldg.y.strut.v_ec_xy, [0, -0.9], atol = 1e-5) |> all

        #off-axis load, forward motion
        kin = KinInit(; h, q_nb = REuler(φ = π/12), v_eb_n = [1e-4, 0, 0] ) |> KinData
        f_ode!(ldg, kin, terrain)
        @test ldg.y.contact.wr_b.F[2] < 0

        #steering
        kin = KinInit(; h, v_eb_n = [10,0,1]) |> KinData
        steering.u.input = 0.5
        f_ode!(ldg, kin, terrain)
        @test ldg.y.strut.v_ec_xy[1] < 10
        @test ldg.y.strut.v_ec_xy[2] < 0

        #braking
        kin = KinInit(; h, q_nb = REuler(φ = 0), v_eb_n = [1e-4, 0, 0] ) |> KinData
        braking.u[] = 1
        f_ode!(ldg, kin, terrain)
        @test ldg.y.contact.μ_max[1] > ldg.y.contact.μ_roll

        #check for f_ode! allocations with wow = true (wow = false does not
        #cover all the code)
        f_ode!(ldg, kin, terrain)
        @test ldg.y.strut.wow == true
        @test @ballocated(f_ode!($ldg, $kin, $terrain)) == 0
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

    #by default LandingGearUnit is initialized with r_bs_b = [0,0,0], so
    #h_s=h_b
    h = h_trn + 0.8
    θ = deg2rad(0)
    φ = deg2rad(0)
    q_nb = REuler(; θ, φ)
    v_eb_n = [0,0,0]
    ω_wb_b = [0,0,0]
    kin = KinInit(; h, v_eb_n, ω_wb_b, q_nb) |> KinData
    f_ode!(ldg, kin, trn)
    f_step!(ldg)
    @show ldg.strut.y.wow
    @show ldg.contact.y.wr_b.F
    @show ldg.contact.y.wr_b.τ
    return

end


end #module