module TestLandingGear

using Test, Plots, UnPack, BenchmarkTools, LinearAlgebra, StaticArrays
using OrdinaryDiffEq, SciMLBase

using Flight
using Flight.Friction: get_μ
using Flight.LandingGear: Rolling, Skidding

export test_landing_gear

function test_landing_gear()
    @testset verbose = true "LandingGear" begin
        test_steering()
        test_braking()
        test_simple_damper()
        test_contact()
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
        @test @ballocated(f_cont!($nb_sys)) == 0
        @test @ballocated(f_disc!($nb_sys)) == 0

        db_sys = System(DirectBraking(η_br = 0.8))
        db_sys.u[] = -2
        @test db_sys.u[] == 0
        db_sys.u[] = 1.3
        @test db_sys.u[] == 1
        db_sys.u[] = 0.5
        @test db_sys.u[] == 0.5
        f_cont!(db_sys)
        @test LandingGear.get_braking_factor(db_sys) == 0.4
        @test @ballocated(f_cont!($db_sys)) == 0
        @test @ballocated(f_disc!($db_sys)) == 0
    end

end

function test_steering()

    @testset verbose = true "Steering" begin

        ns_sys = System(NoSteering())
        @test isnothing(ns_sys.x)
        @test isnothing(ns_sys.u)
        @test isnothing(ns_sys.y)
        @test LandingGear.get_steering_angle(ns_sys) == 0
        @test @ballocated(f_cont!($ns_sys)) == 0
        @test @ballocated(f_disc!($ns_sys)) == 0

        ds_sys = System(DirectSteering(ψ_max = π/8))
        ds_sys.u[] = -2
        @test ds_sys.u[] == -1
        ds_sys.u[] = 1.3
        @test ds_sys.u[] == 1
        ds_sys.u[] = 0.5
        @test ds_sys.u[] == 0.5
        f_cont!(ds_sys)
        @test LandingGear.get_steering_angle(ds_sys) == π/8 * 0.5
        @test @ballocated(f_cont!($ds_sys)) == 0
        @test @ballocated(f_disc!($ds_sys)) == 0
    end

end

function test_simple_damper()

    @testset verbose = true "SimpleDamper" begin
        damper = LandingGear.SimpleDamper()
        @test LandingGear.get_force(damper, -0.1, 0) > 0
        @test LandingGear.get_force(damper, 0, -1) > 0
        @test LandingGear.get_force(damper, damper.ξ_min - 0.1, 0) == 0
    end

end

function test_strut()

    @testset verbose = true "Strut" begin

        damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)
        strut = Strut(l_OsP = 1.0, damper = damper) |> System
        y_default = strut.y

        steering = System(DirectSteering(ψ_max = π/6))
        terrain = HorizontalTerrain()

        #set the initial 2D Location
        l2d = LatLon()
        h_trn = TerrainData(terrain, l2d).altitude

        #wow = false
        kin = Kinematics.Initializer( Ob = GeographicLocation(l2d, h_trn + 1.1)) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        @test strut.y == y_default #strut outputs stay at their defaults

        #wow = true
        Ob = GeographicLocation(l2d, h_trn + 0.9)

        #normal static load
        kin = Kinematics.Initializer(; Ob ) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        @test strut.y.ξ ≈ -0.1 #strut is compressed
        @test strut.y.F ≈ 2500 #pure elastic force, positive along the strut's z axis

        #oblique static load
        kin = Kinematics.Initializer(; Ob, q_nb = REuler(0, 0, π/12)) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        @test strut.y.ξ > -0.1 #strut is less compressed
        @test strut.y.F < 2500 #force is smaller

        #compressing load
        kin = Kinematics.Initializer(; Ob, v_eOb_n = [0,0,1]) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        @test strut.y.ξ_dot < 0 #strut is compressing
        @test strut.y.F ≈ 3500 #damping force is added to elastic force

        #advancing motion and compression
        kin = Kinematics.Initializer(; Ob, v_eOb_n = [10,0,1]) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        @test isapprox.(strut.y.v_eOc_c, [10, 0, 0], atol = 1e-5) |> all

        #lateral motion and compression
        kin = Kinematics.Initializer(; Ob, ω_lb_b = [1,0,0]) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        @test isapprox.(strut.y.v_eOc_c, [0, -0.9, 0], atol = 1e-5) |> all

        #set steering input and update steering system
        steering.u[] = 0.5
        f_cont!(steering)

        #check that steering angle is accounted for
        kin = Kinematics.Initializer(; Ob, v_eOb_n = [10,0,1]) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        @test strut.y.v_eOc_c[1] < 10
        @test strut.y.v_eOc_c[2] < 0

        @test @ballocated(f_cont!($strut, $steering, $terrain, $kin)) == 0
        @test @ballocated(f_disc!($strut)) == 0
    end

end


function test_contact()

    @testset verbose = true "Contact" begin

        #friction parameters
        @test get_μ(Friction.Parameters(Rolling(), DryTarmac), 0.0075) ≈ 0.025
        @test get_μ(Friction.Parameters(Skidding(), DryTarmac), 0.0075) ≈ 0.5
        @test get_μ(Friction.Parameters(Skidding(), WetTarmac), 1e-5) ≈ 0.25
        @test get_μ(Friction.Parameters(Skidding(), IcyTarmac), 10) ≈ 0.025

        contact = LandingGear.Contact() |> System
        @test length(contact.x) == 2
        @test length(contact.u.friction.reset) == 2

        strut = Strut(l_OsP = 1.0) |> System
        steering = System(DirectSteering())
        braking = System(DirectBraking())
        terrain = HorizontalTerrain()

        #set the initial 2D Location
        l2d = LatLon()
        h_trn = TerrainData(terrain, l2d).altitude

        #wow = true
        Ob = GeographicLocation(l2d, h_trn + 0.9)

        #normal static load
        kin = Kinematics.Initializer(; Ob ) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        f_cont!(contact, strut, braking)
        @test contact.u.friction.reset == contact.y.friction.reset == [false, false]
        @test contact.y.μ_eff == [0, 0] #no motion, no effective friction
        @test contact.y.f_c[1:2] == [0, 0] #no effective friction, no tangential force
        @test contact.y.f_c[3] < 0 #ground reaction negative along contact frame z-axis

        kin = Kinematics.Initializer(; Ob, v_eOb_n = [1e-4, 0, 0] ) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        f_cont!(contact, strut, braking)
        @test contact.y.μ_max[1] <= contact.y.μ_roll
        @test contact.y.μ_eff[2] == 0 #no lateral motion, no lateral effective friction
        #low positive longitudinal velocity, small negative longitudinal effective friction
        @test contact.y.μ_eff[1] < 0 && abs(contact.y.μ_eff[1]) < contact.y.μ_max[1]
        @test contact.ẋ[1] > 0 #longitudinal velocity integral should be increasing

        #braking
        kin = Kinematics.Initializer(; Ob, q_nb = REuler(φ = 0), v_eOb_n = [1e-4, 0, 0] ) |> Kinematics.Common
        braking.u[] = 1
        f_cont!(braking)
        f_cont!(strut, steering, terrain, kin)
        f_cont!(contact, strut, braking)
        @test contact.y.μ_max[1] > contact.y.μ_roll

        #non-axial load, forward motion
        kin = Kinematics.Initializer(; Ob, q_nb = REuler(φ = π/12), v_eOb_n = [1e-4, 0, 0] ) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        f_cont!(contact, strut, braking)
        @test contact.y.wr_b.F[2] < 0

        #lateral motion, small velocity
        kin = Kinematics.Initializer(; Ob, q_nb = REuler(φ = 0), v_eOb_n = [0, -1e-4, 0] ) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        f_cont!(contact, strut, braking)
        @test contact.y.wr_b.F[2] > 0
        @test contact.ẋ[2] < 0 #lateral velocity integral should be decreasing

        #lateral motion, large velocity
        kin = Kinematics.Initializer(; Ob, q_nb = REuler(φ = 0), v_eOb_n = [0, -1, 0] ) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin)
        f_cont!(contact, strut, braking)
        @test contact.friction.y.sat[2] == true #large velocity saturates
        @test contact.ẋ[2] == 0 #lateral velocity integral should be decreasing

        #with wow = true, f_disc! should not modify x
        contact.x .= 1
        @test f_disc!(contact) == false
        @test all(contact.x .== 1)

        contact.x .= 0
        @test all(contact.x .== 0) #state has not been modified yet

        #check for f_cont! allocations with wow = true (wow = false exits
        #prematurely)
        @test @ballocated(f_cont!($contact, $strut, $braking)) == 0
        @test @ballocated(f_disc!($contact)) == 0

        #wow = false
        kin_data = Kinematics.Initializer( Ob = GeographicLocation(l2d, h_trn + 1.1)) |> Kinematics.Common
        f_cont!(strut, steering, terrain, kin_data) #update Strut
        f_cont!(contact, strut, braking)
        #the friction regulator's reset input is set and propagated to the outputs
        @test contact.u.friction.reset == contact.y.friction.reset == [true, true]
        @test f_disc!(contact) == false #x was already 0, not modified
        contact.x[1] = 1
        @test f_disc!(contact) == true #now it has
        @test all(contact.x .== 0)

    end

end

function test_landing_gear_unit()

    @testset verbose = true "LandingGearUnit" begin

        ldg = System(LandingGearUnit())

        trn = HorizontalTerrain()
        l2d = LatLon()
        h_trn = Terrain.TerrainData(trn, l2d).altitude

        #wow = true
        kin = Kinematics.Initializer( Ob = GeographicLocation(l2d, h_trn + 0.9), v_eOb_n = [0,1,0]) |> Kinematics.Common

        @test (@ballocated f_cont!($ldg, $kin, $trn)) == 0
        @test @ballocated(f_disc!($ldg)) == 0

    end

end

end #module