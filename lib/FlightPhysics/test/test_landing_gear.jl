module TestLandingGear

using Test
using LinearAlgebra
using BenchmarkTools
using StaticArrays

using FlightCore
using FlightPhysics

using FlightPhysics.LandingGear: Rolling, Skidding, FrictionCoefficients, get_μ

export test_landing_gear

function test_landing_gear()
    @testset verbose = true "LandingGear" begin
        test_steering()
        test_braking()
        test_simple_damper()
        test_landing_gear_unit()
        test_contact_geometry()
    end
end

#terrain surface modelled as a plane through the point at the default 2D
#location and zero orthometric altitude. its inward normal is tilted from the
#local vertical by the slope angle γ about the local north axis, towards the
#east, so the surface rises eastwards
struct SlopedTerrain <: AbstractTerrain
    γ::Float64
end

@no_updates SlopedTerrain

function FlightPhysics.Terrain.SurfaceIntersection(terrain::Model{<:SlopedTerrain},
                            O::Abstract3DPosition,
                            u_e::AbstractVector{<:Real},
                            l_max::Real)

    r_eO_e = Geocentric(O)
    u_e = SVector{3,Float64}(u_e)

    #plane anchor point P0 and inward surface normal
    n_P0_e = NVector()
    r_eP0_e = Geocentric(Geographic(n_P0_e, HOrth(0)))
    kt_P0_e = ltf(n_P0_e)(SVector{3,Float64}(0, sin(terrain.γ), cos(terrain.γ)))

    #ray-plane intersection
    cos_α = kt_P0_e ⋅ u_e
    cos_α > 0 || return SurfaceIntersection() #parallel to or pointing away from the surface
    l = kt_P0_e ⋅ (r_eP0_e - r_eO_e) / cos_α
    0 <= l <= l_max || return SurfaceIntersection() #behind the origin or beyond l_max

    data = TerrainData(P = r_eO_e + l * u_e, kt_P_e = kt_P0_e, surface = DryTarmac)
    SurfaceIntersection(; valid = true, data)

end

function test_braking()

    @testset verbose = true "Braking" begin

        nb_mdl = Model(NoBraking())
        @test isnothing(nb_mdl.x)
        @test isnothing(nb_mdl.u)
        @test LandingGear.get_braking_factor(nb_mdl) == 0
        @test @ballocated(f_ode!($nb_mdl)) == 0
        @test @ballocated(f_step!($nb_mdl)) == 0

        db_mdl = Model(DirectBraking(η_br = 0.8))
        db_mdl.u[] = -2
        @test db_mdl.u[] == 0
        db_mdl.u[] = 1.3
        @test db_mdl.u[] == 1
        db_mdl.u[] = 0.5
        @test db_mdl.u[] == 0.5
        f_ode!(db_mdl)
        @test LandingGear.get_braking_factor(db_mdl) == 0.4
        @test @ballocated(f_ode!($db_mdl)) == 0
        @test @ballocated(f_step!($db_mdl)) == 0
    end

end

function test_steering()

    @testset verbose = true "Steering" begin

        ns_mdl = Model(NoSteering())
        @test isnothing(ns_mdl.x)
        @test isnothing(ns_mdl.u)
        @test LandingGear.get_steering_angle(ns_mdl) == 0
        @test @ballocated(f_ode!($ns_mdl)) == 0
        @test @ballocated(f_step!($ns_mdl)) == 0

        ds_mdl = Model(DirectSteering(ψ_max = π/8))
        ds_mdl.u.engaged = true
        ds_mdl.u.input = -2
        @test ds_mdl.u.input == -1
        ds_mdl.u.input = 1.3
        @test ds_mdl.u.input == 1
        ds_mdl.u.input = 0.5
        @test ds_mdl.u.input == 0.5
        f_ode!(ds_mdl)
        @test LandingGear.get_steering_angle(ds_mdl, 0.54) == π/8 * 0.5
        ds_mdl.u.engaged = false
        f_ode!(ds_mdl)
        @test LandingGear.get_steering_angle(ds_mdl, 1.0) == 1.0
        @test @ballocated(f_ode!($ds_mdl)) == 0
        @test @ballocated(f_step!($ds_mdl)) == 0
    end

end

function test_simple_damper()

    @testset verbose = true "SimpleDamper" begin
        damper = LandingGear.SimpleDamper()
        @test LandingGear.get_force(damper, -0.1, 0) > 0
        @test LandingGear.get_force(damper, 0, -1) > 0
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
            contact = Contact()) |> Model

        (; steering, braking, strut, contact) = ldg

        @test length(contact.x) == 2

        terrain = UniformTerrain() |> Model
        location = NVector()

        #set the initial 2D Location
        h_trn = HOrth(TerrainData(terrain, location))
        h = h_trn + 0.9

        #wow = false
        kin_data = KinInit(; h = h_trn + 1.1) |> KinData
        f_ode!(ldg, terrain, kin_data)

        @test ldg.y.strut.wow === false

        #now, when f_step! executes the friction compensator should be reset
        ldg.x.contact.frc .= 1 #force nonzero value
        @test all(ldg.x.contact.frc .== 1)
        f_step!(ldg)
        @test all(ldg.x.contact.frc .== 0)

        #normal static load
        kin_data = KinInit(; h ) |> KinData
        f_ode!(ldg, terrain, kin_data)
        @test ldg.y.strut.ξ ≈ -0.1 #strut is compressed
        @test ldg.y.strut.F_dmp_zs ≈ 2500 #pure elastic force, positive along the strut's z axis
        @test ldg.y.contact.μ_eff == [0, 0] #no motion, no effective friction
        @test ldg.y.contact.f_c[1:2] == [0, 0] #no effective friction, no tangential force
        @test ldg.y.contact.f_c[3] < 0 #ground reaction negative along contact frame z-axis

        #oblique static load
        kin_data = KinInit(; h, q_nb = REuler(0, 0, π/12)) |> KinData
        f_ode!(ldg, terrain, kin_data)
        @test ldg.y.strut.ξ > -0.1 #strut is less compressed
        @test ldg.y.strut.F_dmp_zs < 2500 #force is smaller

        #axial compressing load
        kin_data = KinInit(; h, v_eb_n = [0,0,1]) |> KinData
        f_ode!(ldg, terrain, kin_data)
        @test ldg.y.strut.ξ_dot < 0 #strut is compressing
        @test ldg.y.strut.F_dmp_zs ≈ 3500 #damping force is added to elastic force

        #low positive longitudinal velocity
        kin_data = KinInit(; h, v_eb_n = [1e-4, 0, 0] ) |> KinData
        f_ode!(ldg, terrain, kin_data)
        @test isapprox.(ldg.y.strut.v_ec_xy, [1e-4, 0], atol = 1e-6) |> all
        @test ldg.y.contact.μ_max[1] <= ldg.y.contact.μ_roll
        @test ldg.y.contact.μ_eff[2] == 0 #no lateral motion, no lateral effective friction
        #longitudinal effective friction should be small and negative
        @test ldg.y.contact.μ_eff[1] < 0 && abs(ldg.y.contact.μ_eff[1]) < ldg.y.contact.μ_max[1]
        @test contact.ẋ[1] < 0 #longitudinal velocity error integral should be increasing

        #low positive lateral velocity
        kin_data = KinInit(; h, q_nb = REuler(φ = 0), v_eb_n = [0, -1e-4, 0] ) |> KinData
        f_ode!(ldg, terrain, kin_data)
        @test ldg.y.contact.wr_b.F[2] > 0
        @test contact.ẋ[2] > 0 #lateral velocity error integral should be increasing

        #large positive lateral velocity
        kin_data = KinInit(; h, q_nb = REuler(φ = 0), v_eb_n = [0, -1, 0] ) |> KinData
        f_ode!(ldg, terrain, kin_data)
        @test contact.frc.y.sat_out[2] == 1 #large velocity saturates
        @test contact.ẋ[2] == 0 #lateral velocity integral should be decreasing

        #advancing motion with compression
        kin_data = KinInit(; h, v_eb_n = [10,0,1]) |> KinData
        f_ode!(ldg, terrain, kin_data)
        @test isapprox.(ldg.y.strut.v_ec_xy, [10, 0], atol = 1e-5) |> all

        #lateral motion with compression
        kin_data = KinInit(; h, ω_wb_b = [1,0,0]) |> KinData
        f_ode!(ldg, terrain, kin_data)
        @test isapprox.(ldg.y.strut.v_ec_xy, [0, -0.9], atol = 1e-5) |> all

        #off-axis load, forward motion
        kin_data = KinInit(; h, q_nb = REuler(φ = π/12), v_eb_n = [1e-4, 0, 0] ) |> KinData
        f_ode!(ldg, terrain, kin_data)
        @test ldg.y.contact.wr_b.F[2] < 0

        #steering
        kin_data = KinInit(; h, v_eb_n = [10,0,1]) |> KinData
        steering.u.input = 0.5
        f_ode!(ldg, terrain, kin_data)
        @test ldg.y.strut.v_ec_xy[1] < 10
        @test ldg.y.strut.v_ec_xy[2] < 0

        #braking
        kin_data = KinInit(; h, q_nb = REuler(φ = 0), v_eb_n = [1e-4, 0, 0] ) |> KinData
        braking.u[] = 1
        f_ode!(ldg, terrain, kin_data)
        @test ldg.y.contact.μ_max[1] > ldg.y.contact.μ_roll

        #check for f_ode! allocations with wow = true (wow = false does not
        #cover all the code)
        f_ode!(ldg, terrain, kin_data)
        @test ldg.y.strut.wow == true
        @test @ballocated(f_ode!($ldg, $terrain, $kin_data)) == 0
        @test @ballocated(f_step!($ldg)) == 0

    end

end

function test_contact_geometry()

    @testset verbose = true "Contact Geometry" begin

        e1 = SVector(1.0, 0, 0)
        e3 = SVector(0.0, 0, 1)

        ldg = LandingGearUnit(;
            strut = Strut(l_0 = 1.0, damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)),
            contact = Contact()) |> Model

        @testset verbose = true "Tilted Strut" begin

            terrain = UniformTerrain() |> Model
            h_trn = HOrth(TerrainData(terrain, NVector()))

            #pitched body, strut origin 0.9 m above the surface
            θ = deg2rad(20)
            kin_data = KinInit(; h = h_trn + 0.9, q_nb = REuler(0, θ, 0)) |> KinData
            f_ode!(ldg, terrain, kin_data)

            (; wow, ξ, α_cs, t_bc) = ldg.y.strut
            @test wow === true
            @test ξ ≈ 0.9/cos(θ) - 1 atol = 1e-9 #the ray is longer than the drop
            @test α_cs ≈ θ atol = 1e-9

            #the contact frame z-axis is the terrain normal, that is, NED down
            @test t_bc.q(e3) ≈ kin_data.q_nb'(e3) atol = 1e-9

            #the contact frame x-axis is orthogonal to it and unit-norm
            @test norm(t_bc.q(e1)) ≈ 1 atol = 1e-9
            @test t_bc.q(e1) ⋅ t_bc.q(e3) ≈ 0 atol = 1e-9

        end

        @testset verbose = true "Sloped Terrain" begin

            γ = deg2rad(15)
            terrain = SlopedTerrain(γ) |> Model

            #level body, strut origin 0.9 m vertically above the plane anchor
            #point, which is therefore the contact point itself
            kin_data = KinInit(; h = HOrth(0.9)) |> KinData
            f_ode!(ldg, terrain, kin_data)

            (; wow, ξ, α_cs, t_bc) = ldg.y.strut
            @test wow === true
            @test ξ ≈ -0.1 atol = 1e-9
            @test α_cs ≈ γ atol = 1e-9

            #body axes are NED axes here, so the contact frame z-axis is the
            #tilted terrain normal
            @test t_bc.q(e3) ≈ SVector(0.0, sin(γ), cos(γ)) atol = 1e-9

            #with the vehicle at rest the ground reaction is purely normal, so
            #it tilts eastwards with the surface
            @test ldg.y.contact.f_c[3] < 0 #reaction along the contact z-axis
            @test ldg.y.contact.F_c[3] < 0 #positive normal force magnitude
            @test ldg.y.contact.wr_b.F[2] < 0

            @test @ballocated(f_ode!($ldg, $terrain, $kin_data)) == 0

        end

        @testset verbose = true "GroundCrash" begin

            terrain = UniformTerrain() |> Model
            h_trn = HOrth(TerrainData(terrain, NVector()))

            #oblique hit: the ray still reaches the surface within l_0, but the
            #strut is nearly tangent to it
            φ = deg2rad(70)
            kin_data = KinInit(; h = h_trn + 0.2, q_nb = REuler(φ = φ)) |> KinData
            f_ode!(ldg, terrain, kin_data)
            @test ldg.y.strut.wow === true
            @test ldg.y.strut.α_cs ≈ φ atol = 1e-9
            @test_throws LandingGear.GroundCrash f_step!(ldg)

            #excessive compression rate
            kin_data = KinInit(; h = h_trn + 0.9, v_eb_n = [0, 0, 11]) |> KinData
            f_ode!(ldg, terrain, kin_data)
            @test_throws LandingGear.GroundCrash f_step!(ldg)

            #inverted attitude: the ray points away from the surface
            kin_data = KinInit(; h = h_trn + 0.5, q_nb = REuler(φ = π)) |> KinData
            f_ode!(ldg, terrain, kin_data)
            @test ldg.y.strut.wow === false
            @test (f_step!(ldg); true) #no contact, no crash

            #excessive damper force: 0.1 m static compression yields 2500 N
            ldg_weak = LandingGearUnit(;
                strut = Strut(l_0 = 1.0, damper = SimpleDamper(k_s = 25000, F_max = 2000))) |> Model
            kin_data = KinInit(; h = h_trn + 0.9) |> KinData
            f_ode!(ldg_weak, terrain, kin_data)
            @test_throws LandingGear.GroundCrash f_step!(ldg_weak)

        end

    end

end


end #module
