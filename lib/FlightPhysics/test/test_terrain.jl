module TestTerrain

using Test
using LinearAlgebra
using BenchmarkTools
using StaticArrays

using FlightCore
using FlightPhysics

export test_terrain

function test_terrain()
    @testset verbose = true "Terrain" begin
        test_surface_intersection()
    end
end

function test_surface_intersection()

    @testset verbose = true "SurfaceIntersection" begin

        h_trn = HOrth(500)
        terrain = HorizontalTerrain(; elevation = h_trn) |> Model

        #pick a location with significant geoid undulation
        loc = LatLon(ϕ = deg2rad(40), λ = deg2rad(-3)) |> NVector
        he_trn = HEllip(h_trn, loc)
        r_srf_e = Cartesian(Geographic(loc, he_trn))[:] #terrain surface point
        dn_e = SVector{3,Float64}(-loc[1], -loc[2], -loc[3]) #geodetic down

        H = 2.0 #ray origin height above the terrain surface
        r_eO_e = Cartesian(Geographic(loc, he_trn + H))[:]

        #vertical ray: hits the surface point right below
        srf = SurfaceIntersection(terrain, r_eO_e, dn_e, 3.0)
        @test srf.valid
        @test norm(srf.r_eP_e - r_srf_e) < 1e-6
        @test srf.kt_P_e ≈ dn_e
        @test srf.surface == DryTarmac

        #the hit point's orthometric altitude must recover the terrain
        #elevation, which confirms the geoid offset is handled
        P = srf.r_eP_e |> Cartesian |> Geographic
        @test HOrth(HEllip(P), NVector(P)) ≈ h_trn atol = 1e-6

        #same ray, insufficient l_max: miss
        @test !SurfaceIntersection(terrain, r_eO_e, dn_e, 1.9).valid

        #tilted ray: hits at l = H / cos(θ)
        θ = deg2rad(30)
        u_e = ltf(loc)(SVector{3,Float64}(sin(θ), 0, cos(θ)))
        srf = SurfaceIntersection(terrain, r_eO_e, u_e, 3.0)
        @test srf.valid
        @test (srf.r_eP_e - r_eO_e) ⋅ u_e ≈ H / cos(θ) atol = 1e-6

        #horizontal ray: miss
        u_hrz_e = ltf(loc)(SVector{3,Float64}(1, 0, 0))
        @test !SurfaceIntersection(terrain, r_eO_e, u_hrz_e, 100.0).valid

        #upward ray: miss
        @test !SurfaceIntersection(terrain, r_eO_e, -dn_e, 3.0).valid

        #origin below the surface, ray pointing down: intersection lies behind
        #the origin, so it must be rejected
        r_low_e = Cartesian(Geographic(loc, he_trn - 1.0))[:]
        @test !SurfaceIntersection(terrain, r_low_e, dn_e, 3.0).valid

        #unbounded ray length must not break the query
        @test SurfaceIntersection(terrain, r_eO_e, dn_e, Inf).valid

        #surface type propagation
        terrain.u[] = IcyTarmac
        @test SurfaceIntersection(terrain, r_eO_e, dn_e, 3.0).surface == IcyTarmac
        terrain.u[] = DryTarmac

        @test @ballocated(SurfaceIntersection($terrain, $r_eO_e, $dn_e, 3.0)) == 0

    end

end

end #module
