module TestTerrain

using Test
using LinearAlgebra
using BenchmarkTools
using StaticArrays

using FlightCore
using FlightPhysics

export test_terrain

struct DummyTerrain <: AbstractTerrain end

function test_terrain()
    @testset verbose = true "Terrain" begin
        test_basics()
        test_uniform_terrain()
    end
end

function test_basics()

    @testset verbose = true "Basics" begin

        #TerrainData defaults must place P on the WGS84 ellipsoid at the
        #default 2D location
        data = TerrainData()
        @test data.P == Geographic(NVector(), HOrth(0))
        @test data.kt_P_e ≈ -NVector()[:]
        @test data.surface == DryTarmac
        @test isbits(data)

        @test HOrth(data) == HOrth(0)

        #default SurfaceIntersection must indicate a miss
        srf = SurfaceIntersection()
        @test !srf.valid
        @test isbits(srf)

        #interface methods must throw for subtypes without an implementation
        trn = DummyTerrain() |> Model
        @test_throws MethodError TerrainData(trn, NVector())
        @test_throws MethodError SurfaceIntersection(trn, Geographic(),
                                                    SVector{3,Float64}(0, 0, 1), 1.0)

    end

end

function test_uniform_terrain()

    @testset verbose = true "UniformTerrain" begin

        #pick a 2D location with significant geoid undulation
        loc = LatLon(ϕ = deg2rad(40), λ = deg2rad(-3)) |> NVector
        dn_e = -loc[:] #geodetic down at loc

        h_trn = HOrth(500)
        terrain = UniformTerrain(; elevation = h_trn) |> Model

        test_terrain_queries(terrain, loc)

        #the hit point's orthometric altitude must recover the terrain
        #elevation, which confirms the geoid offset is handled
        r_eO_e = Cartesian(Geographic(loc, HEllip(h_trn, loc) + 2.0))
        P = SurfaceIntersection(terrain, r_eO_e, dn_e, 3.0).data.P |> Geographic
        @test HOrth(P) ≈ h_trn atol = 1e-6

        #surface type must propagate through both queries
        terrain.u[] = IcyTarmac
        @test TerrainData(terrain, loc).surface == IcyTarmac
        @test SurfaceIntersection(terrain, r_eO_e, dn_e, 3.0).data.surface == IcyTarmac
        terrain.u[] = DryTarmac

        #the datum only affects the Geographic constructor call within each
        #query, so for the ellipsoidal variant we only verify that the
        #elevation flows through unchanged, with no geoid involvement
        h_ell = HEllip(500)
        terrain_ell = UniformTerrain(; elevation = h_ell) |> Model

        @test HEllip(TerrainData(terrain_ell, loc).P) ≈ h_ell

        r_eO_e = Cartesian(Geographic(loc, h_ell + 2.0))
        srf = SurfaceIntersection(terrain_ell, r_eO_e, dn_e, 3.0)
        @test srf.valid
        @test HEllip(srf.data.P) ≈ h_ell atol = 1e-6


    end

end

#TerrainData / SurfaceIntersection query battery. expected results are
#re-derived from the terrain's parameters rather than read back from the
#queries under test
function test_terrain_queries(terrain::Model{<:UniformTerrain}, loc::NVector)

    he_trn = HEllip(terrain.elevation, loc) #surface ellipsoidal altitude at loc
    dn_e = -loc[:] #geodetic down at loc
    r_srf_e = Cartesian(Geographic(loc, he_trn)) #expected surface point

    H = 2.0 #ray origin height above the terrain surface
    r_eO_e = Cartesian(Geographic(loc, he_trn + H))

    #2D query must return the surface point at loc
    data = TerrainData(terrain, loc)
    @test norm(Cartesian(data.P) - r_srf_e) < 1e-6
    @test data.kt_P_e ≈ dn_e

    #vertical ray: hits the surface point right below
    srf = SurfaceIntersection(terrain, r_eO_e, dn_e, 3.0)
    @test srf.valid
    @test norm(srf.data.P - r_srf_e) < 1e-6
    @test srf.data.kt_P_e ≈ dn_e

    #both queries must agree on the same surface point
    @test norm(Cartesian(data.P) - srf.data.P) < 1e-6
    @test data.kt_P_e ≈ srf.data.kt_P_e
    @test data.surface == srf.data.surface

    #same ray, insufficient l_max: miss
    @test !SurfaceIntersection(terrain, r_eO_e, dn_e, 1.9).valid

    #tilted ray: hits at l = H / cos(θ)
    θ = deg2rad(30)
    u_e = ltf(loc)(SVector{3,Float64}(sin(θ), 0, cos(θ)))
    srf = SurfaceIntersection(terrain, r_eO_e, u_e, 3.0)
    @test srf.valid
    @test (srf.data.P - r_eO_e) ⋅ u_e ≈ H / cos(θ) atol = 1e-6

    #horizontal ray: miss
    u_hrz_e = ltf(loc)(SVector{3,Float64}(1, 0, 0))
    @test !SurfaceIntersection(terrain, r_eO_e, u_hrz_e, 100.0).valid

    #upward ray: miss
    @test !SurfaceIntersection(terrain, r_eO_e, -dn_e, 3.0).valid

    #origin below the surface, ray pointing down: intersection lies behind
    #the origin, so it must be rejected
    r_low_e = Cartesian(Geographic(loc, he_trn - 1.0))
    @test !SurfaceIntersection(terrain, r_low_e, dn_e, 3.0).valid

    #unbounded ray length must not break the query
    @test SurfaceIntersection(terrain, r_eO_e, dn_e, Inf).valid

    #zero allocations in both queries
    @test @ballocated(TerrainData($terrain, $loc)) == 0
    @test @ballocated(SurfaceIntersection($terrain, $r_eO_e, $dn_e, 3.0)) == 0

end

end #module
