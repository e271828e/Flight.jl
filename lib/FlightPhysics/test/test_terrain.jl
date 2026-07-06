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
        Δh = 2.0 #ray origin height above the terrain surface

        @testset verbose = true "Orthometric" begin

            h_trn = HOrth(500)
            terrain = UniformTerrain(; elevation = h_trn) |> Model

            P_srf = Cartesian(Geographic(loc, h_trn)) #expected surface point
            O = Cartesian(Geographic(loc, h_trn + Δh))

            #2D query must return the surface point at loc
            data = TerrainData(terrain, loc)
            @test norm(Cartesian(data.P) - P_srf) < 1e-6
            @test data.kt_P_e ≈ dn_e

            #vertical ray: hits the surface point right below
            intersection = SurfaceIntersection(terrain, O, dn_e, 3.0)
            @test intersection.valid
            @test norm(intersection.data.P - P_srf) < 1e-6
            @test intersection.data.kt_P_e ≈ dn_e

            #same ray, insufficient l_max: miss
            @test !SurfaceIntersection(terrain, O, dn_e, 1.9).valid

            #tilted ray: hits at l = Δh / cos(θ)
            θ = deg2rad(30)
            u_e = ltf(loc)(SVector{3,Float64}(sin(θ), 0, cos(θ)))
            intersection = SurfaceIntersection(terrain, O, u_e, 3.0)
            @test intersection.valid
            @test (intersection.data.P - O) ⋅ u_e ≈ Δh / cos(θ) atol = 1e-6

            #horizontal ray: miss
            u_hrz_e = ltf(loc)(SVector{3,Float64}(1, 0, 0))
            @test !SurfaceIntersection(terrain, O, u_hrz_e, 100.0).valid

            #upward ray: miss
            @test !SurfaceIntersection(terrain, O, -dn_e, 3.0).valid

            #origin below the surface, ray pointing down: intersection lies behind
            #the origin, so it must be rejected
            O_low = Cartesian(Geographic(loc, h_trn - Δh))
            @test !SurfaceIntersection(terrain, O_low, dn_e, 3.0).valid

            #unbounded ray length must not break the query
            @test SurfaceIntersection(terrain, O, dn_e, Inf).valid

            #surface type must propagate through both queries
            terrain.u[] = IcyTarmac
            @test TerrainData(terrain, loc).surface == IcyTarmac
            @test SurfaceIntersection(terrain, O, dn_e, 3.0).data.surface == IcyTarmac
            terrain.u[] = DryTarmac

            #zero allocations in both queries
            @test @ballocated(TerrainData($terrain, $loc)) == 0
            @test @ballocated(SurfaceIntersection($terrain, $O, $dn_e, 3.0)) == 0

        end

        @testset verbose = true "Ellipsoidal" begin
            #the datum only affects the Geographic constructor call within each
            #query, so for the ellipsoidal variant we only verify that the
            #elevation flows through unchanged, with no geoid involvement
            h_trn = HEllip(500)
            terrain = UniformTerrain(; elevation = h_trn) |> Model

            P_srf = Cartesian(Geographic(loc, h_trn)) #expected surface point
            O = Cartesian(Geographic(loc, h_trn + Δh))

            @test HEllip(TerrainData(terrain, loc).P) ≈ h_trn

            intersection = SurfaceIntersection(terrain, O, dn_e, 3.0)
            @test intersection.valid
            @test HEllip(intersection.data.P) ≈ h_trn atol = 1e-6

            #zero allocations in both queries
            @test @ballocated(TerrainData($terrain, $loc)) == 0
            @test @ballocated(SurfaceIntersection($terrain, $O, $dn_e, 3.0)) == 0

        end

    end

end


end #module
