module TestLinearization

using Test
using BenchmarkTools
using UnPack
using ComponentArrays

using Flight.FlightCore
using Flight.FlightLib

import Flight.FlightLib.Linearization: LinearizedSS, subsystem, delete_vars

export test_linearization

function test_linearization()

    function build_ss(x0, u0, y0)

        ẋ0 = copy(x0)

        A = x0 * x0'
        B = x0 * u0'
        C = y0 * x0'
        D = y0 * u0'

        return LinearizedSS(ẋ0, x0, u0, y0, A, B, C, D)

    end

    function test_system(ss)

        (; ẋ0, x0, u0, y0, A, B, C, D) = ss

        mdl = Model(ss)

        x = 2mdl.x0
        u = 3mdl.u0
        mdl.x .= x
        mdl.u .= u
        f_ode!(mdl)

        @test all(mdl.ẋ .== ẋ0 + A * (x - x0) + B * (u - u0))
        @test all(mdl.y .== y0 + C * (x - x0) + D * (u - u0))

        @test @ballocated(f_ode!($mdl)) === 0
        @test @ballocated(f_step!($mdl)) === 0
        @test @ballocated(f_periodic!(NoScheduling(), $mdl)) === 0

    end

    @testset verbose = true "Linearization" begin

        x0 = ComponentVector(V = 1.0, q = 0.5, θ = 0.3, α = 5.0)
        u0 = ComponentVector(e = 0.1, a = 0.2)
        y0 = ComponentVector(V = 0.3, q = 0.8, θ = 2.0, α = 3.0, f_z = -9.8)

        @testset verbose = true "Model" begin
            build_ss(x0, u0, y0) |> test_system #with ComponentVectors
            build_ss(x0[:], u0[:], y0[:]) |> test_system #with plain Vectors
        end

        @testset verbose = true "Subsystem" begin

            #with ComponentVectors
            lss = build_ss(x0, u0, y0)
            lss_sub = subsystem(lss; x = (:V, :q), u = (:e,), y = (:V, :q, :f_z))
            @test keys(lss_sub.x0) == (:V, :q)
            @test keys(lss_sub.u0) == (:e,)
            @test keys(lss_sub.y0) == (:V, :q, :f_z)
            @test (
                size(lss_sub.A) == (2, 2) &&
                size(lss_sub.B) == (2, 1) &&
                size(lss_sub.C) == (3, 2) &&
                size(lss_sub.D) == (3, 1)
                )

            lss_sub2 = delete_vars(lss, (:θ, :α, :a))
            @test keys(lss_sub2.x0) == (:V, :q)
            @test keys(lss_sub2.u0) == (:e,)
            @test keys(lss_sub2.y0) == (:V, :q, :f_z)
            @test lss_sub.A == lss_sub2.A
            @test lss_sub.B == lss_sub2.B
            @test lss_sub.C == lss_sub2.C
            @test lss_sub.D == lss_sub2.D

            #with plain Vectors
            lss = build_ss(x0[:], u0[:], y0[:])
            lss_sub = subsystem(lss; x = [1, 3], u = [1], y = [1, 3, 5])
            @test (
                length(lss_sub.x0) == 2 &&
                length(lss_sub.u0) == 1 &&
                length(lss_sub.y0) == 3 &&
                size(lss_sub.A) == (2,2) &&
                size(lss_sub.B) == (2, 1) &&
                size(lss_sub.C) == (3, 2) &&
                size(lss_sub.D) == (3, 1))

        end

    end

end

end #module