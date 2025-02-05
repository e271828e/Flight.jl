module TestDynamics

using Test
using LinearAlgebra
using BenchmarkTools
using StaticArrays

using Flight.FlightCore
using Flight.FlightLib

export test_dynamics

struct TestHarness <: SystemDefinition end

function test_dynamics()

    @testset verbose = true "Dynamics" begin

        th = TestHarness() |> System;
        dyn = System(RigidBodyDynamics())
        kin = System(WA())

        kin_init = KinInit(
            loc = LatLon(0, 0),
            h = HOrth(0),
            ω_wb_b = [0.0, 0.0, 0.0],
            q_nb = REuler(0.0, 0.0, 0.0),
            v_eb_n = [0, 0, 0])

        Systems.init!(kin, kin_init)
        kin_data = KinData(kin)

        #set up a rigid body with unit mass and unit inertia tensor at its
        #center of mass G
        rbd_G = RigidBodyDistribution(1.0, diagm(ones(3)))
        mp_G = MassProperties(rbd_G)

        #let its frame origin Ob be located at G, so that its mass
        #properties are the same in both
        mp_b = mp_G


        @testset verbose = true "Essentials" begin


            #apply a force along each axis with zero torque
            wr_ext_b = Wrench(F = [1, 2, 1])
            hr_b = zeros(SVector{3, Float64})
            actions = Actions(th; mp_b, wr_ext_b, hr_b, kin_data)

            #we expect the angular acceleration to be zero, and linear acceleration
            #added to that of free fall
            f_ode!(dyn, mp_b, kin_data, actions)
            @test all(dyn.ẋ.ω_eb_b .≈ 0)
            @test all(dyn.ẋ.v_eb_b .≈ [1, 2, 1 + gravity(kin_data)])

            f_ode!(dyn, mp_b, kin_data, actions)
            @test all(dyn.y.a_eb_b .≈ [1, 2, 1 + gravity(kin_data)])
            @test all(dyn.y.a_ib_b .≈ [1, 2, 1 + G_n(kin_data)[3]])

            #now let G be located 1 meter ahead from Ob along the x axis
            r_bG_b = [1,0,0]
            t_bG = FrameTransform(r = r_bG_b) #b to G
            mp_b = Dynamics.transform(t_bG, mp_G)

            #apply a 1N force along z_b
            wr_ext_b = Wrench(F = [0, 0, 1], M = zeros(3))
            hr_b = zeros(SVector{3, Float64})
            actions = Actions(th; mp_b, wr_ext_b, hr_b, kin_data)

            #we expect a positive unit angular acceleration around y_b, which
            #will also add to the linear acceleration of frame b
            f_ode!(dyn, mp_b, kin_data, actions)
            @test dyn.ẋ.ω_eb_b[2] .≈ 1
            @test dyn.ẋ.v_eb_b[3] ≈ 2 + gravity(kin_data)

            r_bG_e = kin_data.q_eb'(r_bG_b)
            r_eG_e = Cartesian(kin_data.r_eb_e + r_bG_e)

            @test isapprox(dyn.y.a_eb_b[3], 2 + gravity(r_eG_e), atol = 1e-10)
            @test isapprox(dyn.y.a_ib_b[3], 2 + G_n(r_eG_e)[3], atol = 1e-10)
            @test all(isapprox.(dyn.y.f_G_b, [0, 0, 1], atol = 1e-10))

        end

        @testset verbose = true "Performance" begin

            wr_ext_b = Wrench(F = [1, 2, 1])
            hr_b = zeros(SVector{3, Float64})
            actions = Actions(th; mp_b, wr_ext_b, hr_b, kin_data)

            @test (@ballocated f_ode!($dyn, $mp_b, $kin_data, $actions)) == 0

        end

    end #testset

end



end #module