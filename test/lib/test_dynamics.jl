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
        kin = System(LTF())

        kin_init = KinInit(
            loc = LatLon(0, 0),
            h = HOrth(0),
            ω_lb_b = [0.0, 0.0, 0.0],
            q_nb = REuler(0.0, 0.0, 0.0),
            v_eOb_n = [0, 0, 0])

        Systems.init!(kin, kin_init)
        kin_data = KinData(kin)

        #set up a rigid body with unit mass and unit inertia tensor at its
        #center of mass Gb
        rbd_Gb = RigidBodyDistribution(1.0, diagm(ones(3)))
        mp_Gb = MassProperties(rbd_Gb)

        #let its frame origin Ob be located at Gb, so that its mass
        #properties are the same in both
        mp_Ob = mp_Gb


        @testset verbose = true "Essentials" begin


            #apply a force along each axis with zero torque
            wr_ext_Ob = Wrench(F = [1, 2, 1])
            hr_b = zeros(SVector{3, Float64})
            actions = Actions(th; mp_Ob, wr_ext_Ob, hr_b, kin_data)

            #we expect the angular acceleration to be zero, and linear acceleration
            #added to that of free fall
            f_ode!(dyn, mp_Ob, kin_data, actions)
            @test all(dyn.ẋ.ω_eb_b .≈ 0)
            @test all(dyn.ẋ.v_eOb_b .≈ [1, 2, 1 + gravity(kin_data)])

            f_ode!(dyn, mp_Ob, kin_data, actions)
            @test all(dyn.y.a_eOb_b .≈ [1, 2, 1 + gravity(kin_data)])
            @test all(dyn.y.a_iOb_b .≈ [1, 2, 1 + G_n(kin_data)[3]])

            #now let Gb be located 1 meter ahead from Ob along the x axis
            r_ObGb_b = [1,0,0]
            t_ObGb = FrameTransform(r = r_ObGb_b) #Ob to Gb
            mp_Ob = Dynamics.transform(t_ObGb, mp_Gb)

            #apply a 1N force along z_b
            wr_ext_Ob = Wrench(F = [0, 0, 1], M = zeros(3))
            hr_b = zeros(SVector{3, Float64})
            actions = Actions(th; mp_Ob, wr_ext_Ob, hr_b, kin_data)

            #we expect a positive unit angular acceleration around y_b, which
            #will also add to the linear acceleration at Ob
            f_ode!(dyn, mp_Ob, kin_data, actions)
            @test dyn.ẋ.ω_eb_b[2] .≈ 1
            @test dyn.ẋ.v_eOb_b[3] ≈ 2 + gravity(kin_data)

            r_ObGb_e = kin_data.q_eb'(r_ObGb_b)
            r_eGb_e = Cartesian(kin_data.r_eOb_e + r_ObGb_e)

            @test isapprox(dyn.y.a_eOb_b[3], 2 + gravity(r_eGb_e), atol = 1e-10)
            @test isapprox(dyn.y.a_iOb_b[3], 2 + G_n(r_eGb_e)[3], atol = 1e-10)
            @test all(isapprox.(dyn.y.f_Gb_b, [0, 0, 1], atol = 1e-10))

        end

        @testset verbose = true "Performance" begin

            wr_ext_Ob = Wrench(F = [1, 2, 1])
            hr_b = zeros(SVector{3, Float64})
            actions = Actions(th; mp_Ob, wr_ext_Ob, hr_b, kin_data)

            @test (@ballocated f_ode!($dyn, $mp_Ob, $kin_data, $actions)) == 0

        end

    end #testset

end



end #module