module TestDynamics

using Test
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using UnPack

using Flight.FlightCore
using Flight.FlightLib

export test_dynamics


function test_dynamics()

    @testset verbose = true "Dynamics" begin

        dyn = Model(VehicleDynamics())

        kin_init = KinInit(
            location = LatLon(0.4, 0.5),
            h = HOrth(1000),
            ω_wb_b = [1.0, 1.0, 0.0],
            q_nb = REuler(1.0, 0.5, 0.4),
            v_eb_n = [100, 0, 0])

        kin_data = KinData(kin_init)
        @unpack q_eb, r_eb_e, q_nb = kin_data

        #set up a rigid body distribution with unit mass and unit inertia tensor
        #at its center of mass c
        mp_Σ_c = MassProperties(RigidBodyDistribution(1.0, diagm(ones(3))))

        @testset verbose = true "Essentials" begin

            #start with Ob=Oc
            dyn.u.mp_Σ_b = mp_Σ_c
            dyn.u.wr_Σ_b = Wrench(F = [1, 2, 1])
            @pack! dyn.u = q_eb, r_eb_e

            f_ode!(dyn)
            @test all(dyn.ẋ.ω_eb_b .≈ 0)
            @test all(dyn.ẋ.v_eb_b .≈ [1, 2, 1] +  q_nb'(g_n(Cartesian(r_eb_e))))
            @test all(dyn.y.a_eb_b .≈ dyn.ẋ.v_eb_b)
            @test all(dyn.y.a_ib_b .≈ [1, 2, 1] + q_nb'(G_n(Cartesian(r_eb_e))))

            #now let Oc be located 1 meter ahead from Ob along the x axis
            r_bc_b = [1,0,0]
            t_bc = FrameTransform(r = r_bc_b) #b to c

            dyn.u.mp_Σ_b = t_bc(mp_Σ_c)
            dyn.u.wr_Σ_b = Wrench(F = [0, 0, 1], τ = zeros(3))
            @pack! dyn.u = q_eb, r_eb_e

            #expect a positive unit angular acceleration around y_b
            f_ode!(dyn)
            @test dyn.ẋ.ω_eb_b[2] .≈ 1
            # @test dyn.ẋ.v_eb_b[3] ≈ 2 + q_nb'(g_n(Cartesian(r_eb_e)))[3]

        end

        @testset verbose = true "Performance" begin

            @test (@ballocated f_ode!($dyn)) == 0

        end

    end #testset

end



end #module