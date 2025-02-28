module TestDynamics

using Test
using LinearAlgebra
using BenchmarkTools
using StaticArrays

using Flight.FlightCore
using Flight.FlightLib

export test_dynamics

#define a dummy component that returns the dynamic data at b to be tested
struct TestComponent <: AbstractComponents end

@kwdef mutable struct TestComponentU
    mp_Σ_b::MassProperties = MassProperties()
    wr_Σ_b::Wrench = Wrench()
    ho_Σ_b::SVector{3, Float64} = zeros(3)
end

Systems.U(::TestComponent) = TestComponentU()
Dynamics.get_mp_b(cmp::System{TestComponent}) = cmp.u.mp_Σ_b
Dynamics.get_wr_b(cmp::System{TestComponent}) = cmp.u.wr_Σ_b
Dynamics.get_hr_b(cmp::System{TestComponent}) = cmp.u.ho_Σ_b

function test_dynamics()


    @testset verbose = true "Dynamics" begin

        cmp = System(TestComponent())
        dyn = System(VehicleDynamics())
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
        #center of mass c
        mp_Σ_c = MassProperties(RigidBodyDistribution(1.0, diagm(ones(3))))
        cmp.u.mp_Σ_b = mp_Σ_c

        @testset verbose = true "Essentials" begin

            #start with Ob=Oc
            mp_Σ_b = mp_Σ_c
            cmp.u.mp_Σ_b = mp_Σ_b

            cmp.u.wr_Σ_b = Wrench(F = [1, 2, 1])
            f_ode!(dyn, cmp, kin_data)
            @test all(dyn.ẋ.ω_eb_b .≈ 0)
            @test all(dyn.ẋ.v_eb_b .≈ [1, 2, 1 + gravity(kin_data)])
            @test all(dyn.y.a_eb_b .≈ [1, 2, 1 + gravity(kin_data)])
            @test all(dyn.y.a_ib_b .≈ [1, 2, 1 + G_n(kin_data)[3]])

            #now let Oc be located 1 meter ahead from Ob along the x axis
            r_bc_b = [1,0,0]
            t_bc = FrameTransform(r = r_bc_b) #b to c
            mp_Σ_b = t_bc(mp_Σ_c)
            cmp.u.mp_Σ_b = mp_Σ_b

            cmp.u.wr_Σ_b = Wrench(F = [0, 0, 1], τ = zeros(3))

            #we expect a positive unit angular acceleration around y_b, which
            #will also add to the linear acceleration of frame b
            f_ode!(dyn, cmp, kin_data)
            @test dyn.ẋ.ω_eb_b[2] .≈ 1
            @test dyn.ẋ.v_eb_b[3] ≈ 2 + gravity(kin_data)

        end

        @testset verbose = true "Performance" begin

            @test (@ballocated f_ode!($dyn, $cmp, $kin_data)) == 0

        end

    end #testset

end



end #module