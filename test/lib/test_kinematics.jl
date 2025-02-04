module TestKinematics

using Test
using LinearAlgebra
using BenchmarkTools

using Flight.FlightCore

using Flight.FlightLib.Geodesy
using Flight.FlightLib.Kinematics

export test_kinematics

function test_kinematics()

    @testset verbose = true "Kinematics" begin

        sys_ECEF = System(ECEF())
        sys_WA = System(WA())
        sys_NED = System(NED())

        @testset verbose = true "Performance" begin

            @test (@ballocated f_ode!($sys_ECEF)) == 0
            @test (@ballocated f_ode!($sys_WA)) == 0
            @test (@ballocated f_ode!($sys_NED)) == 0

            sys_WA.x.q_wb[1] = 3 #force renormalization in f_step!
            sys_ECEF.x.q_eb[1] = 3 #force renormalization in f_step!
            @test @ballocated(f_step!($sys_WA)) == 0
            @test @ballocated(f_step!($sys_ECEF)) == 0
            @test @ballocated(f_step!($sys_NED)) == 0

        end

        @testset verbose = true "Initialization" begin

            kin_init = KinInit(
                loc = LatLon(π/3, -π/6),
                h = HOrth(12354),
                ω_wb_b = [0.1, 0.1, -0.2],
                v_eOb_n = [100, 10, -4])

            Systems.init!(sys_ECEF, kin_init)
            Systems.init!(sys_WA, kin_init)
            Systems.init!(sys_NED, kin_init)

            #let the kinematic state propagate to y
            f_ode!(sys_ECEF)
            f_ode!(sys_WA)
            f_ode!(sys_NED)

            #check the initialization yields consistent results between implementations
            @test sys_ECEF.y.q_nb ≈ sys_WA.y.q_nb
            @test sys_ECEF.y.n_e ≈ sys_WA.y.n_e
            @test sys_ECEF.y.h_e ≈ sys_WA.y.h_e
            @test sys_ECEF.y.v_eOb_b ≈ sys_WA.y.v_eOb_b
            @test sys_ECEF.y.ω_eb_b ≈ sys_WA.y.ω_eb_b

            #check that direct KinData initialization is equivalent
            kin_data = KinData(kin_init)
            @test kin_data.q_nb ≈ sys_WA.y.q_nb
            @test kin_data.n_e ≈ sys_WA.y.n_e
            @test kin_data.h_e ≈ sys_WA.y.h_e
            @test kin_data.v_eOb_b ≈ sys_WA.y.v_eOb_b
            @test kin_data.ω_eb_b ≈ sys_WA.y.ω_eb_b

        end

        @testset verbose = true "Simulation" begin

            sim_ECEF = Simulation(sys_ECEF; t_end = 20);
            sim_WA = Simulation(sys_WA; t_end = 20);
            sim_NED = Simulation(sys_NED; t_end = 20);

            Sim.run!(sim_ECEF)
            Sim.run!(sim_WA)
            Sim.run!(sim_NED)

            #### validate WA against ECEF
            @test sys_ECEF.y.q_nb ≈ sys_WA.y.q_nb
            @test sys_ECEF.y.v_eOb_n ≈ sys_WA.y.v_eOb_n
            @test sys_ECEF.y.h_e ≈ sys_WA.y.h_e

            #### validate NED against WA
            @test sys_WA.y.q_nb ≈ sys_NED.y.q_nb
            @test sys_WA.y.v_eOb_n ≈ sys_NED.y.v_eOb_n
            @test sys_WA.y.h_e ≈ sys_NED.y.h_e

        end

    end #testset

end



end #module