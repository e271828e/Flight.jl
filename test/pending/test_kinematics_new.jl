module TestKinematicsNew

using Test
using LinearAlgebra
using BenchmarkTools

using Flight.Systems
using Flight.Sim
using Flight.KinematicsNew
using Flight.Attitude
using Flight.Geodesy

export test_kinematics

function test_kinematics()

    @testset verbose = true "Kinematics" begin

    sys_ECEF = System(ECEF())
    sys_LTF = System(LTF())
    sys_NED = System(NED())

    @test (@ballocated f_cont!($sys_ECEF)) == 0
    @test (@ballocated f_cont!($sys_LTF)) == 0
    @test (@ballocated f_cont!($sys_NED)) == 0

    kin_init = KinematicsNew.InitialCondition(
        Ob = GeographicLocation(LatLon(π/3, -π/6), AltE(12354)),
        ω_lb_b = [0.1, 0.1, -0.2],
        v_eOb_n = [100, 10, -4])

    KinematicsNew.init!(sys_ECEF, kin_init)
    KinematicsNew.init!(sys_LTF, kin_init)
    KinematicsNew.init!(sys_NED, kin_init)

    #let the kinematic state propagate to y
    f_cont!(sys_ECEF)
    f_cont!(sys_LTF)
    f_cont!(sys_NED)

    #check the initialization yields consistent results between implementations
    @test sys_ECEF.y.q_nb ≈ sys_LTF.y.q_nb
    @test sys_ECEF.y.n_e ≈ sys_LTF.y.n_e
    @test sys_ECEF.y.h_e ≈ sys_LTF.y.h_e
    @test sys_ECEF.y.v_eOb_b ≈ sys_LTF.y.v_eOb_b
    @test sys_ECEF.y.ω_eb_b ≈ sys_LTF.y.ω_eb_b

    sim_ECEF = Simulation(sys_ECEF; t_end = 20);
    sim_LTF = Simulation(sys_LTF; t_end = 20);
    sim_NED = Simulation(sys_NED; t_end = 20);

    Sim.run!(sim_ECEF)
    Sim.run!(sim_LTF)
    Sim.run!(sim_NED)

    #### validate LTF against ECEF
    @test sys_ECEF.y.q_nb ≈ sys_LTF.y.q_nb
    @test sys_ECEF.y.v_eOb_n ≈ sys_LTF.y.v_eOb_n
    @test sys_ECEF.y.h_e ≈ sys_LTF.y.h_e

    #### validate NED against LTF
    @test sys_LTF.y.q_nb ≈ sys_NED.y.q_nb
    @test sys_LTF.y.v_eOb_n ≈ sys_NED.y.v_eOb_n
    @test sys_LTF.y.h_e ≈ sys_NED.y.h_e

    end #testset

end



end #module