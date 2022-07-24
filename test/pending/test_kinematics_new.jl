module TestKinematicsNew

using Test
using LinearAlgebra
using BenchmarkTools

using Flight.Systems
using Flight.Sim
using Flight.Kinematics
using Flight.KinematicsNew
using Flight.Attitude
using Flight.Geodesy

export test_kinematics

function test_kinematics()

    sys_LTF = System(KinLTF())
    sys_WA = System(WA())
    sys_NED = System(NED())

    v_eOb_b = [100, 10, -4]
    ω_eb_b = [0.1, 0.1, -0.2]
    h_e = 12354

    sys_LTF.x.vel.v_eOb_b = v_eOb_b
    sys_LTF.x.vel.ω_eb_b = ω_eb_b
    sys_LTF.x.pos.h_e = h_e
    sys_WA.x.vel.v_eOb_b = v_eOb_b
    sys_WA.x.vel.ω_eb_b = ω_eb_b
    sys_WA.x.pos.h_e = h_e
    sys_NED.x.vel.v_eOb_b = v_eOb_b
    sys_NED.x.vel.ω_eb_b = ω_eb_b
    sys_NED.x.pos.h_e = h_e

    # println(@btime f_cont!($sys_LTF))
    # println(@btime f_cont!($sys_WA))
    # println(@btime f_cont!($sys_NED))
    # return

    #### validate WA against LTF

    sim_LTF = Simulation(sys_LTF; t_end = 20);
    sim_WA = Simulation(sys_WA; t_end = 20);

    Sim.run!(sim_LTF)
    Sim.run!(sim_WA)

    println(sys_LTF.y.pos.q_nb ∘ sys_WA.y.pos.q_nb')
    println(sys_LTF.y.pos.h_e - sys_WA.y.pos.h_e)

    #### validate NED against WA

    kin_init = KinematicsNew.InitialCondition(
        Ob = GeographicLocation(LatLon(π/3, -π/6), AltE(12354)),
        v_eOb_n = [100, 10, -4])

    KinematicsNew.init!(sys_WA, kin_init)
    KinematicsNew.init!(sys_NED, kin_init)

    println(sys_WA.y.vel)
    println(sys_NED.y.vel)

    sim_WA = Simulation(sys_WA; t_end = 20);
    sim_NED = Simulation(sys_NED; t_end = 20);

    Sim.run!(sim_WA)
    Sim.run!(sim_NED)

    println(sys_WA.y.pos.q_nb ∘ sys_NED.y.pos.q_nb')
    println(sys_WA.y.pos.h_e - sys_NED.y.pos.h_e)

end



end #module