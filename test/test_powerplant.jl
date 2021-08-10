# module TestPowerplant

using Test
using DifferentialEquations
using Flight.Powerplant

export test_pwp

#integrate over Δt = 1, stopping exactly, and taking as many intermediate
#steps as dictated by the adaptive algorithm
f_benchmark(sys) = (step!(sys, 1, true); reinit!(sys.integrator))

function test_pwp()

    #real time simulation setup, save at all steps
    sys1 = ElectricThrusterSystem(
        method = Heun(),
        adaptive = false,
        dt = 0.01)

    #real time simulation setup, save at specified steps
    sys2 = ElectricThrusterSystem(
        method = Heun(),
        adaptive = false,
        dt = 0.01,
        output_saveat = 0:0.1:1)

    #non-real time simulation setup, save at all steps
    sys4 = ElectricThrusterSystem(
        method = Tsit5(),
        adaptive = true)

    #non-real time simulation setup, save at specified steps (forcing tstops)
    sys3 = ElectricThrusterSystem(
        method = Tsit5(),
        adaptive = true,
        output_saveat = 0:0.1:1)

    @btime f_benchmark($sys1)
    @btime f_benchmark($sys2)
    @btime f_benchmark($sys3)
    @btime f_benchmark($sys4)

end