using Test
using DifferentialEquations
using Flight.Propulsion
using Flight.Model

export test_pwp

#integrate over Δt = 1, stopping exactly, and taking as many intermediate
#steps as dictated by the adaptive algorithm
f_benchmark(mdl) = (reinit!(mdl.integrator); step!(mdl, 1, true))

function test_pwp()

    #NRT setup, save at all adaptive time steps
    mdl1 = ContinuousModel(EThruster())

    #NRT setup, save at specified time steps
    mdl2 = ContinuousModel(EThruster(), y_saveat = 0:0.1:1)

    #RT setup, save at specified time steps (should be multiple of dt to avoid
    #intermediate steps)
    mdl3 = ContinuousModel(EThruster(); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0:0.1:1)

    @btime f_benchmark($mdl1)
    @btime f_benchmark($mdl2)
    @btime f_benchmark($mdl3)

end