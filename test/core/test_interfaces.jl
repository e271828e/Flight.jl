module TestInterfaces

using Test
using BenchmarkTools
using UnPack
using ComponentArrays

using Flight

export test_io

struct TestTarget end

function test_interfaces()

    sys = PICompensator{1}(k_i = 0.5) |> System
    sim = Simulation(sys; t_end = 20, dt = 0.02)
    joy = XBoxController()
    # joysticks = get_connected_joysticks()
    # input = InputManager(sim, joysticks, TestPICompensatorMapping())


    # @sync begin
    #     Threads.@spawn Sim.run!(sim; rate = 1, verbose = true)
    #     Threads.@spawn Sim.run!(input; verbose = true)
    # end

    return sim

end



end #module