using Random
using UnPack

using Flight
using Flight.Components.Discrete: SecondOrder


function demo_stochastic()

    k_u = [1.0, 1.0]
    k_av = [-0.3, -0.3]
    k_ap = [-0.02, -0.02]

    sys = SecondOrder{2}(; k_u, k_av, k_ap) |> System

    rng_init = Xoshiro(0)
    rng_io = Xoshiro(0)
    σ_0 = ComponentVector(v = [0.1, 0.1], p = [1, 1])

    sys_reinit! = let rng_init = rng_init, rng_io = rng_io, σ_0 = σ_0
        function (sys; init_seed = 0, io_seed = 0)
            #set seeds for both initialization and simulation
            @unpack x, u = sys

            #randomize initial state
            Random.seed!(rng_init, init_seed)
            Random.randn!(rng_init, x)
            x.v .*= σ_0.v #scale velocity components with their initial σ
            x.p .*= σ_0.p #scale position components with their initial σ

            #reset the io rng, and assign u for the first step (σ_u is
            #controlled by k_u)
            Random.seed!(rng_io, io_seed)
            Random.randn!(rng_io, u)
        end
    end

    sys_io! = let rng_io = rng_io
        #apply unit noise input, k_u controls σ_u
        function (u, y, t, params)
            Random.randn!(rng_io, u)
        end
    end

    sim = Simulation(sys; t_end = 100, dt = 0.02, Δt = 0.02, sys_reinit!, sys_io!)
    reinit!(sim; init_seed = 0, io_seed = 0) #apply sys_reinit! before running

    return sim

end