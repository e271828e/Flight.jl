module TestStochastic

using Test
using BenchmarkTools
using UnPack
using ComponentArrays
using StaticArrays
using Random
using Statistics

using Flight

export test_stochastic

function test_stochastic()
    @testset verbose = true "Stochastic" begin
        test_ou()
    end
end

function test_ou()

    @testset verbose = true "Ornstein-Uhlenbeck" begin

        T_c = 1.0
        k_w = 1.0

        sys = OrnsteinUhlenbeck{1}(; T_c, k_w) |> System

        @test (@ballocated(f_disc!($sys, 0.1)) == 0)

        σ0 = Stochastic.σ(sys)
        rng_init = Xoshiro(0) #for state initialization
        rng_io = Xoshiro(0) #drives the white noise input

        sys_reinit! = let rng_init = rng_init, rng_io = rng_io
            function (sys; init_seed = 0, io_seed = 0)
                #randomize initial state
                Random.seed!(rng_init, init_seed)
                Random.randn!(rng_init, sys, σ0)

                #set seed for simulation
                Random.seed!(rng_io, io_seed)
            end
        end

        #pass rng_io as additional argument to the System's f_disc!
        sim = Simulation(sys; t_end = 1000, dt = 1, Δt = 1,
                              args_disc = (rng_io,), sys_reinit!)

        #run 1000 trajectories and return their sample variances
        vars = map(1:1000) do seed
            reinit!(sim; init_seed = seed, io_seed = seed)
            Sim.run!(sim)
            th_scalar = collect(get_components(TimeHistory(sim)))[1]
            Statistics.var(th_scalar._data)
        end

        #experimental variance should match the theoretical stationary variance
        @test mean(vars) ≈ Stochastic.σ²(sys)[1] atol = 1e-3

    end

end

end