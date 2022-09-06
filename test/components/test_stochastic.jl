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

        T_c = [1.0]
        k_w = [1.0]

        sys = OrnsteinUhlenbeck{1}(; T_c, k_w) |> System

        #test for allocations
        @test (@ballocated(f_disc!($sys, 0.1)) == 0)

        rng_init = Xoshiro(0)
        rng_io = Xoshiro(0)
        σ_0 = Stochastic.σ(sys) #set the initialization σ equal to the stationary σ

        sys_reinit! = let rng_init = rng_init, rng_io = rng_io, σ_0 = σ_0
            function (sys; init_seed = 0, io_seed = 0)
                #set seeds for both initialization and simulation
                @unpack x, u = sys

                #randomize initial state
                Random.seed!(rng_init, init_seed)
                Random.randn!(rng_init, x)
                x .*= σ_0.x #scale N(0,1) with initialization σ

                #reset the io rng, and assign u for the first step (σ_u is
                #controlled by k_w)
                Random.seed!(rng_io, io_seed)
                Random.randn!(rng_io, u)
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