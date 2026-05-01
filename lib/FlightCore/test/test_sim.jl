module TestSim

using Test
using Logging

using FlightCore

export test_sim

function test_sim()
    @testset verbose = true "Sim" begin

        Logging.disable_logging(Logging.Warn)
        @test_nowarn test_periodic()
        @test_nowarn test_single()
        @test_nowarn test_multirate()
        Logging.disable_logging(Logging.LogLevel(typemin(Int32)))

    end
end

################################################################################
########################### Periodic Dynamics ##################################

struct Periodic <: ModelDefinition end

Modeling.X(::Periodic) = [0.0]
Modeling.Y(::Periodic) = 0.0

@no_ode Periodic
@no_step Periodic

function Modeling.f_periodic!(::NoScheduling, mdl::Model{Periodic})
    mdl.x[1] += 1
    mdl.y = mdl.x[1]
end

function test_periodic()

    #if we set a fixed dt < Δt and adaptive = false, the integrator will take
    #multiple unnecessary steps between periodic updates
    mdl = Periodic() |> Model
    sim = Simulation(mdl; adaptive = false, dt = 0.02, Δt = 1.0)
    step!(sim, 2, true)
    @info sim.integrator.iter
    @info sim.t
    @info mdl.y

    #if we set adaptive = true, it will only take a few intermediate steps
    #for the integrator to extend the proposed dt beyond Δt, on account of ẋ
    #always being 0. from that moment on, it only stops at the periodic update
    #instants
    mdl = Periodic() |> Model
    sim = Simulation(mdl; adaptive = true, dt = 0.02, Δt = 1.0)
    step!(sim, 2, true)
    @info sim.integrator.iter
    @info sim.t
    @info mdl.y

    #here, we set dt = Δt directly, so it only stops at periodic updates right
    #from the start
    mdl = Periodic() |> Model
    sim = Simulation(mdl; dt = 1.0, Δt = 1.0)
    step!(sim, 2, true)
    @info sim.integrator.iter
    @info sim.t
    @info mdl.y

    #setting dt > Δt also works: the integrator will still honor the periodic
    #callback
    mdl = Periodic() |> Model
    sim = Simulation(mdl; Δt = 1.0, dt = 2.0)
    step!(sim, 2, true)
    @info sim.integrator.iter
    @info sim.t
    @info mdl.y

end

################################################################################
################################ Hybrid, Multi-rate ############################

@kwdef struct FirstOrder <: ModelDefinition
    τ::Float64 = 1.0
end

Modeling.X(::FirstOrder) = [0.0]
Modeling.U(::FirstOrder) = Ref(0.0)
Modeling.Y(::FirstOrder) = 0.0

function Modeling.f_ode!(mdl::Model{FirstOrder})
    @info("Called f_ode! with t = $(mdl.t[]), x = $(mdl.x[1]) and y = $(mdl.y)")
    mdl.ẋ .= 1/mdl.τ * (mdl.u[] - mdl.x[1])
    mdl.y = mdl.x[1]
end

function Modeling.f_periodic!(::NoScheduling, mdl::Model{FirstOrder})
    x_new = mdl.x[1] + 0.1
    @info("Called f_periodic! at t = $(mdl.t[]), _n = $(mdl._n[]), _N = $(mdl._N), updating x = $(mdl.x[1]) to x = $(x_new)")
    mdl.x .= x_new
    mdl.y = mdl.x[1]
    # println("Called f_periodic! at t = $(mdl.t[]), got y = $(mdl.y)")
end

@no_step FirstOrder

################################################################################

@kwdef struct Node <: ModelDefinition
    a::FirstOrder = FirstOrder()
    b::Subsampled{FirstOrder} = Subsampled(FirstOrder(), 2)
end

@sm_updates Node

################################################################################

@kwdef struct Root <: ModelDefinition
    a::FirstOrder = FirstOrder()
    b::Subsampled{FirstOrder} = Subsampled(FirstOrder(), 2)
    c::Subsampled{Node} = Subsampled(Node(), 3)
end

@sm_updates Root

function Modeling.f_init!(mdl::Model{Root}, x0::Real = 0.0)
    (mdl.x .= x0)
    # f_periodic!(mdl)
end

################################################################################

function test_single()
    mdl = FirstOrder() |> Model;
    sim = Simulation(mdl; Δt = 1.0, t_end = 30)
    # init!(sim)
    Sim.run!(sim)
    # return TimeSeries(sim)
    return sim
end

function test_multirate()
    mdl = Root() |> Model;
    sim = Simulation(mdl; Δt = 1.0, t_end = 30)
    # init!(sim)
    Sim.run!(sim)
    # return TimeSeries(sim)
    return sim
end

################################################################################
############################# Threading sketches ###############################

#compare the following:
# @time sleep(2)

# @time @async sleep(2)
# @time Threads.@spawn sleep(2)

# wait(@time @async sleep(2))
# wait(@time Threads.@spawn sleep(2))

# @time wait(@async sleep(2))
# @time wait(Threads.@spawn @sleep(2))


function threading_sketch()

    c = Channel{Int}(1)
    @sync begin
        Threads.@spawn begin
            while isopen(c)
                Core.println("Taken $(take!(c))")
            end
        end
        Threads.@spawn begin
            for i in 1:5
                sleep(1)
                @lock c begin
                    if !isready(c)
                        Core.println("Putting $i")
                        put!(c, i)
                    end
                end
            end
            close(c)
            Core.println("Bye")
        end

    end
end

end #module
