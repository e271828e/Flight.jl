module TestSystems

using Test, UnPack, Logging, StructTypes

using Flight.FlightCore

################################################################################
################################ Multi-rate ####################################

@kwdef struct FirstOrder <: SystemDefinition
    τ::Float64 = 1.0
end

Systems.X(::FirstOrder) = [0.0]
Systems.U(::FirstOrder) = Ref(0.0)
Systems.Y(::FirstOrder) = 0.0

function Systems.f_ode!(sys::System{FirstOrder})
    # @info("Called f_ode! with t = $(sys.t[]), x = $(sys.x[1]) and y = $(sys.y)")
    sys.ẋ .= 1/sys.τ * (sys.u[] - sys.x[1])
    sys.y = sys.x[1]
end

function Systems.f_disc!(::NoScheduling, sys::System{FirstOrder})
    x_new = sys.x[1] + 0.1
    @info("Called f_disc! at t = $(sys.t[]), n = $(sys.n[]), updating x = $(sys.x[1]) to x = $(x_new)")
    sys.x .= x_new
    sys.y = sys.x[1]
    # println("Called f_disc! at t = $(sys.t[]), got y = $(sys.y)")
end

@no_step FirstOrder

################################################################################

@kwdef struct Node <: SystemDefinition
    a::FirstOrder = FirstOrder()
    b::Subsampled{FirstOrder} = Subsampled(FirstOrder(), 2)
end

@ss_dynamics Node

################################################################################

@kwdef struct Root <: SystemDefinition
    a::FirstOrder = FirstOrder()
    b::Subsampled{FirstOrder} = Subsampled(FirstOrder(), 2)
    c::Subsampled{Node} = Subsampled(Node(), 3)
end

@ss_dynamics Root

function Systems.init!(sys::System{Root}, x0::Real = 0.0)
    (sys.x .= x0)
    # f_disc!(sys)
end

################################################################################

function test_multirate()
    sys = Root() |> System;
    sim = Simulation(sys; Δt = 1.0, t_end = 30)
    # Sim.init!(sim)
    Sim.run!(sim)
    # return TimeSeries(sim)
    return sim
end


################################################################################
########################### Discrete Dynamics ##################################

struct DiscreteTestComponent <: SystemDefinition end

@kwdef struct DiscreteTestComponentY
    a::Float64 = 0
    b::Float64 = 0
end

Systems.X(::DiscreteTestComponent) = ComponentVector(a = 0.0, b = 0.0)
Systems.Y(::DiscreteTestComponent) = DiscreteTestComponentY()

function Systems.f_disc!(::NoScheduling, sys::System{DiscreteTestComponent})
    sys.x.a += 1
    sys.x.b -= 1
    sys.y = DiscreteTestComponentY(a = x.a, b = x.b)
end

function test_discrete()

    #if we set a fixed dt < Δt and adaptive = false, the integrator may take
    #multiple unnecessary steps between discrete update epochs
    sys = DiscreteTestComponent() |> System
    sim = Simulation(sys, adaptive = false, Δt = 1.0)
    step!(sim, 1, true)
    @show sim.integrator.iter
    @show sim.t
    @show sim.y

    #if we set adaptive = true, it will only take a few intermediate steps
    #before the integrator extends the proposed dt beyond Δt, on account of ẋ
    #always being 0. from that moment on, it only stops at the discrete update
    #epochs
    sys = DiscreteTestComponent() |> System
    sim = Simulation(sys, adaptive = true, Δt = 1.0)
    step!(sim, 1, true)
    @show sim.integrator.iter
    @show sim.t
    @show sim.y

    #here, we set dt = Δt directly, so it only stops at discrete update epochs
    #right from the start
    sys = DiscreteTestComponent() |> System
    sim = Simulation(sys, Δt = 1.0, dt = Δt)
    step!(sim, 1, true)
    @show sim.integrator.iter
    @show sim.t
    @show sim.y

    #setting dt > Δt also works: the integrator will still honor the discrete
    #callback
    sys = DiscreteTestComponent() |> System
    sim = Simulation(sys, Δt = 1.0, dt = 2.0)
    step!(sim, 1, true)
    @show sim.integrator.iter
    @show sim.t
    @show sim.y

end

end
