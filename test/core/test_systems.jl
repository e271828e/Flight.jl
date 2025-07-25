module TestSystems

using Test, UnPack, Logging, StructTypes

using Flight.FlightCore

################################################################################
################################ Multi-rate ####################################

@kwdef struct FirstOrder <: ModelDefinition
    τ::Float64 = 1.0
end

Modeling.X(::FirstOrder) = [0.0]
Modeling.U(::FirstOrder) = Ref(0.0)
Modeling.Y(::FirstOrder) = 0.0

function Modeling.f_ode!(mdl::Model{FirstOrder})
    # @info("Called f_ode! with t = $(mdl.t[]), x = $(mdl.x[1]) and y = $(mdl.y)")
    mdl.ẋ .= 1/mdl.τ * (mdl.u[] - mdl.x[1])
    mdl.y = mdl.x[1]
end

function Modeling.f_disc!(::NoScheduling, mdl::Model{FirstOrder})
    x_new = mdl.x[1] + 0.1
    @info("Called f_disc! at t = $(mdl.t[]), _n = $(mdl._n[]), _N = $(mdl._N), updating x = $(mdl.x[1]) to x = $(x_new)")
    mdl.x .= x_new
    mdl.y = mdl.x[1]
    # println("Called f_disc! at t = $(mdl.t[]), got y = $(mdl.y)")
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

function Modeling.init!(mdl::Model{Root}, x0::Real = 0.0)
    (mdl.x .= x0)
    # f_disc!(mdl)
end

################################################################################

function test_single()
    mdl = FirstOrder() |> Model;
    sim = Simulation(mdl; Δt = 1.0, t_end = 30)
    # Sim.init!(sim)
    Sim.run!(sim)
    # return TimeSeries(sim)
    return sim
end

function test_multirate()
    mdl = Root() |> Model;
    sim = Simulation(mdl; Δt = 1.0, t_end = 30)
    # Sim.init!(sim)
    Sim.run!(sim)
    # return TimeSeries(sim)
    return sim
end


################################################################################
########################### Discrete Dynamics ##################################

struct DiscreteTestComponent <: ModelDefinition end

@kwdef struct DiscreteTestComponentY
    a::Float64 = 0
    b::Float64 = 0
end

Modeling.X(::DiscreteTestComponent) = ComponentVector(a = 0.0, b = 0.0)
Modeling.Y(::DiscreteTestComponent) = DiscreteTestComponentY()

function Modeling.f_disc!(::NoScheduling, mdl::Model{DiscreteTestComponent})
    mdl.x.a += 1
    mdl.x.b -= 1
    mdl.y = DiscreteTestComponentY(a = x.a, b = x.b)
end

function test_discrete()

    #if we set a fixed dt < Δt and adaptive = false, the integrator may take
    #multiple unnecessary steps between discrete update epochs
    mdl = DiscreteTestComponent() |> Model
    sim = Simulation(mdl, adaptive = false, Δt = 1.0)
    step!(sim, 1, true)
    @show sim.integrator.iter
    @show sim.t
    @show sim.y

    #if we set adaptive = true, it will only take a few intermediate steps
    #before the integrator extends the proposed dt beyond Δt, on account of ẋ
    #always being 0. from that moment on, it only stops at the discrete update
    #epochs
    mdl = DiscreteTestComponent() |> Model
    sim = Simulation(mdl, adaptive = true, Δt = 1.0)
    step!(sim, 1, true)
    @show sim.integrator.iter
    @show sim.t
    @show sim.y

    #here, we set dt = Δt directly, so it only stops at discrete update epochs
    #right from the start
    mdl = DiscreteTestComponent() |> Model
    sim = Simulation(mdl, Δt = 1.0, dt = Δt)
    step!(sim, 1, true)
    @show sim.integrator.iter
    @show sim.t
    @show sim.y

    #setting dt > Δt also works: the integrator will still honor the discrete
    #callback
    mdl = DiscreteTestComponent() |> Model
    sim = Simulation(mdl, Δt = 1.0, dt = 2.0)
    step!(sim, 1, true)
    @show sim.integrator.iter
    @show sim.t
    @show sim.y

end

end
