using Flight

#Example of how to define and simulate a purely discrete system

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

function f_dtc()

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
    sim = Simulation(sys, Δt = 1.0, dt = 1.0)
    step!(sim, 1, true)
    @show sim.integrator.iter
    @show sim.t
    @show sim.y

    #setting dt > Δt also works: the integrator will honor the discrete callback
    sys = DiscreteTestComponent() |> System
    sim = Simulation(sys, Δt = 1.0, dt = 2.0)
    step!(sim, 1, true)
    @show sim.integrator.iter
    @show sim.t
    @show sim.y

end
