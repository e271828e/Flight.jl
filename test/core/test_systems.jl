module TestSystems

using Test, UnPack, Logging, StructTypes, JSON3

using Flight.FlightCore

################################################################################

@kwdef struct FirstOrder <: SystemDefinition
    τ::Float64 = 1.0
end

Systems.X(::FirstOrder) = [0.0]
Systems.U(::FirstOrder) = Ref(0.0)
Systems.Y(::FirstOrder) = 0.0

function Systems.f_ode!(sys::System{FirstOrder})
    # @info("Called f_ode! with t = $(sys.t[]), x = $(sys.x[1]) and y = $(sys.y)")
    sys.ẋ .= 1/sys.constants.τ * (sys.u[] - sys.x[1])
    sys.y = sys.x[1]
end

function Systems.f_disc!(::NoScheduling, sys::System{FirstOrder})
    x_new = sys.x[1] + 0.1
    @info("Called f_disc! at t = $(sys.t[]), updating x = $(sys.x[1]) to x = $(x_new)")
    sys.x .= x_new
    sys.y = sys.x[1]
    # println("Called f_disc! at t = $(sys.t[]), got y = $(sys.y)")
end

Systems.init!(sys::System{FirstOrder}, x0::Real = 0.0) = (sys.x .= x0)


@kwdef struct Node <: SystemDefinition
    a::FirstOrder = FirstOrder()
    b::Subsampled{FirstOrder} = Subsampled(FirstOrder(), 2)
end

Systems.init!(sys::System{Node}, x0::Real = 0.0) = (sys.x .= x0)

@kwdef struct Root <: SystemDefinition
    a::FirstOrder = FirstOrder()
    b::Subsampled{FirstOrder} = Subsampled(FirstOrder(), 2)
    c::Subsampled{Node} = Subsampled(Node(), 3)
end

Systems.init!(sys::System{Root}, x0::Real = 0.0) = (sys.x .= x0)


function test_multirate()
    sys = Root() |> System;
    sim = Simulation(sys; Δt = .25)
    Sim.run!(sim)
    return TimeSeries(sim)
end


end
