module TestSim

using Flight.FlightCore

export test_sim_standalone
################################################################################
############################### FirstOrder #####################################

@kwdef struct FirstOrder <: SystemDefinition
    τ::Float64 = 1.0
end

Systems.X(::FirstOrder) = [0.0]
Systems.U(::FirstOrder) = Ref(0.0)
Systems.Y(::FirstOrder) = 0.0

function Systems.f_ode!(sys::System{FirstOrder})
    sys.ẋ .= 1/sys.constants.τ * (sys.u[] - sys.x[1])
    # println("Called f_ode! at t = $(sys.t[])")
end

function Systems.f_disc!(sys::System{FirstOrder}, ::Real)
    println("Called f_disc! at t = $(sys.t[]), got y = $(sys.y)")
end

function Systems.f_step!(sys::System{FirstOrder})
    sys.y = sys.x[1]
    println("Called f_step! at t = $(sys.t[]), y is now $(sys.y)")
end

Systems.init!(sys::System{FirstOrder}, x0::Real = 0.0) = (sys.x .= x0)


function test_sim_standalone()

    sys = FirstOrder() |> System
    sim = Simulation(sys; dt = 0.1, Δt = 1.0, t_end = 5)
    x0 = 1.0
    f_init! = Systems.init!(sys, x0)
    reinit!(sim, f_init!)
    return sim
    # Sim.run!(sim)

end

end #module
