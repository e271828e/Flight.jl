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
    @info("Called f_ode! with t = $(sys.t[]), x = $(sys.x[1]) and y = $(sys.y)")
    sys.ẋ .= 1/sys.constants.τ * (sys.u[] - sys.x[1])
    sys.y = sys.x[1]
end

function Systems.f_step!(sys::System{FirstOrder})
    x_new = sys.x[1] + 1
    # @info("Called f_step! at t = $(sys.t[]) and x = $(sys.x[1]), x updated to $(x_new)")
    # sys.x .= x_new
    #if we want the change in x to be reflected in y at the end of this step
end

function Systems.f_disc!(sys::System{FirstOrder}, ::Real)
    println("Called f_disc! at t = $(sys.t[]), got y = $(sys.y)")
end

Systems.init!(sys::System{FirstOrder}, x0::Real = 0.0) = (sys.x .= x0)


################################################################################
################################# TestSystem ###################################

#discrete system. it simply outputs simulation time and its IO-loopback echo
@kwdef struct TestSystem <: SystemDefinition end

@kwdef mutable struct TestSystemY
    output::Float64 = 0.0
    echo::Float64 = 0.0
end

Systems.U(::TestSystem) = Ref(0.0)
Systems.Y(::TestSystem) = TestSystemY()

function Systems.f_disc!(sys::TestSystem, ::Real)
    @info sys.y
    sys.y = TestSystemY(; output = sys.t[], echo = sys.u[])
end


################################################################################
################################## TestInput ###################################

#the internal Channel for this Device is NOT the one to be used by the SimInput
#wrapper, it is meant to be shared with a TestOutput{T} to emulate a loopback
@kwdef struct TestInput{T} <: InputDevice
    channel::Channel{T} = Channel{T}()
end

IODevices.get_data(input::TestInput) = take!(input.channel)
function Sim.assign_input!(sys::System{TestSystem}, data::Float64, ::DefaultMapping)
    sys.u[] = data
end

Base.Channel(::TestInput{T}, size::Int) where {T} = Channel{T}(size)


################################################################################
################################## TestOutput ##################################

@kwdef struct TestOutput{T} <: OutputDevice
    channel::Channel{T} = Channel{T}()
end

IODevices.handle_data(output::TestOutput{T}, ::T) = put!(output.channel)
function Sim.extract_output(sys::System{TestSystem}, data::Float64, ::DefaultMapping)
    sys.u[] = data
end

Base.Channel(::TestOutput{T}, size::Int) where {T} = Channel{T}(size)





function test_sim_standalone()

    sys = FirstOrder() |> System
    sim = Simulation(sys; dt = 0.1, Δt = 1.0, t_end = 5)
    x0 = 1.0
    f_init! = (sys)->Systems.init!(sys, x0)
    reinit!(sim, f_init!)
    return sim
    # Sim.run!(sim)

end

end #module
