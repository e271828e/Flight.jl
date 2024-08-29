module TestSim

using Test

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


#lesson learned: for a loopback test, it's essential that the input and output
#devices are not mutually locking. otherwise, at least one of them will block
#irrecoverably when the simulation terminates. this coupling may happen for
#example if input and output share a loopback Channel and they make blocking
#put! and take! calls on it.

#to avoid this, it is enough that at least one of them can only block when
#waiting on its SimInterface, but not on the loopback interface. this is
#the case with an UDP loopback, in which the UDPOutput may block when calling
#take! on the SimInterface Channel, but not on its send() call, which is
#nonblocking.


################################################################################
################################# TestSystem ###################################

@kwdef struct TestSystem <: SystemDefinition end

@kwdef mutable struct TestSystemU
    output::Float64 = 0.0
    echo::Float64 = 0.0
end

@kwdef struct TestSystemY
    output::Float64 = 0.0
    echo::Float64 = 0.0
end

Systems.U(::TestSystem) = TestSystemU()
Systems.Y(::TestSystem) = TestSystemY()

function Systems.f_disc!(sys::System{<:TestSystem}, ::Real)
    sys.u.output += 1
    sys.y = TestSystemY(; output = sys.u.output, echo = sys.u.echo)
end


################################################################################
################################## TestInput ###################################

@kwdef struct TestInput{T} <: InputDevice
    channel::Channel{T} = Channel{T}(1)
end

IODevices.init!(::TestInput) = nothing
IODevices.shutdown!(::TestInput) = nothing

#may block on the loopback channel
IODevices.get_data(input::TestInput) = Sim.take!(input.channel)

function Sim.assign_input!(sys::System{TestSystem}, data::Float64, ::DefaultMapping)
    sys.u.echo = data
    # @show sys.u.echo
end

Base.Channel(::TestInput{T}, size::Int) where {T} = Channel{T}(size)


################################################################################
################################## TestOutput ##################################

@kwdef struct TestOutput{T} <: OutputDevice
    channel::Channel{T} = Channel{T}(1)
end

IODevices.init!(::TestOutput) = nothing
IODevices.shutdown!(::TestOutput) = nothing

#will NOT block on the loopback channel
function IODevices.handle_data(output::TestOutput{T}, v::T) where {T}
    Sim.put_nonblocking!(output.channel, v)
end

function Sim.extract_output(sys::System{TestSystem}, ::DefaultMapping)
    # @show sys.y.output
    return sys.y.output
end

Base.Channel(::TestOutput{T}, size::Int) where {T} = Channel{T}(size)

################################################################################


function test_basic_io()

    @testset verbose = true "Channel loopback" begin

        c = Channel{Float64}(1)
        sys = TestSystem() |> System
        sim = Simulation(sys; t_end = 1.0)
        Sim.attach!(sim, TestInput(c))
        Sim.attach!(sim, TestOutput(c))
        input_interface = sim.interfaces[1]
        output_interface = sim.interfaces[2]

        step!(sim)
        @test sim.y.output == 1.0
        @test sim.y.echo == 0.0

        #extracts sys.y.output as a Float64 and puts it into the output
        #SimInterface's channel
        Sim.update!(output_interface, sys)
        @test isready(output_interface.channel)

        #TestOutput takes this value from its parent SimInterface's channel and
        #handles it, which here means putting it into the loopback Channel
        Sim.update!(output_interface)
        @test !isready(output_interface.channel)
        @test isready(output_interface.device.channel)

        #TestInput gets new input data, which here means taking it from the
        #loopback Channel, and then puts it into its parent SimInterface's
        #channel
        Sim.update!(input_interface)
        @test !isready(input_interface.device.channel)
        @test isready(input_interface.channel)

        #reads the pending input from the input SimInterface's channel and assigns
        #it to the System as specified by the SimInterface's mapping
        Sim.update!(input_interface, sys)
        @test !isready(input_interface.channel)

        @test sim.u.echo == 1.0
        step!(sim)
        @test sim.y.output == 2.0
        @test sim.y.echo == 1.0

        #run to completion
        Sim.run!(sim)
        #we should expect the echo value to have increased, but not necessarily
        #to match output-1, because due to the non-blocking SimOutput policy
        #some of the sim.y updates may be skipped

        return sim

    end

end





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
