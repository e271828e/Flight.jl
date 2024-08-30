module TestSim

using Test, Logging

using Flight.FlightCore
using Flight.FlightPhysics

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
    # @info("Called f_ode! with t = $(sys.t[]), x = $(sys.x[1]) and y = $(sys.y)")
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
    # println("Called f_disc! at t = $(sys.t[]), got y = $(sys.y)")
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
    output::UInt8 = 0
    echo::UInt8 = 0
end

@kwdef struct TestSystemY
    output::UInt8 = 0
    echo::UInt8 = 0
end

Systems.U(::TestSystem) = TestSystemU()
Systems.Y(::TestSystem) = TestSystemY()

function Systems.f_disc!(sys::System{<:TestSystem}, ::Real)
    sleep(0.01)
    sys.y = TestSystemY(; output = sys.u.output, echo = sys.u.echo)
end

################################ UDP Loopback ##################################

struct UDPTestMapping <: IOMapping end

function Systems.assign_data!(sys::System{TestSystem}, data::Vector{UInt8},
                            ::UDPInput, ::UDPTestMapping)
    @debug "Got $data"
    sys.u.echo = data[1]
end

#needs to provide a Vector{UInt8}, which is what UDPOutput expects
function Systems.extract_data(sys::System{TestSystem},
                            ::UDPOutput, ::UDPTestMapping)
    data = [sys.y.output]
    @debug "Extracted $data"
    return data
end

function udp_loopback()

    @testset verbose = true "UDP Loopback" begin

        port = 14149
        sys = TestSystem() |> System
        sim = Simulation(sys; t_end = 1.0)
        sim.u.output = UInt8(37)
        Sim.attach!(sim, UDPInput(; port), UDPTestMapping())
        Sim.attach!(sim, UDPOutput(; port), UDPTestMapping())
        Sim.run!(sim)

        #sys.u.output must have propagated to sys.y.output within f_disc!, then
        #to sys.u.echo via loopback, and finally to sys.y.echo within f_disc!
        @test sim.y.output === UInt8(37)

        return sim

    end

end

################################ XPC Loopback ##################################

function Systems.extract_data(::System{TestSystem}, ::XPCClient, ::IOMapping)
    data = KinData() |> XPCPosition
    @debug "Extracted $data"
    return data
end

function xpc_loopback()

    @testset verbose = true "UDP Loopback" begin

        port = 14149
        sys = TestSystem() |> System
        sim = Simulation(sys; t_end = 1.0)
        Sim.attach!(sim, UDPInput(; port), UDPTestMapping())
        Sim.attach!(sim, XPCClient(; port))
        Sim.run!(sim)

        cmd = KinData() |> XPCPosition |> Network.pos_cmd
        #extract_data returns an XPCPosition instance, from which handle_data
        #constructs a position command, which propagates to sys.u.echo via
        #loopback, and finally to sys.y.echo within f_disc!
        @test sim.y.echo === cmd[1]

        return sim

    end

end


################################ JSON Loopback #################################


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
