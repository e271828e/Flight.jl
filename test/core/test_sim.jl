module TestSim

using Test, Logging, StructTypes, JSON3

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


function test_sim_standalone()

    sys = FirstOrder() |> System
    sim = Simulation(sys; dt = 0.1, Δt = 1.0, t_end = 5)
    x0 = 1.0
    f_init! = (sys)->Systems.init!(sys, x0)
    reinit!(sim, f_init!)
    return sim
    # Sim.run!(sim)

end



################################################################################
################################# TestSystem ###################################

#for a loopback test, it's essential that the input and output devices are not
#mutually locking. otherwise, at least one of them will block irrecoverably when
#the simulation terminates. this coupling may happen for example if input and
#output share a loopback Channel and they make blocking put! and take! calls on
#it.

#to avoid this, it is enough that at least one of them can only block when
#waiting on its SimInterface, but not on the loopback interface. this is
#the case with an UDP loopback, in which the UDPOutput may block when calling
#take! on the SimInterface Channel, but not on its send() call, which is
#nonblocking.


@kwdef struct TestSystem <: SystemDefinition end

@kwdef mutable struct TestSystemU
    input::Float64 = 0
end

@kwdef struct TestSystemY
    input::Float64 = 0
end

Systems.U(::TestSystem) = TestSystemU()
Systems.Y(::TestSystem) = TestSystemY()

function Systems.f_disc!(sys::System{<:TestSystem}, ::Real)
    # sleep(0.01)
    sys.y = TestSystemY(; input = sys.u.input)
end

################################ UDP Loopback ##################################

struct UDPTestMapping <: IOMapping end

function Systems.assign_input!(sys::System{TestSystem}, data::Vector{UInt8},
                            ::UDPTestMapping)
    @debug "Got $data"
    sys.u.input = data[1]
end

function Systems.extract_output(::System{TestSystem},
                            ::Type{Vector{UInt8}}, ::UDPTestMapping)
    data = UInt8[37]
    @debug "Extracted $data"
    return data
end

function udp_loopback()

    @testset verbose = true "UDP Loopback" begin

        port = 14149
        sys = TestSystem() |> System
        sim = Simulation(sys; t_end = 1.0)
        Sim.attach!(sim, UDPInput(; port), UDPTestMapping())
        Sim.attach!(sim, UDPOutput(; port), UDPTestMapping())

        # return sim
        Sim.run!(sim)

        #sys.y.output must have propagated to sys.u.input via loopback, and then
        #to sys.y.input within f_disc!
        @test sim.y.input == 37.0

        return sim

    end

end

################################ XPC Loopback ##################################

function Systems.extract_output(::System{TestSystem}, ::Type{XPCPosition}, ::IOMapping)
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
        #extract_output returns an XPCPosition instance, from which handle_data
        #constructs a position command, which propagates to sys.u.input via
        #loopback, and finally to sys.y.input within f_disc!
        @test sim.y.input === Float64(cmd[1])

        return sim

    end

end


################################ JSON Loopback #################################

#declare TestSystemY as immutable for JSON3 parsing
StructTypes.StructType(::Type{TestSystemY}) = StructTypes.Struct()
StructTypes.excludes(::Type{TestSystemY}) = (:input,) #only extract :output

#declare TestSystemU as mutable so that JSON3 can read into it
StructTypes.StructType(::Type{TestSystemU}) = StructTypes.Mutable()

#this doesn't work for switching field values via loopback, because the
#inversion applies both to serializing and deserializing, so it cancels out

# StructTypes.names(::Type{TestSystemU}) = ((:input, :output), (:output,  :input))

struct JSONTestMapping <: IOMapping end

function Systems.extract_output(::System{TestSystem}, ::Type{Vector{UInt8}},
                            ::JSONTestMapping)
    data = (input = 37.0,) |> JSON3.write
    @debug "Extracted $data"
    return Vector{UInt8}(data)
end

function Systems.assign_input!(sys::System{TestSystem}, data::Vector{UInt8},
                            ::JSONTestMapping)
    JSON3.read!(String(data), sys.u)
    @debug "Echo is now $(sys.u.input)"
end

function json_loopback()

    @testset verbose = true "JSON Loopback" begin

        port = 14149
        sys = TestSystem() |> System
        sim = Simulation(sys; t_end = 1.0)
        Sim.attach!(sim, UDPInput(; port), JSONTestMapping())
        Sim.attach!(sim, UDPOutput(; port), JSONTestMapping())
        Sim.run!(sim)

        @test sim.y.input == 37.0

        return sim

    end

end


################################## Joystick ####################################

function Systems.assign_input!(sys::System{TestSystem},
                            data::Joysticks.T16000MData,
                            ::IOMapping)
    sys.u.input = get_axis_value(data, :stick_x)
    @debug "Input $(sys.u.input)"
end

function Systems.assign_input!(sys::System{TestSystem},
                            data::Joysticks.XBoxControllerData,
                            ::IOMapping)
    sys.u.input = get_axis_value(data, :left_stick_x)
    @debug "Input $(sys.u.input)"
end

function joystick_input()

    @testset verbose = true "Joystick Input" begin

        sys = TestSystem() |> System
        sim = Simulation(sys; t_end = 1.0)
        joystick = get_connected_joysticks()[1]
        Sim.attach!(sim, joystick)

        Sim.run_interactive!(sim)

        return sim

    end

end

#REPL:
# with_logger(ConsoleLogger(Logging.Debug)) do
#     sim = TestSim.joystick_input()
# end



end #module
