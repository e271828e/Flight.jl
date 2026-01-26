module TestSim

using Test, Logging, StructTypes, JSON3

using Flight.FlightCore
using Flight.FlightLib

export test_sim

function test_sim()
    @testset verbose = true "Sim" begin

        udp_loopback()
        xp12_loopback()
        json_loopback()

    end
end

################################################################################
############################### Simulation #####################################

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

function Modeling.f_step!(mdl::Model{FirstOrder})
    x_new = mdl.x[1] + 1
    # @info("Called f_step! at t = $(mdl.t[]) and x = $(mdl.x[1]), x updated to $(x_new)")
    # mdl.x .= x_new #if we want the change in x to propagate to y at the end of this step
end

function Modeling.f_periodic!(::NoScheduling, mdl::Model{FirstOrder})
    # println("Called f_periodic! at t = $(mdl.t[]), got y = $(mdl.y)")
end

Modeling.init!(mdl::Model{FirstOrder}, x0::Real = 0.0) = (mdl.x .= x0)


function test_sim_standalone()

    mdl = FirstOrder() |> Model
    sim = Simulation(mdl; dt = 0.1, Δt = 1.0, t_end = 5)
    x0 = 1.0
    init!(sim, x0)
    Sim.run!(sim)
    return sim

end



################################################################################
############################### IO Loopback ####################################

#input and output devices must not be mutually locking. otherwise, at least one
#of them may block irrecoverably when the simulation terminates. this coupling
#may happen for example if input and output share a loopback Channel and they
#make blocking put! and take! calls on it.

#to avoid this, at least one of them should only block when waiting on its
#SimInterface, but not on its external/loopback interface. this is the case with
#an UDP loopback, in which the UDPOutput may block when calling take! on the
#SimInterface Channel, but not on its send() call, which is nonblocking.

@kwdef struct TestSystem <: ModelDefinition end

@kwdef mutable struct TestSystemU
    input::Float64 = 0
end

@kwdef struct TestSystemY
    input::Float64 = 0
end

Modeling.U(::TestSystem) = TestSystemU()
Modeling.Y(::TestSystem) = TestSystemY()

@no_ode TestSystem
@no_step TestSystem

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:TestSystem})
    mdl.y = TestSystemY(; input = mdl.u.input)
end

function GUI.draw(mdl::Model{TestSystem}, label::String = "TestSystem")

    (; input) = mdl.y

    CImGui.Begin(label)

        CImGui.Text("input = $input")

    CImGui.End()

end #function


################################ UDP Loopback ##################################

struct UDPTestMapping <: IOMapping end

function IODevices.assign_input!(mdl::Model{TestSystem},
                            ::UDPTestMapping,
                            data::String)
    # @debug "Got $data"
    mdl.u.input = Vector{UInt8}(data)[1]
    # mdl.u.input = "Hi"
end

function IODevices.extract_output(::Model{TestSystem}, ::UDPTestMapping)
    data = UInt8[37] |> String
    # data = String([0x04]) #EOT character
    # @debug "Extracted $data"
    return data
end

function udp_loopback()

    @testset verbose = true "UDP Loopback" begin

        port = 14141
        mdl = TestSystem() |> Model
        sim = Simulation(mdl; t_end = 1.0)
        Sim.attach!(sim, UDPInput(; port), UDPTestMapping())
        Sim.attach!(sim, UDPOutput(; port), UDPTestMapping())

        #we need to run these paced, otherwise the IO can't keep up
        Sim.run!(sim; pace = 1)

        #mdl.y.output must have propagated to mdl.u.input via loopback, and then
        #to mdl.y.input within f_periodic!
        @test mdl.y.input == 37.0

        return sim

    end

end

################################ XPC Loopback ##################################

function IODevices.extract_output(::Model{TestSystem}, ::XPlane12ControlMapping)
    data = KinData() |> XPlanePose |> Network.xpmsg_set_pose
    return data
end

function xp12_loopback()

    @testset verbose = true "X-Plane 12 Loopback" begin

        port = 14143
        mdl = TestSystem() |> Model
        sim = Simulation(mdl; t_end = 1.0)
        Sim.attach!(sim, UDPInput(; port), UDPTestMapping())
        Sim.attach!(sim, XPlane12Control(; port))

        #we need to run these paced, otherwise the IO can't keep up
        Sim.run!(sim; pace = 1)

        cmd = KinData() |> XPlanePose |> Network.xpmsg_set_pose
        #extract_output returns an XPlanePose instance, from which extract_output
        #constructs a pose command message, which is sent through UDP by
        #handle_data! and finally reaches assign_input! via loopback. the first
        #character is converted to Float64 and assigned to mdl.u.input, and it
        #finally propagates to mdl.y.input within f_periodic!
        @test mdl.y.input === Float64(cmd[1])

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

function IODevices.extract_output(::Model{TestSystem}, ::JSONTestMapping)
    data = (input = 37.0,) |> JSON3.write
    # @info "Extracted $data"
    return data
end

function IODevices.assign_input!(mdl::Model{TestSystem},
                            ::JSONTestMapping,
                            data::String)

    # @info "Got $data"
    JSON3.read!(data, mdl.u)
    # @info "Echo is now $(mdl.u.input)"
end

function json_loopback()

    @testset verbose = true "JSON Loopback" begin

        port = 14142
        mdl = TestSystem() |> Model
        sim = Simulation(mdl; t_end = 1.0)
        Sim.attach!(sim, UDPInput(; port), JSONTestMapping())
        Sim.attach!(sim, UDPOutput(; port), JSONTestMapping())

        #trigger method precompilation
        JSON3.read!(JSON3.write((input = 0.0,)), mdl.u)

        #we need to run these paced, otherwise the IO can't keep up
        Sim.run!(sim; pace = 1)

        @test mdl.y.input == 37.0

        return sim

    end

end


################################## Joystick ####################################

function IODevices.assign_input!(mdl::Model{TestSystem},
                            ::IOMapping,
                            data::Joysticks.T16000MData)
    mdl.u.input = data.axes.stick_x
end

function joystick_input()

    @testset verbose = true "Joystick Input" begin

        mdl = TestSystem() |> Model
        sim = Simulation(mdl; t_end = 10.0)
        joystick = update_connected_joysticks()[1]
        Sim.attach!(sim, joystick)

        Sim.run!(sim; gui = true)

    end

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
