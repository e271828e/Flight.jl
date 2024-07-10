module TestNetworking

using Test
using Sockets, UnPack
using StructTypes, JSON3

using Flight
using Flight.FlightCore.Sim
using Flight.FlightCore.IODevices
using Flight.FlightCore.Networking

export test_json_loopback1

@kwdef struct TestSystem <: SystemDefinition end

@kwdef mutable struct TestSystemU
    input::Float64 = 0.0
end

Systems.U(::TestSystem) = TestSystemU()
Systems.Y(::TestSystem) = 0.0

function Systems.f_disc!(sys::System{TestSystem}, ::Real)
    sys.y = sys.u.input
end

#declare TestSystemU as mutable so that JSON3 can read into it
StructTypes.StructType(::Type{TestSystemU}) = StructTypes.Mutable()

function test_json_loopback1()

    #use simulation time to generate a sinusoidal signal and put its current
    #value into a JSON string. this will be read by the server into the input
    #field of the simulated System's input, which is a TestSystemU struct. note:
    #JSON3 parses only the first JSON message in a string. apparently, it
    #starts parsing and once it finds the closing brace, it discards the rest
    function output_callback(::Sim.Output)::Vector{UInt8}
        value = 2.4
        return Vector{UInt8}("""{"input": $value}""")
    end

    #once TestSystemU has been declared, JSON3 can automatically read a JSON
    #string into one or more of its fields
    function assign_callback!(sys::System{TestSystem}, data::Vector{UInt8}, ::InputMapping)
        JSON3.read!(String(data), sys.u)
        # println(sys.u.input)
        # println(sys.y)
    end

    sys = TestSystem() |> System
    dt = Δt = 0.01
    sim = Simulation(sys; t_end = 1, dt, Δt, save_start = true)

    #trigger compilation of parsing methods before launching the simulation
    JSON3.read!(JSON3.write(sys.u, allow_inf=true), sys.u; allow_inf=true)

    Sim.attach!(sim, UDPClient(; port = 49017, output_callback))
    Sim.attach!(sim, UDPServer(; port = 49017, assign_callback!))
    Sim.run_paced!(sim)

    @test get_data(TimeSeries(sim))[end] == 2.4

end

end #module