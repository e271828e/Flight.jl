module JSONSketch1

using Sockets, UnPack
using StructTypes, JSON3

using Flight
using Flight.FlightCore.Sim
using Flight.FlightCore.IODevices
using Flight.FlightCore.Networking
using Flight.FlightComponents.Control.Discrete: PID, PIDInput

export json_sketch1

#note: JSON3 parses only the first JSON message in a string. apparently, it
#starts parsing and once it finds the closing brace, it discards the rest

#declare PIDInput as mutable so that JSON3 can read into it
StructTypes.StructType(::Type{PIDInput}) = StructTypes.Mutable()

#use simulation time to generate a sinusoidal signal and put its current value
#into a JSON string. this will be read by the server into the input field of the
#simulated System's input, which is a PIDInput struct
function output_callback(sim_out::Sim.Output)::Vector{UInt8}
    freq = 0.5
    value = sin(2π*freq*sim_out.t)
    return Vector{UInt8}("""{"input": $value}""")
end

#once PIDInput has been declared, JSON3 can automatically read a JSON string
#into one or more of its fields
function assign_callback!(sys::System{PID}, data::Vector{UInt8}, ::InputMapping)
    JSON3.read!(String(data), sys.u)
end

function json_sketch1()

    sys = PID() |> System
    sim = Simulation(sys; t_end = 10, dt = 0.02)

    #trigger compilation of parsing methods before launching the simulation
    JSON3.read!(JSON3.write(sys.u, allow_inf=true), sys.u; allow_inf=true)

    Sim.attach!(sim, UDPClient(; port = 49017, output_callback))
    Sim.attach!(sim, UDPServer(; port = 49017, assign_callback!))
    Sim.run_paced!(sim)

end

end #module