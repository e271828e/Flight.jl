module DemoUDP

using Sockets, UnPack
using StructTypes, JSON3

using Flight
using Flight.FlightCore.Sim
using Flight.FlightCore.IODevices
using Flight.FlightComponents.Control.Discrete: PID, PIDInput

export demo_udp

#declare PIDInput as mutable for JSON3 to read into it
StructTypes.StructType(::Type{PIDInput}) = StructTypes.Mutable()

################################################################################
############################### JSON #################################

@kwdef mutable struct JSONReceiver <: InputDevice
    socket::UDPSocket = UDPSocket()
    host::IPv4 = IPv4("127.0.0.1")
    port::Int64 = 49017
    buffer::Vector{UInt8} = zeros(UInt8, 2^16)
end

function IODevices.init!(receiver::JSONReceiver)
    receiver.socket = UDPSocket() #create a new socket on each execution

    @unpack socket, host, port = receiver
    if !bind(socket, host, port; reuseaddr=true)
        @error( "Failed to bind socket to host $host, port $port")
    end
end

function IODevices.assign!(sys::System{PID}, receiver::JSONReceiver, ::DefaultMapping)
    data_str = String(receiver.buffer)
    try
        # JSON3.read(data_str).input
        JSON3.read!(data_str, sys.u)
    catch ex
        println(ex)
    end
end

function IODevices.update!(receiver::JSONReceiver)
    #apparently, there is currently no in place version of recv, so each call
    #allocates a new Vector{UInt8}. this is not ideal, because it will
    #continuously trigger the garbage collector
    receiver.buffer = recv(receiver.socket) #option 1: discard previous contents
    # append!(receiver.buffer, recv(receiver.socket)) #option 2: queue
end

IODevices.shutdown!(receiver::JSONReceiver) = close(receiver.socket)

IODevices.should_close(::JSONReceiver) = false

################################################################################

@kwdef mutable struct JSONSender <: OutputDevice
    socket::UDPSocket = UDPSocket()
    host::IPv4 = IPv4("127.0.0.1")
    port::Int64 = 49017
end

function IODevices.init!(sender::JSONSender)
    sender.socket = UDPSocket() #create a new socket on each execution
end

function IODevices.update!(sender::JSONSender, data::Any)
    @unpack socket, host, port = sender

    send(socket, host, port, """{"input": $(data.t)}""")

    #JSON3 parses only the first JSON message in the string. probably it starts
    #parsing the string and when it finds a closing brace, it discards the rest

    # data_str_first = """{"input": 0.1}"""
    # data_str_second = """{"input": 0.2}""" send(socket, host, port,
    # data_str_first) send(socket, host, port, data_str_second) send(socket,
    # host, port, data_str_first*data_str_second)

end

IODevices.shutdown!(::JSONSender) = nothing
IODevices.should_close(::JSONSender) = false

################################################################################
################################################################################

function demo_udp()

    sys = PID() |> System
    sim = Simulation(sys; t_end = 10, dt = 0.02)

    #trigger compilation of parsing methods before launching the simulation
    JSON3.read!(JSON3.write(sys.u, allow_inf=true), sys.u; allow_inf=true)

    Sim.attach!(sim, JSONSender())
    Sim.attach!(sim, JSONReceiver())
    Sim.run_paced!(sim)

end

end #module