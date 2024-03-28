module Networking

using Sockets
using UnPack
using GLFW
using JSON3

using ..Sim
using ..Systems
using ..IODevices

export UDPSender, UDPReceiver, XPCClient, UDPServer, UDPClient


################################################################################
################################## Output ######################################

mutable struct UDPSender{T <: IPAddr}
    socket::UDPSocket
    host::T
    port::Int
    function UDPSender(; host::T = IPv4("127.0.0.1"),
                         port::Integer = 49017) where {T <: IPAddr}
        new{T}(UDPSocket(), host, port)
    end
end

init!(sender::UDPSender) = (sender.socket = UDPSocket()) #new socket

function Sockets.send(sender::UDPSender, msg)
    @unpack socket, host, port = sender
    send(socket, host, port, msg)
end

shutdown!(::UDPSender) = nothing

################################################################################

struct UDPClient{F1 <: Function, F2 <: Function} <: OutputDevice
    init_callback::F1
    output_callback::F2
    sender::UDPSender
    function UDPClient(; init_callback::F1 = init_callback,
                        output_callback::F2 = output_callback,
                        kwargs...) where {F1 <: Function, F2 <: Function}
        new{F1, F2}(init_callback, output_callback, UDPSender(; kwargs...))
    end
end

function IODevices.init!(client::UDPClient)
    try
        init!(client.sender)
        data = client.init_callback()
        !isempty(data) && send(client.sender, data)
    catch ex
        st = stacktrace(catch_backtrace())
        @warn("Client init failed with $ex in $(st[1])")
    end
end

function IODevices.process_output!(client::UDPClient, output::Sim.Output)
    try
        data = client.output_callback(output)
        !isempty(data) && send(client.sender, data)
    catch ex
        st = stacktrace(catch_backtrace())
        @warn("Client output processing failed with $ex in $(st[1])")
    end
end

IODevices.shutdown!(client::UDPClient) = shutdown!(client.sender)
IODevices.should_close(::UDPClient) = false

init_callback()::Vector{UInt8} = UInt8[]
#output_callback signature; the output is the result of some user-specified
#serialization of simulation output (typically via JSON)
output_callback(::Sim.Output)::Vector{UInt8} = UInt8[]


############################### XPCClient ####################################

struct XPCClient <: OutputDevice
    sender::UDPSender
    function XPCClient(; host::IPAddr = IPv4("127.0.0.1"), port::Integer = 49009)
        new(UDPSender(; host, port))
    end
end

function IODevices.init!(xpc::XPCClient)
    init!(xpc.sender)
    disable_physics!(xpc)
end

 #to be overridden for each System's y
function IODevices.process_output!(xpc::XPCClient, out::Sim.Output)
    set_position!(xpc, out.y)
end

IODevices.should_close(xpc::XPCClient) = false #handled
IODevices.shutdown!(xpc::XPCClient) = shutdown!(xpc.sender)

function set_dref(xpc::XPCClient, dref_id::AbstractString, dref_value::Real)

    #ascii() ensures ASCII data, codeunits returns a CodeUnits object, which
    #behaves similarly to a byte array. this is equivalent to b"text".
    #Vector{UInt8}(dref_id) would also work
    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        UInt8(1),
        Float32(dref_value))

    send(xpc.sender, take!(buffer))
end

function set_dref(xpc::XPCClient, dref_id::AbstractString, dref_value::AbstractVector{<:Real})

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        dref_value |> length |> UInt8,
        Vector{Float32}(dref_value))

    send(xpc.sender, take!(buffer))
end

function disable_physics!(xpc::XPCClient)
    set_dref(xpc, "sim/operation/override/override_planepath", 1)
end

function set_position!(xpc::XPCClient; lat, lon, h_o, psi, theta, phi, aircraft::Integer = 0)

    #all angles must be provided in degrees
    buffer = IOBuffer()
    write(buffer,
        b"POSI\0", UInt8(aircraft),
        Float64(lat), Float64(lon), Float64(h_o),
        Float32(theta), Float32(phi), Float32(psi),
        Float32(-998)) #last one is landing gear (?!)

    send(xpc.sender, take!(buffer))

end



################################################################################
################################## Input ######################################

mutable struct UDPReceiver{T <: IPAddr}
    socket::UDPSocket
    host::T
    port::Int
    function UDPReceiver(; host::T = IPv4("127.0.0.1"),
                           port::Integer = 49017) where {T <: IPAddr}
        new{T}(UDPSocket(), host, port)
    end
end

function init!(receiver::UDPReceiver)
    receiver.socket = UDPSocket() #create a new socket on each initialization
    @unpack socket, host, port = receiver
    if !bind(socket, host, port; reuseaddr=true)
        @error( "Failed to bind socket to host $host, port $port")
    end
end

#apparently, there is currently no in place version of recv, so each call
#allocates a new Vector{UInt8}. this is not ideal, because it will continuously
#trigger the garbage collector
Sockets.recv(receiver::UDPReceiver) = recv(receiver.socket)
shutdown!(receiver::UDPReceiver) = close(receiver.socket)


############################## UDPServer #######################################

mutable struct UDPServer{F <: Function} <: InputDevice
    assign_callback!::F
    receiver::UDPReceiver
    buffer::IOBuffer
    function UDPServer(; assign_callback!::F = assign_callback!,
                        kwargs...) where {F <: Function}
        new{F}(assign_callback!, UDPReceiver(; kwargs...), IOBuffer())
    end
end

IODevices.init!(server::UDPServer) = init!(server.receiver)

function IODevices.update_input!(server::UDPServer)
    #will always write something, because recv blocks until something is
    #actually received
    write(server.buffer, recv(server.receiver))
end

function IODevices.assign_input!(sys::System, server::UDPServer, mapping::InputMapping)
    try
        #because update_input! blocks waiting for an UDP packet and this
        #function is called immediately after update_input!, we can count on
        #input_data not being empty
        input_data = take!(server.buffer)
        server.assign_callback!(sys, input_data, mapping)
    catch ex
        st = stacktrace(catch_backtrace())
        @warn("Server assign callback failed with $ex in $(st[1])")
    end
end

IODevices.shutdown!(server::UDPServer) = shutdown!(server.receiver)
IODevices.should_close(::UDPServer) = false

#server callback signature; the input string is deserialized and used to update
#the simulated System according to the specified InputMapping
assign_callback!(::System, ::Vector{UInt8}, ::InputMapping) = nothing



end #module