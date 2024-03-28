module Networking

using Sockets
using UnPack
using GLFW
using JSON3

using ..Sim
using ..Systems
using ..IODevices

export UDPSender, UDPReceiver, XPCDevice, UDPServer, UDPClient


################################################################################
################################### UDP ########################################

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


################################################################################

struct UDPClient{F <: Function} <: OutputDevice
    callback::F
    sender::UDPSender
    function UDPClient(; callback::F = client_callback,
                        kwargs...) where {F <: Function}
        new{F}(callback, UDPSender(; kwargs...))
    end
end

IODevices.init!(client::UDPClient) = init!(client.sender)

function IODevices.update!(client::UDPClient, data::Sim.Output)
    try
        str = client.callback(data)
        send(client.sender, str)
    catch ex
        st = stacktrace(catch_backtrace())
        @warn("Client update failed with $ex in $(st[1])")
    end
end

IODevices.shutdown!(client::UDPClient) = shutdown!(client.sender)
IODevices.should_close(::UDPClient) = false

#client callback signature; the output string is the result of some
#user-specified serialization of simulation output (typically via JSON)
client_callback(::Sim.Output)::String = ""


################################################################################

mutable struct UDPServer{F <: Function} <: InputDevice
    callback::F
    receiver::UDPReceiver
    buffer::Vector{UInt8}
    function UDPServer(; callback::F = server_callback!,
                        buffer_size::Integer = 1024,
                        kwargs...) where {F <: Function}
        new{F}(callback, UDPReceiver(; kwargs...), zeros(UInt8, buffer_size))
    end
end

IODevices.init!(server::UDPServer) = init!(server.receiver)

#discards previous contents
IODevices.update!(server::UDPServer) = (server.buffer = recv(server.receiver))

#this one queues new data
# IODevices.update!(server::UDPServer) = (append!(server.buffer, recv(server.socket)))

function IODevices.assign!(sys::System, server::UDPServer, mapping::InputMapping)
    try
        server.callback(sys, String(server.buffer), mapping)
    catch ex
        st = stacktrace(catch_backtrace())
        @warn("Server callback failed with $ex in $(st[1])")
    end
end

IODevices.shutdown!(server::UDPServer) = shutdown!(server.receiver)
IODevices.should_close(::UDPServer) = false

#server callback signature; the input string is deserialized and used to update
#the simulated System according to the specified InputMapping
server_callback!(::System, ::String, ::InputMapping) = nothing


################################################################################
############################### XPCDevice ####################################

struct XPCDevice <: OutputDevice
    sender::UDPSender
    function XPCDevice(; host::IPAddr = IPv4("127.0.0.1"), port::Integer = 49009)
        new(UDPSender(; host, port))
    end
end

function set_dref(xpc::XPCDevice, dref_id::AbstractString, dref_value::Real)

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

    send(xpc.sender, buffer.data)
end

function set_dref(xpc::XPCDevice, dref_id::AbstractString, dref_value::AbstractVector{<:Real})

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        dref_value |> length |> UInt8,
        Vector{Float32}(dref_value))

    send(xpc.sender, buffer.data)
end

function disable_physics!(xpc::XPCDevice)
    set_dref(xpc, "sim/operation/override/override_planepath", 1)
end

function set_position!(xpc::XPCDevice; lat, lon, h_o, psi, theta, phi, aircraft::Integer = 0)

    #all angles must be provided in degrees
    buffer = IOBuffer()
    write(buffer,
        b"POSI\0", UInt8(aircraft),
        Float64(lat), Float64(lon), Float64(h_o),
        Float32(theta), Float32(phi), Float32(psi),
        Float32(-998)) #last one is landing gear (?!)

    send(xpc.sender, buffer.data)

end

function IODevices.init!(xpc::XPCDevice)
    init!(xpc.sender)
    disable_physics!(xpc)
end

 #to be overridden for each System's y
IODevices.update!(xpc::XPCDevice, out::Sim.Output) = set_position!(xpc, out.y)

IODevices.should_close(xpc::XPCDevice) = false #handled
IODevices.shutdown!(xpc::XPCDevice) = shutdown!(xpc.sender)

end #module