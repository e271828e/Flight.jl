module Networking

using Sockets
using UnPack
using GLFW

using ..Sim
using ..IODevices

export XPCDevice


################################################################################
############################### UDPDummies #################################

@kwdef mutable struct DummyUDPReceiver <: InputDevice
    socket::UDPSocket = UDPSocket()
    host::IPv4 = IPv4("127.0.0.1")
    port::Int64 = 49017
end

function IODevices.init!(receiver::DummyUDPReceiver)
    receiver.socket = UDPSocket() #create a new socket on each execution

    @unpack socket, host, port = receiver
    if !bind(socket, host, port; reuseaddr=true)
        @error( "Failed to bind socket to host $host, port $port")
    end
end

IODevices.assign!(::Any, ::DummyUDPReceiver, ::DefaultMapping) = nothing

function IODevices.update!(receiver::DummyUDPReceiver)
    data = recv(receiver.socket)
    # println("DummyUDPReceiver got this message: $data")
    println("""{"eng_start": true}""")
end

IODevices.shutdown!(receiver::DummyUDPReceiver) = close(receiver.socket)

IODevices.should_close(::DummyUDPReceiver) = false

################################################################################

@kwdef mutable struct DummyUDPSender <: OutputDevice
    socket::UDPSocket = UDPSocket()
    host::IPv4 = IPv4("127.0.0.1")
    port::Int64 = 49017
end

function IODevices.init!(sender::DummyUDPSender)
    sender.socket = UDPSocket() #create a new socket on each execution
end

function IODevices.update!(sender::DummyUDPSender, ::Any)
    @unpack socket, host, port = sender
    send(socket, host, port, "Hello")
end

IODevices.shutdown!(::DummyUDPSender) = nothing
IODevices.should_close(::DummyUDPSender) = false


################################################################################
############################### XPCDevice ####################################

@kwdef mutable struct XPCDevice <: OutputDevice
    socket::UDPSocket = UDPSocket()
    host::IPv4 = IPv4("127.0.0.1")
    port::Int64 = 49009
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

    send(xpc.socket, xpc.host, xpc.port, buffer.data)
end

function set_dref(xpc::XPCDevice, dref_id::AbstractString, dref_value::AbstractVector{<:Real})

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        dref_value |> length |> UInt8,
        Vector{Float32}(dref_value))

    send(xpc.socket, xpc.host, xpc.port, buffer.data)
end

disable_physics!(xpc::XPCDevice) = set_dref(xpc, "sim/operation/override/override_planepath", 1)

function set_position!(xpc::XPCDevice; lat, lon, h_o, psi, theta, phi, aircraft::Integer = 0)

    #all angles must be provided in degrees
    buffer = IOBuffer()
    write(buffer,
        b"POSI\0", UInt8(aircraft),
        Float64(lat), Float64(lon), Float64(h_o),
        Float32(theta), Float32(phi), Float32(psi),
        Float32(-998)) #last one is landing gear (?!)

    send(xpc.socket, xpc.host, xpc.port, buffer.data)

end
########################### IODevices extensions ###############################

function IODevices.init!(xpc::XPCDevice)
    xpc.socket = UDPSocket() #create a new socket on each execution
    disable_physics!(xpc)
end

 #to be overridden for each System's y
IODevices.update!(xpc::XPCDevice, out::Sim.Output) = set_position!(xpc, out.y)

IODevices.should_close(xpc::XPCDevice) = false #handled
IODevices.shutdown!(::XPCDevice) = nothing

end #module