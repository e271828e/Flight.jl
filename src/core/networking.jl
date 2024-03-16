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
    port::Int64 = 49016
end

function IODevices.init!(receiver::DummyUDPReceiver)
    receiver.socket = UDPSocket() #create a new socket on each execution

    @unpack socket, host, port = receiver
    if !bind(socket, host, port; reuseaddr=true)
        @error( "Failed to bind socket to host $host, port $port")
    end
    # println("DummyUDPReceiver initialized")
end

IODevices.assign!(::Any, ::DummyUDPReceiver, ::DefaultMapping) = nothing

function IODevices.update!(receiver::DummyUDPReceiver)
    data = recv(receiver.socket)
    # println("DummyUDPReceiver got this message: $data")
end

function IODevices.shutdown!(receiver::DummyUDPReceiver)
    # println("DummyUDPReceiver shutting down")
    close(receiver.socket)
end

IODevices.should_close(::DummyUDPReceiver) = false

@kwdef mutable struct DummyUDPSender <: OutputDevice
    socket::UDPSocket = UDPSocket()
    host::IPv4 = IPv4("127.0.0.1")
    port::Int64 = 49016
end

function IODevices.init!(sender::DummyUDPSender)
    sender.socket = UDPSocket() #create a new socket on each execution
    # println("DummyUDPSender initialized")
end

function IODevices.update!(sender::DummyUDPSender, ::Any)
    @unpack socket, host, port = sender
    send(socket, host, port, "Hello")
end

function IODevices.shutdown!(::DummyUDPSender)
    # println("DummyUDPSender shutting down")
end

IODevices.should_close(::DummyUDPSender) = false


################################################################################
############################### XPCDevice ####################################

mutable struct XPCDevice <: OutputDevice
    socket::UDPSocket
    host::IPv4
    port::Int64
    update_interval::Int64
    window::GLFW.Window
    function XPCDevice(; socket = UDPSocket(), host = IPv4("127.0.0.1"), port = 49009, update_interval = 1)
        new(socket, host, port, update_interval) #window uninitialized
    end
end

function set_dref(xp::XPCDevice, dref_id::AbstractString, dref_value::Real)

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

    send(xp.socket, xp.host, xp.port, buffer.data)
end

function set_dref(xp::XPCDevice, dref_id::AbstractString, dref_value::AbstractVector{<:Real})

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        dref_value |> length |> UInt8,
        Vector{Float32}(dref_value))

    send(xp.socket, xp.host, xp.port, buffer.data)
end

disable_physics!(xp::XPCDevice) = set_dref(xp, "sim/operation/override/override_planepath", 1)

function set_position!(xp::XPCDevice; lat, lon, h_o, psi, theta, phi, aircraft::Integer = 0)

    #all angles must be provided in degrees
    buffer = IOBuffer()
    write(buffer,
        b"POSI\0", UInt8(aircraft),
        Float64(lat), Float64(lon), Float64(h_o),
        Float32(theta), Float32(phi), Float32(psi),
        Float32(-998)) #last one is landing gear (?!)

    send(xp.socket, xp.host, xp.port, buffer.data)

end
########################### IODevices extensions ###############################

function IODevices.init!(xp::XPCDevice)
    xp.window = GLFW.CreateWindow(640, 480, "XPCDevice")
    @unpack window, update_interval = xp
    # GLFW.HideWindow(window)
    GLFW.MakeContextCurrent(window)
    GLFW.SwapInterval(update_interval)
    disable_physics!(xp)
end

function IODevices.update!(xp::XPCDevice, out::Sim.Output)
    # GLFW.SwapBuffers(xp.window) #honor the requested update_interval
    set_position!(xp, out.y) #to be overridden for each System's y
    GLFW.PollEvents() #see if we got a shutdown request
end

IODevices.should_close(xp::XPCDevice) = GLFW.WindowShouldClose(xp.window)
IODevices.shutdown!(xp::XPCDevice) = GLFW.DestroyWindow(xp.window)

end #module