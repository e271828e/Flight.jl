module Network

using Sockets
using UnPack
using GLFW
using JSON3

using ..IODevices

export UDPOutput, UDPInput
export XPCClient, XPCPosition

#UDPInput and UDPOutput both use the EOT character as a shutdown request. This
#provides a means to prevent the UDPInput thread from getting stuck in the
#blocking recv call indefinitely. Any source providing data to the Simulation
#via UDPInput should send an EOT character before shutting down. This is done by
#UDPOutput, which avoids issues during loopback tests.

################################################################################
################################# UDInput ######################################

@kwdef mutable struct UDPInput{T, A <: IPAddr} <: InputDevice{T}
    socket::UDPSocket = UDPSocket()
    address::A = IPv4("127.0.0.1")#IP address we'll be listening at
    port::Int = 49017 #port we'll be listening at
    should_close::Bool = false
    function UDPInput(socket::UDPSocket, address::A, port::Integer,
                      should_close::Bool) where {A <: IPAddr}
        new{String, A}(socket, address, port, should_close)
    end
end

function IODevices.init!(device::UDPInput)
    device.socket = UDPSocket() #create a new socket on each initialization
    @unpack socket, address, port = device
    if !bind(socket, address, port; reuseaddr=true)
        @error( "Failed to bind socket to address $address, port $port")
    end
end

IODevices.should_close(device::UDPInput) = device.should_close
IODevices.shutdown!(device::UDPInput) = close(device.socket)

function IODevices.get_data!(device::UDPInput)
    data = recv(device.socket) |> String
    (data === "\x04") && (device.should_close = true) #received EOT character
    return data
end


################################################################################
################################# UDPOutput ####################################

@kwdef mutable struct UDPOutput{T, A <: IPAddr} <: OutputDevice{T}
    socket::UDPSocket = UDPSocket()
    address::A = IPv4("127.0.0.1")#IP address we'll be listening at
    port::Int = 49017 #port we'll be listening at
    function UDPOutput(socket::UDPSocket, address::A, port::Integer) where {A <: IPAddr}
        new{String, A}(socket, address, port)
    end
end

function IODevices.init!(device::UDPOutput)
    device.socket = UDPSocket() #get a new socket on each initialization
end

function IODevices.shutdown!(device::UDPOutput)
    IODevices.handle_data!(device, "\x04") #send EOT character
    close(device.socket)
end

function IODevices.handle_data!(device::UDPOutput, data::String)
    @unpack socket, address, port = device
    # @info "UDPOutput: Sending $(length(data)) bytes"
    # @info "UDPOutput: Sending $data"
    !isempty(data) && send(socket, address, port, data)
end


################################################################################
################################# XPCClient ####################################

@kwdef struct XPCPosition
    ϕ::Float64 = 0.0 #degrees
    λ::Float64 = 0.0 #degrees
    h::Float64 = 0.0 #meters
    ψ::Float32 = 0.0 #degrees
    θ::Float32 = 0.0 #degrees
    φ::Float32 = 0.0 #degrees
    aircraft::UInt8 = 0 #aircraft number
end

struct XPCClient{T, U <: UDPOutput} <: OutputDevice{T}
    udp::U
    function XPCClient(udp::U) where {U <: UDPOutput}
        new{XPCPosition, U}(udp)
    end
end

XPCClient(args...; kwargs...) = XPCClient(UDPOutput(args...; kwargs...))

#disable X-Plane physics
function IODevices.init!(xpc::XPCClient)
    IODevices.init!(xpc.udp)
    IODevices.handle_data!(xpc.udp, dref_cmd(
        "sim/operation/override/override_planepath", 1))
end

IODevices.shutdown!(xpc::XPCClient) = IODevices.shutdown!(xpc.udp)

function IODevices.handle_data!(xpc::XPCClient, data::XPCPosition)
    IODevices.handle_data!(xpc.udp, pos_cmd(data))
end

############################ XPC Command Messages ##############################

#write a scalar or vector value to an arbitrary DREF
function dref_cmd(id::AbstractString, value::Union{Real, AbstractVector{<:Real}})

    #ascii() ensures ASCII data, codeunits returns a CodeUnits object, which
    #behaves similarly to a byte array. this is equivalent to b"text".
    #casting to Vector{UInt8} would also work
    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        id |> length |> UInt8,
        id |> ascii |> codeunits,
        value |> length |> UInt8,
        Float32.(value))

    return String(take!(buffer))
end

#set aircraft position and attitude
function pos_cmd(pos::XPCPosition)

    @unpack ϕ, λ, h, ψ, θ, φ, aircraft = pos

    buffer = IOBuffer(sizehint = 64)
    write(buffer, b"POSI\0", aircraft, ϕ, λ, h, θ, φ, ψ, Float32(-998))

    return String(take!(buffer))

end

end #module