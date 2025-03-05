module Network

using Sockets
using UnPack
using GLFW
using JSON3

using ..IODevices

export UDPOutput, UDPInput
export XP12Client, XP12Pose

#UDPInput and UDPOutput both use the EOT character as a shutdown request. This
#provides a means to prevent the UDPInput thread from getting stuck in the
#blocking recv call indefinitely. Any source providing data to the Simulation
#via UDPInput should send an EOT character before shutting down. This is done by
#UDPOutput, which avoids issues during loopback tests.

################################################################################
################################# UDInput ######################################

@kwdef mutable struct UDPInput{A <: IPAddr} <: InputDevice
    socket::UDPSocket = UDPSocket()
    address::A = IPv4("127.0.0.1")#IP address we'll be listening at
    port::Int = 49017 #port we'll be listening at
    should_close::Bool = false
    function UDPInput(socket::UDPSocket, address::A, port::Integer,
                      should_close::Bool) where {A <: IPAddr}
        new{A}(socket, address, port, should_close)
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

@kwdef mutable struct UDPOutput{A <: IPAddr} <: OutputDevice
    socket::UDPSocket = UDPSocket()
    address::A = IPv4("127.0.0.1")#IP address we'll be listening at
    port::Int = 49017 #port we'll be listening at
    function UDPOutput(socket::UDPSocket, address::A, port::Integer) where {A <: IPAddr}
        new{A}(socket, address, port)
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
################################# XP12Client ###################################

@kwdef struct XP12Pose
    aircraft::Int32 = 0 #aircraft number
    ϕ::Float64 = 47.80433 #degrees
    λ::Float64 = 12.997 #degrees
    h::Float64 = 429.0 #meters
    ψ::Float32 = 157.0 #degrees
    θ::Float32 = 3.7 #degrees
    φ::Float32 = -0.5 #degrees
end

struct XP12Client{U <: UDPOutput} <: OutputDevice
    udp::U
    function XP12Client(udp::U) where {U <: UDPOutput}
        new{U}(udp)
    end
end

function XP12Client(; address = IPv4("127.0.0.1"), port = 49000, kwargs...)
    XP12Client(UDPOutput(; address, port, kwargs...))
end

function IODevices.init!(xpc::XP12Client)
    IODevices.init!(xpc.udp)
    #disable pose updates for aircraft 0
    IODevices.handle_data!(xpc.udp, msg_set_dref(
        "sim/operation/override/override_planepath[0]", 1))
    #override control surface positions (not only visuals!)
    IODevices.handle_data!(xpc.udp, msg_set_dref(
        "sim/operation/override/override_control_surfaces", 1))
end

IODevices.shutdown!(xpc::XP12Client) = IODevices.shutdown!(xpc.udp)

function IODevices.handle_data!(xpc::XP12Client, data::String)
    #give X-Plane some breating room (limit the update rate to 100Hz)
    sleep(0.01)
    IODevices.handle_data!(xpc.udp, data)
    # IODevices.handle_data!(xpc.udp, msg_set_dref(
    #     "sim/flightmodel2/wing/rudder1_deg[10]", 30))
end

############################ XPC Command Messages ##############################

#construct the UDP message to write a scalar or vector value to an arbitrary DREF
function msg_set_dref(dref_id::AbstractString, value::Real)
# function msg_set_dref()

    #ascii() ensures ASCII data, codeunits returns a CodeUnits object, which
    #behaves similarly to a byte array. this is equivalent to b"text".
    #casting to Vector{UInt8} would also work

    #length(dref_id) returns the number of characters in dref_id, not its actual
    #length in bytes. both of these are equal only if dref_id is pure ASCII. to
    #pad the message, we need the actual length in bytes. therefore, we either
    #do ascii(dref_id) (which ensures only ascii characters are present)
    #or accept non-ascii characters and use length(dref_id |> codeunits)
    # dref_id = "sim/operation/override/override_planepath"
    # value = 1

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        Float32(value),
        ascii(dref_id) * "\0", #zero-terminated dref id pure ASCII string
        zeros(UInt8, 499-length(dref_id)) #pad the message to 509 bytes
        )

    return String(take!(buffer))
end

function msg_cmd(cmd_id::AbstractString)

    buffer = IOBuffer()
    write(buffer,
        b"CMND\0",
        ascii(cmd_id) * "\0", #zero-terminated command id pure ASCII string
        zeros(UInt8, 499-length(cmd_id))) #pad the message to 509 bytes

    return String(take!(buffer))
end

#construct the UDP message to set aircraft position and attitude
# function set_xp12pos_msg(pos::XPCPosition)
function msg_set_pose(pose::XP12Pose)

    @unpack aircraft, ϕ, λ, h, ψ, θ, φ = pose

    buffer = IOBuffer()
    write(buffer, b"VEHS\0", aircraft, ϕ, λ, h, ψ, θ, φ)
    # write(buffer, b"VEHX\0", aircraft, ϕ, λ, h, ψ, θ, φ)

    return String(take!(buffer))

end

end #module