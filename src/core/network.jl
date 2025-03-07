module Network

using Sockets
using UnPack
using GLFW
using JSON3

using ..IODevices

export UDPOutput, UDPInput
export XPlane12Output, XPlanePose

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

function IODevices.handle_data!(device::UDPOutput, data::NTuple{N, String}) where N
    foreach(data) do msg
        IODevices.handle_data!(device, msg)
    end
end


################################################################################
################################# XPlane12Output ###################################

struct XPlane12Output{U <: UDPOutput} <: OutputDevice
    udp::U
    function XPlane12Output(udp::U) where {U <: UDPOutput}
        new{U}(udp)
    end
end

function XPlane12Output(; address = IPv4("127.0.0.1"), port = 49000, kwargs...)
    XPlane12Output(UDPOutput(; address, port, kwargs...))
end

function IODevices.init!(xpc::XPlane12Output)

    override_pose = "sim/operation/override/override_planepath[0]"
    override_surf = "sim/operation/override/override_control_surfaces[0]"
    override_prop = "sim/flightmodel2/engines/prop_disc/override[0]"
    override_nws = "sim/operation/override/override_wheel_steer[0]"

    IODevices.init!(xpc.udp)
    msg_tuple = (
        Network.xpmsg_set_dref(override_pose, 1),
        Network.xpmsg_set_dref(override_surf, 1),
        Network.xpmsg_set_dref(override_prop, 1),
        Network.xpmsg_set_dref(override_nws, 1),
    )
    IODevices.handle_data!(xpc.udp, msg_tuple)
end

IODevices.shutdown!(xpc::XPlane12Output) = IODevices.shutdown!(xpc.udp)

function IODevices.handle_data!(xpc::XPlane12Output, data::Union{String, NTuple{N, String}}) where N
    sleep(0.01) #give X-Plane some breathing room
    IODevices.handle_data!(xpc.udp, data)
end

############################### XPlane Messages ################################

@kwdef struct XPlanePose
    ϕ::Float64 = 47.80433 #degrees
    λ::Float64 = 12.997 #degrees
    h::Float64 = 429.0 #meters
    ψ::Float32 = 157.0 #degrees
    θ::Float32 = 3.7 #degrees
    φ::Float32 = -0.5 #degrees
end


#note: length(s::String) returns the number of characters in s, not its actual
#length in bytes, which can be found as length(codeunits(s)) or
#length(Vector{UInt8}(s)). these are equal only if s is pure ASCII

function xpmsg_cmd(cmd_id::AbstractString)

    buffer = IOBuffer()
    write(buffer,
        b"CMND\0",
        ascii(cmd_id) * "\0", #zero-terminated pure ASCII string
        zeros(UInt8, 499-length(cmd_id))) #pad the message to 509 bytes

    return String(take!(buffer))
end

#construct the UDP message to write a scalar value to an arbitrary DREF
function xpmsg_set_dref(dref_id::AbstractString, value::Real)

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        Float32(value),
        ascii(dref_id) * "\0", #zero-terminated pure ASCII string
        zeros(UInt8, 499-length(dref_id)) #pad the message to 509 bytes
        )

    return String(take!(buffer))
end

#construct the UDP message to set aircraft position and attitude
function xpmsg_set_pose(pose::XPlanePose)

    @unpack ϕ, λ, h, ψ, θ, φ = pose
    aircraft = Int32(0) #aircraft index, we only support one for now

    buffer = IOBuffer()
    write(buffer, b"VEHS\0", aircraft, ϕ, λ, h, ψ, θ, φ)
    # write(buffer, b"VEHX\0", aircraft, ϕ, λ, h, ψ, θ, φ) #for X-Plane 11

    return String(take!(buffer))

end

end #module