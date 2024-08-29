module Network

using Sockets
using UnPack
using GLFW
using JSON3

using ..IODevices

export UDPOutput, UDPInput, XPCClient


################################################################################
################################# UDPOutput ####################################

@kwdef mutable struct UDPOutput{T <: IPAddr} <: OutputDevice
    socket::UDPSocket = UDPSocket()
    address::T = IPv4("127.0.0.1") #IP address we'll be sending to
    port::Int = 49017 #port we'll be sending to
end

function IODevices.init!(output::UDPOutput)
    output.socket = UDPSocket() #get a new socket on each initialization
end

function IODevices.handle_data(output::UDPOutput, data::Vector{UInt8})
    try
        !isempty(data) && send(output.socket, output.address, output.port, data)
    catch ex
        st = stacktrace(catch_backtrace())
        @warn("UDPOutput failed with $ex in $(st[1])")
    end
end

IODevices.shutdown!(output::UDPOutput) = close(output.socket)

Base.Channel(::UDPOutput, size::Int) = Channel{Vector{UInt8}}(size)


################################################################################
################################# UDInput ######################################

@kwdef mutable struct UDPInput{T <: IPAddr} <: InputDevice
    socket::UDPSocket = UDPSocket()
    address::T = IPv4("127.0.0.1") #IP address we'll be listening at
    port::Int = 49017 #port we'll be listening at
end

function IODevices.init!(input::UDPInput)
    input.socket = UDPSocket() #create a new socket on each initialization
    @unpack socket, address, port = input
    if !bind(socket, address, port; reuseaddr=true)
        @error( "Failed to bind socket to address $address, port $port")
    end
end

#will either block or write something to the buffer, because recv blocks until
#an UDP packet is actually received
function IODevices.get_data(input::UDPInput)
    #there is currently no in place version of recv, so each call allocates a
    #new Vector{UInt8}. not ideal
    data = recv(input.socket)
    @info "Got $data"
    return data
end

IODevices.shutdown!(input::UDPInput) = close(input.socket)

Base.Channel(::UDPInput, size::Int) = Channel{Vector{UInt8}}(size)

# ############################### XPCClient ####################################

# function XPCClient(; address::IPAddr = IPv4("127.0.0.1"), port::Integer = 49009)
#     UDPClient(; address, port,
#                 init_callback = xpc_init_callback,
#                 output_callback = xpc_output_callback)
# end

# #to be overridden for each System's y
# xpc_output_callback(out::Sim.SimData) = set_xpc_pos!(out.y)
# #disables physics
# xpc_init_callback() = set_xpc_dref("sim/operation/override/override_planepath", 1)

# function set_xpc_dref(dref_id::AbstractString, dref_value::Real)

#     #ascii() ensures ASCII data, codeunits returns a CodeUnits object, which
#     #behaves similarly to a byte array. this is equivalent to b"text".
#     #Vector{UInt8}(dref_id) would also work
#     buffer = IOBuffer()
#     write(buffer,
#         b"DREF\0",
#         dref_id |> length |> UInt8,
#         dref_id |> ascii |> codeunits,
#         UInt8(1),
#         Float32(dref_value))

#     return take!(buffer)
# end

# function set_xpc_dref(dref_id::AbstractString, dref_value::AbstractVector{<:Real})

#     buffer = IOBuffer()
#     write(buffer,
#         b"DREF\0",
#         dref_id |> length |> UInt8,
#         dref_id |> ascii |> codeunits,
#         dref_value |> length |> UInt8,
#         Vector{Float32}(dref_value))

#     return take!(buffer)
# end

# #this needs to be extended in AircraftBase
# function set_xpc_pos!(; lat, lon, h_o, psi, theta, phi, aircraft::Integer = 0)

#     #all angles expected in degrees
#     buffer = IOBuffer()
#     write(buffer,
#         b"POSI\0", UInt8(aircraft),
#         Float64(lat), Float64(lon), Float64(h_o),
#         Float32(theta), Float32(phi), Float32(psi),
#         Float32(-998)) #last one is landing gear (?!)

#     return take!(buffer)

# end




end #module