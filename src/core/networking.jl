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
################################# UDPClient ####################################

mutable struct UDPClient{T <: IPAddr, F1 <: Function, F2 <: Function} <: OutputDevice
    socket::UDPSocket
    address::T #IP address we'll be sending to
    port::Int #port we'll be sending to
    init_callback::F1
    output_callback::F2
    function UDPClient(; address::T = IPv4("127.0.0.1"),
                        port::Integer = 49017,
                        init_callback::F1 = init_callback,
                        output_callback::F2 = output_callback
                        ) where {T <: IPAddr, F1 <: Function, F2 <: Function}
        new{T, F1, F2}(UDPSocket(), address, port, init_callback, output_callback)
    end
end

function IODevices.init!(client::UDPClient)
    client.socket = UDPSocket() #get a new socket on each initialization
    try
        data = client.init_callback()
        !isempty(data) && send(client.socket, client.address, client.port, data)
    catch ex
        st = stacktrace(catch_backtrace())
        @warn("Client init failed with $ex in $(st[1])")
    end
end

function IODevices.process_output!(client::UDPClient, output::Sim.Output)
    try
        data = client.output_callback(output)
        !isempty(data) && send(client.socket, client.address, client.port, data)
    catch ex
        st = stacktrace(catch_backtrace())
        @warn("Client output processing failed with $ex in $(st[1])")
    end
end

IODevices.shutdown!(client::UDPClient) = close(client.socket)
IODevices.should_close(::UDPClient) = false

init_callback()::Vector{UInt8} = UInt8[]
#output_callback signature; the output is the result of some user-specified
#serialization of simulation output (typically via JSON)
output_callback(::Sim.Output)::Vector{UInt8} = UInt8[]


############################### XPCClient ####################################

function XPCClient(; address::IPAddr = IPv4("127.0.0.1"), port::Integer = 49009)
    UDPClient(; address, port,
                init_callback = xpc_init_callback,
                output_callback = xpc_output_callback)
end

#to be overridden for each System's y
xpc_output_callback(out::Sim.Output) = set_xpc_pos!(out.y)
#disables physics
xpc_init_callback() = set_xpc_dref("sim/operation/override/override_planepath", 1)

function set_xpc_dref(dref_id::AbstractString, dref_value::Real)

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

    return take!(buffer)
end

function set_xpc_dref(dref_id::AbstractString, dref_value::AbstractVector{<:Real})

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        dref_value |> length |> UInt8,
        Vector{Float32}(dref_value))

    return take!(buffer)
end

#this needs to be extended in AircraftBase
function set_xpc_pos!(; lat, lon, h_o, psi, theta, phi, aircraft::Integer = 0)

    #all angles expected in degrees
    buffer = IOBuffer()
    write(buffer,
        b"POSI\0", UInt8(aircraft),
        Float64(lat), Float64(lon), Float64(h_o),
        Float32(theta), Float32(phi), Float32(psi),
        Float32(-998)) #last one is landing gear (?!)

    return take!(buffer)

end


################################################################################
################################ UDPServer #####################################

mutable struct UDPServer{T <: IPAddr, F <: Function} <: InputDevice
    buffer::IOBuffer
    socket::UDPSocket
    address::T #IP address we'll be listening at
    port::Int #port we'll be listening at
    assign_callback!::F
    function UDPServer(; address::T = IPv4("127.0.0.1"),
                        port::Integer = 49017,
                        assign_callback!::F = assign_callback!,
                        ) where {T <: IPAddr, F <: Function}
        new{T, F}(IOBuffer(), UDPSocket(), address, port, assign_callback!)
    end
end

function IODevices.init!(server::UDPServer)
    server.socket = UDPSocket() #create a new socket on each initialization
    @unpack socket, address, port = server
    if !bind(socket, address, port; reuseaddr=true)
        @error( "Failed to bind socket to address $address, port $port")
    end
end

#will either block or write something to the buffer, because recv blocks until
#an UDP packet is actually received
function IODevices.update_input!(server::UDPServer)
    #there is currently no in place version of recv, so each call allocates a
    #new Vector{UInt8}. this is not ideal, because it will trigger the garbage
    #collector
    write(server.buffer, recv(server.socket))
end

function IODevices.assign_input!(sys::System, server::UDPServer, mapping::InputMapping)
    try
        #since update_input! blocks waiting for an UDP packet, and this is
        #called immediately afterwards, we can assume input_data is never empty
        data = take!(server.buffer)
        !isempty(data) && server.assign_callback!(sys, data, mapping)
    catch ex
        st = stacktrace(catch_backtrace())
        @warn("Server assign callback failed with $ex in $(st[1])")
    end
end

IODevices.shutdown!(server::UDPServer) = close(server.socket)
IODevices.should_close(::UDPServer) = false

#server callback signature; the input string is deserialized and used to update
#the simulated System according to the specified InputMapping
assign_callback!(::System, ::Vector{UInt8}, ::InputMapping) = nothing



end #module