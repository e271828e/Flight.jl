module XPlane

using Sockets
using UnPack
using GLFW

using ..Sim
using ..IODevices

export XPConnect


################################################################################
############################### XPConnect ####################################

mutable struct XPConnect <: IODevice
    socket::UDPSocket
    host::IPv4
    port::Int64
    update_interval::Int64
    window::GLFW.Window
    function XPConnect(; socket = UDPSocket(), host = IPv4("127.0.0.1"), port = 49009, update_interval = 1)
        new(socket, host, port, update_interval) #window uninitialized
    end
end

function set_dref(xp::XPConnect, dref_id::AbstractString, dref_value::Real)

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

function set_dref(xp::XPConnect, dref_id::AbstractString, dref_value::AbstractVector{<:Real})

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        dref_value |> length |> UInt8,
        Vector{Float32}(dref_value))

    send(xp.socket, xp.host, xp.port, buffer.data)
end

disable_physics!(xp::XPConnect) = set_dref(xp, "sim/operation/override/override_planepath", 1)

function set_position!(xp::XPConnect; lat = -998, lon = -998, h = -998,
                       psi = -998, theta = -998, phi = -998,
                       aircraft::Integer = 0)

    #all angles must be provided in degrees
    buffer = IOBuffer()
    write(buffer,
        b"POSI\0", UInt8(aircraft),
        Float64(lat), Float64(lon), Float64(h),
        Float32(theta), Float32(phi), Float32(psi),
        Float32(-998)) #last one is landing gear (?!)

    send(xp.socket, xp.host, xp.port, buffer.data)

end
########################### IODevices extensions ###############################

function IODevices.init!(xp::XPConnect)
    xp.window = GLFW.CreateWindow(640, 480, "XPConnect")
    @unpack window, update_interval = xp
    # GLFW.HideWindow(window)
    GLFW.MakeContextCurrent(window)
    GLFW.SwapInterval(update_interval)
    disable_physics!(xp)
end

function IODevices.update!(xp::XPConnect, out::Sim.Output)
    GLFW.SwapBuffers(xp.window) #honor the requested update_interval
    set_position!(xp, out.y) #to be overridden for each System's y
    GLFW.PollEvents() #see if we got a shutdown request
end

IODevices.should_close(xp::XPConnect) = GLFW.WindowShouldClose(xp.window)
IODevices.shutdown!(xp::XPConnect) = GLFW.DestroyWindow(xp.window)

end #module