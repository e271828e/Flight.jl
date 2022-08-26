module Output

using Sockets
using UnPack

using Flight.GUI

export XPConnect

export Dashboard



################################################################################
############################### AbstractDevice #################################

abstract type AbstractDevice end

#each device should control its own update rate from within the update! method,
#preferably via calls to blocking functions (such as GLFW.SwapBuffers with
#GLFW.SwapInterval > 0)

init!(device::D) where {D<:AbstractDevice} = MethodError(init!, (device, ))
shutdown!(device::D) where {D<:AbstractDevice} = MethodError(shutdown!, (device, ))
should_close(device::D) where {D<:AbstractDevice} = MethodError(should_close, (device, ))
update!(device::D, data) where {D<:AbstractDevice} = MethodError(update!, (device, data))


################################################################################
############################## Interface #################################

Base.@kwdef struct Interface{D <: AbstractDevice, C <: Channel}
    device::D #underlying IO device
    channel::C #unbuffered channel where the Simulation will put! its output
    ext_shutdown::Bool #whether to observe shutdown requests received by the IO device
end

init!(output::Interface) = init!(output.device)
shutdown!(output::Interface) = shutdown!(output.device)
should_close(output::Interface) = should_close(output.device)

run!(output::Interface; verbose = true) = (Threads.@spawn _run!(output; verbose))

function _run!(io::Interface{D}; verbose = true) where {D}

    @unpack device, channel, ext_shutdown = io

    verbose && println("Output Interface: Starting at thread $(Threads.threadid())...")

    init!(device)

    try

        while true

            data = take!(channel)
            update!(device, data)

            if ext_shutdown && should_close(device)
                println("Output Interface: Shutdown requested")
                break
            end

        end

    catch ex
        if ex isa InvalidStateException
            println("Output Interface: Channel closed")
        else
            println("Output Interface: Error during execution: $ex")
        end

    finally
        println("Output Interface: Exiting...")
        shutdown!(device)
    end

end

################################################################################
############################# Dashboard ##################################

Base.@kwdef struct Dashboard{F <: Function} <: AbstractDevice
    renderer::CImGuiRenderer
    draw::F
end

Dashboard(draw::Function; kwargs...) = Dashboard(CImGuiRenderer(; kwargs...), draw)

function init!(db::Dashboard)
    GUI.init!(db.renderer)
end
shutdown!(db::Dashboard) = GUI.shutdown!(db.renderer)
should_close(db::Dashboard) = GUI.should_close(db.renderer)
update!(db::Dashboard, data) = GUI.render!(db.renderer, db.draw, data)


"""
to be updated as an AbstractDevice for the new Interface
"""


################################################################################
############################### XPConnect ####################################

# Base.@kwdef struct XPConnect <: AbstractDevice
#     socket::UDPSocket = UDPSocket()
#     host::IPv4 = IPv4("127.0.0.1")
#     port::Integer = 49009
# end

# function set_dref(xp::XPConnect, dref_id::AbstractString, dref_value::Real)

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

#     send(xp.socket, xp.host, xp.port, buffer.data)
# end

# function set_dref(xp::XPConnect, dref_id::AbstractString, dref_value::AbstractVector{<:Real})

#     buffer = IOBuffer()
#     write(buffer,
#         b"DREF\0",
#         dref_id |> length |> UInt8,
#         dref_id |> ascii |> codeunits,
#         dref_value |> length |> UInt8,
#         Vector{Float32}(dref_value))

#     send(xp.socket, xp.host, xp.port, buffer.data)
# end

# function display_text(xp::XPConnect, txt::AbstractString, x::Integer = -1, y::Integer = -1)

#     buffer = IOBuffer()
#     txt_ascii = ascii(txt)
#     write(buffer,
#         b"TEXT\0",
#         x |> Int32,
#         y |> Int32,
#         txt_ascii |> length |> UInt8,
#         txt_ascii |> codeunits)

#     send(xp.socket, xp.host, xp.port, buffer.data)
# end

# disable_physics!(xp::XPConnect) = set_dref(xp, "sim/operation/override/override_planepath", 1)

# init!(xp::XPConnect) = disable_physics!(xp)


# function set_position!(xp::XPConnect; lat = -998, lon = -998, h = -998,
#                         psi = -998, theta = -998, phi = -998,
#                         aircraft::Integer = 0)

#     #all angles must be in degrees
#     buffer = IOBuffer()
#     write(buffer,
#         b"POSI\0", UInt8(aircraft),
#         Float64(lat), Float64(lon), Float64(h),
#         Float32(theta), Float32(phi), Float32(psi),
#         Float32(-998)) #last one is landing gear (?!)

#     send(xp.socket, xp.host, xp.port, buffer.data)

# end



end #module