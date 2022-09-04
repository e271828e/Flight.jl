module XPlane

using Sockets
using UnPack

using ..IODevices

export XPConnect


################################################################################
############################### XPConnect ####################################

Base.@kwdef struct XPConnect <: IODevice
    socket::UDPSocket = UDPSocket()
    host::IPv4 = IPv4("127.0.0.1")
    port::Integer = 49009
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

function display_text(xp::XPConnect, txt::AbstractString, x::Integer = -1, y::Integer = -1)

    buffer = IOBuffer()
    txt_ascii = ascii(txt)
    write(buffer,
        b"TEXT\0",
        x |> Int32,
        y |> Int32,
        txt_ascii |> length |> UInt8,
        txt_ascii |> codeunits)

    send(xp.socket, xp.host, xp.port, buffer.data)
end

disable_physics!(xp::XPConnect) = set_dref(xp, "sim/operation/override/override_planepath", 1)

init!(xp::XPConnect) = disable_physics!(xp)


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


end #module