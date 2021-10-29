using Sockets

Base.@kwdef struct XPInterface
    socket::UDPSocket = UDPSocket()
    host::IPv4 = IPv4("127.0.0.1")
    port::Integer = 49009
end

function display_test()

    #XPlane expects native byte-ordering, so no need to convert anything
    xp = XPInterface()

    dref_id = "sim/operation/override/override_planepath"
    dref_value = 1
    set_dref(xp, dref_id, dref_value)
    set_position(xp)

end

function set_dref(xp::XPInterface, dref_id::AbstractString, dref_value::Real)

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

function set_position(xp::XPInterface)
# function set_position(xp::XPInterface, q_nb::Rotation, p::Abstract3DPosition)

    ac = UInt8(0)
    lat = Float64(45)
    lon = Float64(0)
    alt = Float64(10)
    psi = Float32(1)
    theta = Float32(2)
    phi = Float32(3)
    gear = Float32(-998)

    buffer = IOBuffer()
    write(buffer, b"POSI\0", ac, lat, lon, alt, psi, theta, phi, gear)

    send(xp.socket, xp.host, xp.port, buffer.data)

end


function set_dref(dref_id::AbstractString, dref_value::AbstractVector{<:Real})

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        dref_value |> length |> UInt8,
        Vector{Float32}(dref_value))

    send(xp.socket, xp.host, xp.port, buffer.data)
end
    #we could create a method set_dataref for datarefs whose value is an array
    #rather than a scalar. for these, we would cast the value to a Vector{Float32}
