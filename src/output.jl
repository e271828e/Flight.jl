module Output

using Sockets

using Flight.Attitude
using Flight.Geodesy
using Flight.Kinematics

export XPInterface, disable_physics, set_position

Base.@kwdef struct XPInterface
    socket::UDPSocket = UDPSocket()
    host::IPv4 = IPv4("127.0.0.1")
    port::Integer = 49009
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

function set_dref(xp::XPInterface, dref_id::AbstractString, dref_value::AbstractVector{<:Real})

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        dref_value |> length |> UInt8,
        Vector{Float32}(dref_value))

    send(xp.socket, xp.host, xp.port, buffer.data)
end

disable_physics(xp::XPInterface) = set_dref(xp, "sim/operation/override/override_planepath", 1)

function set_position(xp::XPInterface, init::KinInit, aircraft::Integer = 0)
    set_position(xp, KinData(init).pos, aircraft)
end

function set_position(xp::XPInterface, pos::PosData, aircraft::Integer = 0)

    llh = Geographic(pos.ϕ_λ, pos.h_o)
    euler = REuler(pos.q_nb)

    lat = rad2deg(llh.l2d.ϕ)
    lon = rad2deg(llh.l2d.λ)
    alt = llh.alt
    psi = rad2deg(euler.ψ)
    theta = rad2deg(euler.θ)
    phi = rad2deg(euler.φ)
    set_position(xp; lat, lon, alt, psi, theta, phi, aircraft)

end

function set_position(xp::XPInterface; lat = -998, lon = -998, alt = -998,
                        psi = -998, theta = -998, phi = -998,
                        aircraft::Integer = 0)

    #all angles must be in degrees
    buffer = IOBuffer()
    write(buffer,
        b"POSI\0", UInt8(aircraft),
        Float64(lat), Float64(lon), Float64(alt),
        Float32(theta), Float32(phi), Float32(psi),
        Float32(-998)) #last one is landing gear (?!)

    send(xp.socket, xp.host, xp.port, buffer.data)

end

end #module