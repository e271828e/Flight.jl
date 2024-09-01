module JSONFMSSketch

using Sockets, UnPack
using StructTypes, JSON3

using Flight
using Flight.FlightCore.Sim
using Flight.FlightCore.IODevices
using Flight.FlightCore.Network
using Flight.FlightAircraft.C172RPAv1

export json_fms_sketch

function output_callback(sim_data::Sim.SimData)::Vector{UInt8}

    @unpack vehicle, avionics = sim_data.y
    @unpack kinematics, air, components = vehicle

    output = (
        lat = rad2deg(kinematics.ϕ_λ.ϕ),
        lon = rad2deg(kinematics.ϕ_λ.λ),
        h_ellip = Float64(kinematics.h_e),
        h_orth = Float64(kinematics.h_o),
        v_gnd_NED = kinematics.v_eOb_n,
        v_wind_NED = air.v_ew_n,
        h_trn = components.ldg.left.strut.Δh,
        wow = components.ldg.left.strut.wow,
        m_fuel = components.fuel.m_avail,
    )

    json_out = JSON3.write(output)
    return Vector{UInt8}(json_out)
end

#ControllerU is declared as StructTypes.Mutable() in C172RPA.FlightControl, so
#JSON3 can automatically read a JSON string into one or more of its fields
function assign_callback!(sys::System{<:Cessna172RPAv1}, data::Vector{UInt8}, ::IOMapping)

    raw_str = String(data)
    json_dict = raw_str |> JSON3.read |> Dict
    if (:chi_sp in keys(json_dict))
        json_dict[:χ_sp] = json_dict[:chi_sp]
        delete!(json_dict, :chi_sp)
    end
    if (:theta_sp in keys(json_dict))
        json_dict[:θ_sp] = json_dict[:theta_sp]
        delete!(json_dict, :theta_sp)
    end
    clean_str = JSON3.write(json_dict)

    length(clean_str) > 2 && JSON3.read!(clean_str, sys.avionics.fcl.u)
end

function json_fms_sketch(; save::Bool = true)

    h_trn = HOrth(601.55);

    trn = HorizontalTerrain(altitude = h_trn)
    sys = Cessna172RPAv1(LTF(), trn) |> System;
    sim = Simulation(sys; dt = 1/60, Δt = 1/60, t_end = 3600)

    #on ground
    kin_init = KinInit(
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.81);

    # #on air, automatically trimmed by reinit!
    # kin_init = C172.TrimParameters(
    #     Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)), HEllip(1050)))

    #initialize simulated system
    reinit!(sim, kin_init)

    #setup IO devices
    for joystick in get_connected_joysticks()
        Sim.attach!(sim, joystick)
    end
    xpc = XPCClient()
    # xpc = XPCClient(address = IPv4("192.168.1.2"))
    Sim.attach!(sim, xpc)
    Sim.attach!(sim, UDPClient(; address = IPv4("192.168.1.2"), port = 5004, output_callback))
    Sim.attach!(sim, UDPServer(; address = IPv4("192.168.1.3"), port = 2004, assign_callback!))

    #trigger compilation of parsing methods for AvionicsU before launching the
    #simulation
    JSON3.read!(JSON3.write(sys.avionics.fcl.u, allow_inf=true), sys.avionics.fcl.u; allow_inf=true)

    Sim.run_interactive!(sim; pace = 5)

    kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172rpa_v1", "sim_interactive", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172rpa_v1", "sim_interactive", "air"))

    return nothing

end

end #module