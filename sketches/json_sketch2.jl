module JSONSketch2

using Sockets, UnPack
using StructTypes, JSON3

using Flight
using Flight.FlightCore.Sim
using Flight.FlightCore.IODevices
using Flight.FlightCore.Networking
using Flight.FlightAircraft.C172RPAv1

export json_sketch2

function output_callback(sim_out::Sim.Output)::Vector{UInt8}
    freq = 0.1
    φ_sp_max = π/6
    φ_sp = φ_sp_max * sin(2π*freq*sim_out.t)

    #these are all valid empty JSON entities. when passed to JSON3.write, they
    #yield respectively "\"\"", "[]" and "{}", all of length 2
    cmd = ""
    # cmd = []
    # cmd = Dict()

    if sim_out.t > 10
        #these enums will be automatically cast to Ints per the StructTypes
        #methods defined in C172RPAv1
        cmd = (
            vrt_gdc_mode_req = C172RPAv1.vrt_gdc_alt,
            lat_ctl_mode_req = C172RPAv1.lat_φ_β,
            φ_sp = φ_sp,
        )
        # cmd = (
        #     vrt_gdc_mode_req = 1,
        #     lat_ctl_mode_req = 2,
        #     φ_sp = φ_sp,
        # )
    end

    json_cmd = JSON3.write(cmd)
    return Vector{UInt8}(json_cmd)
end

#AvionicsU is declared as StructTypes.Mutable() in C172RPAv1, so JSON3 can
#automatically read a JSON string into one or more of its fields
function assign_callback!(sys::System{<:Cessna172RPAv1}, data::Vector{UInt8}, ::InputMapping)
    #warning! String(data) empties the original data::Vector{UInt8}, so further
    #calls would return an empty string
    str = String(data)

    #if length(str) < 2, it cannot be a valid JSON string. if length(str) == 2
    #it is an empty JSON entity (either string, object or array). instead of
    #this check we could simply call do isempty(JSON3.read(str)) but that would
    #mean parsing the string twice
    length(str) > 2 && JSON3.read!(str, sys.avionics.u)

    # isempty(str) |> println
    # JSON3.read(str) |> isempty |> println
end

function json_sketch2(; save::Bool = true)

    h_trn = HOrth(601.55);

    trn = HorizontalTerrain(altitude = h_trn)
    sys = Cessna172RPAv1(LTF(), trn) |> System;
    sim = Simulation(sys; dt = 1/60, Δt = 1/60, t_end = 20)

    # #on ground
    # kin_init = KinematicInit(
    #     loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
    #     h = h_trn + 1.81);

    #on air, automatically trimmed by reinit!
    kin_init = C172.TrimParameters(
        Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)), HEllip(1050)))

    #initialize simulated system
    reinit!(sim, kin_init)

    #setup IO devices
    for joystick in get_connected_joysticks()
        Sim.attach!(sim, joystick)
    end
    xpc = XPCClient()
    # xpc = XPCClient(address = IPv4("192.168.1.2"))
    Sim.attach!(sim, xpc)
    Sim.attach!(sim, UDPClient(; port = 49017, output_callback))
    Sim.attach!(sim, UDPServer(; port = 49017, assign_callback!))

    #trigger compilation of parsing methods for AvionicsU before launching the
    #simulation
    JSON3.read!(JSON3.write(sys.avionics.u, allow_inf=true), sys.avionics.u; allow_inf=true)

    # return
    Sim.run_paced!(sim)

    kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172rpa_v1", "sim_paced", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172rpa_v1", "sim_paced", "air"))

    return nothing

end

end #module