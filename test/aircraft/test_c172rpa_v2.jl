module TestC172RPAv2

using Test, UnPack, BenchmarkTools, Sockets, JSON3

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

#non-exported stuff
using Flight.FlightLib.Control.Discrete: load_pid_lookup, load_lqr_tracker_lookup
using Flight.FlightAircraft.C172RPA.C172RPAControl: lon_direct, lon_thr_ele, lon_thr_q, lon_thr_θ, lon_thr_EAS, lon_EAS_q, lon_EAS_θ, lon_EAS_clm
using Flight.FlightAircraft.C172RPA.C172RPAControl: lat_direct, lat_p_β, lat_φ_β, lat_χ_β
using Flight.FlightAircraft.C172RPA.C172RPAControl: vrt_gdc_off, vrt_gdc_alt
using Flight.FlightAircraft.C172RPA.C172RPAControl: hor_gdc_off, hor_gdc_line
using Flight.FlightAircraft.C172RPA.C172RPAControl: phase_gnd, phase_air

export test_c172rpa_v2


function test_c172rpa_v2()
    @testset verbose = true "Cessna172 RPAv2" begin

        # test_control_modes()
        # test_guidance_modes()

    end
end

function test_sim(; save::Bool = true)

    h_trn = HOrth(601.55);

    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172RPAv2(WA(), trn) |> System;
    sim = Simulation(ac; t_end = 30)

    #on ground
    initializer = KinInit(
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.81);

    #on air, automatically trimmed
    # initializer = C172.TrimParameters(
    #     Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)), HEllip(1050)))

    reinit!(sim, initializer)

    Sim.run!(sim)

    kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172rpa_v2", "sim", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172rpa_v2", "sim", "air"))

    return nothing

end

function test_sim_interactive(; save::Bool = true)

    h_trn = HOrth(601.55);

    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172RPAv2(WA(), trn) |> System;
    sim = Simulation(ac; dt = 1/60, Δt = 1/60, t_end = 1000)

    #on ground
    initializer = KinInit(
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.81);

    # #on air, automatically trimmed
    # initializer = C172.TrimParameters(
    #     Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)), HEllip(1050)))

    reinit!(sim, initializer)

    for joystick in get_connected_joysticks()
        Sim.attach!(sim, joystick)
    end

    xpc = XPCClient()
    # xpc = XPCClient(address = IPv4("192.168.1.2"))
    Sim.attach!(sim, xpc)
    # Sim.attach!(sim, IODevices.DummyInputDevice())

    Sim.run_interactive!(sim; pace = 1)

    kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172rpa_v2", "sim_interactive", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172rpa_v2", "sim_interactive", "air"))

    return nothing

end


end #module