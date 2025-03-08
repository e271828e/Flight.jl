module TestC172

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

export test_c172

function test_c172()
    @testset verbose = true "Cessna 172" begin
    end
end

# function test_sim(ac::System{<:Cessna172} = Cess; save::Bool = true)
function test_sim(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        h_trn = HOrth(427.2);

        # on ground
        # initializer = KinInit(
        #     loc = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)),
        #     q_nb = REuler(deg2rad(157), 0, 0),
        #     h = h_trn + 1.81);

        # on air, automatically trimmed
        mid_cg_pld = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)
        initializer = C172.TrimParameters(
            Ob = Geographic(LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)), HEllip(650)),
            EAS = 25.0,
            γ_wb_n = 0.0,
            x_fuel = 0.5,
            flaps = 1.0,
            payload = mid_cg_pld)

        trn = HorizontalTerrain(altitude = h_trn)

        ac = Cessna172Xv1(WA(), trn) |> System;

        ac.vehicle.atmosphere.u.v_ew_n .= [0, 0, 0]
        user_callback! = let
            function (ac)
                # u_act = ac.vehicle.components.act.u
                t = ac.t[]
            end
        end

        sim = Simulation(ac; t_end = 30, user_callback!)
        reinit!(sim, initializer)
        Sim.run!(sim)

        # plots = make_plots(sim; Plotting.defaults...)
        kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
        air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
        dyn_plots = make_plots(TimeSeries(sim).vehicle.dynamics; Plotting.defaults...)
        save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172x", "sim", "kin"))
        save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172x", "sim", "air"))
        save && save_plots(dyn_plots, save_folder = joinpath("tmp", "test_c172x", "sim", "dyn"))

    end

end

function test_sim_interactive(; save::Bool = true)

    h_trn = HOrth(427.2);

    # on ground
    initializer = KinInit(
        loc = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)),
        q_nb = REuler(deg2rad(157), 0, 0),
        h = h_trn + 1.81);

    # # on air, automatically trimmed
    # initializer = C172.TrimParameters(
    #     Ob = Geographic(LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)), HEllip(650)))

    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172Xv1(WA(), trn) |> System;

    sim = Simulation(ac; dt = 1/60, Δt = 1/60, t_end = 1000)

    reinit!(sim, initializer)

    for joystick in get_connected_joysticks()
        Sim.attach!(sim, joystick)
    end

    xpc = XPlane12Output()
    # xpc = XPlane12Output(address = IPv4("192.168.1.2"))
    Sim.attach!(sim, xpc)

    Sim.run_interactive!(sim; pace = 1)

    kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172x1", "sim_interactive", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172x1", "sim_interactive", "air"))

    return nothing

end


end