module TestC172

using Test
using UnPack
using BenchmarkTools
using Sockets
using OrdinaryDiffEq: Heun, RK4, Tsit5

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

export test_c172

function test_c172()
    @testset verbose = true "Cessna 172" begin
    end
end

function test_sim(; ac::Cessna172 = Cessna172Sv0(),
                loc::Abstract2DLocation = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)),
                h_trn::HOrth = HOrth(427.2),
                ψ::Real = deg2rad(157),
                situation::Symbol = :ground,
                interactive::Bool = false,
                t_end::Real = 1000,
                save::Bool = true)

    world = SimpleWorld(ac, SimpleAtmosphere(), HorizontalTerrain(altitude = h_trn)) |> System

    if situation === :ground
        initializer = KinInit(; loc, q_nb = REuler(ψ, 0, 0), h = h_trn + 1.81);
    elseif situation === :air
        EAS = 50.0
        flaps = 0.0
        γ_wb_n = 0.0
        x_fuel = 0.5
        payload = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

        initializer = C172.TrimParameters(; Ob = Geographic(loc, HEllip(650)), EAS, γ_wb_n, x_fuel, flaps, payload)
    else
        error("Unknown situation: $situation")
    end

    user_callback! = let
        function (world)
            t = world.t[]
        end
    end

    sim = Simulation(world; algorithm = RK4(), dt = 1/60, t_end, user_callback!)
    # sim = Simulation(world; algorithm = Tsit5(), dt = 1/60, t_end, user_callback!)
    # sim = Simulation(world; algorithm = Heun(), dt = 0.01, t_end, user_callback!)
    reinit!(sim, initializer)

    if interactive
        xp = XPlane12Output(address = IPv4("127.0.0.1"), port = 49000)
        for joystick in get_connected_joysticks()
            Sim.attach!(sim, joystick)
        end
        Sim.attach!(sim, xp)
        Sim.run_interactive!(sim)
    else
        Sim.run!(sim)
    end

    kin_plots = make_plots(TimeSeries(sim).ac.vehicle.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).ac.vehicle.airflow; Plotting.defaults...)
    dyn_plots = make_plots(TimeSeries(sim).ac.vehicle.dynamics; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172", "sim", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172", "sim", "air"))
    save && save_plots(dyn_plots, save_folder = joinpath("tmp", "test_c172", "sim", "dyn"))

end



end