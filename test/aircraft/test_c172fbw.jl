module TestC172FBW

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

export test_c172fbw

function test_c172fbw()
    @testset verbose = true "Cessna172FBW" begin

        test_trimming()
        test_linearization()
        test_system_methods()
        test_sim(save = false)

    end
end

function test_trimming()

    @testset verbose = true "Trimming" begin

        vehicle = System(C172FBW.Vehicle())
        trim_params = C172.TrimParameters()
        state = C172.TrimState()

        f_target = C172.get_f_target(vehicle, trim_params)

        @test @ballocated($f_target($state)) === 0

        success, _ = Systems.init!(vehicle, trim_params)

        @test success

    end #testset

end #function

function test_linearization()

    #needs actual tests

    @testset verbose = true "Linearization" begin
        ss = Cessna172FBW(NED()) |> System |> Control.Continuous.LinearizedSS

    end #testset

end #function

function test_system_methods()

        @testset verbose = true "System Methods" begin

            trn = HorizontalTerrain()
            loc = NVector()
            trn_data = TerrainData(trn, loc)
            kin_init = KinInit( h = trn_data.altitude + 1.8);

            ac_WA = Cessna172FBW(WA(), trn) |> System;
            ac_ECEF = Cessna172FBW(ECEF(), trn) |> System;
            ac_NED = Cessna172FBW(NED(), trn) |> System;

            Systems.init!(ac_WA, kin_init)
            Systems.init!(ac_ECEF, kin_init)
            Systems.init!(ac_NED, kin_init)

            f_ode!(ac_WA)
            #make sure we are on the ground to ensure landing gear code coverage
            @test ac_WA.y.vehicle.components.ldg.left.strut.wow == true

            #all three kinematics implementations must be supported, no allocations
            @test @ballocated(f_ode!($ac_WA)) == 0
            @test @ballocated(f_step!($ac_WA)) == 0
            @test @ballocated(f_disc!($ac_WA)) == 0

            @test @ballocated(f_ode!($ac_ECEF)) == 0
            @test @ballocated(f_step!($ac_ECEF)) == 0
            @test @ballocated(f_disc!($ac_ECEF)) == 0

            @test @ballocated(f_ode!($ac_NED)) == 0
            @test @ballocated(f_step!($ac_NED)) == 0
            @test @ballocated(f_disc!($ac_NED)) == 0

        end

    return nothing

end

function test_sim(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        ac = Cessna172FBW() |> System;

        mid_cg_pld = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

        ac.vehicle.atmosphere.u.v_ew_n .= [0, 0, 0]

        trim_params = C172.TrimParameters(
        Ob = Geographic(LatLon(), HOrth(1000)),
        EAS = 25.0,
        γ_wb_n = 0.0,
        x_fuel = 0.5,
        flaps = 1.0,
        payload = mid_cg_pld)

        exit_flag, trim_state = Systems.init!(ac, trim_params)
        @test exit_flag === true

        user_callback! = let

            function (ac)

                u_act = ac.vehicle.components.act.u
                t = ac.t[]

                # u_act.throttle_cmd = 0.2
                # u_act.aileron_cmd = (t < 5 ? 0.25 : 0.0)
                # u_act.elevator_cmd = 0.0
                # u_act.rudder_cmd = 0.1
                # u_act.aileron_cmd = 0.1
                # u_act.brake_left = 1
                # u_act.brake_right = 1

            end
        end

        sim = Simulation(ac; t_end = 30, user_callback!, adaptive = true)
        Sim.run!(sim)

        # plots = make_plots(sim; Plotting.defaults...)
        kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
        air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
        dyn_plots = make_plots(TimeSeries(sim).vehicle.dynamics; Plotting.defaults...)
        save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172fbw", "sim", "kin"))
        save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172fbw", "sim", "air"))
        save && save_plots(dyn_plots, save_folder = joinpath("tmp", "test_c172fbw", "sim", "dyn"))

    end

end


function test_sim_interactive(; save::Bool = true)

    h_trn = HOrth(427.2);

    # # on ground
    # initializer = KinInit(
    #     loc = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)),
    #     q_nb = REuler(deg2rad(157), 0, 0),
    #     h = h_trn + 1.81);

    # on air, automatically trimmed
    initializer = C172.TrimParameters(
        Ob = Geographic(LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)), HEllip(650)))

    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172FBW(WA(), trn) |> System;

    sim = Simulation(ac; dt = 1/60, Δt = 1/60, t_end = 1000)

    reinit!(sim, initializer)

    for joystick in get_connected_joysticks()
        Sim.attach!(sim, joystick)
    end

    xpc = XP12Client()
    # xpc = XP12Client(address = IPv4("192.168.1.2"))
    Sim.attach!(sim, xpc)

    Sim.run_interactive!(sim; pace = 1)

    kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172rpa_v1", "sim_interactive", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172rpa_v1", "sim_interactive", "air"))

    return nothing

end




end #module