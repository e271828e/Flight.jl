module TestC172RPA

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

export test_c172rpa

function test_c172rpa()
    @testset verbose = true "Cessna 172RPA" begin
        test_trimming()
        test_linearization()
        test_system_methods()
        test_sim(save = false)
    end
end

function test_trimming()

    @testset verbose = true "Trimming" begin

        trim_params = C172.TrimParameters()
        state = C172.TrimState()

        #test on vehicle
        vehicle = System(C172RPA.Vehicle())
        f_target = C172.get_f_target(vehicle, trim_params)
        @test @ballocated($f_target($state)) === 0
        success, _ = Systems.init!(vehicle, trim_params)
        @test success

        #test on direct variant
        ac = System(Cessna172RPA())
        success, _ = Systems.init!(ac, trim_params)
        @test success

    end #testset

end #function

function test_linearization()

    @testset verbose = true "Linearization" begin

        ss = Cessna172RPA(NED()) |> System |> Control.Continuous.LinearizedSS

    end #testset

end #function


function test_system_methods()

        @testset verbose = true "System Methods" begin

            trn = HorizontalTerrain()
            loc = NVector()
            trn_data = TerrainData(trn, loc)
            kin_init = KinInit( h = trn_data.altitude + 1.8);

            ac_WA = System(Cessna172RPA(WA(), trn));
            ac_ECEF = System(Cessna172RPA(ECEF(), trn));
            ac_NED = System(Cessna172RPA(NED(), trn));

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

        ac = Cessna172RPA() |> System;

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

            end
        end

        sim = Simulation(ac; t_end = 30, user_callback!, adaptive = true)
        Sim.run!(sim)

        # plots = make_plots(sim; Plotting.defaults...)
        kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
        air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
        rb_plots = make_plots(TimeSeries(sim).vehicle.accelerations; Plotting.defaults...)
        save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172rpa", "sim", "kin"))
        save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172rpa", "sim", "air"))
        save && save_plots(rb_plots, save_folder = joinpath("tmp", "test_c172rpa", "sim", "accelerations"))

    end

end


function test_sim_interactive(; save::Bool = true)

    h_trn = HOrth(601.55);

    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172RPA(WA(), trn) |> System

    kin_init = KinInit(
        v_eb_n = [0, 0, 0],
        ω_wb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.0),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0);

    Systems.init!(ac, kin_init)

    sim = Simulation(ac; dt = 0.02, Δt = 0.02, t_end = 300)

    for joystick in get_connected_joysticks()
        Sim.attach!(sim, joystick)
    end

    xpc = XPCClient()
    # xpc = XPCClient(address = IPv4("192.168.1.2"))
    Sim.attach!(sim, xpc)

    Sim.run_interactive!(sim)

    plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
    # plots = make_plots(TimeSeries(sim); Plotting.defaults...)
    save && save_plots(plots, save_folder = joinpath("tmp", "test_c172rpa", "sim_interactive"))

    return nothing

end

end #module