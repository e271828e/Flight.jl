module TestC172RBase

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight.FlightCore
using Flight.FlightCore.Sim
using Flight.FlightCore.Visualization

using Flight.FlightPhysics
using Flight.FlightComponents

using Flight.FlightAircraft.AircraftBase
using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172R
using Flight.FlightAircraft.C172RBase

export test_c172r_base


function test_c172r_base()
    @testset verbose = true "Cessna172RBase" begin

        test_system_methods()
        test_sim(save = false)

    end
end

function test_system_methods()

        @testset verbose = true "System Methods" begin

            trn = HorizontalTerrain()
            loc = NVector()
            trn_data = TerrainData(trn, loc)
            kin_init = KinematicInit( h = trn_data.altitude + 1.8);

            ac_LTF = System(Cessna172RBase(LTF(), trn));
            ac_ECEF = System(Cessna172RBase(ECEF(), trn));
            ac_NED = System(Cessna172RBase(NED(), trn));

            Systems.init!(ac_LTF, kin_init)
            Systems.init!(ac_ECEF, kin_init)
            Systems.init!(ac_NED, kin_init)

            f_ode!(ac_LTF)
            #make sure we are on the ground to ensure landing gear code coverage
            @test ac_LTF.y.physics.platform.ldg.left.strut.wow == true

            #all three kinematics implementations must be supported, no allocations
            @test @ballocated(f_ode!($ac_LTF)) == 0
            @test @ballocated(f_step!($ac_LTF)) == 0
            @test @ballocated(f_disc!($ac_LTF, 1)) == 0

            @test @ballocated(f_ode!($ac_ECEF)) == 0
            @test @ballocated(f_step!($ac_ECEF)) == 0
            @test @ballocated(f_disc!($ac_ECEF, 1)) == 0

            @test @ballocated(f_ode!($ac_NED)) == 0
            @test @ballocated(f_step!($ac_NED)) == 0
            @test @ballocated(f_disc!($ac_NED, 1)) == 0

        end

    return nothing

end

function test_sim(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        ac = Cessna172RBase() |> System;

        mid_cg_pld = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

        ac.physics.atmosphere.u.v_ew_n .= [0, 0, 0]

        trim_params = C172.TrimParameters(
        Ob = Geographic(LatLon(), HOrth(1000)),
        EAS = 25.0,
        γ_wOb_n = 0.0,
        x_fuel = 0.5,
        flaps = 1.0,
        payload = mid_cg_pld)

        exit_flag, trim_state = trim!(ac, trim_params)
        @test exit_flag === true

        sys_io! = let

            function (ac)

                u_act = ac.physics.platform.act.u
                t = ac.t[]

            end
        end

        sim = Simulation(ac; t_end = 30, sys_io!, adaptive = true)
        Sim.run!(sim)

        # plots = make_plots(sim; Plotting.defaults...)
        kin_plots = make_plots(TimeSeries(sim).physics.kinematics; Plotting.defaults...)
        air_plots = make_plots(TimeSeries(sim).physics.air; Plotting.defaults...)
        rb_plots = make_plots(TimeSeries(sim).physics.dynamics; Plotting.defaults...)
        save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172r_base", "sim", "kin"))
        save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172r_base", "sim", "air"))
        save && save_plots(rb_plots, save_folder = joinpath("tmp", "test_c172r_base", "sim", "dynamics"))

    end

end


function test_sim_paced(; save::Bool = true)

    h_trn = HOrth(601.55);

    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172RBase(LTF(), trn) |> System

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.0),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0);

    Systems.init!(ac, kin_init)

    sim = Simulation(ac; dt = 0.02, Δt = 0.02, t_end = 300)

    interfaces = Vector{IODevices.Interface}()
    for joystick in get_connected_joysticks()
        push!(interfaces, attach_io!(sim, joystick))
    end

    # xp = XPCDevice()
    xp = XPCDevice(host = IPv4("192.168.1.2"))
    push!(interfaces, attach_io!(sim, xp))

    @sync begin
        for interface in interfaces
            Threads.@spawn IODevices.start!(interface)
        end
        Threads.@spawn Sim.run_paced!(sim; pace = 1)
    end

    plots = make_plots(TimeSeries(sim).physics.kinematics; Plotting.defaults...)
    # plots = make_plots(TimeSeries(sim); Plotting.defaults...)
    save && save_plots(plots, save_folder = joinpath("tmp", "test_c172r_base", "sim_paced"))

    return nothing

end




end #module