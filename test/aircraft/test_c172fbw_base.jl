module TestC172FBWBase

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight.FlightCore.Systems
using Flight.FlightCore.Sim
using Flight.FlightCore.Plotting
using Flight.FlightCore.IODevices
using Flight.FlightCore.Joysticks
using Flight.FlightCore.XPC

using Flight.FlightPhysics.Attitude
using Flight.FlightPhysics.Geodesy
using Flight.FlightPhysics.Kinematics
using Flight.FlightPhysics.Environment

using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.World

using Flight.FlightAircraft.C172FBWBase

export test_c172fbw_base


function test_c172fbw_base()
    @testset verbose = true "Cessna172FBWBase" begin

        test_system_methods()
        test_sim(save = false)

    end
end

function test_system_methods()

        @testset verbose = true "System Methods" begin

            env = SimpleEnvironment() |> System

            loc = NVector()
            trn_data = TerrainData(env.trn, loc)
            kin_init = KinematicInit( h = trn_data.altitude + 1.8);

            ac_LTF = System(Cessna172FBWBase(LTF()));
            ac_ECEF = System(Cessna172FBWBase(ECEF()));
            ac_NED = System(Cessna172FBWBase(NED()));

            init_kinematics!(ac_LTF, kin_init)
            init_kinematics!(ac_ECEF, kin_init)
            init_kinematics!(ac_NED, kin_init)

            f_ode!(ac_LTF, env)
            #make sure we are on the ground to ensure landing gear code coverage
            @test ac_LTF.y.physics.airframe.ldg.left.strut.wow == true

            #all three kinematics implementations must be supported, no allocations
            @test @ballocated(f_ode!($ac_LTF, $env)) == 0
            @test @ballocated(f_step!($ac_LTF)) == 0
            @test @ballocated(f_disc!($ac_LTF, 1, $env)) == 0

            @test @ballocated(f_ode!($ac_ECEF, $env)) == 0
            @test @ballocated(f_step!($ac_ECEF)) == 0
            @test @ballocated(f_disc!($ac_ECEF, 1, $env)) == 0

            @test @ballocated(f_ode!($ac_NED, $env)) == 0
            @test @ballocated(f_step!($ac_NED)) == 0
            @test @ballocated(f_disc!($ac_NED, 1, $env)) == 0

        end

    return nothing

end

function test_sim(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        h_trn = HOrth(608.55);

        ac = Cessna172FBWBase();
        env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
        world = SimpleWorld(ac, env) |> System;

        kin_init = KinematicInit(
            v_eOb_n = [30, 0, 0],
            ω_lb_b = [0, 0, 0],
            q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
            loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
            h = h_trn + 1.9 + 2200.5);

        init_kinematics!(world, kin_init)

        world.ac.physics.airframe.act.u.eng_start = true #engine start switch on
        world.env.atm.wind.u.v_ew_n .= [0, 0, 0]

        sys_io! = let

            function (world)

                u_act = world.ac.physics.airframe.act.u
                t = world.t[]

                u_act.throttle_cmd = 0.2
                u_act.aileron_cmd = (t < 5 ? 0.25 : 0.0)
                u_act.elevator_cmd = 0.0
                u_act.rudder_cmd = 0.0
                u_act.brake_left = 1
                u_act.brake_right = 1

            end
        end

        sim = Simulation(world; t_end = 300, sys_io!, adaptive = true)
        Sim.run!(sim, verbose = true)

        # plots = make_plots(sim; Plotting.defaults...)
        plots = make_plots(TimeHistory(sim).ac.physics.kinematics; Plotting.defaults...)
        save && save_plots(plots, save_folder = joinpath("tmp", "test_c172r_base", "sim"))

    end

end


function test_sim_paced(; save::Bool = true)

    h_trn = HOrth(601.55);

    ac = Cessna172FBWBase();
    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
    world = SimpleWorld(ac, env) |> System;

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.3),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0);

    init_kinematics!(world, kin_init)

    sim = Simulation(world; dt = 0.02, Δt = 0.02, t_end = 300)

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
        Threads.@spawn Sim.run_paced!(sim; pace = 1, verbose = true)
    end

    plots = make_plots(TimeHistory(sim).ac.physics.kinematics; Plotting.defaults...)
    # plots = make_plots(TimeHistory(sim); Plotting.defaults...)
    save && save_plots(plots, save_folder = joinpath("tmp", "test_c172r_base", "sim_paced"))

    return nothing

end




end #module