module TestC172Rv1

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight

export test_c172rv1


function test_c172rv1()
    @testset verbose = true "Cessna172Rv1" begin

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

            ac = System(Cessna172Rv1());

            init_kinematics!(ac, kin_init)

            f_ode!(ac, env) #make sure we're on the ground
            @test ac.y.airframe.ldg.left.strut.wow == true

            @test @ballocated(f_ode!($ac, $env)) == 0
            @test @ballocated(f_step!($ac)) == 0
            @test @ballocated(f_disc!($ac, 0.02, $env)) == 0

        end

    return nothing

end

function test_sim(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        h_trn = HOrth(608.55);

        ac = Cessna172Rv1();
        env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
        world = SimpleWorld(ac, env) |> System;

        kin_init = KinematicInit(
            v_eOb_n = [30, 0, 0],
            ω_lb_b = [0, 0, 0],
            q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
            loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
            h = h_trn + 1.9 + 2200.5);

        init_kinematics!(world, kin_init)

        world.u.ac.avionics.eng_start = true #engine start switch on
        world.u.env.atm.wind.v_ew_n .= [0, 0, 0]

        sys_io! = let

            function (u, s, y, t, params)

                u.ac.avionics.throttle = 0.2
                u.ac.avionics.aileron = (t < 5 ? 0.25 : 0.0)
                u.ac.avionics.elevator = 0.0
                u.ac.avionics.rudder = 0.0
                u.ac.avionics.brake_left = 1
                u.ac.avionics.brake_right = 1

            end
        end

        sim = Simulation(world; Δt = 0.02, t_end = 300, sys_io!, adaptive = false)
        Sim.run!(sim, verbose = true)

        # plots = make_plots(sim; Plotting.defaults...)
        plots = make_plots(TimeHistory(sim).ac.kinematics; Plotting.defaults...)
        save && save_plots(plots, save_folder = joinpath("tmp", "sim_test"))

        # return sim
        return world

    end

end

function test_sim_paced(; save::Bool = true)

    h_trn = HOrth(601.55);

    ac = Cessna172Rv1();
    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
    world = SimpleWorld(ac, env) |> System;

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.3),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0);

    init_kinematics!(world, kin_init)

    sim = Simulation(world; Δt = 0.02, t_end = 300)

    interfaces = Vector{IODevices.Interface}()
    for joystick in get_connected_joysticks()
        push!(interfaces, attach_io!(sim, joystick))
    end

    # xp = XPConnect()
    xp = XPConnect(host = IPv4("192.168.1.2"))
    push!(interfaces, attach_io!(sim, xp))

    @sync begin
        for interface in interfaces
            Threads.@spawn IODevices.start!(interface)
        end
        Threads.@spawn Sim.run_paced!(sim; rate = 1, verbose = true)
    end

    plots = make_plots(TimeHistory(sim).ac.kinematics; Plotting.defaults...)
    # plots = make_plots(TimeHistory(sim); Plotting.defaults...)
    save && save_plots(plots, save_folder = joinpath("tmp", "paced_sim_test"))

    return nothing

end


end #module