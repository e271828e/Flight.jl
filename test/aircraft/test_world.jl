module TestWorld

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight

export test_world

function test_world()
    @testset verbose = true "World" begin

        @testset verbose = true "Performance" begin test_system() end
        @testset verbose = true "Simulation" begin test_sim(save = false) end

    end
end

function test_system()

    h_trn = HOrth(608.55);

    ac = Cessna172R();
    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
    world = SimpleWorld(ac, env) |> System

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.8 + 0);

    init_kinematics!(world, kin_init)

    #make sure we're on the ground
    f_ode!(world)
    @test world.ac.y.airframe.ldg.left.strut.wow == true

    @test @ballocated(f_ode!($world)) == 0
    @test @ballocated(f_step!($world)) == 0
    @test @ballocated(f_disc!($world, 0.2)) == 0

    return nothing

end

function test_sim(; save::Bool = true)

    h_trn = HOrth(608.55);

    ac = Cessna172R();
    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
    world = SimpleWorld(ac, env) |> System;

    kin_init = KinematicInit(
        v_eOb_n = [30, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 2200.5);

    init_kinematics!(world, kin_init)

    world.ac.u.avionics.eng_start = true #engine start switch on
    world.env.atm.u.wind.v_ew_n .= [0, 0, 0]

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

    sim = Simulation(world; t_end = 300, sys_io!, adaptive = true)
    Sim.run!(sim, verbose = true)

    # plots = make_plots(sim; Plotting.defaults...)
    plots = make_plots(TimeHistory(sim).ac.kinematics; Plotting.defaults...)
    save && save_plots(plots, save_folder = joinpath("tmp", "sim_test"))

    # return sim

end



function test_sim_paced(; save::Bool = true)

    h_trn = HOrth(601.55);

    ac = Cessna172R();
    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
    world = SimpleWorld(ac, env) |> System;

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.3),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0);

    init_kinematics!(world, kin_init)

    sim = Simulation(world; t_end = 300)

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