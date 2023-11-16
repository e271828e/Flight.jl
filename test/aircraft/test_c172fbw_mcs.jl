module TestC172FBWMCS

using Test, UnPack, BenchmarkTools, Sockets

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
using Flight.FlightComponents.Piston
using Flight.FlightComponents.World

using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172FBWMCS

export test_c172fbw_mcs


function test_c172fbw_mcs()
    @testset verbose = true "Cessna172FBW MCS" begin

        test_system_methods()
        test_sim(save = false)

    end
end

function test_system_methods()

        @testset verbose = true "System Methods" begin

            env = SimpleEnvironment() |> System

            loc = NVector()
            trn_data = TerrainData(env.trn, loc)
            kin_init_gnd = KinematicInit( h = trn_data.altitude + 1.8);
            kin_init_air = KinematicInit( h = trn_data.altitude + 1000);

            ac = System(Cessna172FBWMCS());

            #to exercise all airframe functionality, including landing gear, we
            #need to be on the ground with the engine running
            init_kinematics!(ac, kin_init_gnd)
            ac.avionics.u.inceptors.eng_start = true #engine start switch on
            f_disc!(ac, 0.02, env) #run avionics for the engine start signal to propagate
            f_ode!(ac, env)
            f_step!(ac)
            f_ode!(ac, env)
            f_step!(ac)
            @test ac.y.physics.airframe.ldg.left.strut.wow == true
            @test ac.y.physics.airframe.pwp.engine.state === Piston.eng_starting

            @test @ballocated(f_ode!($ac, $env)) == 0
            @test @ballocated(f_step!($ac)) == 0
            @test @ballocated(f_disc!($ac, 0.02, $env)) == 0

            #now we put the aircraft in flight
            init_kinematics!(ac, kin_init_air)
            f_ode!(ac, env)
            @test ac.y.physics.airframe.ldg.left.strut.wow == false
            @test @ballocated(f_ode!($ac, $env)) == 0
            @test @ballocated(f_step!($ac)) == 0

            #testing the different avionics modes for allocations is a bit more
            #involved
        end

    return nothing

end

function test_cas(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        world = SimpleWorld(Cessna172FBWMCS()) |> System;

        design_condition = C172.TrimParameters(
            Ob = Geographic(LatLon(), HOrth(1000)),
            EAS = 40.0,
            γ_wOb_n = 0.0,
            x_fuel = 0.5,
            flaps = 0.0,
            payload = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50))

        exit_flag, trim_state = trim!(design_condition, trim_params)
        @test exit_flag === true

        sys_io! = let

            function (world)

                t = world.t[]

            end
        end

        sim = Simulation(world; dt = 0.01, Δt = 0.01, t_end = 60, sys_io!, adaptive = false)
        # sim = Simulation(world; dt = 0.01, Δt = 0.01, t_end = 60, adaptive = false)
        Sim.run!(sim, verbose = true)

        # plots = make_plots(sim; Plotting.defaults...)
        kin_plots = make_plots(TimeHistory(sim).ac.physics.kinematics; Plotting.defaults...)
        air_plots = make_plots(TimeHistory(sim).ac.physics.air; Plotting.defaults...)
        rb_plots = make_plots(TimeHistory(sim).ac.physics.rigidbody; Plotting.defaults...)
        save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172fbw_mcs", "cas", "kin"))
        save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172fbw_mcs", "cas", "air"))
        save && save_plots(rb_plots, save_folder = joinpath("tmp", "test_c172fbw_mcs", "cas", "rigidbody"))

        return nothing

    end
end

function test_sim(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        world = SimpleWorld(Cessna172FBWMCS()) |> System;

        mid_cg_pld = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

        trim_params = C172.TrimParameters(
        Ob = Geographic(LatLon(), HOrth(1000)),
        EAS = 55.0,
        γ_wOb_n = 0.0,
        x_fuel = 0.5,
        flaps = 0.0,
        payload = mid_cg_pld)

        exit_flag, trim_state = trim!(world, trim_params)
        @test exit_flag === true

        sys_io! = let

            function (world)
                t = world.t[]
            end
        end

        # sim = Simulation(world; dt = 0.01, Δt = 0.01, t_end = 60, sys_io!, adaptive = false)
        sim = Simulation(world; dt = 0.01, Δt = 0.01, t_end = 60, adaptive = false)
        Sim.run!(sim, verbose = true)

        # plots = make_plots(sim; Plotting.defaults...)
        kin_plots = make_plots(TimeHistory(sim).ac.physics.kinematics; Plotting.defaults...)
        air_plots = make_plots(TimeHistory(sim).ac.physics.air; Plotting.defaults...)
        rb_plots = make_plots(TimeHistory(sim).ac.physics.rigidbody; Plotting.defaults...)
        save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172fbw_mcs", "sim", "kin"))
        save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172fbw_mcs", "sim", "air"))
        save && save_plots(rb_plots, save_folder = joinpath("tmp", "test_c172fbw_mcs", "sim", "rigidbody"))

        return nothing

    end
end


function test_sim_paced(; save::Bool = true)

    h_trn = HOrth(601.55);

    ac = Cessna172FBWMCS();
    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
    world = SimpleWorld(ac, env) |> System;

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.3),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0);

    init_kinematics!(world, kin_init)

    sim = Simulation(world; dt = 0.01, Δt = 0.01, t_end = 600)

    interfaces = Vector{IODevices.Interface}()
    for joystick in get_connected_joysticks()
        push!(interfaces, attach_io!(sim, joystick))
    end

    xp = XPCDevice()
    # xp = XPCDevice(host = IPv4("192.168.1.2"))
    push!(interfaces, attach_io!(sim, xp))

    @sync begin
        for interface in interfaces
            Threads.@spawn IODevices.start!(interface)
        end
        Threads.@spawn Sim.run_paced!(sim; pace = 1, verbose = true)
    end

    kin_plots = make_plots(TimeHistory(sim).ac.physics.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeHistory(sim).ac.physics.air; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172fbw_mcs", "sim_paced", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172fbw_mcs", "sim_paced", "air"))

    return nothing

end


end #module