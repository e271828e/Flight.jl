module TestC172FBWCAS

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

using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172FBWCAS

export test_c172fbw_cas


function test_c172fbw_cas()
    @testset verbose = true "Cessna172FBW CAS" begin

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

            ac = System(Cessna172FBWCAS());

            #to exercise all airframe functionality, including landing gear, we
            #need to be on the ground with the engine running
            init_kinematics!(ac, kin_init_gnd)
            ac.avionics.u.physical.eng_start = true #engine start switch on
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
            u_physical = ac.avionics.u.physical
            u_digital = ac.avionics.u.digital

            #we start by testing semiautomatic modes. enabling the outermost
            #loop for each control channel will prompt execution of all inner
            #control loops
            u_digital.lon_mode_sel = C172FBWCAS.lon_mode_semi
            u_digital.lat_mode_sel = C172FBWCAS.lat_mode_semi

            u_digital.throttle_mode_sel = C172FBWCAS.airspeed_throttle_mode
            u_digital.roll_mode_sel = C172FBWCAS.bank_angle_mode
            u_digital.pitch_mode_sel = C172FBWCAS.climb_rate_mode
            u_digital.yaw_mode_sel = C172FBWCAS.sideslip_mode
            u_digital.TAS_dmd = 40
            u_digital.φ_dmd = 0.1
            u_digital.c_dmd = 1
            u_physical.yaw_input = 0.02

            f_disc!(ac, 0.02, env)

            y_mod = ac.avionics.y.moding
            y_act = ac.avionics.y.actuation

            @test y_mod.flight_phase == C172FBWCAS.phase_air
            @test y_mod.throttle_mode === C172FBWCAS.airspeed_throttle_mode
            @test y_mod.roll_mode === C172FBWCAS.bank_angle_mode
            @test y_mod.pitch_mode === C172FBWCAS.climb_rate_mode
            @test y_mod.yaw_mode === C172FBWCAS.sideslip_mode

            #the demands should have propagated through the control loops to the
            #actuators
            @test Float64(y_act.throttle_cmd) > 0
            @test Float64(y_act.aileron_cmd) > 0
            @test Float64(y_act.elevator_cmd) > 0
            @test Float64(y_act.rudder_cmd) > 0

            #now all outermost loops are active, test for allocations
            @test @ballocated(f_disc!($ac, 0.02, $env)) == 0

        end

    return nothing

end

function test_sim(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        h_trn = HOrth(608.55);

        world = SimpleWorld(Cessna172FBWCAS()) |> System;

        u_pld = world.ac.physics.airframe.pld.u
        u_pld.m_pilot = 75
        u_pld.m_copilot = 75
        u_pld.m_lpass = 0
        u_pld.m_rpass = 0
        u_pld.m_baggage = 100

        world.env.atm.wind.u.v_ew_n .= [0, 0, 0]

        trim_params = C172FBW.TrimParameters(
            Ob = Geographic(LatLon(), h_trn + 1000),
            ψ_nb = 0.0,
            TAS = 60.0,
            γ_wOb_n = 0.0,
            ψ_lb_dot = 0.0,
            θ_lb_dot = 0.0,
            β_a = 0.0,
            fuel_load = 1,
            mixture = 0.5,
            flaps = 0.0)

        exit_flag, trim_state = trim!(world, trim_params)
        @test exit_flag === true

        # return world.ac.avionics.y.actuation

        sys_io! = let

            function (world)

                t = world.t[]

                u_physical = world.ac.avionics.u.physical
                u_digital = world.ac.avionics.u.digital

                if 4 < t < 5
                    # world.env.atm.wind.u.v_ew_n[3] = 5
                    u_digital.pitch_mode_sel = C172FBWCAS.pitch_rate_mode
                    # u_digital.pitch_mode_sel = C172FBWCAS.pitch_angle_mode
                    u_physical.pitch_input = 0.5
                    # u_digital.θ_dmd = 0.1
                    # u_digital.c_dmd = π/6

                    # u_digital.throttle_mode_sel = C172FBWCAS.direct_throttle_mode
                    # u_digital.throttle_mode_sel = C172FBWCAS.airspeed_throttle_mode
                    # u_physical.throttle = 1
                    # u_digital.TAS_dmd = 50

                    # u_digital.roll_mode_sel = C172FBWCAS.course_angle_mode
                    # u_digital.roll_mode_sel = C172FBWCAS.bank_angle_mode
                    # u_physical.roll_input = 0.1
                    # u_digital.φ_dmd = 0
                    # u_digital.χ_dmd = π

                    # u_digital.yaw_mode_sel = C172FBWCAS.sideslip_mode
                    # u_physical.yaw_input = 0
                else
                    u_physical.pitch_input = 0

                end
            end
        end

        sim = Simulation(world; dt = 0.01, Δt = 0.01, t_end = 10, sys_io!, adaptive = false)
        Sim.run!(sim, verbose = true)

        # plots = make_plots(sim; Plotting.defaults...)
        kin_plots = make_plots(TimeHistory(sim).ac.physics.kinematics; Plotting.defaults...)
        air_plots = make_plots(TimeHistory(sim).ac.physics.air; Plotting.defaults...)
        save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172fbw_cas", "sim", "kin"))
        save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172fbw_cas", "sim", "air"))

        return nothing

    end
end


function test_sim_paced(; save::Bool = true)

    h_trn = HOrth(601.55);

    ac = Cessna172FBWCAS();
    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
    world = SimpleWorld(ac, env) |> System;

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.3),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0);

    init_kinematics!(world, kin_init)

    sim = Simulation(world; dt = 0.01, Δt = 0.02, t_end = 300)

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
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172fbw_cas", "sim_paced", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172fbw_cas", "sim_paced", "air"))

    return nothing

end


end #module