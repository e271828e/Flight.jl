module TestC172CAS

using Test, UnPack, BenchmarkTools, Sockets

using Flight.FlightCore
using Flight.FlightCore.Sim
using Flight.FlightCore.Visualization

using Flight.FlightPhysics

using Flight.FlightComponents

using Flight.FlightAircraft.AircraftBase
using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172CAS

export test_c172cas


function test_c172cas()
    @testset verbose = true "Cessna172FBW CAS" begin

        test_system_methods()
        test_sim(save = false)

    end
end

function test_system_methods()

        @testset verbose = true "System Methods" begin

            trn = HorizontalTerrain()
            loc = NVector()
            trn_data = TerrainData(trn, loc)
            kin_init_gnd = KinematicInit( h = trn_data.altitude + 1.8);
            kin_init_air = KinematicInit( h = trn_data.altitude + 1000);

            ac = System(Cessna172CAS());

            #to exercise all airframe functionality, including landing gear, we
            #need to be on the ground with the engine running
            Systems.init!(ac, kin_init_gnd)
            ac.avionics.u.inceptors.eng_start = true #engine start switch on
            f_disc!(ac, 0.02) #run avionics for the engine start signal to propagate
            f_ode!(ac)
            f_step!(ac)
            f_ode!(ac)
            f_step!(ac)
            @test ac.y.physics.airframe.ldg.left.strut.wow == true
            @test ac.y.physics.airframe.pwp.engine.state === Piston.eng_starting

            @test @ballocated(f_ode!($ac)) == 0
            @test @ballocated(f_step!($ac)) == 0
            @test @ballocated(f_disc!($ac, 0.02)) == 0

            #now we put the aircraft in flight
            Systems.init!(ac, kin_init_air)
            f_ode!(ac)
            @test ac.y.physics.airframe.ldg.left.strut.wow == false
            @test @ballocated(f_ode!($ac)) == 0
            @test @ballocated(f_step!($ac)) == 0

            #testing the different avionics modes for allocations is a bit more
            #involved
            u_inceptors = ac.avionics.u.inceptors
            u_digital = ac.avionics.u.digital

            #we start by testing semiautomatic modes. enabling the outermost
            #loop for each control channel will prompt execution of all inner
            #control loops
            u_digital.lon_mode_sel = C172CAS.lon_mode_semi
            u_digital.lat_mode_sel = C172CAS.lat_mode_semi

            u_digital.throttle_mode_sel = C172CAS.EAS_throttle_mode
            u_digital.roll_mode_sel = C172CAS.course_angle_mode
            u_digital.pitch_mode_sel = C172CAS.climb_rate_mode
            # u_digital.yaw_mode_sel = C172CAS.sideslip_mode
            u_digital.EAS_dmd = 40
            u_digital.χ_dmd = 0.1
            u_digital.c_dmd = 1
            # u_inceptors.yaw_input = 0.02

            f_disc!(ac, 0.02)

            y_mod = ac.avionics.y.moding

            @test y_mod.flight_phase == C172CAS.phase_air
            @test y_mod.throttle_mode === C172CAS.EAS_throttle_mode
            @test y_mod.roll_mode === C172CAS.course_angle_mode
            @test y_mod.pitch_mode === C172CAS.climb_rate_mode
            # @test y_mod.yaw_mode === C172CAS.sideslip_mode

            #the demands should have propagated through the control loops
            @test Float64(ac.y.avionics.throttle_ctl.thr_cmd) > 0
            @test Float64(ac.y.avionics.roll_ctl.a_cmd) > 0
            @test Float64(ac.y.avionics.pitch_ctl.e_cmd) > 0
            # @test Float64(y_act.rudder_cmd) > 0

            #now all outermost loops are active, test for allocations
            @test @ballocated(f_disc!($ac, 0.02)) == 0

        end

    return nothing

end

function test_sim(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        ac = Cessna172CAS() |> System;

        mid_cg_pld = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

        ac.physics.atmosphere.u.v_ew_n .= [0, 0, 0]

        trim_params = C172.TrimParameters(
        Ob = Geographic(LatLon(), HOrth(1000)),
        EAS = 55.0,
        γ_wOb_n = 0.0,
        x_fuel = 0.5,
        flaps = 0.0,
        payload = mid_cg_pld)

        exit_flag, trim_state = trim!(ac, trim_params)
        @test exit_flag === true

        sys_io! = let

            function (ac)

                t = ac.t[]

                u_inceptors = ac.avionics.u.inceptors
                u_digital = ac.avionics.u.digital

                # @show ac.avionics.y.actuation.throttle_cmd

                # u_digital.throttle_mode_sel = C172CAS.direct_throttle_mode
                # u_digital.throttle_mode_sel = C172CAS.EAS_throttle_mode
                # u_inceptors.throttle = 1
                # u_digital.EAS_dmd = 25

                # u_digital.pitch_mode_sel = C172CAS.pitch_rate_mode
                # u_digital.pitch_mode_sel = C172CAS.pitch_angle_mode
                # u_digital.pitch_mode_sel = C172CAS.climb_rate_mode
                # u_digital.pitch_mode_sel = C172CAS.EAS_pitch_mode
                # @show ac.avionics.y.moding.pitch_mode
                # u_inceptors.pitch_input = 0.005
                # u_digital.θ_dmd = 0.01π
                u_digital.c_dmd = 0.0

                u_digital.EAS_dmd = 30.0
                u_inceptors.flaps = 1
                u_digital.lon_mode_sel = C172CAS.lon_mode_alt
                u_digital.h_dmd = 1000.0
                u_digital.h_ref = C172CAS.orthometric

                # u_digital.roll_mode_sel = C172CAS.roll_rate_mode
                u_digital.roll_mode_sel = C172CAS.bank_angle_mode
                # u_digital.roll_mode_sel = C172CAS.course_angle_mode
                # u_inceptors.roll_input = 0.0
                u_digital.φ_dmd = π/6
                # u_digital.χ_dmd = π

                # u_digital.yaw_mode_sel = C172CAS.sideslip_mode
                # u_digital.yaw_mode_sel = C172CAS.direct_rudder_mode
                # u_inceptors.yaw_input = 0.1


                if 0 < t <= 10
                    # u_inceptors.roll_input = .0
                    # u_inceptors.pitch_input = .0
                    # u_inceptors.yaw_input = .01
                elseif 10 < t < 45
                    # u_inceptors.roll_input = 0.0
                    # u_inceptors.pitch_input = 0.0
                    # u_inceptors.yaw_input = 1
                    # u_digital.EAS_dmd = 45
                else
                    # f = 0.5
                    # u_digital.EAS_dmd = 50
                end
            end
        end

        sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 60, sys_io!, adaptive = false)
        # sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 60, adaptive = false)
        Sim.run!(sim)

        # plots = make_plots(sim; Plotting.defaults...)
        kin_plots = make_plots(TimeSeries(sim).physics.kinematics; Plotting.defaults...)
        air_plots = make_plots(TimeSeries(sim).physics.air; Plotting.defaults...)
        rb_plots = make_plots(TimeSeries(sim).physics.dynamics; Plotting.defaults...)
        save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172cas", "sim", "kin"))
        save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172cas", "sim", "air"))
        save && save_plots(rb_plots, save_folder = joinpath("tmp", "test_c172cas", "sim", "dynamics"))

        return nothing

    end
end


function test_sim_paced(; save::Bool = true)

    h_trn = HOrth(601.55);

    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172CAS(LTF(), trn) |> System;

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.0),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0);

    Systems.init!(ac, kin_init)

    sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 600)

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
        Threads.@spawn Sim.run_paced!(sim; pace = 1)
    end

    kin_plots = make_plots(TimeSeries(sim).physics.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).physics.air; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172cas", "sim_paced", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172cas", "sim_paced", "air"))

    return nothing

end

end #module