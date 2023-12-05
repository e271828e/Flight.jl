module TestC172MCS

using Test, UnPack, BenchmarkTools, Sockets

using Flight.FlightCore
using Flight.FlightCore.Sim
using Flight.FlightCore.Visualization

using Flight.FlightPhysics
using Flight.FlightComponents

using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172MCS

export test_c172_mcs


function test_c172_mcs()
    @testset verbose = true "Cessna172 MCS" begin

        test_system_methods()
        test_sim(save = false)

    end
end

function terrain_init()

    #HorizontalTerrain altitude is defined as HOrth. However, KinematicInit
    #works internally with HEllip. Because the difference between orthometric
    #and ellipsoidal altitude depends on geographic location, if we create two
    #different kinematic initializers with the same terrain-referenced
    #orthometric altitudes but at different locations, we will get two different
    #ellipsoidal altitudes. However, this should lead to the same
    #terrain-relative heights and therefore to the same behavior.

    h_trn = HOrth(0)
    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172MCS(LTF(), trn) |> System;

    kin_init_1 = KinematicInit(
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9);

    kin_init_2 = KinematicInit( h = TerrainData(trn).altitude + 1.9);

    @show kin_init_1
    @show kin_init_2

    Systems.init!(ac, kin_init_1)
    sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 600)
    step!(sim, 10, true)
    kin_plots = make_plots(TimeHistory(sim).physics.kinematics; Plotting.defaults...)
    save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172_mcs", "trn_init", "1"))

    Systems.init!(ac, kin_init_2)
    sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 600)
    step!(sim, 10, true)
    kin_plots = make_plots(TimeHistory(sim).physics.kinematics; Plotting.defaults...)
    save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172_mcs", "trn_init", "2"))

end

function test_avionics()

    data_folder = joinpath(dirname(dirname(@__DIR__)),
        normpath("src/aircraft/c172/c172fbw/variants/mcs/data"))

    y_kin(ac::System{<:Cessna172MCS}) = ac.y.physics.kinematics
    y_air(ac::System{<:Cessna172MCS}) = ac.y.physics.air

    @testset verbose = true "Avionics" begin

    h_trn = HOrth(0)
    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172MCS(LTF(), trn) |> System;
    design_point = C172.TrimParameters()
    av = ac.avionics

    ############################################################################
    ############################## Ground ######################################

    @testset verbose = true "Ground" begin

    #we don't really need to provide a specific sys_init! function, because
    #sys_init! defaults to Systems.init!, which for Aircraft has methods
    #accepting both a Kinematics.Initializer and an AbstractTrimParameters
    sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 600)
    kin_init_gnd = KinematicInit( h = TerrainData(trn).altitude + 1.9);
    reinit!(sim, kin_init_gnd)

    step!(sim, 1, true)

    @test ac.avionics.y.flight_phase === C172MCS.phase_gnd

    #set arbitrary control and guidance modes
    av.u.ver_gdc_mode_req = C172MCS.ver_gdc_alt
    av.u.hor_gdc_mode_req = C172MCS.hor_gdc_line
    av.u.lon_ctl_mode_req = C172MCS.lon_EAS_clm
    av.u.lat_ctl_mode_req = C172MCS.lat_p_β
    av.u.throttle_input = 0.1
    av.u.roll_input = 0.2
    av.u.pitch_input = 0.3
    av.u.yaw_input = 0.4
    step!(sim, 1, true)

    @test av.y.flight.ver_gdc_mode === C172MCS.ver_gdc_off
    @test av.y.flight.hor_gdc_mode === C172MCS.hor_gdc_off
    @test av.y.flight.lon_ctl_mode === C172MCS.lon_SAS_off
    @test av.y.flight.lat_ctl_mode === C172MCS.lat_SAS_off
    @test av.y.actuation.throttle_cmd == 0.1
    @test av.y.actuation.aileron_cmd == 0.2
    @test av.y.actuation.elevator_cmd == 0.3
    @test av.y.actuation.rudder_cmd == 0.4

    # @test @ballocated(f_ode!($ac)) == 0
    # @test @ballocated(f_step!($ac)) == 0
    # @test @ballocated(f_disc!($ac, 0.01)) == 0

    end #testset

    ############################################################################
    ################################# Air ######################################

    @testset verbose = true "Air" begin

    #create a simulation with trim! as System reinitializer
    sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 600)

    #put the aircraft in its nominal design point
    reinit!(sim, design_point)
    y_kin_trim = y_kin(ac)

    ############################### direct control #############################
    reinit!(sim, design_point)
    step!(sim, 0.01, true)
    @test av.y.flight.lon_ctl_mode === C172MCS.lon_SAS_off
    @test av.y.flight.lat_ctl_mode === C172MCS.lat_SAS_off

    #with direct surface control, trim state must be initially preserved
    step!(sim, 10, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b, y_kin_trim.ω_lb_b; atol = 1e-5))
    @test all(isapprox.(y_kin(ac).v_eOb_b, y_kin_trim.v_eOb_b; atol = 1e-2))

    # @test @ballocated(f_disc!($ac, 0.01)) == 0

    ############################ thr+ele SAS mode ##############################
    #we test the longitudinal SAS first, because we want to test the lateral
    #modes with it enabled
    reinit!(sim, design_point)
    av.u.lon_ctl_mode_req = C172MCS.lon_thr_ele
    step!(sim, 0.01, true)
    @test av.y.flight.lon_ctl_mode === C172MCS.lon_thr_ele

    #check the correct parameters are loaded and assigned to the controller
    te2te_lookup = C172MCS.load_lqr_tracker_lookup(joinpath(data_folder, "te2te_lookup.h5"))
    C_fwd = te2te_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
    @test all(isapprox.(av.y.flight.lon_ctl.te2te_lqr.C_fwd, C_fwd; atol = 1e-6))

    #with thr+ele SAS active, trim state must be preserved for longer
    step!(sim, 30, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b[2], y_kin_trim.ω_lb_b[2]; atol = 1e-5))
    @test all(isapprox.(y_kin(ac).v_eOb_b[1], y_kin_trim.v_eOb_b[1]; atol = 1e-2))

    # @test @ballocated(f_disc!($ac, 0.01)) == 0

    ################################ φ + β #####################################
    reinit!(sim, design_point)
    av.u.lon_ctl_mode_req = C172MCS.lon_thr_ele
    av.u.lat_ctl_mode_req = C172MCS.lat_φ_β
    step!(sim, 0.01, true)
    @test av.y.flight.lat_ctl_mode === C172MCS.lat_φ_β

    #check the correct parameters are loaded and assigned to the controller
    φβ2ar_lookup = C172MCS.load_lqr_tracker_lookup(joinpath(data_folder, "φβ2ar_lookup.h5"))
    C_fwd = φβ2ar_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
    @test all(isapprox.(av.y.flight.lat_ctl.φβ2ar_lqr.C_fwd, C_fwd; atol = 1e-6))

    #with setpoints matching their trim values, the control mode must activate
    #without transients
    step!(sim, 1, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b[2], y_kin_trim.ω_lb_b[2]; atol = 1e-5))
    @test all(isapprox.(y_kin(ac).v_eOb_b[1], y_kin_trim.v_eOb_b[1]; atol = 1e-2))

    #correct tracking
    av.u.φ_sp = 0.1π
    av.u.β_sf = 1.0
    av.u.yaw_input = 0.05
    step!(sim, 5, true)
    @test av.flight.lat_ctl.u.φ_sp != 0
    @test av.flight.lat_ctl.u.β_sp != 0
    @test isapprox(av.u.φ_sp, y_kin(ac).e_nb.φ; atol = 1e-3)
    @test isapprox(Float64(av.u.yaw_input), y_air(ac).β_b; atol = 1e-3)

    # @test @ballocated(f_disc!($ac, 0.01)) == 0

    ################################ p + β ######################################
    reinit!(sim, design_point)
    av.u.lon_ctl_mode_req = C172MCS.lon_thr_ele
    av.u.lat_ctl_mode_req = C172MCS.lat_p_β
    step!(sim, 0.01, true)
    @test av.y.flight.lat_ctl_mode === C172MCS.lat_p_β

    #check the correct parameters are loaded and assigned to the controllers
    φβ2ar_lookup = C172MCS.load_lqr_tracker_lookup(joinpath(data_folder, "φβ2ar_lookup.h5"))
    C_fwd = φβ2ar_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
    @test all(isapprox.(av.y.flight.lat_ctl.φβ2ar_lqr.C_fwd, C_fwd; atol = 1e-6))

    #the control mode must activate without transients
    step!(sim, 1, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b[2], y_kin_trim.ω_lb_b[2]; atol = 1e-5))
    @test all(isapprox.(y_kin(ac).v_eOb_b[1], y_kin_trim.v_eOb_b[1]; atol = 1e-2))

    #the controller must keep trim values in steady state
    step!(sim, 10, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b[2], y_kin_trim.ω_lb_b[2]; atol = 1e-5))
    @test all(isapprox.(y_kin(ac).v_eOb_b[1], y_kin_trim.v_eOb_b[1]; atol = 1e-2))

    av.u.p_sf = 1.0
    av.u.β_sf = 1.0
    av.u.roll_input = 0.01
    av.u.yaw_input = 0.05
    step!(sim, 10, true)
    @test av.flight.lat_ctl.u.p_sp != 0
    @test av.flight.lat_ctl.u.β_sp != 0
    @test isapprox(Float64(av.u.roll_input), y_kin(ac).ω_lb_b[1]; atol = 1e-3)
    @test isapprox(Float64(av.u.yaw_input), y_air(ac).β_b; atol = 1e-3)

    # @test @ballocated(f_disc!($ac, 0.01)) == 0

    ################################ χ + β #####################################
    reinit!(sim, design_point)
    av.u.lon_ctl_mode_req = C172MCS.lon_thr_ele
    av.u.lat_ctl_mode_req = C172MCS.lat_χ_β
    step!(sim, 0.01, true)
    @test av.y.flight.lat_ctl_mode === C172MCS.lat_χ_β

    #check the correct parameters are loaded and assigned to the controller
    χ2φ_lookup = C172MCS.load_pid_lookup(joinpath(data_folder, "χ2φ_lookup.h5"))
    k_p = χ2φ_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
    @test all(isapprox.(av.y.flight.lat_ctl.χ2φ_pid.k_p, k_p; atol = 1e-6))

    #with setpoints matching their trim values, the control mode must activate
    #without transients
    step!(sim, 1, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b[2], y_kin_trim.ω_lb_b[2]; atol = 1e-5))
    @test all(isapprox.(y_kin(ac).v_eOb_b[1], y_kin_trim.v_eOb_b[1]; atol = 1e-2))

    #correct tracking
    av.u.χ_sp = π/2
    av.u.β_sf = 1.0
    av.u.yaw_input = 0.0
    step!(sim, 29, true)
    @test av.flight.lat_ctl.u.χ_sp != 0
    @test isapprox(av.u.χ_sp, y_kin(ac).χ_gnd; atol = 1e-2)
    # @test isapprox(Float64(av.u.yaw_input), y_air(ac).β_b; atol = 1e-3)

    #correct tracking with 10m/s of crosswind
    ac.physics.atmosphere.u.v_ew_n[1] = 10
    step!(sim, 10, true)
    @test isapprox(av.u.χ_sp, y_kin(ac).χ_gnd; atol = 1e-2)
    ac.physics.atmosphere.u.v_ew_n[1] = 0

    # @test @ballocated(f_disc!($ac, 0.01)) == 0


    ############################################################################
    #the rest of longitudinal modes we test with lateral SAS enabled

    ################################# thr+q ####################################
    reinit!(sim, design_point)

    av.u.lon_ctl_mode_req = C172MCS.lon_thr_q
    av.u.lat_ctl_mode_req = C172MCS.lat_p_β
    step!(sim, 0.01, true)
    @test av.y.flight.lon_ctl_mode === C172MCS.lon_thr_q

    #check the correct parameters are loaded and assigned to the controller
    q2e_lookup = C172MCS.load_pid_lookup(joinpath(data_folder, "q2e_lookup.h5"))
    k_p = q2e_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
    @test all(isapprox.(av.y.flight.lon_ctl.q2e_pid.k_p, k_p; atol = 1e-6))

    #when trim setpoints are kept, the control mode must activate without
    #transients
    step!(sim, 1, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b[2], y_kin_trim.ω_lb_b[2]; atol = 1e-5))
    @test all(isapprox.(y_kin(ac).v_eOb_b[1], y_kin_trim.v_eOb_b[1]; atol = 1e-2))

    #correct tracking
    av.u.q_sf = 1.0
    av.u.pitch_input = 0.01
    av.u.throttle_input = 1
    step!(sim, 10, true)

    @test av.flight.lon_ctl.u.q_sp != 0
    @test isapprox(av.flight.lon_ctl.u.q_sp, y_kin(ac).ω_lb_b[2]; atol = 1e-3)

    #note: throttle_cmd != throttle_input, because we have a SAS in between!

    # @test @ballocated(f_disc!($ac, 0.01)) == 0

    ################################ EAS + q ###################################
    reinit!(sim, design_point)

    av.u.lon_ctl_mode_req = C172MCS.lon_EAS_q
    av.u.lat_ctl_mode_req = C172MCS.lat_p_β
    step!(sim, 0.01, true)
    @test av.y.flight.lon_ctl_mode === C172MCS.lon_EAS_q

    #check the correct parameters are loaded and assigned to the controller
    v2t_lookup = C172MCS.load_pid_lookup(joinpath(data_folder, "v2t_lookup.h5"))
    k_p = v2t_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
    @test all(isapprox.(av.y.flight.lon_ctl.v2t_pid.k_p, k_p; atol = 1e-6))

    #when trim setpoints are kept, the control mode must activate without
    #transients
    step!(sim, 1, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b[2], y_kin_trim.ω_lb_b[2]; atol = 1e-5))
    @test all(isapprox.(y_kin(ac).v_eOb_b[1], y_kin_trim.v_eOb_b[1]; atol = 1e-2))

    #correct tracking
    av.u.pitch_input = 0.0
    av.u.EAS_sp = 45
    step!(sim, 60, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b[2], y_kin_trim.ω_lb_b[2]; atol = 1e-5))
    @test all(isapprox.(y_air(ac).EAS, av.u.EAS_sp; atol = 1e-1))

    # @test @ballocated(f_disc!($ac, 0.01)) == 0


    ############################ EAS + climb rate ##############################
    reinit!(sim, design_point)

    av.u.lon_ctl_mode_req = C172MCS.lon_EAS_clm
    av.u.lat_ctl_mode_req = C172MCS.lat_p_β
    step!(sim, 0.01, true)
    @test av.y.flight.lon_ctl_mode === C172MCS.lon_EAS_clm

    #check the correct parameters are loaded and assigned to the controller
    vc2te_lookup = C172MCS.load_lqr_tracker_lookup(joinpath(data_folder, "vc2te_lookup.h5"))
    C_fwd = vc2te_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
    @test all(isapprox.(av.y.flight.lon_ctl.vc2te_lqr.C_fwd, C_fwd; atol = 1e-6))

    #when trim setpoints are kept, the control mode must activate without
    #transients
    step!(sim, 1, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b[2], y_kin_trim.ω_lb_b[2]; atol = 1e-5))
    @test all(isapprox.(y_kin(ac).v_eOb_b[1], y_kin_trim.v_eOb_b[1]; atol = 1e-2))

    #correct tracking
    av.u.EAS_sp = 45
    av.u.clm_sp = 2
    step!(sim, 20, true)
    @test all(isapprox.(y_kin(ac).v_eOb_n[3], -av.u.clm_sp; atol = 1e-2))
    @test all(isapprox.(y_air(ac).EAS, av.u.EAS_sp; atol = 1e-1))

    # @test @ballocated(f_disc!($ac, 0.01)) == 0


    ############################## EAS + throttle ##############################
    reinit!(sim, design_point)

    av.u.lon_ctl_mode_req = C172MCS.lon_EAS_thr
    av.u.lat_ctl_mode_req = C172MCS.lat_p_β
    step!(sim, 0.01, true)
    @test av.y.flight.lon_ctl_mode === C172MCS.lon_EAS_thr

    #check the correct parameters are loaded and assigned to the controller
    vt2te_lookup = C172MCS.load_lqr_tracker_lookup(joinpath(data_folder, "vt2te_lookup.h5"))
    C_fwd = vt2te_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
    @test all(isapprox.(av.y.flight.lon_ctl.vt2te_lqr.C_fwd, C_fwd; atol = 1e-6))

    #when trim setpoints are kept, the control mode must activate without
    #transients
    step!(sim, 1, true)
    @test all(isapprox.(y_kin(ac).ω_lb_b[2], y_kin_trim.ω_lb_b[2]; atol = 1e-5))
    @test all(isapprox.(y_kin(ac).v_eOb_b[1], y_kin_trim.v_eOb_b[1]; atol = 1e-2))

    #correct tracking
    av.u.EAS_sp = 45
    av.u.throttle_input = 1
    step!(sim, 30, true)
    @test all(isapprox.(ac.y.physics.airframe.act.throttle_cmd, av.flight.u.throttle_sp; atol = 1e-2))
    @test all(isapprox.(y_air(ac).EAS, av.u.EAS_sp; atol = 1e-1))

    return
    # @test @ballocated(f_disc!($ac, 0.01)) == 0


    kin_plots = make_plots(TimeHistory(sim).physics.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeHistory(sim).physics.air; Plotting.defaults...)
    save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172_mcs", "avionics", "kin"))
    save_plots(air_plots, save_folder = joinpath("tmp", "test_c172_mcs", "avionics", "air"))

    #note: throttle_cmd != throttle_input, because we have a SAS in between!

    # @test @ballocated(f_disc!($ac, 0.01)) == 0



    end #testset

    end #testset

    return


end


function test_cas(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        ac = Cessna172MCS() |> System;
        design_condition = C172.TrimParameters()

        exit_flag, trim_state = trim!(ac, design_condition)

        @test exit_flag === true

        sys_io! = let

            function (ac)

                t = ac.t[]

                u_inceptors = ac.avionics.u.inceptors
                u_digital = ac.avionics.u.digital

                u_digital.lon_mode_sel = C172MCS.lon_θ_EAS

                # u_digital.lon_mode_sel = C172MCS.lon_q_EAS
                # if 5 < t < 15
                #     u_inceptors.pitch_input = 0.001
                # elseif 15 < t < 25
                #     u_inceptors.pitch_input = -0.001
                # else
                #     u_inceptors.pitch_input = 0
                # end

            end
        end

        sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 30, sys_io!, adaptive = false)
        # sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 60, adaptive = false)
        Sim.run!(sim, verbose = true)

        # plots = make_plots(sim; Plotting.defaults...)
        kin_plots = make_plots(TimeHistory(sim).physics.kinematics; Plotting.defaults...)
        air_plots = make_plots(TimeHistory(sim).physics.air; Plotting.defaults...)
        rb_plots = make_plots(TimeHistory(sim).physics.rigidbody; Plotting.defaults...)
        save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172_mcs", "cas", "kin"))
        save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172_mcs", "cas", "air"))
        save && save_plots(rb_plots, save_folder = joinpath("tmp", "test_c172_mcs", "cas", "rigidbody"))

        return nothing

    end
end




end #module