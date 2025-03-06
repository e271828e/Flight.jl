module TestC172FBWv1

using Test, UnPack, BenchmarkTools, Sockets

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

#non-exported stuff
using Flight.FlightLib.Control.Discrete: load_pid_lookup, load_lqr_tracker_lookup
using Flight.FlightAircraft.C172FBW.C172FBWControl: lon_direct, lon_thr_ele, lon_thr_q, lon_thr_θ, lon_thr_EAS, lon_EAS_q, lon_EAS_θ, lon_EAS_clm
using Flight.FlightAircraft.C172FBW.C172FBWControl: lat_direct, lat_p_β, lat_φ_β, lat_χ_β
using Flight.FlightAircraft.C172FBW.C172FBWControl: vrt_gdc_off, vrt_gdc_alt
using Flight.FlightAircraft.C172FBW.C172FBWControl: hor_gdc_off, hor_gdc_line
using Flight.FlightAircraft.C172FBW.C172FBWControl: phase_gnd, phase_air

export test_c172fbw_v1


function test_c172fbw_v1()
    @testset verbose = true "Cessna172 FBWv1" begin

        test_control_modes()
        test_guidance_modes()

    end
end

y_kin(ac::System{<:Cessna172FBWv1}) = ac.y.vehicle.kinematics
y_air(ac::System{<:Cessna172FBWv1}) = ac.y.vehicle.air
y_aero(ac::System{<:Cessna172FBWv1}) = ac.y.vehicle.components.aero

function test_control_modes()

    data_folder = joinpath(dirname(dirname(@__DIR__)),
        normpath("src/aircraft/c172/c172fbw/control/data"))

    @testset verbose = true "Control Modes" begin

    h_trn = HOrth(0)
    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172FBWv1(WA(), trn) |> System;
    ctl = ac.avionics.ctl

    init_gnd = KinInit( h = TerrainData(trn).altitude + 1.9);
    init_air = C172.TrimParameters()

    #we don't really need to provide a specific sys_init! function, because
    #sys_init! defaults to Systems.init!, which for Aircraft has methods
    #accepting both a Kinematics.Initializer and an AbstractTrimParameters
    dt = Δt = 0.01
    sim = Simulation(ac; dt, Δt, t_end = 600)

    ############################################################################
    ############################## Ground ######################################

    @testset verbose = true "Ground" begin

    reinit!(sim, init_gnd)

    @test ac.y.avionics.ctl.flight_phase === phase_gnd

    #set arbitrary control and guidance modes
    ctl.u.vrt_gdc_mode_req = vrt_gdc_alt
    ctl.u.hor_gdc_mode_req = hor_gdc_line
    ctl.u.lon_ctl_mode_req = lon_EAS_clm
    ctl.u.lat_ctl_mode_req = lat_p_β
    ctl.u.throttle_sp_input = 0.1
    ctl.u.aileron_sp_input = 0.2
    ctl.u.elevator_sp_input = 0.3
    ctl.u.rudder_sp_input = 0.4

    step!(sim, ctl.Δt, true)

    @test ctl.y.flight_phase === phase_gnd

    #the mode requests are overridden due to phase_gnd
    @test ctl.y.vrt_gdc_mode === vrt_gdc_off
    @test ctl.y.hor_gdc_mode === hor_gdc_off
    @test ctl.y.lon_ctl_mode === lon_direct
    @test ctl.y.lat_ctl_mode === lat_direct

    #control laws outputs must have propagated to actuator inputs (not yet to
    #their outputs, that requires a subsequent call to f_ode!)
    @test ac.vehicle.components.act.throttle.u[] == 0.1
    @test ac.vehicle.components.act.aileron.u[] == 0.2
    @test ac.vehicle.components.act.elevator.u[] == 0.3
    @test ac.vehicle.components.act.rudder.u[] == 0.4

    #must reset scheduling counter before standalone calls to f_disc!, but
    #without calling Sim.reinit! so that the controller state is preserved
    # ac.n[] = 0
    # @test @ballocated(f_ode!($ac)) == 0
    # @test @ballocated(f_step!($ac)) == 0
    # @test @ballocated(f_disc!($ac)) == 0

    end #testset

    ############################################################################
    ################################# Air ######################################

    @testset verbose = true "Air" begin

    #put the aircraft in its nominal design point
    reinit!(sim, init_air)
    y_kin_trim = y_kin(ac)

    ############################### direct control #############################

    @testset verbose = true "lon_direct + lat_direct" begin

        reinit!(sim, init_air)
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_direct
        @test ctl.y.lat_ctl_mode === lat_direct

        #with direct surface control, trim state must be initially preserved
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b, y_kin_trim.ω_wb_b; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b, y_kin_trim.v_eb_b; atol = 1e-2))

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0

    end #testset

    ############################ thr+ele SAS mode ##############################

    @testset verbose = true "lon_thr_ele" begin

        #we test the longitudinal SAS first, because we want to test the lateral
        #modes with it enabled
        reinit!(sim, init_air)
        ctl.u.lon_ctl_mode_req = lon_thr_ele
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_thr_ele

        #check the correct parameters are loaded and assigned to the controller
        te2te_lookup = load_lqr_tracker_lookup(joinpath(data_folder, "te2te_lookup.h5"))
        C_fwd = te2te_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
        @test all(isapprox.(ctl.y.lon_ctl.te2te_lqr.C_fwd, C_fwd; atol = 1e-6))

        #with thr+ele SAS active, trim state must be preserved for longer
        step!(sim, 30, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0

        end #testset

    ################################ φ + β #####################################

    @testset verbose = true "lat_φ_β" begin

        reinit!(sim, init_air)
        ctl.u.lon_ctl_mode_req = lon_thr_ele
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lat_ctl_mode === lat_φ_β

        #check the correct parameters are loaded and assigned to the controller
        φβ2ar_lookup = load_lqr_tracker_lookup(joinpath(data_folder, "φβ2ar_lookup.h5"))
        C_fwd = φβ2ar_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
        @test all(isapprox.(ctl.y.lat_ctl.φβ2ar_lqr.C_fwd, C_fwd; atol = 1e-6))

        #with setpoints matching their trim values, the control mode must activate
        #without transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_sp = π/12
        ctl.u.β_sp = deg2rad(3)
        step!(sim, 10, true)
        @test isapprox(ctl.u.φ_sp, y_kin(ac).e_nb.φ; atol = 1e-3)
        @test isapprox(Float64(ctl.u.β_sp), y_aero(ac).β; atol = 1e-3)

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0

    end

    ################################ p + β #####################################

    @testset verbose = true "lat_p_β" begin

        reinit!(sim, init_air)
        ctl.u.lon_ctl_mode_req = lon_thr_ele
        ctl.u.lat_ctl_mode_req = lat_p_β
        step!(sim, Δt, true)
        @test ctl.y.lat_ctl_mode === lat_p_β

        #check the correct parameters are loaded and assigned to the controllers
        φβ2ar_lookup = load_lqr_tracker_lookup(joinpath(data_folder, "φβ2ar_lookup.h5"))
        C_fwd = φβ2ar_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
        @test all(isapprox.(ctl.y.lat_ctl.φβ2ar_lqr.C_fwd, C_fwd; atol = 1e-6))

        #the control mode must activate without transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #the controller must keep trim values in steady state
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        ctl.u.p_sp = 0.02
        ctl.u.β_sp = deg2rad(3)
        step!(sim, 10, true)
        @test isapprox(Float64(ctl.u.p_sp), y_kin(ac).ω_wb_b[1]; atol = 1e-3)
        @test isapprox(ctl.u.β_sp, y_aero(ac).β; atol = 1e-3)

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0

    end

    ################################ χ + β #####################################

    @testset verbose = true "lat_χ_β" begin

        reinit!(sim, init_air)
        ctl.u.lon_ctl_mode_req = lon_thr_ele
        ctl.u.lat_ctl_mode_req = lat_χ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lat_ctl_mode === lat_χ_β

        #check the correct parameters are loaded and assigned to the controller
        χ2φ_lookup = load_pid_lookup(joinpath(data_folder, "χ2φ_lookup.h5"))
        k_p = χ2φ_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
        @test all(isapprox.(ctl.y.lat_ctl.χ2φ_pid.k_p, k_p; atol = 1e-6))

        #with setpoints matching their trim values, the control mode must activate
        #without transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking
        ctl.u.χ_sp = π/2
        step!(sim, 29, true)
        @test ctl.lat_ctl.u.χ_sp != 0
        @test isapprox(ctl.u.χ_sp, y_kin(ac).χ_gnd; atol = 1e-2)
        # @test isapprox(Float64(ctl.u.yaw_input), y_aero(ac).β; atol = 1e-3)

        #correct tracking with 10m/s of crosswind
        ac.vehicle.atmosphere.u.v_ew_n[1] = 10
        step!(sim, 10, true)
        @test isapprox(ctl.u.χ_sp, y_kin(ac).χ_gnd; atol = 1e-2)
        ac.vehicle.atmosphere.u.v_ew_n[1] = 0

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0

    end

    ############################################################################

    #now we proceed to test the remaining longitudinal modes we test with
    #lateral p + β mode enabled

    ############################### lon_thr_q ##################################

    @testset verbose = true "lon_thr_q" begin

        reinit!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_thr_q
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_thr_q

        #check the correct parameters are loaded and assigned to the controller
        q2e_lookup = load_pid_lookup(joinpath(data_folder, "q2e_lookup.h5"))
        k_p = q2e_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
        @test all(isapprox.(ctl.y.lon_ctl.q2e_pid.k_p, k_p; atol = 1e-6))

        #when trim setpoints are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_sp = π/12
        ctl.u.q_sp = 0.01
        step!(sim, 10, true)

        @test ctl.lon_ctl.u.q_sp != 0
        @test isapprox(ctl.lon_ctl.u.q_sp, y_kin(ac).ω_wb_b[2]; atol = 1e-3)
        @test isapprox(Float64(ac.y.vehicle.components.act.throttle.cmd),
                        Float64(ctl.u.throttle_sp_input + ctl.u.throttle_sp_offset); atol = 1e-3)

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0


    end


    ############################## lon_thr_θ ###################################

    @testset verbose = true "lon_thr_θ" begin

        reinit!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_thr_θ
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_thr_θ

        #when trim setpoints are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_sp = π/6
        ctl.u.θ_sp = deg2rad(5)
        step!(sim, 10, true)
        @test isapprox(y_kin(ac).e_nb.θ, ctl.u.θ_sp; atol = 1e-4)

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0


    end


    ################################ lon_thr_EAS ###############################

    @testset verbose = true "lon_thr_EAS" begin

        reinit!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_thr_EAS
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_thr_EAS

        #check the correct parameters are loaded and assigned to the controller
        v2θ_lookup = load_pid_lookup(joinpath(data_folder, "v2θ_lookup.h5"))
        k_p = v2θ_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
        @test all(isapprox.(ctl.y.lon_ctl.v2θ_pid.k_p, k_p; atol = 1e-6))

        #when trim setpoints are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_sp = π/6
        ctl.u.EAS_sp = 45
        step!(sim, 30, true)
        @test all(isapprox.(y_air(ac).EAS, ctl.u.EAS_sp; atol = 1e-1))

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0


    end

    ################################ lon_EAS_q #################################

    @testset verbose = true "lon_EAS_q" begin

        reinit!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_EAS_q
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_EAS_q

        #check the correct parameters are loaded and assigned to v2t, the q
        #tracker is shared with other modes
        v2t_lookup = load_pid_lookup(joinpath(data_folder, "v2t_lookup.h5"))
        k_p = v2t_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
        @test all(isapprox.(ctl.y.lon_ctl.v2t_pid.k_p, k_p; atol = 1e-6))

        #when trim setpoints are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking
        ctl.u.q_sp = -0.01
        step!(sim, 10, true)
        ctl.u.q_sp = 0.005
        step!(sim, 10, true)
        ctl.u.q_sp = 0.0
        step!(sim, 20, true)

        @test isapprox(ctl.lon_ctl.u.q_sp, y_kin(ac).ω_wb_b[2]; atol = 1e-3)
        @test all(isapprox.(y_air(ac).EAS, ctl.u.EAS_sp; atol = 1e-1))

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0

    end


    ################################ lon_EAS_q #################################

    @testset verbose = true "lon_EAS_θ" begin

        reinit!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_EAS_θ
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_EAS_θ

        #when trim setpoints are kept, the control mode must activate without
        #transients
        step!(sim, 0.1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_sp = π/6
        ctl.u.θ_sp = deg2rad(3)
        step!(sim, 10, true)
        ctl.u.θ_sp = -deg2rad(3)
        step!(sim, 60, true)

        @test isapprox(ctl.lon_ctl.u.θ_sp, y_kin(ac).e_nb.θ; atol = 1e-3)
        @test all(isapprox.(y_air(ac).EAS, ctl.u.EAS_sp; atol = 1e-1))

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0

    end

    ############################## lon_EAS_clm #################################

    @testset verbose = true "lon_EAS_clm" begin

        reinit!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_EAS_clm
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_EAS_clm

        #check the correct parameters are loaded and assigned to the controller
        vc2te_lookup = load_lqr_tracker_lookup(joinpath(data_folder, "vc2te_lookup.h5"))
        C_fwd = vc2te_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
        @test all(isapprox.(ctl.y.lon_ctl.vc2te_lqr.C_fwd, C_fwd; atol = 1e-6))

        #when trim setpoints are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_sp = π/6
        ctl.u.EAS_sp = 45
        ctl.u.clm_sp = 2
        step!(sim, 30, true)
        @test all(isapprox.(y_kin(ac).v_eb_n[3], -ctl.u.clm_sp; atol = 1e-1))
        @test all(isapprox.(y_air(ac).EAS, ctl.u.EAS_sp; atol = 1e-1))

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        # ac.n[] = 0
        # @test @ballocated(f_disc!($ac)) == 0

        # kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
        # air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
        # save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172_fbw_v1", "avionics", "kin"))
        # save_plots(air_plots, save_folder = joinpath("tmp", "test_c172_fbw_v1", "avionics", "air"))
        # return TimeSeries(sim)


    end #testset

    end #testset

    end #testset

end #function


function test_guidance_modes()

    @testset verbose = true "Guidance Modes" begin

    h_trn = HOrth(0)
    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172FBWv1(WA(), trn) |> System;
    ctl = ac.avionics.ctl
    init_air = C172.TrimParameters()

    dt = Δt = 0.01
    sim = Simulation(ac; dt, Δt, t_end = 600)

    @testset verbose = true "Altitude Guidance" begin

        reinit!(sim, init_air)
        y_kin_trim = y_kin(ac)

        ctl.u.vrt_gdc_mode_req = vrt_gdc_alt
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.vrt_gdc_mode === vrt_gdc_alt
        @test ctl.y.lon_ctl_mode === lon_EAS_clm

        #when trim setpoints are kept, the guidance mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #all tests while turning
        ctl.u.φ_sp = π/12

        ctl.u.h_sp = y_kin_trim.h_e + 100
        step!(sim, 1, true)
        @test ctl.y.lon_ctl_mode === lon_thr_EAS
        step!(sim, 60, true) #altitude is captured
        @test ctl.y.lon_ctl_mode === lon_EAS_clm
        @test isapprox.(y_kin(ac).h_e - ctl.u.h_sp, 0.0; atol = 1e-1)

        #setpoint changes within the current threshold do not prompt a mode change
        ctl.u.h_sp = y_kin(ac).h_e - ctl.alt_gdc.s.h_thr / 2
        step!(sim, 1, true)
        @test ctl.y.lon_ctl_mode === lon_EAS_clm
        step!(sim, 30, true) #altitude is captured
        @test isapprox.(y_kin(ac).h_e - ctl.u.h_sp, 0.0; atol = 1e-1)

        ctl.u.h_sp = y_kin_trim.h_e - 100
        step!(sim, 1, true)
        @test ctl.y.lon_ctl_mode === lon_thr_EAS
        step!(sim, 80, true) #altitude is captured
        @test ctl.y.lon_ctl_mode === lon_EAS_clm
        @test isapprox.(y_kin(ac).h_e - ctl.u.h_sp, 0.0; atol = 1e-1)

        @test ctl.y.lon_ctl_mode === lon_EAS_clm

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        ac.n[] = 0
        @test @ballocated(f_disc!($ac)) == 0

        ctl.u.h_sp = y_kin_trim.h_e + 100
        step!(sim, 1, true)
        @test ctl.y.lon_ctl_mode === lon_thr_EAS
        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        ac.n[] = 0
        @test @ballocated(f_disc!($ac)) == 0

        # kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
        # air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
        # save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172_fbw_v1", "avionics", "kin"))
        # save_plots(air_plots, save_folder = joinpath("tmp", "test_c172_fbw_v1", "avionics", "air"))
        # return TimeSeries(sim)

    end

    end #testset

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
    ac = Cessna172FBWv1(WA(), trn) |> System;

    sim = Simulation(ac; dt = 1/60, Δt = 1/60, t_end = 1000)
    reinit!(sim, initializer)

    for joystick in get_connected_joysticks()
        Sim.attach!(sim, joystick)
    end

    xpc = XPlaneOutput()
    # xpc = XPlaneOutput(address = IPv4("192.168.1.2"))
    Sim.attach!(sim, xpc)

    Sim.run_interactive!(sim)

    kin_plots = make_plots(TimeSeries(sim).vehicle.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).vehicle.air; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172fbw_v1", "sim_interactive", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172fbw_v1", "sim_interactive", "air"))

    return nothing

end


end #module