module TestC172Xv1

using Test, UnPack, BenchmarkTools, Sockets, JSON3

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

#non-exported stuff
using Flight.FlightLib.Control.Discrete: load_pid_lookup, load_lqr_tracker_lookup
using Flight.FlightAircraft.C172X.C172XControl: ModeControlLon, ModeControlLat
using Flight.FlightAircraft.C172X.C172XControl: ModeGuidanceLon, ModeGuidanceLat
using Flight.FlightAircraft.C172X.C172XControl: FlightPhase

export test_c172x1


function test_c172x1()
    @testset verbose = true "Cessna 172Xv1" begin

        test_control_modes()
        test_guidance_modes()

    end
end

y_kin(aircraft::Model{<:Cessna172Xv1}) = aircraft.y.vehicle.kinematics
y_air(aircraft::Model{<:Cessna172Xv1}) = aircraft.y.vehicle.airflow
y_aero(aircraft::Model{<:Cessna172Xv1}) = aircraft.y.vehicle.systems.aero

function test_control_modes()

    data_folder = joinpath(dirname(dirname(dirname(@__DIR__))),
        normpath("src/aircraft/c172/c172x/control/data"))

    @testset verbose = true "Control Modes" begin

    h_trn = HOrth(0.0)
    world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

    aircraft = world.aircraft
    act = aircraft.vehicle.systems.act
    ctl = aircraft.avionics.ctl

    init_gnd = C172.Init(KinInit( h = h_trn + 1.9))
    init_air = C172.TrimParameters()

    dt = Δt = 0.01

    sim = Simulation(world; dt, Δt, t_end = 600)

    ############################################################################
    ############################## Ground ######################################

    @testset verbose = true "Ground" begin

    Sim.init!(sim, init_gnd)

    #set arbitrary control and guidance modes
    ctl.u.mode_gdc_lon_req = ModeGuidanceLon.alt
    ctl.u.mode_gdc_lat_req = ModeGuidanceLat.seg
    ctl.u.mode_ctl_lon_req = ModeControlLon.EAS_clm
    ctl.u.mode_ctl_lat_req = ModeControlLat.p_β
    ctl.u.throttle_axis = 0.1
    ctl.u.aileron_axis = 0.2
    ctl.u.elevator_axis = 0.3
    ctl.u.rudder_axis = 0.4

    #step for one controller sample period
    step!(sim, ctl.Δt, true)

    #make sure we're on the ground
    @test ctl.y.flight_phase === FlightPhase.gnd

    #the mode requests are overridden due to FlightPhase.gnd
    @test ctl.y.mode_gdc_lon === ModeGuidanceLon.off
    @test ctl.y.mode_gdc_lat === ModeGuidanceLat.off
    @test ctl.y.mode_ctl_lon === ModeControlLon.direct
    @test ctl.y.mode_ctl_lat === ModeControlLat.direct

    #control laws outputs must have propagated to actuator inputs (not yet to
    #their outputs, that requires a subsequent call to f_ode!)
    @test act.throttle.u[] == 0.1
    @test act.aileron.u[] == 0.2
    @test act.elevator.u[] == 0.3
    @test act.rudder.u[] == 0.4

    #test for allocation. f_periodic! will need further tests in other control modes
    @test @ballocated(f_ode!($world)) == 0
    @test @ballocated(f_step!($world)) == 0
    @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end #testset


    # ############################################################################
    # ################################# Air ######################################

    @testset verbose = true "Air" begin

    #put the aircraft in its nominal design point
    Sim.init!(sim, init_air)
    y_kin_trim = y_kin(aircraft)

    ############################### direct control #############################

    @testset verbose = true "ModeControlLon.direct + ModeControlLat.direct" begin

        Sim.init!(sim, init_air)
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.direct
        @test ctl.y.mode_ctl_lat === ModeControlLat.direct

        #with direct surface control, trim state must be initially preserved
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b, y_kin_trim.ω_wb_b; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b, y_kin_trim.v_eb_b; atol = 1e-2))

        #test for allocations in the controller's current state
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end #testset

    ############################ thr+ele SAS mode ##############################

    @testset verbose = true "ModeControlLon.sas" begin

        #we test the longitudinal SAS first, because we want to test the lateral
        #modes with it enabled
        Sim.init!(sim, init_air)
        ctl.u.mode_ctl_lon_req = ModeControlLon.sas
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.sas

        #check the correct parameters are loaded and assigned to the controller
        te2te_lookup = load_lqr_tracker_lookup(joinpath(data_folder, "e2e_lookup.h5"))
        C_fwd = te2te_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).C_fwd
        @test all(isapprox.(ctl.y.ctl_lon.e2e_lqr.C_fwd, C_fwd; atol = 1e-6))

        #a small initial transient when engaging the SAS is acceptable
        #once active, trim equilibrium must be preserved for longer
        step!(sim, 30, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end #testset

    ################################ ail_rud SAS ###############################

    @testset verbose = true "ModeControlLat.sas" begin

        Sim.init!(sim, init_air)
        ctl.u.mode_ctl_lon_req = ModeControlLon.sas
        ctl.u.mode_ctl_lat_req = ModeControlLat.sas
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lat === ModeControlLat.sas

        #check the correct parameters are loaded and assigned to the controller
        ar2ar_lookup = load_lqr_tracker_lookup(joinpath(data_folder, "ar2ar_lookup.h5"))
        C_fwd = ar2ar_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).C_fwd
        @test all(isapprox.(ctl.y.ctl_lat.ar2ar_lqr.C_fwd, C_fwd; atol = 1e-6))

        #with ail+rud SAS active, trim state must be preserved for longer
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[1], y_kin_trim.ω_wb_b[1]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end


    ################################ φ + β #####################################

    @testset verbose = true "ModeControlLat.φ_β" begin

        Sim.init!(sim, init_air)
        ctl.u.mode_ctl_lon_req = ModeControlLon.sas
        ctl.u.mode_ctl_lat_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lat === ModeControlLat.φ_β

        #check the correct parameters are loaded and assigned to the controller
        φβ2ar_lookup = load_lqr_tracker_lookup(joinpath(data_folder, "φβ2ar_lookup.h5"))
        C_fwd = φβ2ar_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).C_fwd
        @test all(isapprox.(ctl.y.ctl_lat.φβ2ar_lqr.C_fwd, C_fwd; atol = 1e-6))

        #a small initial transient when engaging the SAS is acceptable
        #once active, trim equilibrium must be preserved
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[1]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/12
        ctl.u.β_ref = deg2rad(3)
        step!(sim, 10, true)
        @test isapprox(ctl.u.φ_ref, y_kin(aircraft).e_nb.φ; atol = 1e-3)
        @test isapprox(Float64(ctl.u.β_ref), y_aero(aircraft).β; atol = 1e-3)

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ################################ p + β #####################################

    @testset verbose = true "ModeControlLat.p_β" begin

        Sim.init!(sim, init_air)

        #enable SAS first and let the small initial transient die out
        ctl.u.mode_ctl_lon_req = ModeControlLon.sas
        ctl.u.mode_ctl_lat_req = ModeControlLat.sas
        step!(sim, 1, true)

        ctl.u.mode_ctl_lat_req = ModeControlLat.p_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lat === ModeControlLat.p_β

        #check the correct parameters are loaded and assigned to the controller
        p2φ_lookup = load_pid_lookup(joinpath(data_folder, "p2φ_lookup.h5"))
        k_p = p2φ_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.ctl_lat.p2φ_pid.k_p, k_p; atol = 1e-6))

        #the control mode must activate without transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #the controller must keep trim values in steady state
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        ctl.u.p_ref = 0.02
        ctl.u.β_ref = deg2rad(3)
        step!(sim, 10, true)
        @test isapprox(Float64(ctl.u.p_ref), y_kin(aircraft).ω_wb_b[1]; atol = 1e-3)
        @test isapprox(ctl.u.β_ref, y_aero(aircraft).β; atol = 1e-3)

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end


    ################################ χ + β #####################################

    @testset verbose = true "ModeControlLat.χ_β" begin

        Sim.init!(sim, init_air)

        #enable SAS first and let the small initial transient die out
        ctl.u.mode_ctl_lon_req = ModeControlLon.sas
        ctl.u.mode_ctl_lat_req = ModeControlLat.sas
        step!(sim, 1, true)

        ctl.u.mode_ctl_lat_req = ModeControlLat.χ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lat === ModeControlLat.χ_β

        #check the correct parameters are loaded and assigned to the controller
        χ2φ_lookup = load_pid_lookup(joinpath(data_folder, "χ2φ_lookup.h5"))
        k_p = χ2φ_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.ctl_lat.χ2φ_pid.k_p, k_p; atol = 1e-6))

        #with reference values matching their trim values, the control mode must activate
        #without transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking
        ctl.u.χ_ref = π/2
        step!(sim, 29, true)
        @test ctl.ctl_lat.u.χ_ref != 0
        @test isapprox(ctl.u.χ_ref, y_kin(aircraft).χ_gnd; atol = 1e-2)
        # @test isapprox(Float64(ctl.u.yaw_axis), y_aero(aircraft).β; atol = 1e-3)

        #correct tracking with 10m/s of crosswind (N, current heading is E)
        world.atmosphere.wind.u.N = 10
        step!(sim, 10, true)
        @test isapprox(ctl.u.χ_ref, y_kin(aircraft).χ_gnd; atol = 1e-2)
        world.atmosphere.wind.u.N = 0

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ############################################################################

    #now we proceed to test the remaining longitudinal modes with lateral p + β
    #mode enabled

    ############################### lon_thr_q ##################################

    @testset verbose = true "ModeControlLon.thr_q" begin

        Sim.init!(sim, init_air)

        ctl.u.mode_ctl_lon_req = ModeControlLon.thr_q
        ctl.u.mode_ctl_lat_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.thr_q

        #check the correct parameters are loaded and assigned to the controller
        q2e_lookup = load_pid_lookup(joinpath(data_folder, "q2e_lookup.h5"))
        k_p = q2e_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.ctl_lon.q2e_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/12
        ctl.u.q_ref = 0.01
        step!(sim, 10, true)

        @test ctl.ctl_lon.u.q_ref != 0
        @test isapprox(ctl.ctl_lon.u.q_ref, y_kin(aircraft).ω_wb_b[2]; atol = 1e-3)
        @test isapprox(Float64(aircraft.y.vehicle.systems.act.throttle.cmd),
                        Float64(ctl.u.throttle_axis + ctl.u.throttle_offset); atol = 1e-3)

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ############################## lon_thr_θ ###################################

    @testset verbose = true "ModeControlLon.thr_θ" begin

        Sim.init!(sim, init_air)

        ctl.u.mode_ctl_lon_req = ModeControlLon.thr_θ
        ctl.u.mode_ctl_lat_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.thr_θ

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/6
        ctl.u.θ_ref = deg2rad(5)
        step!(sim, 10, true)
        @test isapprox(y_kin(aircraft).e_nb.θ, ctl.u.θ_ref; atol = 1e-4)

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end


    ################################ lon_thr_EAS ###############################

    @testset verbose = true "ModeControlLon.thr_EAS" begin

        Sim.init!(sim, init_air)

        ctl.u.mode_ctl_lon_req = ModeControlLon.thr_EAS
        ctl.u.mode_ctl_lat_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.thr_EAS

        #check the correct parameters are loaded and assigned to the controller
        v2θ_lookup = load_pid_lookup(joinpath(data_folder, "v2θ_lookup.h5"))
        k_p = v2θ_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.ctl_lon.v2θ_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/6
        ctl.u.EAS_ref = 45
        step!(sim, 30, true)
        @test all(isapprox.(y_air(aircraft).EAS, ctl.u.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ################################ lon_EAS_q #################################

    @testset verbose = true "ModeControlLon.EAS_q" begin

        Sim.init!(sim, init_air)

        ctl.u.mode_ctl_lon_req = ModeControlLon.EAS_q
        ctl.u.mode_ctl_lat_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_q

        #check the correct parameters are loaded and assigned to v2t, the q
        #tracker is shared with other modes
        v2t_lookup = load_pid_lookup(joinpath(data_folder, "v2t_lookup.h5"))
        k_p = v2t_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.ctl_lon.v2t_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking
        ctl.u.q_ref = -0.01
        step!(sim, 10, true)
        ctl.u.q_ref = 0.005
        step!(sim, 10, true)
        ctl.u.q_ref = 0.0
        step!(sim, 20, true)

        @test isapprox(ctl.ctl_lon.u.q_ref, y_kin(aircraft).ω_wb_b[2]; atol = 1e-3)
        @test all(isapprox.(y_air(aircraft).EAS, ctl.u.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end


    ################################ lon_EAS_q #################################

    @testset verbose = true "ModeControlLon.EAS_θ" begin

        Sim.init!(sim, init_air)

        ctl.u.mode_ctl_lon_req = ModeControlLon.EAS_θ
        ctl.u.mode_ctl_lat_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_θ

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 0.1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/6
        ctl.u.θ_ref = deg2rad(3)
        step!(sim, 10, true)
        ctl.u.θ_ref = -deg2rad(3)
        step!(sim, 60, true)

        @test isapprox(ctl.ctl_lon.u.θ_ref, y_kin(aircraft).e_nb.θ; atol = 1e-3)
        @test all(isapprox.(y_air(aircraft).EAS, ctl.u.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ############################## lon_EAS_clm #################################

    @testset verbose = true "ModeControlLon.EAS_clm" begin

        Sim.init!(sim, init_air)

        ctl.u.mode_ctl_lon_req = ModeControlLon.EAS_clm
        ctl.u.mode_ctl_lat_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm

        #check the correct parameters are loaded and assigned to the controller
        c2θ_lookup = load_pid_lookup(joinpath(data_folder, "c2θ_lookup.h5"))
        k_p = c2θ_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.ctl_lon.c2θ_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/6
        ctl.u.EAS_ref = 45
        ctl.u.clm_ref = 2
        step!(sim, 30, true)
        @test all(isapprox.(y_kin(aircraft).v_eb_n[3], -ctl.u.clm_ref; atol = 1e-1))
        @test all(isapprox.(y_air(aircraft).EAS, ctl.u.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

        # return sim

    end #testset

    end #testset

    end #testset

end #function


function test_guidance_modes()

    @testset verbose = true "Guidance Modes" begin

    h_trn = HOrth(0.0)
    world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

    aircraft = world.aircraft
    ctl = aircraft.avionics.ctl

    init_air = C172.TrimParameters()
    dt = Δt = 0.01

    sim = Simulation(world; dt, Δt, t_end = 600)

    @testset verbose = true "ModeGuidanceLon.alt" begin

        Sim.init!(sim, init_air)
        y_kin_trim = y_kin(aircraft)

        ctl.u.mode_gdc_lon_req = ModeGuidanceLon.alt
        ctl.u.mode_ctl_lat_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_gdc_lon === ModeGuidanceLon.alt
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm

        #when trim reference values are kept, the guidance mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #all tests while turning
        ctl.u.φ_ref = π/12

        ctl.u.h_target = y_kin_trim.h_e + 100
        step!(sim, 1, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.thr_EAS
        step!(sim, 60, true) #altitude is captured
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm
        @test isapprox.(y_kin(aircraft).h_e - HEllip(ctl.u.h_target), 0.0; atol = 1e-1)

        #reference changes within the current threshold do not prompt a mode change
        ctl.u.h_target = y_kin(aircraft).h_e - ctl.gdc_lon_alt.constants.h_thr / 2
        step!(sim, 1, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm
        step!(sim, 30, true) #altitude is captured
        @test isapprox.(y_kin(aircraft).h_e - HEllip(ctl.u.h_target), 0.0; atol = 1e-1)

        ctl.u.h_target = y_kin_trim.h_e - 100
        step!(sim, 1, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.thr_EAS
        step!(sim, 80, true) #altitude is captured
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm
        @test isapprox.(y_kin(aircraft).h_e - HEllip(ctl.u.h_target), 0.0; atol = 1e-1)

        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm

        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

        ctl.u.h_target = y_kin_trim.h_e + 100
        step!(sim, 1, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.thr_EAS

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

        # return TimeSeries(sim)

    end

    end #testset

end


############################### JSON Loopback Test #############################

struct JSONTestMapping <: IOMapping end

function IODevices.extract_output(mdl::Model{<:SimpleWorld}, ::JSONTestMapping)
    freq = 0.1
    φ_ref_max = π/6
    φ_ref = φ_ref_max * sin(2π*freq*mdl.t[])

    #these are all valid empty JSON entities. when passed to JSON3.write, they
    #yield respectively "\"\"", "[]" and "{}", all of length 2
    cmd = ""
    # cmd = []
    # cmd = Dict()

    if mdl.t[] > 5
        #these enums will be automatically cast to Ints per the StructTypes
        #methods defined in C172X.C172XControl
        cmd = (
            mode_gdc_lon_req = ModeGuidanceLon.alt,
            mode_ctl_lat_req = ModeControlLat.φ_β,
            φ_ref = φ_ref,
        )

        #therefore, these would also work
        # cmd = (
        #     mode_gdc_lon_req = 1,
        #     mode_ctl_lat_req = 2,
        #     φ_ref = φ_ref,
        # )
    end

    return JSON3.write(cmd)
end

function IODevices.assign_input!(world::Model{<:SimpleWorld}, ::JSONTestMapping, data::String)
    #ControllerU is declared as StructTypes.Mutable() in C172X.C172XControl,
    #so JSON3 can automatically read a JSON string into one or more of its
    #fields

    #caution: String(data) empties the original data::Vector{UInt8}, so
    #further calls would return an empty string
    str = String(data)

    #if length(str) < 2, it cannot be a valid JSON string. if length(str) == 2
    #it is an empty JSON entity (either string, object or array). instead of
    #this check we could simply call isempty(JSON3.read(str)) but that would
    #mean parsing the string twice
    length(str) > 2 && JSON3.read!(str, world.aircraft.avionics.ctl.u)

    # isempty(str) |> println
    # JSON3.read(str) |> isempty |> println
end

function test_json_loopback(; save::Bool = true)


    h_trn = HOrth(427.2);
    world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

    sim = Simulation(world; dt = 1/60, Δt = 1/60, t_end = 30)

    #on air, automatically trimmed by reinit!
    initializer = C172.TrimParameters(
        Ob = Geographic(LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)), HEllip(650)))

    #initialize simulated system
    Sim.init!(sim, initializer)

    #the loopback interface must not share its port with the XPlane12Control!
    Sim.attach!(sim, XPlane12Control(; port = 49000))
    Sim.attach!(sim, UDPInput(; port = 49017), JSONTestMapping())
    Sim.attach!(sim, UDPOutput(; port = 49017), JSONTestMapping())

    #trigger compilation of parsing methods for AvionicsU before launching the
    #simulation
    JSON3.read!(JSON3.write(world.aircraft.avionics.ctl.u, allow_inf=true), world.aircraft.avionics.ctl.u; allow_inf=true)

    Sim.run_interactive!(sim)

    save && save_plots(TimeSeries(sim).aircraft.vehicle.kinematics,
                        normpath("tmp/plots/test_c172x1/test_json_loopback/kin");
                        Plotting.defaults...)

    return nothing

end

end #module