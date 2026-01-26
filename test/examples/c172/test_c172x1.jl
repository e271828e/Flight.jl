module TestC172Xv1

using Test, BenchmarkTools, Sockets, JSON3, Logging

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightExamples

using Flight.FlightExamples.C172: is_on_gnd
using Flight.FlightLib.Control.Discrete: build_lookup_pid, build_lookup_lqr
using Flight.FlightExamples.C172X.C172XControl: ModeControlLon, ModeControlLat, AltTrackingState

export test_c172x1

y_kin(aircraft::Model{<:Cessna172Xv1}) = aircraft.y.vehicle.kinematics
y_air(aircraft::Model{<:Cessna172Xv1}) = aircraft.y.vehicle.airflow
y_aero(aircraft::Model{<:Cessna172Xv1}) = aircraft.y.vehicle.systems.aero


function test_c172x1(; alloc::Bool = true)

    @testset verbose = true "Cessna 172Xv1" begin

    data_folder = joinpath(dirname(dirname(dirname(@__DIR__))),
        normpath("src/examples/c172/c172x/control/data"))

    h_trn = HOrth(0.0)
    world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

    aircraft = world.aircraft
    act = aircraft.vehicle.systems.act
    ctl = aircraft.avionics

    init_gnd = C172.Init(KinInit( h = h_trn + 1.9))

    #using the default TrimParameters() is crucial for these tests. this ensures
    #we are exactly at one of the design points, with exactly computed
    #controller parameters, rather than ones interpolated from the lookup
    #tables. this of course would not be a problem performance-wise, but it is
    #when assessing whether the SAS loops exactly respect the trim condition
    init_air = C172.TrimParameters()

    dt = Δt = 0.01

    sim = Simulation(world; dt, Δt, t_end = 600)

    #verify parent-child input linkage
    @test ctl.u.lon === ctl.lon.u
    @test ctl.u.lat === ctl.lat.u

    ############################################################################
    ############################## Ground ######################################

    @testset verbose = true "Ground" begin

    init!(sim, init_gnd)

    #set arbitrary control modes
    ctl.u.lon.mode_req = ModeControlLon.EAS_clm
    ctl.u.lat.mode_req = ModeControlLat.p_β
    ctl.u.lon.throttle_axis = 0.1
    ctl.u.lon.elevator_axis = 0.3
    ctl.u.lat.aileron_axis = 0.2
    ctl.u.lat.rudder_axis = 0.4

    #step for one controller sample period
    step!(sim, ctl.Δt, true)

    #make sure we're on the ground
    @test is_on_gnd(aircraft.vehicle)

    #the mode requests are overridden due to wow = gnd
    @test ctl.y.lon.mode === ModeControlLon.direct
    @test ctl.y.lat.mode === ModeControlLat.direct

    #control laws outputs must have propagated to actuator inputs (not yet to
    #their outputs, that requires a subsequent call to f_ode!)
    @test act.throttle.u[] == 0.1
    @test act.elevator.u[] == 0.3
    @test act.rudder.u[] == 0.4
    @test act.aileron.u[] == 0.2

    if alloc
        @test @ballocated(f_ode!($world)) == 0
        @test @ballocated(f_step!($world)) == 0
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0
    end

    end #testset

    # ############################################################################
    # ################################# Air ######################################

    # @testset verbose = true "Air" begin

    #put the aircraft at its nominal design point and save its output for later
    init!(sim, init_air)
    y_kin_trim = y_kin(aircraft)

    ############################### direct control #############################

    @testset verbose = true "ModeControlLon.direct + ModeControlLat.direct" begin

        init!(sim, init_air)
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon.mode === ModeControlLon.direct
        @test ctl.y.lat.mode === ModeControlLat.direct

        #with direct surface control, trim state must be initially preserved
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b, y_kin_trim.ω_wb_b; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b, y_kin_trim.v_eb_b; atol = 1e-2))

        #test for allocations in the controller's current state
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end #testset

    ############################ thr+ele SAS mode ##############################

    @testset verbose = true "ModeControlLon.sas" begin

        #we test the longitudinal SAS first, because we want to test the lateral
        #modes with it enabled
        init!(sim, init_air)
        ctl.u.lon.mode_req = ModeControlLon.sas
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon.mode === ModeControlLon.sas

        #check the correct parameters are loaded and assigned to the controller
        te2te_lookup = build_lookup_lqr(joinpath(data_folder, "te2te.h5"))
        K_fwd = te2te_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).K_fwd
        @test all(isapprox.(ctl.y.lon.te2te_lqr.K_fwd, K_fwd; atol = 1e-6))

        #a small initial transient when engaging the SAS is acceptable
        #once active, trim equilibrium must be preserved for longer
        step!(sim, 30, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end #testset

    ################################ ail_rud SAS ###############################

    @testset verbose = true "ModeControlLat.sas" begin

        init!(sim, init_air)
        ctl.u.lon.mode_req = ModeControlLon.sas
        ctl.u.lat.mode_req = ModeControlLat.sas
        step!(sim, ctl.Δt, true)
        @test ctl.y.lat.mode === ModeControlLat.sas

        #check the correct parameters are loaded and assigned to the controller
        ar2ar_lookup = build_lookup_lqr(joinpath(data_folder, "ar2ar.h5"))
        K_fwd = ar2ar_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).K_fwd
        @test all(isapprox.(ctl.y.lat.ar2ar_lqr.K_fwd, K_fwd; atol = 1e-6))

        #with ail+rud SAS active, trim state must be preserved for longer
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[1], y_kin_trim.ω_wb_b[1]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ################################ φ + β #####################################

    @testset verbose = true "ModeControlLat.φ_β" begin

        init!(sim, init_air)
        ctl.u.lon.mode_req = ModeControlLon.sas
        ctl.u.lat.mode_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lat.mode === ModeControlLat.φ_β

        #check the correct parameters are loaded and assigned to the controller
        φβ2ar_lookup = build_lookup_lqr(joinpath(data_folder, "φβ2ar.h5"))
        K_fwd = φβ2ar_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).K_fwd
        @test all(isapprox.(ctl.y.lat.φβ2ar_lqr.K_fwd, K_fwd; atol = 1e-6))

        #a small initial transient when engaging the SAS is acceptable
        #once active, trim equilibrium must be preserved
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[1]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.lat.φ_ref = π/12
        ctl.u.lat.β_ref = deg2rad(3)
        step!(sim, 10, true)
        @test isapprox(ctl.u.lat.φ_ref, y_kin(aircraft).e_nb.φ; atol = 1e-3)
        @test isapprox(Float64(ctl.u.lat.β_ref), y_aero(aircraft).β; atol = 1e-3)

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ################################ p + β #####################################

    @testset verbose = true "ModeControlLat.p_β" begin

        init!(sim, init_air)

        #enable SAS first and let the small initial transient die out
        ctl.u.lon.mode_req = ModeControlLon.sas
        ctl.u.lat.mode_req = ModeControlLat.sas
        step!(sim, 1, true)

        ctl.u.lat.mode_req = ModeControlLat.p_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lat.mode === ModeControlLat.p_β

        #check the correct parameters are loaded and assigned to the controller
        p2φ_lookup = build_lookup_pid(joinpath(data_folder, "p2φ.h5"))
        k_p = p2φ_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.lat.p2φ_pid.k_p, k_p; atol = 1e-6))

        #the control mode must activate without transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #the controller must keep trim values in steady state
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        ctl.u.lat.p_ref = 0.02
        ctl.u.lat.β_ref = deg2rad(3)
        step!(sim, 10, true)
        @test isapprox(Float64(ctl.u.lat.p_ref), y_kin(aircraft).ω_wb_b[1]; atol = 1e-3)
        @test isapprox(ctl.u.lat.β_ref, y_aero(aircraft).β; atol = 1e-3)

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ################################ χ + β #####################################

    @testset verbose = true "ModeControlLat.χ_β" begin

        init!(sim, init_air)

        #enable SAS first and let the small initial transient die out
        ctl.u.lon.mode_req = ModeControlLon.sas
        ctl.u.lat.mode_req = ModeControlLat.sas
        step!(sim, 1, true)

        ctl.u.lat.mode_req = ModeControlLat.χ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lat.mode === ModeControlLat.χ_β

        #check the correct parameters are loaded and assigned to the controller
        χ2φ_lookup = build_lookup_pid(joinpath(data_folder, "χ2φ.h5"))
        k_p = χ2φ_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.lat.χ2φ_pid.k_p, k_p; atol = 1e-6))

        #with reference values matching their trim values, the control mode must activate
        #without transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking
        ctl.u.lat.χ_ref = π/2
        step!(sim, 29, true)
        @test ctl.lat.u.χ_ref != 0
        @test isapprox(ctl.u.lat.χ_ref, y_kin(aircraft).χ_gnd; atol = 1e-2)
        # @test isapprox(Float64(ctl.u.yaw_axis), y_aero(aircraft).β; atol = 1e-3)

        #correct tracking with 10m/s of crosswind (N, current heading is E)
        world.atmosphere.wind.u.N = 10
        step!(sim, 10, true)
        @test isapprox(ctl.u.lat.χ_ref, y_kin(aircraft).χ_gnd; atol = 1e-2)
        world.atmosphere.wind.u.N = 0

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ############################################################################

    #test the remaining longitudinal modes with lateral p + β mode enabled

    ############################### lon_thr_q ##################################

    @testset verbose = true "ModeControlLon.thr_q" begin

        init!(sim, init_air)

        ctl.u.lon.mode_req = ModeControlLon.thr_q
        ctl.u.lat.mode_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon.mode === ModeControlLon.thr_q

        #check the correct parameters are loaded and assigned to the controller
        q2e_lookup = build_lookup_pid(joinpath(data_folder, "q2e.h5"))
        k_p = q2e_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.lon.q2e_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.lat.φ_ref = π/12
        ctl.u.lon.q_ref = 0.01
        step!(sim, 10, true)

        @test ctl.u.lon.q_ref != 0
        @test isapprox(ctl.u.lon.q_ref, y_kin(aircraft).ω_wb_b[2]; atol = 1e-3)
        @test isapprox(Float64(aircraft.y.vehicle.systems.act.throttle.cmd),
                        Float64(ctl.u.lon.throttle_axis + ctl.u.lon.throttle_offset); atol = 1e-3)

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ############################## lon_thr_θ ###################################

    @testset verbose = true "ModeControlLon.thr_θ" begin

        init!(sim, init_air)

        ctl.u.lon.mode_req = ModeControlLon.thr_θ
        ctl.u.lat.mode_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon.mode === ModeControlLon.thr_θ

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.lat.φ_ref = π/6
        ctl.u.lon.θ_ref = deg2rad(5)
        step!(sim, 10, true)
        @test isapprox(y_kin(aircraft).e_nb.θ, ctl.u.lon.θ_ref; atol = 1e-4)

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ################################ lon_thr_EAS ###############################

    @testset verbose = true "ModeControlLon.thr_EAS" begin

        init!(sim, init_air)

        ctl.u.lon.mode_req = ModeControlLon.thr_EAS
        ctl.u.lat.mode_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon.mode === ModeControlLon.thr_EAS

        #check the correct parameters are loaded and assigned to the controller
        tv2te_lookup = build_lookup_lqr(joinpath(data_folder, "tv2te.h5"))
        K_fwd = tv2te_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).K_fwd
        @test all(isapprox.(ctl.y.lon.tv2te_lqr.K_fwd, K_fwd; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.lat.φ_ref = π/6
        ctl.u.lon.EAS_ref = 45
        step!(sim, 30, true)
        @test all(isapprox.(y_air(aircraft).EAS, ctl.u.lon.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ################################ lon_EAS_q #################################

    @testset verbose = true "ModeControlLon.EAS_q" begin

        init!(sim, init_air)

        ctl.u.lon.mode_req = ModeControlLon.EAS_q
        ctl.u.lat.mode_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon.mode === ModeControlLon.EAS_q

        #check the correct parameters are loaded and assigned to v2t, the q
        #tracker is shared with other modes
        v2t_lookup = build_lookup_pid(joinpath(data_folder, "v2t.h5"))
        k_p = v2t_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.lon.v2t_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking
        ctl.u.lon.q_ref = -0.005
        step!(sim, 20, true)
        @test isapprox(ctl.u.lon.q_ref, y_kin(aircraft).ω_wb_b[2]; atol = 1e-3)
        @test all(isapprox.(y_air(aircraft).EAS, ctl.u.lon.EAS_ref; atol = 1)) #allow 1m/s error

        ctl.u.lon.q_ref = 0.005
        step!(sim, 20, true)
        @test isapprox(ctl.u.lon.q_ref, y_kin(aircraft).ω_wb_b[2]; atol = 1e-3)
        @test all(isapprox.(y_air(aircraft).EAS, ctl.u.lon.EAS_ref; atol = 1)) #allow 1m/s error

        ctl.u.lon.q_ref = 0.0
        step!(sim, 20, true)
        @test isapprox(ctl.u.lon.q_ref, y_kin(aircraft).ω_wb_b[2]; atol = 1e-3)
        @test all(isapprox.(y_air(aircraft).EAS, ctl.u.lon.EAS_ref; atol = 1)) #allow 1m/s error

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ################################ lon_EAS_q #################################

    @testset verbose = true "ModeControlLon.EAS_θ" begin

        init!(sim, init_air)

        ctl.u.lon.mode_req = ModeControlLon.EAS_θ
        ctl.u.lat.mode_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon.mode === ModeControlLon.EAS_θ

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 0.1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.lat.φ_ref = π/6
        ctl.u.lon.θ_ref = deg2rad(3)
        step!(sim, 10, true)
        ctl.u.lon.θ_ref = -deg2rad(3)
        step!(sim, 60, true)

        @test isapprox(ctl.u.lon.θ_ref, y_kin(aircraft).e_nb.θ; atol = 1e-3)
        @test all(isapprox.(y_air(aircraft).EAS, ctl.u.lon.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    ############################## lon_EAS_clm #################################

    @testset verbose = true "ModeControlLon.EAS_clm" begin

        init!(sim, init_air)

        ctl.u.lon.mode_req = ModeControlLon.EAS_clm
        ctl.u.lat.mode_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon.mode === ModeControlLon.EAS_clm

        #check the correct parameters are loaded and assigned to the controller
        c2θ_lookup = build_lookup_pid(joinpath(data_folder, "c2θ.h5"))
        k_p = c2θ_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).k_p
        @test all(isapprox.(ctl.y.lon.c2θ_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.lat.φ_ref = π/6
        ctl.u.lon.EAS_ref = 45
        ctl.u.lon.clm_ref = 2
        step!(sim, 30, true)
        @test all(isapprox.(y_kin(aircraft).v_eb_n[3], -ctl.u.lon.clm_ref; atol = 1e-1))
        @test all(isapprox.(y_air(aircraft).EAS, ctl.u.lon.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

        # return sim

    end #testset

    @testset verbose = true "ModeControlLon.alt" begin

        init!(sim, init_air)

        ctl.u.lon.mode_req = ModeControlLon.EAS_alt
        ctl.u.lat.mode_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)

        #check the correct parameters are loaded and assigned to the controller
        vh2te_lookup = build_lookup_lqr(joinpath(data_folder, "vh2te.h5"))
        K_fwd = vh2te_lookup(y_air(aircraft).EAS, Float64(y_kin(aircraft).h_e)).K_fwd
        @test all(isapprox.(ctl.y.lon.vh2te_lqr.K_fwd, K_fwd; atol = 1e-6))

        #h_ref should have been initialized at its trim value, so the initial
        #altitude tracking state should be hold
        @test ctl.y.lon.h_state === AltTrackingState.hold
        @test ctl.y.lon.mode === ModeControlLon.EAS_alt

        #when trim reference values are kept, the guidance mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #all tests while turning
        ctl.u.lat.φ_ref = π/12

        ctl.u.lon.h_ref = y_kin_trim.h_e + 100
        step!(sim, 1, true)
        @test ctl.y.lon.h_state === AltTrackingState.acquire
        @test ctl.y.lon.mode === ModeControlLon.thr_EAS

        step!(sim, 60, true) #altitude is captured
        @test ctl.y.lon.h_state === AltTrackingState.hold
        @test isapprox.(y_kin(aircraft).h_e - ctl.u.lon.h_ref, 0.0; atol = 1e-1)

        #reference changes within the current threshold do not prompt a mode change
        ctl.u.lon.h_ref = y_kin(aircraft).h_e - ctl.lon.parameters.h_thr / 2
        step!(sim, 1, true)
        @test ctl.y.lon.h_state === AltTrackingState.hold
        step!(sim, 30, true) #altitude is captured
        @test isapprox.(y_kin(aircraft).h_e - ctl.u.lon.h_ref, 0.0; atol = 1e-1)

        ctl.u.lon.h_ref = y_kin_trim.h_e - 100
        step!(sim, 1, true)
        @test ctl.y.lon.h_state === AltTrackingState.acquire
        step!(sim, 80, true) #altitude is captured
        @test ctl.y.lon.h_state === AltTrackingState.hold
        @test isapprox.(y_kin(aircraft).h_e - ctl.u.lon.h_ref, 0.0; atol = 1e-1)

        #test for allocations in EAS_alt mode
        @test ctl.y.lon.mode === ModeControlLon.EAS_alt
        alloc && @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    end #testset

end #function



#for interactive ControlLaws validation
function test_interactive(init::C172.TrimParameters = C172.TrimParameters())

    world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain()) |> Model
    sim = Simulation(world; dt = 0.01, t_end = 1000)
    init!(sim, init)
    Sim.run!(sim; gui = true)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/test_c172x1/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/test_c172x1/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/test_c172x1/dyn"); Plotting.defaults...)
    # save_plots(TimeSeries(sim).aircraft.vehicle.dynamics; Plotting.defaults...)

    return sim

end


end #module