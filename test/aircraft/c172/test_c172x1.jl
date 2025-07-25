module TestC172Xv1

using Test, UnPack, BenchmarkTools, Sockets, JSON3

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

#non-exported stuff
using Flight.FlightLib.Control.Discrete: load_pid_lookup, load_lqr_tracker_lookup
using Flight.FlightAircraft.C172X.C172XControl: lon_direct, lon_sas, lon_thr_q, lon_thr_θ, lon_thr_EAS, lon_EAS_q, lon_EAS_θ, lon_EAS_clm
using Flight.FlightAircraft.C172X.C172XControl: lat_direct, lat_sas, lat_p_β, lat_φ_β, lat_χ_β
using Flight.FlightAircraft.C172X.C172XControl: vrt_gdc_off, vrt_gdc_alt
using Flight.FlightAircraft.C172X.C172XControl: hor_gdc_off, hor_gdc_line
using Flight.FlightAircraft.C172X.C172XControl: phase_gnd, phase_air

export test_c172x1


function test_c172x1()
    @testset verbose = true "Cessna 172Xv1" begin

        test_control_modes()
        test_guidance_modes()

    end
end

y_kin(ac::Model{<:Cessna172Xv1}) = ac.y.vehicle.kinematics
y_air(ac::Model{<:Cessna172Xv1}) = ac.y.vehicle.airflow
y_aero(ac::Model{<:Cessna172Xv1}) = ac.y.vehicle.systems.aero

function test_control_modes()

    data_folder = joinpath(dirname(dirname(dirname(@__DIR__))),
        normpath("src/aircraft/c172/c172x/control/data"))

    @testset verbose = true "Control Modes" begin

    h_trn = HOrth(0.0)
    world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

    ac = world.ac
    act = ac.vehicle.systems.act
    ctl = ac.avionics.ctl

    init_gnd = C172.Init(KinInit( h = h_trn + 1.9))
    init_air = C172.TrimParameters()

    dt = Δt = 0.01

    sim = Simulation(world; dt, Δt, t_end = 600)

    ############################################################################
    ############################## Ground ######################################

    @testset verbose = true "Ground" begin

    Sim.init!(sim, init_gnd)

    #set arbitrary control and guidance modes
    ctl.u.vrt_gdc_mode_req = vrt_gdc_alt
    ctl.u.hor_gdc_mode_req = hor_gdc_line
    ctl.u.lon_ctl_mode_req = lon_EAS_clm
    ctl.u.lat_ctl_mode_req = lat_p_β
    ctl.u.throttle_axis = 0.1
    ctl.u.aileron_axis = 0.2
    ctl.u.elevator_axis = 0.3
    ctl.u.rudder_axis = 0.4

    #step for one controller sample period
    step!(sim, ctl.Δt, true)

    #make sure we're on the ground
    @test ctl.y.flight_phase === phase_gnd

    #the mode requests are overridden due to phase_gnd
    @test ctl.y.vrt_gdc_mode === vrt_gdc_off
    @test ctl.y.hor_gdc_mode === hor_gdc_off
    @test ctl.y.lon_ctl_mode === lon_direct
    @test ctl.y.lat_ctl_mode === lat_direct

    #control laws outputs must have propagated to actuator inputs (not yet to
    #their outputs, that requires a subsequent call to f_ode!)
    @test act.throttle.u[] == 0.1
    @test act.aileron.u[] == 0.2
    @test act.elevator.u[] == 0.3
    @test act.rudder.u[] == 0.4

    #test for allocation. f_disc! will need further tests in other control modes
    @test @ballocated(f_ode!($world)) == 0
    @test @ballocated(f_step!($world)) == 0
    @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end #testset


    # ############################################################################
    # ################################# Air ######################################

    @testset verbose = true "Air" begin

    #put the aircraft in its nominal design point
    Sim.init!(sim, init_air)
    y_kin_trim = y_kin(ac)

    ############################### direct control #############################

    @testset verbose = true "lon_direct + lat_direct" begin

        Sim.init!(sim, init_air)
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_direct
        @test ctl.y.lat_ctl_mode === lat_direct

        #with direct surface control, trim state must be initially preserved
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b, y_kin_trim.ω_wb_b; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b, y_kin_trim.v_eb_b; atol = 1e-2))

        #test for allocations in the controller's current state
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end #testset

    ############################ thr+ele SAS mode ##############################

    @testset verbose = true "lon_sas" begin

        #we test the longitudinal SAS first, because we want to test the lateral
        #modes with it enabled
        Sim.init!(sim, init_air)
        ctl.u.lon_ctl_mode_req = lon_sas
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_sas

        #check the correct parameters are loaded and assigned to the controller
        te2te_lookup = load_lqr_tracker_lookup(joinpath(data_folder, "e2e_lookup.h5"))
        C_fwd = te2te_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
        @test all(isapprox.(ctl.y.lon_ctl.e2e_lqr.C_fwd, C_fwd; atol = 1e-6))

        #a small initial transient when engaging the SAS is acceptable
        #once active, trim equilibrium must be preserved for longer
        step!(sim, 30, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end #testset

    ################################ ail_rud SAS ###############################

    @testset verbose = true "lat_sas" begin

        Sim.init!(sim, init_air)
        ctl.u.lon_ctl_mode_req = lon_sas
        ctl.u.lat_ctl_mode_req = lat_sas
        step!(sim, ctl.Δt, true)
        @test ctl.y.lat_ctl_mode === lat_sas

        #check the correct parameters are loaded and assigned to the controller
        ar2ar_lookup = load_lqr_tracker_lookup(joinpath(data_folder, "ar2ar_lookup.h5"))
        C_fwd = ar2ar_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
        @test all(isapprox.(ctl.y.lat_ctl.ar2ar_lqr.C_fwd, C_fwd; atol = 1e-6))

        #with ail+rud SAS active, trim state must be preserved for longer
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[1], y_kin_trim.ω_wb_b[1]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end


    ################################ φ + β #####################################

    @testset verbose = true "lat_φ_β" begin

        Sim.init!(sim, init_air)
        ctl.u.lon_ctl_mode_req = lon_sas
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lat_ctl_mode === lat_φ_β

        #check the correct parameters are loaded and assigned to the controller
        φβ2ar_lookup = load_lqr_tracker_lookup(joinpath(data_folder, "φβ2ar_lookup.h5"))
        C_fwd = φβ2ar_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).C_fwd
        @test all(isapprox.(ctl.y.lat_ctl.φβ2ar_lqr.C_fwd, C_fwd; atol = 1e-6))

        #a small initial transient when engaging the SAS is acceptable
        #once active, trim equilibrium must be preserved
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[1]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/12
        ctl.u.β_ref = deg2rad(3)
        step!(sim, 10, true)
        @test isapprox(ctl.u.φ_ref, y_kin(ac).e_nb.φ; atol = 1e-3)
        @test isapprox(Float64(ctl.u.β_ref), y_aero(ac).β; atol = 1e-3)

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end

    ################################ p + β #####################################

    @testset verbose = true "lat_p_β" begin

        Sim.init!(sim, init_air)

        #enable SAS first and let the small initial transient die out
        ctl.u.lon_ctl_mode_req = lon_sas
        ctl.u.lat_ctl_mode_req = lat_sas
        step!(sim, 1, true)

        ctl.u.lat_ctl_mode_req = lat_p_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lat_ctl_mode === lat_p_β

        #check the correct parameters are loaded and assigned to the controller
        p2φ_lookup = load_pid_lookup(joinpath(data_folder, "p2φ_lookup.h5"))
        k_p = p2φ_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
        @test all(isapprox.(ctl.y.lat_ctl.p2φ_pid.k_p, k_p; atol = 1e-6))

        #the control mode must activate without transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #the controller must keep trim values in steady state
        step!(sim, 10, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        ctl.u.p_ref = 0.02
        ctl.u.β_ref = deg2rad(3)
        step!(sim, 10, true)
        @test isapprox(Float64(ctl.u.p_ref), y_kin(ac).ω_wb_b[1]; atol = 1e-3)
        @test isapprox(ctl.u.β_ref, y_aero(ac).β; atol = 1e-3)

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end


    ################################ χ + β #####################################

    @testset verbose = true "lat_χ_β" begin

        Sim.init!(sim, init_air)

        #enable SAS first and let the small initial transient die out
        ctl.u.lon_ctl_mode_req = lon_sas
        ctl.u.lat_ctl_mode_req = lat_sas
        step!(sim, 1, true)

        ctl.u.lat_ctl_mode_req = lat_χ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lat_ctl_mode === lat_χ_β

        #check the correct parameters are loaded and assigned to the controller
        χ2φ_lookup = load_pid_lookup(joinpath(data_folder, "χ2φ_lookup.h5"))
        k_p = χ2φ_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
        @test all(isapprox.(ctl.y.lat_ctl.χ2φ_pid.k_p, k_p; atol = 1e-6))

        #with reference values matching their trim values, the control mode must activate
        #without transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking
        ctl.u.χ_ref = π/2
        step!(sim, 29, true)
        @test ctl.lat_ctl.u.χ_ref != 0
        @test isapprox(ctl.u.χ_ref, y_kin(ac).χ_gnd; atol = 1e-2)
        # @test isapprox(Float64(ctl.u.yaw_axis), y_aero(ac).β; atol = 1e-3)

        #correct tracking with 10m/s of crosswind (N, current heading is E)
        world.atm.wind.u.N = 10
        step!(sim, 10, true)
        @test isapprox(ctl.u.χ_ref, y_kin(ac).χ_gnd; atol = 1e-2)
        world.atm.wind.u.N = 0

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end

    ############################################################################

    #now we proceed to test the remaining longitudinal modes with lateral p + β
    #mode enabled

    ############################### lon_thr_q ##################################

    @testset verbose = true "lon_thr_q" begin

        Sim.init!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_thr_q
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_thr_q

        #check the correct parameters are loaded and assigned to the controller
        q2e_lookup = load_pid_lookup(joinpath(data_folder, "q2e_lookup.h5"))
        k_p = q2e_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
        @test all(isapprox.(ctl.y.lon_ctl.q2e_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/12
        ctl.u.q_ref = 0.01
        step!(sim, 10, true)

        @test ctl.lon_ctl.u.q_ref != 0
        @test isapprox(ctl.lon_ctl.u.q_ref, y_kin(ac).ω_wb_b[2]; atol = 1e-3)
        @test isapprox(Float64(ac.y.vehicle.systems.act.throttle.cmd),
                        Float64(ctl.u.throttle_axis + ctl.u.throttle_offset); atol = 1e-3)

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end

    ############################## lon_thr_θ ###################################

    @testset verbose = true "lon_thr_θ" begin

        Sim.init!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_thr_θ
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_thr_θ

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/6
        ctl.u.θ_ref = deg2rad(5)
        step!(sim, 10, true)
        @test isapprox(y_kin(ac).e_nb.θ, ctl.u.θ_ref; atol = 1e-4)

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end


    ################################ lon_thr_EAS ###############################

    @testset verbose = true "lon_thr_EAS" begin

        Sim.init!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_thr_EAS
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_thr_EAS

        #check the correct parameters are loaded and assigned to the controller
        v2θ_lookup = load_pid_lookup(joinpath(data_folder, "v2θ_lookup.h5"))
        k_p = v2θ_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
        @test all(isapprox.(ctl.y.lon_ctl.v2θ_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/6
        ctl.u.EAS_ref = 45
        step!(sim, 30, true)
        @test all(isapprox.(y_air(ac).EAS, ctl.u.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end

    ################################ lon_EAS_q #################################

    @testset verbose = true "lon_EAS_q" begin

        Sim.init!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_EAS_q
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_EAS_q

        #check the correct parameters are loaded and assigned to v2t, the q
        #tracker is shared with other modes
        v2t_lookup = load_pid_lookup(joinpath(data_folder, "v2t_lookup.h5"))
        k_p = v2t_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
        @test all(isapprox.(ctl.y.lon_ctl.v2t_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking
        ctl.u.q_ref = -0.01
        step!(sim, 10, true)
        ctl.u.q_ref = 0.005
        step!(sim, 10, true)
        ctl.u.q_ref = 0.0
        step!(sim, 20, true)

        @test isapprox(ctl.lon_ctl.u.q_ref, y_kin(ac).ω_wb_b[2]; atol = 1e-3)
        @test all(isapprox.(y_air(ac).EAS, ctl.u.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end


    ################################ lon_EAS_q #################################

    @testset verbose = true "lon_EAS_θ" begin

        Sim.init!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_EAS_θ
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_EAS_θ

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 0.1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/6
        ctl.u.θ_ref = deg2rad(3)
        step!(sim, 10, true)
        ctl.u.θ_ref = -deg2rad(3)
        step!(sim, 60, true)

        @test isapprox(ctl.lon_ctl.u.θ_ref, y_kin(ac).e_nb.θ; atol = 1e-3)
        @test all(isapprox.(y_air(ac).EAS, ctl.u.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

    end

    ############################## lon_EAS_clm #################################

    @testset verbose = true "lon_EAS_clm" begin

        Sim.init!(sim, init_air)

        ctl.u.lon_ctl_mode_req = lon_EAS_clm
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.lon_ctl_mode === lon_EAS_clm

        #check the correct parameters are loaded and assigned to the controller
        c2θ_lookup = load_pid_lookup(joinpath(data_folder, "c2θ_lookup.h5"))
        k_p = c2θ_lookup(y_air(ac).EAS, Float64(y_kin(ac).h_e)).k_p
        @test all(isapprox.(ctl.y.lon_ctl.c2θ_pid.k_p, k_p; atol = 1e-6))

        #when trim reference values are kept, the control mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #correct tracking while turning
        ctl.u.φ_ref = π/6
        ctl.u.EAS_ref = 45
        ctl.u.clm_ref = 2
        step!(sim, 30, true)
        @test all(isapprox.(y_kin(ac).v_eb_n[3], -ctl.u.clm_ref; atol = 1e-1))
        @test all(isapprox.(y_air(ac).EAS, ctl.u.EAS_ref; atol = 1e-1))

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

        # return sim

    end #testset

    end #testset

    end #testset

end #function


function test_guidance_modes()

    @testset verbose = true "Guidance Modes" begin

    h_trn = HOrth(0.0)
    world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

    ac = world.ac
    ctl = ac.avionics.ctl

    init_air = C172.TrimParameters()
    dt = Δt = 0.01

    sim = Simulation(world; dt, Δt, t_end = 600)

    @testset verbose = true "Altitude Guidance" begin

        Sim.init!(sim, init_air)
        y_kin_trim = y_kin(ac)

        ctl.u.vrt_gdc_mode_req = vrt_gdc_alt
        ctl.u.lat_ctl_mode_req = lat_φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.vrt_gdc_mode === vrt_gdc_alt
        @test ctl.y.lon_ctl_mode === lon_EAS_clm

        #when trim reference values are kept, the guidance mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(ac).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(ac).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #all tests while turning
        ctl.u.φ_ref = π/12

        ctl.u.h_ref = y_kin_trim.h_e + 100
        step!(sim, 1, true)
        @test ctl.y.lon_ctl_mode === lon_thr_EAS
        step!(sim, 60, true) #altitude is captured
        @test ctl.y.lon_ctl_mode === lon_EAS_clm
        @test isapprox.(y_kin(ac).h_e - HEllip(ctl.u.h_ref), 0.0; atol = 1e-1)

        #reference changes within the current threshold do not prompt a mode change
        ctl.u.h_ref = y_kin(ac).h_e - ctl.alt_gdc.s.h_thr / 2
        step!(sim, 1, true)
        @test ctl.y.lon_ctl_mode === lon_EAS_clm
        step!(sim, 30, true) #altitude is captured
        @test isapprox.(y_kin(ac).h_e - HEllip(ctl.u.h_ref), 0.0; atol = 1e-1)

        ctl.u.h_ref = y_kin_trim.h_e - 100
        step!(sim, 1, true)
        @test ctl.y.lon_ctl_mode === lon_thr_EAS
        step!(sim, 80, true) #altitude is captured
        @test ctl.y.lon_ctl_mode === lon_EAS_clm
        @test isapprox.(y_kin(ac).h_e - HEllip(ctl.u.h_ref), 0.0; atol = 1e-1)

        @test ctl.y.lon_ctl_mode === lon_EAS_clm

        #must reset scheduling counter before standalone calls to f_disc!, but
        #without calling Sim.reinit! so that the controller state is preserved
        world._n[] = 0
        @test @ballocated(f_disc!($world)) == 0

        ctl.u.h_ref = y_kin_trim.h_e + 100
        step!(sim, 1, true)
        @test ctl.y.lon_ctl_mode === lon_thr_EAS

        #test for allocations in the current control mode
        @test @ballocated(f_disc!(NoScheduling(), $world)) == 0

        # return TimeSeries(sim)

    end

    end #testset

end

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
            vrt_gdc_mode_req = vrt_gdc_alt,
            lat_ctl_mode_req = lat_φ_β,
            φ_ref = φ_ref,
        )

        #therefore, these would also work
        # cmd = (
        #     vrt_gdc_mode_req = 1,
        #     lat_ctl_mode_req = 2,
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
    length(str) > 2 && JSON3.read!(str, world.ac.avionics.ctl.u)

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
    JSON3.read!(JSON3.write(world.ac.avionics.ctl.u, allow_inf=true), world.ac.avionics.ctl.u; allow_inf=true)

    Sim.run_interactive!(sim)

    save && save_plots(TimeSeries(sim).ac.vehicle.kinematics,
                        normpath("tmp/plots/test_c172x1/test_json_loopback/kin");
                        Plotting.defaults...)

    return nothing

end


function test_sim(; save::Bool = true)

    h_trn = HOrth(427.2);

    # on ground
    # initializer = KinInit(
    #     loc = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)),
    #     q_nb = REuler(deg2rad(157), 0, 0),
    #     h = h_trn + 1.81) |> C172.Init

    # on air, automatically trimmed
    initializer = C172.TrimParameters(
        Ob = Geographic(LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)), HEllip(650)))

    world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

    sim = Simulation(world; t_end = 30)
    Sim.init!(sim, initializer)

    Sim.run!(sim)

    save && save_plots(TimeSeries(sim).ac.vehicle.kinematics,
                        normpath("tmp/plots/test_c172x1/test_sim/kin");
                        Plotting.defaults...)

    return nothing

end

function test_sim_interactive(; save::Bool = true)

    h_trn = HOrth(427.2);

    # on ground
    initializer = C172.Init(KinInit(
        loc = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)),
        q_nb = REuler(deg2rad(157), 0, 0),
        h = h_trn + C172.Δh_to_gnd))

    # # on air, automatically trimmed
    # initializer = C172.TrimParameters(
    #     Ob = Geographic(LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)), HEllip(650)))

    trn = HorizontalTerrain(h_trn)
    ac = Cessna172Xv1(WA(), trn) |> Model;

    sim = Simulation(ac; dt = 1/60, Δt = 1/60, t_end = 1000)

    Sim.init!(sim, initializer)

    for joystick in update_connected_joysticks()
        Sim.attach!(sim, joystick)
    end

    xpc = XPlane12Control()
    # xpc = XPlane12Control(address = IPv4("192.168.1.2"))
    Sim.attach!(sim, xpc)

    Sim.run_interactive!(sim; pace = 1)

    save && save_plots(TimeSeries(sim).ac.vehicle.kinematics,
                        normpath("tmp/plots/test_c172x1/test_sim_interactive/kin");
                        Plotting.defaults...)

    return nothing

end



end #module