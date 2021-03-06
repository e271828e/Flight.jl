module TestPiston

using Test
using Plots
using UnPack
using BenchmarkTools
using LinearAlgebra
using SciMLBase

using Flight
using Flight.Piston: Engine, Thruster, MagicFuelSupply
using Flight.Piston: inHg2Pa, ft2m, h2δ, p2δ, ft2m, compute_π_ISA_pow
using Flight.Piston: eng_off, eng_starting, eng_running
using Flight.Air
using Flight.Air: p_std, T_std

export test_piston

function test_piston()
    test_engine_dataset()
    test_engine_dynamics()
    test_thruster_dynamics()
end

function test_engine_dataset()
    n_stall = 0.15
    n_cutoff = 1.4
    dataset = Piston.generate_dataset(; n_stall, n_cutoff)
    ω_rated = 2700
    P_rated = 200

    @testset verbose = true "EngineDataset" begin

        @testset verbose = true "δ_wot" begin

            let δ_wot = dataset.δ_wot
                #these graphs have been retouched, so allow more leeway here
                @test δ_wot(1800/ω_rated, inHg2Pa(20)/p_std) ≈ (9500 |> ft2m |> h2δ) atol = 0.1
                @test δ_wot(2700/ω_rated, inHg2Pa(22)/p_std) ≈ (7000 |> ft2m |> h2δ) atol = 0.1
                @test δ_wot(2100/ω_rated, inHg2Pa(16)/p_std) ≈ (15250 |> ft2m |> h2δ) atol = 0.1
                @test δ_wot(2300/ω_rated, inHg2Pa(12)/p_std) ≈ (22000 |> ft2m |> h2δ) atol = 0.1
            end

        end #testset

        @testset verbose = true "π_std" begin

            let π_std = dataset.π_std
                @test π_std(1800/ω_rated, inHg2Pa(20)/p_std) * P_rated ≈ 71 atol = 1
                @test π_std(2050/ω_rated, inHg2Pa(24)/p_std) * P_rated ≈ 113 atol = 1
                @test π_std(2400/ω_rated, inHg2Pa(17)/p_std) * P_rated ≈ 85 atol = 1
                @test π_std(2400/ω_rated, inHg2Pa(28.8)/p_std) * P_rated ≈ 176 atol = 1
            #test some values here
            end

        end #testset

        @testset verbose = true "π_wot" begin

            let π_wot = dataset.π_wot
                #these graphs have been retouched, so allow more leeway here
                @test π_wot(1800/ω_rated, 3e3 |> ft2m |> h2δ) * P_rated ≈ 108 atol = 3
                @test π_wot(2300/ω_rated, 2.4e3 |> ft2m |> h2δ) * P_rated ≈ 153 atol = 3
                @test π_wot(2500/ω_rated, 10e3 |> ft2m |> h2δ) * P_rated ≈ 129 atol = 3
                @test π_wot(2000/ω_rated, 20e3 |> ft2m |> h2δ) * P_rated ≈ 65 atol = 3

            end

        end #testset

        @testset verbose = true "π_ISA_pow" begin

            π_ISA_pow = let dataset = dataset
                (n, μ, δ) -> compute_π_ISA_pow(dataset, n, μ, δ)
            end

            #at n_stall and below, power is zero regardless of MAP value
            @test π_ISA_pow(n_stall, 0, 1) ≈ 0
            @test π_ISA_pow(n_stall, dataset.μ_wot(n_stall, 1), 1) ≈ 0
            @test π_ISA_pow(0.5*n_stall, 0.5, 1) ≈ 0

            #as soon as n rises above n_stall, power starts increasing with MAP
            @test π_ISA_pow(1.5*n_stall, 0.5, 1) > π_ISA_pow(1.5*n_stall, 0.3, 1)

            #sanity checks against IO360 performance charts
            @test 71 <  π_ISA_pow(1800/ω_rated, inHg2Pa(20)/p_std, 3e3 |> ft2m |> h2δ) * P_rated < 84
            @test 131 < π_ISA_pow(2310/ω_rated, inHg2Pa(23.6)/p_std, 2.4e3 |> ft2m |> h2δ) * P_rated < 139
            @test 102 < π_ISA_pow(2500/ω_rated, inHg2Pa(18)/p_std, 10e3 |> ft2m |> h2δ) * P_rated < 119

        end #testset

    end #testset

end #function


function test_engine_dynamics()

    @testset verbose = true "EngineDynamics" begin

        kin = Kinematics.Initializer(v_eOb_n = [50, 0, 0]) |> Kinematics.Common
        atm = Atmosphere() |> System
        air = AirflowData(kin, atm)
        eng = Engine() |> System
        fuel = System(MagicFuelSupply())

        ω = 100.0
        y_init = eng.y
        f_cont!(eng, air, ω)
        @test eng.y != y_init #y must have been updated

        ω = 0.0
        eng.d.state = eng_off
        f_cont!(eng, air, ω)
        @test eng.y.M == 0

        eng.u.start = true
        f_disc!(eng, fuel, ω)
        @test eng.d.state == eng_starting

        ω = 0.9eng.idle.params.ω_target
        f_disc!(eng, fuel, ω) #with ω <= ω_target, engine won't leave the starting state
        @test eng.d.state == eng_starting
        f_cont!(eng, air, ω)
        @test eng.y.M > 0 #it should output the starter torque

        ω = 1.1eng.idle.params.ω_target
        f_disc!(eng, fuel, ω) #engine should start now
        @test eng.d.state == eng_running
        f_cont!(eng, air, ω)

        #engine will not generate torque because the idle controller's state is
        #initialized to 0, and ω is currently above ω_target, so μ_ratio_idle is
        #set to 0, and therefore idle power is also 0
        @test eng.y.M == 0

        #if we give it some throttle, we should get output power
        eng.u.thr = 0.1
        f_cont!(eng, air, ω)
        @test eng.y.M > 0

        #commanded shutdown
        eng.d.state = eng_running
        eng.u.shutdown = true
        f_disc!(eng, fuel, ω) #engine should start now
        eng.u.shutdown = false
        @test eng.d.state == eng_off

        #stall shutdown
        eng.d.state = eng_running
        ω = 0.95eng.params.ω_stall
        f_disc!(eng, fuel, ω) #engine should start now
        @test eng.d.state == eng_off
        ω = 1.1eng.idle.params.ω_target
        eng.d.state = eng_running

        #without fuel, the engine should shut down
        fuel.u[] = false
        f_disc!(eng, fuel, ω)
        @test eng.d.state == eng_off

        #and then fail to start, even above the required speed
        eng.u.start = true
        f_disc!(eng, fuel, ω)
        @test eng.d.state == eng_starting
        f_disc!(eng, fuel, ω)
        @test eng.d.state != eng_running

        #when fuel is available, the engine starts
        fuel.u[] = true
        f_disc!(eng, fuel, ω)
        @test eng.d.state == eng_running

        @test @ballocated(f_cont!($eng, $air, $ω)) == 0
        @test @ballocated(f_disc!($eng, $fuel, $ω)) == 0

        return eng

    end #testset

end #function

function test_thruster_dynamics()

    @testset verbose = true "NewThrusterDynamics" begin

        #initialize auxiliary elements
        kin = Kinematics.Initializer(v_eOb_n = [0, 0, 0]) |> Kinematics.Common
        atm = Atmosphere() |> System
        air = AirflowData(kin, atm)
        fuel = System(MagicFuelSupply())
        thr = Thruster() |> System
        sim = Simulation(thr, args_c = (air, kin), args_d = (fuel,), t_end = 100)

        sim.u.engine.start = true

        #take two steps for the start command to take effect
        step!(sim)
        step!(sim)
        @test sim.y.engine.state === Piston.eng_starting

        #give it a few seconds to get to stable idle RPMs
        step!(sim, 5, true)
        @test sim.y.engine.state === Piston.eng_running
        @test sim.y.engine.ω ≈ thr.engine.idle.params.ω_target atol = 1
        sim.u.engine.start = false

        # @test sim.y.transmission.ΔP ≈ 0 atol = 1e-8

        #thruster should be pushing
        @test get_wr_b(sim.sys).F[1] > 0
        #and receiving a CCW torque
        @test get_wr_b(sim.sys).M[1] < 0

        #give it some throttle and see the RPMs increase
        sim.u.engine.thr = 1
        step!(sim, 5, true)
        @test sim.y.engine.ω > 2thr.engine.idle.params.ω_target

        #back to idle, make sure the idle controller kicks back in and idle RPMs
        #stabilize around their target value
        sim.u.engine.thr = 0
        step!(sim, 5, true)
        @test sim.y.engine.ω ≈ thr.engine.idle.params.ω_target atol = 1


        ########## Propeller sense and transmission gear ratio mismatch ########

        #coupling a CCW propeller with a positive gear ratio (non-inverting)
        #transmission, leads to the propeller turning in the wrong sense

        @test_throws AssertionError Thruster(
            propeller = Propeller(sense = Propellers.CCW),
            n = 1)


        ################### Variable pitch CCW thruster ########################

        #CCW propeller should be coupled with negative gear ratio transmission
        thr = Thruster(
            propeller = Propeller(pitch = VariablePitch(), sense = Propellers.CCW),
            n = -1
        ) |> System

        sim = Simulation(thr, args_c = (air, kin), args_d = (fuel,), t_end = 100)

        sim.u.propeller[] = 0
        sim.u.engine.start = true
        step!(sim, 5, true) #give it a few seconds to get to stable idle RPMs
        sim.u.engine.start = false
        @test sim.y.engine.state === Piston.eng_running

        @test sim.y.propeller.ω ≈ -sim.y.engine.ω
        @test get_wr_b(sim.sys).F[1] > 0
        @test get_wr_b(sim.sys).M[1] > 0 #CW opposing torque

        #change propeller pitch and check that the idle controller raises the
        #idle manifold pressure to hold the target idle RPMs
        sim.u.propeller[] = 0.1
        step!(sim, 5, true)
        @test sim.y.engine.ω ≈ thr.engine.idle.params.ω_target atol = 1

        #with the throttle well above idle, the idle controller will not be
        #holding RPMs anymore
        sim.u.engine.thr = 0.5
        step!(sim, 5, true)

        #now, increasing pitch gives a higher thrust coefficient, but also a
        #higher torque coefficient, which drives down RPMs, which in turn
        #reduces absolute thrust. in general, whether absolute thrust increases
        #or decreases with propeller pitch depends on operating conditions.
        ω_tmp = sim.y.engine.ω
        sim.u.propeller[] = 0.2
        step!(sim, 5, true)
        @test sim.y.engine.ω < ω_tmp

        #starved engine shuts down
        fuel.u[] = false
        step!(sim, 1, true)
        @test sim.y.engine.state == Piston.eng_off
        step!(sim, 5, true)

        #after a few seconds the engine should have stopped completely due to
        #friction
        @test sim.y.engine.ω ≈ 0.0 atol = 1e-10

        #start the engine again to test for allocations
        sim.u.engine.start = true
        step!(sim, 5, true)
        sim.u.engine.start = false

        @test @ballocated(f_cont!($thr, $air, $kin)) == 0
        @test @ballocated(f_disc!($thr, $fuel)) == 0

        sim = Simulation(thr, args_c = (air, kin), args_d = (fuel,), t_end = 100, y_save_on = false)
        sim.u.engine.start = true
        step!(sim, 5, true)
        sim.u.engine.start = false
        # @show @ballocated(step!($sim, 0.1, true))

        # return sim

    end #testset

end #function

function plot_dataset(; plot_settings...)

    dataset = Engine().dataset

    n_plot = range(0, 1.5, length = 100)
    δ_plot = range(1, 0, length = 100)
    @show μ_plot = range(0.1p_std, inHg2Pa(30), length = 10)/p_std

    π_std_plot = [dataset.π_std(n, μ) for (n, μ) in Iterators.product(n_plot, μ_plot)]
    # plot(μ_plot, π_std_plot'; plot_settings...)
    plot(n_plot, π_std_plot; plot_settings...)

    # π_wot_plot = [dataset.π_wot(n,p) for (n,p) in Iterators.product(n_plot, δ_plot)]
    # plot(δ_plot, π_wot_plot'; plot_settings...)
    # plot(n_plot, π_wot_plot; plot_settings...)


end #function

end #module