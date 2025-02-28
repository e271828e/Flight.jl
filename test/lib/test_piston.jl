module TestPiston

using Test
using UnPack
using BenchmarkTools
using LinearAlgebra

using Flight.FlightCore
using Flight.FlightLib

#non-exported stuff
using Flight.FlightLib.Air: p_std, T_std
using Flight.FlightLib.Piston: PistonEngine, PistonThruster
using Flight.FlightLib.Piston: inHg2Pa, ft2m, h2δ, p2δ, ft2m, compute_π_ISA_pow
using Flight.FlightLib.Piston: eng_off, eng_starting, eng_running

export test_piston

################################################################################
################################################################################

@kwdef struct TestHarness{T <: PistonThruster} <: SystemDefinition
    thruster::T = PistonThruster()
end

@kwdef mutable struct TestHarnessU
    air_data::AirData = AirData()
    kin_data::KinData = KinData()
    fuel_available::Bool = true
end

Systems.U(::TestHarness) = TestHarnessU()
Systems.Y(sd::TestHarness) = (thruster = Systems.Y(sd.thruster),)

function Systems.f_ode!(harness::System{<:TestHarness})
    @unpack air_data, kin_data, fuel_available = harness.u
    f_ode!(harness.thruster, air_data, kin_data)
    Systems.update_y!(harness)
    #alternatively, harness.y = (thruster = thruster.y,)
end

function Systems.f_step!(harness::System{<:TestHarness})
    @unpack fuel_available = harness.u
    f_step!(harness.thruster, fuel_available)
end

Systems.f_disc!(::NoScheduling, ::System{<:TestHarness}) = nothing

################################################################################

function test_piston()
    @testset verbose = true "Piston" begin
        test_engine_lookup()
        test_engine_response()
        test_thruster_response()
    end
end

function test_engine_lookup()
    n_stall = 0.15
    n_cutoff = 1.4
    lookup = Piston.generate_lookup(; n_stall, n_cutoff)
    ω_rated = 2700
    P_rated = 200

    @testset verbose = true "Engine Lookup" begin

        @testset verbose = true "δ_wot" begin

            let δ_wot = lookup.δ_wot
                #these graphs have been retouched, so allow more leeway here
                @test δ_wot(1800/ω_rated, inHg2Pa(20)/p_std) ≈ (9500 |> ft2m |> h2δ) atol = 0.1
                @test δ_wot(2700/ω_rated, inHg2Pa(22)/p_std) ≈ (7000 |> ft2m |> h2δ) atol = 0.1
                @test δ_wot(2100/ω_rated, inHg2Pa(16)/p_std) ≈ (15250 |> ft2m |> h2δ) atol = 0.1
                @test δ_wot(2300/ω_rated, inHg2Pa(12)/p_std) ≈ (22000 |> ft2m |> h2δ) atol = 0.1
            end

        end #testset

        @testset verbose = true "π_std" begin

            let π_std = lookup.π_std
                @test π_std(1800/ω_rated, inHg2Pa(20)/p_std) * P_rated ≈ 71 atol = 1
                @test π_std(2050/ω_rated, inHg2Pa(24)/p_std) * P_rated ≈ 113 atol = 1
                @test π_std(2400/ω_rated, inHg2Pa(17)/p_std) * P_rated ≈ 85 atol = 1
                @test π_std(2400/ω_rated, inHg2Pa(28.8)/p_std) * P_rated ≈ 176 atol = 1
            #test some values here
            end

        end #testset

        @testset verbose = true "π_wot" begin

            let π_wot = lookup.π_wot
                #these graphs have been retouched, so allow more leeway here
                @test π_wot(1800/ω_rated, 3e3 |> ft2m |> h2δ) * P_rated ≈ 108 atol = 3
                @test π_wot(2300/ω_rated, 2.4e3 |> ft2m |> h2δ) * P_rated ≈ 153 atol = 3
                @test π_wot(2500/ω_rated, 10e3 |> ft2m |> h2δ) * P_rated ≈ 129 atol = 3
                @test π_wot(2000/ω_rated, 20e3 |> ft2m |> h2δ) * P_rated ≈ 65 atol = 3

            end

        end #testset

        @testset verbose = true "π_ISA_pow" begin

            π_ISA_pow = let lookup = lookup
                (n, μ, δ) -> compute_π_ISA_pow(lookup, n, μ, δ)
            end

            #at n_stall and below, power is zero regardless of MAP value
            @test π_ISA_pow(n_stall, 0, 1) ≈ 0
            @test π_ISA_pow(n_stall, lookup.μ_wot(n_stall, 1), 1) ≈ 0
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


function test_engine_response()

    @testset verbose = true "Engine Response" begin

        kin = KinInit(h = HEllip(), v_eb_n = [50, 0, 0]) |> KinData
        atm = AtmData()
        air = AirData(kin, atm)
        eng = PistonEngine() |> System

        eng.u.τ_load = -10
        eng.u.J_load = 0.1

        eng.x.ω = 100.0
        y_init = eng.y
        f_ode!(eng, air)
        @test eng.y != y_init #y must have been updated

        eng.x.ω = 0.0
        eng.s.state = eng_off
        f_ode!(eng, air)
        @test eng.y.τ_shaft == 0

        eng.u.start = true
        f_step!(eng)
        @test eng.s.state == eng_starting

        eng.x.ω = 0.9eng.constants.ω_idle
        f_step!(eng) #with ω <= ω_idle, engine won't leave the starting state
        @test eng.s.state == eng_starting
        f_ode!(eng, air)
        @test eng.y.τ_shaft > 0 #it should output the starter torque

        eng.x.ω = 1.1eng.constants.ω_idle
        f_step!(eng) #engine should start now
        @test eng.s.state == eng_running
        f_ode!(eng, air)

        #if we give it some throttle, we should get output power
        eng.u.throttle = 0.1
        f_ode!(eng, air)
        @test eng.y.τ_shaft > 0

        #commanded stop
        eng.s.state = eng_running
        eng.u.stop = true
        f_step!(eng) #engine should start now
        eng.u.stop = false
        @test eng.s.state == eng_off

        #stall stop
        eng.s.state = eng_running
        eng.x.ω = 0.95eng.constants.ω_stall
        f_step!(eng)
        @test eng.s.state == eng_off
        eng.x.ω = 1.1eng.constants.ω_idle
        eng.s.state = eng_running

        #without fuel, the engine should shut down
        f_step!(eng, false)
        @test eng.s.state == eng_off

        #and then fail to start, even above the required speed
        eng.u.start = true
        f_step!(eng, false)
        @test eng.s.state == eng_starting
        f_step!(eng, false)
        @test eng.s.state != eng_running

        #when fuel is available, the engine starts
        f_step!(eng)
        @test eng.s.state == eng_running

        @test @ballocated(f_ode!($eng, $air)) == 0
        @test @ballocated(f_step!($eng, true)) == 0

        return

    end #testset

end #function


function test_thruster_response()

    @testset verbose = true "Thruster Dynamics" begin

        #initialize auxiliary elements
        hrn = TestHarness() |> System

        sim = Simulation(hrn, t_end = 100)

        hrn.thruster.engine.u.start = true

        #take two steps for the start command to take effect
        step!(sim)
        step!(sim)
        @test sim.y.thruster.engine.state === eng_starting

        #give it a few seconds to get to stable idle RPMs
        step!(sim, 10, true)
        @test sim.y.thruster.engine.state === eng_running
        @test sim.y.thruster.engine.ω ≈ hrn.thruster.engine.constants.ω_idle atol = 1
        hrn.thruster.engine.u.start = false

        #thruster should be pushing
        @test get_wr_b(sim.sys.thruster).F[1] > 0
        #and receiving a CCW torque
        @test get_wr_b(sim.sys.thruster).τ[1] < 0


        #give it some throttle and see the RPMs increase
        hrn.thruster.engine.u.throttle = 1
        step!(sim, 5, true)
        @test sim.y.thruster.engine.ω > 2hrn.thruster.engine.constants.ω_idle

        #back to idle, make sure the idle controller kicks back in and idle RPMs
        #stabilize around their target value
        hrn.thruster.engine.u.throttle = 0
        step!(sim, 5, true)
        @test sim.y.thruster.engine.ω ≈ hrn.thruster.engine.constants.ω_idle atol = 1

        ########## Propeller sense and transmission gear ratio mismatch ########

        #coupling a CCW propeller with a positive gear ratio (non-inverting)
        #transmission, leads to the propeller turning in the wrong sense

        @test_throws AssertionError PistonThruster(
            propeller = Propeller(sense = Propellers.CCW),
            gear_ratio = 1)

        ################### Variable pitch CCW thruster ########################

        #CCW propeller should be coupled with negative gear ratio transmission
        hrn = PistonThruster(
            propeller = Propeller(Propellers.Lookup(Δβ_range = range(0, 0.2, length = 11)); sense = Propellers.CCW),
            gear_ratio = -1) |> TestHarness |> System

        sim = Simulation(hrn, t_end = 100)

        hrn.thruster.propeller.u[] = 0
        hrn.thruster.engine.u.start = true
        step!(sim, 5, true) #give it a few seconds to get to stable idle RPMs
        hrn.thruster.engine.u.start = false
        @test hrn.y.thruster.engine.state === eng_running

        @test hrn.y.thruster.propeller.ω ≈ -sim.y.thruster.engine.ω
        @test get_wr_b(sim.sys.thruster).F[1] > 0
        @test get_wr_b(sim.sys.thruster).τ[1] > 0 #CW opposing torque

        #change propeller pitch and check that the idle controller raises the
        #idle manifold pressure to hold the target idle RPMs
        hrn.thruster.propeller.u[] = 0.1
        step!(sim, 5, true)
        @test sim.y.thruster.engine.ω ≈ hrn.thruster.engine.constants.ω_idle atol = 1

        #with the throttle well above idle, the idle controller will not be
        #holding RPMs anymore
        hrn.thruster.engine.u.throttle = 0.5
        step!(sim, 5, true)

        #now, increasing pitch gives a higher thrust coefficient, but also a
        #higher torque coefficient, which drives down RPMs, which in turn
        #reduces absolute thrust. in general, whether absolute thrust increases
        #or decreases with propeller pitch depends on operating conditions.
        ω_tmp = sim.y.thruster.engine.ω
        hrn.thruster.propeller.u[] = 0.2
        step!(sim, 5, true)
        @test sim.y.thruster.engine.ω < ω_tmp

        #starved engine shuts down
        hrn.u.fuel_available = false
        step!(sim, 1, true)
        @test sim.y.thruster.engine.state == eng_off
        step!(sim, 5, true)

        #after a few seconds the engine should have stopped completely due to
        #friction
        @test sim.y.thruster.engine.ω ≈ 0.0 atol = 1e-10

        #start the engine again to test for allocations
        hrn.u.fuel_available = true
        hrn.thruster.engine.u.start = true
        step!(sim, 5, true)
        hrn.thruster.engine.u.start = false

        @test sim.y.thruster.engine.state === eng_running

        @test @ballocated(f_ode!($hrn)) == 0
        @test @ballocated(f_step!($hrn)) == 0

        # sim = Simulation(thr, args_ode = (air, kin), args_step = (fuel,), t_end = 100, save_on = false)
        # sim.u.engine.start = true
        # step!(sim, 5, true)
        # sim.u.engine.start = false
        # @show @ballocated(step!($sim, 0.1, true))


    end #testset

end #function

function plot_lookup(; plot_settings...)

    lookup = PistonEngine().lookup

    n_plot = range(0, 1.5, length = 100)
    δ_plot = range(1, 0, length = 100)
    @show μ_plot = range(0.1p_std, inHg2Pa(30), length = 10)/p_std

    π_std_plot = [lookup.π_std(n, μ) for (n, μ) in Iterators.product(n_plot, μ_plot)]
    # plot(μ_plot, π_std_plot'; plot_settings...)
    plot(n_plot, π_std_plot; plot_settings...)

    # π_wot_plot = [lookup.π_wot(n,p) for (n,p) in Iterators.product(n_plot, δ_plot)]
    # plot(δ_plot, π_wot_plot'; plot_settings...)
    # plot(n_plot, π_wot_plot; plot_settings...)


end #function

end #module