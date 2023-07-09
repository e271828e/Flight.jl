module TestC172Rv0

using Test
using UnPack
using BenchmarkTools
using Sockets

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

using Flight.FlightAircraft.Aircraft
using Flight.FlightAircraft.World

using Flight.FlightAircraft.C172R

export test_c172rv0


function test_c172rv0()
    @testset verbose = true "Cessna172Rv0" begin

        test_system_methods()
        test_trimming()
        test_sim(save = false)

    end
end

function test_system_methods()

        @testset verbose = true "System Methods" begin

            env = SimpleEnvironment() |> System

            loc = NVector()
            trn_data = TerrainData(env.trn, loc)
            kin_init = KinematicInit( h = trn_data.altitude + 1.8);

            ac_LTF = System(Cessna172Rv0(LTF()));
            ac_ECEF = System(Cessna172Rv0(ECEF()));
            ac_NED = System(Cessna172Rv0(NED()));

            init_kinematics!(ac_LTF, kin_init)
            init_kinematics!(ac_ECEF, kin_init)
            init_kinematics!(ac_NED, kin_init)

            f_ode!(ac_LTF, env)
            #make sure we are on the ground to ensure landing gear code coverage
            @test ac_LTF.y.airframe.ldg.left.strut.wow == true

            #all three kinematics implementations must be supported, no allocations
            @test @ballocated(f_ode!($ac_LTF, $env)) == 0
            @test @ballocated(f_step!($ac_LTF)) == 0
            @test @ballocated(f_disc!($ac_LTF, 1, $env)) == 0

            @test @ballocated(f_ode!($ac_ECEF, $env)) == 0
            @test @ballocated(f_step!($ac_ECEF)) == 0
            @test @ballocated(f_disc!($ac_ECEF, 1, $env)) == 0

            @test @ballocated(f_ode!($ac_NED, $env)) == 0
            @test @ballocated(f_step!($ac_NED)) == 0
            @test @ballocated(f_disc!($ac_NED, 1, $env)) == 0

        end

    return nothing

end


function test_trimming()

    @testset verbose = true "Trimming" begin

    @testset verbose = true "θ Constraint" begin

        #precompute v_wOb_b
        α_a = 0.15
        β_a = -0.11
        TAS = 100
        v_wOa_a = Atmosphere.get_velocity_vector(TAS, α_a, β_a)
        v_wOb_b = C172Rv0.Airframe.f_ba.q(v_wOa_a)

        #set γ_wOb_n and φ_nb arbitrarily and compute θ_nb
        γ_wOb_n = -0.07 #set arbitrarily
        ψ_nb = 0.3 #inconsequential
        φ_nb = 0.7
        θ_nb = C172Rv0.θ_constraint(; v_wOb_b, γ_wOb_n, φ_nb)

        #then construct e_nb, transform v_wOb_b to v_wOb_n, recompute γ_wOb_n
        #and check it matches the original value
        e_nb = REuler(ψ_nb, θ_nb, φ_nb)
        v_wOb_n = e_nb(v_wOb_b)
        γ_wOb_n_test = Attitude.inclination(v_wOb_n)

        @test γ_wOb_n_test ≈ γ_wOb_n
        @test @ballocated(C172Rv0.θ_constraint(; v_wOb_b = $v_wOb_b, γ_wOb_n = $γ_wOb_n, φ_nb = $φ_nb)) === 0

    end

    @testset verbose = true "Assignment" begin

        ac = System(Cessna172Rv0())

        state = C172Rv0.TrimState(;
            α_a = 0.08, φ_nb = 0.3, n_eng = 0.8,
            throttle = 0.61, aileron = 0.01, elevator = -0.025, rudder = 0.0)

        params = C172Rv0.TrimParameters(;
            loc = LatLon(), h = HOrth(1000),
            ψ_nb = 0.2, TAS = 40.0, γ_wOb_n = 0.0, ψ_lb_dot = 0.2, θ_lb_dot = 0.2,
            β_a = 0.3, fuel = 0.5, mixture = 0.5, flaps = 0.0)

        env = SimpleEnvironment() |> System
        # env.atm.u.wind.v_ew_n = [4, 2, 4]

        C172Rv0.assign!(ac, env, params, state)

        e_lb = e_nb = REuler(ac.y.kinematics.q_nb)
        v_wOb_n = e_nb(ac.y.air.v_wOb_b)

        @test e_nb.φ ≈ state.φ_nb
        @test ac.y.airframe.aero.α ≈ state.α_a
        @test ac.y.airframe.pwp.engine.ω == state.n_eng * ac.airframe.pwp.engine.params.ω_rated
        @test ac.u.airframe.act.throttle == state.throttle
        @test ac.u.airframe.act.aileron_trim == state.aileron
        @test ac.u.airframe.act.elevator_trim == state.elevator
        @test ac.u.airframe.act.rudder_trim == state.rudder

        @test e_nb.ψ ≈ params.ψ_nb
        @test Attitude.inclination(v_wOb_n) ≈ params.γ_wOb_n atol = 1e-12
        @test ac.y.kinematics.common.ω_lb_b ≈ Attitude.ω(e_lb, [params.ψ_lb_dot, params.θ_lb_dot, 0])
        @test ac.y.air.TAS ≈ params.TAS
        @test ac.y.airframe.aero.β ≈ params.β_a
        @test ac.x.airframe.fuel[1] == params.fuel
        @test ac.u.airframe.act.mixture == params.mixture
        @test ac.u.airframe.act.flaps == params.flaps

        #setting α_filt = α and β_filt = β should have zeroed their derivatives
        @test ac.ẋ.airframe.aero.α_filt ≈ 0.0 atol = 1e-12
        @test ac.ẋ.airframe.aero.β_filt ≈ 0.0 atol = 1e-12

        @test (@ballocated C172Rv0.assign!($ac, $env, $params, $state))===0

    end

    @testset verbose = true "Optimization" begin

        ac = System(Cessna172Rv0())
        env = System(SimpleEnvironment())
        trim_params = C172Rv0.TrimParameters()
        state = C172Rv0.TrimState()

        f_target = C172Rv0.get_target_function(ac, env, trim_params)

        @test @ballocated($f_target($state)) === 0

        success, _ = C172Rv0.trim!(ac; env, trim_params)

        @test success

    end

    end #testset

end #function


function test_sim(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        h_trn = HOrth(608.55);

        ac = Cessna172Rv0();
        env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
        world = SimpleWorld(ac, env) |> System;

        kin_init = KinematicInit(
            v_eOb_n = [30, 0, 0],
            ω_lb_b = [0, 0, 0],
            q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
            loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
            h = h_trn + 1.9 + 2200.5);

        init_kinematics!(world, kin_init)

        world.u.ac.airframe.act.eng_start = true #engine start switch on
        world.u.env.atm.wind.v_ew_n .= [0, 0, 0]

        sys_io! = let

            function (u, s, y, t, params)

                u.ac.airframe.act.throttle = 0.2
                u.ac.airframe.act.aileron = (t < 5 ? 0.25 : 0.0)
                u.ac.airframe.act.elevator = 0.0
                u.ac.airframe.act.rudder = 0.0
                u.ac.airframe.act.brake_left = 1
                u.ac.airframe.act.brake_right = 1

            end
        end

        sim = Simulation(world; t_end = 300, sys_io!, adaptive = true)
        Sim.run!(sim, verbose = true)

        # plots = make_plots(sim; Plotting.defaults...)
        plots = make_plots(TimeHistory(sim).ac.kinematics; Plotting.defaults...)
        save && save_plots(plots, save_folder = joinpath("tmp", "sim_test"))

        # return sim
        return world

    end

end


function test_sim_paced(; save::Bool = true)

    h_trn = HOrth(601.55);

    ac = Cessna172Rv0();
    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
    world = SimpleWorld(ac, env) |> System;

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.3),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0);

    init_kinematics!(world, kin_init)

    sim = Simulation(world; t_end = 300)

    interfaces = Vector{IODevices.Interface}()
    for joystick in get_connected_joysticks()
        push!(interfaces, attach_io!(sim, joystick))
    end

    # xp = XPCInterface()
    xp = XPCInterface(host = IPv4("192.168.1.2"))
    push!(interfaces, attach_io!(sim, xp))

    @sync begin
        for interface in interfaces
            Threads.@spawn IODevices.start!(interface)
        end
        Threads.@spawn Sim.run_paced!(sim; rate = 1, verbose = true)
    end

    plots = make_plots(TimeHistory(sim).ac.kinematics; Plotting.defaults...)
    # plots = make_plots(TimeHistory(sim); Plotting.defaults...)
    save && save_plots(plots, save_folder = joinpath("tmp", "paced_sim_test"))

    return nothing

end




end #module