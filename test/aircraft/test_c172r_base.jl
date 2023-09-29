module TestC172RBase

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

using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.World

using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172RBase

export test_c172r_base


function test_c172r_base()
    @testset verbose = true "Cessna172RBase" begin

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

            ac_LTF = System(Cessna172RBase(LTF()));
            ac_ECEF = System(Cessna172RBase(ECEF()));
            ac_NED = System(Cessna172RBase(NED()));

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

    @testset verbose = true "Assignment" begin

        ac = System(Cessna172RBase())

        state = C172RBase.TrimState(;
            α_a = 0.08, φ_nb = 0.3, n_eng = 0.8,
            throttle = 0.61, aileron = 0.01, elevator = -0.025, rudder = 0.0)

        params = C172RBase.TrimParameters(;
            loc = LatLon(), h = HOrth(1000),
            ψ_nb = 0.2, TAS = 40.0, γ_wOb_n = 0.0, ψ_lb_dot = 0.2, θ_lb_dot = 0.2,
            β_a = 0.3, fuel = 0.5, mixture = 0.5, flaps = 0.0)

        env = SimpleEnvironment() |> System
        # env.atm.u.wind.v_ew_n = [4, 2, 4]

        C172RBase.assign!(ac, env, params, state)

        e_lb = e_nb = ac.y.kinematics.e_nb
        v_wOb_n = e_nb(ac.y.air.v_wOb_b)

        @test e_nb.φ ≈ state.φ_nb
        @test ac.y.airframe.aero.α ≈ state.α_a
        @test ac.y.airframe.pwp.engine.ω == state.n_eng * ac.airframe.pwp.engine.params.ω_rated
        @test ac.airframe.act.u.throttle == state.throttle
        @test ac.airframe.act.u.aileron_offset == state.aileron
        @test ac.airframe.act.u.elevator_offset == state.elevator
        @test ac.airframe.act.u.rudder_offset == state.rudder

        @test e_nb.ψ ≈ params.ψ_nb
        @test Attitude.inclination(v_wOb_n) ≈ params.γ_wOb_n atol = 1e-12
        @test ac.y.kinematics.common.ω_lb_b ≈ Attitude.ω(e_lb, [params.ψ_lb_dot, params.θ_lb_dot, 0])
        @test ac.y.air.TAS ≈ params.TAS
        @test ac.y.airframe.aero.β ≈ params.β_a
        @test ac.x.airframe.fuel[1] == params.fuel
        @test ac.airframe.act.u.mixture == params.mixture
        @test ac.airframe.act.u.flaps == params.flaps

        #setting α_filt = α and β_filt = β should have zeroed their derivatives
        @test ac.ẋ.airframe.aero.α_filt ≈ 0.0 atol = 1e-12
        @test ac.ẋ.airframe.aero.β_filt ≈ 0.0 atol = 1e-12

        @test (@ballocated C172RBase.assign!($ac, $env, $params, $state))===0

    end

    @testset verbose = true "Optimization" begin

        ac = System(Cessna172RBase())
        env = System(SimpleEnvironment())
        trim_params = C172RBase.TrimParameters()
        state = C172RBase.TrimState()

        f_target = C172RBase.get_target_function(ac, env, trim_params)

        @test @ballocated($f_target($state)) === 0

        success, _ = trim!(ac; env, trim_params)

        @test success

    end

    end #testset

end #function


function test_sim(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        h_trn = HOrth(608.55);

        ac = Cessna172RBase();
        env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
        world = SimpleWorld(ac, env) |> System;

        kin_init = KinematicInit(
            v_eOb_n = [30, 0, 0],
            ω_lb_b = [0, 0, 0],
            q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
            loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
            h = h_trn + 1.9 + 2200.5);

        init_kinematics!(world, kin_init)

        world.ac.airframe.act.u.eng_start = true #engine start switch on
        world.env.atm.wind.u.v_ew_n .= [0, 0, 0]

        sys_io! = let

            function (world)

                u_act = world.ac.airframe.act.u
                t = world.t[]

                u_act.throttle = 0.2
                u_act.aileron = (t < 5 ? 0.25 : 0.0)
                u_act.elevator = 0.0
                u_act.rudder = 0.0
                u_act.brake_left = 1
                u_act.brake_right = 1

            end
        end

        sim = Simulation(world; t_end = 300, sys_io!, adaptive = true)
        Sim.run!(sim, verbose = true)

        # plots = make_plots(sim; Plotting.defaults...)
        plots = make_plots(TimeHistory(sim).ac.kinematics; Plotting.defaults...)
        save && save_plots(plots, save_folder = joinpath("tmp", "test_c172r_base", "sim"))

        # return sim
        # return world

    end

end


function test_sim_paced(; save::Bool = true)

    h_trn = HOrth(601.55);

    ac = Cessna172RBase();
    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
    world = SimpleWorld(ac, env) |> System;

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.3),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0);

    init_kinematics!(world, kin_init)

    sim = Simulation(world; dt = 0.02, Δt = 0.02, t_end = 300)

    interfaces = Vector{IODevices.Interface}()
    for joystick in get_connected_joysticks()
        push!(interfaces, attach_io!(sim, joystick))
    end

    # xp = XPCDevice()
    xp = XPCDevice(host = IPv4("192.168.1.2"))
    push!(interfaces, attach_io!(sim, xp))

    @sync begin
        for interface in interfaces
            Threads.@spawn IODevices.start!(interface)
        end
        Threads.@spawn Sim.run_paced!(sim; pace = 1, verbose = true)
    end

    plots = make_plots(TimeHistory(sim).ac.kinematics; Plotting.defaults...)
    # plots = make_plots(TimeHistory(sim); Plotting.defaults...)
    save && save_plots(plots, save_folder = joinpath("tmp", "test_c172r_base", "sim_paced"))

    return nothing

end




end #module