module TestC172R

using Test
using UnPack
using BenchmarkTools

using Flight

export test_c172r

function test_c172r()
    @testset verbose = true "Cessna172R" begin

        @testset verbose = true "Performance" begin test_system() end
        @testset verbose = true "Simulation" begin test_sim_nrt(save = false) end
        @testset verbose = true "Trimming" begin test_trimming() end

    end
end

function test_system()

        trn = HorizontalTerrain();
        atm = System(Atmosphere());

        ac_LTF = System(Cessna172R(LTF()));
        ac_ECEF = System(Cessna172R(ECEF()));
        ac_NED = System(Cessna172R(NED()));

        #all three kinematics implementations must be supported, no allocations
        @test @ballocated(f_cont!($ac_LTF, $atm, $trn)) == 0
        @test @ballocated(f_disc!($ac_LTF)) == 0

        @test @ballocated(f_cont!($ac_ECEF, $atm, $trn)) == 0
        @test @ballocated(f_disc!($ac_ECEF)) == 0

        @test @ballocated(f_cont!($ac_NED, $atm, $trn)) == 0
        @test @ballocated(f_disc!($ac_NED)) == 0

    return nothing

end

function get_input_callback()

    inputs = init_joysticks() |> values |> collect
    # inputs = [XBoxController(),]

    let inputs = inputs

        function (u) #anonymous
            for input in inputs
                Input.update!(input)
                Input.assign!(u.avionics, input)
            end
        end

    end

end

function get_output_callback()

    # outputs = [XPInterface(host = IPv4("192.168.1.2"))] #Parsec
    outputs = [XPInterface(),]
    Output.init!.(outputs)

    let outputs = outputs

        function (y) #anonymous
            for output in outputs
                Output.update!(output, y.kinematics.pos)
            end
        end

    end

end

function test_sim_nrt(; save::Bool = true)

    h_trn = HOrth(608.55);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(Atmosphere());
    ac = System(Cessna172R());
    kin_init = KinematicInit(
        v_eOb_n = [30, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 2200.5);

    Aircraft.init!(ac, kin_init)
    ac.u.avionics.eng_start = true #engine start switch on

    atm.u.wind.v_ew_n[1] = 0

    callback! = let

        function (u, t, y, params)

            ac.u.avionics.throttle = 0.2
            ac.u.avionics.yoke_Δx = (t < 5 ? 0.2 : 0.0)
            ac.u.avionics.yoke_Δy = 0.0
            ac.u.avionics.pedals = 0.0
            ac.u.avionics.brake_left = 0
            ac.u.avionics.brake_right = 0

        end
    end

    sim = Simulation(ac; args_c = (atm, trn), t_end = 150, sim_callback = callback!)

    Sim.run!(sim)
    plots = make_plots(sim; Plotting.defaults...)
    save ? save_plots(plots, save_folder = joinpath("tmp", "nrt_sim_test")) : nothing

    return sim

end

function test_sim_rt(; save::Bool = true)

    h_trn = HOrth(608.55);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(Atmosphere());
    ac = System(Cessna172R());
    kin_init = KinematicInit(
        v_eOb_n = [30, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 100.5);

    Aircraft.init!(ac, kin_init)
    ac.u.avionics.eng_start = true #engine start switch on
    ac.u.avionics.throttle = 1

    callback! = let input_callback! = get_input_callback(), output_callback = get_output_callback()

        function (u, t, y, params)
            input_callback!(u)
            output_callback(y)
        end
    end

    sim = Simulation(ac;
        args_c = (atm, trn),
        t_end = 120,
        sim_callback = callback!,
        realtime = true,
        )


    Sim.run!(sim)
    plots = make_plots(sim; Plotting.defaults...)
    save ? save_plots(plots, save_folder = joinpath("tmp", "rt_sim_test")) : nothing

    return sim

end

function test_trimming()

    @testset verbose = true "θ Constraint" begin

        #precompute v_wOb_b
        α_a = 0.15
        β_a = -0.11
        TAS = 100
        v_wOa_a = Air.get_velocity_vector(TAS, α_a, β_a)
        v_wOb_b = C172R.f_ba.q(v_wOa_a)

        #set γ_wOb_n and φ_nb arbitrarily and compute θ_nb
        γ_wOb_n = -0.07 #set arbitrarily
        ψ_nb = 0.3 #inconsequential
        φ_nb = 0.7
        θ_nb = C172R.Trim.θ_constraint(; v_wOb_b, γ_wOb_n, φ_nb)

        #then construct e_nb, transform v_wOb_b to v_wOb_n, recompute γ_wOb_n
        #and check it matches the original value
        e_nb = REuler(ψ_nb, θ_nb, φ_nb)
        v_wOb_n = e_nb(v_wOb_b)
        γ_wOb_n_test = Attitude.inclination(v_wOb_n)

        @test γ_wOb_n_test ≈ γ_wOb_n
        @test @ballocated(C172R.Trim.θ_constraint(; v_wOb_b = $v_wOb_b, γ_wOb_n = $γ_wOb_n, φ_nb = $φ_nb)) === 0

    end

    @testset verbose = true "Assignment" begin

        ac = System(Cessna172R())

        state = C172R.Trim.State(;
            α_a = 0.08, φ_nb = 0.3, n_eng = 0.8,
            throttle = 0.61, yoke_x = 0.01, yoke_y = -0.025, pedals = 0.0)

        params = C172R.Trim.Parameters(;
            loc = LatLon(), h = HOrth(1000),
            ψ_nb = 0.2, TAS = 40.0, γ_wOb_n = 0.0, ψ_lb_dot = 0.2, β_a = 0.3,
            fuel = 0.5, mixture = 0.5, flaps = 0.0)

        atm = System(Atmosphere());
        # atm.u.wind.v_ew_n = [4, 2, 4]

        trn = HorizontalTerrain()

        C172R.Trim.assign!(ac, atm, trn, state, params)

        e_lb = e_nb = REuler(ac.y.kinematics.q_nb)
        v_wOb_n = e_nb(ac.y.airflow.v_wOb_b)

        @test e_nb.φ ≈ state.φ_nb
        @test ac.y.airframe.aero.α ≈ state.α_a
        @test ac.y.airframe.pwp.engine.ω == state.n_eng * ac.airframe.pwp.engine.params.ω_rated
        @test ac.u.avionics.throttle == state.throttle
        @test ac.u.avionics.yoke_x == state.yoke_x
        @test ac.u.avionics.yoke_y == state.yoke_y
        @test ac.u.avionics.pedals == state.pedals

        @test e_nb.ψ ≈ params.ψ_nb
        @test Attitude.inclination(v_wOb_n) ≈ params.γ_wOb_n atol = 1e-12
        @test ac.y.kinematics.common.ω_lb_b ≈ Attitude.ω(e_lb, [params.ψ_lb_dot, 0, 0])
        @test ac.y.airflow.TAS ≈ params.TAS
        @test ac.y.airframe.aero.β ≈ params.β_a
        @test ac.x.airframe.fuel[1] == params.fuel
        @test ac.u.avionics.mixture == params.mixture
        @test ac.u.avionics.flaps == params.flaps

        #setting α_filt = α and β_filt = β should have zeroed their derivatives
        @test ac.ẋ.airframe.aero.α_filt ≈ 0.0 atol = 1e-12
        @test ac.ẋ.airframe.aero.β_filt ≈ 0.0 atol = 1e-12

        @test (@ballocated C172R.Trim.assign!($ac, $atm, $trn, $state, $params))===0

    end

    @testset verbose = true "Optimization" begin

        ac = System(Cessna172R())
        atm = System(Atmosphere())
        trn = HorizontalTerrain() #zero orthometric altitude
        state = C172R.Trim.State()
        params = C172R.Trim.Parameters()

        f_target = C172R.Trim.get_target_function(ac, atm, trn, params)

        @assert @ballocated($f_target($state)) === 0

    end
end





end #module