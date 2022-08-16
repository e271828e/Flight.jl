using Flight
using OrdinaryDiffEq

function drop_test(; save::Bool = true)

    h_trn = HOrth(608.55);

    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn)) |> System
    ac = System(Cessna172R());
    kin_init = KinematicInit(
        v_eOb_n = [30, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 0.5);

    Aircraft.init!(ac, kin_init)
    ac.u.avionics.eng_start = true #engine start switch on
    ac.u.avionics.brake_left = 1.
    ac.u.avionics.brake_right = 1.

    sim = Simulation(ac; args_ode = (env, ), algorithm = Tsit5(), t_end = 25, adaptive = true)
    Sim.run!(sim, verbose = true)
    # plots = make_plots(sim; Plotting.defaults...)
    plots = make_plots(TimeHistory(sim); Plotting.defaults...)
    save ? save_plots(plots, save_folder = joinpath("tmp", "drop_test", "Tsit5")) : nothing

    return sim

end