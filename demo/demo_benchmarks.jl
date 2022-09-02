using Flight
using OrdinaryDiffEq

function benchmark_pi()

    sys = PICompensator{1}() |> System
    sys_init! = (s) -> nothing #avoid reinit! warnings
    sim = Simulation(sys; t_end = 100, sys_init!)
    b = @benchmarkable Sim.run!($sim) setup=reinit!($sim)
    return b

end

function benchmark_ac()

    ac = System(Cessna172R())
    env = System(SimpleEnvironment())
    kin_init = KinematicInit( v_eOb_n = [30, 0, 0], h = HOrth(1.8 + 2000))
    ac.u.avionics.eng_start = true

    sys_init! = let kin_init = kin_init
        ac -> Aircraft.init!(ac, kin_init)
    end

    sim = Simulation(ac; args_ode = (env,), t_end = 100, adaptive = true, sys_init!)
    Sim.run!(sim)
    plots = make_plots(TimeHistory(sim).kinematics; Plotting.defaults...)
    save_plots(plots, save_folder = joinpath("tmp", "ac_benchmark"))

    b = @benchmarkable Sim.run!($sim) setup=reinit!($sim)
    return b

end