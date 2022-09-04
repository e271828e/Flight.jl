using Flight
using OrdinaryDiffEq

function benchmark_pi()

    sys = PICompensator{1}() |> System
    sys_reinit! = (s) -> nothing #avoid reinit! warnings
    sys_reinit!(sys) #not called automatically by the Simulation constructor
    sim = Simulation(sys; t_end = 100, sys_reinit!)
    b = @benchmarkable Sim.run!($sim) setup=reinit!($sim)
    return b

end

function benchmark_ac()

    ac = System(Cessna172R())
    env = System(SimpleEnvironment())
    kin_init = KinematicInit( v_eOb_n = [30, 0, 0], h = HOrth(1.8 + 2000))

    sys_reinit! = let kin_init = kin_init
        function (ac)
            Aircraft.init!(ac, kin_init)
            ac.u.avionics.eng_start = true
        end
    end

    sys_reinit!(ac) #not called automatically by the Simulation constructor
    sim = Simulation(ac; args_ode = (env,), t_end = 100, adaptive = true, sys_reinit!)
    Sim.run!(sim)
    plots = make_plots(TimeHistory(sim).kinematics; Plotting.defaults...)
    save_plots(plots, save_folder = joinpath("tmp", "ac_benchmark"))

    b = @benchmarkable Sim.run!($sim) setup=reinit!($sim)
    return b

end