using Flight

function benchmark_pi()

    sys = PICompensator{1}() |> System
    sys_init! = (s) -> nothing #inhibit warnings
    sim = Simulation(sys; t_end = 1, sys_init!)
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

    sim = Simulation(ac; args_ode = (env,), t_end = 100, adaptive = false, sys_init!)
    Sim.run!(sim)
    plots = make_plots(TimeHistory(sim).kinematics; Plotting.defaults...)
    save_plots(plots, save_folder = joinpath("tmp", "ac_benchmark"))

    b = @benchmarkable Sim.run!($sim) setup=reinit!($sim)
    return b

end

struct TestPICompensatorMapping <: AbstractInputMapping end

function Input.assign!(u::Essentials.PICompensatorU{1},
                       joystick::Joystick{XBoxController},
                       ::TestPICompensatorMapping)

    u.input .= get_axis_value(joystick, :right_analog_y)
    u.reset .= Input.was_released(joystick, :button_A)
    u.sat_enable .⊻= Input.was_released(joystick, :button_Y)

end

function test_input_pi()

    sys = PICompensator{1}(k_i = 0.5) |> System
    sim = Simulation(sys; t_end = 30)
    joysticks = get_connected_joysticks()
    input = InputManager(sim, joysticks, TestPICompensatorMapping())

    @sync begin
        Threads.@spawn Sim.run!(input; verbose = true)
        Threads.@spawn Sim.run!(sim; rate = 2, verbose = true)
    end

end