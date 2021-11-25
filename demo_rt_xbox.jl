using Flight

using Base.Iterators
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using BenchmarkTools
using GLFW


################################## AC2 ######################################

function demo_rt_xbox()

    println("Generalize this for non-rt")

    h_trn = AltOrth(609.5);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(AtmosphereDescriptor());
    ac = System(C172Aircraft());
    ac_mdl = Model(ac, (trn, atm); dt = 0.02, t_end = 30, adaptive = false, solver = RK4(), y_saveat = 0.02);

    kin_init = KinInit(v_eOb_b = [10, 0, 0],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = π, θ = 0.05, φ = 0.00),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.531818), λ = deg2rad(-3.574862)),
                                        h_trn + 0.85 + 0.0));

    atm.u.wind.v_ew_n[1] = 0
    ac.u.throttle = 0
    println("Map the throttle correctly!")

    init!(ac, kin_init)
    #if the model was instantiated in advance, we need this to change its internal state vector!
    reinit!(ac_mdl, ac.x)

    xp = XPInterface()
    disable_physics(xp)
    set_position(xp, kin_init)

    output_div = 1

    #add any joysticks connected before GLFW was imported
    init_joysticks()

    # #enable detection of further joystick connections or disconnections
    # GLFW.SetJoystickCallback(joystick_callback)
    # GLFW.PollEvents() #check for newly connected joysticks

    # window = GLFW.CreateWindow(640, 480, "GLFW Callback Test")
    # GLFW.MakeContextCurrent(window)

    t_wall = time()

    # for i in take(ac_mdl.integrator, 5)
    for i in ac_mdl.integrator

        #when using the integrator as an iterator, we don't need to step, it is done
        #automatically at the beginning of each iteration. therefore, the dt we need
        #for pacing the simulation is not the proposed dt for the next step, but the
        #one in the step the integrator has already taken.

        #retrieve the dt just taken by the integrator
        dt = ac_mdl.dt

        #compute the wall time corresponding to the newly updated simulation
        t_wall_next = t_wall + dt

        #busy wait until the wall time catches up
        while (time() < t_wall_next) end
        # println(time()-t_wall_next)

        t_wall = t_wall_next
        # println(t_wall - t_wall_0)

        if ac_mdl.success_iter % output_div == 0
            set_position(xp, ac_mdl.sys.y.kin.pos)
        end

        #try to catch potential joystick disconnections
        # GLFW.PollEvents()
        for joystick in values(connected_joysticks)
            update_joystick(joystick)
            Aircraft.assign_joystick_inputs!(ac, joystick)
            # println(joystick.axes.data)
            # println(ac.u.pedals)
            # println(ac.u.throttle)
            println(ac.u.brake_left)
            println(ac.u.brake_right)
        end

        # Swap front and back buffers
        GLFW.SwapBuffers(window)

        #when GLFW is implemented, build the model with t_end = Inf and use this to break
        # if !GLFW.WindowShouldClose(window)
        #     break
        # end

    end

    # plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12)
    # plots(ac_mdl; save_path = joinpath("tmp", "plots_demo_rt"), plot_settings...)


end