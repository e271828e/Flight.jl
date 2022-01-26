using Flight

using Base.Iterators
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using Sockets
using BenchmarkTools
using GLFW


function demo_rt_xbox()

    println("Generalize this for non-rt")

    h_trn = AltOrth(609.5);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(AtmosphereDescriptor());
    ac = System(BeaverDescriptor());
    ac_mdl = Model(ac, (trn, atm); dt = 0.02, t_end = 90, adaptive = false, solver = RK4(), y_saveat = 0.02);

    kin_init = KinInit(v_eOb_b = [0, 0, 0],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = 0, θ = 0.07π, φ = 0.00),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
                                        h_trn + 2.5 - 0.15 + 0.0));

    atm.u.wind.v_ew_n[1] = 0
    ac.u.controls.throttle = 0

    init!(ac, kin_init)
    #since the model was instantiated in advance, we need this to change its
    #internal initial condition
    reinit!(ac_mdl, ac.x)

    # xp = XPInterface()
    xp = XPInterface(host = IPv4("192.168.1.2"))
    disable_physics(xp)
    set_position(xp, kin_init)

    output_div = 1

    #add any joysticks connected before GLFW was imported
    init_joysticks()

    #enable detection of further joystick connections or disconnections
    GLFW.SetJoystickCallback(joystick_callback)
    GLFW.PollEvents() #check for newly connected joysticks

    # window = GLFW.CreateWindow(640, 480, "GLFW Callback Test")
    # GLFW.MakeContextCurrent(window)

    t_wall = time()

    # error()

    # for i in take(ac_mdl.integrator, 5)
    for i in ac_mdl.integrator

        #the integrator steps automatically at the beginning of each iteration

        #retrieve the dt just taken by the integrator
        dt = ac_mdl.dt

        #compute the wall time epoch corresponding to the simulation time epoch
        #we just reached
        t_wall_next = t_wall + dt

        #busy wait while wall time catches up
        while (time() < t_wall_next) end
        # println(time()-t_wall_next)

        t_wall = t_wall_next

        if ac_mdl.success_iter % output_div == 0
            set_position(xp, ac_mdl.sys.y.kinematics.pos)
        end

        for joystick in values(connected_joysticks)

            update_joystick(joystick)
            Input.assign_joystick_inputs!(ac, joystick)

            if ac_mdl.success_iter % 40 == 0
                pwf(ac.u.controls)
                # ac.y.kinematics.pos.h_o
                # println(ac.u.yoke_x_trim.val, " ",  ac.u.yoke_x.val, " ", ac.subsystems.afm.subsystems.aero.u.δa.val)
                # println(ac.u.yoke_y_trim.val, " ",  ac.u.yoke_y.val, " ", ac.subsystems.afm.subsystems.aero.u.δe.val)
                # println(ac.subsystems.afm.subsystems.aero.u.δr.val)
                # println(ac.subsystems.afm.subsystems.aero.u.δf.val)
            end

        end

        # Swap front and back buffers
        # GLFW.SwapBuffers(window)

        #when GLFW is implemented, build the model with t_end = Inf and use this to break
        # if !GLFW.WindowShouldClose(window)
        #     break
        # end

    end

    plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12)
    plots(ac_mdl; save_path = joinpath("tmp", "plots_demo_rt"), plot_settings...)


end