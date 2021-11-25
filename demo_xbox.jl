using GLFW

includet("src/input.jl")

function demo_xbox()

    window = GLFW.CreateWindow(640, 480, "GLFW Callback Test")
    GLFW.MakeContextCurrent(window)

    # Input callbacks
    GLFW.SetJoystickCallback(joystick_callback)
    GLFW.SetKeyCallback(window, (_, key, scancode, action, mods) -> begin
        name = GLFW.GetKeyName(key, scancode)
        if name === nothing
            println("scancode $scancode ", action)
        else
            println("key $name ", action)
        end
    end)

    #add any joysticks connected before GLFW was imported
    init_joysticks()
    #add any joysticks connected since GLFW was imported
    GLFW.PollEvents()

    t_start = 0
    t_end = 10
    dt = 1

    t_wall = time()

    t = t_start
    while t < t_end

        #this executes the callbacks for all the events received since the last call
        #to PollEvents
        GLFW.PollEvents()
        #update joystick inputs for callbacks
        # cursor_pos = GLFW.GetCursorPos(window) #placeholder for joystick reads
        for joystick in values(connected_joysticks)
            update(joystick)
            println(interface.axes.data)
        end
        # println(cursor_pos)

        #take a simulation step
        #step!()

        #update simulation time
        t += dt

        #update wall time
        t_wall_next = t_wall + dt
        t_wall = t_wall_next

        # Swap front and back buffers
        GLFW.SwapBuffers(window)

        #busy wait until the end of this frame
        while (time() < t_wall) & !GLFW.WindowShouldClose(window)

            #nothing to do here. to increase responsiveness we could take the chance
            #to keep polling asynchronous events, but this increases processor load
            #(why? this is already a busy wait)
            # GLFW.PollEvents()

        end
        # println(time()-t_wall)


    end

    GLFW.DestroyWindow(window);

end
