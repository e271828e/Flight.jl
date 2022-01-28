using Flight
using GLFW

function test_xbox()

    #detect joysticks already connected
    init_joysticks()

    #enable detection of further joystick connections or disconnections
    GLFW.SetJoystickCallback(joystick_callback)
    GLFW.PollEvents() #check for newly connected joysticks

    run_time = 10.0
    dt = 0.5

    t_wall_0 = time()
    t_wall = t_wall_0
    while t_wall - t_wall_0 < run_time

        #compute the wall time corresponding to the newly updated simulation
        t_wall_next = t_wall + dt

        #busy wait until the wall time catches up
        while (time() < t_wall_next) end

        t_wall = t_wall_next

        #try to catch potential joystick disconnections
        GLFW.PollEvents()
        for joystick in values(connected_joysticks)
            Input.update!(joystick)
            println(get_axis_value(joystick))
            println()
            println(get_button_state(joystick, :button_A))
            println(get_button_change(joystick, :button_A))
        end

    end

end