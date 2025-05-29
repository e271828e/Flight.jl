**Installation notes on separate documentation section. As usual, recommend using Revise**

Explain what XPlane12Output does.


The ```FlightCore.Joysticks``` module provides a joystick interface layer via [SDL2_jll]
(https://github.com/JuliaBinaryWrappers/SDL2_jll.jl). Currently, Thrustmaster's T.16000M is the only
supported model. Adding support for other joysticks is relatively straightforward, but we are not be
getting into that in this example. If you happen to have a T.16000M available, you can plug it in
now and do the following:

Otherwise, don't worry. For this example, the SAS and autopilot modes will do most of the flying for
us, so you can easily control the aircraft through the GUI.

!!! note "X-Plane UDP Settings"

    By default, the ```XPlane12Output``` constructor assumes X-Plane 12 is running on your local machine
    and listening on port 49000. If this is not the case, you should provide the appropriate IP address
    and port via keyword arguments. For example:

    ```@repl
    using Flight
    using Sockets #bring IPv4 into scope
    xp = XPlane12Output(address = IPv4("192.168.1.2"), port = 49001)
    ```

!!! note

    The X-Plane 12 demo is time-limited to 15 minutes. After that, an on-screen pop-up will appear. To
    get rid of it, you can restart X-Plane and click on *Resume Last Flight*. You will then need to
    abort the Julia simulation and run it again (we will see how to do this shortly).
