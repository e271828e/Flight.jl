module Robot2DDemos

using Plots, LaTeXStrings

using FlightCore
using FlightPhysics
using FlightApps
using FlightApps.Robot2D: Robot, InitParameters, mode_m, mode_v, mode_η

using FlightCore.GUI

export robot2d_sim

function robot2d_sim(gui::Bool = false)

    #instantiate the Robot (Vehicle and Controller) and set up a Simulation
    mdl = Model(Robot())
    sim = Simulation(mdl; t_end = 100, dt = 0.01, Δt = 0.02,
                            gui_settings = GUI.Settings(; theme = GUI.ColorTheme.light))

    #use the default initialization
    init_params = InitParameters()
    init!(sim, init_params)
    # mdl.controller.u.mode = Robot2D.mode_η
    # mdl.controller.u.η_ref = 5
    run!(sim; gui)
    ts = TimeSeries(sim)
    plot(ts.controller.v_ref; plot_title = "Velocity", label = "Command", ylabel=L"$v \ (m / s)$") |> display
    plot!(ts.vehicle.v; label = "Response")

end

end
