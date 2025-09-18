
using Flight
using ControlSystems, RobustAndOptimalControl, Plots, LaTeXStrings

function quickexample1()

    #create nonlinear simulation model
    world = SimpleWorld(aircraft = Cessna172Sv0(NED())) |> Model;
    sim = Simulation(world; t_end = 10)

    #define trim condition
    trim_params = C172.TrimParameters()

    #initialize simulation at trim condition
    Sim.init!(sim, trim_params)

    #advance the simulation for one second from trimmed condition
    Sim.step!(sim, 1, true)
    #apply 10% elevator increment
    world.aircraft.vehicle.systems.act.u.elevator += 0.1
    #run to completion
    Sim.run!(sim)

    #extract pitch rate from simulation results
    ts = TimeSeries(sim);
    _, q_nonlinear, _ = get_components(ts.aircraft.vehicle.kinematics.ω_wb_b)

    #linearize aircraft around trim condition and convert to NamedStateSpace
    nss = linearize!(world.aircraft, trim_params) |> named_ss
    #extract elevator to pitch rate SISO system
    e2q = nss[:q, :elevator]
    #simulate a 0.1 step input at t=1
    y, t, _, _ = lsim(e2q, (x, t)->[0.1]*(t>=1), 0:0.01:10)
    q_linear = vec(y)

    plot(q_nonlinear; plot_title = "Pitch Rate", label = "Nonlinear", ylabel=L"$q \ (rad/s)$") |> display
    plot!(t, q_linear, label = "Linear")

end


function quickexample2()

    world = SimpleWorld(aircraft = Cessna172Sv0(NED())) |> Model;
    sim = Simulation(world; t_end = 10)
    trim_params = C172.TrimParameters()
    Sim.init!(sim, trim_params)

    Sim.step!(sim, 1, true)
    world.aircraft.vehicle.systems.act.u.elevator += 0.1
    Sim.run!(sim)
    ts = TimeSeries(sim);
    θ_nonlinear = ts.aircraft.vehicle.kinematics.e_nb.θ

    lss = linearize!(world.aircraft, trim_params)
    nss = named_ss(lss)
    e2θ = nss[:θ, :elevator]
    y, t, _, _ = lsim(e2θ, (x, t)->[0.1]*(t>=1), 0:0.01:10)

    #linear model simulation returns the perturbation Δθ around trim condition
    Δθ_linear = vec(y)
    #retrieve trim θ value from linearized system
    θ_trim = lss.y0[:θ]
    #compute total θ
    θ_linear = θ_trim .+ Δθ_linear

    plot(θ_nonlinear; plot_title = "Pitch Angle", label = "Nonlinear", ylabel=L"$\theta \ (rad)$") |> display
    plot!(t, θ_linear; label = "Linear")

end


function quickexample3()

    world = SimpleWorld(aircraft = Cessna172Xv2()) |> Model;
    sim = Simulation(world; t_end = 600)
    trim_params = C172.TrimParameters()

    Sim.init!(sim, trim_params)

    #set longitudinal control mode to EAS + climb rate
    world.aircraft.avionics.ctl.lon.u.mode_req = C172XControl.ModeControlLon.EAS_clm
    #set lateral control mode to bank angle + sideslip
    world.aircraft.avionics.ctl.lat.u.mode_req = C172XControl.ModeControlLat.φ_β
    #set climb rate reference to 2.0 m/s
    world.aircraft.avionics.ctl.lon.u.clm_ref = 2.0
    #set bank angle reference to 30 degrees
    world.aircraft.avionics.ctl.lat.u.φ_ref = deg2rad(30)
    #set wind to 1 m/s North
    world.atmosphere.wind.u.N = 1.0

    Sim.run!(sim)
    ts = TimeSeries(sim)
    #generate kinematics data plots
    kin_plots = make_plots(ts.aircraft.vehicle.kinematics)
    display(kin_plots[:Ob_t3d])

end
