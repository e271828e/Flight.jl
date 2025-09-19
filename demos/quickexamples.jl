
using Flight
using ControlSystems, RobustAndOptimalControl, Plots, LaTeXStrings

function quickexample1()

    #set up nonlinear simulation
    world = Model(SimpleWorld(aircraft = Cessna172Sv0(NED())))
    sim = Simulation(world; t_end = 10)

    #initialize at default trim condition
    trim_params = C172.TrimParameters()
    init!(sim, trim_params)

    #advance simulation for one second from trim condition
    step!(sim, 1, true)
    #apply 10% elevator increment
    world.aircraft.vehicle.systems.act.u.elevator += 0.1
    #run to completion
    run!(sim)
    #extract pitch rate from simulation results
    ts = TimeSeries(sim);
    _, q_nonlinear, _ = get_components(ts.aircraft.vehicle.kinematics.ω_wb_b)

    #extract aircraft submodel and linearize it around the trim condition
    lss = linearize!(world.aircraft, trim_params)
    #convert to NamedStateSpace
    nss = named_ss(lss)
    #extract elevator to pitch rate linear SISO system
    e2q = nss[:q, :elevator]
    #simulate a 0.1 step input at t=1
    y, t, _, _ = lsim(e2q, (x, t)->[0.1]*(t>=1), 0:0.01:10)
    q_linear = vec(y)

    plot(q_nonlinear; plot_title = "Pitch Rate", label = "Nonlinear", ylabel=L"$q \ (rad/s)$") |> display
    plot!(t, q_linear, label = "Linear")

end


function quickexample2()

    #set up nonlinear simulation
    world = SimpleWorld(aircraft = Cessna172Sv0(NED())) |> Model
    sim = Simulation(world; t_end = 10)

    #initialize at default trim condition
    trim_params = C172.TrimParameters()
    init!(sim, trim_params)

    #advance simulation 1 second from trim condition
    step!(sim, 1, true)
    #apply 10% elevator increment
    world.aircraft.vehicle.systems.act.u.elevator += 0.1
    #run to completion
    run!(sim)
    #extract pitch angle TimeSeries from simulation results
    ts = TimeSeries(sim);
    θ_nonlinear = ts.aircraft.vehicle.kinematics.e_nb.θ

    #extract aircraft submodel and linearize it around the trim condition
    lss = linearize!(world.aircraft, trim_params)
    #convert to NamedStateSpace
    nss = named_ss(lss)
    #extract elevator to pitch angle linear SISO system
    e2θ = nss[:θ, :elevator]
    #simulate a 0.1 step input at t=1
    y, t, _, _ = lsim(e2θ, (x, t)->[0.1]*(t>=1), 0:0.01:10)

    #get perturbation Δθ around trim condition from linear simulation results
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

    init!(sim, trim_params)

    #set longitudinal control mode to airspeed + climb rate
    world.aircraft.avionics.ctl.lon.u.mode_req = C172XControl.ModeControlLon.EAS_clm
    #set lateral control mode to bank angle + sideslip
    world.aircraft.avionics.ctl.lat.u.mode_req = C172XControl.ModeControlLat.φ_β
    #set climb rate reference to 2.0 m/s
    world.aircraft.avionics.ctl.lon.u.clm_ref = 2.0
    #set bank angle reference to 30 degrees
    world.aircraft.avionics.ctl.lat.u.φ_ref = deg2rad(30)
    #set wind to 1 m/s North
    world.atmosphere.wind.u.N = 1.0

    run!(sim)
    ts = TimeSeries(sim)
    #generate kinematics data plots
    kin_plots = make_plots(ts.aircraft.vehicle.kinematics)
    display(kin_plots[:Ob_t3d])

end
