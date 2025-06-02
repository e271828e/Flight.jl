using Sockets
using Flight

function ex01(; ac::Cessna172 = Cessna172Xv1(),
                situation::Symbol = :ground,
                xp12_address = IPv4("127.0.0.1"),
                xp12_port = 49000,
                save::Bool = true)

    #3D position and geographic heading for Salzburg airport (LOWS) runway 15
    p_LOWS15 = Geographic(LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)), HOrth(427.2))
    ψ_LOWS15 = deg2rad(157)

    #horizontal terrain model with altitude matching LOWS runway 15
    trn = HorizontalTerrain(altitude = HOrth(p_LOWS15))

    #default atmospheric model
    atm = SimpleAtmosphere()

    #define world and build System for simulation
    world = SimpleWorld(ac, atm, trn) |> System

    if situation === :ground
        #initial condition specified through aircraft frame kinematics
        initializer = KinInit(;
            q_nb = REuler(ψ_LOWS15, 0, 0), #attitude with respect to NED frame
            loc = LatLon(p_LOWS15), #2D location
            h = HOrth(p_LOWS15) + C172.Δh_to_gnd, #altitude
            ω_wb_b = zeros(3), #angular velocity
            v_eb_n = zeros(3), #velocity
            )
    elseif situation === :air
        #initial condition specified by trim parameters
        initializer = C172.TrimParameters(;
            Ob = Geographic(LatLon(p_LOWS15), HEllip(650)), #3D position
            EAS = 50.0, #equivalent airspeed
            ψ_nb = ψ_LOWS15, #geographic heading
            ψ_wb_dot = 0.0, #turn rate
            flaps = 0.0, #flap setting
            γ_wb_n = 0.0, #flight path angle
            x_fuel = 0.5, #available fuel fraction
        )
    else
        error("Unknown situation: $situation")
    end

    #create a Simulation with specified integration step size and stop time
    sim = Simulation(world; dt = 0.01, t_end = 1000)

    Sim.init!(sim, initializer)

    xp = XPlane12Control(address = xp12_address, port = xp12_port)
    Sim.attach!(sim, xp)
    for joystick in update_connected_joysticks()
        isa(joystick, Joysticks.T16000M) && Sim.attach!(sim, joystick)
    end

    Sim.run_interactive!(sim)

    kin_plots = make_plots(TimeSeries(sim).ac.vehicle.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).ac.vehicle.airflow; Plotting.defaults...)
    dyn_plots = make_plots(TimeSeries(sim).ac.vehicle.dynamics; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "plots", "ex01", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "plots", "ex01", "air"))
    save && save_plots(dyn_plots, save_folder = joinpath("tmp", "plots", "ex01", "dyn"))

    return sim

end