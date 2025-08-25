using Sockets
using Flight

function tutorial01(; aircraft::Cessna172 = Cessna172Yv1(),
                situation::Symbol = :ground,
                xp12_address = IPv4("127.0.0.1"),
                xp12_port = 49000,
                )

    #3D position and geographic heading for Salzburg airport (LOWS) runway 15
    loc_LOWS15 = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997))
    h_LOWS15 = HOrth(427.2)
    ψ_LOWS15 = deg2rad(157)

    #horizontal terrain model with elevation matching LOWS runway 15
    terrain = HorizontalTerrain(h_LOWS15)

    #default atmospheric model
    atmosphere = SimpleAtmosphere()

    #define world and build Model for simulation
    world = SimpleWorld(aircraft, atmosphere, terrain) |> Model

    if situation === :ground
        #initial condition specified through aircraft frame kinematics
        initializer = KinInit(;
            location = loc_LOWS15, #2D location
            h = h_LOWS15 + C172.Δh_to_gnd, #altitude
            q_nb = REuler(ψ_LOWS15, 0, 0), #attitude with respect to NED frame
            ω_wb_b = zeros(3), #angular velocity
            v_eb_n = zeros(3), #velocity
            ) |> C172.Init

    elseif situation === :air
        #initial condition specified by trim parameters
        initializer = C172.TrimParameters(;
            Ob = Geographic(loc_LOWS15, h_LOWS15 + 500), #500 m above LOWS runway 15
            EAS = 50.0, #equivalent airspeed
            ψ_nb = ψ_LOWS15, #geographic heading
            γ_wb_n = 0.0, #wind-relative flight path angle
            ψ_wb_dot = 0.0, #turn rate
            flaps = 0.0, #flap setting
            fuel_load = 0.5, #available fuel fraction
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

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/tutorial01/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/tutorial01/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/tutorial01/dyn"); Plotting.defaults...)
    # save_plots(TimeSeries(sim).aircraft.vehicle.dynamics; Plotting.defaults...)

    return sim

end