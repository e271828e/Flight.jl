using Sockets
using Flight
using OrdinaryDiffEq: Heun, RK4

using Flight.FlightAircraft.C172Y.C172YControl: ModeControlLon, ModeControlLat,
    ModeGuidanceLon, ModeGuidanceLat, FlightPhase

function tutorial02a()

    #3D position and geographic heading for Salzburg airport (LOWS) runway 15
    loc_LOWS15 = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997))
    h_LOWS15 = HOrth(427.2)
    ψ_LOWS15 = deg2rad(157)

    aircraft = Cessna172Yv1()
    terrain = HorizontalTerrain(h_LOWS15)
    atmosphere = SimpleAtmosphere()

    mdl = SimpleWorld(aircraft, atmosphere, terrain) |> Model

    initializer = C172.TrimParameters(;
        Ob = Geographic(loc_LOWS15, h_LOWS15 + 500), #500 m above LOWS runway 15
        EAS = 50.0, #equivalent airspeed
        ψ_nb = ψ_LOWS15, #geographic heading
        γ_wb_n = 0.0, #wind-relative flight path angle
        ψ_wb_dot = 0.0, #turn rate
        flaps = 0.0, #flap setting
        fuel_load = 0.5, #available fuel fraction
    )

    Modeling.init!(mdl, initializer)

    sim = Simulation(mdl; dt = 0.02, t_end = 60)

    Sim.step!(sim, 5)
    @show mdl.aircraft.avionics.ctl.u.elevator_offset = 0.1
    Sim.step!(sim, 2)
    @show mdl.aircraft.avionics.ctl.u.elevator_offset = -0.1
    Sim.step!(sim, 2)
    @show mdl.aircraft.avionics.ctl.u.elevator_offset = 0
    Sim.run!(sim)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/tutorial02/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/tutorial02/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/tutorial02/dyn"); Plotting.defaults...)

    return sim

end

function tutorial02b()

    #3D position and geographic heading for Salzburg airport (LOWS) runway 15
    loc_LOWS15 = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997))
    h_LOWS15 = HOrth(427.2)
    ψ_LOWS15 = deg2rad(157)

    aircraft = Cessna172Yv1()
    terrain = HorizontalTerrain(h_LOWS15)
    atmosphere = SimpleAtmosphere()

    mdl = SimpleWorld(aircraft, atmosphere, terrain) |> Model

    initializer = C172.TrimParameters(;
        Ob = Geographic(loc_LOWS15, h_LOWS15 + 500), #500 m above LOWS runway 15
        EAS = 50.0, #equivalent airspeed
        ψ_nb = ψ_LOWS15, #geographic heading
        γ_wb_n = 0.0, #wind-relative flight path angle
        ψ_wb_dot = 0.0, #turn rate
        flaps = 0.0, #flap setting
        fuel_load = 0.5, #available fuel fraction
    )

    user_callback! = function(mdl::Model)
        t = mdl.t[]
        u_ctl = mdl.aircraft.avionics.ctl.u
        if 5 <= t < 7
            u_ctl.elevator_offset = 0.1
        elseif 7 <= t < 9
            u_ctl.elevator_offset = -0.1
        else
            u_ctl.elevator_offset = 0
        end
    end

    sim = Simulation(mdl; dt = 0.02, t_end = 60, user_callback!)
    Sim.init!(sim, initializer)
    Sim.run!(sim)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/tutorial02/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/tutorial02/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/tutorial02/dyn"); Plotting.defaults...)

    return sim

end

function tutorial02c()

    #3D position and geographic heading for Salzburg airport (LOWS) runway 15
    loc_LOWS15 = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997))
    h_LOWS15 = HOrth(427.2)
    ψ_LOWS15 = deg2rad(157)

    aircraft = Cessna172Yv1()
    terrain = HorizontalTerrain(h_LOWS15)
    atmosphere = SimpleAtmosphere()

    mdl = SimpleWorld(aircraft, atmosphere, terrain) |> Model

    initializer = KinInit(;
        location = loc_LOWS15, #2D location
        h = h_LOWS15 + C172.Δh_to_gnd, #altitude
        q_nb = REuler(ψ_LOWS15, 0, 0), #attitude with respect to NED frame
        ω_wb_b = zeros(3), #angular velocity
        v_eb_n = zeros(3), #velocity
        ) |> C172.Init

    Modeling.init!(mdl, initializer)

    phase = Ref(:startup)
    user_callback! = let phase = phase, u = mdl.aircraft.avionics.ctl.u
        function (mdl::Model{<:SimpleWorld})
            if phase[] === :startup
                u.eng_start = true
                if mdl.y.aircraft.vehicle.systems.pwp.engine.state === Piston.EngineState.running
                    println("Engine started")
                    phase[] = :takeoff
                end
            end
        end
    end

    #create a Simulation with specified integration step size and stop time
    sim = Simulation(mdl; user_callback!, dt = 0.02, t_end = 100)

    Sim.step!(sim, 2, true)

    # save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/tutorial02/kin"); Plotting.defaults..., linewidth = 2,)
    # save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/tutorial02/air"); Plotting.defaults...)
    # save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/tutorial02/dyn"); Plotting.defaults...)

    return sim

end

function tutorial02d()

    #3D position and geographic heading for Salzburg airport (LOWS) runway 15
    loc_LOWS15 = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997))
    h_LOWS15 = HOrth(427.2)
    ψ_LOWS15 = deg2rad(157)

    aircraft = Cessna172Yv1()
    terrain = HorizontalTerrain(h_LOWS15)
    atmosphere = SimpleAtmosphere()

    mdl = SimpleWorld(aircraft, atmosphere, terrain) |> Model

    initializer = C172.TrimParameters(;
        Ob = Geographic(loc_LOWS15, h_LOWS15 + 500), #500 m above LOWS runway 15
        EAS = 50.0, #equivalent airspeed
        ψ_nb = ψ_LOWS15, #geographic heading
        γ_wb_n = 0.0, #wind-relative flight path angle
        ψ_wb_dot = 0.0, #turn rate
        flaps = 0.0, #flap setting
        fuel_load = 0.5, #available fuel fraction
    )

    #inspect model hierarchy here.
    Modeling.init!(mdl, initializer)

    sim = Simulation(mdl; dt = 0.02, t_end = 2000)

    Sim.step!(sim, 10, true)
    @show mdl.aircraft.avionics.ctl.u.lon_ctl_mode_req = lon_thr_EAS
    @show mdl.aircraft.avionics.ctl.u.throttle_axis = 1.0
    @show mdl.aircraft.avionics.ctl.u.EAS_ref = 40.0

    @show mdl.aircraft.avionics.ctl.u.lat_ctl_mode_req = lat_φ_β
    @show mdl.aircraft.avionics.ctl.u.φ_ref = deg2rad(30.0)
    @show mdl.aircraft.avionics.ctl.u.β_ref = deg2rad(0.0)

    # Sim.step!(sim, 60, true)
    Sim.step!(sim, 60)
    # Sim.run!(sim)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/tutorial02/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/tutorial02/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/tutorial02/dyn"); Plotting.defaults...)

    return sim

end