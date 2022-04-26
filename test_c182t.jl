using Flight

using Base.Iterators
using Dates
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using Sockets
using BenchmarkTools
using GLFW

function aerodynamics_test()

    h_trn = AltOrth(601.3);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(AtmosphereDescriptor());
    aero = System(C182T.Aero());
    pwp = System(C182T.Pwp());
    kin_init = KinInit(v_eOb_n = [40, 0, 0],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.0),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
                                        h_trn - 0 + 102.0));

    kin_data = KinData(kin_init)
    air_data = AirData(kin_data, atm)

    aero.u.e = 0.0
    aero.u.a = 0.0
    aero.u.r = 0.0
    aero.u.f = 0.0

    f_cont!(aero, pwp, air_data, kin_data, trn)

    aero.y |> pwf

    f_disc!(aero)

    aero.d |> pwf

end

function forward_drop_test()

    h_trn = AltOrth(608.55);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(AtmosphereDescriptor());
    ac = System(C182TDescriptor());
    kin_init = KinInit(v_eOb_n = [30, 0, 0],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = 0, θ = 0.2, φ = 0.0),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
                                        h_trn + 1.9 + 2220.5));

    Aircraft.init!(ac, kin_init)
    # #if the model was instantiated before setting the system's initial
    # condition, we need this to update the model's initial condition
    # reinit!(mdl, ac.x)

    ac.u.avionics.pedals = -.5
    ac.u.avionics.brake_left = 1
    ac.u.avionics.brake_right = 1
    ac.u.avionics.throttle = 0.3

    sim = SimulationRun(
        # model = Model(ac, (trn, atm); t_end = 3600, adaptive = true),
        model = Model(ac, (trn, atm); t_end = 30, adaptive = false, solver = RK4(), dt = 0.02, y_saveat = 0.02),
        # outputs = [XPInterface(host = IPv4("192.168.1.2"))], #Parsec
        outputs = [XPInterface()], #localhost
        realtime = false,
        )

    Simulation.run!(sim)
    # plots = make_plots(sim.model)
    # save_plots(plots, save_folder = joinpath("tmp", "drop_test", Dates.format(now(), "yyyy_mm_dd_HHMMSS")))

    return sim

end


function free_flight()

    h_trn = AltOrth(608.55);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(AtmosphereDescriptor());
    ac = System(C182TDescriptor());
    kin_init = KinInit(v_eOb_n = [0, 0, 0],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
                                        h_trn + 2.0 - 0.15 + 00.0));

    Aircraft.init!(ac, kin_init)

    # #if the model was instantiated before setting the system's initial
    # condition, we need this to update the model's initial condition
    # reinit!(mdl, ac.x)

    sim = SimulationRun(
        # model = Model(ac, (trn, atm); t_end = 3600, adaptive = true),
        model = Model(ac, (trn, atm); t_end = 100, adaptive = false, solver = RK4(), dt = 0.02, y_saveat = 0.02),
        inputs = init_joysticks() |> values |> collect,
        outputs = [XPInterface(host = IPv4("192.168.1.2"))], #Parsec
        # outputs = [XPInterface()], #localhost
        realtime = true,
        )

    Simulation.run!(sim)
    plots = make_plots(sim.model)
    save_plots(plots, save_folder = joinpath("tmp", "free_flight", Dates.format(now(), "yyyy_mm_dd_HHMMSS")))

    return sim

end