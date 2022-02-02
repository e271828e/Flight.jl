using Flight

using Base.Iterators
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
    kin_init = KinInit(v_eOb_b = [40, 0, 0],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.0),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
                                        h_trn - 0 + 2.0));

    kin_data = KinData(kin_init)
    air_data = AirData(kin_data, atm)

    aero.u.e = 0.2
    aero.u.a = 0.0
    aero.u.r = 0.0
    aero.u.f = 0.333

    f_cont!(aero, pwp, air_data, kin_data, trn)

    aero.y |> pwf

    f_disc!(aero)

    aero.d |> pwf

end

function forward_drop_test()

    h_trn = AltOrth(601.3);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(AtmosphereDescriptor());
    ac = System(C182TDescriptor());
    kin_init = KinInit(v_eOb_b = [0, 0, 0],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.15),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
                                        h_trn + 2.0 - 0.15 + 0.0));

    Aircraft.init!(ac, kin_init)

    # #if the model was instantiated before setting the system's initial
    # condition, we need this to update the model's initial condition
    # reinit!(mdl, ac.x)

    sim = SimulationRun(
        # model = Model(ac, (trn, atm); t_end = 3600, adaptive = true),
        model = Model(ac, (trn, atm); t_end = 10, adaptive = false, solver = RK4(), dt = 0.02, y_saveat = 0.02),
        outputs = [XPInterface(host = IPv4("192.168.1.2"))], #Parsec
        # outputs = [XPInterface()], #localhost
        realtime = true,
        plot_enable = false,
        plot_path = joinpath("tmp", "c182t", "forward_drop")
        )

    Simulation.run!(sim)

end


function free_flight()

    h_trn = AltOrth(601.3);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(AtmosphereDescriptor());
    ac = System(C182TDescriptor());
    kin_init = KinInit(v_eOb_b = [0, 0, 0],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.15),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
                                        h_trn + 2.0 - 0.15 + 0.0));

    Aircraft.init!(ac, kin_init)

    # #if the model was instantiated before setting the system's initial
    # condition, we need this to update the model's initial condition
    # reinit!(mdl, ac.x)

    sim = SimulationRun(
        # model = Model(ac, (trn, atm); t_end = 3600, adaptive = true),
        model = Model(ac, (trn, atm); t_end = 60, adaptive = false, solver = RK4(), dt = 0.02, y_saveat = 0.02),
        inputs = init_joysticks() |> values |> collect,
        outputs = [XPInterface(host = IPv4("192.168.1.2"))], #Parsec
        # outputs = [XPInterface()], #localhost
        realtime = true,
        plot_enable = false,
        plot_path = joinpath("tmp", "beaver", "free_flight")
        )

    Simulation.run!(sim)

end