using Flight

function static_ground_test()

    h_trn = AltOrth(607.7);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(AtmosphereDescriptor());
    ac = System(C182TDescriptor());
    ac_mdl = Model(ac, (trn, atm); dt = 0.02, t_end = 90, adaptive = false, solver = RK4(), y_saveat = 0.02);

    kin_init = KinInit(v_eOb_b = [0, 0, 0],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = 0, θ = 0.07π, φ = 0.00),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
                                        h_trn + 2.5 - 0.15 + 0.0));

    init!(ac, kin_init)

end
