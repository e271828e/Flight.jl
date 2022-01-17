using Flight

h_trn = AltOrth(811);
trn = HorizontalTerrain(altitude = h_trn);
atm = System(AtmosphereDescriptor());

aero = System(C172.Aero());
pwp = System(C172.Pwp());

kin = KinInit(
    v_eOb_b = [50, 0, 5],
    ω_lb_b = [0, 0, 0],
    q_nb = REuler(ψ = 0, θ = 0.0, φ = 0),
    Ob = Geographic(alt = h_trn + 0.85)) |> KinData;

air = AirData(kin, atm);

aero.u.e = 0
aero.u.a = 0.
aero.u.r = 0.
aero.u.f = 0

f_cont!(aero, pwp, air, kin, trn);

aero.y |> pwf