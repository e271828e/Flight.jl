using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using BenchmarkTools

function ac_benchmark()

    h_trn = AltOrth()
    trn = HorizontalTerrain(altitude = h_trn)
    atm = System(AtmosphereDescriptor())
    ac = System(C172Aircraft());

    kin_init = KinInit( Ob = Geographic(LatLon(), h_trn + 0.9));

    init!(ac, kin_init)

    b = @benchmark f_cont!($ac, $trn, $atm)
    @show median(b)

    ac_mdl = Model(ac, (trn, atm); dt = 0.01, adaptive = false, solver = Heun(), y_saveat = 0.1);

    b = @benchmark step!($ac_mdl) setup=(reinit!($ac_mdl))
    @show median(b)

    b = @benchmark step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl))
    @show median(b)

    ac_mdl = Model(ac, (trn, atm); dt = 0.02, adaptive = false, solver = RK4(), y_saveat = 0.1);

    b = @benchmark step!($ac_mdl) setup=(reinit!($ac_mdl))
    @show median(b)

    b = @benchmark step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl))
    @show median(b)

    ac_mdl = Model(ac, (trn, atm)); #save at all steps

    b = @benchmark step!($ac_mdl) setup=(reinit!($ac_mdl))
    @show median(b)

    b = @benchmark step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl))
    @show median(b)

end

function ac_test02()

    h_trn = AltOrth(811)

    trn = HorizontalTerrain(altitude = h_trn)
    atm = System(AtmosphereDescriptor())
    ac = System(C172Aircraft());

    kin_init = KinInit(v_eOb_b = [1, 0, 0.1],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.531818), λ = deg2rad(-3.574862)),
                                        h_trn + 0.85));

    init!(ac, kin_init)
    ac.u.throttle = 0
    ac.u.brake_left = 0
    ac.u.brake_right = 0
    atm.u.wind.v_ew_n[1] = 0

    # f_cont!(ac, trn, atm)

    ac_mdl = Model(ac, (trn, atm); dt = 0.01, adaptive = false, solver = Heun());
    # reinit!(ac_mdl, tf = 20)
    step!(ac_mdl, 5, true)
    ac.u.pedals = 1
    step!(ac_mdl, 5, true)

    plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12)
    plots(ac_mdl; save_path = joinpath("tmp", "plots_test02"), plot_settings...)

    return nothing

end


function ac_test03()


end