using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using BenchmarkTools

function ac_test01()

    trn = DummyTerrain()
    atm_sys = System(AtmosphereCmp())
    ac = AircraftBase(
        kin = KinLTF(),
        mass = TunableMass(),
        aero = SimpleDrag(),
        pwp = AirframeGroup((
            left = EThruster(motor = ElectricMotor(α = CW)),
            right = EThruster(motor = ElectricMotor(α = CCW)))),
    )
    ac_sys = System(ac);
    f_cont!(ac_sys, trn, atm_sys);
    y_ac = ac_sys.y

    ac_sys = System(ac); #should remake the system, because it sets the Model's initial condition upon creation
    ac_mdl = Model(ac_sys, (trn, atm_sys); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0.1);
    ac_sys.subsystems.pwp.u.left.throttle = 1 #the same
    b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)

    ac_sys = System(ac); #should remake the system, because it sets the Model's initial condition upon creation
    ac_mdl = Model(ac_sys, (trn, atm_sys));
    ac_sys.subsystems.pwp.u.left.throttle = 1 #the same
    b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)

    #kinematics torture test: do not try this without an adaptive integrator! the
    #"thing" shoots southward from the North pole. if you're looking exactly south,
    #your heading is ill-defined, so it chatters between π and -π at different
    #integration step. but the quaternion representation doesn't care. bank angle
    #wraps around -π twice longitude is not defined at the North pole, so it jumps
    #instantaneously from its random initialization value (0) to around π. it
    #doesn't chatter because due to Coriolis acceleration, it does not remain at π.

    #reinitialize at the North Pole, set end time and run to completion
    ac_sys.subsystems.pwp.u.left.throttle = 0. #the same
    ac_sys.subsystems.pwp.u.right.throttle = 0. #the same
    atm_sys.u.wind.v_ew_n[1] = -5
    reinit!(ac_mdl)
    x0 = copy(ac_mdl.x)
    init!(x0.kin, KinInit(v_eOb_b = [10, 0, 0], Ob = Geographic(LatLon(ϕ = 0), AltOrth(1000))))
    reinit!(ac_mdl, x0, tf = 10)
    solve!(ac_mdl)

    #this should be set at startup.jl
    plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12)

    plots(ac_mdl; save_path = joinpath("tmp", "plots_drag"), plot_settings...)

end

function ac_test02()

    h_trn = AltOrth(811)

    trn = HorizontalTerrain(altitude = h_trn)

    atm_sys = System(AtmosphereCmp())
    atm_sys.u.wind.v_ew_n[1] = 0

    ac = C172Aircraft()
    ac_sys = System(ac);

    kin_init = KinInit(v_eOb_b = [5, 0, 0],
                        ω_lb_b = [0, 0, 0],
                        q_nb = REuler(ψ = 0, θ = 0, φ = 0),
                        Ob = Geographic(LatLon(ϕ = deg2rad(40.531818), λ = deg2rad(-3.574862)),
                                        h_trn + 1));
    init!(ac_sys.x.kin, kin_init)
    ac_sys.subsystems.pwp.u.left.throttle = 1
    ac_sys.subsystems.pwp.u.right.throttle = 1
    ac_sys.subsystems.ldg.u.center.steering[] = 0
    ac_sys.subsystems.ldg.u.left.braking[] = 1
    ac_sys.subsystems.ldg.u.right.braking[] = 1

    @btime f_cont!($ac_sys, $trn, $atm_sys);

    ac_mdl = Model(ac_sys, (trn, atm_sys); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0.02);
    # b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)

    # reinit!(ac_mdl)
    # step!(ac_mdl, 5, true)
    # ac_sys.subsystems.ldg.u.left.braking[] = 0
    # ac_sys.subsystems.ldg.u.right.braking[] = 0
    # step!(ac_mdl, 5, true)
    solve!(ac_mdl)
    @show ac_mdl.sys.t

    plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12)

    plots(ac_mdl; save_path = joinpath("tmp", "plots_test02"), plot_settings...)

    # reinit!(ac_mdl)
    # x0 = copy(ac_mdl.x)
    # init!(x0.kin, KinInit(v_eOb_b = [10, 0, 0], Ob = Geographic(LatLon(ϕ = 0), AltOrth(1000))))
    # reinit!(ac_mdl, x0, tf = 10)
    # solve!(ac_mdl)

end