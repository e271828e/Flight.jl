using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using BenchmarkTools


plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12)

trn = DummyTerrainModel()
atm_sys = System(AtmosphereCmp())

ac = TestAircraft(
    mass = ConstantMass(),
    ctrl = NoControlMapping(),
    aero = NullAirframeComponent(),
    pwp = NullAirframeComponent(),
    ldg = NullAirframeComponent(),
    srf = NullAirframeComponent()
)

ac_sys = System(ac); #should remake the system, because it sets the Model's initial condition upon creation
init!(ac_sys.x.kin, KinInit(v_eOb_b = [10, 0, 0], Ob = Geographic(LatLon(ϕ = 0), AltOrth(1000))))
ac_mdl = Model(ac_sys, (trn, atm_sys); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0.1);
solve!(ac_mdl)
plots(ac_mdl; save_path = joinpath("tmp", "plots_ac_01"), plot_settings...)

x0 = copy(ac_mdl.x)
init!(x0.kin, KinInit(v_eOb_b = [10, 0, 0], ω_lb_b = [0.5, 0.5, 0], Ob = Geographic(LatLon(ϕ = 0), AltOrth(1000))))
reinit!(ac_mdl, x0)
solve!(ac_mdl)
plots(ac_mdl; save_path = joinpath("tmp", "plots_ac_02"), plot_settings...)

################################## AC2 ######################################

ac = TestAircraft(
    kin = KinLTF(),
    mass = ConstantMass(),
    aero = SimpleDrag(),
    pwp = AirframeGroup((
        left = EThruster(motor = ElectricMotor(α = CW)),
        right = EThruster(motor = ElectricMotor(α = CCW)))),
)

atm_sys.u.wind.v_ew_n[1] = -5
ac_sys = System(ac);
init!(ac_sys.x.kin, KinInit(v_eOb_b = [10, 0, 0], Ob = Geographic(LatLon(ϕ = 0), AltOrth(1000))))
ac_mdl = Model(ac_sys, (trn, atm_sys); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0.1);
solve!(ac_mdl)
plots(ac_mdl; save_path = joinpath("tmp", "plots_ac_03"), plot_settings...)

ac_sys.subsystems.pwp.u.left.throttle = .1
ac_sys.subsystems.pwp.u.right.throttle = .1
reinit!(ac_mdl)
solve!(ac_mdl)
plots(ac_mdl; save_path = joinpath("tmp", "plots_ac_04"), plot_settings...)

ac_sys.subsystems.pwp.u.left.throttle = 0
ac_sys.subsystems.pwp.u.right.throttle = .2
reinit!(ac_mdl)
solve!(ac_mdl)
plots(ac_mdl; save_path = joinpath("tmp", "plots_ac_05"), plot_settings...)

b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)

# ac_sys = System(ac); #should remake the system, because it sets the Model's initial condition upon creation
# ac_mdl = Model(ac_sys, (trn, atm_sys));
# ac_sys.subsystems.pwp.u.left.throttle = 1 #the same
# b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)

# #kinematics torture test: do not try this without an adaptive integrator! the
# #"thing" shoots southward from the North pole. if you're looking exactly south,
# #your heading is ill-defined, so it chatters between π and -π at different
# #integration step. but the quaternion representation doesn't care. bank angle
# #wraps around -π twice longitude is not defined at the North pole, so it jumps
# #instantaneously from its random initialization value (0) to around π. it
# #doesn't chatter because due to Coriolis acceleration, it does not remain at π.

# #reinitialize at the North Pole, set end time and run to completion
# ac_sys.subsystems.pwp.u.left.throttle = 0. #the same
# ac_sys.subsystems.pwp.u.right.throttle = 0. #the same
# atm_sys.u.wind.v_ew_n[1] = -5
# reinit!(ac_mdl)
# x0 = copy(ac_mdl.x)
# init!(x0.kin, KinInit(v_eOb_b = [10, 0, 0], Ob = Geographic(LatLon(ϕ = 0), AltOrth(1000))))
# reinit!(ac_mdl, x0, tf = 10)
# solve!(ac_mdl)

# #this should be set at startup.jl
# plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12)

# plots(ac_mdl; save_path = joinpath("tmp", "plots_drag"), plot_settings...)