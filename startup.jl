using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra


air = Y(AirData())
thr = EThruster();
thr_sys = ContinuousSystem(thr);
thr_mdl = ContinuousModel(thr_sys, (air,))
SciMLBase.step!(thr_mdl)

g = AirframeComponentGroup(left = EThruster(), right = EThruster());
g_sys = ContinuousSystem(g);
g_mdl = ContinuousModel(g_sys, (air,))

trn = DummyTerrainModel()
atm = DummyAtmosphericModel()
ac = TestAircraft();

ac_sys = ContinuousSystem(ac);
ac_mdl = ContinuousModel(ac_sys, (trn, atm));
ac_mdl.sys.subsystems.pwp.u.left.throttle = 1 #the same
step!(ac_mdl)
b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)

ac_sys = ContinuousSystem(ac);
ac_mdl = ContinuousModel(ac_sys, (trn, atm); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0:0.1:100);
ac_mdl.sys.subsystems.pwp.u.left.throttle = 1 #the same
b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)

# mdl.u.pwp.left.throttle = 0
# mdl.u.pwp.right.throttle = 1
# # step!(mdl)

# mdl = Model.ContinuousModel(EThruster())