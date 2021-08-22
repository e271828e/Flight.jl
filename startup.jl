using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra


air = AirDataY()

thr = EThruster();
thr_sys = HybridSystem(thr);
thr_mdl = HybridModel(thr_sys, (air,))
step!(thr_mdl)

g = AirframeComponentGroup(left = EThruster(), right = EThruster());
g_sys = HybridSystem(g);
g_mdl = HybridModel(g_sys, (air,))
step!(g_mdl)

trn = DummyTerrainModel()
atm = DummyAtmosphericModel()
ac = TestAircraft();

ac_sys = HybridSystem(ac);
ac_mdl = HybridModel(ac_sys, (trn, atm));
ac_mdl.sys.subsystems.pwp.u.left.throttle = 1 #the same
b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)

ac_sys = HybridSystem(ac);
ac_mdl = HybridModel(ac_sys, (trn, atm); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0:0.1:100);
ac_mdl.sys.subsystems.pwp.u.left.throttle = 1 #the same
b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)
