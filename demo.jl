using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using BenchmarkTools

using StructArrays

using Flight.Plotting

air = AirY()

thr = EThruster();
thr_sys = HybridSystem(thr);
y_thr = f_cont!(thr_sys, air);
thr_mdl = HybridModel(thr_sys, (air,))
thr_mdl.u.throttle = 1
step!(thr_mdl, 10, true)

t = thr_mdl.log.t
y_sa = StructArray(thr_mdl.log.saveval)
thplot(t, y_sa.hr_b)

g = ACGroup(left = EThruster(), right = EThruster());
g_sys = HybridSystem(g);
y_g = f_cont!(g_sys, air);
get_wr_b(y_g)
g_mdl = HybridModel(g_sys, (air,))
step!(g_mdl)

x_kin = x0(Kin())
dx_pos = copy(x_kin.pos)
y_kin = f_kin!(dx_pos, x_kin)
y_acc = AccY()
y_air = AirY()
y_ac = Aircraft.TestAircraftY(y_kin, y_acc, y_air, y_g)

trn = DummyTerrainModel()
atm = DummyAtmosphericModel()
ac = TestAircraft();

ac_sys = HybridSystem(ac);
y_ac = f_cont!(ac_sys, trn, atm);

ac_mdl = HybridModel(ac_sys, (trn, atm));
ac_mdl.sys.subsystems.pwp.u.left.throttle = 1 #the same
b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)

#this could set at startup.jl
plot_settings = (linewidth=2, plot_titlefontfamily="Times Roman", plot_titlefontsize = 20, legendfontsize=14, )
plots(ac_mdl; plot_settings...)

ac_mdl = HybridModel(ac_sys, (trn, atm); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0:0.1:100);
ac_mdl.sys.subsystems.pwp.u.left.throttle = 1 #the same
b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)