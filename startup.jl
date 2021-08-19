using Flight
using OrdinaryDiffEq
using LinearAlgebra


pa = ParametricAircraft();
trn = Aircraft.DummyTerrainModel();
atm = Aircraft.DummyAtmosphericModel();

mdl = Model.ContinuousModel(pa, (trn, atm); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0:0.1:100);
# mdl = Model.ContinuousModel(pa, (trn, atm); dt = 0.02, adaptive = false, method = RK4(), y_saveat = 0:0.1:100);
# mdl = ContinuousModel(pa);

mdl.u.pwp.left.throttle = 0
mdl.u.pwp.right.throttle = 1
# step!(mdl)

b = @benchmarkable step!($mdl,1, true) setup=(reinit!($mdl)); run(b)
# b = @benchmarkable step!($mdl,10, true) setup=(reinit!($mdl)); run(b)