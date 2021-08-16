using Flight
using OrdinaryDiffEq
using LinearAlgebra


pa = ParametricAircraft();
# mdl = ContinuousModel(pa);

mdl = Model.ContinuousModel(pa; dt = 0.01, adaptive = false, method = Heun(),
y_saveat = 0:0.1:100);

mdl.u.pwp.left.throttle = 0.1
mdl.u.pwp.right.throttle = 0.1
# step!(mdl)

b = @benchmarkable step!($mdl,1, true) setup=(reinit!($mdl)); run(b)