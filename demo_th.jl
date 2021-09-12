using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using BenchmarkTools

using StructArrays

using Flight.Kinematics
using Flight.Plotting

air = AirData()

thr = EThruster();
thr_sys = HybridSystem(thr);
f_cont!(thr_sys, air);
y_thr = thr_sys.y
thr_mdl = HybridModel(thr_sys, (air,))
thr_mdl.u.throttle = 1
step!(thr_mdl, 10, true)

g = ACGroup(left = EThruster(), right = EThruster());
g_sys = HybridSystem(g);
f_cont!(g_sys, air);
y_g = g_sys.y
g_mdl = HybridModel(g_sys, (air,))
g_mdl.u.left.throttle = 1
g_mdl.u.right.throttle = 1
step!(g_mdl)
