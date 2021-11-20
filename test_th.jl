using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using BenchmarkTools

using StructArrays

using Flight.Kinematics
using Flight.Plotting

atm_sys = System(AtmosphereCmp())
p = Geographic(alt = AltOrth(13000))
x_kin = init_x0(KinInit(v_eOb_b = [300, 10, 30], Ob = p, q_nb = Ry(0.01)))
dx_pos = copy(x_kin.pos)
kin_data = f_kin!(dx_pos, x_kin)
air_data = AirData(kin_data, atm_sys)



thr = EThruster();
thr_sys = System(thr);
f_cont!(thr_sys, kin_data, air_data);
y_thr = thr_sys.y
thr_mdl = Model(thr_sys, (kin_data, air_data,))
thr_mdl.u.throttle = 1
step!(thr_mdl, 10, true)

g = AirframeGroup(left = EThruster(), right = EThruster());
g_sys = System(g);
f_cont!(g_sys, kin_data, air_data);
y_g = g_sys.y
g_mdl = Model(g_sys, (kin_data, air,))
g_mdl.u.left.throttle = 1
g_mdl.u.right.throttle = 1
step!(g_mdl, 10, true)
