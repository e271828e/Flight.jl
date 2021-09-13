using Flight.Kinematics
using Flight.Atmosphere
using Flight.Airdata
using Flight.Utils

atm_sys = HybridSystem(AtmosphereCmp())
p = Geographic(alt = AltOrth(10000))
x_kin = get_x0(KinInit(v_eOb_b = [300, 10, 30], Ob = p, q_nb = Ry(0.01)))
dx_pos = copy(x_kin.pos)
y_kin = f_kin!(dx_pos, x_kin)

atm_sys.u.wind.v_ew_n[1] = 50
atm_sys.u.isa_.T_sl += 10 #ISA+10 day

AirData(atm_sys, y_kin) |> pwf