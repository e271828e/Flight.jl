using Flight.Geodesy
using Flight.Kinematics
using Flight.Air
using Flight.Air

function test_airflow()

    atm_sys = System(Atmosphere())
    p = Geographic(h = HOrth(13000))
    x_kin = init_x(KinLTF(), KinInit(v_eOb_n = [300, 10, 0], Ob = p, q_nb = Ry(0.01)))
    dx_pos = copy(x_kin.pos)
    kin_data = f_kin!(dx_pos, x_kin)

    atm_sys.u.wind.v_ew_n[1] = 50
    atm_sys.u.air.T_sl += 10 #ISA+10 day

    b = @benchmarkable AirflowData($kin_data, $atm_sys); run(b)

end