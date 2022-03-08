using Flight.Geodesy
using Flight.Kinematics
using Flight.Atmosphere
using Flight.Airdata

function test_airdata()

    atm_sys = System(AtmosphereDescriptor())
    p = Geographic(alt = AltOrth(13000))
    x_kin = init_x(KinLTF(), KinInit(v_eOb_b = [300, 10, 30], Ob = p, q_nb = Ry(0.01)))
    dx_pos = copy(x_kin.pos)
    kin_data = f_kin!(dx_pos, x_kin)

    atm_sys.u.wind.v_ew_n[1] = 50
    atm_sys.u.static.T_sl += 10 #ISA+10 day

    b = @benchmarkable AirData($kin_data, $atm_sys); run(b)

end