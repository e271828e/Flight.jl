using Flight

function test_aero()

    #because EAS is altitude independent, the stall EAS will be roughly the same
    #at all altitudes
    h_trn = HOrth(-50); #avoid terrain interference when trimming at SL
    trn = HorizontalTerrain(altitude = h_trn)
    ac = Cessna172Sv0(WA(), trn) |> System;

    mid_cg_pld = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

    trim_params = C172.TrimParameters(
        Ob = Geographic(LatLon(), HOrth(0000)),
        EAS = 25.0,
        γ_wb_n = 0.0,
        x_fuel = 0.5,
        flaps = 1.0,
        payload = mid_cg_pld,
        )

    exit_flag, trim_state = Systems.init!(ac, trim_params)

    return trim_state

end