using HDF5
using Interpolations

function generate_data(fname = "src/aircraft/c172r/aero_dataset.h5")

    h5open(fname, "w") do fid

        ############################## C_D ##################################

        create_group(fid, "C_D")
        C_D = fid["C_D"]

        C_D["zero"] = 0.027

        create_group(C_D, "δe")
        C_D["δe"]["δe"] = [-1.0 0.0 1.0] |> vec
        C_D["δe"]["data"] = [ 0.06 0 0.06] |> vec

        create_group(C_D, "β")
        C_D["β"]["β"] = [-1.0 0.0 1.0] |> vec
        C_D["β"]["data"] = [ 0.17 0 0.17] |> vec

        create_group(C_D, "ge")
        C_D["ge"]["Δh_nd"] = [ 0.0000 0.1000 0.1500 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 1.1000 ] |> vec
        C_D["ge"]["data"] = [ 0.4800 0.5150 0.6290 0.7090 0.8150 0.8820 0.9280 0.9620 0.9880 1.0000 1.0000 1.0000 1.0000 ] |> vec

        create_group(C_D, "δf")
        C_D["δf"]["δf"] = deg2rad.([ 0.0000	10.0000 20.0000 30.0000 ]) |> vec
        C_D["δf"]["data"] = [ 0.0000 0.0070 0.0120 0.0180 ] |> vec

        create_group(C_D, "α_δf")
        C_D["α_δf"]["α"] = [ -0.0873 -0.0698 -0.0524 -0.0349 -0.0175 0.0000	0.0175	0.0349	0.0524	0.0698	0.0873 0.1047	0.1222	0.1396	0.1571	0.1745	0.1920	0.2094	0.2269	0.2443	0.2618	0.2793	0.2967	0.3142	0.3316	0.3491] |> vec
        C_D["α_δf"]["δf"] = deg2rad.([ 0.0000	10.0000	20.0000	30.0000 ]) |> vec
        C_D["α_δf"]["data"] = [ 0.0041	0.0000	0.0005	0.0014
                                0.0013	0.0004	0.0025	0.0041
                                0.0001	0.0023	0.0059	0.0084
                                0.0003	0.0057	0.0108	0.0141
                                0.0020	0.0105	0.0172	0.0212
                                0.0052	0.0168	0.0251	0.0299
                                0.0099	0.0248	0.0346	0.0402
                                0.0162	0.0342	0.0457	0.0521
                                0.0240	0.0452	0.0583	0.0655
                                0.0334	0.0577	0.0724	0.0804
                                0.0442	0.0718	0.0881	0.0968
                                0.0566	0.0874	0.1053	0.1148
                                0.0706	0.1045	0.1240	0.1343
                                0.0860	0.1232	0.1442	0.1554
                                0.0962	0.1353	0.1573	0.1690
                                0.1069	0.1479	0.1708	0.1830
                                0.1180	0.1610	0.1849	0.1975
                                0.1298	0.1746	0.1995	0.2126
                                0.1424	0.1892	0.2151	0.2286
                                0.1565	0.2054	0.2323	0.2464
                                0.1727	0.2240	0.2521	0.2667
                                0.1782	0.2302	0.2587	0.2735
                                0.1716	0.2227	0.2507	0.2653
                                0.1618	0.2115	0.2388	0.2531
                                0.1475	0.1951	0.2214	0.2351
                                0.1097	0.1512	0.1744	0.1866
        ]

        ############################## C_Y ##################################

        create_group(fid, "C_Y")
        C_Y = fid["C_Y"]

        C_Y["δr"] = 0.1870
        C_Y["δa"] = 0.0

        create_group(C_Y, "β_δf")
        C_Y["β_δf"]["β"] = [-0.3490 0 0.3490] |> vec
        C_Y["β_δf"]["δf"] = deg2rad.([0 30]) |> vec
        C_Y["β_δf"]["data"] = [
                            0.1370	0.1060
                            0.0000	0.0000
                            -0.1370	-0.1060
        ]
        create_group(C_Y, "p")
        C_Y["p"]["α"] = [0.0 0.094] |> vec
        C_Y["p"]["δf"] = deg2rad.([0 30]) |> vec
        C_Y["p"]["data"] = [
                            -0.0750	-0.1610
                            -0.1450	-0.2310
        ]
        create_group(C_Y, "r")
        C_Y["r"]["α"] = [0.0 0.094] |> vec
        C_Y["r"]["δf"] = deg2rad.([0 30]) |> vec
        C_Y["r"]["data"] = [
                            0.2140	0.1620
                            0.2670	0.2150
        ]


        ############################### C_L #################################

        create_group(fid, "C_L")
        C_L = fid["C_L"]

        C_L["δe"] = 0.4300
        C_L["q"] = 3.900
        C_L["α_dot"] = 1.700

        create_group(C_L, "ge")
        C_L["ge"]["Δh_nd"] = [ 0.0000 0.1000 0.1500 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 1.1000 ] |> vec
        C_L["ge"]["data"] = [ 1.2030 1.1270 1.0900 1.0730 1.0460 1.0550 1.0190 1.0130 1.0080 1.0060 1.0030 1.0020 1.0000 ] |> vec

        create_group(C_L, "α")
        C_L["α"]["α"] = [ -0.0900 0.0000	0.0900	0.1000	0.1200	0.1400	0.1600	0.1700	0.1900	0.2100	0.2400	0.2600	0.2800	0.3000	0.3200	0.3400	0.3600	] |> vec
        C_L["α"]["stall"] = [0.0 1.0] |> vec
        C_L["α"]["data"] = [-0.2200	-0.2200
                           	0.2500	0.2500
                           	0.7300	0.7300
                           	0.8300	0.7800
                           	0.9200	0.7900
                           	1.0200	0.8100
                           	1.0800	0.8200
                           	1.1300	0.8300
                           	1.1900	0.8500
                           	1.2500	0.8600
                           	1.3500	0.8800
                           	1.4400	0.9000
                           	1.4700	0.9200
                           	1.4300	0.9500
                           	1.3800	0.9900
                           	1.3000	1.0500
                           	1.1500	1.1500
        ]

        create_group(C_L, "δf")
        C_L["δf"]["δf"] = deg2rad.([ 0.0000	10.0000 20.0000 30.0000 ]) |> vec
        C_L["δf"]["data"] = [ 0.0000 0.2 0.3 0.35] |> vec


        ############################### C_l #################################

        create_group(fid, "C_l")
        C_l = fid["C_l"]

        C_l["δa"] = 0.229
        C_l["δr"] = 0.0147
        C_l["β"] = -0.09226
        C_l["p"] = -0.4840

        create_group(C_l, "r")
        C_l["r"]["α"] = [0.0 0.094] |> vec
        C_l["r"]["δf"] = deg2rad.([0 30]) |> vec
        C_l["r"]["data"] = [
                            0.0798	0.1246
                            0.1869	0.2317
        ]

        ############################# C_m ###################################

        create_group(fid, "C_m")
        C_m = fid["C_m"]

        C_m["zero"] = 0.100
        C_m["δe"] = -1.1220
        C_m["α"] = -1.8000
        C_m["q"] = -12.400
        C_m["α_dot"] = -7.2700

        create_group(C_m, "δf")
        C_m["δf"]["δf"] = deg2rad.([0 10 20 30]) |> vec
        C_m["δf"]["data"] = [ 0.0000 -0.0654 -0.0981 -0.1140 ] |> vec

        ############################# C_n ###################################

        create_group(fid, "C_n")
        C_n = fid["C_n"]

        C_n["δr"] = -0.0430
        C_n["δa"] = -0.0053
        C_n["β"] = 0.05874
        C_n["p"] = -0.0278
        C_n["r"] = -0.0937

    end

end

function load_aero_data( fname = "src/aircraft/c172r/aero_dataset.h5")

    fid = h5open(fname, "r")

    gr_C_D = fid["C_D"]
    gr_C_Y = fid["C_Y"]
    gr_C_L = fid["C_L"]
    gr_C_l = fid["C_l"]
    gr_C_m = fid["C_m"]
    gr_C_n = fid["C_n"]

    C_D = (
        z = gr_C_D["zero"] |> read,
        β = LinearInterpolation(gr_C_D["β"]["β"] |> read, gr_C_D["β"]["data"] |> read, extrapolation_bc = Flat()),
        δe = LinearInterpolation(gr_C_D["δe"]["δe"] |> read, gr_C_D["δe"]["data"] |> read, extrapolation_bc = Flat()),
        δf = LinearInterpolation(gr_C_D["δf"]["δf"] |> read, gr_C_D["δf"]["data"] |> read, extrapolation_bc = Flat()),
        α_δf = LinearInterpolation((gr_C_D["α_δf"]["α"] |> read,  gr_C_D["α_δf"]["δf"] |> read), gr_C_D["α_δf"]["data"] |> read, extrapolation_bc = Flat()),
        ge = LinearInterpolation(gr_C_D["ge"]["Δh_nd"] |> read, gr_C_D["ge"]["data"] |> read, extrapolation_bc = Flat())
    )

    C_Y = (
        δr = gr_C_Y["δr"] |> read,
        δa = gr_C_Y["δa"] |> read,
        β_δf = LinearInterpolation((gr_C_Y["β_δf"]["β"] |> read,  gr_C_Y["β_δf"]["δf"] |> read), gr_C_Y["β_δf"]["data"] |> read, extrapolation_bc = Flat()),
        p = LinearInterpolation((gr_C_Y["p"]["α"] |> read,  gr_C_Y["p"]["δf"] |> read), gr_C_Y["p"]["data"] |> read, extrapolation_bc = Flat()),
        r = LinearInterpolation((gr_C_Y["r"]["α"] |> read,  gr_C_Y["r"]["δf"] |> read), gr_C_Y["r"]["data"] |> read, extrapolation_bc = Flat()),
    )

    C_L = (
        δe = gr_C_L["δe"] |> read,
        q = gr_C_L["q"] |> read,
        α_dot = gr_C_L["α_dot"] |> read,
        α = LinearInterpolation((gr_C_L["α"]["α"] |> read,  gr_C_L["α"]["stall"] |> read), gr_C_L["α"]["data"] |> read, extrapolation_bc = Flat()),
        δf = LinearInterpolation(gr_C_L["δf"]["δf"] |> read, gr_C_L["δf"]["data"] |> read, extrapolation_bc = Flat()),
        ge = LinearInterpolation(gr_C_L["ge"]["Δh_nd"] |> read, gr_C_L["ge"]["data"] |> read, extrapolation_bc = Flat())
    )

    C_l = (
        δa = gr_C_l["δa"] |> read,
        δr = gr_C_l["δr"] |> read,
        β = gr_C_l["β"] |> read,
        p = gr_C_l["p"] |> read,
        r = LinearInterpolation((gr_C_l["r"]["α"] |> read,  gr_C_l["r"]["δf"] |> read), gr_C_l["r"]["data"] |> read, extrapolation_bc = Flat()),
    )

    C_m = (
        z = gr_C_m["zero"] |> read,
        δe = gr_C_m["δe"] |> read,
        α = gr_C_m["α"] |> read,
        q = gr_C_m["q"] |> read,
        α_dot = gr_C_m["α_dot"] |> read,
        δf = LinearInterpolation(gr_C_m["δf"]["δf"] |> read, gr_C_m["δf"]["data"] |> read, extrapolation_bc = Flat()),
    )

    C_n = (
        δr = gr_C_n["δr"] |> read,
        δa = gr_C_n["δa"] |> read,
        β = gr_C_n["β"] |> read,
        p = gr_C_n["p"] |> read,
        r = gr_C_n["r"] |> read,
    )

    close(fid)

    return (C_D = C_D, C_Y = C_Y, C_L = C_L, C_l = C_l, C_m = C_m, C_n = C_n)

end
