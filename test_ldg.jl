using StaticArrays
using BenchmarkTools
using Flight

function test_steering()

    no_steering_sys = System(NoSteering())
    set_steering_input(no_steering_sys, -1)
    f_cont!(no_steering_sys)
    get_steering_angle(no_steering_sys)

    direct_steering_sys = System(DirectSteering())
    set_steering_input(direct_steering_sys, -1)
    f_cont!(direct_steering_sys)
    try
        set_steering_input(direct_steering_sys, 2)
    catch AssertionError
        print("OK")
    end
    get_steering_angle(direct_steering_sys)

end

function test_braking()

    no_braking_sys = System(NoBraking())
    set_braking_input(no_braking_sys, 0)
    f_cont!(no_braking_sys)
    get_braking_coefficient(no_braking_sys)

    direct_braking_sys = System(DirectBraking())
    set_braking_input(direct_braking_sys, 0.1)
    f_cont!(direct_braking_sys)
    try
        set_braking_input(direct_braking_sys, -.5)
    catch AssertionError
        print("OK")
    end
    get_braking_coefficient(direct_braking_sys)

end

function test_damper()

    damper = LandingGear.SimpleDamper()
    LandingGear.get_damper_force(damper, -0.1, 0) > 0
    LandingGear.get_damper_force(damper, 0, -1) > 0
    LandingGear.get_damper_force(damper, damper.ξ_min - 0.1, 0) == 0

end

function test_strut()

    trn = HorizontalTerrain()

    strut_sys = System(Strut())

    #set the initial 2D Location
    l2d = LatLon()

    #get terrain altitude
    h_trn = Terrain.get_terrain_data(trn, l2d).altitude

    #wow = false
    kin_init = KinInit(
        Ob = Geographic(l2d, h_trn + 10), #initialize 1m above the terrain
        q_nb = REuler(0, 0, 0),
        v_eOb_b = [0,0,0]
    )
    kin_data = KinData(kin_init)

    f_cont!(strut_sys, kin_data, trn, 0.0)
    @btime f_cont!($strut_sys, $kin_data, $trn, 0.0)

    #wow = true
    kin_init = KinInit(
        Ob = Geographic(l2d, h_trn + 1 + 0*cos(π/6)),
        q_nb = REuler(0, 0, 0*π/6),
        v_eOb_b = [0,0,0],
        ω_lb_b = [1,0,0]
    )
    kin_data = KinData(kin_init)

    f_cont!(strut_sys, kin_data, trn, 0.0)
    @btime f_cont!($strut_sys, $kin_data, $trn, 0.0)

end

# function test_contact()

#     roll = LandingGear.Rolling()
#     skid = LandingGear.Skidding(wet = LandingGear.StaticDynamic(0.2, 0.1))

#     dry = Terrain.DryTarmac; wet = Terrain.WetTarmac; icy = Terrain.IcyTarmac
#     LandingGear.get_μ(roll, dry)
#     LandingGear.get_μ(roll, wet)

#     LandingGear.get_μ(skid, dry)
#     LandingGear.get_μ(skid, wet)
#     LandingGear.get_μ(skid, icy)

#     LandingGear.get_μ(roll, dry, 0)
#     LandingGear.get_μ(roll, dry, 0.5)
#     LandingGear.get_μ(roll, dry, 1)

#     LandingGear.get_μ(skid, icy, 0.1)

#     contact = Contact()
#     LandingGear.ContactY()

#     get_x0(contact)
#     get_y0(contact)

#     contact_sys = System(contact)

#     f_cont!(contact_sys, true, wet, 1, @SVector [0, 0, 0])
#     @btime f_cont!($contact_sys, true, $wet, 1, @SVector [0, 0, 0])
#     @btime f_cont!($contact_sys, false, $wet, 1, @SVector [0, 0, 0])

# end


# ldg_sys = System(ldg);

# f_cont!(ldg_sys, kin_data, trn)
# @btime f_cont!($ldg_sys, $kin_data, $trn)

# end