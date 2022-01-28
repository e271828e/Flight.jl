using StaticArrays
using OrdinaryDiffEq
using SciMLBase
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

    steering_sys = System(DirectSteering())
    strut_sys = System(Strut())

    terrain = HorizontalTerrain()

    #set the initial 2D Location
    l2d = LatLon()

    #get terrain altitude
    h_trn = Terrain.get_terrain_data(terrain, l2d).altitude

    #wow = false
    kin_init = KinInit(
        Ob = Geographic(l2d, h_trn + 10),
        q_nb = REuler(0, 0, 0),
        v_eOb_b = [0,0,0]
    )
    kin_data = KinData(kin_init)

    f_cont!(strut_sys, steering_sys, terrain, kin_data)
    @btime f_cont!($strut_sys, $steering_sys, $terrain, $kin_data)

    #wow = true
    kin_init = KinInit(
        Ob = Geographic(l2d, h_trn + 0*cos(π/6)),
        q_nb = REuler(0, 0, 0*π/6),
        v_eOb_b = [0,0,0],
        ω_lb_b = [1,0,0]
    )
    kin_data = KinData(kin_init)

    f_cont!(strut_sys, steering_sys, terrain, kin_data)
    @btime f_cont!($strut_sys, $steering_sys, $terrain, $kin_data)

end

function test_friction()

    roll = LandingGear.Rolling()
    skid = LandingGear.Skidding(wet = LandingGear.StaticDynamic(0.2, 0.1))

    dry = Terrain.DryTarmac; wet = Terrain.WetTarmac; icy = Terrain.IcyTarmac
    LandingGear.get_μ(roll, dry)
    LandingGear.get_μ(roll, wet)

    LandingGear.get_μ(skid, dry)
    LandingGear.get_μ(skid, wet)
    LandingGear.get_μ(skid, icy)

    LandingGear.get_μ(roll, dry, 0)
    LandingGear.get_μ(roll, dry, 0.5)
    LandingGear.get_μ(roll, dry, 1)

    LandingGear.get_μ(skid, icy, 0.1)

end

function test_ldg_unit()

    ldg_sys = System(LandingGearUnit())

    terrain = HorizontalTerrain()

    #set the initial 2D Location
    l2d = LatLon()

    #get terrain altitude
    h_trn = Terrain.get_terrain_data(terrain, l2d).altitude

    #wow = true
    kin_init = KinInit(
        Ob = Geographic(l2d, h_trn + 0.9),
        q_nb = REuler(0, 0, 0),
        v_eOb_b = [0,1,0]
    )
    kinematics = KinData(kin_init)

    @btime f_cont!($ldg_sys, $kinematics, $terrain)
    f_cont!(ldg_sys, kinematics, terrain)

end

#test on-aircraft integration
function test_ldg_ac()

    atm_sys = System(AtmosphereDescriptor());

    terrain = HorizontalTerrain()

    l2d = LatLon()
    h_trn = Terrain.get_terrain_data(terrain, l2d).altitude

    #wow = true
    kin_init = KinInit(
        Ob = Geographic(l2d, h_trn - 0.1),
        q_nb = REuler(0, 0, 0),
        v_eOb_b = [0,1,0]
    )

    ac_sys = System(AircraftBase());

    Aircraft.init!(ac_sys.x.kin, kin_init)
    ac_sys.x.kin.pos
    ac_sys.x.kin.vel

    f_cont!(ac_sys, terrain, atm_sys);
    @show ac_sys.ẋ.kin.vel
    println(ac_sys.y.ldg)
    @btime f_cont!($ac_sys, $terrain, $atm_sys)

    ac_mdl = Model(ac_sys, (terrain, atm_sys); dt = 0.01, adaptive = false, solver = Heun(), y_saveat = 0.1);
    b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)

    plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12)

    plots(ac_mdl; save_path = joinpath("tmp", "test_ldg_ac"), plot_settings...)

end