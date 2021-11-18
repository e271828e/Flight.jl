using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using BenchmarkTools


#first run the test in isolation!
function test_ldg_ac()

    atm_sys = System(AtmosphereCmp());

    terrain = HorizontalTerrain()

    #set initial 2D Location
    l2d = LatLon()

    #get terrain altitude
    h_trn = Terrain.get_terrain_data(terrain, l2d).altitude

    #wow = true
    kin_init = KinInit(
        Ob = Geographic(l2d, h_trn + 0.9),
        q_nb = REuler(0, 0, 0),
        v_eOb_b = [0,1,0]
    )

    ac = TestAircraft(
        kin = KinLTF(),
        mass = ConstantMass(),
        aero = SimpleDrag(),
        ldg = LandingGearUnit(steering = DirectSteering(), braking = DirectBraking()),
        pwp = AirframeGroup((
            left = EThruster(motor = ElectricMotor(α = CW)),
            right = EThruster(motor = ElectricMotor(α = CCW)))),
    );
    ac_sys = System(ac);

    init!(ac_sys.x.kin, kin_init)
    ac_sys.x.kin.pos
    ac_sys.x.kin.vel

    f_cont!(ac_sys, terrain, atm_sys);
    @show ac_sys.ẋ.kin.vel
    println(ac_sys.y.ldg)
    # @btime f_cont!($ac_sys, $trn, $atm_sys)
    # # y_ac = ac_sys.y
    # ac_mdl = Model(ac_sys, (trn, atm_sys); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0.1);
    # b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)
    # return ac_mdl

end