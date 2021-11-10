using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using BenchmarkTools


#first run the test in isolation!

trn = HorizontalTerrain();
atm_sys = System(AtmosphereCmp());
ac = TestAircraft(
    kin = KinLTF(),
    mass = ConstantMass(),
    aero = SimpleDrag(),
    ldg = LandingGearLeg(steering = DirectSteering(), braking = DirectBraking()),
    pwp = AirframeGroup((
        left = EThruster(motor = ElectricMotor(α = CW)),
        right = EThruster(motor = ElectricMotor(α = CCW)))),
);
ac_sys = System(ac);
f_cont!(ac_sys, trn, atm_sys);
y_ac = ac_sys.y
ac_mdl = Model(ac_sys, (trn, atm_sys); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0.1);
b = @benchmarkable step!($ac_mdl, 1, true) setup=(reinit!($ac_mdl)); run(b)