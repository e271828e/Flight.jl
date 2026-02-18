using Revise

includet(normpath("test_quaternions.jl")); using .TestQuaternions
includet(normpath("test_attitude.jl")); using .TestAttitude
includet(normpath("test_geodesy.jl")); using .TestGeodesy
includet(normpath("test_kinematics.jl")); using .TestKinematics
includet(normpath("test_dynamics.jl")); using .TestDynamics
includet(normpath("test_linearization.jl")); using .TestLinearization
includet(normpath("test_control.jl")); using .TestControl
includet(normpath("test_propellers.jl")); using .TestPropellers
includet(normpath("test_piston.jl")); using .TestPiston
includet(normpath("test_landing_gear.jl")); using .TestLandingGear
includet(normpath("test_aircraft_base.jl")); using .TestAircraftBase
includet(normpath("test_world.jl")); using .TestWorld

test_quaternions()
test_attitude()
test_geodesy()
test_kinematics()
test_dynamics()
test_linearization()
test_control()
test_propellers()
test_piston()
test_landing_gear()
test_aircraft_base()
test_world()