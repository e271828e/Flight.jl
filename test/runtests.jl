using Flight

include(normpath("physics/test_quaternions.jl")); using .TestQuaternions
include(normpath("physics/test_attitude.jl")); using .TestAttitude
include(normpath("physics/test_geodesy.jl")); using .TestGeodesy
include(normpath("physics/test_kinematics.jl")); using .TestKinematics

include(normpath("aircraft/test_control.jl")); using .TestControl
include(normpath("aircraft/test_propellers.jl")); using .TestPropellers
include(normpath("aircraft/test_landing_gear.jl")); using .TestLandingGear
include(normpath("aircraft/test_piston.jl")); using .TestPiston
include(normpath("aircraft/test_c172r.jl")); using .TestC172R
include(normpath("aircraft/test_world.jl")); using .TestWorld

test_quaternions()
test_attitude()
test_geodesy()
test_kinematics()

test_control()
test_propellers()
test_landing_gear()
test_piston()
test_c172r()
test_world()