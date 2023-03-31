using Flight

include("physics/test_quaternions.jl"); using .TestQuaternions
include("physics/test_attitude.jl"); using .TestAttitude
include("physics/test_geodesy.jl"); using .TestGeodesy
include("physics/test_kinematics.jl"); using .TestKinematics

include("aircraft/test_control.jl"); using .TestControl
include("aircraft/test_propellers.jl"); using .TestPropellers
include("aircraft/test_landing_gear.jl"); using .TestLandingGear
include("aircraft/test_piston.jl"); using .TestPiston

include("aircraft/test_c172r.jl"); using .TestC172R
include("aircraft/test_world.jl"); using .TestWorld

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