using Flight

include("physics/test_quaternions.jl"); using .TestQuaternions
include("physics/test_attitude.jl"); using .TestAttitude
include("physics/test_geodesy.jl"); using .TestGeodesy
include("physics/test_kinematics.jl"); using .TestKinematics

include("components/test_general.jl"); using .TestGeneral
include("components/test_propellers.jl"); using .TestPropellers
include("components/test_landing_gear.jl"); using .TestLandingGear
include("components/test_piston.jl"); using .TestPiston

include("aircraft/test_c172r.jl"); using .TestC172R

test_quaternions()
test_attitude()
test_geodesy()
test_kinematics()

test_general()
test_propellers()
test_landing_gear()
test_piston()

test_c172r()