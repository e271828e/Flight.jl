using Flight

include(normpath("lib/test_quaternions.jl")); using .TestQuaternions
include(normpath("lib/test_attitude.jl")); using .TestAttitude
include(normpath("lib/test_geodesy.jl")); using .TestGeodesy
include(normpath("lib/test_kinematics.jl")); using .TestKinematics
include(normpath("lib/test_dynamics.jl")); using .TestDynamics
include(normpath("lib/test_control.jl")); using .TestControl
include(normpath("lib/test_propellers.jl")); using .TestPropellers
include(normpath("lib/test_piston.jl")); using .TestPiston
include(normpath("lib/test_landing_gear.jl")); using .TestLandingGear
include(normpath("lib/test_aircraft_base.jl")); using .TestAircraft

include(normpath("aircraft/c172/test_c172s.jl")); using .TestC172S

include(normpath("aircraft/c172/test_c172x.jl")); using .TestC172X
include(normpath("aircraft/c172/test_c172x1.jl")); using .TestC172Xv1

test_quaternions()
test_attitude()
test_geodesy()
test_kinematics()
test_dynamics()
test_control()
test_propellers()
test_piston()
test_landing_gear()
test_aircraft_base()

test_c172s()

test_c172x()
test_c172x1()