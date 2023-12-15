using Flight

include(normpath("physics/test_quaternions.jl")); using .TestQuaternions
include(normpath("physics/test_attitude.jl")); using .TestAttitude
include(normpath("physics/test_geodesy.jl")); using .TestGeodesy
include(normpath("physics/test_kinematics.jl")); using .TestKinematics

include(normpath("components/test_control.jl")); using .TestControl
include(normpath("components/test_propellers.jl")); using .TestPropellers
include(normpath("components/test_landing_gear.jl")); using .TestLandingGear
include(normpath("components/test_piston.jl")); using .TestPiston

include(normpath("aircraft/test_c172r.jl")); using .TestC172R
include(normpath("aircraft/test_c172r_base.jl")); using .TestC172RBase

include(normpath("aircraft/test_aircraft_base.jl")); using .TestAircraft
include(normpath("aircraft/test_c172fbw.jl")); using .TestC172FBW
include(normpath("aircraft/test_c172fbw_base.jl")); using .TestC172FBWBase
include(normpath("aircraft/test_c172cas.jl")); using .TestC172CAS
include(normpath("aircraft/test_c172mcs.jl")); using .TestC172MCS

test_quaternions()
test_attitude()
test_geodesy()
test_kinematics()

test_control()
test_propellers()
test_landing_gear()
test_piston()

test_aircraft_base()

test_c172r()
test_c172r_base()

test_c172fbw()
test_c172fbw_base()
test_c172cas()
test_c172mcs()