using Flight

include(normpath("physics/test_quaternions.jl")); using .TestQuaternions
include(normpath("physics/test_attitude.jl")); using .TestAttitude
include(normpath("physics/test_geodesy.jl")); using .TestGeodesy
include(normpath("physics/test_kinematics.jl")); using .TestKinematics

include(normpath("components/test_control.jl")); using .TestControl
include(normpath("components/test_propellers.jl")); using .TestPropellers
include(normpath("components/test_piston.jl")); using .TestPiston
include(normpath("components/test_landing_gear.jl")); using .TestLandingGear

include(normpath("aircraft/test_aircraft_base.jl")); using .TestAircraft

include(normpath("aircraft/test_c172r.jl")); using .TestC172R

include(normpath("aircraft/test_c172fbw.jl")); using .TestC172FBW
include(normpath("aircraft/test_c172fbw_v1.jl")); using .TestC172FBWv1

include(normpath("aircraft/test_c172rpa.jl")); using .TestC172RPA
include(normpath("aircraft/test_c172rpa_v1.jl")); using .TestC172RPAv1

# test_quaternions()
# test_attitude()
# test_geodesy()
# test_kinematics()

# test_control()
# test_propellers()
# test_piston()
# test_landing_gear()

test_aircraft_base()

test_c172r()

test_c172fbw()
test_c172fbw_v1()

test_c172rpa()
test_c172rpa_v1()