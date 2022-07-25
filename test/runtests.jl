using Flight
using Test

include("test_quaternions.jl")
include("test_attitude.jl")
include("test_geodesy.jl")
include("test_kinematics.jl")
include("test_propellers.jl")
include("test_friction.jl")
include("test_piston.jl")
include("test_landing_gear.jl")
include("test_integration.jl")

using .TestQuaternions
using .TestAttitude
using .TestGeodesy
using .TestKinematics
using .TestPropellers
using .TestFriction
using .TestPiston
using .TestLandingGear
using .TestIntegration

test_quaternions()
test_attitude()
test_geodesy()
test_kinematics()
test_propellers()
test_friction()
test_landing_gear()
test_piston()
test_integration()