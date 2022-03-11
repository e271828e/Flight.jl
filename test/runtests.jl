using Flight
using Test

include("test_quaternions.jl")
include("test_attitude.jl")
include("test_geodesy.jl")
include("test_propellers.jl")
include("test_piston.jl")

using .TestQuaternions
using .TestAttitude
using .TestGeodesy
using .TestPropellers
using .TestPiston

test_quaternions()
test_attitude()
test_geodesy()
test_propellers()
test_piston()
