using Flight
using Test

include("test_quaternions.jl")
include("test_attitude.jl")
include("test_geodesy.jl")

using .TestQuaternions
using .TestAttitude
using .TestGeodesy

test_quaternions()
test_attitude()
test_geodesy()
