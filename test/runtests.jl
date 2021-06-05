using Flight
using Test

include("test_quaternions.jl")
include("test_attitude.jl")
include("test_wgs84.jl")
#using TestQuaternions tries to load package TestQuaternions
using .TestQuaternions
using .TestAttitude
using .TestWGS84

# test_quaternions()
# test_attitude()
test_wgs84()
