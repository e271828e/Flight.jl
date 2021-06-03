using Flight
using Test

include("test_quaternions.jl")
include("test_attitude.jl")
#using TestQuaternions tries to load package TestQuaternions
using .TestQuaternions
using .TestAttitude

# test_quaternions()
test_attitude()
