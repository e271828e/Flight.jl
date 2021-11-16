using Flight
using Test

include("test_quaternions.jl")
include("test_attitude.jl")
include("test_geodesy.jl")
# include("test_kinematics.jl")

using .TestQuaternions
using .TestAttitude
using .TestGeodesy
# using .TestKinematics

# test_quaternions()
# test_attitude()
test_geodesy()
# test_kinematics()
