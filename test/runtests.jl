using Flight
using Test

include("test_quaternions.jl")
include("test_attitude.jl")
include("test_wgs84.jl")
include("test_lbv.jl")

using .TestQuaternions
using .TestAttitude
using .TestWGS84
using .TestLBV

# test_quaternions()
# test_attitude()
test_wgs84()
# test_lbv()
