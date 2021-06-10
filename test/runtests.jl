using Flight
using Test

include("test_quaternions.jl")
include("test_attitude.jl")
include("test_wgs84.jl")
# include("test_lbv.jl")
include("test_lbv2.jl")

using .TestQuaternions
using .TestAttitude
using .TestWGS84
# using .TestLBV
using .TestLBV2

# test_quaternions()
# test_attitude()
# test_wgs84()
# test_lbv()
test_lbv2()
