using Flight
using Test

include("test_quaternions.jl")
#using TestQuaternions tries to load package TestQuaternions
using .TestQuaternions

@testset "Flight.jl" begin
    test_quaternions()
    # Write your tests here.
end