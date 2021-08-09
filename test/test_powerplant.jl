# module TestPowerplant

using Test
using DifferentialEquations
using Flight.Powerplant

export test_pwp

function test_pwp()

    sys = ElectricThrusterSystem()

    #integrate over Δt = 5, stopping exactly, and taking as many intermediate
    #steps as dictated by the adaptive algorithm
    step!(sys, 5, true)
    #test integration from a starting time to

end