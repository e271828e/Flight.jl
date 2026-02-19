using Revise
using Test

includet(normpath("test_gui.jl")); using .TestGUI

# @testset verbose = true "FlightCore" begin
    test_gui()
# end