using Revise
using Test

includet(normpath("test_gui.jl")); using .TestGUI
includet(normpath("test_sim.jl")); using .TestSim
includet(normpath("test_network.jl")); using .TestNetwork

test_gui()
test_sim()
test_network()
