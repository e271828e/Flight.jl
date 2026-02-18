using Revise

includet(normpath("c172/test_c172s.jl")); using .TestC172S
includet(normpath("c172/test_c172x.jl")); using .TestC172X
includet(normpath("c172/test_c172x1.jl")); using .TestC172Xv1
includet(normpath("c172/test_c172x2.jl")); using .TestC172Xv2
includet(normpath("robot2d/test_robot2d.jl")); using .TestRobot2D

test_c172s()
test_c172x()
test_c172x1()
test_c172x2()
test_robot2d()