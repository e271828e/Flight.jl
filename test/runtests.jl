using Flight
using Revise

includet(normpath("core/test_gui.jl")); using .TestGUI

includet(normpath("lib/test_quaternions.jl")); using .TestQuaternions
includet(normpath("lib/test_attitude.jl")); using .TestAttitude
includet(normpath("lib/test_geodesy.jl")); using .TestGeodesy
includet(normpath("lib/test_kinematics.jl")); using .TestKinematics
includet(normpath("lib/test_dynamics.jl")); using .TestDynamics
includet(normpath("lib/test_linearization.jl")); using .TestLinearization
includet(normpath("lib/test_control.jl")); using .TestControl
includet(normpath("lib/test_propellers.jl")); using .TestPropellers
includet(normpath("lib/test_piston.jl")); using .TestPiston
includet(normpath("lib/test_landing_gear.jl")); using .TestLandingGear
includet(normpath("lib/test_aircraft_base.jl")); using .TestAircraftBase
includet(normpath("lib/test_world.jl")); using .TestWorld

includet(normpath("apps/c172/test_c172s.jl")); using .TestC172S
includet(normpath("apps/c172/test_c172x.jl")); using .TestC172X
includet(normpath("apps/c172/test_c172x1.jl")); using .TestC172Xv1
includet(normpath("apps/c172/test_c172x2.jl")); using .TestC172Xv2
includet(normpath("apps/robot2d/test_robot2d.jl")); using .TestRobot2D

test_gui()
test_quaternions()
test_attitude()
test_geodesy()
test_kinematics()
test_dynamics()
test_linearization()
test_control()
test_propellers()
test_piston()
test_landing_gear()
test_aircraft_base()
test_world()

test_c172s()
test_c172x()
test_c172x1()
test_c172x2()
test_robot2d()