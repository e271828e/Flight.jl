module Flight

using Reexport
@reexport using BenchmarkTools

#general
include("utils.jl")
include("modeling.jl")
include("plotting.jl")
include("input.jl")
include("output.jl")
include("simulation.jl")

#math
include("quaternions.jl")
include("attitude.jl")

#environment & dynamics
include("geodesy.jl")
include("terrain.jl")
include("atmosphere.jl")
include("kinematics.jl")
include("dynamics.jl")
include("friction.jl")

# #aircraft components
include("airdata.jl")
include("propulsion.jl")
include("landinggear.jl")
include("propellers.jl")
include("piston.jl")

# # #aircraft
include("aircraft.jl")
include("c172r.jl")

@reexport using .Utils
@reexport using .Modeling
@reexport using .Plotting
@reexport using .Simulation

@reexport using .Input
@reexport using .Output

@reexport using .Quaternions
@reexport using .Attitude

@reexport using .Geodesy
@reexport using .Terrain
@reexport using .Atmosphere
@reexport using .Kinematics
@reexport using .Dynamics
@reexport using .Friction

@reexport using .Airdata
@reexport using .Propulsion
@reexport using .LandingGear
@reexport using .Propellers
@reexport using .Piston

@reexport using .Aircraft
@reexport using .C172R

export ftest

ftest() = println("Welcome to Flight")



end
