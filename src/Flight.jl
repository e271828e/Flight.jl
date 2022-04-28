module Flight

using Reexport
@reexport using BenchmarkTools

#general
include("utils.jl")
include("systems.jl")
include("modeling.jl")
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
include("air.jl")
include("dynamics.jl")
include("friction.jl")

# #aircraft components
include("electrics.jl")
include("landinggear.jl")
include("propellers.jl")
include("piston.jl")

# # #aircraft
include("aircraft.jl")
include("c172r.jl")

include("plotting.jl")

@reexport using .Utils
@reexport using .Systems
@reexport using .Modeling
@reexport using .Simulation

@reexport using .Input
@reexport using .Output

@reexport using .Quaternions
@reexport using .Attitude

@reexport using .Geodesy
@reexport using .Terrain
@reexport using .Atmosphere
@reexport using .Kinematics
@reexport using .Air
@reexport using .Dynamics
@reexport using .Friction

@reexport using .Electrics
@reexport using .LandingGear
@reexport using .Propellers
@reexport using .Piston

@reexport using .Aircraft
@reexport using .C172R

@reexport using .Plotting


end
