module Flight

using Reexport
@reexport using BenchmarkTools

#general
include("misc.jl")
include("plotting.jl")
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
include("dynamics.jl")

#aircraft components
include("airdata.jl")
include("propulsion.jl")
include("landinggear.jl")
include("propellers.jl")
include("piston.jl")
include("aerodynamics.jl")

# #aircraft
include("aircraft.jl")
include("c182t.jl")

@reexport using .Misc
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

@reexport using .Airdata
@reexport using .Propulsion
@reexport using .LandingGear
@reexport using .Propellers
@reexport using .Piston
@reexport using .Aerodynamics

@reexport using .Aircraft
@reexport using .C182T

export ftest

ftest() = println("Welcome to Flight")

end
