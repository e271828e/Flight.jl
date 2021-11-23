module Flight

using Reexport
@reexport using BenchmarkTools

#general
include("utils.jl")
include("modeling.jl")
include("plotting.jl")

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
include("components.jl")
include("airdata.jl")
include("propulsion.jl")
include("landinggear.jl")
include("aerodynamics.jl")

# #aircraft
include("aircraft.jl")
include("c172.jl")

#output interfaces
include("output.jl")

@reexport using .Utils
@reexport using .Modeling
@reexport using .Plotting

@reexport using .Quaternions
@reexport using .Attitude

@reexport using .Geodesy
@reexport using .Terrain
@reexport using .Atmosphere

@reexport using .Kinematics
@reexport using .Dynamics

@reexport using .Components
@reexport using .Airdata
@reexport using .Propulsion
@reexport using .LandingGear
@reexport using .Aerodynamics

@reexport using .Aircraft
@reexport using .C172

@reexport using .Output

export ftest

ftest() = println("Welcome to Flight")

end
