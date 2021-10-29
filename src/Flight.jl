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
include("kinematics.jl")
include("terrain.jl")
include("atmosphere.jl")
include("dynamics.jl")

#airframe
include("airframe.jl")
include("airdata.jl")
include("propulsion.jl")
include("landinggear.jl")
include("aerodynamics.jl")
include("statemachine.jl")

#aircraft
include("aircraft.jl")


@reexport using .Utils
@reexport using .ModelingTools
@reexport using .Plotting

@reexport using .Quaternions
@reexport using .Attitude

@reexport using .Geodesy
@reexport using .Terrain
@reexport using .Atmosphere
@reexport using .Kinematics
@reexport using .Dynamics

@reexport using .Airframe
@reexport using .Airdata
@reexport using .Propulsion
@reexport using .Aerodynamics
@reexport using .StateMachine
@reexport using .LandingGear

@reexport using .Aircraft

export ftest

ftest() = println("Welcome to Flight")

end
