module Flight

using Reexport
@reexport using BenchmarkTools

#general
include("utils.jl")
include("plotting.jl")

#math
include("quaternions.jl")
include("attitude.jl")

#systems and models
include("system.jl")
include("model.jl")

#environment & dynamics
include("geodesy.jl")
include("kinematics.jl")
include("terrain.jl")
include("atmosphere.jl")
include("dynamics.jl")

# # #airframe
include("airframe.jl")
include("airdata.jl")
include("propulsion.jl")
# # include("landinggear.jl")
include("statemachine.jl")

# # #aircraft
include("aircraft.jl")


@reexport using .Utils
@reexport using .Plotting

@reexport using .Quaternions
@reexport using .Attitude

@reexport using .System
@reexport using .Model

@reexport using .Geodesy
@reexport using .Kinematics
@reexport using .Terrain
@reexport using .Atmosphere
@reexport using .Dynamics

@reexport using .Airframe
@reexport using .Airdata
@reexport using .Propulsion
@reexport using .StateMachine
# # @reexport using .LandingGear

@reexport using .Aircraft

export ftest

ftest() = println("Welcome to Flight")

end
