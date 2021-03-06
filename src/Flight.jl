module Flight

using Reexport
@reexport using BenchmarkTools

include("core/utils.jl")
include("core/systems.jl")
include("core/sim.jl")
include("core/input.jl")
include("core/output.jl")
include("core/plotting.jl")

include("physics/quaternions.jl")
include("physics/attitude.jl")
include("physics/geodesy.jl")
include("physics/kinematics.jl")
include("physics/air.jl")
include("physics/terrain.jl")
include("physics/dynamics.jl")
include("physics/friction.jl")

include("aircraft/aircraft.jl")
include("aircraft/landinggear.jl")
include("aircraft/propellers.jl")
include("aircraft/piston.jl")
include("aircraft/electrics.jl")

include("aircraft/c172r/c172r.jl")

@reexport using .Utils
@reexport using .Systems
@reexport using .Sim
@reexport using .Input
@reexport using .Output
@reexport using .Plotting

@reexport using .Quaternions
@reexport using .Attitude
@reexport using .Geodesy
@reexport using .Kinematics
@reexport using .Air
@reexport using .Terrain
@reexport using .Dynamics
@reexport using .Friction

@reexport using .Aircraft
@reexport using .Electrics
@reexport using .LandingGear
@reexport using .Propellers
@reexport using .Piston

@reexport using .C172R



end
