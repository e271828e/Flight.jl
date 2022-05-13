module Flight

using Reexport
@reexport using BenchmarkTools

include("core/utils.jl")
include("core/systems.jl")
include("core/sim.jl")
include("core/input.jl")
include("core/output.jl")
include("core/plotting.jl")

include("common/quaternions.jl")
include("common/attitude.jl")
include("common/geodesy.jl")
include("common/kinematics.jl")
include("common/air.jl")
include("common/terrain.jl")
include("common/dynamics.jl")
include("common/friction.jl")

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

@reexport using .Electrics
@reexport using .LandingGear
@reexport using .Propellers
@reexport using .Piston

@reexport using .Aircraft
@reexport using .C172R



end
