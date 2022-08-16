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
include("physics/rigidbody.jl")

include("components/common/common.jl")

include("components/environment/atmosphere.jl")
include("components/environment/terrain.jl")
include("components/environment/environment.jl")

include("components/aircraft/landinggear.jl")
include("components/aircraft/propellers.jl")
include("components/aircraft/piston.jl")
include("components/aircraft/aircraft.jl")
include("components/aircraft/c172r/c172r.jl")

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
@reexport using .RigidBody

@reexport using .Common

@reexport using .Atmosphere
@reexport using .Terrain
@reexport using .Environment

@reexport using .LandingGear
@reexport using .Propellers
@reexport using .Piston
@reexport using .Aircraft
@reexport using .C172R



end
