module FlightPhysics

using Reexport

include("quaternions.jl"); @reexport using .Quaternions
include("attitude.jl"); @reexport using .Attitude
include("geodesy.jl"); @reexport using .Geodesy
include("kinematics.jl"); @reexport using .Kinematics
include("dynamics.jl"); @reexport using .Dynamics
include("atmosphere.jl"); @reexport using .Atmosphere
include("terrain.jl"); @reexport using .Terrain

end