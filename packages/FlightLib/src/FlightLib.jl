module FlightLib

using Reexport

include("quaternions.jl"); @reexport using .Quaternions
include("attitude.jl"); @reexport using .Attitude
include("geodesy.jl"); @reexport using .Geodesy
include("kinematics.jl"); @reexport using .Kinematics
include("dynamics.jl"); @reexport using .Dynamics
include("atmosphere.jl"); @reexport using .Atmosphere
include("terrain.jl"); @reexport using .Terrain
include("linearization.jl"); @reexport using .Linearization
include("control.jl"); @reexport using .Control
include("propellers.jl"); @reexport using .Propellers
include("piston.jl"); @reexport using .Piston
include("landinggear.jl"); @reexport using .LandingGear
include("aircraftbase.jl"); @reexport using .AircraftBase
include("world.jl"); @reexport using .World

end