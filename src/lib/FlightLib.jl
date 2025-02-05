module FlightLib

using Reexport

include("quaternions.jl"); @reexport using .Quaternions
include("attitude.jl"); @reexport using .Attitude
include("geodesy.jl"); @reexport using .Geodesy
include("kinematics.jl"); @reexport using .Kinematics
include("dynamics.jl"); @reexport using .Dynamics
include("air.jl"); @reexport using .Air
include("terrain.jl"); @reexport using .Terrain

include("control.jl"); @reexport using .Control
include("propellers.jl"); @reexport using .Propellers
include("piston.jl"); @reexport using .Piston
include("landinggear.jl"); @reexport using .LandingGear
include("aircraftbase.jl"); @reexport using .AircraftBase

include("srukf.jl"); @reexport using .SRUKF

end