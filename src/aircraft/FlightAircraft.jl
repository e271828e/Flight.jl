module FlightAircraft

using Reexport

include("control.jl"); @reexport using .Control
include("landinggear.jl"); @reexport using .LandingGear
include("propellers.jl"); @reexport using .Propellers
include("piston.jl"); @reexport using .Piston
include("aircraft.jl"); @reexport using .Aircraft
include("world.jl"); @reexport using .World

include(normpath("c172r/c172r.jl")); @reexport using .C172R

end