module FlightAircraft

using Reexport

include("control.jl"); using .Control
include("landinggear.jl"); using .LandingGear
include("propellers.jl"); using .Propellers
include("piston.jl"); using .Piston
include("aircraft.jl"); using .Aircraft
include("world.jl"); using .World

include(normpath("c172r/c172r.jl")); @reexport using .C172R

end