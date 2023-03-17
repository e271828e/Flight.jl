module FlightAircraft

using Reexport

include("control.jl"); @reexport using .Control
include("landinggear.jl"); @reexport using .LandingGear
include("propellers.jl"); @reexport using .Propellers
include("piston.jl"); @reexport using .Piston
include("template.jl"); @reexport using .Template
include("c172r/C172R.jl"); @reexport using .C172R

end