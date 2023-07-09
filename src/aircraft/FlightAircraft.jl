module FlightAircraft

using Reexport

include(normpath("c172r/c172r.jl")); @reexport using .C172R

end