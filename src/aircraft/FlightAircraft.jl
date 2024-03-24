module FlightAircraft

using Reexport

include("aircraftbase.jl"); @reexport using .AircraftBase
include(normpath("c172/c172.jl")); @reexport using .C172

end