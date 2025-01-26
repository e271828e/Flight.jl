module Flight

using Reexport

export FlightCore, FlightLib, FlightAircraft

include(normpath("core/FlightCore.jl")); @reexport using .FlightCore
include(normpath("lib/FlightLib.jl")); @reexport using .FlightLib
include(normpath("aircraft/FlightAircraft.jl")); @reexport using .FlightAircraft


end
