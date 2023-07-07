module Flight

using Reexport

include(normpath("core/FlightCore.jl")); @reexport using .FlightCore
include(normpath("physics/FlightPhysics.jl")); @reexport using .FlightPhysics
include(normpath("aircraft/FlightAircraft.jl")); @reexport using .FlightAircraft

end
