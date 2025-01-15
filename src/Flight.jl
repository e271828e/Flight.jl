module Flight

using Reexport

export FlightCore, FlightPhysics, FlightComponents, FlightAircraft, FlightNavigation

include(normpath("core/FlightCore.jl")); @reexport using .FlightCore
include(normpath("physics/FlightPhysics.jl")); @reexport using .FlightPhysics
include(normpath("components/FlightComponents.jl")); @reexport using .FlightComponents
include(normpath("aircraft/FlightAircraft.jl")); @reexport using .FlightAircraft
include(normpath("navigation/FlightNavigation.jl")); @reexport using .FlightNavigation


end
