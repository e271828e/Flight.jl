module Flight

using Reexport

export FlightCore, FlightLib, FlightApps

include(normpath("core/FlightCore.jl")); @reexport using .FlightCore
include(normpath("lib/FlightLib.jl")); @reexport using .FlightLib
include(normpath("apps/FlightApps.jl")); @reexport using .FlightApps


end
