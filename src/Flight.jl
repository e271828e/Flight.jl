module Flight

using Reexport
@reexport using BenchmarkTools

include("core/FlightCore.jl"); @reexport using .FlightCore
include("physics/FlightPhysics.jl"); @reexport using .FlightPhysics
include("components/FlightComponents.jl"); @reexport using .FlightComponents
include("aircraft/FlightAircraft.jl"); @reexport using .FlightAircraft

end
