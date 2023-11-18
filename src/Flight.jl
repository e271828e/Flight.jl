module Flight

using Reexport

include(normpath("core/FlightCore.jl")); @reexport using .FlightCore
include(normpath("physics/FlightPhysics.jl")); @reexport using .FlightPhysics
include(normpath("components/FlightComponents.jl")); @reexport using .FlightComponents
include(normpath("aircraft/FlightAircraft.jl")); @reexport using .FlightAircraft

#only reexport some essential, frequently used stuff
@reexport using .FlightCore.Systems
@reexport using .FlightCore.Sim

@reexport using .FlightPhysics.Attitude
@reexport using .FlightPhysics.Kinematics
@reexport using .FlightPhysics.Geodesy
@reexport using .FlightPhysics.Atmosphere

@reexport using .FlightComponents.Control
@reexport using .FlightComponents.Aircraft

end
