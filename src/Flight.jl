module Flight

using Reexport

include(normpath("core/FlightCore.jl")); @reexport using .FlightCore
include(normpath("physics/FlightPhysics.jl")); @reexport using .FlightPhysics
include(normpath("components/FlightComponents.jl")); @reexport using .FlightComponents
include(normpath("aircraft/FlightAircraft.jl")); @reexport using .FlightAircraft

#we only reexport some essential modules
@reexport using .FlightCore.Systems
@reexport using .FlightCore.Sim

@reexport using .FlightPhysics.Attitude
@reexport using .FlightPhysics.Kinematics
@reexport using .FlightPhysics.Geodesy
@reexport using .FlightPhysics.Environment: SimpleEnvironment

@reexport using .FlightComponents.Aircraft: trim!, linearize!
@reexport using .FlightComponents.World: SimpleWorld

end
