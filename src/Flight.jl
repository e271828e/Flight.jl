module Flight

using Reexport

export FlightCore, FlightLib, FlightExamples

include(normpath("core/FlightCore.jl")); @reexport using .FlightCore
include(normpath("lib/FlightLib.jl")); @reexport using .FlightLib
include(normpath("examples/FlightExamples.jl")); @reexport using .FlightExamples


end
