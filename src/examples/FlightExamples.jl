module FlightExamples

using Reexport

include(normpath("c172/c172.jl")); @reexport using .C172

end