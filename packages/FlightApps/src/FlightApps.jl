module FlightApps

using Reexport

include(normpath("c172/c172.jl")); @reexport using .C172
include(normpath("robot2d/robot2d.jl")); @reexport using .Robot2D

end