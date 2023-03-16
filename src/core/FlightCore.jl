module FlightCore

using Reexport

include("iodevices.jl"); @reexport using .IODevices
include("gui.jl"); @reexport using .GUI
include("utils.jl"); @reexport using .Utils
include("systems.jl"); @reexport using .Systems
include("sim.jl"); @reexport using .Sim
include("plotting.jl"); @reexport using .Plotting
include("joysticks.jl"); @reexport using .Joysticks
include("xplane.jl"); @reexport using .XPlane

end