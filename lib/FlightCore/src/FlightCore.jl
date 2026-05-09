module FlightCore

using Reexport
@reexport using About: about

include("modeling.jl"); @reexport using .Modeling
include("iodevices.jl"); @reexport using .IODevices
include("network.jl"); @reexport using .Network
include("joysticks.jl"); @reexport using .Joysticks
include("gui.jl"); @reexport using .GUI
include("sim.jl"); @reexport using .Sim
include("plotting.jl"); @reexport using .Plotting

end
