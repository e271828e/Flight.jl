module Components

using Reexport

include("general/continuous.jl"); using .Continuous #avoid name clash with ControlSystems
include("general/discrete.jl"); using .Discrete
include("environment/atmosphere.jl"); @reexport using .Atmosphere
include("environment/terrain.jl"); @reexport using .Terrain
include("environment/environment.jl"); @reexport using .Environment
include("aircraft/landinggear.jl"); @reexport using .LandingGear
include("aircraft/propellers.jl"); @reexport using .Propellers
include("aircraft/piston.jl"); @reexport using .Piston

end