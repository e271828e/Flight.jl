module Flight

using Reexport

#(sub) modules
# include("lbv.jl")
# @reexport using .LBV

include("quaternions.jl")
include("attitude.jl")
include("wgs84.jl") #many, many constants here. should keep it a module
include("system.jl")
include("model.jl")
include("airdata.jl")
include("kinematics.jl")
include("airframe.jl")
include("powerplant.jl")
include("aircraft.jl")

@reexport using .Quaternions
@reexport using .Attitude
@reexport using .WGS84
@reexport using .System
@reexport using .Model
@reexport using .Airdata
@reexport using .Kinematics
@reexport using .Airframe
@reexport using .Powerplant
@reexport using .Aircraft

println("REMINDER: Set normalization = false")

export ftest

ftest() = println("Welcome to Flight")

end
