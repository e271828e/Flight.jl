module Flight

using Reexport

#(sub) modules
# include("lbv.jl")
# @reexport using .LBV

include("quaternions.jl")
@reexport using .Quaternions

include("attitude.jl")
@reexport using .Attitude

include("wgs84.jl") #there are many, many constants here. should keep it a module
@reexport using .WGS84

include("system.jl")
@reexport using .System

include("airdata.jl")
@reexport using .AirData

include("kinematics.jl")
@reexport using .Kinematics

include("airframe.jl")
@reexport using .Airframe

include("powerplant.jl")
@reexport using .Powerplant



export ftest

function ftest()
    println("Welcome to Flight")
end
# Write your package code here.

end
