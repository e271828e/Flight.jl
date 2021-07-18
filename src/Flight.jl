module Flight

include("quaternions.jl")
include("attitude.jl")
include("wgs84.jl")
include("lbv.jl")
include("kinematics.jl")
include("airframe.jl")
# include("system.jl")

using Reexport
@reexport using .Quaternions
@reexport using .Attitude
@reexport using .WGS84
@reexport using .LBV
@reexport using .Kinematics
@reexport using .Airframe
# @reexport using .System

export ftest

function ftest()
    println("Welcome to Flight")
end
# Write your package code here.

end
