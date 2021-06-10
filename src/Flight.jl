module Flight

include("quaternions.jl")
include("attitude.jl")
include("wgs84.jl")
# include("lbv.jl")
include("lbv2.jl")
# include("../test/test_lbv.jl")

using Reexport
@reexport using .Quaternions
@reexport using .Attitude
@reexport using .WGS84
# @reexport using .LabelledBlockVector
@reexport using .LBV

export ftest

function ftest()
    println("Welcome to Flight")
end
# Write your package code here.

end
