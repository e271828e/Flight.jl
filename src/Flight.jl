module Flight

include("quaternions.jl")
include("attitude.jl")

using Reexport
@reexport using .Quaternions
@reexport using .Attitude

export ftest

function ftest()
    println("Welcome to Flight")
end
# Write your package code here.

end
