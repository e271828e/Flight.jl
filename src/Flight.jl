module Flight

include("quaternions.jl")
include("attitude.jl")
include("sketches/statevector.jl")

using Reexport
@reexport using .Quaternions
@reexport using .Attitude
@reexport using .StateVector

export ftest

function ftest()
    println("Welcome to Flight")
end
# Write your package code here.

end
