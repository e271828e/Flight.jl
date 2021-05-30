module Flight

include("quaternions.jl")
include("attitude.jl")
# include("sketches/statevector.jl")
# include("sketches/statevectorpbv.jl")

using Reexport
@reexport using .Quaternions
@reexport using .Attitude
# @reexport using .StateVector
# @reexport using .StateVectorPBV

export ftest

function ftest()
    println("Welcome to Flight")
end
# Write your package code here.

end
