module C172R

using Reexport

include("airframe.jl"); using .Airframe
include(normpath("variants/c172rv0.jl")); @reexport using .C172Rv0
include(normpath("variants/c172rv1.jl")); @reexport using .C172Rv1
include(normpath("variants/c172rv2.jl")); @reexport using .C172Rv2

end