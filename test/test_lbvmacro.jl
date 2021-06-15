module TestLBV


########################### submodule Rbd ###################

module Rbd
using Flight.LBV
export XRbd

LBV.@lbv_macro XRbd LBVLeaf{4}

end #submodule

#################### submodule Ldg ########################

module Ldg
using Flight.LBV
export XLdg

LBV.@lbv_macro XLdg LBVLeaf{6}
# println(LBV.XRbd)

end #submodule

#################### submodule Aircraft ########################

module Aircraft
using Flight.LBV
using ..Rbd #needed to access XRbd
using ..Ldg #needed to access XLdg
export XAircraft

length(XRbd)
println(typeof(XRbd))
LBV.@lbv_macro XAircraft XRbd

end #submodule

#################### Tests ########################

using Flight.LBV
using Reexport
@reexport using .Rbd
@reexport using .Ldg
@reexport using .Aircraft

export test_lbvmacro

function test_lbvmacro()


end


end #module