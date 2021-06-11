module TestLBV


########################### submodule Rbd ###################

module Rbd
using Flight.LBV
export XRbd

register_LBVNode(:XRbd, (:att, :vel, :pos), (LBVLeaf{4}, LBVLeaf{3}, LBVLeaf{4})) |> eval

end

#################### end submodule Rbd ########################

#################### submodule Ldg ########################

module Ldg
using Flight.LBV
export XLdg

register_LBVNode(:XLdg, (:μReg,), (LBVLeaf{2},)) |> eval

end

#################### end submodule Ldg ########################

#################### submodule Aircraft ########################

module Aircraft
using Flight.LBV
using ..Rbd #needed to access XRbd
using ..Ldg #needed to access XLdg
export XAircraft

register_LBVNode(:XAircraft, (:rbd, :ldg, :pwp), (XRbd, XLdg, LBVLeaf{4})) |> eval

end #submodule

#################### end submodule Aircraft ########################

using Flight.LBV
using Reexport
@reexport using .Rbd
@reexport using .Ldg
@reexport using .Aircraft

export test_lbv

function test_lbv()
    x = XAircraft()
    x_rbd = x.rbd
    x_rbd[4] = 0
    @show x
    v = XAircraft(view(x,:))
    v.pwp .= 3
    @show x
    y = similar(x)
    @. y = x + sin(x) + 2x
    @show y


end


end #module