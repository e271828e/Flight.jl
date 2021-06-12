module TestLBV


########################### submodule Rbd ###################

module Rbd
using Flight.LBV
export XRbd

register_LBVNode(:XRbd, (:att, :vel, :pos), (LBVLeaf{4}, LBVLeaf{3}, LBVLeaf{4})) |> eval

end #submodule

#################### submodule Ldg ########################

module Ldg
using Flight.LBV
export XLdg

register_LBVNode(:XLdg, (:μReg,), (LBVLeaf{2},)) |> eval

end #submodule

#################### submodule Aircraft ########################

module Aircraft
using Flight.LBV
using ..Rbd #needed to access XRbd
using ..Ldg #needed to access XLdg
export XAircraft

register_LBVNode(:XAircraft, (:rbd, :ldg, :pwp), (XRbd, XLdg, LBVLeaf{4})) |> eval

end #submodule

#################### Tests ########################

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
    #broadcasting with LBVLeafs of mixed parametric subtypes
    xv_rbd = x.pwp .+ v.pwp
    #broadcasting with LBVNodes of mixed parametric subtypes
    xv = x .+ v
    y = similar(x)
    @. y = x + sin(x) + 2x
    @show y

    # z1 = XAircraft{Complex}() #this yields undefs, which break operations,
    # need to initialize to zeros
    z1 = XAircraft(zeros(Complex, length(XAircraft)))
    z2 = XAircraft(ones(Int64, length(XAircraft)))
    z = exp.(z1) + z2
    @show z


end


end #module