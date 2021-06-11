module TestLBV2


########################### submodule Rbd ###################

module Rbd
using Flight.LBV
export XRbd

# const XRbd = LBVNode{:XRbd} #UnionAll, equivalent to LBVBlock{:XRbd, D} where {D}

#this defines the method for the UnionAll LBVNode{:XRbd,D} where D, and all its
#specific subtypes such as LBVNode{:XRbd, }
# LBV.descriptor(::Type{<:LBVNode{:XRbd}}) = (att = LBVLeaf{4}, vel = LBVLeaf{3}, pos = LBVLeaf{4})

#generated code
# function Base.length(::Type{<:LBVNode{:XRbd}})
#     #queries the descriptor with the UnionAll as input
#     return sum(length.(values(descriptor(LBVNode{:XRbd}))))
# # end

register_node(:XRbd, (:att, :vel, :pos), (LBVLeaf{4}, LBVLeaf{3}, LBVLeaf{4})) |> eval

end

#################### end submodule Rbd ########################

#################### submodule Ldg ########################

module Ldg
using Flight.LBV
export XLdg

register_node(:XLdg, (:μReg,), (LBVLeaf{2},)) |> eval

end

#################### end submodule Ldg ########################

#################### submodule Aircraft ########################

module Aircraft
using Flight.LBV
using ..Rbd #needed to access XRbd
using ..Ldg #needed to access XLdg
export XAircraft

register_node(:XAircraft, (:rbd, :ldg, :pwp), (XRbd, XLdg, LBVLeaf{4})) |> eval

end #submodule

#################### end submodule Aircraft ########################

using Flight.LBV
using Reexport
@reexport using .Rbd
@reexport using .Ldg
@reexport using .Aircraft

export test_lbv2

function test_lbv2()
    x = XAircraft()
    x_rbd = x.rbd
    x_rbd[4] = 0
    @show x
    v = XAircraft(view(x,:))
    v.pwp .= 3
    @show x


end


end #module