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

const XLdg = LBVNode{:XLdg}
LBV.descriptor(::Type{<:LBVNode{:XLdg}}) = (μReg = LBVLeaf{2},)
#register_LBVNode(:XLdg) |> eval
#this accesses the descriptor internally, and runs from there

#generated code
function Base.length(::Type{<:LBVNode{:XLdg}})
    #queries the descriptor with the UnionAll as input
    return sum(length.(values(descriptor(LBVNode{:XLdg}))))
end

# eval(register_type(XLdg_desc))

end

#################### end submodule Ldg ########################

#################### submodule Aircraft ########################

module Aircraft
using Flight.LBV
using ..Rbd #needed to access XRbd
using ..Ldg #needed to access XLdg
export XAircraft

const XAircraft = LBVNode{:XAircraft}
LBV.descriptor(::Type{<:LBVNode{:XAircraft}}) = (rbd = XRbd, ldg = XLdg, pwp = LBVLeaf{4})
function Base.length(::Type{<:LBVNode{:XAircraft}})
    #queries the descriptor with the UnionAll as input
    return sum(length.(values(descriptor(LBVNode{:XAircraft}))))
end

# eval(register_type(XAircraft_desc))

end #submodule

#################### end submodule Aircraft ########################

using Flight.LBV
using Reexport
@reexport using .Rbd
@reexport using .Ldg
@reexport using .Aircraft

export test_lbv2

function test_lbv2()
end


end #module