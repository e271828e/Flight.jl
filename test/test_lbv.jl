module TestLBV


########################### submodule Rbd ###################

module Rbd
using Flight.LabelledBlockVector
export XRbd_desc, XRbd

const XRbd_desc = Node{:XRbd}((att = Leaf(4), vel = Leaf(3), pos = Leaf(4)))

eval(register_type(XRbd_desc))

# @generated function Base.getindex(x::XRbd, ::Val{s}) where {s}
#     #within the @generated function body, x is a type, but s is a Symbol, since
#     #it is extracted from a type parameter.
#     Core.println("Generated function getindex parsed for type $x, symbol $s")
#     child = XRbd_desc[s]
#     brange = XRbd_child_ranges[s]
#     # error("Consider the case where brange is nothing, no method should be
#     # generated in that case")
#     if brange === nothing #zero-length child
#         return :(Vector{eltype(getfield(x,:data))}[])
#     end
#     if isa(child, Leaf)
#         return :(view(getfield(x,:data), $brange))
#     else #<: Node
#         btype = block_type(child)
#         return :($btype(view(getfield(x,:data), $brange)))
#     end
# end

end

#################### end submodule Rbd ########################

#################### submodule Ldg ########################

module Ldg
using Flight.LabelledBlockVector
export XLdg_desc

const XLdg_desc = Leaf(4)
eval(register_type(XLdg_desc))

end

#################### end submodule Ldg ########################

#################### submodule Aircraft ########################

module Aircraft
using Flight.LabelledBlockVector
using ..Rbd #needed to access XRbd
using ..Ldg #needed to access XLdg
export XAircraft, XAircraft_desc

const XAircraft_desc = Node{:XAircraft}((rbd = XRbd_desc, ldg = XLdg_desc))

eval(register_type(XAircraft_desc))

end #submodule

#################### end submodule Aircraft ########################

using Flight.LabelledBlockVector
using Reexport
@reexport using .Rbd
@reexport using .Aircraft

export test_lbv

function test_lbv()
    # println(methods(descriptor))
    x_rbd = XRbd(rand(11))
    @show x_rbd
    x_aircraft = XAircraft(rand(length(XAircraft)))
    @show x_aircraft
    x_rbd_view = XRbd(x_aircraft.rbd)
    @show x_rbd_view
    x_rbd_view .= 0
    @show x_aircraft
end



end #module