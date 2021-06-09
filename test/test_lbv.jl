module TestLBV


########################### submodule Rbd ###################

module Rbd
using Flight.LabelledBlockVector
export XRbd_desc, XRbd

#VERY IMPORTANT: here, we are extending the descriptor() function originally
#defined in the Flight.LabelledBlockVector module (by the way, this syntax works
#because Flight is a package). therefore, whenever we do "using
#Flight.LabelledBlockVector", this extension method for XRbd will be seen.
#however, if for instance we do "include("/src/lbv.jl")"; using
#.LabelledBlockVectors in the REPL, this is effectively a different module!!
#therefore, the method extensions previously added to Flight.LabelledBlockVector
#will not be part of it!!
#

const XRbd_desc = Node{:XRbd}((att = Leaf(4), vel = Leaf(3), pos = Leaf(4)))

const XRbd_block_length = block_length(XRbd_desc)
const XRbd_child_ranges = child_ranges(XRbd_desc)

eval(register_type(XRbd_desc))

#para no hacerlo demasiado complicado, quiza seria mejor construir estas
#directamente como expresiones, en vez de generated functions. porque un quote
#que contenga una generated function puede ser jodido. es meta-metaprogramming
# function generate_getindex()

#     return quote
#         @generated function Base.getindex(x::XRbXRbds}) where {s}
#             #within the @generated function body, x is a type, but s is a Symbol, since
#             #it is extracted from a type parameter.
#             Core.println("Generated function getindex parsed for type $x, symbol $s")
#             child = XRbd_desc[s]
#             brange = XRbd_child_ranges[s]
#             if isa(child, Leaf)
#                 return :(view(getfield(x,:data), $brange))
#             else #<: Node
#                 btype = block_type(child)
#                 return :($btype(view(getfield(x,:data), $brange)))
#             end
#         end
#     end

# end
# eval(generate_getindex())

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

#the notation x.att .= 4 calls getindex, but x.att = ones(4) calls setindex!, so
#we need both.
# @generated function Base.setindex!(x::XRbd, v, ::Val{s}) where {s}
#     Core.println("Generated function setindex! parsed for type $x, symbol $s")
#     brange = XRbd_child_ranges[s]
#     :(setindex!(getfield(x, :data), v, $brange))
# end

# struct XRbdStyle{D} <: Broadcast.AbstractArrayStyle{1} end
# XRbdStyle{D}(::Val{1}) where {D} = XRbdStyle{D}()
# Base.BroadcastStyle(::Type{XRbd{D}}) where {D} = XRbdStyle{D}()
# function Base.similar(::Broadcast.Broadcasted{XRbdStyle{D}}, ::Type{ElType}) where {D, ElType}
#     # println("Called similar for XRbd with type parameter $D")
#     similar(XRbd{D})
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

#customizar la representacion para que aparezcan los nombres de los child blocks

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