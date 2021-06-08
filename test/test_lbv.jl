module TestLBV


########################### submodule Rbd ###################

module Rbd
using Flight.LabelledBlockVector
export XRbd

#VERY IMPORTANT: here, we are extending the descriptor() function originally
#defined in the Flight.LabelledBlockVector module (by the way, this syntax works
#because Flight is a package). therefore, whenever we do "using
#Flight.LabelledBlockVector", this extension method for XRbd will be seen.
#however, if for instance we do "include("/src/lbv.jl")"; using
#.LabelledBlockVectors in the REPL, this is effectively a different module!!
#therefore, the method extensions previously added to Flight.LabelledBlockVector
#will not be part of it!!
#

const XRbd_descriptor = Node{:XRbd}((att = Leaf(4), vel = Leaf(3), pos = Leaf(4)))
#otros ejemplos
# const XStatelessSystem_desc = Empty
# const XAnotherSystem_descriptor = Node{:XAnotherSystem}((a = Leaf{3}, b = XStatelessSystem_desc))

#this could be helped by these macros
# @LBV XRbd (att = Leaf(4), vel = Leaf(3), pos = Leaf(4))
# @LBV XStatelessSystem ()
# @LBV XStatelessSystem
# @LBV XAnotherSystem (a = Leaf(3), b = XStatelessSystem_desc)
#the macros call the code generation functions

#now all that follows is generated as an expression within a function. this
#function accepts simply a LBVDescriptor as its only input.

#these stay inside the generating function
const XRbd_block_length = block_length(XRbd_descriptor)
const XRbd_block_ranges = block_ranges(XRbd_descriptor)
println(XRbd_block_length)
println(XRbd_block_ranges)


#now we know the name of the LBlock subtype from the Node label, its length and
#the block_ranges
struct XRbd{D} <: LBlock{D}
    data::D
    function XRbd{D}(data::D) where {D}
        @assert length(data) == XRbd_block_length "Expected an input of length $XRbd_block_length"
        new{D}(data)
    end
end

XRbd(data::D) where {D<:AbstractVector{Float64}} = XRbd{D}(data)

#before extending anything make sure this node subtype has not been defined already
@assert !hasmethod(LabelledBlockVector.descriptor, (Type{XRbd},)) "XRbd subtype already defined"
#with "using", we need the LabelledBlockVector qualifier to extend. see:
LabelledBlockVector.descriptor(::Type{T}) where {T <: XRbd} = XRbd_descriptor
#this descriptor method extension allows us to simply export XRbd and access its
#descriptor from another module without explicitly exporting the descriptor as
#well. anyone who is using both Rbd and Flight.LabelledBlockArrays can get the
#XRbd_descriptor by simply calling descriptor(XRbd). this will not work if
#lbv.jl is included locally with then "using .LabelledBlockArrays"

@assert hasmethod(LabelledBlockVector.descriptor, (Type{XRbd},)) "XRbd failed"

#this syntax accomodates all XRbd subtypes, which include
#XRbd{Vector{Float64}}, XRbd{SubArray{Float64,...}}, etc, but also the
#unqualified XRbd itself! it is useful for similar(::Type{LBlock})
Base.length(::Type{T}) where {T <: XRbd} = XRbd_block_length
Base.similar(::Type{T}) where {T <: XRbd} = XRbd(Vector{Float64}(undef, XRbd_block_length))
Base.size(::XRbd) = (XRbd_block_length,)

#this syntax does not accomodate XRbd without type parameters, only those
#parametric subtypes that are qualified with ANY parameter. but SOME parameter.
# Base.length(::Type{XRbd{D}}) where {D} = XRbd_block_length

#para no hacerlo demasiado complicado, quiza seria mejor construir estas
#directamente como expresiones, en vez de generated functions. porque un quote
#que contenga una generated function puede ser jodido. es meta-metaprogramming
# function generate_getindex()

#     return quote
#         @generated function Base.getindex(x::XRbXRbds}) where {s}
#             #within the @generated function body, x is a type, but s is a Symbol, since
#             #it is extracted from a type parameter.
#             Core.println("Generated function getindex parsed for type $x, symbol $s")
#             child = XRbd_descriptor[s]
#             brange = XRbd_block_ranges[s]
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

@generated function Base.getindex(x::XRbd, ::Val{s}) where {s}
    #within the @generated function body, x is a type, but s is a Symbol, since
    #it is extracted from a type parameter.
    Core.println("Generated function getindex parsed for type $x, symbol $s")
    child = XRbd_descriptor[s]
    brange = XRbd_block_ranges[s]
    # error("Consider the case where brange is nothing, no method should be
    # generated in that case")
    if isa(child, Leaf)
        return :(view(getfield(x,:data), $brange))
    else #<: Node
        btype = block_type(child)
        return :($btype(view(getfield(x,:data), $brange)))
    end
end

#the notation x.att .= 4 calls getindex, but x.att = ones(4) calls setindex!, so
#we need both.
@generated function Base.setindex!(x::XRbd, v, ::Val{s}) where {s}
    Core.println("Generated function setindex! parsed for type $x, symbol $s")
    brange = XRbd_block_ranges[s]
    :(setindex!(getfield(x, :data), v, $brange))
end

struct XRbdStyle{D} <: Broadcast.AbstractArrayStyle{1} end
XRbdStyle{D}(::Val{1}) where {D} = XRbdStyle{D}()
Base.BroadcastStyle(::Type{XRbd{D}}) where {D} = XRbdStyle{D}()
function Base.similar(::Broadcast.Broadcasted{XRbdStyle{D}}, ::Type{ElType}) where {D, ElType}
    # println("Called similar for XRbd with type parameter $D")
    similar(XRbd{D})
end

# #it is much faster to perform basic operations on the underlying data than
# #broadcasting. Broadcasting should be used only as a fallback for generic
# #functions
Base.@propagate_inbounds Base.:(+)(x1::XRbd, x2::XRbd) = XRbd(getfield(x1,:data) + getfield(x2,:data))
Base.@propagate_inbounds Base.:(-)(x1::XRbd, x2::XRbd) = XRbd(getfield(x1,:data) + getfield(x2,:data))
Base.@propagate_inbounds Base.:(*)(x::XRbd, a::Real) = XRbd(a * getfield(x,:data))
Base.@propagate_inbounds Base.:(*)(a::Real, x::XRbd) = x * a

end

#################### end submodule Rbd ########################

#################### submodule Ldg ########################

module Ldg
using Flight.LabelledBlockVector

struct LBlockLdg{D} <: LBlock{D}
    data::D
end

end

#################### end submodule Ldg ########################

#################### submodule Aircraft ########################

module Aircraft
using Flight.LabelledBlockVector
using ..Rbd #needed to access XRbd
using ..Ldg #needed to access XLdg
export XAircraft

const XAircraft_descriptor = Node{:XAircraft}((rbd = descriptor(XRbd), ldg = Leaf(4), pwp = Empty))
if hasmethod(LabelledBlockVector.descriptor, (Type{XRbd},))
    println("Good, XRbd subtype already defined")
end

const XAircraft_block_length = block_length(XAircraft_descriptor)
const XAircraft_block_ranges = block_ranges(XAircraft_descriptor)
println(XAircraft_block_length)
println(XAircraft_block_ranges)


#customizar la representacion para que aparezcan los nombres de los child blocks

#now we know the name of the LBlock subtype from the Node label, its length and
#the block_ranges
struct XAircraft{D} <: LBlock{D}
    data::D
    function XAircraft{D}(data::D) where {D}
        @assert length(data) == XAircraft_block_length "Expected an input of length $XAircraft_block_length"
        new{D}(data)
    end
end

XAircraft(data::D) where {D<:AbstractVector{Float64}} = XAircraft{D}(data)

end #submodule

#################### end submodule Aircraft ########################

using Flight.LabelledBlockVector
using Reexport
@reexport using .Rbd
@reexport using .Aircraft

export test_lbv

function test_lbv()
    println(methods(descriptor))
    descriptor(XRbd)
    x = XRbd(rand(11))
    return x
    # println(Rbd.XRbd <: LBlock{D} where {D})
    # LabelledBlockVector.blockranges(Rbd.XRbd{Vector{Float64}})
    # XRbd_data = rand(length(Rbd.XRbd{Vector{Float64}}))
    # x_rbd = Rbd.XRbd(view(XRbd_data, :))
    # y_rbd = Rbd.XRbd(XRbd_data)
    # println(typeof(x_rbd))
    # println(typeof(y_rbd))

    # aircraft_data = rand(length(Aircraft.LBlockAircraft{Vector{Float64}}))
    # x_aircraft = Aircraft.LBlockAircraft(aircraft_data)
    # @show x_aircraft
    # # @show x
    # # y = LBlock{:aicraft}(view(rand(16), :))
    # # @show x[Val(:att)]
    # # y[Val(:)]
end



end #module