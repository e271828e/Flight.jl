#when using "using", we need the LabelledBlockVector qualifier to extend. see:
#https://docs.julialang.org/en/v1/manual/modules/#using-and-import-with-specific-identifiers,-and-adding-methods

#customizar la representacion para que aparezcan los nombres de los child blocks


module TestLBV

module Rbd
using Flight.LabelledBlockVector


#LBlock subtype must be defined from a macro as follows
# @LBV NodeRbd (att = Leaf{4}, vel = Leaf{3}, pos = Leaf{4})

#
#bueno, en realidad para lo anterior tengo que haber elegido un nombre para el
#LBlock subtype. pero eso para mas adelante, cuando genere una macro. tendra una
#pinta como:

#primero, definir el descriptor, que NO tiene que ser un method, sino un
#NamedTuple de Symbols a Union{Leaf, Node} (probar como restringir esto, y si
#acaso definir un alias const NodeDescriptor). suponiendo que existen funciones
#Base.length() para todos los types que forman parte del descriptor, basta una
#llamada a precompute_lengths para obtener la longitud de antemano
const rbd_descriptor = Node{:LBlockRbd}((att = Leaf(4), vel = Leaf(3), pos = Leaf(4)))

#now all that follows is generated as an expression within a function
const rbd_length = block_length(rbd_descriptor)
const rbd_block_ranges = block_ranges(rbd_descriptor)
# println(rbd_length)
# println(rbd_block_ranges)

#now we know the name of the LBlock subtype from the Node label, its length and
#the block_ranges
struct LBlockRbd{D} <: LBlock{D}
    data::D
    function LBlockRbd{D}(data::D) where {D}
        @assert length(data) == rbd_length "Expected an input of length $rbd_length"
        new{D}(data)
    end
end

LBlockRbd(data::D) where {D<:AbstractVector{Float64}} = LBlockRbd{D}(data)

#before extending anything make sure this node subtype has not been defined already
@assert !hasmethod(Base.size, (Type{LBlockRbd},)) "LBlockRbd subtype already defined"

#this syntax accomodates all LBlockRbd subtypes, which include
#LBlockRbd{Vector{Float64}}, LBlockRbd{SubArray{Float64,...}}, etc, but also the
#unqualified LBlockRbd itself! it is useful for similar(::Type{LBlock})
Base.length(::Type{T}) where {T <: LBlockRbd} = rbd_length
Base.similar(::Type{T}) where {T <: LBlockRbd} = LBlockRbd(Vector{Float64}(undef, rbd_length))
Base.size(::LBlockRbd) = (rbd_length,)

#this syntax does not accomodate LBlockRbd without type parameters, only those
#parametric subtypes that are qualified with ANY parameter. but SOME parameter.
# Base.length(::Type{LBlockRbd{D}}) where {D} = rbd_length

#para no hacerlo demasiado complicado, quiza seria mejor construir estas
#directamente como expresiones, en vez de generated functions. porque un quote
#que contenga una generated function puede ser jodido. es meta-metaprogramming
# function generate_getindex()

#     return quote
#         @generated function Base.getindex(x::LBlockRbd, ::Val{s}) where {s}
#             #within the @generated function body, x is a type, but s is a Symbol, since
#             #it is extracted from a type parameter.
#             Core.println("Generated function getindex parsed for type $x, symbol $s")
#             child = rbd_descriptor[s]
#             brange = rbd_block_ranges[s]
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
@generated function Base.getindex(x::LBlockRbd, ::Val{s}) where {s}
    #within the @generated function body, x is a type, but s is a Symbol, since
    #it is extracted from a type parameter.
    Core.println("Generated function getindex parsed for type $x, symbol $s")
    child = rbd_descriptor[s]
    brange = rbd_block_ranges[s]
    if isa(child, Leaf)
        return :(view(getfield(x,:data), $brange))
    else #<: Node
        btype = block_type(child)
        return :($btype(view(getfield(x,:data), $brange)))
    end
end

#the notation x.att .= 4 calls getindex, but x.att = ones(4) calls setindex!, so
#we need both
@generated function Base.setindex!(x::LBlockRbd, v, ::Val{s}) where {s}
    Core.println("Generated function setindex! parsed for type $x, symbol $s")
    brange = rbd_block_ranges[s]
    :(setindex!(getfield(x, :data), v, $brange))
end

struct LBlockRbdStyle{D} <: Broadcast.AbstractArrayStyle{1} end
LBlockRbdStyle{D}(::Val{1}) where {D} = LBlockRbdStyle{D}()
Base.BroadcastStyle(::Type{LBlockRbd{D}}) where {D} = LBlockRbdStyle{D}()
function Base.similar(::Broadcast.Broadcasted{LBlockRbdStyle{D}}, ::Type{ElType}) where {D, ElType}
    # println("Called similar for LBlockRbd with type parameter $D")
    similar(LBlockRbd{D})
end

# #it is much faster to perform basic operations on the underlying data than
# #broadcasting. Broadcasting should be used only as a fallback for generic
# #functions
Base.@propagate_inbounds Base.:(+)(x1::LBlockRbd, x2::LBlockRbd) = LBlockRbd(getfield(x1,:data) + getfield(x2,:data))
Base.@propagate_inbounds Base.:(-)(x1::LBlockRbd, x2::LBlockRbd) = LBlockRbd(getfield(x1,:data) + getfield(x2,:data))
Base.@propagate_inbounds Base.:(*)(x::LBlockRbd, a::Real) = LBlockRbd(a * getfield(x,:data))
Base.@propagate_inbounds Base.:(*)(a::Real, x::LBlockRbd) = x * a

end

#pero una vez definidos los LBVDescriptors, que exporta cada module? Ya no puede
#ser un tipo, tiene que ser un objeto. Una instancia de un Node concreto. O
#varias distintas.





module Ldg
using Flight.LabelledBlockVector

struct LBlockLdg{D} <: LBlock{D}
    data::D
end


end

module Aircraft
using Flight.LabelledBlockVector
using ..Rbd #needed to access LBlockRbd
using ..Ldg #needed to access

struct LBlockAircraft{D} <: LBlock{D}
    data::D
end
#need to specify type parameters so that parameter. pero en
#realidad, puesto que
aircraft_desc = (rbd = Rbd.LBlockRbd{Vector{Float64}}, ldg = Ldg.LBlockLdg{Vector{Float64}})
# function LabelledBlockVector.descriptor(::Type{LBlockAircraft{D}}) where {D}
#     (rbd = Rbd.LBlockRbd{D}, ldg = Ldg.LBlockLdg{D})
# end

end
#
#STEP 2:
#READ HANDS ON DESIGN PATTERNS and understand how to interpolate symbols into
#expressions. either a macro or eval could work. no, IT MUST BE EVAL!

using Flight.LabelledBlockVector
using .Rbd, .Ldg, .Aircraft
export test_lbv

function test_lbv()
    # println(Rbd.LBlockRbd <: LBlock{D} where {D})
    # LabelledBlockVector.blockranges(Rbd.LBlockRbd{Vector{Float64}})
    # rbd_data = rand(length(Rbd.LBlockRbd{Vector{Float64}}))
    # x_rbd = Rbd.LBlockRbd(view(rbd_data, :))
    # y_rbd = Rbd.LBlockRbd(rbd_data)
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