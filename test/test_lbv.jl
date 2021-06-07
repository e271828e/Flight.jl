#when using "using", we need the LabelledBlockVector qualifier to extend. see:
#https://docs.julialang.org/en/v1/manual/modules/#using-and-import-with-specific-identifiers,-and-adding-methods
#maybe create a macro to make this a bit more convenient


#problem: to construct the LBV hierarchy and the associated methods, the LBV
#functions need to have access from the LBV module to the descriptor methods for
#each of the Node parametric types contained in the hierarchy. this can only be
#achieved (AFAIK) by extending the descriptor() method originally defined in the
#LBV module. otherwise, these descriptor methods will not be in scope of the LBV
#methods that require them.

# however, if each module using the LBV module defines a method
#descriptor(::Type{Node{T}}) for its own parameter T, then there is the
#potential for one module overwritting another module's descriptor if they
#happen to choose the same parameter T. this is not solved by making Node
#abstract and requiring each module to define its own concrete subtype, because
#similary there is the possibility of module defining different implementations
#of descriptor(::Type{MyNodeSubtype}), because both think that MyNodeSubtype is
#being defined by nobody else.

#OPTION 1: create and export a macro from LBV that generates code function that
#does everything locally. for example, in module Aircraft I would do:
# @LBV :aircraft (rbd = Node{:rbd}, ldg = Node{:ldg})
# #this generates the following code
# descriptor(::Type{Node{:aircraft}}) = (rbd = Node{:rbd}, ldg = Node{:ldg})
# blocklengths(::Type{Node{:aircraft}}) = length.(values(descriptor(Node{:aircraft})))
# totallength(::Type{Node{:aircraft}}) = sum(blocklengths(Node{:aircraft}))
# function get_methods(::Type{Node{:aircraft}})
#     #getindex and setindex from block identifiers
# end
#now, for the above code to run, it needs access to
#descriptor(::Type{Node{:rbd}}) and descriptor(::Type{Node{:ldg}}). unless these
#are explictly requested by the Aircraft module by writing "using Rbd" and
#"using Ldg", they will not be available. this requires in turn that Rbd and Ldg
#export their respective descriptor() methods.

#OPTION 2: extend the descriptor method in LabelledBlockVector, but do so with a
#macro that checks if a method with the same signature already exists, and if
#so, raise an error. this is the simplest one, and the most secure. each module
#can then export its descriptor method

#OPTION 3: use a macro to define an alias const MyNodeType for Node{:T}. Then
#create a unique type parameter (mangled) T to avoid clashes. if some other
#module creates the same MyNodeType, hopefully the constants will clash. the
#problem is that behind the scenes, we have Node{:Tmangled}, so we will need to
#also define show methods that display MyNodeType. each module exports its own
#NodeType, so that I can do
# @LBV XAircraft (rbd = XRbd, ldg = XLdg)

#OPTION 4: define Node as an abstract type. use "global" to prevent multiple
#definitions of the same method

#OPTION 5: the descriptor should not be methods, but types!! that is what should
#be exported. what I need to do is pass the descriptor itself to the functions
#in LBV

#customizar la representacion para que aparezcan los nombres de los child blocks

#crear un type Descriptor que contiene un NTuple{Symbol} y un NTuple{DataType}

module TestLBV

module Rbd
using Flight.LabelledBlockVector

#supongamos que yo tengo un descriptor DescRbd

#Leaf es distinto cualitativamente de Node. Leaf no puede alojar data. Es solo
#un descriptor. Necesito un equivalente para Node. Basicamente, me defino un
#AbstractDescriptor, del que Leaf{S} es un subtype y NodeDescriptor{L} es otro.
#L sera el nombre que queremos que tenga el Node subtype. Dentro, tendra 2
#Tuples, uno de

#Node subtype must be defined from a macro as follows
@LBV NodeRbd (att = Leaf{4}, vel = Leaf{3}, pos = Leaf{4})

#
#bueno, en realidad para lo anterior tengo que haber elegido un nombre para el
#Node subtype. pero eso para mas adelante, cuando genere una macro. tendra una
#pinta como:

#primero, definir el descriptor, que NO tiene que ser un method, sino un
#NamedTuple de Symbols a Union{Leaf, Node} (probar como restringir esto, y si
#acaso definir un alias const NodeDescriptor). suponiendo que existen funciones
#Base.length() para todos los types que forman parte del descriptor, basta una
#llamada a precompute_lengths para obtener la longitud de antemano
const rbd_descriptor = (att = Leaf{4}, vel = Leaf{3}, pos = Leaf{4})
const rbd_length = LabelledBlockVector.precompute_length(rbd_descriptor)
const rbd_block_ranges = LabelledBlockVector.precompute_block_ranges(rbd_descriptor)
println(rbd_length)
println(rbd_block_ranges)

#ahora, con este valor numerico ya conocido, puedo definir directamente el Node
#subtype con un inner constructor que haga un assert de longitud al input data:
struct NodeRbd{D} <: Node{D}
    data::D
    function NodeRbd{D}(data::D) where {D}
        @assert length(data) == rbd_length "Expected an input of length $rbd_length"
        new{D}(data)
    end
end
#un outer constructor para comodidad
NodeRbd(data::D) where {D<:AbstractVector{Float64}} = NodeRbd{D}(data)

#before extending anything make sure this node subtype has not been defined already
if !hasmethod(Base.length, (Type{NodeRbd},))
    println("OK, NodeRbd has not been defined yet")
end

#this syntax takes care of all NodeRbd subtypes, which include
#NodeRbd{Vector{Float64}}, NodeRbd{SubArray{Float64,...}}, etc, but also the
#unqualified NodeRbd itself! the essential purpose of this method is to allow
#other Nodes to include this Node type in their own descriptors and precompute
#their lengths. make sure the let block is required, it probably isn't
Base.length(::Type{T}) where {T <: NodeRbd} = rbd_length

if hasmethod(Base.length, (Type{NodeRbd},))
    println("NodeRbd registered")
end
#this syntax does not accomodate NodeRbd without type parameters, only those
#parametric subtypes that are qualified with ANY parameter. but SOME parameter.
# function Base.length(::Type{NodeRbd{D}}) where {D}
#     let rbd_length = rbd_length
#         return rbd_length
#     end
# end

#y esto, que me lo pide el AbstractArray interface, y tambien es conocido de
#antemano.
Base.size(::NodeRbd) = (rbd_length,)
#aunque haria falta un Base.length(::NodeRbd{D} where {D}), pero eso me lo da
#indirectamente Base.size

@generated function Base.getindex(x::NodeRbd, ::Val{s}) where {s}
    #within the @generated function body, x is a type, but s is a Symbol, since
    #it is extracted from a type parameter.
    Core.println("Generated function getindex parsed for type $x, symbol $s")
    block_type = rbd_descriptor[s]
    block_range = rbd_block_ranges[s]
    if block_type <: Leaf
        return :(view(getfield(x,:data), $block_range))
    else #<: Node
        return :($block_type(view(getfield(x,:data), $block_range)))
    end
end

Base.setindex!(x::NodeRbd, v, ::Val{s}) where {s} = (x[Val(s)] .= v)

#despues definir los getindex[Val(s)] y setindex, getproperty y setproperty

#TERMINAR ESTO Y REPETIRLO PARA LOS OTROS
#SOLO CUANDO TENGA LOS TRES PROBADOS introducir las generating functions

#
# eval(generate_getindex_sym(NodeRbd, :att))
# eval(generate_getindex_sym(NodeRbd, :vel))
# eval(generate_getindex_sym(NodeRbd, :pos))

#pero claro, si defino Base.length para un NodeRbd{D}, eso me obliga a definirlo
#con un parametro D en cualquier descriptor que haga uso de el. aunque puesto
#que esto lo va a procesar una macro, podria prescindir de ello



end

module Ldg
using Flight.LabelledBlockVector

struct NodeLdg{D} <: Node{D}
    data::D
end


end

module Aircraft
using Flight.LabelledBlockVector
using ..Rbd #needed to access NodeRbd
using ..Ldg #needed to access

struct NodeAircraft{D} <: Node{D}
    data::D
end
#need to specify type parameters so that parameter. pero en
#realidad, puesto que
aircraft_desc = (rbd = Rbd.NodeRbd{Vector{Float64}}, ldg = Ldg.NodeLdg{Vector{Float64}})
# function LabelledBlockVector.descriptor(::Type{NodeAircraft{D}}) where {D}
#     (rbd = Rbd.NodeRbd{D}, ldg = Ldg.NodeLdg{D})
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
    # println(Rbd.NodeRbd <: Node{D} where {D})
    # LabelledBlockVector.blockranges(Rbd.NodeRbd{Vector{Float64}})
    # rbd_data = rand(length(Rbd.NodeRbd{Vector{Float64}}))
    # x_rbd = Rbd.NodeRbd(view(rbd_data, :))
    # y_rbd = Rbd.NodeRbd(rbd_data)
    # println(typeof(x_rbd))
    # println(typeof(y_rbd))

    # aircraft_data = rand(length(Aircraft.NodeAircraft{Vector{Float64}}))
    # x_aircraft = Aircraft.NodeAircraft(aircraft_data)
    # @show x_aircraft
    # # @show x
    # # y = Node{:aicraft}(view(rand(16), :))
    # # @show x[Val(:att)]
    # # y[Val(:)]
end

end #module