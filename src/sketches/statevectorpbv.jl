module StateVectorPBV

using BlockArrays, StaticArrays
import BlockArrays: axes, viewblock, getblock

#additions to export are not tracked by Revise!
export Node, Leaf

#setting default values directly:
# https://mauro3.github.io/Parameters.jl/v0.9/manual.html

#customizar la representacion para que aparezcan los nombres de los child blocks

#a pesar de todo esto, sigue siendo 20 veces mas lento que el PseudoBlockArray
#desnudo. no entiendo por que. en cualquier caso... importa tanto? es decir, si
#incluso aunque opere con mis arrays ineficientes en la integracion numerica, la
#mayor parte del tiempo de calculo se va a ir internamente en la evaluacion de
#las f's. que van a devolver x_dots.

# Leaf{K} = MVector{K, Float64} where {K}
struct Leaf{S} end
Base.length(::Type{Leaf{S}}) where {S} = S

struct Node{S} <: AbstractBlockVector{Float64}
    data::PseudoBlockVector{Float64}
    Node{S}(data::PseudoBlockVector{Float64}) where {S} = (#=println("Inner const");=# new{S}(data))
end

function Node{S}(data::AbstractVector{Float64}) where {S}
    # println("Abstract Vector const")
    # println(collect(blocklengths(Node{S})))
    # @assert length(data) == length(Node{S})
    Node{S}(PseudoBlockVector(data, blocklengths(Node{S})))
end

Node{S}() where {S} = Node{S}(PseudoBlockVector{Float64}(undef, blocklengths(Node{S})))

#AbstractBlockArray interface
axes(x::Node) = axes(getfield(x,:data))
viewblock(x::Node, block) = viewblock(getfield(x, :data), block)

descriptor(::Type{Node}) = error("To be implemented by each Node type")



#Array Interface ###############

# blocklengths(::Type{NodeRbd}) = [4, 3, 3]
# Base.length(::Type{NodeRbd}) = 10

######################### MAKE GENERATED #############################
blocklengths(::Type{Node{S}}) where {S} = collect(length.(values(descriptor(Node{S}))))
Base.length(::Type{Node{S}}) where {S} = sum(blocklengths(Node{S}))

#AbstractArray interface
Base.@propagate_inbounds Base.getindex(x::Node, i::Integer)::Float64 = getindex(getfield(x,:data), i)
Base.@propagate_inbounds Base.getindex(x::Node, i::Colon)::Vector{Float64} = getindex(getfield(x,:data),  i)
Base.@propagate_inbounds Base.getindex(x::Node, i::AbstractUnitRange)::Vector{Float64} = getindex(getfield(x,:data),  i)

Base.@propagate_inbounds Base.setindex!(x::Node, v, i::Integer) = setindex!(getfield(x,:data), v, i)
Base.@propagate_inbounds Base.setindex!(x::Node, v, i::Colon) = setindex!(getfield(x,:data), v, i)
Base.@propagate_inbounds Base.setindex!(x::Node, v, i::AbstractUnitRange) = setindex!(getfield(x,:data), v, i)

#it is much faster to perform basic operations on the underlying PBV than
#broadcasting them. Broadcasting should be used only as a fallback for generic
#functions
Base.:(+)(x1::Node{S}, x2::Node{S}) where {S} = (#=println("No broadcast");=# Node{S}(getfield(x1,:data) + getfield(x2,:data)))
Base.:(-)(x1::Node{S}, x2::Node{S}) where {S} = (#=println("No broadcast");=# Node{S}(getfield(x1,:data) - getfield(x2,:data)))
Base.:(*)(x::Node{S}, a::Real) where {S} = (#=println("No broadcast");=# Node{S}(a * getfield(x,:data)))
Base.:(*)(a::Real, x::Node{S}) where {S} = (#=println("No broadcast");=# x * a)
#=
Base.@propagate_inbounds Base.setindex!(x::Node, v, ::Val{s}) where {s} = (x[Val(s)] .= v)

@generated function Base.getindex(x::Node, ::Val{s}) where {s}
#within the @generated function body, x is a type, but s is a Symbol (since it
#is extracted from a type parameter).
    Core.println("Generated function getindex parsed for type $x")
    blocknumber = findfirst(i->i==s, keys(descriptor(x)))
    blocktype = descriptor(x)[s]
    #getblock is called on the data field and thus dispatches to the BlockArrays
    #method
    if blocktype <: Node
        return :($blocktype(view(getfield(x,:data), Block($blocknumber)))) #enforce return type for stability
    else #Leaf
        return :(view(getfield(x,:data), Block($blocknumber)))
    end
end

Base.getproperty(x::Node, s::Symbol) = getindex(x, Val(s))
Base.setproperty!(x::Node, s::Symbol, v) = setindex!(x, v, Val(s))

=#

struct NodeStyle{S} <: Broadcast.AbstractArrayStyle{1} end
NodeStyle{S}(::Val{1}) where {S} = NodeStyle{S}()
Base.BroadcastStyle(::Type{Node{S}}) where {S} = NodeStyle{S}()
function Base.similar(bc::Broadcast.Broadcasted{NodeStyle{S}}, ::Type{ElType}) where {S, ElType}
    println("Called similar bc")
    pbv = PseudoBlockVector{Float64}(undef, axes(bc))
    Node{S}(pbv)
end
Base.similar(::Type{Node{S}}) where {S} = Node{S}(PseudoBlockVector{Float64}(undef, blocklengths(Node{S})))
#
Base.dataids(x::Node{S}) where {S} = Base.dataids(x.data)



const NodeRbd = Node{:rbd}
descriptor(::Type{NodeRbd}) = (att = Leaf{4}, vel = Leaf{3}, pos = Leaf{3})




end


# #CONVERT TO @GENERATED
# function allocate_x!(sys::NodeSystem{S}, x::AbstractVector)
#     children_xlengths = [xlength(child) for child in sys.children]
#     #@assert lengths are correct
#     sys.x = PseudoBlockArray(x, children_xlengths)
#     for (i, child) in enumerate(sys.children)
#         allocate_x!(child, view(sys.x, Block(i)))
#     end
# end

# function allocate_x!(sys::LeafSystem{S}, x::AbstractVector)
#     sys.x = x
# end

#en el momento de crear un sistema, el constructor de sys asigna a sys.x un
#Node{S} del tipo apropiado. si tiene child subsystems, hara
#assign_subsystem_states!(sys)

#ahora, en este caso, lo que necesito hacer es, si aircraft tiene como child
#subsystems rbd y ldg:
#function assign_subsystem_states!(aircraft)
#aircraft.rbd.x = aircraft.x[:rbd]
#assign_subsystem_states!(aircraft.rbd)
#aircraft.ldg.x = aircraft.x[:ldg]
#assign_subsystem_states!(aircraft.ldg)
#end

#es decir, como getitem devuelve ya block views, esto consiste simplemente en
#asignar al x de cada subsystem el block view correspondiente del parent
#subsystem. a continuacion, tengo que hacer lo mismo de forma recursiva con cada
#child

#el problema de hacer assign_subsystem_states en el constructor es que, al ser
#una funcion recursiva, habra multiples pasadas por un mismo subsistema. por
#ejemplo, att sera asignado cuando construya rbd. pero cuando construya
#aircraft, asignara a rbd, y este a su vez recursivamente a att. esto en la
#practica no es un problema, porque solo va a ocurrir al crear la estructura del
#sistema, no en simulacion. pero si es poco estetico, se puede evitar haciendo
#que esta asignacion ocurra no en el propio constructor, sino con la funcion
#initialize!(sys). y la primera vez que haga step tambien se puede llamar a
#initialize state

#al final, lo que tenia en Python no estaba tan mal!!!!!! se trata de replicarlo
#con el enfoque Julia, nada mas. pero el concepto de asignar views a diferentes
#subsistemas, usar pseudoblockarrays, usar @generated functions son las
#novedades, y lo que va a hacer que vaya rapido.





#problema: si defino el x_aircraft como PseudoBlockArray, para que al cambiar
#los x_rbd y x_ldg me cambien los valores del x_aircraft, necesito asignar
#x_rbd = view(x_aircraft, Block(1)), x_ldg = view(x_aircraft, Block(2)). pero
#claro, para eso necesito assign_blocks!((x_rbd, x_ldg), x_aircraft). yo
#internamente compruebo que x_rbd tiene el tamano y el tipo correcto
