module LabelledBlockVector

using BlockArrays, StaticArrays
import BlockArrays: axes, viewblock, getblock

#additions to export are not tracked by Revise!
export Node, Leaf, descriptor, generate_getindex_sym

#setting default values directly:
# https://mauro3.github.io/Parameters.jl/v0.9/manual.html

#customizar la representacion para que aparezcan los nombres de los child blocks

#since Leaf is only a dummy type for dispatching and computing lengths, every
#system must have a Node as its StateVector.
#if it has only one Leaf, the descriptor will be descriptor(Node{:simplesys}) =
#(singleleaf = Leaf{3})

# Leaf{K} = MVector{K, Float64} where {K}
struct Leaf{S} end
Base.length(::Type{Leaf{S}}) where {S} = S

# alternative definition, which enables different data types (vectors and views)
abstract type Node{D<:AbstractVector{Float64}} <: AbstractVector{Float64} end

Node(data::D) where {D<:AbstractVector{Float64}} = Node{D}(data)


#Array Interface ###############

function precompute_length(desc::NamedTuple)
    sum(length.(values(desc)))
end

function precompute_block_ranges(desc::NamedTuple)
    block_lengths = length.(values(desc))
    block_ranges = Vector{UnitRange{Int}}(undef, length(desc))
    offset = 0
    for (i,l) in enumerate(block_lengths)
        block_ranges[i] = (1 + offset) : (l + offset)
        offset += l
    end
    return NamedTuple{keys(desc)}(Tuple(block_ranges))
end

# function precompute_blockranges(::Type{T}) where {D, T <: Node{D}}
#     # println("Called blockranges with $T")
#     lengths = blocklengths(T)
#     blockranges = Vector{UnitRange{Int}}(undef, length(lengths))
#     offset = 0
#     for (i,l) in enumerate(lengths)
#         blockranges[i] = (1 + offset) : (l + offset)
#         offset += l
#     end
#     return blockranges
# end

# function precompute_blocklengths(desc::NamedTuple)
#     length.(values(desc))
# end

#GENERATE AT PARSE TIME for each Node subtype. replace with a constant output
#for each subtype.
#ASSERT THE DATA PASSED TO THE CONSTRUCTOR IS OF THE APPROPRIATE LENGTH
# function Base.length(::Type{T}) where {D, T<:Node{D}}
#     println("Called length with $T")
#     sum(blocklengths(T))
# end
# Base.size(x::Node) = size(getfield(x,:data))

#probably unnecessary
# function blockoffsets(::Type{Node{S}}) where {S}
#     offsets = Vector{Int}([])
#     current_offset = 1
#     for l in blocklengths(Node{S})
#         append!(offsets, current_offset)
#         current_offset += l
#     end
#     return offsets
# end



#aqui devolvere una llamada a Node{blocktype}(view(getfield(x,:data)),
#blockrange), donde blocktype y blockrange los obtengo del descriptor,
#recursivamente si hace falta! si !blocktype <: Node, entonces es Leaf. y
#entonces devuelvo una view sin mas. probar. pensar si blockranges debe devolver
#directamente un NamedTuple para no tener que andar buscando aqui el blocknumber
#en realidad, esta funcion debe devolver un Vector{Expr}, cada uno
#correspondiente a un getindex. despues, hago for ex in v eval(ex) end
# function generate_getindex_sym(::Type{Node{L}}, s::Symbol) where L
#     println("Getting descriptor for Node{:$L}")
#     d = descriptor(Node{L})
#     block_length = length(d[s])
#     println("Generating getindex for type Node{$L}, Val($s)")
#     println("Replace this with a true")
#     type_par = QuoteNode(L)
#     sym = QuoteNode(s)
#     return :(Base.getindex(x::Node{$type_par}, ::Val{$sym}) = getindex(getfield(x,:data), 1:$block_length))
# end

# https://stackoverflow.com/questions/48675081/is-it-possible-to-implement-a-type-factory-in-julia-without-using-eval
# https://stackoverflow.com/questions/39385408/julia-macros-with-multiple-return-expressions
#anadir despues getproperty(x, s) = getindex(x, Val(s))

#AbstractArray interface
#add Base.@propagate_inbounds!!!!!!!
Base.@propagate_inbounds Base.getindex(x::Node, i) = getindex(getfield(x,:data), i)
Base.@propagate_inbounds Base.setindex!(x::Node, v, i) = setindex!(getfield(x,:data), v, i)

#it is much faster to perform basic operations on the underlying PBV than
#broadcasting them. Broadcasting should be used only as a fallback for generic
#functions
Base.:(+)(x1::T, x2::T) where {T<:Node} = (#=println("No broadcast");=# T(getfield(x1,:data) + getfield(x2,:data)))
Base.:(-)(x1::T, x2::T) where {T<:Node} = (#=println("No broadcast");=# T(getfield(x1,:data) - getfield(x2,:data)))
Base.:(*)(x::T, a::Real) where {T<:Node} = (#=println("No broadcast");=# T(a * getfield(x,:data)))
Base.:(*)(a::Real, x::T) where {T<:Node} = (#=println("No broadcast");=# x * a)
#=
Base.@propagate_inbounds Base.setindex!(x::Node, v, ::Val{s}) where {s} = (x[Val(s)] .= v)


Base.getproperty(x::Node, s::Symbol) = getindex(x, Val(s))
Base.setproperty!(x::Node, s::Symbol, v) = setindex!(x, v, Val(s))

=#

# struct NodeStyle{S} <: Broadcast.AbstractArrayStyle{1} end
# NodeStyle{S}(::Val{1}) where {S} = NodeStyle{S}()
# Base.BroadcastStyle(::Type{Node{S}}) where {S} = NodeStyle{S}()
# function Base.similar(bc::Broadcast.Broadcasted{NodeStyle{S}}, ::Type{ElType}) where {S, ElType}
#     # println("Called similar bc")
#     pbv = PseudoBlockVector{Float64}(undef, axes(bc))
#     Node{S}(pbv)
# end
# Base.similar(::Type{Node{S}}) where {S} = Node{S}(PseudoBlockVector{Float64}(undef, blocklengths(Node{S})))
# #
# Base.dataids(x::Node{S}) where {S} = Base.dataids(x.data)







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
