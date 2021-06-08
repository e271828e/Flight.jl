module LabelledBlockVector

using Base: Unsigned
using BlockArrays, StaticArrays
import BlockArrays: axes, viewblock, getblock

#additions to export are not tracked by Revise!
export Node, Leaf, Empty, LBlock, descriptor, block_type, block_length, block_ranges

#ONCE THE CODE GENERATING FUNCTIONS ARE HERE, block_type, block_length,
#block_ranges will not have to be exported

#setting default values directly:
# https://mauro3.github.io/Parameters.jl/v0.9/manual.html


abstract type LBVDescriptor end

struct Leaf{S} <: LBVDescriptor end
#a Leaf only holds length information. since labels are contained in each
#Descriptor's parent, and Leaf is parent to no other Descriptor, it does not
#need to hold a label
Leaf(S::Integer) = (@assert S>=0; Leaf{S}())
block_length(::Leaf{S}) where {S} =  S
const Empty = Leaf(0)

struct Node{T} <: LBVDescriptor #T: type of the associated block, N: number of children
    children::NamedTuple
end
Node(label:: Symbol, nt) = Node{label}(nt)

block_type(::Node{T}) where {T} = T

block_length(n::Node) = sum(block_length.(values(n.children)))

function block_ranges(n::Node)
    child_lengths = block_length.(values(n.children))
    blk_ranges = Vector{Union{Nothing, UnitRange{Int}}}(undef, length(n.children))
    offset = 0
    for (i, l) in enumerate(child_lengths)
        blk_ranges[i] = (l > 0 ? ((1 + offset):(l + offset)) : nothing)
        offset += 0
    end
    return NamedTuple{keys(n.children)}(Tuple(blk_ranges))
end

Base.getindex(n::Node, s::Symbol) = getindex(n.children, s)


#type parameter D enables different underlying data types (vectors and views)
abstract type LBlock{D<:AbstractVector{Float64}} <: AbstractVector{Float64} end

descriptor(::Type{T}) where {T<:LBlock} = error("Must be implemented by concrete subtypes")

#create custom show methods. these can use the descriptor to know which blocks
#we have, their block ranges, etc!

#Array Interface ###############

# https://stackoverflow.com/questions/48675081/is-it-possible-to-implement-a-type-factory-in-julia-without-using-eval
# https://stackoverflow.com/questions/39385408/julia-macros-with-multiple-return-expressions


#Abstract types must not make any assumptions about data fields. either define
#getdata methods, or move these to each specific subtype

Base.@propagate_inbounds Base.getindex(x::LBlock, i) = getindex(getfield(x,:data), i)
Base.@propagate_inbounds Base.setindex!(x::LBlock, v, i) = setindex!(getfield(x,:data), v, i)
Base.@propagate_inbounds Base.getindex(::LBlock, ::Val{s}) where {s} = error("Must be defined by concrete subtype")
Base.@propagate_inbounds Base.setindex!(::LBlock, v, ::Val{s}) where {s} = error("Must be defined by concrete subtype")
Base.@propagate_inbounds Base.getproperty(x::LBlock, s::Symbol) = getindex(x, Val(s))
Base.@propagate_inbounds Base.setproperty!(x::LBlock, s::Symbol, v) = setindex!(x, v, Val(s))

function generate_constructor(d::LBVDescriptor) #just a demo
    fname = :getindex
    lbv_type = block_type(d)
    println(lbv_type)
    ex = quote
        @generated function Base.getindex(x::$(lbv_type), ::Val{s}) where {s}
            #within the @generated function body, x is a type, but s is a Symbol, since
            #it is extracted from a type parameter.
            Core.println("Generated function $($fname) parsed for type $x, symbol $s")
            child = rbd_descriptor[s]
            brange = rbd_block_ranges[s]
            if isa(child, Leaf)
                return :(view(getfield(x,:data), $brange))
            else #<: Node
                btype = block_type(child)
                return :($btype(view(getfield(x,:data), $brange)))
            end
        end
    end
    return Base.remove_linenums!(ex)
end


end



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
