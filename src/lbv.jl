module LabelledBlockVector

using Base: Unsigned
using BlockArrays, StaticArrays
import BlockArrays: axes, viewblock, getblock

#additions to export are not tracked by Revise!
export Node, Leaf, Empty, LBlock#, descriptor
export block_type, block_length, child_ranges
export descriptor, is_registered, register_type, generate_constructors, generate_array_basics
# export @LBV

#ONCE THE CODE GENERATING FUNCTIONS ARE HERE, block_type, block_length,
#child_ranges will not have to be exported

#setting default values directly:
# https://mauro3.github.io/Parameters.jl/v0.9/manual.html


abstract type LBVDescriptor end

#convenience subtype to indicate a zero length LBV
struct Empty <: LBVDescriptor end
block_length(::Empty) = 0

struct Leaf{S} <: LBVDescriptor end
Leaf(S::Integer) = (@assert S>=0; Leaf{S}())
block_length(::Leaf{S}) where {S} =  S

struct Node{T} <: LBVDescriptor #T: type of the associated block, N: number of children
    children::NamedTuple
end
Node(label:: Symbol, nt) = Node{label}(nt)
block_length(n::Node) = sum(block_length.(values(n.children)))
block_type(::Node{T}) where {T} = T

function child_ranges(n::Node)
    child_lengths = block_length.(values(n.children))
    blk_ranges = Vector{Union{Nothing, UnitRange{Int}}}(undef, length(n.children))
    offset = 0
    for (i, l) in enumerate(child_lengths)
        blk_ranges[i] = (l > 0 ? ((1 + offset):(l + offset)) : nothing)
        offset += l
    end
    return NamedTuple{keys(n.children)}(Tuple(blk_ranges))
end

Base.getindex(n::Node, s::Symbol) = getindex(n.children, s)


#type parameter D enables different underlying data types (vectors and views)
abstract type LBlock{D<:AbstractVector{Float64}} <: AbstractVector{Float64} end

is_registered(::Type{<:Any}) = false
descriptor(::Type{<:Any}) = error("Must be implemented by concrete LBlock subtypes")

#create custom show methods. these can use the descriptor to know which blocks
#we have, their block ranges, etc!

#Array Interface ###############

# https://stackoverflow.com/questions/48675081/is-it-possible-to-implement-a-type-factory-in-julia-without-using-eval
# https://stackoverflow.com/questions/39385408/julia-macros-with-multiple-return-expressions


#Abstract types must not make any assumptions about data fields. either define
#getdata methods, or move these to each specific subtype

Base.@propagate_inbounds Base.getindex(::LBlock, ::Val{s}) where {s} = error("Must be defined by concrete subtype")
Base.@propagate_inbounds Base.setindex!(::LBlock, v, ::Val{s}) where {s} = error("Must be defined by concrete subtype")


# register_type(::Leaf) = :()

# function register_type(desc::Node)
#     btype = block_type(desc)
#     if is_registered(Val(btype))
#         println("WARNING: LBV type :$btype already registered")
#     end
#     btype_qn = QuoteNode(btype) #this is required if we want the quotes to show!
#     #DEFER THIS UNTIL ALL METHODS HAVE BEEN CREATED
#     ex = quote
#         LabelledBlockVector.is_registered(::Val{$btype_qn}) = true
#         println("Registered LBV subtype ", string($btype_qn))
#         # println("Registered LBV subtype ", string($btype))
#     end
#     return ex
# end


register_type(::Leaf) = :() #a leaf does not require
function register_type(desc::Node) #just a demo

    blength = block_length(desc)
    ranges = child_ranges(desc)
    btype = block_type(desc)
    btype_qn = QuoteNode(btype)

    blength == 0 && return :() #ignore zero-length blocks

    ex = Expr(:block) #equivalent to ex = quote end

    ex_const = quote

        println("Generating constructors for type $($btype_qn)... ")
        println("Block length: $($blength)")
        println("Block children: $(typeof($desc.children))")
        println("Child ranges: $(($ranges))")

        struct $btype{D} <: LBlock{D}
            data::D
            function $btype{D}(data::D) where {D}
                @assert length(data) == $blength "Expected an input of length $($blength)"
                new{D}(data)
            end
        end
        $btype(data::D) where {D<:AbstractVector{Float64}} = $btype{D}(data)
        $btype() = $btype(Vector{Float64}(undef, $blength))

        if is_registered($btype)
            println("\nWARNING: Type $($btype_qn) already registered")
        end

        LabelledBlockVector.is_registered(::Type{<:$btype}) = true
        LabelledBlockVector.descriptor(::Type{T} where {T<:$btype}) = Node{$btype_qn}($(desc.children))

    end

    push!(ex.args, ex_const)

    ex_array = quote

        println("Generating AbstractArray interface for type $($btype_qn)... ")
        # Base.length(::Type{T}) where {T <: $btype} = $blength

        Base.similar(::Type{<:$btype}) = $btype(Vector{Float64}(undef, $blength))
        Base.similar(::$btype) = similar($btype)
        Base.size(::$btype) = ($blength,)
        Base.getindex(x::$btype, i) = getindex(getfield(x,:data), i)
        Base.setindex!(x::$btype, v, i) = setindex!(getfield(x,:data), v, i)

    end

    push!(ex.args, ex_array)

    for (label, child, _range) in zip(keys(desc.children), desc.children, ranges)
        label_qn = QuoteNode(label)
        if _range === nothing #zero-length child
            ex3 = quote
                Base.getindex(x::$btype, ::Val{$label_qn}) = Vector{eltype(getfield(x,:data))}[]
            end
        elseif isa(child, Leaf)
            ex3 = quote
                Base.getindex(x::$btype, ::Val{$label_qn}) = view(getfield(x,:data), $_range)
            end
        elseif isa(child, Node)
            child_type = block_type(child)
            ex3 = quote
                Base.getindex(x::$btype, ::Val{$label_qn}) = $child_type(view(getfield(x,:data), $_range))
            end
        else
            error("Unexpected child type")
        end
        push!(ex.args, ex3)
    end

    return Base.remove_linenums!(ex)
    # return ex
end

# generate_array_basics(::Leaf) = :() #a leaf does not require
# function generate_array_basics(desc::Node)
#     btype = block_type(desc)
#     blength = block_length(desc)
#     btype_qn = QuoteNode(btype)
#     # println("Registering array basics for LBV type $btype")
#     ex = quote
#         # @assert !hasmethod(Base.size, ($btype,)) "Type "*string($btype)*" already registered"
#         #this syntax accomodates all XRbd subtypes, which include
#         #XRbd{Vector{Float64}}, XRbd{SubArray{Float64,...}}, etc, but also the
#         #unqualified XRbd itself! it is useful for similar(::Type{LBlock})
#         # Base.length(::Type{T}) where {T <: $btype} = $blength
#         Base.similar(::Type{<:$btype}) = $btype(Vector{Float64}(undef, $blength))
#         Base.similar(::$btype) = similar($btype)
#         Base.size(::$btype) = ($blength,)
#         Base.getindex(x::$btype, i) = getindex(getfield(x,:data), i)
#         Base.setindex!(x::$btype, v, i) = setindex!(getfield(x,:data), v, i)
#         # @assert hasmethod(Base.size, ($btype,)) "Failed to register type" * string($btype)
#         println("Generated array basics for type ", string($btype_qn))
#         println("Generated array basics for type ", string($btype))
#     end
#     return Base.remove_linenums!(ex)
#     # return ex
# end

#=
macro LBV(descriptor_symbol)
    # return :(generate_constructors)
    println(typeof(descriptor_symbol))
    println("$descriptor_symbol")
    println("$descriptor_symbol")
    #this inserts an expression in the source code to be evaluated later, after
    #macro parsing time. this expression is a call to a function that itself
    #returns an expression, and the evaluation of that expression

    #this returns an expression within an expression. when evaluated on macro
    #execution, the result is an expression! it will then be executed when the
    #module is compiled
    # ex = QuoteNode(:(eval(generate_constructors($descriptor_symbol))))
    ex = :(eval(generate_constructors($descriptor_symbol)))
    #if this works, i can also generate code to register the descriptor method
    #in LabelledBlockVector as before
    return ex

    # return :( :( eval(sin($($(esc(descriptor_symbol))))) ) )
end
=#

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
