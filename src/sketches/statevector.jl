module StateVector

using BlockArrays, StaticArrays
import BlockArrays: axes, viewblock, getblock

#additions to export are not tracked by Revise!
export Node, lbv_demo

#setting default values directly:
# https://mauro3.github.io/Parameters.jl/v0.9/manual.html

#customizar la representacion para que aparezcan los nombres de los child blocks

#when operating with a complete Node, use [:] to output a Vector{Float64} view
#and ensure type stability.

Leaf{K} = MVector{K, Float64} where {K}

struct Node{S} <: AbstractBlockVector{Float64}
    data::BlockVector{Float64}
    Node{S}(data::BlockVector{Float64}) where {S} = new{S}(data)
end

#as an alternative to the runtime block type checks or conversions, use
#childtypes(Node{S}) to obtain the exact ordered block types for Node{S},
#and then use a macro to build the constructor with that specific input argument
#type. take the tuple of DataTypes, extract them as expressions/symbols, and put
#them inside a Tuple{}. simply doing typeof(makechildren(Node{S})) will not
#work, because it requires the constructor, which we are making right now!
#also, this constructor must be defined as an OUTER constructor, with an
#implicit constructor taking only data as input, but which should not be called
#directly. to prevent calling it, we could implement a single argument version
#that throws an error, and then another with a second dummy argument that is
#used by the type-correct outer constructor to call it.
#principle: if blocks are provided, they must be of the exact right type.

#the above approach implies that no single outer constructor exists for a
#generic Node. the constructor for Node{:specific} will be defined wherever the
#Node{:specific} type is defined by the macro. what happens when one tries to
#define such a type outside the scope of a module? nothing, the descriptor
#function and constructor are generated in that scope. but what about the
#generic type Node itself? will it be visible to the macro? yes if both Node and
#the macro are exported together by the StateVector module
function Node{S}(blocks::NTuple{N, Union{Node, Leaf}} where {N}) where {S}
    # println("Called block assembly constructor for $S")
    # println("Converting inputs of type $(typeof.(blocks)) to $properblocktypes")
    blocks = Tuple(convert(T, b) for (T, b) in zip(childtypes(Node{S}), blocks))
    # @assert typeof.(blocks) == childtypes(Node{S})
    data = mortar(collect(blocks))
    Node{S}(data)
end

descriptor(::Type{Node}) = error("To be implemented by each Node type")
childtypes(T::Type{Node{S}} where{S}) = values(descriptor(T))

#when the blocks are not supplied directly, we still need to be able to
#construct a Node{S} instance that respects the specific types at each node of
#its underlying hierarchy. this is done by having a similar() method that
#dispatches on type rather than on an instance, and calling it recursively from
#top to bottom. any Node hierarchy will be made of Node nodes and MVectors at
#the leafs. MVectors already have a similar() method dispatching on type. to
#achieve the above, we need another one for Node, which will implement the
#recursion
makechildren(::Type{Node{S}}) where {S} = (println("Creating children for Node{:$S}..."); similar.(childtypes(Node{S})))
Base.similar(::Type{Node{S}}) where {S} = Node{S}(makechildren(Node{S}))
Base.similar(::Node{S}) where {S} = similar(Node{S})

Node{S}(v::AbstractVector) where {S} = (x = similar(Node{S}) ; x .= v)
Node{S}() where {S} = (x = similar(Node{S}) ; x .= 0)

Base.convert(::Type{Node{S}}, x::Node{S}) where {S} = (println("Conversion unnecessary for $(typeof(x))"); x)
Base.convert(::Type{Node{S}}, v::AbstractVector) where {S} = (println("Converting $(typeof(v)) to Node{$S}"); Node{S}(v))
Base.copy(x::Node{S}) where {S} = Node{S}(Tuple(copy.(blocks(getfield(x,:data)))))

#AbstractBlockArray interface
axes(x::Node) = axes(getfield(x,:data))
viewblock(x::Node, block) = viewblock(getfield(x, :data), block)

#AbstractArray interface
Base.getindex(x::Node, i::Integer)::Float64 = getindex(getfield(x,:data), i)
Base.getindex(x::Node, i::Colon)::Vector{Float64} = getindex(getfield(x,:data),  i)
Base.getindex(x::Node, i::AbstractUnitRange)::Vector{Float64} = getindex(getfield(x,:data),  i)
@generated function Base.getindex(x::Node, ::Val{s}) where {s}
#within the @generated function body, x is a type, but s is a Symbol (since it
#is extracted from a type parameter).
    Core.println("Generated function getindex parsed for type $x")
    blocknumber = findfirst(i->i==s, keys(descriptor(x)))
    blocktype = descriptor(x)[s]
    #getblock is called on the data field and thus dispatches to the BlockArrays method
    :(getblock(getfield(x,:data), $blocknumber)::$blocktype) #enforce return type for stability
end

Base.setindex!(x::Node, v, i::Integer) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::Node, v, i::Colon) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::Node, v, i::AbstractUnitRange) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::Node, v, ::Val{s}) where {s} = (x[Val(s)] .= v)

Base.getproperty(x::Node, s::Symbol) = getindex(x, Val(s))
Base.setproperty!(x::Node, s::Symbol, v) = setindex!(x, v, Val(s))

########### BROADCASTING #############

#if implemented, generate at compile time
# Base.length(::Type{Node{S}}) where {S} = sum(length.(childtypes(Node{S})))
#Base.length(::Node{S}) where {S} = length(Node{S})

# broadcast con un solo argumento
# al efectuar la llamada a broadcast, como solo tengo un argumento, se llama a
# BroadcastStyle() con el tipo de ese argumento.

# Julia coge una broadcasted expression y genera un arbol donde cada nodo es un
# objeto de tipo Broadcasted que contiene la operacion y los operandos. por
# defecto, lo que hace a continuacion es recorrer el arbol determinando el
# BroadcastStyle resultante de cada nodo. Esto lo hace examinando cogiendo por
# ejemplo dos tipos que operan entre si, y llamando a Base.BroadcastStyle(t1,
# t2). el resultado de esta llamada es el tipo Base.Broadcast.BroadcastStyle
# dominante entre ambas expresiones. asi

#como determina Julia el BroadcastStyle resultante de una determinada expresion?
#va recorriendo

#cuando acaba, el resultado es una variable de tipo Broadcasted{WinningStyle}
# de ese  , tengo un veredicto sobre cual debe
#ser el output de la operacion. ese veredicto es un style. cuando se ha
#determinado el ResultingStyle a partir del argumento o argumentos, se hace
#dispatch a similar() con ResultingStyle como argumento. Y ese similar en
#principio debe devolver Y en mi caso, si el
#veredicto es que el ResultingStyle es mi NodeBlockStyle, pues

#como hacer que mi NodeBlockStyle gane frente a otros como AbstractArrayStyle
#cuando la operacion es binaria? pues me defino methods con dos argumentos para
#BroadcastStyle()


#sidenote: writing an expression will all operators dotted, such as x .= x .+ 1,
#means that broadcast! will be called instead of broadcast. as a result, the
#operation will be done in place: instead of allocating a new temporary array
#x_tmp (via a call to similar()) to operate on and copying the result to x,
#Julia reuses the original memory for x. for Node{S}, this is very, very
#interesting, because similar(::Type{Node{S}}) is not efficient: it is a
#recursive procedure, which builds each sub-block from top to bottom. if we can
#spare that work on state vector operations, much better! this will not work if
#any of the operators are not dotted. for example, x .= 3 * x will not be done
#in place. it must be x .= 3.*x

struct NodeBlockStyle{S} <: Broadcast.AbstractArrayStyle{1} end
NodeBlockStyle{S}(::Val{1}) where {S} = NodeBlockStyle{S}()

#as a result of this definition, when a Node{S} is cast together with any other
#AbstractArray, the winning BroadcastStyle is always NodeBlockStyle{S}. it also
#handles the case of broadcast unary operations (see binary broadcasting rules
#in doc to see why), and binary operations with the same Node type parameter.
#however, it pleasantly fails to determine a winner when two Nodes with
#different type parameters are broadcast together in a binary operation
Base.BroadcastStyle(::Type{Node{S}}) where {S} = NodeBlockStyle{S}()


#once the winning BroadcastStyle is determined to be a NodeBlockStyle{S}, this
#defines the resulting data structure that must be produced
Base.similar(::Broadcast.Broadcasted{NodeBlockStyle{S}}, ::Type{ElType}) where {S, ElType} = similar(Node{S})


################ THIS GOES IN A TEST MODULE ###########################

#create a macro that, given a type parameter, generates the descriptor
#function dynamically
const NodeRbd = Node{:rbd}
const NodeLdg = Node{:ldg}
const NodePwp = Node{:pwp}
const NodeAircraft = Node{:aircraft}
const NodeFakeRbd = Node{:fake_rbd}

descriptor(::Type{NodeRbd}) = (att = Leaf{4}, vel = Leaf{3}, pos = Leaf{3})
descriptor(::Type{NodeFakeRbd}) = descriptor(NodeRbd)
descriptor(::Type{NodeLdg}) = (nlg = Leaf{3}, mlg = Leaf{3})
descriptor(::Type{NodePwp}) = (left = Leaf{2}, right = Leaf{2})
descriptor(::Type{NodeAircraft}) = (rbd = NodeRbd, ldg = NodeLdg, pwp = NodePwp)

function lbv_demo()
    att = Leaf{4}(rand(4))
    vel = Leaf{3}(2ones(3))
    pos = Leaf{3}(3ones(3))
    pos = Leaf{3}(3ones(3))
    x_rbd = NodeRbd((att, vel, pos))

    nlg = Leaf{3}(5ones(3))
    mlg = Leaf{3}(2ones(3))
    x_ldg = NodeLdg((nlg, mlg))

    x_pwp = NodePwp(rand(4))

    x_aircraft = NodeAircraft((x_rbd, x_ldg, x_pwp))
    # display(x_aircraft)
    x_rbd_retrieved = x_aircraft[Val(:rbd)]
    x_rbd_retrieved[:] .= 111
    x_aircraft.ldg .= 222
    display(x_aircraft)

    return (
        x_aircraft,
        x_aircraft[Val(:rbd)],
        x_aircraft[Val(:ldg)][:],
        x_aircraft[1:4],
        x_aircraft.rbd.att,
        x_aircraft.ldg[:],
        x_aircraft.ldg.nlg,
    )


end

end