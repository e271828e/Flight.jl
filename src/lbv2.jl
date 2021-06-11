module LBV

export LBVLeaf, LBVNode
export is_registered, descriptor
export register_node

abstract type AbstractLBV{D<:AbstractVector{Float64}} <: AbstractVector{Float64} end

struct LBVLeaf{S,D} <: AbstractLBV{D}
    data::D
    function LBVLeaf{S,D}(data) where {S,D}
        data_length = length(data)
        block_length = length(length(LBVLeaf{S,D}))
        @assert length(data) == length(LBVLeaf{S,D}) "Got input array of length $data_length, expected $block_length"
        new{L,D}(data)
    end
end
#convenience constructor
##################### AVOID FOR GETINDEX, call inner constructor directly with
#the known size
LBVLeaf(data::D) where {D} = LBVLeaf{length(data),D}(data)

Base.length(::Type{<:LBVLeaf{S}}) where {S} = S
#########################THESE ARE PENDING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Base.similar(::btype) = similar($btype)
Base.size(::LBVLeaf{S}) where {S} = (S,)
Base.getindex(x::LBVLeaf, i) = getindex(getfield(x,:data), i)
Base.setindex!(x::LBVLeaf, v, i) = setindex!(getfield(x,:data), v, i)
# Base.similar(::Type{<:$btype}) = $btype(Vector{Float64}(undef, $blength))
# Base.length(::Type{<:$btype}) = $blength



struct LBVNode{L,D} <: AbstractLBV{D}
    data::D
    function LBVNode{L,D}(data) where {L,D}
        data_length = length(data)
        block_length = length(LBVNode{L,D})
        @assert data_length == block_length "Got input array of length $data_length, expected $block_length"
        new{L,D}(data)
    end
end
LBVNode{L}(data::D) where {L,D} = LBVNode{L,D}(data)
LBVNode{L}(input::LBVNode) where {L} = LBVNode{L}(input.data) #conversion between equal length LBVNodes
LBVNode{L}() where {L} = LBVNode{L}(Vector{Float64}(undef, length(LBVNode{L})))

# is_registered(::Type{LBVNode{T,D}} where {T,D}) = false
is_registered(::Type{<:LBVNode}) = false
descriptor(::Type{<:LBVNode}) = error("To be implemented for each type parameter")
descriptor(::T) where {T<:LBVNode}= descriptor(T)

#########################THESE ARE PENDING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Base.similar(::btype) = similar($btype)
Base.size(x::LBVNode) = size(getfield(x,:data))
Base.getindex(x::LBVNode, i) = getindex(getfield(x,:data), i)
Base.setindex!(x::LBVNode, v, i) = setindex!(getfield(x,:data), v, i)
# Base.similar(::Type{<:$btype}) = $btype(Vector{Float64}(undef, $blength))
# Base.length(::Type{<:$btype}) = $blength


#code_generation
function register_node(typepar::Symbol, child_labels::NTuple{N, Symbol},
    child_types::NTuple{N, Any}) where {N}

    @assert all([t <: Union{LBVLeaf, LBVNode} for t in child_types]) "Invalid child type"

    #the descriptor is also be constructed as the right hand side of the
    #descriptor method in the generated code, but as an expression to be
    #evaluated
    desc = NamedTuple{child_labels}(child_types)
    node_length = sum(length.(values(desc)))
    println("Generating code for $typepar = LBVNode{$(QuoteNode(typepar))}...")
    println("Node length: $node_length")

    ex = quote

        if is_registered(LBVNode{$(QuoteNode(typepar))})
            println("\nWARNING: Type $($(QuoteNode(typepar))) already registered")
        end

        #define a convenient shorthand for exporting and method signatures
        const $typepar = LBVNode{$(QuoteNode(typepar))}

        LBV.is_registered(::Type{<:$typepar}) = true
        LBV.descriptor(::Type{<:$typepar}) = NamedTuple{$child_labels}($child_types)

        Base.length(::Type{<:$typepar}) = $node_length

    end

    return Base.remove_linenums!(ex)
    # return ex
end



#to be called only at code generation, not at runtime
function child_ranges(::Type{<:LBVNode{L}}) where {L}
    desc = descriptor(LBVNode{L})
    ranges = Vector{UnitRange{Int}}(undef, length(desc))
    offset = 0
    for (i, l) in enumerate(length.(values(desc)))
        ranges[i] = (1 + offset):(l + offset)
        offset += l
    end
    return NamedTuple{keys(desc)}(Tuple(ranges))
end


# macro register_descriptor(typepar, desc)
#     # println(typeof(typepar))
#     ex = quote
#         $typepar
#         # const typepar = LBVNode{$(QuoteNode($typepar))}
#         LBV.descriptor(::Type{<:LBVNode{$(QuoteNode(typepar))}}) = $desc
#     end
#     return QuoteNode(ex)
# end

# function register_type()

#la expresion que devuelve una macro se evalua inmediatamente, antes de compile
#time. pero claro, si la envuelvo en un QuoteNode, el resultado de esa
#evaluacion seria a su vez una expresion. que todavia quedaria por evaluar. esto
#es lo mismo que estaba consiguiendo con las code generation functions. pero es
#que las expresiones que devuelven esas funciones las estaba envolviendo en un
#eval() para que se ejecutaran en compile time. aqui tendria que hacer lo mismo
#con la expresion no evaluada devuelva por la macro. y por tanto, no gano nada
#respecto de hacerlo con funciones salvo que la macro absorbe el descriptor
#directamente como expresion y me ahorra
# yo puedo hacer que la macro devuelva una expresion sin evaluarla, igual que lo
# hacia con la funcion. pero claro, luego tendria que anadir el eval en el
# momento adecuado

end

# Que LNode sea un tipo concreto, y que este parametrizado por un Symbol que
# nos permita inferir su metadata. Por ejemplo, tendriamos:
# is_registered(::Type{LNode{T}} where {T}) = false
# is_registered(::Type{LNode{:XRbd}}) = true metadata(::Type{LBlock}) = (att =
# Leaf{3},

# La siguiente cuestion es: donde esta almacenado el metadata en este caso? como
# ahora ya no necesito generated functions, puedo usar un method descriptor()
# registrado en LabelledNodeVector, que haga dispatch con el type parameter de
# LNode.

# LabelledBV.desc(::Type{<:LBV} = error("Not implemented") const XRbd =
# LBV{:XRbd} LabelledBV.desc(::Type{XRbd}) = (att = Leaf{4}, vel = Leaf{3}, pos
# = Leaf{4}) const XAircraft = LBV{:XAircraft}
# LabelledBV.desc(::Type{XAircraft}) = (att = Leaf{4}, vel = Leaf{3}, pos =
# Leaf{4})

# Puedo tener un AbstractLBV con subtypes LBVNode{T,D} y LBVLeaf{S,D}. Quiza
# habra algunos methods que antes me definia ad hoc para cada subtipo que
# seguramente se puedan convertir en genericos. Por ejemplo, el constructor? Lo
# unico especifico es la llamada a length. Pero como length lo voy a definir con
# code generation, podria hacer en el inner constructor: function
# LBVNode{L,D}(data::D) where {L,D} @assert length(data) == length(LBVNode{L,D})
# new{L,D}(data) end

# Asi no tengo que estar pendiente de tener descriptors y block types, solo
# block types. A la funcion de code generation le paso simplemente el type
# parameter, y el rhs del descriptor method en formato expression (para no tener
# que reconstruirlo luego): generate_code(:XRbd, :(att = Leaf{4}, vel = Leaf{3},
# pos = Leaf{4}) |> eval

# seguir el mismo camino que hasta ahora: crear las instancias concretas y luego
# generalizar los methods con code generation

# Ojo: aquí Leaf sí que va a almacenar datos. Representa simplemente un bloque
# genérico que por no tener children, no necesita un Symbol como type parameter,
# sólo un Int. Y tampoco necesita, claro, tener methods generados dinamicamente,
# serán genéricos. Pero es importante tener Leafs, porque se pueden usar para
# declarar los tipos de State e Input vectors de Systems, y actúan como
# contenedores válidos de las views de un parent subsystem. Por ejemplo, si un
# muReg define un x::Leaf{2}, cuando vaya a ensamblar el descriptor para generar
# código para el parent system Ldg, lo que haré será definir en el descriptor
# para construir xldg un campo reg=typeof(sys.reg.x). Ese type será Leaf, que es
# un Abstract LBV, que es lo que debo aceptar como válido en un descriptor. Si
# el descriptor va en formato Expr, tendré que interpolar en él.
