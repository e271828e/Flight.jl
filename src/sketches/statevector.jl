module StateVector

using BlockArrays, StaticArrays
import BlockArrays: axes, viewblock, getblock, setblock!

#additions to export are not tracked by Revise!
export Node, lbv_demo

#setting default values directly:
# https://mauro3.github.io/Parameters.jl/v0.9/manual.html

#problem: since the keys of a named tuple are part of the type itself,
#generating a NamedTuple for this, unless these keys can be inferred and
#replaced by constants by the compiler, they will lead to type instability
# https://discourse.julialang.org/t/named-tuple-constructor-type-unstable/23461
#here they suggest using either Val or a Dict
#apparently, Val will work as long as it is hardcoded somewhere or can be
#inferred directly from some type parameter:
#https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type

#so, the first question is... what the fuck do i need to implement for a
#concrete subtype of the AbstractBlockArray interface? by wading through the
#source code in the BlockArrays repo it is very, very hard to tell. however, a
#reasonable answer is: since BlockArray is a concrete subtype of
#AbstractBlockArray, the methods required must have been implemented by
#BlockArray. going to blockarray.jl, where BlockArray is defined, one finds
#sections "AbstractBlockArray interface", "AbstractArray interface" and
#"Indexing". if we reimplement this methods and simply forward the calls to the
#BlockVector field, all should be fine turns out the rest of methods are already
#implemented either by AbstractBlockArray or AbstractArray and do what they need

#####################################################
# about integrating both Labels and BlockTypes in the type parameter itself...
# it is not possible. Julia doc says:

# Both abstract and concrete types can be parameterized by other types. They can
# also be parameterized by symbols, by values of any type for which isbits
# returns true (essentially, things like numbers and bools that are stored like
# C types or structs with no pointers to other objects), and also by tuples
# thereof. Type parameters may be omitted when they do not need to be referenced
# or restricted.


#customizar la representacion para que aparezcan los nombres de los child blocks

#when operating with a complete Node, use [:] to output a Vector{Float64} view
#and ensure type stability.

#in Julia, x[:] =... should be replaced by x .= ...
#https://julialang.org/blog/2017/01/moredots/

#todo: need a way to restrict the specification of the descriptor. ONLY CAN HAVE
#Node OR LEAF!! Maybe create a new type called NodeDescriptor, whose constructor
#only accepts Symbols and Union{Node, Leaf}
#only if i do this can i relax the constructor argument types from Union{Leaf,
#Node} to AbstractVector. otherwise, someone could call the block assembly
#constructor with a regular Vector{Float64}, which has no similar() method


#TODO: currently, broadcasting externally breaks the types
#investigate:
#https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
#i may find that after implementing that BroadcastStyle stuff i no longer need
#to call similar() inside operator overloads

Leaf{K} = MVector{K, Float64} where {K}

struct Node{S} <: AbstractBlockVector{Float64}
    data::BlockVector{Float64}
    function Node{S}(blocks::NTuple{N, Union{Node, Leaf}} where {N}) where {S}
        println("Called block assembly constructor for $S")
        properblocktypes=getblocktypes(Node{S})
        println("Converting inputs of type $(typeof.(blocks)) to $properblocktypes")
        converted_blocks = Tuple(convert(T, b) for (T, b) in zip(properblocktypes, blocks))
        data = mortar(collect(converted_blocks))
        new{S}(data)
    end
end

getdescriptor(::Type{Node}) = error("To be implemented by each Node type")
getblocktypes(T::Type{Node{S}} where{S}) = values(getdescriptor(T))

#the goal is to construct a new Node instance that respects the specific types at
#each node of the Node hierarchy. we do this recursively from top to bottom, by
#assemblying new instances of the types indicated by the Node descriptor at each
#level. a way to do this is have a similar() method that dispatches on type
#rather than on an instance, and call it recursively from the top. any Node
#hierarchy will be made of Node nodes and MVectors at the leafs. MVectors already
#have a similar() method dispatching on type. to achieve the above, we need
#another one for Node, which will implement the recursion
Base.similar(::Type{Node{S}}) where {S} = Node{S}(similar.(getblocktypes(Node{S})))
Base.similar(::Node{S}) where {S} = similar(Node{S})

Node{S}(v::AbstractVector) where {S} = (x = similar(Node{S}) ; x .= v)
Node{S}() where {S} = (x = similar(Node{S}) ; x .= 0)

Base.convert(::Type{Node{S}}, v::Node{S}) where{S} = (println("No need to convert $S"); v) #do nothing!
Base.convert(::Type{Node{S}}, v::AbstractVector) where {S} = (println("Converting $(typeof(v)) to Node{$S}"); Node{S}(v))

#TO COMPUTE AT COMPILE TIME!!!!!!!!!!!!!!!!!!!
Base.length(::Type{Node{S}}) where {S} = sum(length.(getblocktypes(Node{S})))
Base.copy(x::Node{S}) where {S} = Node{S}(Tuple(copy.(blocks(getfield(x,:data)))))
#Base.length(::Node{S}) where {S} = length(Node{S})

#AbstractBlockArray interface
axes(x::Node) = axes(getfield(x,:data))
viewblock(x::Node, block) = viewblock(getfield(x, :data), block)

#AbstractArray interface
Base.getindex(x::Node, i::Integer)::Float64 = getindex(getfield(x,:data), i)
Base.getindex(x::Node, i::Colon)::Vector{Float64} = getindex(getfield(x,:data),  i)
Base.getindex(x::Node, i::AbstractUnitRange)::Vector{Float64} = getindex(getfield(x,:data),  i)

#within the @generated function body, x is a type, and s is a Symbol (it is
#extracted from a type parameter). we could do Node{T} where {T} and T would
#be also a Symbol type parameter, but we can dispatch directly on typeof(x)
@generated function Base.getindex(x::Node, ::Val{s}) where {s}
    Core.println("Generated function getindex parsed for type $x")
    blocknumber = findfirst(i->i==s, keys(getdescriptor(x)))
    blocktype = getdescriptor(x)[s]
    #getblock dispatches to BlockArrays
    :(getblock(getfield(x,:data), $blocknumber)::$blocktype) #enforce return type for stability
end

Base.setindex!(x::Node, v, i::Integer) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::Node, v, i::Colon) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::Node, v, i::AbstractUnitRange) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::Node, v, ::Val{s}) where {s} = (x[Val(s)] .= v)

Base.getproperty(x::Node, s::Symbol) = getindex(x, Val(s))
Base.setproperty!(x::Node, s::Symbol, v) = setindex!(x, v, Val(s))

#these non-broadcast operators preserve the Node hierarchy
Base.:(+)(x1::Node{S}, x2::Union{Real, Node{S}}) where {S} = (y = similar(Node{S}) ; y .= x1 .+ x2)
Base.:(-)(x1::Node{S}, x2::Union{Real, Node{S}}) where {S} = (y = similar(Node{S}) ; y .= x1 .- x2)
Base.:(*)(x1::Node{S}, a::Real) where {S} = (y = similar(Node{S}) ; y .= x1 .* a)
Base.:(*)(a::Real, x1::Node{S}) where {S} = x1 * a
Base.:(/)(x1::Node{S}, a::Real) where {S} = (y = similar(Node{S}) ; y .= x1 ./ a)


################ THIS GOES IN A TEST MODULE ###########################

const NodeRbd = Node{:rbd}
const NodeLdg = Node{:ldg}
const NodePwp = Node{:pwp}
const NodeAircraft = Node{:aircraft}

getdescriptor(::Type{NodeRbd}) = (att = Leaf{4}, vel = Leaf{3}, pos = Leaf{3})
getdescriptor(::Type{NodeLdg}) = (nlg = Leaf{3}, mlg = Leaf{3})
getdescriptor(::Type{NodePwp}) = (left = Leaf{2}, right = Leaf{2})
getdescriptor(::Type{NodeAircraft}) = (rbd = NodeRbd, ldg = NodeLdg, pwp = NodePwp)

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
    display(x_aircraft)
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




















#ahora, la cuestion es, si yo parametrizo Node solo con el type parameter Labels,
#puedo tener algun conflicto si yo por ejemplo defino un Node para LandingGear
#que sea Node{(:right, :left)} y otro para PowerPlant que tambien sea
#Node{(:right, left)}, pero que cada uno tenga bloques de tamano distinto? En
#teoria podria ser un problema, porque siendo los mismos tipos, pero teniendo
#los mismos tamanos, la funcion getindex(x,:right) deberia ser distinta en ambos
#casos, porque los bloques son de distintos tamanos.

#Pero no deberia. Porque si tengo dos tipos con el mismo nombre en dos modulos,
#cada uno estara cualificado con su #modulo, si no Julia declarara ambiguedad. O
#sea, yo tendre Ldg.Node{...} y Pwp.Node{...}. Para cada uno, en sus respectivos
#modulos, se habran generado localmente funciones distintas. Pero hacer un toy
#example por si acaso, con dos tipos triviales del mismo nombre en dos modulos
#distintos, importar los modulos en un cierto codigo y ver que pasa

#gran problema: parece que en cuanto hago un mortar de AbstractBlockArrays
#que no son BlockArrays, se jode el type stability de los elementos que
#extraigo. un arreglo es contener el type instability declarando los tipos tras
#un getindex

#y es que es normal!!! si yo tengo un Node que esta formado por blocks, cada uno
#de los cuales es un tipo, y por tanto cada vez que llamo a setindex con un
#Symbol como argumento devuelvo un block de tipo distinto, el compilador no
#puede saber exactamente que tipo de retorno se va a encontrar en cada caso. es
#un caso claro de type instability: el tipo de retorno depende del argumento de
#la funcion. no tengo muy claro el impacto de esto en el codigo generado. es
#decir, yo podria perfectamente implementar un concrete subtype de
#AbstractBlockVector{Float64} que almacenara una serie de strings, y que a
#traves de una logica perversa, cuando llamo a sus getproperty o setproperty,
#devuelvan alguno de esos strings en vez de un Float64. y eso hace ademas que la
#cosa se propague. aunque yo tras la primera llamada a getproperty devuelva un
#subtipo de AbstractBlockVector{Float64}, como el compilador no lo sabe,
#cualquier llamada posterior a getproperty para ese objeto volvera a ser
#imprededcible conclusion: a menos que

# se podria pensar que una alternativa para lograr type stability es preservar
# el return type real (y con el la estructura anidada de todo lo que hay por
# debajo) para get_property y forzar el return type a Vector{Float64} para
# getindex, de manera que entonces se podria hacer algo como x[:rbd] para
# obtener un Node{(:att, :vel)} y algo como x.rbd[:att] para obtener el
# Vector{Float64} subyacente, y que como la ultima llamada es a getindex, el
# compilador sabria que el resultado es un Vector{Float64}. pero no funciona,
# porque lo que el compilador no puede saber es que la primera llamada x[:rbd]
# realmente va a devolver un Node, por lo que no sabe que la getindex a la que se
# invoca despues corresponde realmente a Node, y por tanto podria no ser esa a la
# que le hemos restringido el return type


end #module
