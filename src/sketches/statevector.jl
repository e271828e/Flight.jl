module StateVector

using BlockArrays, StaticArrays
import BlockArrays: axes, viewblock, getblock

#additions to export are not tracked by Revise. instead, until Julia is
#restarted, we need to do Flight.LandingGear.ldg_demo
export XLeaf, LBV, xleaf_demo, lbv_demo

#setting default values directly:
# https://mauro3.github.io/Parameters.jl/v0.9/manual.html

struct XLeaf <: AbstractBlockVector{Float64}
    data::BlockVector{Float64}
    metadata::@NamedTuple{att::Int, vel::Int}
    #supress all other constructors to ensure initial types and values are correct
    function XLeaf()
        data = mortar([rand(1), ones(2)])
        #improvement: generate metadata automatically from the tuple of Symbols,
        #simply adding their position automatically and putting everything in a
        #Dict or NamedTuple?
        metadata = (att = 1, vel = 2)
        new(data, metadata)
    end
end

#so, the first question is... what the fuck do i need to implement for a
#concrete subtype of the AbstractBlockArray interface? by wading through the
#source code in the BlockArrays repo it is very, very hard to tell. however, a
#reasonable answer is: since BlockArray is a concrete subtype of
#AbstractBlockArray, the methods required must have been implemented by
#BlockArray. going to blockarray.jl, where BlockArray is defined, one finds
#sections "AbstractBlockArray interface", "AbstractArray interface" and
#"Indexing". if we reimplement this methods and simply forward the calls to the
#BlockVector field, all should be fine
#turns out the rest of methods are already implemented either by
#AbstractBlockArray or AbstractArray and do what they need

#these few could be forwarded automatically with a macro

#AbstractBlockArray interface
axes(x::XLeaf) = axes(x.data)
viewblock(x::XLeaf, block) = viewblock(x.data, block)
Base.getindex(x::XLeaf, blockindex::BlockIndex{1}) = getindex(x.data, blockindex) #not essential
#AbstractArray interface
Base.getindex(x::XLeaf, i::Integer) = getindex(x.data, i)
Base.setindex!(x::XLeaf, v, i::Integer) = setindex!(x.data, v, i)


function Base.getindex(x::XLeaf, s::Symbol)
    return x.data[Block(x.metadata[s])]
end

function Base.setindex!(x::XLeaf, v, s::Symbol)
    x.data[Block(x.metadata[s])] = v
end

#next, define XNode <: AbstractBlockVector{Float64} it must define the same
#methods as XLeaf, but its constructor, instead of setting directly data and
#metadata, it must receive an N tuple of AbstractBlockVector{Float64} and
#another of Symbols/Strings. alternatively, we could receive a NamedTuple, see
#what's better for type stability.

#problem: since the keys of a named tuple are part of the type itself,
#generating a NamedTuple for this, unless these keys can be inferred and
#replaced by constants by the compiler, they will lead to type instability
# https://discourse.julialang.org/t/named-tuple-constructor-type-unstable/23461
#here they suggest using either Val or a Dict
#apparently, Val will work as long as it is hardcoded somewhere or can be
#inferred directly from some type parameter:
#https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type

# For instance, we get two instances of XLeaf() and two identifiers (:mlg, :nlg).
#from the two XLeaf instances, we build another one using mortar and store it in
#the data field. for the identifiers, we generate a metadata NamedTuple or Dict
#mapping each now, the good news is that each block of this XNode is itself
#either a XNode or a Leaf. and this means that it holds both its own sub-blocks
#and their idents! so this solves the chained sub-block access by name!

#now, thinking about the commonalities, the only difference between XNode and
#XLeaf is where data and metadata come from. in XLeaf, it is hardcoded in the
#constructor. in XNode they are given as arguments. other than that, the way
#data and metadata are stored and accessed are the same. so it is not out of the
#question to create a single BlockStateVector <: AbstractBlockVector. then,
#there is no need for subclassing it any further. Each module that defines a
#specific BlockStateVector Leaf simply defines a method that assembles it and
#returns an instance of it! the only thing to care for is type stability: ensure
#that the compiler has the information required. see above. either Dict or Val.
#but probably a solution could be: leave the type of metadata open in the struct
#declaration, then make sure to pass a hardcoded NamedTuple on each call to the
#BlockStateVector constructor. if we use Symbols to generate the NT
#automatically within the constructor, it will be type unstable, because the
#constructor argument values (the Symbols) determine the NT type parameters. A
#solution would be to use Val. to enforce having the same number of input blocks
#and input labels, we could define them as NTuple{N, Symbol}, NTuple{N,
#AbstractBlockVector} arguments

#with all those pieces in place, it may be time to return to...
# http://www.stochasticlifestyle.com/zero-cost-abstractions-in-julia-indexing-vectors-by-name-with-labelledarrays/
#... adapting it for my  own case
# i would only need LVector{Syms}, because T is always Float64 and A is always
# BlockVector. no need to parameterize anything else. with this approach, i
# would have to define an LVector{Syms} specific to each Node or Leaf. it may be
# more efficient, but it is more cumbersome and less flexible

#DO NOT USE StaticArrays as blocks! their data cannot be referenced inside BlockVector

struct LBV{Labels} <: AbstractBlockVector{Float64}
    data::BlockVector{Float64}
    #this should never be called with a single element. the raison d'etre for
    #a LBV is assembling multiple blocks of other LBV or Vector{Float64}s!
    function LBV{Labels}(blocks::NTuple{N, Union{LBV, MVector{K, Float64} where {K}, Vector{Float64}}} where{N}) where {Labels}
        println("Called block assembly constructor with $(typeof(blocks))")
        data = mortar(collect(blocks))
        new{Labels}(data)
    end
end

axes(x::LBV) = axes(getfield(x,:data))
viewblock(x::LBV, block) = viewblock(getfield(x, :data), block)

Base.getindex(x::LBV, i::Integer)::Float64 = getindex(getfield(x,:data), i)
Base.getindex(x::LBV, i::Colon)::Vector{Float64} = getindex(getfield(x,:data),  i)
Base.getindex(x::LBV, i::UnitRange{Int})::Vector{Float64} = getindex(getfield(x,:data),  i)
# Base.getindex(x::LBV, blockindex::BlockIndex{1})::Float64 = getindex(getfield(x,:data), blockindex) #not essential

Base.setindex!(x::LBV, v, i::Integer) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::LBV, v, i::Colon) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::LBV, v, i::UnitRange{Int}) = setindex!(getfield(x,:data), v, i)
# Base.setindex!(x::LBV, v, i::Union{Integer, UnitRange{Int}, Colon}) = setindex!(getfield(x,:data), v, i)

#one step further is to tell the compiler exactly WHICH type of LBV will that be
#how? well, if i say that whenever i say x_aircraft[:rbd] the result will be of
#the type LBV_Rbd = LBV{RbdLabels}, that's what i need


const AircraftLabels = (:rbd, :ldg)
const AircraftSV = LBV{AircraftLabels}

const RbdLabels = (:att, :vel)
const RbdSV = LBV{RbdLabels}

const LdgLabels = (:nlg, :mlg)
const LdgSV = LBV{LdgLabels}

#problema: Si AttSV y VelSV no son sino alias para un Vector ordinario, no puede
#haber confusion en los methods getindex y setindex con Val? NO! porque en
#realidad nunca voy a llamar a esos methods si no es para un LBV! porque los
#leafs no van a tener Labels. y por tanto no tiene sentido para ellos. asi que
#basicamente no me molesto en definir AttSV y VelSV. todo lo que sea un leaf,
#Vector{Float64} y sea acbo


getblock(x::LBV, k::Integer) = getfield(x,:data)[Block(k)]

Base.getindex(x::AircraftSV, ::Val{:rbd})::RbdSV = getblock(x, 1)
Base.getindex(x::AircraftSV, ::Val{:ldg})::LdgSV = getblock(x, 2)
Base.getproperty(x::AircraftSV, s::Symbol) = getindex(x, Val(s))

Base.setindex!(x::AircraftSV, v, ::Val{:rbd}) = setindex!(getblock(x, 1), v, :)


Base.getindex(x::RbdSV, ::Val{:att})::Vector{Float64} = getblock(x, 1)
Base.getindex(x::RbdSV, ::Val{:vel})::Vector{Float64} = getblock(x, 2)
Base.getproperty(x::RbdSV, s::Symbol) = getindex(x, Val(s))

Base.getindex(x::LdgSV, ::Val{:nlg})::Vector{Float64} = getblock(x, 1)
Base.getindex(x::LdgSV, ::Val{:mlg})::Vector{Float64} = getblock(x, 2)
Base.getproperty(x::LdgSV, s::Symbol) = getindex(x, Val(s))

# #all these need to be replaced by generated functions. from the value type S,
# which is a #Symbol, it needs to find: which block number corresponds to it
# (this can be #done with findfirst on the tuple of AircraftLabels), and which
# LBV specific subtype corresponds to it (for example, to :rbd #corresponds
# LBV_Rbd = LBV{RbdLabels}, to :ldg corresponds LBV_Ldg, etc

# @generated function Base.getindex(x::LBV{CHANGETHISAircraftLabels}, ::Val{S}) where {S<:Symbol}
#     blocknumber = getblocknumber(x, S) #this is a generated function so x is a type!
#     blocktype = getblocktype(x, S)
#     :(getblock($blocknumber)::$blocktype) #return it and enforce type
# end

# #must create functions getblocknumber & getblocktype that dispatch on the LBV
# specific subtype as the first argument, and on the block label s on the
# second. note that these functions must take a Type as a first argument, not a
# value of a specified type. this is important, because within the generated
# functions, variable names actually denote types. The blocktype can also be a
# MVector, or whatever I want, as long as its constructor accepts the output of
# getblock

#the key to type stability is that these functions are called at compile time,
# not at runtime, if they were called at runtime, since their output determines
# the types of the extracted blocks... type instability!


#THE NEXT STEP: PARAMETRIC TYPE!! Nope. YAGNI.


function lbv_demo()
    att = Vector{Float64}(rand(4))
    att = MVector{4, Float64}(rand(4))
    vel = MVector{3, Float64}(2ones(3))
    x_rbd = RbdSV((att, vel))

    nlg = 5ones(3)
    mlg = 2ones(3)
    x_ldg = LdgSV((nlg, mlg))

    x_aircraft = AircraftSV((x_rbd, x_ldg))
    display(x_aircraft)
    x_rbd_retrieved = x_aircraft[Val(:rbd)]
    x_rbd_retrieved[:] .= 0
    display(x_aircraft)

    return (
        x_aircraft,
        x_aircraft[Val(:rbd)],
        x_aircraft[Val(:ldg)][:],
        x_aircraft[1:4],
        x_aircraft.rbd.att,
        x_aircraft.ldg[:] .+ 2
    )


end
#ahora, la cuestion es, si yo parametrizo LBV solo con el type parameter Labels,
#puedo tener algun conflicto si yo por ejemplo defino un LBV para LandingGear
#que sea LBV{(:right, :left)} y otro para PowerPlant que tambien sea
#LBV{(:right, left)}, pero que cada uno tenga bloques de tamano distinto? En
#teoria podria ser un problema, porque siendo los mismos tipos, pero teniendo
#los mismos tamanos, la funcion getindex(x,:right) deberia ser distinta en ambos
#casos, porque los bloques son de distintos tamanos.

#Pero no deberia. Porque si tengo dos tipos con el mismo nombre en dos modulos,
#cada uno estara cualificado con su #modulo, si no Julia declarara ambiguedad. O
#sea, yo tendre Ldg.LBV{...} y Pwp.LBV{...}. Para cada uno, en sus respectivos
#modulos, se habran generado localmente funciones distintas. Pero hacer un toy
#example por si acaso, con dos tipos triviales del mismo nombre en dos modulos
#distintos, importar los modulos en un cierto codigo y ver que pasa

#gran problema: parece que en cuanto hago un mortar de AbstractBlockArrays
#que no son BlockArrays, se jode el type stability de los elementos que
#extraigo. un arreglo es contener el type instability declarando los tipos tras
#un getindex

#y es que es normal!!! si yo tengo un LBV que esta formado por blocks, cada uno
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
# obtener un LBV{(:att, :vel)} y algo como x.rbd[:att] para obtener el
# Vector{Float64} subyacente, y que como la ultima llamada es a getindex, el
# compilador sabria que el resultado es un Vector{Float64}. pero no funciona,
# porque lo que el compilador no puede saber es que la primera llamada x[:rbd]
# realmente va a devolver un LBV, por lo que no sabe que la getindex a la que se
# invoca despues corresponde realmente a LBV, y por tanto podria no ser esa a la
# que le hemos restringido el return type



function xleaf_demo()
    # x_nlg = XLeaf()
    # x_mlg = XLeaf()
    # x = mortar([x_nlg, x_mlg])
    # x[Block(2)][:att] = [π]
    # x[Block(1)][:vel] = [11, 92.4]
    # x1 = x[Block(1)][:att][1]
    # display(x)
end


end #module
