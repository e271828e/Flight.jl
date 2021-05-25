module StateVector

using BlockArrays, StaticArrays
import BlockArrays: axes, viewblock, getblock

#additions to export are not tracked by Revise!
export LBV, lbv_demo

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

#with all those pieces in place, it may be time to return to...
# http://www.stochasticlifestyle.com/zero-cost-abstractions-in-julia-indexing-vectors-by-name-with-labelledarrays/
#... adapting it for my  own case
# i would only need LVector{Syms}, because T is always Float64 and A is always
# BlockVector. no need to parameterize anything else. with this approach, i
# would have to define an LVector{Syms} specific to each Node or Leaf. it may be
# more efficient, but it is more cumbersome and less flexible

#DO NOT USE StaticArrays as blocks! their data cannot be referenced inside BlockVector

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



# about integrating both Labels and BlockTypes in the type parameter itself...
# it is not possible. Julia doc says:

# Both abstract and concrete types can be parameterized by other types. They can
# also be parameterized by symbols, by values of any type for which isbits
# returns true (essentially, things like numbers and bools that are stored like
# C types or structs with no pointers to other objects), and also by tuples
# thereof. Type parameters may be omitted when they do not need to be referenced
# or restricted.

# so, i can use a tuple of symbols representing the labels for the different
# blocks. but nothing more. so, maybe a cleaner solution would be to simply use
# a Symbol identifier for the LBV that helps with dispatch. and then do what i
# was doing: define a block descriptor for each LBV i want to define. then
# dispatch to the correct getblocktype and getblocknumber using that Symbol
# identifier
# what i can do is:
#for Aircraft:
# getblockdescriptor(::LBV{:Aircraft}) = (:rbd = LBV{:Rbd}, :ldg = LBV{:Ldg})
#and now, generic functions:
# getblocknumber(T, s::Symbol) = findfirst(i->i==s, keys(getblockdescriptor(T)))
# getblocktype(T, s::Symbol) = getblockdescriptor(T)[s]

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

#AbstractBlockArray interface
axes(x::LBV) = axes(getfield(x,:data))
viewblock(x::LBV, block) = viewblock(getfield(x, :data), block)

#AbstractArray interface
Base.getindex(x::LBV, i::Integer)::Float64 = getindex(getfield(x,:data), i)
Base.getindex(x::LBV, i::Colon)::Vector{Float64} = getindex(getfield(x,:data),  i)
Base.getindex(x::LBV, i::UnitRange{Int})::Vector{Float64} = getindex(getfield(x,:data),  i)
# Base.getindex(x::LBV, blockindex::BlockIndex{1})::Float64 = getindex(getfield(x,:data), blockindex) #not essential

Base.setindex!(x::LBV, v, i::Integer) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::LBV, v, i::Colon) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::LBV, v, i::UnitRange{Int}) = setindex!(getfield(x,:data), v, i)
# Base.setindex!(x::LBV, v, i::Union{Integer, UnitRange{Int}, Colon}) = setindex!(getfield(x,:data), v, i)

getblock(x::LBV, k::Integer) = getfield(x,:data)[Block(k)]


const RbdLabels = (:att, :vel, :pos)
const LBVRbd = LBV{RbdLabels}

const LdgLabels = (:nlg, :mlg)
const LBVLdg = LBV{LdgLabels}

const AircraftLabels = (:rbd, :ldg)
const LBVAircraft = LBV{AircraftLabels}

# #all these need to be replaced by generated functions. from the value type S,
# which is a #Symbol, it needs to find: which block number corresponds to it
# (this can be #done with findfirst on the tuple of AircraftLabels), and which
# LBV specific subtype corresponds to it (for example, to :rbd #corresponds
# LBV_Rbd = LBV{RbdLabels}, to :ldg corresponds LBV_Ldg, etc

#the key to type stability is that these functions are called at compile time,
# not at runtime, if they were called at runtime, since their output determines
# the types of the extracted blocks... type instability!


# Base.getindex(x::LBVAircraft, ::Val{:rbd})::LBVRbd = getblock(x, 1)
# Base.getindex(x::LBVAircraft, ::Val{:ldg})::LBVLdg = getblock(x, 2)
# Base.getproperty(x::LBVAircraft, s::Symbol) = getindex(x, Val(s))

# Base.getindex(x::LBVRbd, ::Val{:att})::Vector{Float64} = getblock(x, 1)
# Base.getindex(x::LBVRbd, ::Val{:vel})::Vector{Float64} = getblock(x, 2)
# Base.getproperty(x::LBVRbd, s::Symbol) = getindex(x, Val(s))

# Base.getindex(x::LBVLdg, ::Val{:nlg})::Vector{Float64} = getblock(x, 1)
# Base.getindex(x::LBVLdg, ::Val{:mlg})::Vector{Float64} = getblock(x, 2)
# Base.getproperty(x::LBVLdg, s::Symbol) = getindex(x, Val(s))

#these aren't necessary! once the slices or the complete block are retrieved,
#the assignment changes them. but, since they are actually a reference to the
#data held by the parent, this changes that as well.
# Base.setindex!(x::LBVAircraft, v, ::Val{:rbd}) = setindex!(getblock(x, 1), v, :)
# Base.setproperty!(x::LBVAircraft, v, s::Symbol) = setindex!(x, v, Val(s))

Aircraft_block_descriptor = (rbd = LBVRbd, ldg = LBVLdg)
getblocknumber(::Type{LBVAircraft}, s::Symbol) = findfirst(i->i==s, keys(Aircraft_block_descriptor))
getblocktype(::Type{LBVAircraft}, s::Symbol) = Aircraft_block_descriptor[s]

Rbd_block_descriptor = (att = MVector{4, Float64}, vel = MVector{3, Float64}, pos = MVector{3, Float64})
getblocknumber(::Type{LBVRbd}, s::Symbol) = findfirst(i->i==s, keys(Rbd_block_descriptor))
getblocktype(::Type{LBVRbd}, s::Symbol) = Rbd_block_descriptor[s]

Ldg_block_descriptor = (nlg = Vector{Float64}, mlg = Vector{Float64})
getblocknumber(::Type{LBVLdg}, s::Symbol) = findfirst(i->i==s, keys(Ldg_block_descriptor))
getblocktype(::Type{LBVLdg}, s::Symbol) = Ldg_block_descriptor[s]

@generated function Base.getindex(x::LBV, ::Val{s}) where {s}
    #within the @generated function body, x is a type, and s is a Symbol (it is
    #extracted from a type parameter). we could do LBV{T} where {T} and T would
    #be also a Symbol type parameter, but we can dispatch directly on typeof(x)
    Core.println("Generated function getindex parsed for type $x")
    blocknumber = getblocknumber(x, s) #this is a generated function so x is a type!
    blocktype = getblocktype(x, s)
    :(getblock(getfield(x,:data), $blocknumber)::$blocktype) #return it and enforce type
end
Base.getproperty(x::LBV, s::Symbol) = getindex(x, Val(s))


function lbv_demo()
    att = MVector{4, Float64}(rand(4))
    vel = MVector{3, Float64}(2ones(3))
    pos = MVector{3, Float64}(3ones(3))
    x_rbd = LBVRbd((att, vel, pos))

    nlg = 5ones(3)
    mlg = 2ones(3)
    x_ldg = LBVLdg((nlg, mlg))

    x_aircraft = LBVAircraft((x_rbd, x_ldg))
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
        # x_aircraft.ldg.nlg,
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


end #module
