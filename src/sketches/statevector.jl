module StateVector

using BlockArrays, StaticArrays
import BlockArrays: axes, viewblock, getblock, setblock!

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

# #o quiza:
# function LBV{:Aircraft}(blocks::values(getdescriptor(LBV{:Aircraft})))
#     LBV{:Aircraft}(blocks)
# end
# #o genericamente:
# function LBV{S}(blocks::values(getdescriptor(LBV{S}))) where {S<:Symbol}
#     LBV{S}(blocks)
# end
# #y el siguiente inner:
# function LBV{S}(blocks) where {S<:Symbol}
#     println("Called block assembly constructor with $(typeof(blocks))")
#     data = mortar(collect(blocks))
#     new{S}(data)
# end

#llegados a este punto, la cosa empieza a oler a que seria interesante generar
#los constructors con metaprogramming. yo defino una macro a la que le paso el
#LBV Identifier, los BlockLabels y los BlockTypes, y el me genera un outer
#constructor adecuado al que enchufarle un Tuple de blocks. ese a su vez se los
#pasa al inner constructor


#CUSTOMIZAR LA REPRESENTACION PARA QUE APAREZCAN LOS NOMBRES DE LOS CHILD BLOCKS



struct LBV{S} <: AbstractBlockVector{Float64}
    data::BlockVector{Float64}
    function LBV{S}(blocks::NTuple{N, Union{LBV, MVector{K, Float64} where {K}, Vector{Float64}}} where{N}) where {S}
        println("Called block assembly constructor with $(typeof(blocks))")
        #ensure the input block types match exactly those prescribed by the
        #parametric type descriptor. TRY TO DO THIS IN COMPILE TIME. USE MACROS
        #TO GENERATE OUTER CONSTRUCTOR
        @assert typeof.(blocks) == values(getdescriptor(LBV{S}))
        data = mortar(collect(blocks))
        new{S}(data)
    end
end

#AbstractBlockArray interface
axes(x::LBV) = axes(getfield(x,:data))
viewblock(x::LBV, block) = viewblock(getfield(x, :data), block)

#AbstractArray interface
Base.getindex(x::LBV, i::Integer)::Float64 = getindex(getfield(x,:data), i)
Base.getindex(x::LBV, i::Colon)::Vector{Float64} = getindex(getfield(x,:data),  i)
Base.getindex(x::LBV, i::AbstractUnitRange)::Vector{Float64} = getindex(getfield(x,:data),  i)
# Base.getindex(x::LBV, blockindex::BlockIndex{1})::Float64 = getindex(getfield(x,:data), blockindex) #not essential

Base.setindex!(x::LBV, v, i::Integer) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::LBV, v, i::Colon) = setindex!(getfield(x,:data), v, i)
Base.setindex!(x::LBV, v, i::AbstractUnitRange) = setindex!(getfield(x,:data), v, i)
# Base.setindex!(x::LBV, v, i::Union{Integer, UnitRange{Int}, Colon}) = setindex!(getfield(x,:data), v, i)

#new
getdescriptor(::Type{LBV}) = error("To be implemented by each LBV parametric type")
getblocknumber(T::Type{LBV{S}} where{S}, s::Symbol) = findfirst(i->i==s, keys(getdescriptor(T)))
getblocktype(T::Type{LBV{S}} where{S}, s::Symbol) = getdescriptor(T)[s]

#within the @generated function body, x is a type, and s is a Symbol (it is
#extracted from a type parameter). we could do LBV{T} where {T} and T would
#be also a Symbol type parameter, but we can dispatch directly on typeof(x)
@generated function Base.getindex(x::LBV, ::Val{s}) where {s}
    Core.println("Generated function getindex parsed for type $x")
    blocknumber = getblocknumber(x, s) #this is a generated function so x is a type!
    blocktype = getblocktype(x, s)
    #getblock dispatches to BlockArrays
    :(getblock(getfield(x,:data), $blocknumber)::$blocktype) #enforce return type for stability
end
#the block may be either a MVector or a LBV. therefore, we cannot access the
#data field for the latter. instead, we need to use the implementation agnostic
#colon notation (we want to set the whole block). first, we leverage the already
#available @generated getindex to get a reference to the block
Base.setindex!(x::LBV, v, ::Val{s}) where {s} = (x[Val(s)][:] = v)

Base.getproperty(x::LBV, s::Symbol) = getindex(x, Val(s))
Base.setproperty!(x::LBV, s::Symbol, v) = setindex!(x, v, Val(s))


################ THIS GOES IN A TEST MODULE ###########################


#maybe slurp input blocks to LBV constructor, more convenient for the caller

const LBVRbd = LBV{:rbd}
const LBVLdg = LBV{:ldg}
const LBVAircraft = LBV{:aircraft}

getdescriptor(::Type{LBVRbd}) = (att = MVector{4, Float64}, vel = MVector{3, Float64}, pos = MVector{3, Float64})
getdescriptor(::Type{LBVLdg}) = (nlg = Vector{Float64}, mlg = Vector{Float64})
getdescriptor(::Type{LBVAircraft}) = (rbd = LBVRbd, ldg = LBVLdg)

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
        x_aircraft.ldg.nlg,
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
