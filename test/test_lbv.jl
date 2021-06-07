#when using "using", we need the LabelledBlockVector qualifier to extend. see:
#https://docs.julialang.org/en/v1/manual/modules/#using-and-import-with-specific-identifiers,-and-adding-methods
#maybe create a macro to make this a bit more convenient


#problem: to construct the LBV hierarchy and the associated methods, the LBV
#functions need to have access from the LBV module to the descriptor methods for
#each of the Node parametric types contained in the hierarchy. this can only be
#achieved (AFAIK) by extending the descriptor() method originally defined in the
#LBV module. otherwise, these descriptor methods will not be in scope of the LBV
#methods that require them.

# however, if each module using the LBV module defines a method
#descriptor(::Type{Node{T}}) for its own parameter T, then there is the
#potential for one module overwritting another module's descriptor if they
#happen to choose the same parameter T. this is not solved by making Node
#abstract and requiring each module to define its own concrete subtype, because
#similary there is the possibility of module defining different implementations
#of descriptor(::Type{MyNodeSubtype}), because both think that MyNodeSubtype is
#being defined by nobody else.

#OPTION 1: create and export a macro from LBV that generates code function that
#does everything locally. for example, in module Aircraft I would do:
# @LBV :aircraft (rbd = Node{:rbd}, ldg = Node{:ldg})
# #this generates the following code
# descriptor(::Type{Node{:aircraft}}) = (rbd = Node{:rbd}, ldg = Node{:ldg})
# blocklengths(::Type{Node{:aircraft}}) = length.(values(descriptor(Node{:aircraft})))
# totallength(::Type{Node{:aircraft}}) = sum(blocklengths(Node{:aircraft}))
# function get_methods(::Type{Node{:aircraft}})
#     #getindex and setindex from block identifiers
# end
#now, for the above code to run, it needs access to
#descriptor(::Type{Node{:rbd}}) and descriptor(::Type{Node{:ldg}}). unless these
#are explictly requested by the Aircraft module by writing "using Rbd" and
#"using Ldg", they will not be available. this requires in turn that Rbd and Ldg
#export their respective descriptor() methods.

#OPTION 2: extend the descriptor method in LabelledBlockVector, but do so with a
#macro that checks if a method with the same signature already exists, and if
#so, raise an error. this is the simplest one, and the most secure. each module
#can then export its descriptor method

#OPTION 3: use a macro to define an alias const MyNodeType for Node{:T}. Then
#create a unique type parameter (mangled) T to avoid clashes. if some other
#module creates the same MyNodeType, hopefully the constants will clash. the
#problem is that behind the scenes, we have Node{:Tmangled}, so we will need to
#also define show methods that display MyNodeType. each module exports its own
#NodeType, so that I can do
# @LBV XAircraft (rbd = XRbd, ldg = XLdg)

#OPTION 4: define Node as an abstract type. use "global" to prevent multiple
#definitions of the same method

#OPTION 5: the descriptor should not be methods, but types!! that is what should
#be exported. what I need to do is pass the descriptor itself to the functions
#in LBV

# Pero cuidado: me estoy preocupando mucho de los conflicting names con los
# descriptors, cuando en última instancia voy a tener extensions a getindex y
# setindex en Base para los diferentes type parameters o concrete LBV subtypes.
# El problema NO es sólo conflicto en LBV, sino en Base! Y esto no voy a poder
# evitarlo! Sólo puedo asegurarme de que no se ha definido otro Node con el
# mismo type parameter. Comprobar methods.

# En LabelledArrays no es posible definir LArrays con el mismo type parameter y
# distinto comportamiento porque realmente el type parameter contiene toda la
# especificidad de ese subtipo. Si el type parameter es igual, el subtipo es por
# definicion igual. Aqui no es el caso.
module TestLBV

module Rbd
using Flight.LabelledBlockVector


LabelledBlockVector.descriptor(::Type{Node{:rbd}}) = (att = Leaf{4}, vel = Leaf{3}, pos = Leaf{3})
eval(generate_getindex_sym(Node{:rbd}, :att))
eval(generate_getindex_sym(Node{:rbd}, :vel))
eval(generate_getindex_sym(Node{:rbd}, :pos))

end

module Ldg
using Flight.LabelledBlockVector
export descriptor

LabelledBlockVector.descriptor(::Type{Node{:ldg}}) = (mlg = Leaf{3}, nlg = Leaf{3})
# LabelledBlockVector.descriptor(::Type{Node{:rbd}}) = (mlg = Leaf{3}, nlg = Leaf{3})

end

module Aircraft
using Flight.LabelledBlockVector

#this module should not redefine LabelledBlockVector.descriptor(Node{:rbd}).
#that is the definition of type piracy. the problem is i cannot know if another
#module has extended LabelledBV with the descriptor for the same type parameter
#that i am extending as well. alternativa: definir Node como abstract type (casi
#todos los metodos son agnosticos respecto de los datos), y definir un
#AircraftNode{D}, donde data::D y D<:AV{F64}es el unico type parameter. lo que
#si deberia definir es un method getdata para sustituir a getfield(x,:data).
#este getdata debe ser overloaded por cada concrete subtype. con una macro puedo
#generar la definicion struct, el descriptor la funcion getdata y ya esta. no
#necesito mas
LabelledBlockVector.descriptor(::Type{Node{:aircraft}}) = (rbd = Node{:rbd}, ldg = Node{:ldg})

end
#
#STEP 2:
#READ HANDS ON DESIGN PATTERNS and understand how to interpolate symbols into
#expressions. either a macro or eval could work. no, IT MUST BE EVAL!

using Flight.LabelledBlockVector
using .Aircraft
export test_lbv

function test_lbv()
    data = rand(length(Node{:rbd}))
    x = Node{:rbd}(view(data, :))
    # y = Node{:aicraft}(view(rand(16), :))
    @show x[Val(:att)]
    # y[Val(:)]
end

end #module