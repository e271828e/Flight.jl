# module TestLBV

# using .LabelledBlockVector: descriptor
# using Flight.LabelledBlockVector

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

#an alternative would be to export a macro that generates and executes a
#function that does everything locally: define the descriptor, from the
# @LBV :rbd (att = Leaf{4}, vel = Leaf{3}, pos = Leaf{3})
# @LBV :aircract (rbd = Node{:rbd}, ldg = Node{:ldg})
# #this generates the following code
# descriptor(::)

#descriptor compute the block lengths (this requires )

module Rbd
using Flight.LabelledBlockVector

LabelledBlockVector.descriptor(::Type{Node{:rbd}}) = (att = Leaf{4}, vel = Leaf{3}, pos = Leaf{3})
eval(generate_getindex_sym(Node{:rbd}, :att))
eval(generate_getindex_sym(Node{:rbd}, :vel))

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

function test_lbv()
    data = rand(length(Node{:rbd}))
    x = Node{:rbd}(view(data, :))
    x[Val(:att)]
end

# end #module