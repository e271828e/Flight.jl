# module TestLBV

using Test
using LinearAlgebra

includet("../src/lbv.jl")
using .LabelledBlockVector
# using .LabelledBlockVector: descriptor
# using Flight.LabelledBlockVector

#when using "using", we need the LabelledBlockVector qualifier to extend. see:
#https://docs.julialang.org/en/v1/manual/modules/#using-and-import-with-specific-identifiers,-and-adding-methods
#maybe create a macro to make this a bit more convenient

LabelledBlockVector.descriptor(::Type{Node{:rbd}}) = (att = Leaf{4}, vel = Leaf{3}, pos = Leaf{3})
LabelledBlockVector.descriptor(::Type{Node{:ldg}}) = (mlg = Leaf{3}, nlg = Leaf{3})
LabelledBlockVector.descriptor(::Type{Node{:aircraft}}) = (rbd = Node{:rbd}, ldg = Node{:ldg})
#
#STEP 2:
#READ HANDS ON DESIGN PATTERNS and understand how to interpolate symbols into
#expressions. either a macro or eval could work. no, IT MUST BE EVAL!

eval(generate_getindex_sym(Node{:rbd}, :att))
eval(generate_getindex_sym(Node{:rbd}, :vel))


function test_lbv()
    data = rand(length(Node{:rbd}))
    x = Node{:rbd}(view(data, :))
    x[Val(:att)]
end

# end #module