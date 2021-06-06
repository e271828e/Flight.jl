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

macro make_constructor(nodetype)
    println(typeof(nodetype))
    println(:($nodetype))
    d = descriptor(Node{:($nodetype)})
    println(d)
    return quote
        $(esc(Node)){$nodetype}() = println("Dummy function")
    end
end



function generate_getindex(::Type{Node{:rbd}}, s::Symbol)
    d = descriptor(Node{:rbd})
    block_length = length(d[s])
    println("Generating getindex for Val($s)")
    sym = QuoteNode(s)
    return :(Base.getindex(x::Node{:rbd}, ::Val{$sym}) = getindex(getfield(x,:data), 1:$block_length))
end

eval(generate_getindex(Node{:rbd}, :att))
eval(generate_getindex(Node{:rbd}, :vel))




# function Base.getindex(x::Node, ::Val{:att})
#     return view(getfield(x,:data), Block(1))
# end

# function Base.getindex(x::Node, ::Val{:vel})
#     return view(getfield(x,:data), Block(2))
# end

function test_lbv()
    x_rbd = Node{:rbd}()
    x_ldg = Node{:ldg}()
    x_aircraft = Node{:aircraft}()
    return x_aircraft
end

# end #module