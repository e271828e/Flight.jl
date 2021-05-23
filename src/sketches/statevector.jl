module StateVector

using BlockArrays#, ArrayLayouts
import BlockArrays: axes, viewblock

#additions to export are not tracked by Revise. instead, until Julia is
#restarted, we need to do Flight.LandingGear.ldg_demo
export XLeaf, sv_demo

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
Base.getindex(x::XLeaf, blockindex::BlockIndex) = getindex(x.data, blockindex) #not essential
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
#BlockStateVector constructor. if we use Symbols to generate the NT automatically
#within the constructor, it will be type unstable, because the constructor
#argument values (the Symbols) determine the NT type parameters. A solution
#would be to use Val

#with all those pieces in place, it may be time to return to...
# http://www.stochasticlifestyle.com/zero-cost-abstractions-in-julia-indexing-vectors-by-name-with-labelledarrays/
#... adapting it for my  own case
# i would only need LVector{Syms}, because T is always Float64 and A is always
# BlockVector. no need to parameterize anything else. with this approach, i
# would have to define an LVector{Syms} specific to each Node or Leaf. it may be
# more efficient, but it is more cumbersome and less flexible


function sv_demo()
    x_nlg = XLeaf()
    x_mlg = XLeaf()
    x = mortar([x_nlg, x_mlg])
    x[Block(2)][:att] = [π]
    x[Block(1)][:vel] = [11, 92.4]
    display(x)
end

end #module
