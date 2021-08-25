module LBV

using StaticArrays

export AbstractLBV, LBVLeaf, LBVNode
export @define_node

abstract type AbstractLBV{T, D<:AbstractVector{T}} <: AbstractVector{T} end

#note: it is still unclear whether there is any benefit in allowing zero-length
#blocks. it may allow treating systems more homogeneously, regardless of whether
#they have inputs/states or not. when this becomes clear, remove support for
#zero length blocks if it does not afford anything

########################### LBVLeaf ############################

struct LBVLeaf{L, T, D <: AbstractVector{T}} <: AbstractLBV{T, D}
    data::D
    function LBVLeaf{L,T,D}(data) where {L,T,D} #for some reason, you can't restrict type parameters here!
        if !(L > 0)
            throw(ArgumentError("LBVLeaf length must be positive"))
        end
        if length(data) != length(LBVLeaf{L,T,D})
            throw(ArgumentError("Got input length $(length(data)), expected $(length(LBVLeaf{L,T,D}))"))
        end
        new{L,T,D}(data)
    end
end
LBVLeaf{L,T}(data::D) where {L,T,D} = LBVLeaf{L,T,D}(data) #T/D inconsistencies caught by the inner const
LBVLeaf{L}(data::D) where {L,D} = LBVLeaf{L,eltype(D),D}(data)
LBVLeaf(data::D) where {D} = LBVLeaf{length(data),eltype(data),typeof(data)}(data) #avoid for efficiency

LBVLeaf{L,T}(::UndefInitializer) where {L,T} = LBVLeaf{L,T}(MVector{L,T}(undef))
LBVLeaf{L}(::UndefInitializer) where {L} = LBVLeaf{L,Float64}(undef)

LBVLeaf{L,T}() where {L,T} = LBVLeaf{L,T}(undef)
LBVLeaf{L}() where {L} = LBVLeaf{L,Float64}()

Base.similar(::Type{<:LBVLeaf{L,T}}) where {L,T} = LBVLeaf{L,T}(MVector{L,T}(undef))
Base.similar(::LBVLeaf{L,T}) where {L,T} = similar(LBVLeaf{L,T})
Base.copy(x::LBVLeaf{L}) where {L} = LBVLeaf{L}(copy(getfield(x, :data)))
# Base.convert(::Type{<:LBVLeaf{L}}, v::AbstractVector) where {L} = LBVLeaf{L}(v)

###### Abstract Array #######

Base.length(::Type{<:LBVLeaf{L}}) where {L} = L
Base.size(::LBVLeaf{L}) where {L} = (L,)
Base.getindex(x::LBVLeaf, i) = getindex(getfield(x,:data), i)
Base.setindex!(x::LBVLeaf, v, i) = setindex!(getfield(x,:data), v, i)

####### Custom Broadcasting #######

struct LBVLeafStyle{L,T} <: Broadcast.AbstractArrayStyle{1} end

LBVLeafStyle{L,T}(::Val{1}) where {L,T} = LBVLeafStyle{L,T}()
Base.BroadcastStyle(::Type{LBVLeaf{L,T,D}}) where {L,T,D} = LBVLeafStyle{L,T}()

function Base.similar(::Broadcast.Broadcasted{LBVLeafStyle{L,T}}, ::Type{ElType}) where {L,T,ElType}
    similar(LBVLeaf{L,T})
end

function Base.BroadcastStyle(::LBVLeafStyle{L,T1}, ::LBVLeafStyle{L,T2}) where {L,T1,T2}
    LBVLeafStyle{L,promote_type(T1, T2)}()
end


########################### LBVNode ############################

struct LBVNode{S, T, D <: AbstractVector{T}} <: AbstractLBV{T, D} #S: identifier
    data::D
    function LBVNode{S,T,D}(data::D) where {S,T,D} #for some reason, you can't restrict type parameters here
        assert_symbol(S)
        if length(data) != length(LBVNode{S,T,D})
            throw(ArgumentError("Got input length $(length(data)), expected $(length(LBVNode{S,T,D}))"))
        end
        new{S,T,D}(data)
    end
end
@generated assert_symbol(x) = x<:Symbol ? nothing : :(throw(TypeError(:LBVNode, "inner constructor", Symbol, $x)))

LBVNode{S,T}(data::D) where {S,T,D} = LBVNode{S,T,D}(data) #T/D inconsistencies caught by the inner const
LBVNode{S}(data::D) where {S,D} = LBVNode{S,eltype(D),D}(data)

LBVNode{S,T}(::UndefInitializer) where {S,T} = LBVNode{S,T}(MVector{length(LBVNode{S,T}) ,T}(undef))
LBVNode{S}(::UndefInitializer) where {S} = LBVNode{S,Float64}(undef)

LBVNode{S,T}() where {S,T} = LBVNode{S,T}(undef)
LBVNode{S}() where {S} = LBVNode{S,Float64}()

####### Abstract Array #############

Base.size(x::LBVNode) = size(getfield(x,:data))
Base.getindex(x::LBVNode, i) = getindex(getfield(x,:data), i)
Base.setindex!(x::LBVNode, v, i) = setindex!(getfield(x,:data), v, i)

Base.similar(::Type{<:LBVNode{S,T}}) where {S,T} = LBVNode{S,T}(MVector{length(LBVNode{S,T}),T}(undef))
Base.similar(::LBVNode{S,T}) where {S,T} = similar(LBVNode{S,T})
Base.copy(x::LBVNode{S}) where {S} = LBVNode{S}(copy(getfield(x, :data)))
# Base.convert(::Type{<:LBVNode{S}}, v::AbstractVector) where {S} = LBVNode{S}(v)

####### Custom Broadcasting #######

struct LBVNodeStyle{S,T} <: Broadcast.AbstractArrayStyle{1} end

LBVNodeStyle{S,T}(::Val{1}) where {S,T} = LBVNodeStyle{S,T}()
Base.BroadcastStyle(::Type{LBVNode{S,T,D}}) where {S,T,D} = LBVNodeStyle{S,T}()

function Base.similar(::Broadcast.Broadcasted{LBVNodeStyle{S,T}}, ::Type{ElType}) where {S,T,ElType}
    similar(LBVNode{S,T})
end

function Base.BroadcastStyle(::LBVNodeStyle{S,T1}, ::LBVNodeStyle{S,T2}) where {S,T1,T2}
    LBVNodeStyle{S,promote_type(T1, T2)}()
end


######### Code Generation #########

Base.getproperty(x::LBVNode, s::Symbol) = getproperty(x, Val(s))
Base.setproperty!(x::LBVNode, s::Symbol, v) = setproperty!(x, Val(s), v)

#define and register a LBVNode
macro define_node(name, children)

    ex = quote

        println("Generating code for $($(QuoteNode(name))) = ",
                "LBVNode{:$($(QuoteNode(name)))}...")

        if isdefined($__module__, $(QuoteNode(name)))
            println("WARNING: Type $($(QuoteNode(name))) already defined ",
            "in module $($(QuoteNode(__module__)))")
        end

        #use a let block to create a local scope so that:

        #1) variables captured in methods do not inhabit the global module
        #   scope. this avoids the performance hit from capturing global
        #   variables in a closure. the global keyword is required to make
        #   methods and variables defined within the let block visible in the
        #   global module scope.

        #2) macro hygiene is respected (we have escaped the whole quote)

        let children = $(children)

            global const $(name) = LBVNode{$(QuoteNode(name))}

            println(children)
            child_labels = keys(children)
            child_types = values(children)
            @assert all([t <: Union{LBVLeaf, LBVNode} for t in child_types]) "Invalid child type"
            node_length = sum(length.(child_types))

            global Base.propertynames(::$(name), private::Bool = false) = child_label
            global Base.length(::Type{<:$(name)}) = node_length
            global get_children(::Type{<:$(name)}) = children
            global get_children(::$(name)) = children

            offset = 0
            for (child_label, child_type) in zip(child_labels, child_types)

                child_length = length(child_type)
                child_range = (1 + offset):(child_length + offset)
                offset += child_length
                println("Generating helper methods for child $child_label, range $child_range")

                global function Base.getproperty(x::$(name), ::Val{child_label})
                    return (child_type)(view(getfield(x,:data), child_range))
                end
                global function Base.setproperty!(x::$(name), ::Val{child_label}, v)
                    return setindex!(getfield(x, :data), v, child_range)
                end

            end

            global function Base.show(io::IO, ::MIME"text/plain", x::$(name))

                println(io, "$node_length-element ", typeof(x), " with blocks:")
                for (child_label, child_type) in zip(child_labels, child_types)
                    print(io, "\t", child_label, ": ", getproperty(x, child_label)[:])
                    child_label == child_labels[end] ? break : println()
                end
            end

        end

    end
    return Base.remove_linenums!(esc(ex)) #esc: escapes the whole expression
    # return esc(ex)
end


end #module