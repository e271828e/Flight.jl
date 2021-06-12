module LBV

export LBVLeaf, LBVNode
export is_registered, descriptor
export register_LBVNode

abstract type AbstractLBV{D<:AbstractVector{Float64}} <: AbstractVector{Float64} end

########################### LBVLeaf ############################

struct LBVLeaf{L, D <: AbstractVector{Float64}} <: AbstractLBV{D}
    data::D
    #need to implement a constructor with explicit type parameter for extracting
    #and reconstructing Leaf children blocks from the types stored in the
    #descriptor
    function LBVLeaf{L}(data::D) where {L, D} #for some reason, you can't restrict type parameters here!
        if !(L > 0)
            throw(ArgumentError("LBVLeaf length must be positive"))
        end
        if length(data) != length(LBVLeaf{L,D})
            throw(ArgumentError("Got input length $(length(data)), expected $(length(LBVLeaf{L,D}))"))
        end
        new{L, D}(data)
    end
end
LBVLeaf{L}() where {L} = LBVLeaf{L}(Vector{Float64}(undef, L))
LBVLeaf{L}(x::LBVLeaf{L}) where {L} = x
LBVLeaf(x::LBVLeaf) = x
LBVLeaf(data::D) where {D} = LBVLeaf{length(data)}(data)

###### Abstract Array #######

Base.length(::Type{<:LBVLeaf{L}}) where {L} = L
Base.size(::LBVLeaf{L}) where {L} = (L,)
Base.getindex(x::LBVLeaf, i) = getindex(getfield(x,:data), i)
Base.setindex!(x::LBVLeaf, v, i) = setindex!(getfield(x,:data), v, i)
Base.similar(::Type{LBVLeaf{L,D}}) where {L,D} = LBVLeaf{L}(Vector{eltype(D)}(undef, L))
Base.similar(::LBVLeaf{L,D}) where {L,D} = similar(LBVLeaf{L,D})
Base.copy(x::LBVLeaf) = LBVLeaf(copy(getfield(x, :data)))

####### Custom Broadcasting #######

struct LBVLeafStyle{L,D} <: Broadcast.AbstractArrayStyle{1} end

LBVLeafStyle{L,D}(::Val{1}) where {L,D} = LBVLeafStyle{L,D}()
Base.BroadcastStyle(::Type{LBVLeaf{L,D}}) where {L,D} = LBVLeafStyle{L,D}()

function Base.similar(::Broadcast.Broadcasted{LBVLeafStyle{L, D}}, ::Type{ElType}) where {L, D, ElType}
    similar(LBVLeaf{L, D})
end

function Base.BroadcastStyle(::LBVLeafStyle{L,D1}, ::LBVLeafStyle{L,D2}) where {L,D1,D2}
    LBVLeafStyle{L,Vector{promote_type(eltype(D1), eltype(D2))}}()
end


########################### LBVNode ############################

struct LBVNode{S, D <: AbstractVector{Float64}} <: AbstractLBV{D} #S: identifier
    data::D
    function LBVNode{S}(data::D) where {S, D} #for some reason, you can't restrict type parameters here
        assert_symbol(S)
        if length(data) != length(LBVNode{S,D})
            throw(ArgumentError("Got input length $(length(data)), expected $(length(LBVNode{S,D}))"))
        end
        new{S,D}(data)
    end
end
@generated assert_symbol(x) = x<:Symbol ? nothing : :(throw(TypeError(:LBVNode, "inner constructor", Symbol, $x)))
LBVNode{S}(x::LBVNode{S}) where {S} = x
LBVNode{S}(x::LBVNode) where {S} = LBVNode{S}(getfield(x,:data)) #from a LBVNode with different label but equal length
LBVNode{S}() where {S} = LBVNode{S}(Vector{Float64}(undef, length(LBVNode{S})))
LBVNode(x::LBVNode) = x

is_registered(::Type{<:LBVNode}) = false
descriptor(::Type{<:LBVNode}) = error("To be implemented for each type parameter")
descriptor(::T) where {T<:LBVNode}= descriptor(T)

####### Abstract Array #############
################################### Base.@propagate_inbounds for Val(s methods)
###and @_propagate_inbounds_meta aand getindex for array of symbols
Base.size(x::LBVNode) = size(getfield(x,:data))
Base.getindex(x::LBVNode, i) = getindex(getfield(x,:data), i)
Base.setindex!(x::LBVNode, v, i) = setindex!(getfield(x,:data), v, i)
Base.similar(::Type{LBVNode{S, D}}) where {S, D} = LBVNode{S}(Vector{eltype(D)}(undef, length(LBVNode{S, D})))
Base.similar(::LBVNode{S, D}) where {S, D} = similar(LBVNode{S, D})
Base.copy(x::LBVNode{S}) where {S} = LBVNode{S}(copy(getfield(x, :data)))

Base.getproperty(x::LBVNode, s::Symbol) = getindex(x, Val(s))
Base.setproperty!(x::LBVNode, s::Symbol, v) = setindex!(x, v, Val(s))

####### Custom Broadcasting #######

struct LBVNodeStyle{S,D} <: Broadcast.AbstractArrayStyle{1} end
LBVNodeStyle{S, D}(::Val{1}) where {S, D} = LBVNodeStyle{S, D}()
Base.BroadcastStyle(::Type{LBVNode{S, D}}) where {S, D} = LBVNodeStyle{S, D}()

function Base.similar(::Broadcast.Broadcasted{LBVNodeStyle{S, D}}, ::Type{ElType}) where {S, D, ElType}
    similar(LBVNode{S, D})
end

function Base.BroadcastStyle(::LBVNodeStyle{S,D1}, ::LBVNodeStyle{S,D2}) where {S,D1,D2}
    LBVNodeStyle{S,Vector{promote_type(eltype(D1), eltype(D2))}}()
end

######### Code Generation #########

function register_LBVNode(typepar::Symbol, child_labels::NTuple{N, Symbol},
    child_types::NTuple{N, Any}) where {N}

    @assert all([t <: Union{LBVLeaf, LBVNode} for t in child_types]) "Invalid child type"

    #the descriptor is also constructed as the right hand side of the descriptor
    #method in the generated code, but as an expression
    desc = NamedTuple{child_labels}(child_types)
    node_length = sum(length.(values(desc)))
    println("Generating code for $typepar = LBVNode{$(QuoteNode(typepar))}...")

    ex = Expr(:block) #equivalent to ex = quote end

    push!(ex.args, quote

        if is_registered(LBVNode{$(QuoteNode(typepar))})
            println("WARNING: Type $($(QuoteNode(typepar))) already registered")
        end

        #shorthand for exporting and more convenient method signatures
        const $typepar = LBVNode{$(QuoteNode(typepar))}

        LBV.is_registered(::Type{<:$typepar}) = true
        LBV.descriptor(::Type{<:$typepar}) = NamedTuple{$child_labels}($child_types)
        Base.length(::Type{<:$typepar}) = $node_length

    end)

    offset = 0
    for (child_label, child_type) in zip(keys(desc), values(desc))

        child_length = length(child_type)
        child_range = (1 + offset):(child_length + offset)
        offset += child_length

        push!(ex.args, quote

            println("Generating getindex and setindex Symbol methods for child "*
            ":$($(QuoteNode(child_label))), range $($child_range)... ")
            Base.getindex(x::$typepar, ::Val{$(QuoteNode(child_label))}) = ($child_type)(view(getfield(x,:data), $child_range))
            Base.setindex!(x::$typepar, v, ::Val{$(QuoteNode(child_label))}) = setindex!(getfield(x, :data), v, $child_range)

        end)

    end

    # return Base.remove_linenums!(ex)
    return ex
end

end