module LabelledBlockVector

export Node, Leaf, Empty, LBlock
export descriptor, is_registered, register_type



########################## LBV Descriptor #################

abstract type LBVDescriptor end

#convenience subtype to indicate a zero length LBV
struct Empty <: LBVDescriptor end
block_length(::Empty) = 0

struct Leaf{S} <: LBVDescriptor end
Leaf(S::Integer) = (@assert S>0 "Leaf cannot have zero length"; Leaf{S}())
block_length(::Leaf{S}) where {S} =  S

struct Node{T} <: LBVDescriptor #T: type of the associated block, N: number of children
    children::NamedTuple
end
Node(label:: Symbol, nt) = Node{label}(nt)
block_length(n::Node) = sum(block_length.(values(n.children)))
block_type(::Node{T}) where {T} = T

function child_ranges(n::Node)
    child_lengths = block_length.(values(n.children))
    blk_ranges = Vector{UnitRange{Int}}(undef, length(n.children))
    offset = 0
    for (i, l) in enumerate(child_lengths)
        blk_ranges[i] = (1 + offset):(l + offset)
        offset += l
    end
    return NamedTuple{keys(n.children)}(Tuple(blk_ranges))
end

Base.getindex(n::Node, s::Symbol) = getindex(n.children, s)



######################## LBlock ############################

#type parameter D enables different underlying data types (vectors and views)
abstract type LBlock{D<:AbstractVector{Float64}} <: AbstractVector{Float64} end

is_registered(::Type{<:Any}) = false
descriptor(::Type{<:Any}) = error("Must be implemented by concrete LBlock subtypes")

#create custom show methods. these can use the descriptor to know which blocks
#we have, their block ranges, etc!

register_type(::Leaf) = :() #a leaf does not require specific code
function register_type(desc::Node)

    btype = block_type(desc)
    blength = block_length(desc)
    ranges = child_ranges(desc)

    ex = Expr(:block) #equivalent to ex = quote end

    ex_const = quote

        println("Generating constructors for type $($(QuoteNode(btype)))... ")
        println("Length: $($blength)")
        println("Children: $(typeof($desc.children))")

        struct $btype{D} <: LBlock{D}
            data::D
            function $btype{D}(data::D) where {D}
                @assert length(data) == $blength "Expected an input of length $($blength)"
                new{D}(data)
            end
        end
        $btype(input::D) where {D<:AbstractVector{Float64}} = $btype{D}(input)
        $btype(input::$btype) = input
        $btype() = $btype(Vector{Float64}(undef, $blength))

        if is_registered($btype)
            println("\nWARNING: Type $($(QuoteNode(btype))) already registered")
        end

        LabelledBlockVector.is_registered(::Type{<:$btype}) = true
        LabelledBlockVector.descriptor(::Type{T} where {T<:$btype}) = Node{$(QuoteNode(btype))}($(desc.children))

    end

    ex_array = quote

        println("Generating AbstractArray interface for type $($(QuoteNode(btype)))... ")
        Base.similar(::$btype) = similar($btype)
        Base.size(::$btype) = ($blength,)
        Base.getindex(x::$btype, i) = getindex(getfield(x,:data), i)
        Base.setindex!(x::$btype, v, i) = setindex!(getfield(x,:data), v, i)
        Base.similar(::Type{<:$btype}) = $btype(Vector{Float64}(undef, $blength))
        Base.length(::Type{<:$btype}) = $blength

    end

    ex_getsetindex = quote end

    for (child_label, child_desc, _range) in zip(keys(desc.children), desc.children, ranges)

        push!(ex_getsetindex.args, quote
            println("Generating getindex and setindex Symbol methods for type $($(QuoteNode(btype))), child $($(QuoteNode(child_label))), range $($_range)... ")
        end)

        if isa(child_desc, Leaf)
            push!(ex_getsetindex.args, quote
                Base.getindex(x::$btype, ::Val{$(QuoteNode(child_label))}) = view(getfield(x,:data), $_range)
                Base.setindex!(x::$btype, v, ::Val{$(QuoteNode(child_label))}) = setindex!(getfield(x, :data), v, $_range)
            end)
        else #isa(child_desc, Node)
            push!(ex_getsetindex.args, quote
                Base.getindex(x::$btype, ::Val{$(QuoteNode(child_label))}) = $(block_type(child_desc))(view(getfield(x,:data), $_range))
                Base.setindex!(x::$btype, v, ::Val{$(QuoteNode(child_label))}) = setindex!(getfield(x, :data), v, $_range)
            end)
        end

    end

    ex_getsetproperty = quote

        println("Generating getproperty and setproperty Symbol methods for type $($(QuoteNode(btype)))... ")
        Base.getproperty(x::$btype, s::Symbol) = getindex(x, Val(s))
        Base.setproperty!(x::$btype, s::Symbol, v) = setindex!(x, v, Val(s))

    end

    bstyle = (Meta.parse("$(btype)Style"))
    ex_broadcast = quote
        struct $bstyle{D} <: Broadcast.AbstractArrayStyle{1} end
        $bstyle{D}(::Val{1}) where {D} = $bstyle{D}()
        Base.BroadcastStyle(::Type{$btype{D}}) where {D} = $bstyle{D}()
        function Base.similar(::Broadcast.Broadcasted{$bstyle{D}}, ::Type{ElType}) where {D, ElType}
            similar($btype{D})
        end
    end

    push!(ex.args, ex_const, ex_array, ex_getsetindex, ex_getsetproperty, ex_broadcast)

    return Base.remove_linenums!(ex)
    # return ex
end


end

#when a System has children and it is initialized, it checks which of its
#children's descriptors are empty. for those that are not, it assigns to their
#state vector a view of its own state vector, as self.children[:child_id].x =
#self.x.child_id (which is a view, as defined by the Symbol getindex method).
#this should propagate recursively down to the leaf subsystems

#en el momento de crear un sistema, si se le pasan  child subsystems, el
#constructor hara assign_subsystem_states!(sys). cada uno de esos child
#subsystems ya estara creado y tendra su propio x local. entonces:

#function assign_subsystem_states!(aircraft)
#aircraft.rbd.x = aircraft.x[:rbd]
#assign_subsystem_states!(aircraft.rbd)
#aircraft.ldg.x = aircraft.x[:ldg]
#assign_subsystem_states!(aircraft.ldg)
#end

#el problema de hacer assign_subsystem_states en el constructor es que, al ser
#una funcion recursiva, habra multiples pasadas por un mismo subsistema. por
#ejemplo, att sera asignado cuando construya rbd. pero cuando construya
#aircraft, asignara a rbd, y este a su vez recursivamente a att. esto en la
#practica no es un problema, porque solo va a ocurrir al crear la estructura del
#sistema, no en simulacion. pero si es poco estetico, se puede evitar haciendo
#que esta asignacion ocurra no en el propio constructor, sino con la funcion
#initialize!(sys). y la primera vez que haga step tambien se puede llamar a
#initialize state