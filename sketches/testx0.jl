using Flight
using ComponentArrays

#this works conceptually, but to be done properly it requires a tree. a simple
#Node type with label, data, parent and children with a getproperty overload

function define_x0()
    aaa = nothing
    aab = ones(3)
    aba = ones(2)
    abb = ones(4)
    abc = nothing
    baa = ones(3)
    bab = nothing
    bba = nothing
    bbb = nothing
    bca = fill(2, 3)
    bcb = nothing
    bcc = fill(-5, 2)
    aa = (a = aaa, b = aab)
    ab = (a = aba, b = abb, c = abc)
    ba = (a = baa, b = bab)
    bb = (a = bba, b = bbb)
    bc = (a = bca, b = bcb, c = bcc)
    a = (a = aa, b=ab)
    b = (a = ba, b=bb, c=bc)
    root = (a=a, b=b)
    return root
end

function build_dict(node::Union{NamedTuple, Dict})
    # d = Dict{Symbol, Union{Dict, AbstractVector{<:Real}, Nothing}}()
    d = Dict{Symbol, Any}()
    for k in keys(node)
        d[k] = build_dict(node[k])
    end
    return d
end
build_dict(node::Union{AbstractVector{<:Real}, Nothing}) = node

function assemble_x0(node::Dict)
    blocks = Dict{Symbol, AbstractVector{<:Real}}()
    for label in keys(node)
        x_child = assemble_x0(node[label])
        if x_child !== nothing
            blocks[label] = x_child
        end
    end
    return (!isempty(blocks) ? ComponentVector(blocks) : nothing)
end
assemble_x0(node::AbstractVector) = node
assemble_x0(::Nothing) = nothing


function assign_x0!(node::Dict, x::ComponentVector)
    for label in keys(node)
            println("Called, label $label")
        if label in keys(x)
            assign_x0!(node[label], view(x,label))
            node[label] = view(x, label)
        end
    end
end
assign_x0!(::AbstractVector, ::AbstractVector) = nothing

x0_tree = define_x0() |> build_dict
x0 = assemble_x0(x0_tree)
assign_x0!(x0_tree, x0)
