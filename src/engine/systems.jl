module Systems

using ComponentArrays
using DataStructures
using AbstractTrees

export Component, System
export SystemẊ, SystemX, SystemY, SystemU, SystemS
export init_ẋ, init_x, init_y, init_u, init_s
export f_ode!, f_step!, f_disc!, update_y!


################################################################################
############################## Component ################################

abstract type Component end

function DataStructures.OrderedDict(g::Component)
    fields = propertynames(g)
    values = map(λ -> getproperty(g, λ), fields)
    OrderedDict(k => v for (k, v) in zip(fields, values))
end

################################################################################
############################## SystemTrait #####################################

abstract type SystemTrait end

struct SystemẊ <: SystemTrait end
struct SystemX <: SystemTrait end
struct SystemY <: SystemTrait end
struct SystemU <: SystemTrait end
struct SystemS <: SystemTrait end


################################################################################
################################### System #####################################

const XType = Union{Nothing, AbstractVector{Float64}}

#needs the C type parameter for dispatch, the rest for type stability
#must be mutable to allow y updates
mutable struct System{C <: Component, X <: XType, Y, U, D, P, S}
    ẋ::X #continuous state derivative
    x::X #continuous state
    y::Y #output
    u::U #control input
    s::D #discrete state
    t::Base.RefValue{Float64} #allows implicit propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

#default trait initializer. if the descriptor has any Component fields of
#its own, these are considered children and traits are (recursively) initialized
#from them
function init(trait::Union{SystemX, SystemY, SystemU, SystemS}, cmp::Component)
    #get those fields that are themselves Components
    children = filter(p -> isa(p.second, Component), OrderedDict(cmp))
    #build an OrderedDict with the initialized traits for each of those
    trait_dict = OrderedDict(k => init(trait, v) for (k, v) in pairs(children))
    #forward it to the OrderedDict initializers
    init(trait, trait_dict)
end

#fallback method for state vector derivative initialization
function init(::SystemẊ, cmp::Component)
    x = init(SystemX(), cmp)
    !isnothing(x) ? x |> zero : nothing
end

#initialize traits from OrderedDict
function init(::SystemTrait, dict::OrderedDict)
    filter!(p -> !isnothing(p.second), dict) #drop Nothing entries
    isempty(dict) && return nothing
    if all(v -> isa(v, AbstractVector), values(dict))
        return ComponentVector(dict)
    else
        return NamedTuple(dict)
    end
end

#shorthands (do not extend these)
init_ẋ(cmp::Component) = init(SystemẊ(), cmp)
init_x(cmp::Component) = init(SystemX(), cmp)
init_y(cmp::Component) = init(SystemY(), cmp)
init_u(cmp::Component) = init(SystemU(), cmp)
init_s(cmp::Component) = init(SystemS(), cmp)

#need a custom getproperty to address the following scenario: we have a System a
#with children b and c. if neither b and c have inputs, u = init(a, SystemU())
#will return nothing. when the System constructor for a tries to retrieve a.u.b
#and a.u.c to pass them along as inputs for subsystems b and c, it will be
#accessing fields b and c of a nothing variable.
function maybe_getproperty(input, label)
    !isnothing(input) && (label in propertynames(input)) ? getproperty(input, label) : nothing
end

function System(cmp::Component,
                ẋ = init_ẋ(cmp), x = init_x(cmp), y = init_y(cmp),
                u = init_u(cmp), s = init_s(cmp), t = Ref(0.0))

    child_names = filter(p -> (p.second isa Component), OrderedDict(cmp)) |> keys |> Tuple
    child_systems = (System(map((λ)->maybe_getproperty(λ, name), (cmp, ẋ, x, y, u, s))..., t) for name in child_names) |> Tuple
    subsystems = NamedTuple{child_names}(child_systems)

    params = NamedTuple(n=>getfield(cmp,n) for n in propertynames(cmp) if !(n in child_names))
    params = (!isempty(params) ? params : nothing)

    System{map(typeof, (cmp, x, y, u, s, params, subsystems))...}(
                         ẋ, x, y, u, s, t, params, subsystems)

end


Base.getproperty(sys::System, name::Symbol) = getproperty(sys, Val(name))
Base.setproperty!(sys::System, name::Symbol, value) = setproperty!(sys, Val(name), value)

@generated function Base.getproperty(sys::System, ::Val{S}) where {S}
    if S ∈ fieldnames(System)
        return :(getfield(sys, $(QuoteNode(S))))
    else
        return :(getfield(getfield(sys, :subsystems), $(QuoteNode(S))))
    end
end

#disallow setting any System field other than y to avoid breaking the references
#with its subsystems' fields
@generated function Base.setproperty!(sys::System, ::Val{S}, value) where {S}
    if S === :y
        return :(setfield!(sys, $(QuoteNode(S)), value))
    else
        return :(error("A System's $S cannot be reassigned, only mutated in place"))
    end
end

########################## f_ode! f_step! f_disc! ##############################

#f_ode! must update sys.ẋ, compute and reassign sys.y, then return nothing

#f_step! and f_disc! are allowed to modify a System's u, s and x. if they do so,
#they must return true, otherwise false

#caution: if a System subtype defines a f_ode!, f_disc! or f_step! method, but
#its interface is incorrectly specified, dispatch will silently revert to the
#fallback. this has potential to cause subtle bugs. a test should be used to
#confirm that the desired method is indeed being dispatched to!

#fallback method for node Systems. tries calling f_ode! on all subsystems with
#the same arguments provided to the parent System, then assembles a NamedTuple
#from the subsystems' outputs. override as required.
@inline @generated function (f_ode!(sys::System{C, X, Y, U, D, P, S}, args...)
                            where {C<:Component, X <: XType, Y, U, D, P, S})

    # Core.println("@generated f_ode! called for type $C")
    ex_main = Expr(:block)

    #call f_ode! on each subsystem
    ex_calls = Expr(:block)
    for label in fieldnames(S)
        push!(ex_calls.args,
            :(f_ode!(sys.subsystems[$(QuoteNode(label))], args...)))
    end

    ex_update_y = :(update_y!(sys))

    push!(ex_main.args, ex_calls)
    push!(ex_main.args, ex_update_y)
    push!(ex_main.args, :(return nothing))

    return ex_main

end

#fallback method for node Systems. tries calling f_step! on all subsystems with
#the same arguments provided to the parent System, then ORs their outputs.
#override as required
@inline @generated function (f_step!(sys::System{C, X, Y, U, D, P, S}, args...)
                            where {C<:Component, X <: XType, Y, U, D, P, S})

    # Core.println("Generated f_step! called for type $C")
    # Core.println()

    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))

    #call f_step! on each subsystem
    for label in fieldnames(S)
        #we need all f_step! calls executed, so we can't just chain them with ||
        push!(ex.args,
            :(x_mod = x_mod || f_step!(sys.subsystems[$(QuoteNode(label))], args...)))
    end
    return ex

end

#fallback method for node Systems. tries calling f_disc! on all subsystems with
#the same arguments provided to the parent System, then ORs their outputs.
#override as required.
@inline @generated function (f_disc!(sys::System{C, X, Y, U, D, P, S}, Δt, args...)
                            where {C<:Component, X <: XType, Y, U, D, P, S})

    # Core.println("@generated f_disc! called for $C")
    # Core.println()

    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))

    #call f_step! on each subsystem
    for label in fieldnames(S)
        #we need all f_step! calls executed, so we can't just chain them with ||
        push!(ex.args,
            :(x_mod = x_mod || f_disc!(sys.subsystems[$(QuoteNode(label))], Δt, args...)))
    end
    return ex

end

#fallback method for updating a System's output. it assembles the outputs from
#its subsystems into a NamedTuple, then assigns it to the System's y field
# @inline function (update_y!(sys::System{C, X, Y})
#     where {C<:Component, X, Y <: Nothing})
# end

@inline function (update_y!(sys::System{C, X, Y})
    where {C<:Component, X, Y})
end

@inline @generated function (update_y!(sys::System{C, X, Y})
    where {C<:Component, X, Y <: NamedTuple{L, M}} where {L, M})

    #L contains the field names of those subsystems which have outputs. retrieve
    #the y's of those subsystems and assemble them into a NamedTuple, which will
    #have the same type as Y

    #initialize main expression
    ex_main = Expr(:block)

    #build a tuple expression with subsystem outputs
    ex_ss_outputs = Expr(:tuple) #tuple expression for children's outputs
    for label in L
        push!(ex_ss_outputs.args,
            :(sys.subsystems[$(QuoteNode(label))].y))
    end

    #build a NamedTuple from the subsystem's labels and the constructed tuple
    ex_nt = Expr(:call, Expr(:curly, NamedTuple, L), ex_ss_outputs)

    #assign the result to the parent system's y
    ex_assign = Expr(:(=), :(sys.y), ex_nt)

    push!(ex_main.args, ex_assign)
    push!(ex_main.args, :(return nothing))

    return ex_main

end


################################################################################
############################## Visualization ###################################

Base.@kwdef struct SystemTreeNode
    label::Symbol = :root
    type::DataType #Component type
    function SystemTreeNode(label::Symbol, type::DataType)
        @assert (type <: Component) && (!isabstracttype(type))
        new(label, type)
    end
end

SystemTreeNode(::Type{C}) where {C<:Component} = SystemTreeNode(type = C)

function AbstractTrees.children(node::SystemTreeNode)
    return [SystemTreeNode(name, type) for (name, type) in zip(
            fieldnames(node.type), fieldtypes(node.type))
            if type <: Component]
end

function AbstractTrees.printnode(io::IO, node::SystemTreeNode)
    print(io, ":"*string(node.label)*" ($(node.type))")
end

AbstractTrees.print_tree(cmp::Type{C}; kwargs...) where {C<:Component} =
    print_tree(SystemTreeNode(cmp); kwargs...)

AbstractTrees.print_tree(::C; kwargs...) where {C<:Component} = print_tree(C; kwargs...)
AbstractTrees.print_tree(::System{C}; kwargs...) where {C} = print_tree(C; kwargs...)

end #module