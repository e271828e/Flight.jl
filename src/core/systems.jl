module Systems

using ComponentArrays
using DataStructures
using AbstractTrees

using ..GUI

export SystemDefinition, SystemTrait, System
export f_ode!, f_step!, f_disc!, update_y!


################################################################################
############################## SystemDefinition ################################

abstract type SystemDefinition end

function Base.NamedTuple(sd::SystemDefinition)
    ks = propertynames(sd)
    vs = map(λ -> getproperty(sd, λ), ks)
    NamedTuple{ks}(vs)
end

function delete_nothings(nt::NamedTuple)
    valid_keys = filter(k -> !isnothing(getproperty(nt, k)), keys(nt))
    valid_values = map(λ -> getproperty(nt, λ), valid_keys)
    NamedTuple{valid_keys}(valid_values)
end

################################################################################
############################## SystemTrait #####################################

abstract type SystemTrait end

struct Ẋ <: SystemTrait end
struct X <: SystemTrait end
struct Y <: SystemTrait end
struct U <: SystemTrait end
struct S <: SystemTrait end

#default trait constructors
(::Union{Type{U}, Type{S}})(::SystemDefinition) = nothing

function (trait::Union{Type{Ẋ}, Type{X}, Type{Y}})(sd::SystemDefinition)
    #get those fields that are themselves SystemDefinitions
    children_keys = filter(propertynames(sd)) do name
        getproperty(sd, name) isa SystemDefinition
    end
    children = map(λ -> getproperty(sd, λ), children_keys)
    children_traits = map(child-> trait(child), children)
    nt = NamedTuple{children_keys}(children_traits)
    return trait(nt)
end

#from NamedTuple
function (::Type{<:SystemTrait})(nt::NamedTuple)
    filtered_nt = delete_nothings(nt)
    isempty(filtered_nt) && return nothing
    if all(v -> isa(v, AbstractVector), values(filtered_nt)) #x and ẋ
        return ComponentVector(filtered_nt)
    else #u, s and y
        return filtered_nt
    end
end

#y must always be a NamedTuple, even if all subsystem's y are StaticArrays;
#otherwise update_y! does not work
function Y(nt::NamedTuple)
    filtered_nt = delete_nothings(nt)
    return !isempty(filtered_nt) ? filtered_nt : nothing
end

#fallback method for state vector derivative initialization
function Ẋ(sd::SystemDefinition)
    x = X(sd)
    return (!isnothing(x) ? (x |> zero) : nothing)
end


################################################################################
################################### System #####################################

const XType = Union{Nothing, AbstractVector{Float64}}

#needs the SD type parameter for dispatch, the rest for type stability
#must be mutable to allow y updates
mutable struct System{SD <: SystemDefinition, X <: XType, Y, U, S, P, B}
    ẋ::X #continuous state derivative
    x::X #continuous state
    y::Y #output
    u::U #control input
    s::S #discrete state
    t::Base.RefValue{Float64} #allows implicit propagation of t updates down the subsystem hierarchy
    constants::P
    subsystems::B
end


function System(sd::SystemDefinition,
                ẋ = Ẋ(sd), x = X(sd), y = Y(sd),
                u = U(sd), s = S(sd), t = Ref(0.0))

    #construct subsystems from those fields of sd which are themselves
    #SystemDefinitions
    children_names = filter(propertynames(sd)) do name
        getproperty(sd, name) isa SystemDefinition
    end

    children_systems = map(children_names) do child_name

        child_definition = getproperty(sd, child_name)

        child_properties = map((ẋ, x, y, u, s), (Ẋ, X, Y, U, S)) do parent, trait
            if !isnothing(parent) && (child_name in propertynames(parent))
                getproperty(parent, child_name)
            else
                trait(child_definition)
            end
        end

        System(child_definition, child_properties..., t)

    end

    subsystems = NamedTuple{children_names}(children_systems)

    #the remaining fields of the SystemDefinition are saved as parameters
    constants = NamedTuple(n=>getfield(sd, n) for n in propertynames(sd) if !(n in children_names))
    constants = (!isempty(constants) ? constants : nothing)

    sys = System{map(typeof, (sd, x, y, u, s, constants, subsystems))...}(
                    ẋ, x, y, u, s, t, constants, subsystems)

    init!(sys)

    return sys

end

init!(::System) = nothing

function reset!(sys::System)
    foreach(sys.subsystems) do ss
        Systems.reset!(ss)
    end
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

######################## f_ode! f_step! f_disc! ############################

#f_ode! must update sys.ẋ, compute and reassign sys.y, then return nothing

#f_step! and f_disc! are allowed to modify a System's u, s and x. if they do
#modify x, they must return true, otherwise false

#caution: if a System subtype defines a f_ode!, f_disc! or f_step! method, but
#its interface is incorrectly specified, dispatch will silently revert to the
#fallback. this has potential to cause subtle bugs. a test using @which should
#be used to confirm that the desired method is indeed being dispatched to!

#fallback method for node Systems. tries calling f_ode! on all subsystems with
#the same arguments provided to the parent System, then assembles a NamedTuple
#from the subsystems' outputs. override as required.
@inline function (f_ode!(sys::System{SD, X, Y, U, S, P, B}, args...)
                where {SD <: SystemDefinition, X <: XType, Y, U, S, P, B})

    map(ss-> f_ode!(ss, args...), values(sys.subsystems))
    update_y!(sys)

end

#fallback method for node Systems. tries calling f_disc! on all subsystems with
#the same arguments provided to the parent System. also updates y, since f_disc!
#is where discrete Systems are expected to update their output. override as
#required.
@inline function (f_disc!(sys::System{SD, X, Y, U, S, P, B}, Δt, args...)
                    where {SD <: SystemDefinition, X <: XType, Y, U, S, P, B})

    map(ss-> f_disc!(ss, Δt, args...), values(sys.subsystems))
    update_y!(sys)

end

#fallback method for node Systems. tries calling f_step! on all subsystems with
#the same arguments provided to the parent System. does NOT update y. override
#as required
@inline function (f_step!(sys::System{SD, X, Y, U, S, P, B}, args...)
                    where {SD <: SystemDefinition, X <: XType, Y, U, S, P, B})

    map(ss-> f_step!(ss, args...), values(sys.subsystems))

end

@inline function (update_y!(sys::System{SD, X, Y})
    where {SD <: SystemDefinition, X, Y})
    # println("For Anything")
    # println("Hi")
end

# @inline function (update_y!(sys::System{SD, X, Nothing})
#     where {SD <: SystemDefinition, X})
#     println("Ho")
# end

#fallback method for updating a System's NamedTuple output. it assembles the
#outputs from its subsystems into a NamedTuple, then assigns it to the System's
#y field
@inline function (update_y!(sys::System{SD, X, Y})
    where {SD <: SystemDefinition, X, Y <: NamedTuple{L, M}} where {L, M})

    # println("For NT")
    #the keys of NamedTuple sys.y identify those subsystems with non-null
    #outputs; retrieve their updated ys and assemble them into a NamedTuple of
    #the same type
    ys = map(id -> getproperty(sys.subsystems[id], :y), L)
    sys.y = NamedTuple{L}(ys)

end


################################################################################
############################### Inspection #####################################

@kwdef struct SystemTreeNode
    label::Symbol = :root
    type::DataType #SystemDefinition type
    function SystemTreeNode(label::Symbol, type::DataType)
        @assert (type <: SystemDefinition) && (!isabstracttype(type))
        new(label, type)
    end
end

SystemTreeNode(::Type{SD}) where {SD <: SystemDefinition} = SystemTreeNode(type = SD)

function AbstractTrees.children(node::SystemTreeNode)
    return [SystemTreeNode(name, type) for (name, type) in zip(
            fieldnames(node.type), fieldtypes(node.type))
            if type <: SystemDefinition]
end

function AbstractTrees.printnode(io::IO, node::SystemTreeNode)
    print(io, ":"*string(node.label)*" ($(node.type))")
end

AbstractTrees.print_tree(sd::Type{SD}; kwargs...) where {SD <: SystemDefinition} =
    print_tree(SystemTreeNode(sd); kwargs...)

AbstractTrees.print_tree(::SD; kwargs...) where {SD <: SystemDefinition} = print_tree(SD; kwargs...)
AbstractTrees.print_tree(::System{SD}; kwargs...) where {SD} = print_tree(SD; kwargs...)

end #module