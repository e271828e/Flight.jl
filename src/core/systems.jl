module Systems

using ComponentArrays
using DataStructures
using AbstractTrees

using ..GUI
using ..IODevices

export SystemDefinition, SystemTrait, System
export f_ode!, f_step!, f_disc!, assemble_y!


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
    else #u, s
        return filtered_nt
    end
end

#y must always be a NamedTuple, even if all subsystem's y are StaticArrays;
#otherwise assemble_y! does not work
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

#D type parameter is needed for dispatching, the rest for type stability.
#mutability required for y updates; reassigning any other field is disallowed to
#avoid breaking the references with the subsystem hierarchy

#storing t as a RefValue allows implicit propagation of t updates down the
#subsystem hierarchy

mutable struct System{D <: SystemDefinition, Y, U, X, S, C, B}
    y::Y #output
    const u::U #input
    const ẋ::X #continuous state derivative
    const x::X #continuous state
    const s::S #discrete state
    const t::Base.RefValue{Float64} #simulation time
    const constants::C
    const subsystems::B
end


function System(sd::SystemDefinition,
                y = Y(sd), u = U(sd), ẋ = Ẋ(sd), x = X(sd), s = S(sd), t = Ref(0.0))

    #construct subsystems from those fields of sd which are themselves
    #SystemDefinitions
    children_names = filter(propertynames(sd)) do name
        getproperty(sd, name) isa SystemDefinition
    end

    children_systems = map(children_names) do child_name

        child_definition = getproperty(sd, child_name)

        child_properties = map((y, u, ẋ, x, s), (Y, U, Ẋ, X, S)) do parent, trait
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

    sys = System{map(typeof, (sd, y, u, x, s, constants, subsystems))...}(
                    y, u, ẋ, x, s, t, constants, subsystems)

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

@generated function Base.getproperty(sys::System, ::Val{S}) where {S}
    if S ∈ fieldnames(System)
        return :(getfield(sys, $(QuoteNode(S))))
    else
        return :(getfield(getfield(sys, :subsystems), $(QuoteNode(S))))
    end
end


################################################################################
######################## f_ode! f_step! f_disc! ################################

#f_ode! must update sys.ẋ, compute and reassign sys.y, then return nothing

#f_step! and f_disc! are allowed to modify a System's u, s and x. if they do
#modify x, they must return true, otherwise false

#caution: if a System subtype defines a f_ode!, f_disc! or f_step! method, but
#its interface is incorrectly specified, dispatch will silently revert to the
#fallback. this has potential to cause subtle bugs. a test using @which should
#be used to confirm that the desired method is indeed being dispatched to!

#note: for some reason, map and foreach tend to allocate as the recursion
#reaches further down the subsystem hierarchy


############################ Fallback methods ##################################

#tries calling f_ode! on all subsystems with the same arguments provided to the
#parent System, then assembles a NamedTuple from the subsystems' outputs.
#override as required.
@inline function f_ode!(sys::System, args...)
    # sys.subsystems |> keys |> println
    for ss in sys.subsystems
        f_ode!(ss, args...)
    end
    assemble_y!(sys)
    return nothing

end

#tries calling f_disc! on all subsystems with the same arguments provided to the
#parent System. also updates y, since f_disc! is where discrete Systems are
#expected to update their output. override as required.
@inline function f_disc!(sys::System, Δt, args...)
    # sys.subsystems |> keys |> println
    for ss in sys.subsystems
        f_disc!(ss, Δt, args...)
    end
    assemble_y!(sys)
    return nothing

end

#tries calling f_step! on all subsystems with the same arguments provided to the
#parent System. does NOT update y. override as required
@inline function f_step!(sys::System, args...)
    for ss in sys.subsystems
        f_step!(ss, args...)
    end
    return nothing

end

#fallback method for node Systems with NamedTuple output
@inline function (assemble_y!(sys::System{D, Y})
    where {D <: SystemDefinition, Y <: NamedTuple{L, M}} where {L, M})

    ys = map(id -> getproperty(sys.subsystems[id], :y), L)
    sys.y = NamedTuple{L}(ys)

end

@inline function assemble_y!(sys::System{D}) where {D}
    if !isempty(sys.subsystems)
        error("An assemble_y! method must be explicitly implemented for node "*
        "Systems with an output type other than NamedTuple, $D doesn't have one")
    end
end


################################################################################
#################################### I/O #######################################

#to add support for an InputDevice, the System should extend this function for
#the data type produced by that InputDevice. additional customization is
#possible by dispatching on a specific IOMapping subtype
function assign_input!(sys::System, data::Any, mapping::IOMapping)
    MethodError(assign_input!, (sys, data, mapping)) |> throw
end

assign_input!(::System, ::Nothing, ::IOMapping) = nothing

#to add support for an OutputDevice, the System should extend this function for
#the Type expected by that OutputDevice. additional customization is possible
#by dispatching on a specific IOMapping subtype
function extract_output(sys::System, type::Type{<:Any}, mapping::IOMapping)
    MethodError(extract_output, (sys, type, mapping)) |> throw
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

SystemTreeNode(::Type{D}) where {D <: SystemDefinition} = SystemTreeNode(type = D)

function AbstractTrees.children(node::SystemTreeNode)
    return [SystemTreeNode(name, type) for (name, type) in zip(
            fieldnames(node.type), fieldtypes(node.type))
            if type <: SystemDefinition]
end

function AbstractTrees.printnode(io::IO, node::SystemTreeNode)
    print(io, ":"*string(node.label)*" ($(node.type))")
end

AbstractTrees.print_tree(sd::Type{D}; kwargs...) where {D <: SystemDefinition} =
    print_tree(SystemTreeNode(sd); kwargs...)

AbstractTrees.print_tree(::D; kwargs...) where {D <: SystemDefinition} = print_tree(D; kwargs...)
AbstractTrees.print_tree(::System{D}; kwargs...) where {D} = print_tree(D; kwargs...)

end #module