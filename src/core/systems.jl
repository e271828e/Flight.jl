module Systems

using ComponentArrays
using DataStructures
using AbstractTrees

using ..IODevices

export SystemDefinition, SystemTrait, System
export Subsampled, Scheduling, NoScheduling
export f_ode!, f_step!, f_disc!, update_y!
export @no_cont, @no_disc, @no_step, @no_dynamics
export @ss_cont, @ss_disc, @ss_step, @ss_dynamics



################################################################################
############################## SystemTrait #####################################

abstract type SystemTrait end

struct Ẋ <: SystemTrait end
struct X <: SystemTrait end
struct Y <: SystemTrait end
struct U <: SystemTrait end
struct S <: SystemTrait end


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
#otherwise update_y! does not work
function Y(nt::NamedTuple)
    filtered_nt = delete_nothings(nt)
    return !isempty(filtered_nt) ? filtered_nt : nothing
end


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

########################## Default Trait Constructors ##########################

(::Union{Type{U}, Type{S}})(::SystemDefinition) = nothing

#returns NamedTuple
function (trait::Union{Type{X}, Type{Y}})(sd::SystemDefinition)
    children_keys = filter(propertynames(sd)) do name
        getproperty(sd, name) isa SystemDefinition
    end
    children_definitions = map(λ -> getproperty(sd, λ), children_keys)
    children_traits = map(child-> trait(child), children_definitions)
    nt = NamedTuple{children_keys}(children_traits)
    return trait(nt)
end

#fallback method for state vector derivative
function Ẋ(sd::SystemDefinition)
    x = X(sd)
    return (!isnothing(x) ? (x |> zero) : nothing)
end

################################################################################
################################# Subsampled ###################################

struct Subsampled{D <: SystemDefinition} <: SystemDefinition
    sd::D
    K::Int #parent-relative discrete sampling period multiplier
end

Subsampled(sd::SD) where {SD <: SystemDefinition} = Subsampled{SD}(sd, 1)

(trait::Union{Type{U}, Type{S}})(ss::Subsampled) = trait(ss.sd)
(trait::Union{Type{X}, Type{Y}})(ss::Subsampled) = trait(ss.sd)
(trait::Type{Ẋ})(ss::Subsampled) = Ẋ(ss.sd)

################################################################################
################################### System #####################################

#D type parameter is needed for dispatching, the rest for type stability.
#mutability only required for y updates; reassigning any other field is
#disallowed to avoid breaking the references with the subsystem hierarchy

#storing t, n and Δt_root as RefValues allows implicit propagation of updates
#down the subsystem hierarchy

mutable struct System{D <: SystemDefinition, Y, U, X, S, C, B}
    y::Y #output
    const u::U #input
    const ẋ::X #continuous state derivative
    const x::X #continuous state
    const s::S #discrete state
    const N::Int #discrete sampling period multipler
    const Δt_root::Base.RefValue{Float64} #root system discrete sampling period
    const t::Base.RefValue{Float64} #simulation time
    const n::Base.RefValue{Int} #simulation discrete iteration counter
    const constants::C
    const subsystems::B
end

function System(sd::SystemDefinition,
                y = Y(sd), u = U(sd), ẋ = Ẋ(sd), x = X(sd), s = S(sd),
                N = 1, Δt_root = Ref(1.0), t = Ref(0.0), n = Ref(0))

    #construct subsystems from SystemDefinition fields
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

        System(child_definition, child_properties..., N, Δt_root, t, n)

    end

    subsystems = NamedTuple{children_names}(children_systems)

    #the remaining fields of the SystemDefinition are saved as parameters
    constants = NamedTuple(n=>getfield(sd, n) for n in propertynames(sd) if !(n in children_names))
    constants = (!isempty(constants) ? constants : nothing)

    sys = System{map(typeof, (sd, y, u, x, s, constants, subsystems))...}(
                    y, u, ẋ, x, s, N, Δt_root, t, n, constants, subsystems)

    init!(sys)

    return sys

end

function System(ss::Subsampled,
                y = Y(ss), u = U(ss), ẋ = Ẋ(ss), x = X(ss), s = S(ss),
                N = 1, Δt_root = Ref(1.0), t = Ref(0.0), n = Ref(0))
    System(ss.sd, y, u, ẋ, x, s, N * ss.K, Δt_root, t, n)
end


init!(::System, args...) = nothing

function reset!(sys::System)
    foreach(sys.subsystems) do ss
        Systems.reset!(ss)
    end
end

Base.getproperty(sys::System, name::Symbol) = getproperty(sys, Val(name))

@generated function Base.getproperty(sys::System, ::Val{S}) where {S}
    if S ∈ fieldnames(System)
        return :(getfield(sys, $(QuoteNode(S))))
    elseif S === :Δt
        return :(getfield(sys, :Δt_root)[] * getfield(sys, :N))
    else
        return :(getfield(getfield(sys, :subsystems), $(QuoteNode(S))))
    end
end


################################################################################
########################## System Update Methods ###############################

abstract type MaybeSchedule end
struct Schedule <: MaybeSchedule end
struct NoScheduling <: MaybeSchedule end

function f_ode!(sys::System, args...)
    MethodError(f_ode!, (sys, args...)) |> throw
end

@inline f_disc!(sys::System, args...) = f_disc!(Schedule(), sys, args...)

function f_disc!(::Schedule, sys::System, args...)
    (sys.n[] % sys.N == 0) && f_disc!(NoScheduling(), sys, args...)
end

function f_disc!(sch::NoScheduling, sys::System, args...)
    MethodError(f_disc!, (sch, sys, args...)) |> throw
end

function f_step!(sys::System, args...)
    MethodError(f_step!, (sys, args...)) |> throw
end

#default for Systems with no output
update_y!(::System{<:SystemDefinition, Nothing}) = nothing

#default for Systems with NamedTuple output
@inline function (update_y!(sys::System{D, Y})
    where {D <: SystemDefinition, Y <: NamedTuple{L, M}} where {L, M})

    ys = map(id -> getproperty(sys.subsystems[id], :y), L)
    sys.y = NamedTuple{L}(ys)
end

#else, we require an explicit implementation
@inline function update_y!(sys::System{D}) where {D}
    if !isempty(sys.subsystems)
        error("An update_y! method must be explicitly implemented for node "*
        "Systems with an output type other than NamedTuple, $D doesn't have one")
    end
end


########################## Convenience Macros ##################################

#no continuous dynamics
macro no_cont(sd)
    esc(:(Systems.f_ode!(::System{<:($sd)}, args...) = nothing))
end

#no discrete dynamics
macro no_disc(sd)
    esc(:(Systems.f_disc!(::NoScheduling, ::System{<:($sd)}, args...) = nothing))
end

#no post-step update
macro no_step(sd)
    esc(:(Systems.f_step!(::System{<:($sd)}, args...) = nothing))
end

macro no_dynamics(sd)
    esc(quote @no_cont $sd; @no_step $sd; @no_disc $sd end)
end

#recursive fallbacks: apply the call to all subsystems and updates output
macro ss_cont(sd)
    esc(quote
        @inline function Systems.f_ode!(sys::System{<:($sd)}, args...)
            for ss in sys.subsystems
                f_ode!(ss, args...)
            end
            update_y!(sys)
            return nothing
        end
    end)
end

macro ss_disc(sd)
    esc(quote
        @inline function Systems.f_disc!(::NoScheduling, sys::System{<:($sd)}, args...)
            for ss in sys.subsystems
                f_disc!(ss, args...)
            end
            update_y!(sys)
            return nothing
        end
    end)
end

macro ss_step(sd)
    esc(quote
        @inline function Systems.f_step!(sys::System{<:($sd)}, args...)
            for ss in sys.subsystems
                f_step!(ss, args...)
            end
            # update_y!(sys) #output update not essential here
            return nothing
        end
    end)
end

macro ss_dynamics(sd)
    esc(quote @ss_cont $sd; @ss_step $sd; @ss_disc $sd end)
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
function extract_output(sys::System, device::OutputDevice, mapping::IOMapping)
    MethodError(extract_output, (sys, device, mapping)) |> throw
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