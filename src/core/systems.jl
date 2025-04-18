module Systems

using ComponentArrays
using DataStructures
using AbstractTrees

using ..IODevices

export SystemDefinition, System
export Subsampled, Scheduling, NoScheduling
export f_ode!, f_step!, f_disc!, update_output!
export @no_ode, @no_disc, @no_step, @no_dynamics
export @ss_ode, @ss_disc, @ss_step, @ss_dynamics


################################################################################
############################## SystemTrait #####################################

abstract type SystemTrait end

struct Ẋ <: SystemTrait end #continuous state derivative
struct X <: SystemTrait end #continuous state
struct Y <: SystemTrait end #output
struct U <: SystemTrait end #input
struct S <: SystemTrait end #discrete state


################################################################################
############################## SystemDefinition ################################

abstract type SystemDefinition end

#################### Default SystemTrait Constructors #######################

(::Union{Type{U}, Type{S}})(::SystemDefinition) = nothing

function (Trait::Union{Type{X}, Type{Y}})(sd::D) where {D <: SystemDefinition}

    children_names = filter(fieldnames(D)) do name
        getfield(sd, name) isa SystemDefinition
    end
    children_fields = map(λ -> getfield(sd, λ), children_names)
    children_traits = map(child-> Trait(child), children_fields)
    nt = NamedTuple{children_names}(children_traits)

    nonempty_names = filter(k -> !isnothing(getproperty(nt, k)), keys(nt))
    nonempty_traits = map(λ -> getproperty(nt, λ), nonempty_names)
    filtered_nt = NamedTuple{nonempty_names}(nonempty_traits)

    if isempty(filtered_nt)
        return nothing
    elseif all(v -> isa(v, AbstractVector), values(filtered_nt)) #x and ẋ
        return ComponentVector(filtered_nt)
    else #y
        return filtered_nt
    end
end

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

Subsampled(sd::D) where {D <: SystemDefinition} = Subsampled{D}(sd, 1)

(Trait::Union{Type{U}, Type{S}})(ss::Subsampled) = Trait(ss.sd)
(Trait::Union{Type{X}, Type{Y}})(ss::Subsampled) = Trait(ss.sd)
(Trait::Type{Ẋ})(ss::Subsampled) = Ẋ(ss.sd)

################################################################################
################################### System #####################################

#type parameter D is needed for dispatching, the rest for type stability.
#System's mutability only required for y updates; reassigning any other field is
#disallowed to avoid breaking the references in the subsystem hierarchy

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

function System(sd::D,
                y = Y(sd), u = U(sd), ẋ = Ẋ(sd), x = X(sd), s = S(sd),
                N = 1, Δt_root = Ref(1.0), t = Ref(0.0), n = Ref(0)) where {D <: SystemDefinition}

    if !isbits(y)
        @warn "The output defined for $D is not an isbits type.
        For performance reasons, it is highly advisable to use a concrete,
        immutable, stack-allocated type as System output"
    end

    sd_fieldnames = fieldnames(D)

    foreach(sd_fieldnames) do name
        name ∉ (fieldnames(System)..., :Δt) || error(
            "Identifier $name cannot be used in a SystemDefinition field")
    end

    #construct subsystems from sd's SystemDefinition fields
    children_names = filter(sd_fieldnames) do name
        getfield(sd, name) isa SystemDefinition
    end

    children_systems = map(children_names) do child_name

        child_definition = getproperty(sd, child_name)

        child_traits = map((y, u, ẋ, x, s), (Y, U, Ẋ, X, S)) do parent_trait, Trait
            if child_name in propertynames(parent_trait)
                getproperty(parent_trait, child_name)
            else
                Trait(child_definition)
            end
        end

        System(child_definition, child_traits..., N, Δt_root, t, n)

    end

    subsystems = NamedTuple{children_names}(children_systems)

    #the remaining sd fields are stored as constants
    constants = NamedTuple(n=>getfield(sd, n) for n in sd_fieldnames if !(n in children_names))
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

################################################################################

Base.getproperty(sys::System, name::Symbol) = getproperty(sys, Val(name))

function Base.propertynames(sys::S) where {S <: System}
    (fieldnames(S)...,
    propertynames(sys.constants)...,
    propertynames(sys.subsystems)...,
    :Δt,
    )
end

#no clashes can occur between constant and subsystem names. thus, for
#convenience we allow direct access to both
@generated function Base.getproperty(sys::System{D, Y, U, X, S, C, B}, ::Val{P}
        ) where {D, Y, U, X, S, C, B, P}
    if P ∈ fieldnames(System)
        return :(getfield(sys, $(QuoteNode(P))))
    elseif P ∈ fieldnames(B)
        return :(getfield(getfield(sys, :subsystems), $(QuoteNode(P))))
    elseif P ∈ fieldnames(C)
        return :(getfield(getfield(sys, :constants), $(QuoteNode(P))))
    elseif P === :Δt #actual discrete sampling period
        return :(getfield(sys, :Δt_root)[] * getfield(sys, :N))
    else
        return :(error("Failed to retrieve property $P from System"))
    end
end


################################################################################
########################## System Update Methods ###############################

#note: for a root System, these methods must be implemented without additional
#arguments

abstract type MaybeScheduling end
struct Scheduling <: MaybeScheduling end
struct NoScheduling <: MaybeScheduling end


init!(::System, args...; kwargs...) = nothing

function reset!(sys::System)
    foreach(sys.subsystems) do ss
        Systems.reset!(ss)
    end
end
#continuous dynamics, to be extended by Systems
function f_ode!(sys::System, args...)
    MethodError(f_ode!, (sys, args...)) |> throw
end

#post-step update, to be extended by Systems
function f_step!(sys::System, args...)
    MethodError(f_step!, (sys, args...)) |> throw
end

#unscheduled discrete update, to be extended by Systems
function f_disc!(sch::NoScheduling, sys::System, args...)
    MethodError(f_disc!, (sch, sys, args...)) |> throw
end

#scheduled discrete update, to be called (not extended!) by Systems
@inline f_disc!(sys::System, args...) = f_disc!(Scheduling(), sys, args...)

function f_disc!(::Scheduling, sys::System, args...)
    (sys.n[] % sys.N == 0) && f_disc!(NoScheduling(), sys, args...)
end


########################## Convenience Methods #################################

#output update fallback, may be used by node Systems with NamedTuple output (the
#default for node Systems)
@inline function (update_output!(sys::System{D, Y})
    where {D <: SystemDefinition, Y <: NamedTuple{L, M}} where {L, M})

    ys = map(id -> getfield(getfield(getfield(sys, :subsystems), id), :y), L)
    sys.y = NamedTuple{L}(ys)
end

#any other output type needs custom implementation
@inline function update_output!(sys::System{D}) where {D}
    if !isempty(getfield(sys, :subsystems))
        error("An update_output! method must be explicitly implemented for node "*
        "Systems with an output type other than NamedTuple; $D doesn't have one")
    end
end

update_output!(::System{<:SystemDefinition, Nothing}) = nothing


########################## Convenience Macros ##################################

#no continuous dynamics
macro no_ode(sd)
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

#no dynamics at all
macro no_dynamics(sd)
    esc(quote @no_ode $sd; @no_step $sd; @no_disc $sd end)
end

#recursive fallbacks
macro ss_ode(sd)
    esc(quote
        @inline function Systems.f_ode!(sys::System{<:($sd)}, args...)
            for ss in sys.subsystems
                f_ode!(ss, args...)
            end
            update_output!(sys)
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
            update_output!(sys)
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
            # update_output!(sys) #output update not essential here
            return nothing
        end
    end)
end

macro ss_dynamics(sd)
    esc(quote @ss_ode $sd; @ss_step $sd; @ss_disc $sd end)
end


################################################################################
############################### Inspection #####################################

#prevent monstrous (nested, parameterized) types from flooding the REPL
function Base.show(io::IO, ::MIME"text/plain", x::Union{Type{<:System}, Type{<:SystemDefinition}})
    str = sprint(show, x)
    maxlen = 100
    length(str) > maxlen ? print(io, first(str, maxlen), "...") : print(io, str)
end

function Base.show(io::IO, ::MIME"text/plain", x::Union{System, SystemDefinition})
    str = sprint(show, x)
    maxlen = 200
    length(str) > maxlen ? print(io, first(str, maxlen), "...") : print(io, str)
end

#inspect hierarchy with AbstractTrees.print_tree
AbstractTrees.children(node::System) = node.subsystems
AbstractTrees.printnode(io::IO, ::System{D}) where {D} = show(io, MIME"text/plain"(), D)

end #module