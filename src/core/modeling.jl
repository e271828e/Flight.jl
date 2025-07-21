module Modeling

using ComponentArrays
using DataStructures
using AbstractTrees

using ..IODevices

export ModelDefinition, Model
export Subsampled, Scheduling, NoScheduling
export f_ode!, f_step!, f_disc!, update_output!
export @no_ode, @no_disc, @no_step, @no_updates
export @ss_ode, @ss_disc, @ss_step, @ss_dynamics


################################################################################
############################## ModelTrait #####################################

abstract type ModelTrait end

struct Ẋ <: ModelTrait end #continuous state derivative
struct X <: ModelTrait end #continuous state
struct Y <: ModelTrait end #output
struct U <: ModelTrait end #input
struct S <: ModelTrait end #discrete state


################################################################################
############################## ModelDefinition ################################

abstract type ModelDefinition end

#################### Default ModelTrait Constructors #######################

(::Union{Type{U}, Type{S}})(::ModelDefinition) = nothing

function (Trait::Union{Type{X}, Type{Y}})(md::D) where {D <: ModelDefinition}

    children_names = filter(fieldnames(D)) do name
        getfield(md, name) isa ModelDefinition
    end
    children_fields = map(λ -> getfield(md, λ), children_names)
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

function Ẋ(md::ModelDefinition)
    x = X(md)
    return (!isnothing(x) ? (x |> zero) : nothing)
end

################################################################################
################################# Subsampled ###################################

struct Subsampled{D <: ModelDefinition} <: ModelDefinition
    md::D
    K::Int #parent-relative discrete sampling period multiplier
end

Subsampled(md::D) where {D <: ModelDefinition} = Subsampled{D}(md, 1)

(Trait::Union{Type{U}, Type{S}})(ss::Subsampled) = Trait(ss.md)
(Trait::Union{Type{X}, Type{Y}})(ss::Subsampled) = Trait(ss.md)
(Trait::Type{Ẋ})(ss::Subsampled) = Ẋ(ss.md)

################################################################################
################################### Model #####################################

#type parameter D is needed for dispatching, the rest for type stability.
#Model's mutability only required for y updates; reassigning any other field is
#disallowed to avoid breaking the references in the submodel hierarchy

#storing t, n and Δt_root as RefValues allows implicit propagation of updates
#down the submodel hierarchy

mutable struct Model{D <: ModelDefinition, Y, U, X, S, C, B}
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
    const submodels::B
end

function Model(md::D,
                y = Y(md), u = U(md), ẋ = Ẋ(md), x = X(md), s = S(md),
                N = 1, Δt_root = Ref(1.0), t = Ref(0.0), n = Ref(0)) where {D <: ModelDefinition}

    if !isbits(y)
        @warn "The output defined for $D is not an isbits type.
        For performance reasons, it is strongly recommended to use a concrete,
        immutable, stack-allocated type as Model output"
    end

    sd_fieldnames = fieldnames(D)

    foreach(sd_fieldnames) do name
        name ∉ (fieldnames(Model)..., :Δt) || error(
            "Identifier $name cannot be used in a ModelDefinition field")
    end

    #construct submodels from md's ModelDefinition fields
    children_names = filter(sd_fieldnames) do name
        getfield(md, name) isa ModelDefinition
    end

    children_models = map(children_names) do child_name

        child_definition = getproperty(md, child_name)

        child_traits = map((y, u, ẋ, x, s), (Y, U, Ẋ, X, S)) do parent_trait, Trait
            if child_name in propertynames(parent_trait)
                getproperty(parent_trait, child_name)
            else
                Trait(child_definition)
            end
        end

        Model(child_definition, child_traits..., N, Δt_root, t, n)

    end

    submodels = NamedTuple{children_names}(children_models)

    #the remaining md fields are stored as constants
    constants = NamedTuple(n=>getfield(md, n) for n in sd_fieldnames if !(n in children_names))
    constants = (!isempty(constants) ? constants : nothing)

    mdl = Model{map(typeof, (md, y, u, x, s, constants, submodels))...}(
                    y, u, ẋ, x, s, N, Δt_root, t, n, constants, submodels)

    init!(mdl)

    return mdl

end

function Model(ss::Subsampled,
                y = Y(ss), u = U(ss), ẋ = Ẋ(ss), x = X(ss), s = S(ss),
                N = 1, Δt_root = Ref(1.0), t = Ref(0.0), n = Ref(0))
    Model(ss.md, y, u, ẋ, x, s, N * ss.K, Δt_root, t, n)
end

################################################################################

Base.getproperty(mdl::Model, name::Symbol) = getproperty(mdl, Val(name))

function Base.propertynames(mdl::M) where {M <: Model}
    (fieldnames(S)...,
    propertynames(mdl.constants)...,
    propertynames(mdl.submodels)...,
    :Δt,
    )
end

#no clashes can occur between constant and submodel names. thus, for
#convenience we allow direct access to both
@generated function Base.getproperty(mdl::Model{D, Y, U, X, S, C, B}, ::Val{P}
        ) where {D, Y, U, X, S, C, B, P}
    if P ∈ fieldnames(Model)
        return :(getfield(mdl, $(QuoteNode(P))))
    elseif P ∈ fieldnames(B)
        return :(getfield(getfield(mdl, :submodels), $(QuoteNode(P))))
    elseif P ∈ fieldnames(C)
        return :(getfield(getfield(mdl, :constants), $(QuoteNode(P))))
    elseif P === :Δt #actual discrete sampling period
        return :(getfield(mdl, :Δt_root)[] * getfield(mdl, :N))
    else
        return :(error("Failed to retrieve property $P from Model"))
    end
end


################################################################################
########################## Model Update Methods ###############################

#note: for a root Model, these methods must be implemented without additional
#arguments

abstract type MaybeScheduling end
struct Scheduling <: MaybeScheduling end
struct NoScheduling <: MaybeScheduling end


#fallback for the single-argument call (used in constructor)
init!(::Model) = nothing

#continuous dynamics, to be extended by Models
function f_ode!(mdl::Model, args...)
    MethodError(f_ode!, (mdl, args...)) |> throw
end

#post-step update, to be extended by Models
function f_step!(mdl::Model, args...)
    MethodError(f_step!, (mdl, args...)) |> throw
end

#unscheduled discrete update, to be extended by Models
function f_disc!(sch::NoScheduling, mdl::Model, args...)
    MethodError(f_disc!, (sch, mdl, args...)) |> throw
end

#scheduled discrete update, to be called (not extended!) by Models
@inline f_disc!(mdl::Model, args...) = f_disc!(Scheduling(), mdl, args...)

function f_disc!(::Scheduling, mdl::Model, args...)
    (mdl.n[] % mdl.N == 0) && f_disc!(NoScheduling(), mdl, args...)
end


########################## Convenience Methods #################################

#output update fallback, may be used by node Models with NamedTuple output (the
#default for node Models)
@inline function (update_output!(mdl::Model{D, Y})
    where {D <: ModelDefinition, Y <: NamedTuple{L, M}} where {L, M})

    ys = map(id -> getfield(getfield(getfield(mdl, :submodels), id), :y), L)
    mdl.y = NamedTuple{L}(ys)
    nothing
end

#any other output type needs custom implementation
@inline function update_output!(mdl::Model{D}) where {D}
    if !isempty(getfield(mdl, :submodels))
        error("An update_output! method must be explicitly implemented for node "*
        "Models with an output type other than NamedTuple; $D doesn't have one")
    end
end

update_output!(::Model{<:ModelDefinition, Nothing}) = nothing


########################## Convenience Macros ##################################

#no continuous dynamics
macro no_ode(md)
    esc(:(Modeling.f_ode!(::Model{<:($md)}, args...) = nothing))
end

#no discrete dynamics
macro no_disc(md)
    esc(:(Modeling.f_disc!(::NoScheduling, ::Model{<:($md)}, args...) = nothing))
end

#no post-step update
macro no_step(md)
    esc(:(Modeling.f_step!(::Model{<:($md)}, args...) = nothing))
end

#no dynamics at all
macro no_updates(md)
    esc(quote @no_ode $md; @no_step $md; @no_disc $md end)
end

#recursive fallbacks
macro ss_ode(md)
    esc(quote
        @inline function Modeling.f_ode!(mdl::Model{<:($md)}, args...)
            for ss in mdl.submodels
                f_ode!(ss, args...)
            end
            update_output!(mdl)
            return nothing
        end
    end)
end

macro ss_disc(md)
    esc(quote
        @inline function Modeling.f_disc!(::NoScheduling, mdl::Model{<:($md)}, args...)
            for ss in mdl.submodels
                f_disc!(ss, args...)
            end
            update_output!(mdl)
            return nothing
        end
    end)
end

macro ss_step(md)
    esc(quote
        @inline function Modeling.f_step!(mdl::Model{<:($md)}, args...)
            for ss in mdl.submodels
                f_step!(ss, args...)
            end
            update_output!(mdl) #output update not usually required here
            return nothing
        end
    end)
end

macro ss_dynamics(md)
    esc(quote @ss_ode $md; @ss_step $md; @ss_disc $md end)
end


################################################################################
############################### Inspection #####################################

#prevent monstrous type signatures from flooding the REPL
function Base.show(io::IO, ::MIME"text/plain", x::Union{Type{<:Model}, Type{<:ModelDefinition}})
    str = sprint(show, x)
    maxlen = 100
    length(str) > maxlen ? print(io, first(str, maxlen), "...") : print(io, str)
end

function Base.show(io::IO, ::MIME"text/plain", x::Union{Model, ModelDefinition})
    str = sprint(show, x)
    maxlen = 200
    length(str) > maxlen ? print(io, first(str, maxlen), "...") : print(io, str)
end

#enable hierarchy inspection with AbstractTrees.print_tree
AbstractTrees.children(node::Model) = node.submodels
AbstractTrees.printnode(io::IO, ::Model{D}) where {D} = show(io, MIME"text/plain"(), D)

end #module