module Modeling

using ComponentArrays
using DataStructures
using AbstractTrees

using ..IODevices

export ModelDefinition, Model
export Subsampled, Scheduling, NoScheduling
export f_ode!, f_step!, f_periodic!, f_output!, init!
export @no_ode, @no_step, @no_periodic, @no_updates
export @sm_ode, @sm_step, @sm_periodic, @sm_updates


################################################################################
############################## ModelDescriptor #################################

abstract type ModelDescriptor end

struct Ẋ <: ModelDescriptor end #continuous state derivative
struct X <: ModelDescriptor end #continuous state
struct Y <: ModelDescriptor end #output
struct U <: ModelDescriptor end #input
struct S <: ModelDescriptor end #discrete state


################################################################################
############################## ModelDefinition #################################

abstract type ModelDefinition end

#################### Default ModelDescriptor Constructors ######################

function get_children_properties(T::Type{<:ModelDescriptor}, md::D) where {D <: ModelDefinition}

    children_names = filter(fieldnames(D)) do name
        getfield(md, name) isa ModelDefinition
    end
    children_definitions = map(λ -> getfield(md, λ), children_names)
    children_properties = map(child-> T(child), children_definitions)
    nt = NamedTuple{children_names}(children_properties)

    something_names = filter(k -> !isnothing(getproperty(nt, k)), keys(nt))
    something_properties = map(λ -> getproperty(nt, λ), something_names)
    return NamedTuple{something_names}(something_properties)

end

(::Union{Type{U}, Type{S}})(::ModelDefinition) = nothing

function (T::Type{X})(md::D) where {D <: ModelDefinition}
    filtered_nt = get_children_properties(T, md)
    return (isempty(filtered_nt) ? nothing : ComponentVector(filtered_nt))
end

function (T::Type{Y})(md::D) where {D <: ModelDefinition}
    filtered_nt = get_children_properties(T, md)
    return (isempty(filtered_nt) ? nothing : filtered_nt)
end

function Ẋ(md::ModelDefinition)
    x = X(md)
    return (!isnothing(x) ? (x |> zero) : nothing)
end

################################################################################
################################# Subsampled ###################################

struct Subsampled{D <: ModelDefinition} <: ModelDefinition
    md::D
    K::Int #parent-relative sampling period multiplier
end

Subsampled(md::D) where {D <: ModelDefinition} = Subsampled{D}(md, 1)

(T::Union{Type{U}, Type{S}})(ss::Subsampled) = T(ss.md)
(::Type{X})(ss::Subsampled) = X(ss.md)
(::Type{Ẋ})(ss::Subsampled) = Ẋ(ss.md)
(::Type{Y})(ss::Subsampled) = Y(ss.md)

################################################################################
################################### Model #####################################

# - type parameter D is needed for dispatching, the rest for type stability.
# - Model's mutability only required for y updates; reassigning any other field
# is disallowed to avoid breaking the references in the submodel hierarchy.
# - t, _n and _Δt_root are stored as RefValues to allow implicit propagation of
# updates down the submodel hierarchy

mutable struct Model{D <: ModelDefinition, Y, U, X, S, P, B}
    y::Y #output
    const u::U #input
    const ẋ::X #continuous state derivative
    const x::X #continuous state
    const s::S #discrete state
    const t::Base.RefValue{Float64} #simulation time
    const parameters::P
    const submodels::B
    const _Δt_root::Base.RefValue{Float64} #root model sampling period
    const _n::Base.RefValue{Int} #simulation periodic update counter
    const _N::Int #root model-relative sampling period multipler
end

function Model(md::D,
                y = Y(md), u = U(md), ẋ = Ẋ(md), x = X(md), s = S(md), t = Ref(0.0),
                _Δt_root = Ref(1.0), _n = Ref(0), _N = 1) where {D <: ModelDefinition}

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

        child_properties = map((y, u, ẋ, x, s), (Y, U, Ẋ, X, S)) do parent_property, T
            if child_name in propertynames(parent_property)
                getproperty(parent_property, child_name)
            else
                T(child_definition)
            end
        end

        Model(child_definition, child_properties..., t, _Δt_root, _n, _N)

    end

    submodels = NamedTuple{children_names}(children_models)

    #the remaining md fields are stored as parameters
    parameters = NamedTuple(c=>getfield(md, c) for c in sd_fieldnames if !(c in children_names))

    mdl = Model{map(typeof, (md, y, u, x, s, parameters, submodels))...}(
                    y, u, ẋ, x, s, t, parameters, submodels, _Δt_root, _n, _N)

    init!(mdl)

    return mdl

end

function Model(ss::Subsampled,
                y = Y(ss), u = U(ss), ẋ = Ẋ(ss), x = X(ss), s = S(ss), t = Ref(0.0),
                _Δt_root = Ref(1.0), _n = Ref(0), _N = 1)
    Model(ss.md, y, u, ẋ, x, s, t, _Δt_root, _n, _N * ss.K)
end

################################################################################

function Base.propertynames(mdl::Model)
    parameter_names = isnothing(mdl.parameters) ? tuple() : keys(mdl.parameters)
    submodel_names = isnothing(mdl.submodels) ? tuple() : keys(mdl.submodels)
    (:x, :ẋ, :s, :u, :y, :t, :Δt, :parameters, parameter_names..., :submodels, submodel_names...)
end

Base.getproperty(mdl::Model, name::Symbol) = getproperty(mdl, Val(name))

#no clashes can occur between constant and submodel names. thus, for
#convenience we allow direct access to both
@generated function Base.getproperty(mdl::Model{D, Y, U, X, S, P, B}, ::Val{Property}
        ) where {D, Y, U, X, S, P, B, Property}
    if Property ∈ fieldnames(Model)
        return :(getfield(mdl, $(QuoteNode(Property))))
    elseif Property ∈ fieldnames(B)
        return :(getfield(getfield(mdl, :submodels), $(QuoteNode(Property))))
    elseif Property ∈ fieldnames(P)
        return :(getfield(getfield(mdl, :parameters), $(QuoteNode(Property))))
    elseif Property === :Δt #actual sampling period
        return :(getfield(mdl, :_Δt_root)[] * getfield(mdl, :_N))
    else
        return :(error("Failed to retrieve property $Property from Model"))
    end
end


################################################################################
########################## Model Update Methods ###############################

#notes:
#-a root Model these must implement these without additional arguments
#-all method implementations MUST return nothing to avoid type instability

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

#unscheduled periodic update, to be extended by Models
function f_periodic!(sch::NoScheduling, mdl::Model, args...)
    MethodError(f_periodic!, (sch, mdl, args...)) |> throw
end

#scheduled periodic update, to be called (not extended!) by Models
@inline function f_periodic!(mdl::Model, args...)
    f_periodic!(Scheduling(), mdl, args...)
    return nothing
end

#scheduled periodic update
function f_periodic!(::Scheduling, mdl::Model, args...)
    (mdl._n[] % mdl._N == 0) && f_periodic!(NoScheduling(), mdl, args...)
    return nothing
end

#generic output update fallback
@inline function (f_output!(mdl::Model{D, Y}) where {D <: ModelDefinition, Y})

    submodels = mdl.submodels
    ys = map(id -> getfield(getfield(submodels, id), :y), keys(submodels))
    mdl.y = Y(ys...)
    nothing
end

#output update fallback for node Models with NamedTuple output (the default)
@inline function (f_output!(mdl::Model{D, Y})
    where {D <: ModelDefinition, Y <: NamedTuple{L, M}} where {L, M})

    ys = map(id -> getfield(getfield(getfield(mdl, :submodels), id), :y), L)
    mdl.y = NamedTuple{L}(ys)
    nothing
end

#output update fallback for node Models with Nothing output
f_output!(::Model{<:ModelDefinition, Nothing}) = nothing


########################## Convenience Macros ##################################

#no continuous dynamics
macro no_ode(md)
    esc(:(Modeling.f_ode!(::Model{<:($md)}, args...) = nothing))
end

#no post-step update
macro no_step(md)
    esc(:(Modeling.f_step!(::Model{<:($md)}, args...) = nothing))
end

#no periodic dynamics
macro no_periodic(md)
    esc(:(Modeling.f_periodic!(::NoScheduling, ::Model{<:($md)}, args...) = nothing))
end

#no dynamics whatsoever
macro no_updates(md)
    esc(quote @no_ode $md; @no_step $md; @no_periodic $md end)
end

#recursive fallbacks
macro sm_ode(md)
    esc(quote
        @inline function Modeling.f_ode!(mdl::Model{<:($md)}, args...)
            for ss in mdl.submodels
                f_ode!(ss, args...)
            end
            f_output!(mdl)
            return nothing
        end
    end)
end

macro sm_step(md)
    esc(quote
        @inline function Modeling.f_step!(mdl::Model{<:($md)}, args...)
            for ss in mdl.submodels
                f_step!(ss, args...)
            end
            f_output!(mdl) #output update not usually required, consider removal
            return nothing
        end
    end)
end

macro sm_periodic(md)
    esc(quote
        @inline function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:($md)}, args...)
            for ss in mdl.submodels
                f_periodic!(ss, args...)
            end
            f_output!(mdl)
            return nothing
        end
    end)
end


macro sm_updates(md)
    esc(quote @sm_ode $md; @sm_step $md; @sm_periodic $md end)
end


################################################################################
############################### Inspection #####################################

#we need to prevent monstrous type signatures from flooding the REPL

#only print the truncated type
function Base.show(io::IO, ::D) where {D <: ModelDefinition}
    maxlen = 100
    md_str = sprint(show, D)
    md_str = (length(md_str) < maxlen ? md_str : first(md_str, maxlen) * "...")
    write(io, md_str * "(...)")
end

#only print the truncated ModelDefinition type parameter
function Base.show(io::IO, ::Model{D}) where {D <: ModelDefinition}
    maxlen = 100 #maximum length for ModelDefinition type parameter
    md_str = sprint(show, D)
    md_str = (length(md_str) < maxlen ? md_str : first(md_str, maxlen) * "...")
    write(io, "Model{" * md_str * "}(...)")
end

#enable hierarchy inspection with AbstractTrees.print_tree
AbstractTrees.children(node::Model) = node.submodels

#print the Model hierarchy showing the Model at each node
function AbstractTrees.print_tree(mdl::Model, s::Symbol; kw...)
    print_tree(stdout, mdl, s; kw...)
end

#print the Model hierarchy showing property s at each node
function AbstractTrees.print_tree(io::IO, mdl::Model, s::Symbol; kw...)
    print_tree((io::IO, mdl::Model)->print(io, getproperty(mdl, s)), io, mdl; kw...)
end

end #module