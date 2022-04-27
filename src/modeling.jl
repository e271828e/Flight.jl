module Modeling

using Dates
using UnPack
using ComponentArrays, StructArrays, RecursiveArrayTools
using SciMLBase: ODEProblem, u_modified!, init as init_problem
using OrdinaryDiffEq: ODEIntegrator, Tsit5
using DiffEqCallbacks: SavingCallback, DiscreteCallback, CallbackSet, SavedValues

import AbstractTrees: children, printnode, print_tree
import SciMLBase: step!, solve!, reinit!, get_proposed_dt
import DataStructures: OrderedDict

export f_cont!, f_disc!, step!
export SystemDescriptor, SystemGroupDescriptor, NullSystemDescriptor, System, Model, TimeHistory
export SystemẊ, SystemX, SystemY, SystemU, SystemD


################################################################################
############################## SystemDescriptor ################################

abstract type SystemDescriptor end #anything from which we can build a System

abstract type SystemData end

struct SystemẊ <: SystemData end
struct SystemX <: SystemData end
struct SystemY <: SystemData end
struct SystemU <: SystemData end
struct SystemD <: SystemData end

################################################################################
################################### System #####################################

#need the T type parameter for dispatch, the rest for type stability. since
#Systems are only meant to be instantiated during initialization, making
#them mutable does not hurt performance (no runtime heap allocations)
mutable struct System{T <: SystemDescriptor, X, Y, U, D, P, S}
    ẋ::X #continuous dynamics state vector derivative
    x::X #continuous dynamics state vector
    y::Y #output
    u::U #control input
    d::D #discrete dynamics state
    t::Base.RefValue{Float64} #allows implicit propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

function OrderedDict(g::SystemDescriptor)
    fields = propertynames(g)
    values = map(λ -> getproperty(g, λ), fields)
    OrderedDict(k => v for (k, v) in zip(fields, values))
end

#fallback method for state vector derivative initialization
function init(d::SystemDescriptor, ::SystemẊ)
    x = init(d, SystemX())
    !isnothing(x) ? x |> similar |> zero : nothing
end

#if the SystemDescriptor has no SystemDescriptor children, and does not override
#this method, it will initialize all its traits to nothing
function init(desc::SystemDescriptor, trait::Union{SystemX, SystemY, SystemU, SystemD})

    child_descriptors = filter(p -> isa(p.second, SystemDescriptor), OrderedDict(desc))

    dict = OrderedDict()
    for (name, child) in child_descriptors
        child_data = init(child, trait)
        !isnothing(child_data) ? dict[name] = child_data : nothing
    end

    init(dict, trait)

end

function init(dict::OrderedDict, ::Union{SystemX})
    !isempty(dict) ? ComponentVector(dict) : nothing
end

function init(dict::OrderedDict, ::Union{SystemY, SystemU, SystemD})
    !isempty(dict) ? NamedTuple{Tuple(keys(dict))}(values(dict)) : nothing
end

#suppose we have a System a with children b and c. the System constructor will
#try to retrieve views a.x.b and a.x.c and assign them as state vectors b.x and
#c.x. but if b and c have no continuous states, x = init(a, SystemX()) will
#return nothing. so the constructor will try to retrieve fields from a nothing
#variable. this function handles this scenario
function maybe_getproperty(input, label)
    !isnothing(input) && (label in propertynames(input)) ? getproperty(input, label) : nothing
end

function System(desc::SystemDescriptor,
                ẋ = init(desc, SystemẊ()), x = init(desc, SystemX()),
                y = init(desc, SystemY()), u = init(desc, SystemU()),
                d = init(desc, SystemD()), t = Ref(0.0))

    child_names = filter(p -> (p.second isa SystemDescriptor), OrderedDict(desc)) |> keys |> Tuple
    child_systems = (System(map((λ)->maybe_getproperty(λ, name), (desc, ẋ, x, y, u, d))..., t) for name in child_names) |> Tuple
    subsystems = NamedTuple{child_names}(child_systems)

    params = NamedTuple(n=>getfield(desc,n) for n in propertynames(desc) if !(n in child_names))
    params = (!isempty(params) ? params : nothing)

    System{map(typeof, (desc, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)

end

#f_disc! is free to modify a Hybrid system's discrete state, control inputs and
#continuous state. if it modifies the latter, it must return true, false
#otherwise. no fallbacks are provided for safety reasons: if the intended
#f_cont! or f_disc! implementations for the System have the wrong interface, the
#dispatch will silently revert to the fallback, which does nothing and may not
#be obvious at all.

f_cont!(sys::System, args...) = MethodError(f_cont!, (sys, args...)) |> throw
(f_disc!(sys::System, args...)::Bool) = MethodError(f_disc!, (sys, args...)) |> throw

Base.getproperty(sys::System, s::Symbol) = getproperty(sys, Val(s))

@generated function Base.getproperty(sys::System, ::Val{S}) where {S}
    if S ∈ fieldnames(System)
        return :(getfield(sys, $(QuoteNode(S))))
    else
        return :(getfield(getfield(sys, :subsystems), $(QuoteNode(S))))
    end
end


@inline function (assemble_y!(sys::System{T, X, Y})
    where {T<:SystemDescriptor, X, Y <: Nothing})
end

@inline @generated function (assemble_y!(sys::System{T, X, Y})
    where {T<:SystemDescriptor, X, Y <: NamedTuple{L, M}} where {L, M})

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

################################ NullSystem ################################

struct NullSystemDescriptor <: SystemDescriptor end

@inline f_cont!(::System{NullSystemDescriptor}, args...) = nothing
@inline (f_disc!(::System{NullSystemDescriptor}, args...)::Bool) = false

######################### SystemGroupDescriptors ############################

#abstract supertype for any SystemDescriptor grouping several children
#SystemDescriptors with common f_cont! and f_disc! interfaces. it provides
#automatically generated methods for these functions.

abstract type SystemGroupDescriptor <: SystemDescriptor end

#default implementation calls f_cont! on all Group subsystems with the same
#arguments provided to the parent System, then builds a NamedTuple with the
#subsystem outputs. can be overridden as required.
@inline @generated function (f_cont!(sys::System{T, X, Y, U, D, P, S}, args...)
    where {T<:SystemGroupDescriptor, X, Y, U, D, P, S})

    # Core.println("Generated function called")
    ex_main = Expr(:block)

    #call f_cont! on each subsystem
    ex_calls = Expr(:block)
    for label in fieldnames(S)
        push!(ex_calls.args,
            :(f_cont!(sys.subsystems[$(QuoteNode(label))], args...)))
    end

    ex_assemble_y = :(assemble_y!(sys))

    push!(ex_main.args, ex_calls)
    push!(ex_main.args, ex_assemble_y)
    push!(ex_main.args, :(return nothing))

    return ex_main

end

#default implementation calls f_disc! on all Node subsystems with the same
#arguments provided to the parent Node's System, then ORs their outputs.
#can be overridden as required
@inline @generated function (f_disc!(sys::System{T, X, Y, U, D, P, S}, args...)
    where {T<:SystemGroupDescriptor, X, Y, U, D, P, S})

    # Core.print("Generated function called")
    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))

    #call f_disc! on each subsystem
    for label in fieldnames(S)
        #we need all f_disc! calls executed, so | must be used instead of ||
        push!(ex.args,
            :(x_mod = x_mod | f_disc!(sys.subsystems[$(QuoteNode(label))], args...)))
    end
    return ex

end

################################################################################
############################## Visualization ###################################

Base.@kwdef struct SystemTreeNode
    label::Symbol = :root
    type::DataType #SystemDescriptor type
    function SystemTreeNode(label::Symbol, type::DataType)
        @assert (type <: SystemDescriptor) && (!isabstracttype(type))
        new(label, type)
    end
end

SystemTreeNode(::T) where {T<:SystemDescriptor} = SystemTreeNode(type = T)
SystemTreeNode(::System{D}) where {D} = SystemTreeNode(type = D)

function children(node::SystemTreeNode)
    return [SystemTreeNode(name, type) for (name, type) in zip(
            fieldnames(node.type), fieldtypes(node.type))
            if type <: SystemDescriptor]
end

function printnode(io::IO, node::SystemTreeNode)
    print(io, ":"*string(node.label)*" ($(node.type))")
end

print_tree(desc::SystemDescriptor) = print_tree(SystemTreeNode(desc))
print_tree(sys::System) = print_tree(SystemTreeNode(sys))


# function AbstractTrees.children(x::ComponentVector)
#     c = []
#     for k in keys(x)
#         k isa Symbol ? push!(c, getproperty(x, k)) : nothing
#     end

#     # vector = [getproperty(x, k) for k in keys(x) if keys(x) isa Tuple]
#     return c
# end

# function AbstractTrees.printnode(io::IO, x::ComponentVector)
#     nothing
# end



################################################################################
################################## Model #######################################

#in this design, the t and x fields of m.sys behave only as temporary
#storage for f_cont! and f_disc! calls, so we have no guarantees about their
#status after a certain step. the only valid sources for t and x at any
#given moment is the integrator's t and u
struct Model{S <: System, I <: ODEIntegrator, L <: SavedValues}
    sys::S
    integrator::I
    log::L

    function Model(
        sys, #System
        args_c::Tuple = (), #externally supplied arguments to f_cont!
        args_d::Tuple = (); #externally supplied arguments to f_disc!
        solver = Tsit5(), #
        t_start = 0.0,
        t_end = 10.0,
        y_saveat = Float64[],
        save_on = false,
        integrator_kwargs...)

        #save_on is set to false because we are not usually interested in the
        #naked System's state vector. everything we need should be available in
        #the System's output struct saved by the SavingCallback
        saveat_arr = (y_saveat isa Real ? (t_start:y_saveat:t_end) : y_saveat)

        params = (sys = sys, args_c = args_c, args_d = args_d)

        log = SavedValues(Float64, typeof(sys.y))

        dcb = DiscreteCallback((u, t, integrator)->true, f_dcb!)
        scb = SavingCallback(f_scb, log, saveat = saveat_arr)
        cb_set = CallbackSet(dcb, scb)

        #use the current System's x value as initial condition. need a copy,
        #otherwise a reference is used and the initial value is overwritten
        x0 = copy(sys.x)
        problem = ODEProblem{true}(f_update!, x0, (t_start, t_end), params)
        integrator = init_problem(problem, solver; callback = cb_set, save_on, integrator_kwargs...)
        new{typeof(sys), typeof(integrator), typeof(log)}(sys, integrator, log)
    end
end

Base.getproperty(mdl::Model, s::Symbol) = getproperty(mdl, Val(s))

@generated function Base.getproperty(mdl::Model, ::Val{S}) where {S}
    if S === :t
        return :(getproperty(getfield(mdl, :sys), $(QuoteNode(S)))[])
    elseif S ∈ fieldnames(System)
        return :(getproperty(getfield(mdl, :sys), $(QuoteNode(S))))
    else
        return :(getfield(mdl, $(QuoteNode(S))))
    end
end

#neither f_update!, f_dcb! nor f_scb should allocate

#integrator ODE function. basically, a wrapper around the System's continuous
#dynamics function. as a side effect, modifies the System's ẋ and y.
function f_update!(u̇, u, p, t)

    @unpack sys, args_c = p

    sys.x .= u #assign current integrator solution to System continuous state
    sys.t[] = t #ditto for time

    f_cont!(sys, args_c...) #call continuous dynamics, updates sys.ẋ and sys.y

    u̇ .= sys.ẋ #update the integrator's derivative

    return nothing

end

#DiscreteCallback function, called on every integration step. brings the
#System's internal x and y up to date with the last integrator's solution step,
#then calls the System's discrete dynamics. as it is, if the System's y
#depends on x or d, and these are modified by f_disc!, the change will not be
#reflected on y until the following integration step
function f_dcb!(integrator)

    @unpack t, u, p = integrator
    @unpack sys, args_c, args_d = p

    sys.x .= u #assign the updated integrator's state to the System's continuous state
    sys.t[] = t #ditto for time

    #at this point sys.y and sys.ẋ hold the values from the last solver
    #evaluation of f_cont!, not the one corresponding to the updated x. with x
    #up to date, we can now compute the final sys.y and sys.ẋ for this epoch
    f_cont!(sys, args_c...)

    #with the System's ẋ and y up to date, call the discrete dynamics function
    x_modified = f_disc!(sys, args_d...)

    u .= sys.x #assign the (potentially modified) sys.x back to the integrator

    u_modified!(integrator, x_modified)

end

#SavingCallback function, gets called at the end of each step after f_disc!
f_scb(x, t, integrator) = deepcopy(integrator.p.sys.y)

step!(mdl::Model, args...) = step!(mdl.integrator, args...)

#makes no sense for a simulation
# solve!(mdl::Model) = solve!(mdl.integrator)

get_proposed_dt(mdl::Model) = get_proposed_dt(mdl.integrator)

function reinit!(mdl::Model, args...; kwargs...)

    #for an ODEIntegrator, the optional args... is simply a new initial
    #condition. if not specified, the original initial condition is used.
    reinit!(mdl.integrator, args...; kwargs...)

    #apparently, reinit! calls f_update!, so all System fields are updated in
    #the process, no need to do it manually

    resize!(mdl.log.t, 1)
    resize!(mdl.log.saveval, 1)
    return nothing
end


################################################################################
################################## TimeHistory #######################################

mutable struct TimeHistory{T}
    _t::Vector{Float64}
    _data::Vector{T}
    function TimeHistory(t::AbstractVector{Float64}, y::AbstractVector{T}) where {T}
        @assert length(t) == length(y)
        new{T}(t, y)
    end
end

TimeHistory(t::Real, y) = TimeHistory([Float64(t)], [y])

function TimeHistory(t::AbstractVector{<:Real}, M::Matrix{<:Real})
    #each Matrix column interpreted as one Vector value
    TimeHistory(t, [M[:, i] for i in 1:size(M,2)])
end

TimeHistory(mdl::Model) = TimeHistory(mdl.log.t, mdl.log.saveval)

Base.length(th::TimeHistory) = length(th._t)

function Base.getproperty(th::TimeHistory, s::Symbol)
    t = getfield(th, :_t)
    y = getfield(th, :_data)
    if s === :_t
        return t
    elseif s === :_data
        return y
    else
        return TimeHistory(t, getproperty(StructArray(y), s))
    end
end

Base.getindex(th::TimeHistory, i) = TimeHistory(th._t[i], th._data[i])

#for inspection
get_child_names(::T) where {T <: TimeHistory} = get_child_names(T)
get_child_names(::Type{TimeHistory{T}}) where {T} = fieldnames(T)

#could be rewritten as @generated to avoid allocation if needed
function get_scalar_components(th::TimeHistory{<:AbstractVector{T}}) where {T<:Real}
    [TimeHistory(th._t, y) for y in th._data |> StructArray |> StructArrays.components]
end



end #module