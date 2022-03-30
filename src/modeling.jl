module Modeling

using Dates
using UnPack
using ComponentArrays, RecursiveArrayTools
using SciMLBase: ODEProblem, u_modified!, init as init_problem
using OrdinaryDiffEq: ODEIntegrator, Tsit5
using DiffEqCallbacks: SavingCallback, DiscreteCallback, CallbackSet, SavedValues

import SciMLBase: step!, solve!, reinit!, get_proposed_dt
import DataStructures: OrderedDict
import Flight.Plotting: plots

export f_cont!, f_disc!
export SystemDescriptor, SystemGroupDescriptor, NullSystemDescriptor, System, Model
export SystemẊ, SystemX, SystemY, SystemU, SystemD


############################# SystemDescriptor ############################

abstract type SystemDescriptor end #anything from which we can build a System

abstract type SystemData end

struct SystemẊ <: SystemData end
struct SystemX <: SystemData end
struct SystemY <: SystemData end
struct SystemU <: SystemData end
struct SystemD <: SystemData end

############################# System ############################

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

############################# Model ############################

#in this design, the t and x fields of m.sys behave only as temporary
#storage for f_cont! and f_disc! calls, so we have no guarantees about their
#status after a certain step. the only valid sources for t and x at any
#given moment is the integrator's t and u
struct Model{S <: System, I <: ODEIntegrator, L <: SavedValues}
    sys::S
    integrator::I
    log::L

    function Model(sys, args_c::Tuple = (), args_d::Tuple = ();
        solver = Tsit5(), t_start = 0.0, t_end = 10.0, y_saveat = Float64[],
        save_on = false, int_kwargs...)

        #save_on is set to false because we are not usually interested in saving
        #the naked state vector. the output saved by the SavingCallback is all
        #we need for insight
        saveat_arr = (y_saveat isa Real ? (t_start:y_saveat:t_end) : y_saveat)

        params = (sys = sys, args_c = args_c, args_d = args_d)

        log = SavedValues(Float64, typeof(sys.y))

        dcb = DiscreteCallback((u, t, integrator)->true, f_dcb!)
        scb = SavingCallback(f_scb, log, saveat = saveat_arr)
        cb_set = CallbackSet(dcb, scb)

        x0 = sys.x #the integrator creates its own copy
        problem = ODEProblem{true}(f_update!, x0, (t_start, t_end), params)
        integrator = init_problem(problem, solver; callback = cb_set, save_on, int_kwargs...)
        new{typeof(sys), typeof(integrator), typeof(log)}(sys, integrator, log)
    end
end

#function barriers: the System and arguments to f_cont! are first extracted from
#integrator parameters, then used as an argument in the call to the actual
#update & callback functions, forcing the compiler to specialize; accesing sys.x
#and sys.ẋ directly instead causes type instability

#function barrier for integrator update
f_update!(ẋ, x, p, t) = f_update!(ẋ, x, t, p.sys, p.args_c)

#in-place integrator update function
function f_update!(ẋ::X, x::X, t::Real, sys::System{T,X}, args_c) where {T, X}
    sys.x .= x
    sys.t[] = t
    f_cont!(sys, args_c...) #updates sys.ẋ and sys.y
    ẋ .= sys.ẋ
    return nothing
end

#function barrier for discrete callback
function f_dcb!(integrator)
    x = integrator.u; t = integrator.t; p = integrator.p
    x_modified = f_dcb!(x, t, p.sys, p.args_c, p.args_d)
    u_modified!(integrator, x_modified)
end

#DiscreteCallback function (called on every integration step). this callback
#brings the System's internal x and y up to date with the last integrator's
#solution step, then executes the System's discrete update function
function f_dcb!(x::X, t::Real, sys::System{T,X}, args_c, args_d) where {T,X}

    sys.x .= x #assign the updated integrator's state to the system's local continuous state
    sys.t[] = t #ditto for time

    #at this point sys.y and sys.ẋ hold the values from the last solver evaluation of
    #f_cont!, not the one corresponding to the updated x. with x up to date, we
    #can now compute the correct sys.y sys.ẋ for this epoch
    f_cont!(sys, args_c...) #updates sys.y, but leaves sys.x unmodified

    #with the system's outputs up to date, call the discrete update function
    x_modified = f_disc!(sys, args_d...) #this may modify sys.x
    x .= sys.x #assign the (potentially modified) sys.x back to the integrator

    #note: as it is, if the System's y depends on x or d, and these are modified
    #by f_disc!, the change will not be reflected on y until the following
    #integration step

    return x_modified
end

#function barrier for saving callback
f_scb(x, t, integrator) = f_scb(x, t, integrator.p.sys, integrator.p.args_c)

#SavingCallback function, this gets called at the end of each step after f_disc!
function f_scb(::X, ::Real, sys::System{T,X}, args_c) where {T,X}
    return deepcopy(sys.y)
end

function Base.getproperty(m::Model, s::Symbol)
    if s === :t
        return m.integrator.t
    elseif s === :x
        return m.integrator.u
    elseif s === :y
        return m.sys.y
    elseif s === :u
        return m.sys.u
    elseif s ∈ (:sys, :integrator, :log)
        return getfield(m, s)
    else
        return getproperty(m.integrator, s)
    end
end

step!(m::Model, args...) = step!(m.integrator, args...)

solve!(m::Model) = solve!(m.integrator)

get_proposed_dt(m::Model) = get_proposed_dt(m.integrator)

function reinit!(m::Model, args...; kwargs...)

    #for an ODEIntegrator, the optional args... is simply a new initial
    #condition. if not specified, the original initial condition is used
    reinit!(m.integrator, args...; kwargs...)

    #grab the updated t and x from the integrator (in case they were reset by
    #the input arguments). this is not strictly necessary, since they are merely
    #buffers. just for consistency.
    m.sys.t[] = m.integrator.t
    m.sys.x .= m.integrator.u

    resize!(m.log.t, 1)
    resize!(m.log.saveval, 1)
    return nothing
end

function plots(mdl::Model; mode::Symbol = :basic,
    save_path::Union{String,Nothing} = nothing, kwargs...)
    #generate default path tmp/plots/current_date
    save_path = (save_path === nothing ?
        joinpath("tmp", Dates.format(now(), "yyyy_mm_dd_HHMMSS")) : save_path)
    mkpath(save_path)
    plots(mdl.log.t, mdl.log.saveval; mode, save_path, kwargs...)
end

#the following causes type instability and kills performance:
# function f_update!(ẋ, x, p, t)
    # @unpack sys, args_c = p
    # sys.x .= x
    # sys.t[] = t
    # f_cont!(sys, args_c...)
    # ẋ = sys.ẋ
# end

# the reason seems to be that having sys stored in p obfuscates type inference.
# when unpacking sys, the compiler can no longer tell its type, and therefore
# has no knowledge of the types of sys.x, sys.dx, sys.y and sys.t. since these
# are being assigned to and read from, the type instability kills performance.

# this can be fixed by storing the x, dx and y fields of sys directly as entries
# of p. this probably fixes their types during construction, so when they are
# accessed later in the closure, the type instability is no longer an issue.

# however, this is redundant! we already have x, dx, y and t inside of sys. a
# more elegant alternative is simply to use a function barrier, first extract
# sys, then call another function using it as an argument. this forces the
# compiler to infer its type, and therefore it specializes the time-critical
# assignment statements to their actual types.

###########################################################################

#to try:

#in System, define and extend f_branch!

# #individual Component
# f_branch!(y, dx, x, u, t, sys, args...) = f_branch!(Val(has_input(sys)), y, dx, x, u, t, args...)
# f_branch!(::Val{true}, y, dx, x, u, t, sys, args...) = f_cont!(y, dx, x, u, t, sys, args...)
# f_cont!(::HasInput, y, dx, x ,u, t, sys, args...) = f_cont!(y, dx, x, u, t, sys, args...)
# f_cont!(::HasNoInput, y, dx, x, u, t, sys, args...) = f_cont!(y, dx, x, t, sys, args...)

# #for a AirframeGroup
# f_cont!(MaybeInput(S), MaybeOutput(S), y, dx, x, u, t, sys, args...)
# f_cont!(::HasInput, ::HasOutput, y, dx, x ,u, t, sys, args...)
# #now, this method needs to consider the possibility for each component that it
# #may have or not Input or Output. so it must do
# for (label, component) in zip(keys(C), values(C))
#     if MaybeInput(typeof(component)) #need tocheck, because if it has no input, u[label] will not exist!
#         f_cont!(y_cmp, dx_cmp, x_cmp, u_cmp, t, cmp, args...)
#     else
#         f_cont!(y_cmp, dx_cmp, x_cmp, t, cmp, args...)
#     end
# end

end