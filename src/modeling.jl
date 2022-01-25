module Modeling

using Dates
using UnPack
using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks
using ComponentArrays, RecursiveArrayTools

export f_cont!, f_disc!
export SystemDescriptor, SystemGroupDescriptor, NullSystemDescriptor, System, Model


############################# SystemDescriptor ############################

abstract type SystemDescriptor end #anything from which we can build a System

#every SystemDescriptor's init_x must return an AbstractVector{<:Real}, even if
#it's inherently discrete and its f_cont! does nothing. this is ugly but ensures
#composability of Systems without the hassle of dealing automatically with empty
#state vector blocks, which is worsened by the need to assign views from the
#root System's state vector to its children System's state vectors
init_x(::SystemDescriptor) = [0.0]
init_y(::SystemDescriptor) = nothing #sytems are not required to have outputs
init_u(::SystemDescriptor) = nothing #sytems are not required to have control inputs
init_d(::SystemDescriptor) = nothing #systems are not required to have discrete states

############################# System ############################

#need the T type parameter for dispatch, the rest for type stability making
#System mutable does not hurt performance, because Systems should only
#instantiated upon initialization. therefore, no runtime heap allocations
mutable struct System{T, X <: AbstractVector{<:Float64},
                    Y, U, D, P, S}
    ẋ::X #continuous state vector derivative
    x::X #continuous state vector (to be used as a buffer for f_cont! evaluation)
    y::Y #output state
    u::U #control inputs
    d::D #discrete state
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

function System(c::SystemDescriptor, ẋ = init_x(c), x = init_x(c),
                    y = init_y(c), u = init_u(c), d = init_d(c), t = Ref(0.0))

    params = c #assign the system descriptor itself as a System parameter
    subsystems = nothing
    System{map(typeof, (c, x, y, u, d, params, subsystems))...}(
                                    ẋ, x, y, u, d, t, params, subsystems)
end

#f_disc! is free to modify a Hybrid system's discrete state, control inputs and
#continuous state. if it modifies the latter, it must return true, false
#otherwise. no fallbacks are provided for safety reasons: if the intended
#f_cont! or f_disc! implementations for the System have the wrong interface, the
#dispatch will silently revert to the fallback, which does nothing and may not
#be obvious at all.
f_cont!(::System, args...) = MethodError(f_cont!, args) |> throw
(f_disc!(::System, args...)::Bool) = MethodError(f_disc!, args) |> throw

######################### NullSystem ############################

struct NullSystemDescriptor <: SystemDescriptor end

@inline f_cont!(::System{NullSystemDescriptor}, args...) = nothing
@inline (f_disc!(::System{NullSystemDescriptor}, args...)::Bool) = false


######################## SystemGroup #############################

#abstract type providing convenience methods for composite Systems

abstract type SystemGroupDescriptor <: SystemDescriptor end

Base.keys(g::SystemGroupDescriptor) = propertynames(g)
Base.values(g::SystemGroupDescriptor) = map(λ -> getproperty(g, λ), keys(g))

init_x(g::SystemGroupDescriptor) = NamedTuple{keys(g)}(init_x.(values(g))) |> ComponentVector
init_y(g::SystemGroupDescriptor) = NamedTuple{keys(g)}(init_y.(values(g)))
init_u(g::SystemGroupDescriptor) = NamedTuple{keys(g)}(init_u.(values(g)))
init_d(g::SystemGroupDescriptor) = NamedTuple{keys(g)}(init_d.(values(g)))

function System(g::SystemGroupDescriptor, ẋ = init_x(g), x = init_x(g),
                    y = init_y(g), u = init_u(g), d = init_d(g), t = Ref(0.0))

    ss_list = Vector{System}()
    ss_labels = keys(g)

    for label in ss_labels
        push!(ss_list, System(map((λ)->getproperty(λ, label), (g, ẋ, x, y, u, d))..., t))
    end

    params = nothing
    subsystems = NamedTuple{ss_labels}(ss_list)

    System{map(typeof, (g, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)

end

#default implementation calls f_cont! on all Node subsystems with the same
#arguments provided to the parent Node's System, then builds the NamedTuple
#can be overridden for specific SystemGroupDescriptor subtypes if needed
@inline @generated function f_cont!(sys::System{T}, args...) where {T<:SystemGroupDescriptor}

    # Core.print("Generated function called")
    ex_main = Expr(:block)

    #call f_cont! on each subsystem
    ex_calls = Expr(:block)
    ss_labels = fieldnames(T)
    for label in ss_labels
        push!(ex_calls.args,
            :(f_cont!(sys.subsystems[$(QuoteNode(label))], args...)))
    end

    #retrieve the y from each subsystem and build a tuple from them
    ex_tuple = Expr(:tuple)
    for label in ss_labels
        push!(ex_tuple.args,
            :(sys.subsystems[$(QuoteNode(label))].y))
    end

    #build a NamedTuple from the subsystem's labels and the constructed tuple
    ex_y = Expr(:call, Expr(:curly, NamedTuple, ss_labels), ex_tuple)

    #assign the result to the parent system's y
    ex_assign_y = Expr(:(=), :(sys.y), ex_y)

    #push everything into the main block expression
    push!(ex_main.args, ex_calls)
    push!(ex_main.args, ex_assign_y)
    push!(ex_main.args, :(return nothing))

    return ex_main

end

#default implementation calls f_disc! on all Node subsystems with the same
#arguments provided to the parent Node's System, then ORs their outputs.
#can be overridden for specific SystemGroupDescriptor subtypes if needed
@inline @generated function (f_disc!(sys::System{T}, args...)::Bool
    ) where {T<:SystemGroupDescriptor}

    # Core.print("Generated function called")
    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))
    for label in fieldnames(T)
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
struct Model{S <: System,
                   I <: OrdinaryDiffEq.ODEIntegrator,
                   L <: SavedValues}

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

        # x0 = copy(sys.x) #not needed, the integrator creates its own copy
        x0 = sys.x
        problem = ODEProblem{true}(f_update!, x0, (t_start, t_end), params)
        integrator = init(problem, solver; callback = cb_set, save_on, int_kwargs...)
        new{typeof(sys), typeof(integrator), typeof(log)}(sys, integrator, log)
    end
end

#these functions are better defined outside the constructor; closures seem to
#have some overhead (?)

#function barriers: the System is first extracted from integrator.p,
#then used as an argument in the call to the actual update & callback functions,
#forcing the compiler to specialize for the specific System subtype;
#accesing sys.x and sys.ẋ directly instead causes type instability
f_update!(ẋ, x, p, t) = f_update!(ẋ, x, t, p.sys, p.args_c)
f_scb(x, t, integrator) = f_scb(x, t, integrator.p.sys, integrator.p.args_c)
function f_dcb!(integrator)
    x = integrator.u; t = integrator.t; p = integrator.p
    x_modified = f_dcb!(x, t, p.sys, p.args_c, p.args_d)
    u_modified!(integrator, x_modified)
end

#in-place integrator update function
function f_update!(ẋ::X, x::X, t::Real, sys::System{T,X}, args_c) where {T, X}

    sys.x .= x
    sys.t[] = t
    f_cont!(sys, args_c...) #updates sys.ẋ and sys.y
    ẋ .= sys.ẋ

    return nothing
end

#DiscreteCallback function (called on every integration step). in addition to
#specific functionality implemented by f_disc!, this callback ensures that
#the System's internal x and y are up to date with the last integrator's
#solution step
function f_dcb!(x::X, t::Real, sys::System{T,X}, args_c, args_d) where {T,X}

    sys.x .= x #assign the updated integrator's state to the system's local continuous state
    sys.t[] = t #ditto for time

    #at this point sys.y holds the output from the last solver evaluation of
    #f_cont!, not the one corresponding to the updated x. with x up to date, we
    #can now compute the correct sys.y for this epoch
    f_cont!(sys, args_c...) #updates sys.y, but leaves sys.x unmodified
    #this could be commented for additional performance, at the cost of
    #a slight output inaccuracy at integration epochs

    #with the system's outputs up to date, call the discrete update function
    x_modified = f_disc!(sys, args_d...) #this may modify sys.x
    x .= sys.x #assign the (potentially modified) sys.x back to the integrator

    #note that if x is modified by f_disc!, it will not be reflected in sys.y
    #until the next call to f_dcb!

    return x_modified
end

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

SciMLBase.step!(m::Model, args...) = step!(m.integrator, args...)

SciMLBase.solve!(m::Model) = solve!(m.integrator)

SciMLBase.get_proposed_dt(m::Model) = get_proposed_dt(m.integrator)

function SciMLBase.reinit!(m::Model, args...; kwargs...)

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

#the following causes type instability and destroys performance:
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