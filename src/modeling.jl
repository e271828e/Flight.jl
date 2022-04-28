module Modeling

using UnPack
using StructArrays
using SciMLBase: ODEProblem, u_modified!, init as init_problem
using OrdinaryDiffEq: ODEIntegrator, RK4
using DiffEqCallbacks: SavingCallback, DiscreteCallback, CallbackSet, SavedValues

using ..Systems

import SciMLBase: step!, solve!, reinit!, get_proposed_dt

export Model, TimeHistory


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
        solver = RK4(), #
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