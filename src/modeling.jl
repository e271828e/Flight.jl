module ModelingTools

using Dates
using UnPack
using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks
using ComponentArrays, RecursiveArrayTools
using Flight.Utils

export get_x0, get_y0, get_u0, get_d0, f_cont!, f_disc!
export AbstractComponent, System, Model


############################# AbstractComponent ############################

abstract type AbstractComponent end #anything from which we can build a System

#every AbstractComponent's get_x0 must return an AbstractVector{<:Real}, even if
#its inherently discrete and its f_cont! does nothing. this is ugly but ensures
#composability of HybridSystems without the hassle of dealing automatically with
#empty state vector blocks, which is magnified by the need to assign views from
#the root System's state vector to each children in its hierarchy
get_x0(::AbstractComponent) = [0.0]
get_y0(::AbstractComponent) = nothing #sytems are not required to have outputs
get_u0(::AbstractComponent) = nothing #sytems are not required to have control inputs
get_d0(::AbstractComponent) = nothing #systems are not required to have discrete states



############################# System ############################

#need the C type parameter for dispatch, the rest for type stability
#making System mutable does not hurt performance, because System instances are
#only instantiated upon initialization, so no runtime heap allocations
mutable struct System{C, X <: AbstractVector{<:Float64},
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

function System(c::AbstractComponent, ẋ = get_x0(c), x = get_x0(c),
                    y = get_y0(c), u = get_u0(c), d = get_d0(c), t = Ref(0.0))
    params = c #assign the component descriptor itself as a System parameter
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
f_cont!(::S, args...) where {S<:System} = no_extend_error(f_cont!, S)
(f_disc!(::S, args...)::Bool) where {S<:System} = no_extend_error(f_disc!, S)



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
        method = Tsit5(), t_start = 0.0, t_end = 10.0, y_saveat = Float64[],
        save_on = false, int_kwargs...)

        #save_on is set to false because we are not usually interested in saving
        #the naked state vector. the output saved by the SavingCallback is all
        #we need for insight
        saveat_arr = (y_saveat isa Real ? (t_start:y_saveat:t_end) : y_saveat)

        params = (sys = sys, args_c = args_c, args_d = args_d)

        # y₀ = f_cont!(sys, args_c...)
        log = SavedValues(Float64, typeof(sys.y))

        dcb = DiscreteCallback((u, t, integrator)->true, f_dcb!)
        scb = SavingCallback(f_scb, log, saveat = saveat_arr)
        cb_set = CallbackSet(dcb, scb)

        x0 = copy(sys.x)
        problem = ODEProblem{true}(f_update!, x0, (t_start, t_end), params)
        integrator = init(problem, method; callback = cb_set, save_on, int_kwargs...)
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
    x_modified = f_dcb!(x, t, p.sys, p.args_d)
    u_modified!(integrator, x_modified)
end

#in-place integrator update function
function f_update!(ẋ::X, x::X, t::Real, sys::System{C,X}, args_c) where {C, X}
    sys.x .= x
    sys.t[] = t
    f_cont!(sys, args_c...) #updates sys.ẋ and sys.y
    ẋ .= sys.ẋ
    return nothing
end

#DiscreteCallback function (called on every integration step). among other
#potential uses provided by f_disc!, this callback ensures that the System's
#internal x is up to date with the integrator's last solution value
function f_dcb!(x::X, t::Real, sys::System{C,X}, args_d) where {C,X}
    sys.x .= x #assign the integrator's state to the system's local continuous state
    sys.t[] = t #ditto for time
    x_modified = f_disc!(sys, args_d...)
    x .= sys.x #assign the (potentially) modified continuous state back to the integrator
    return x_modified
end

#SavingCallback function
function f_scb(x::X, t::Real, sys::System{C,X}, args_c) where {C,X}
    sys.x .= x
    sys.t[] = t
    f_cont!(sys, args_c...)
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