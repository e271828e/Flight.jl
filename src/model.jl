module Model

using Dates
using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks, RecursiveArrayTools
using UnPack

using Flight.System
import Flight.Plotting: plots

export HybridModel

abstract type AbstractModel end

############### HybridModel #####################

#in this design, the t and x fields of m.sys behave only as temporary
#storage for f_cont! and f_disc! calls, so we have no guarantees about their
#status after a certain step. the only valid sources for t and x at any
#given moment is the integrator's t and u
struct HybridModel{S <: HybridSystem, I <: OrdinaryDiffEq.ODEIntegrator, L <: SavedValues} <: AbstractModel

    sys::S
    integrator::I
    log::L

    function HybridModel(sys, args_c::Tuple = (), args_d::Tuple = ();
        method = Tsit5(), t_start = 0.0, t_end = 10.0, y_saveat = Float64[],
        int_kwargs...)

        params = (sys = sys, args_c = args_c, args_d = args_d)

        y₀ = f_cont!(sys, args_c...)
        log = SavedValues(Float64, typeof(y₀))

        dcb = DiscreteCallback((u, t, integrator)->true, f_dcb!)
        scb = SavingCallback(f_scb, log, saveat = y_saveat)
        cb_set = CallbackSet(dcb, scb)

        x0 = copy(sys.x)
        problem = ODEProblem{true}(f_update!, x0, (t_start, t_end), params)
        integrator = init(problem, method; callback = cb_set, save_everystep = false, int_kwargs...)
        new{typeof(sys), typeof(integrator), typeof(log)}(sys, integrator, log)
    end
end

#these functions are better defined outside the constructor; apparently closures
#are not as efficient

#function barriers: the HybridSystem is first extracted from integrator.p,
#then used as an argument in the call to the actual update & callback functions,
#forcing the compiler to specialize for the specific HybridSystem subtype;
#accesing sys.x and sys.ẋ directly without this causes type instability
f_update!(ẋ, x, p, t) = f_update!(ẋ, x, t, p.sys, p.args_c)
f_scb(x, t, integrator) = f_scb(x, t, integrator.p.sys, integrator.p.args_c)
function f_dcb!(integrator)
    x = integrator.u; t = integrator.t; p = integrator.p
    x_modified = f_dcb!(x, t, p.sys, p.args_d)
    u_modified!(integrator, x_modified)
end

#in-place integrator update function
function f_update!(ẋ::X, x::X, t::Real, sys::HybridSystem{C,X}, args_c) where {C, X}
    sys.x .= x
    sys.t[] = t
    f_cont!(sys, args_c...) #updates sys.ẋ and sys.y
    ẋ .= sys.ẋ
    return nothing
end

#DiscreteCallback function (called on every integration step)
function f_dcb!(x::X, t::Real, sys::HybridSystem{C,X}, args_d) where {C,X}
    sys.x .= x #assign the integrator's state to the system's local continuous state
    sys.t[] = t #ditto for time
    x_modified = f_disc!(sys, args_d...)
    x .= sys.x #assign the (potentially) modified continuous state back to the integrator
    # println(x_modified)
    return x_modified
end

#SavingCallback function
function f_scb(x::X, t::Real, sys::HybridSystem{C,X}, args_c) where {C,X}
    sys.x .= x
    sys.t[] = t
    return f_cont!(sys, args_c...)
end


function Base.getproperty(m::HybridModel, s::Symbol)
    if s === :t
        return m.integrator.t
    elseif s === :x
        return m.integrator.u
    elseif s === :u
        return m.sys.u
    elseif s in (:sys, :integrator, :log)
        return getfield(m, s)
    else
        return getproperty(m.integrator, s)
    end
end

SciMLBase.step!(m::HybridModel, args...) = step!(m.integrator, args...)

function SciMLBase.reinit!(m::HybridModel, args...)
    reinit!(m.integrator, args...)

    #restore HybridSystems internal time and state to their original values
    #(which were stored by the integrator upon construction, and now supplied by
    #it)
    m.sys.t[] = m.integrator.t
    m.sys.x .= m.integrator.u

    resize!(m.log.t, 1)
    resize!(m.log.saveval, 1)
    return nothing
end

function plots(mdl::HybridModel; mode::Symbol = :basic, save_path::Union{String,Nothing} = nothing, kwargs...)
    #generate default path tmp/plots/current_date
    save_path = (save_path === nothing ? joinpath("tmp", Dates.format(now(), "yyyy_mm_dd_HHMMSS")) : save_path)
    mkpath(save_path)
    plots(mdl.log.t, mdl.log.saveval; mode, save_path, kwargs...)
end

# #delegate recursive plotting to the simulated System's AbstractY rplot method
# function System.rplot(log::DiffEqCallbacks.SavedValues, args...)
#     isempty(log.t) ? error("Can't plot an empty log") : rplot(log.t, log.saveval, args...)
# end

#the following causes type instability and destroys performance:
# function f_update!(ẋ, x, p, t)
    # @unpack sys, args_c = p
    # sys.x .= x
    # sys.t[] = t
    # f_cont!(sys, args_c...)
    # ẋ = sys.ẋ
# end

#the reason seems to be that having sys stored in p obfuscates type
#inference. when unpacking sys, the compiler can no longer tell its
#type, and therefore has no knowledge of the types of sys.x, sys.dx,
#sys.y and sys.t. since these are being assigned to and read from,
#the type instability kills performance.

#this can be fixed by storing the x, dx and y fields of sys directly
#as entries of p. this probably fixes their types during
#construction, so when they are accessed later in the closure, the
#type instability is no longer an issue.

#however, this is redundant! we already have x, dx, y and t inside
#of sys. a more elegant alternative is simply to use a function
#barrier, first extract sys, then call another function using it as
#an argument. this forces the compiler to infer its type, and
#therefore it specializes the time-critical assignment statements to
#their actual types.
end