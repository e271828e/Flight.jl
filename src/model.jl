module Model

using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks, RecursiveArrayTools
using UnPack

using Flight.System
# import Flight.System: plotlog

export ContinuousModel

#consider a HybridModel{C}, which generalizes ContinuousModel{C}, providing:
#1) an array xd of discrete states
#2) a Discrete callback called on every integration step, which implements the
#   difference equation xd1 = f(xd0, x0, ud0, u0, t0). this callback should set
#   u_modified! to false

abstract type AbstractModel end

############### ContinuousModel #####################

#a ContinuousModel holds a System without a discrete state vector (xd) or
#discrete input vector (ud), but it still provides a discrete function to be
#called on each integration step for bookkeeping (for example, quaternion
#renormalization, )

#in this design, the t and x fields of m.sys behave only as temporary
#storage for f_cont! and f_disc! calls, so we have no guarantees about their
#status after a certain step. the only valid sources for t and x at any
#given moment is the integrator's t and u
struct ContinuousModel{S, I <: OrdinaryDiffEq.ODEIntegrator, L <: SavedValues} <: AbstractModel

    sys::S
    integrator::I#just for annotation purposes
    log::L

    function ContinuousModel(sys, args_c::Tuple = (), args_d::Tuple = ();
        method = Tsit5(), t_start = 0.0, t_end = 10.0, y_saveat = Float64[],
        int_kwargs...)

        params = (sys = sys, args_c = args_c, args_d = args_d)

        log = SavedValues(Float64, typeof(sys.y))

        dcb = DiscreteCallback((u, t, integrator)->true, f_dcb!)
        scb = SavingCallback(f_scb, log, saveat = y_saveat)
        cb_set = CallbackSet(dcb, scb)

        problem = ODEProblem{true}(f_update!, copy(sys.x), (t_start, t_end), params)
        integrator = init(problem, method; callback = cb_set, save_everystep = false, int_kwargs...)
        new{typeof(sys), typeof(integrator), typeof(log)}(sys, integrator, log)
    end
end

#these functions are better defined outside the constructor; apparently closures
#are not as efficient

#function barriers: the ContinuousSystem is first extracted from integrator.p,
#then used as an argument in the call to the actual update & callback functions,
#forcing the compiler to specialize those
f_update!(ẋ, x, p, t) = f_update!(ẋ, x, t, p.sys, p.args_c)
f_scb(x, t, integrator) = f_scb(x, t, integrator.p.sys, integrator.p.args_c)
function f_dcb!(integrator)
    x = integrator.u; t = integrator.t; p = integrator.p
    x_modified = f_dcb!(x, t, p.sys, p.args_d)
    u_modified!(integrator, x_modified)
end

#in-place integrator update function
function f_update!(ẋ::X, x::X, t::Real, sys::ContinuousSystem{C,X}, args_c) where {C, X}
    sys.x .= x
    sys.t[] = t
    f_cont!(sys, args_c...) #updates sys.ẋ and sys.y
    ẋ .= sys.ẋ
end

#DiscreteCallback function (called on every integration step)
function f_dcb!(x::X, t::Real, sys::ContinuousSystem{C,X}, args_d) where {C,X}
    sys.x .= x #assign the integrator's state to the system's local continuous state
    sys.t[] = t #ditto for time
    x_modified = f_disc!(sys, args_d...)
    x .= sys.x #assign the (potentially) modified continuous state back to the integrator
    # println(x_modified)
    return x_modified
end

#SavingCallback function
function f_scb(x::X, t::Real, sys::ContinuousSystem{C,X}, args_c) where {C,X}
    sys.x .= x
    sys.t[] = t
    f_cont!(sys, args_c...) #updates sys.ẋ and sys.y
    return copy(sys.y)
end


function Base.getproperty(m::ContinuousModel, s::Symbol)
    if s === :t
        return m.integrator.t
    elseif s === :x
        return m.integrator.u
    elseif s === :u
        return m.sys.u
    elseif s === :y #should only retrieve this after a saving step
        return m.sys.y
    elseif s in (:sys, :integrator, :log)
        return getfield(m, s)
    else
        return getproperty(m.integrator, s)
    end
end

SciMLBase.step!(m::ContinuousModel, args...) = step!(m.integrator, args...)

function SciMLBase.reinit!(m::ContinuousModel, args...)
    reinit!(m.integrator, args...)

    #restore ContinuousSystems internal time and state to their original values
    #(which were stored by the integrator upon construction, and now supplied by
    #it)
    m.sys.t[] = m.integrator.t
    m.sys.x .= m.integrator.u

    resize!(m.log.t, 1)
    resize!(m.log.saveval, 1)
    return nothing
end

#replace this with the appropriate Plot recipes, etc
#careful with overloading plot without importing. see how it is done properly in
#Plots
function plotlog(mdl::ContinuousModel)

    t = mdl

    #this produces a Vector of Y(aircraft) with length(saveval) which can be
    #indexed as a matrix
    v = VectorOfArray(mdl.log.saveval)

    #this produces an actual ComponentMatrix of size [length(Y(aircraft)),
    #length(saveval)] matrix. the awesome part is that its row dimension
    #preserves the Axes metadata in the original Y(aircraft). its column
    #dimension has a Flat Axis.
    m = convert(Array, v)

    #after this preprocessing the model's System and each sub-component in its
    #hierarchy can recursively extract each block and delegate its plots to
    #its children
    log = (t = mdl.log.t, y = m)
    plotlog(log, mdl.sys)

end

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