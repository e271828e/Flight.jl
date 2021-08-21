module Model

using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks, RecursiveArrayTools
using UnPack

using Flight.System
# import Flight.System: plotlog

export ContinuousModel

#consider a HybridModel{C}, which generalizes ContinuousModel{C}, providing:
#1) an array xd of discrete states and an array ud of discrete inputs. these
#   should be also stored in the integrator parameters. ud, like u, is modified
#   externally
#2) a Discrete callback called on every integration step, which implements the
#   difference equation xd1 = f(xd0, x0, ud0, u0, t0). this callback should set
#   u_modified! to false
#3) an optional Iterative callback called periodically or after a certain number
#   of integration steps to handle numerical errors (this is the way to
#   implement quaternion renormalization)

############### ContinuousModel #####################

abstract type AbstractModel end

#a ContinuousModel holds a System without a discrete state vector (xd) or
#discrete input vector (ud), but it still provides a discrete function to be
#called on each integration step for bookkeeping (for example, quaternion
#renormalization, )
struct ContinuousModel{I <: OrdinaryDiffEq.ODEIntegrator, L <: SavedValues} <: AbstractModel
# <:OrdinaryDiffEq.ODEIntegrator
# <:SavedValues
    sys::ContinuousSystem
    integrator::I#just for annotation purposes
    log::L

    function ContinuousModel(sys, args_c::Tuple = (), args_d::Tuple = ();
        method = Tsit5(), t_start = 0.0, t_end = 10.0, y_saveat = Float64[],
        int_kwargs...)

        #pass the y cache for f_update! to have somewhere to write to, then
        #throw it away. what matters in this call is the update to the ẋ passed
        #by the integrator
        function f_update!(ẋ, x, p, t)
            @unpack sys, args_c = p
            sys.x .= x
            sys.t[] = t
            f_cont!(sys, args_c...) #updates sys.ẋ and sys.y
            ẋ .= sys.ẋ
            return nothing
        end

        #the dummy ẋ cache is passed for f_update! to have somewhere to write
        #to without clobbering the integrator's du, then it is thrown away. copy
        #and output the updated y
        function f_save(x, t, integrator)
            @unpack sys, args_c = integrator.p
            sys.x .= x
            sys.t[] = t
            f_cont!(sys, args_c...) #updates sys.ẋ and sys.y
            return copy(sys.y)
        end

        function f_dcb!(integrator)
            @unpack sys, args_d = integrator.p
            modified_x = f_disc!(sys, args_d...)
            u_modified!(integrator, modified_x)
        end

        params = (sys = sys, args_c = args_c, args_d = args_d)
        log = SavedValues(Float64, typeof(sys.y))

        dcb = DiscreteCallback((u, t, integrator)->true, f_dcb!)
        scb = SavingCallback(f_save, log, saveat = y_saveat)
        cb_set = CallbackSet(dcb, scb)

        problem = ODEProblem{true}(f_update!, copy(sys.x), (t_start, t_end), params)
        integrator = init(problem, method; callback = cb_set, save_everystep = false, int_kwargs...)
        new{typeof(integrator), typeof(log)}(sys, integrator, log)
    end
end


function Base.getproperty(m::ContinuousModel, s::Symbol)
    #sys also has stored t and x, but since they are written on every call to
    #f_update! we have no guarantees about their status after a certain step
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

end