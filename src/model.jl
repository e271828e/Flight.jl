module Model

using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks, RecursiveArrayTools
using UnPack

using Flight.System
import Flight.System: plotlog

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
struct ContinuousModel{S, I <: OrdinaryDiffEq.ODEIntegrator, L <: SavedValues} <: AbstractModel
# <:OrdinaryDiffEq.ODEIntegrator
# <:SavedValues
    integrator::I#just for annotation purposes
    log::L

    function ContinuousModel(sys::S; x₀ = X(sys), u₀ = U(sys), data = D(sys),
        method = Tsit5(), t_start = 0.0, t_end = 10.0, y_saveat = Float64[],
        kwargs...) where {S<:AbstractSystem}

        #pass the y cache for f_update! to have somewhere to write to, then
        #throw it away. what matters in this call is the update to the ẋ passed
        #by the integrator
        function f_update!(ẋ, x, p, t)
            @unpack y_tmp, u, data, sys = p
            f_cont!(y_tmp, ẋ, x, u, t, data, sys) #throw away y
        end

        function f_dcb!(integrator)
            @unpack u, data, sys = integrator.p
            t = integrator.t
            x = integrator.u
            modified_x = f_disc!(x, u, t, data, sys)
            u_modified!(integrator, modified_x)
        end

        #the dummy ẋ cache is passed for f_update! to have somewhere to write
        #to without clobbering the integrator's du, then it is thrown away. copy
        #and output the updated y
        function f_save(x, t, integrator)
            @unpack y_tmp, ẋ_tmp, u, data, sys = integrator.p
            f_cont!(y_tmp, ẋ_tmp, x, u, t, data, sys)
            return copy(y_tmp)
        end

        params = (u = u₀, y_tmp = Y(sys), ẋ_tmp = X(sys), data = data, sys = sys)
        log = SavedValues(Float64, typeof(params.y_tmp))

        dcb = DiscreteCallback((u, t, integrator)->true, f_dcb!)
        scb = SavingCallback(f_save, log, saveat = y_saveat)
        cb_set = CallbackSet(dcb, scb)

        problem = ODEProblem{true}(f_update!, x₀, (t_start, t_end), params)
        integrator = init(problem, method; callback = cb_set, save_everystep = false, kwargs...)
        new{S, typeof(integrator), typeof(log)}(integrator, log)
    end
end


function Base.getproperty(m::ContinuousModel, s::Symbol)
    if s === :t
        return m.integrator.t
    elseif s === :x
        return m.integrator.u
    elseif s in (:integrator, :log)
        return getfield(m, s)
    elseif s in (:u, :ẋ_tmp, :y_tmp, :data, :sys)
        return getproperty(m.integrator.p, s)
    else
        return getproperty(m.integrator, s)
    end
end

# Base.getproperty(m::ContinuousModel, s::Symbol) = getproperty(m, Val(s))

# Base.getproperty(m::ContinuousModel, ::Val{:integrator}) = getfield(m, :integrator)
# Base.getproperty(m::ContinuousModel, ::Val{:log}) = getfield(m, :log)

# #forward everything else to the integrator...
# Base.getproperty(m::ContinuousModel, ::Val{S}) where {S} = getproperty(getfield(m, :integrator), S)

# #...except for x (because DiffEqs calls the state u, instead of x)
# Base.getproperty(m::ContinuousModel, ::Val{:x}) = m.integrator.u #state vector

# Base.getproperty(m::ContinuousModel, ::Val{:u}) = m.integrator.p.u #input vector
# Base.getproperty(m::ContinuousModel, ::Val{:ẋ}) = m.integrator.p.ẋ #ẋ cache
# Base.getproperty(m::ContinuousModel, ::Val{:y}) = m.integrator.p.y #y cache
# Base.getproperty(m::ContinuousModel, ::Val{:data}) = m.integrator.p.data #external data cache
# Base.getproperty(m::ContinuousModel, ::Val{:sys}) = m.integrator.p.sys


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