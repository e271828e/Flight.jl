module Model

using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks
using UnPack

using Flight.System

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

struct ContinuousModel{S} <: AbstractModel

    integrator::OrdinaryDiffEq.ODEIntegrator #just for annotation purposes
    log::SavedValues

    function ContinuousModel(sys::S; x₀ = X(sys), u₀ = U(sys), data = D(sys),
        method = Tsit5(), t_start = 0.0, t_end = 10.0, y_saveat = Float64[],
        kwargs...) where {S<:AbstractSystem}

        #pass the y cache for f_update! to have somewhere to write to, then
        #throw it away. what matters in this call is the update to the ẋ passed
        #by the integrator
        function f_step!(ẋ, x, p, t)
            @unpack u, data, sys = p
            f_output!(ẋ, x, u, t, data, sys) #throw away y
        end

        #the dummy ẋ cache is passed for f_update! to have somewhere to write
        #to without clobbering the integrator's du, then it is thrown away. copy
        #and output the updated y
        function f_save(x, t, integrator)
            @unpack ẋ_dummy, u, data, sys = integrator.p
            y = f_output!(ẋ_dummy, x, u, t, data, sys)
            return y
        end

        y₀, ẋ₀ = output_init(x₀, u₀, t_start, data, sys)

        params = (u = u₀, ẋ_dummy = ẋ₀, data = data, sys = sys)
        log = SavedValues(Float64, typeof(y₀))
        scb = SavingCallback(f_save, log, saveat = y_saveat)

        problem = ODEProblem{true}(f_step!, x₀, (t_start, t_end), params)
        integrator = init(problem, method; callback = scb, save_everystep = false, kwargs...)
        new{S}(integrator, log)
    end
end

Base.getproperty(m::ContinuousModel, s::Symbol) = getproperty(m, Val(s))

Base.getproperty(m::ContinuousModel, ::Val{:system}) = m.integrator.p.sys
Base.getproperty(m::ContinuousModel, ::Val{:integrator}) = getfield(m, :integrator)
Base.getproperty(m::ContinuousModel, ::Val{:log}) = getfield(m, :log)

#forward everything else to the integrator...
Base.getproperty(m::ContinuousModel, ::Val{S}) where {S} = getproperty(getfield(m, :integrator), S)

#...except for x and u (because DiffEqs calls the state u, instead of x)
Base.getproperty(m::ContinuousModel, ::Val{:x}) = m.integrator.u #state vector
Base.getproperty(m::ContinuousModel, ::Val{:u}) = m.integrator.p.u #input vector

Base.getproperty(m::ContinuousModel, ::Val{:ẋ}) = m.integrator.p.ẋ #ẋ cache
Base.getproperty(m::ContinuousModel, ::Val{:y}) = m.integrator.p.y #y cache
Base.getproperty(m::ContinuousModel, ::Val{:data}) = m.integrator.p.data #external data cache

SciMLBase.step!(m::ContinuousModel, args...) = step!(m.integrator, args...)
SciMLBase.reinit!(m::ContinuousModel, args...) = reinit!(m.integrator, args...)


end