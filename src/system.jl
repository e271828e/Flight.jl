module System

using DifferentialEquations

export SystemDescriptor, ContinuousSystem

abstract type SystemDescriptor end

struct ContinuousSystem{D}

    integrator::OrdinaryDiffEq.ODEIntegrator #just for annotation purposes
    log::SavedValues
    function ContinuousSystem(d::D, x₀, u₀, f_update!, f_output; method = Tsit5(), kwargs...) where {D<:SystemDescriptor}

        params = (u = u₀, d = d)
        y₀ = f_output(x₀, u₀, 0, d)
        log = SavedValues(Float64, typeof(y₀))
        scb = SavingCallback(f_output, log)

        problem = ODEProblem{true}(f_update!, x₀, (0, Inf), params)
        integrator = init(problem, method; callback = scb, save_everystep = false, kwargs...)
        new{D}(integrator, log)
    end
end

Base.getproperty(sys::ContinuousSystem, s::Symbol) = getproperty(sys, Val(s))

Base.getproperty(sys::ContinuousSystem, ::Val{:descriptor}) = sys.integrator.p.descriptor
Base.getproperty(sys::ContinuousSystem, ::Val{:integrator}) = getfield(sys, :integrator)
Base.getproperty(sys::ContinuousSystem, ::Val{:log}) = getfield(sys, :log)

#forward everything else to the integrator...
Base.getproperty(sys::ContinuousSystem, ::Val{S}) where {S} = getproperty(getfield(sys, :integrator), S)

#...except for x and u (because DiffEqs calls the state u, instead of x)
Base.getproperty(sys::ContinuousSystem, ::Val{:u}) = sys.integrator.p.u #input vector
Base.getproperty(sys::ContinuousSystem, ::Val{:x}) = sys.integrator.u #state vector

# f_update!(ẋ, x, p, t) = f_update!(ẋ, x, p.u, t, p.d)
# f_output(x, t, integrator) = f_output(x, integrator.p.u, t, integrator.p.d)

DifferentialEquations.step!(sys::ContinuousSystem) = step!(sys.integrator)
DifferentialEquations.step!(sys::ContinuousSystem, dt::Real) = step!(sys.integrator, dt, stop_at_tdt = true)

end