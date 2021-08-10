module System

using DifferentialEquations

export init_x, init_u

abstract type Descriptor end

init_x(::Descriptor) = nothing
init_u(::Descriptor) = nothing

struct Continuous{D}

    integrator::OrdinaryDiffEq.ODEIntegrator #just for annotation purposes
    log::SavedValues
    function Continuous(d::D, x₀, u₀, f_update!, f_output;
                            t_start = 0.0,
                            t_end = 10.0,
                            method = Tsit5(),
                            output_saveat = Float64[],
                            kwargs...) where {D<:Descriptor}

        f_step!(ẋ, x, p, t) = f_update!(ẋ, x, p.u, t, p.d)
        f_save(x, t, integrator) = f_output(x, integrator.p.u, t, integrator.p.d)

        params = (u = u₀, d = d)
        y₀ = f_output(x₀, u₀, t_start, d)
        log = SavedValues(Float64, typeof(y₀))
        scb = SavingCallback(f_save, log, saveat = output_saveat) #ADD A FLAG TO DISABLE SAVING OPTIONALLY, IT REDUCES ALLOCATIONS

        problem = ODEProblem{true}(f_step!, x₀, (t_start, t_end), params)
        integrator = init(problem, method; callback = scb, save_everystep = false, kwargs...)
        new{D}(integrator, log)
    end
end

Base.getproperty(sys::Continuous, s::Symbol) = getproperty(sys, Val(s))

Base.getproperty(sys::Continuous, ::Val{:descriptor}) = sys.integrator.p.d
Base.getproperty(sys::Continuous, ::Val{:integrator}) = getfield(sys, :integrator)
Base.getproperty(sys::Continuous, ::Val{:log}) = getfield(sys, :log)

#forward everything else to the integrator...
Base.getproperty(sys::Continuous, ::Val{S}) where {S} = getproperty(getfield(sys, :integrator), S)

#...except for x and u (because DiffEqs calls the state u, instead of x)
Base.getproperty(sys::Continuous, ::Val{:u}) = sys.integrator.p.u #input vector
Base.getproperty(sys::Continuous, ::Val{:x}) = sys.integrator.u #state vector

DifferentialEquations.step!(sys::Continuous, args...) = step!(sys.integrator, args...)
DifferentialEquations.reinit!(sys::Continuous, args...) = reinit!(sys.integrator, args...)

end