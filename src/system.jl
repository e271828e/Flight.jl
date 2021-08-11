module System

using DifferentialEquations
using UnPack

export x_template, u_template, init_output, f_output!

abstract type Descriptor end

x_template(::Type{<:Descriptor}) = error("To be overridden by each Descriptor subtype")
u_template(::Type{<:Descriptor}) = error("To be overridden by each Descriptor subtype")
x_template(::D) where {D<:Descriptor} = x_template(D)
u_template(::D) where {D<:Descriptor} = u_template(D)
# init_output(x, u, t, ::Descriptor) = error("To be extended by each Descriptor subtype")
# f_output!(out, x, u, t, ::Descriptor) = error("To be extended by each Descriptor subtype")
init_output(::Any, ::Any, ::Real, ::Descriptor) = error("To be extended by each Descriptor subtype")
f_output!(::Any, ::Any, ::Any, ::Real, ::Descriptor) = error("To be extended by each Descriptor subtype")

struct Continuous{D}

    integrator::OrdinaryDiffEq.ODEIntegrator #just for annotation purposes
    log::SavedValues
    function Continuous(d::D; x₀ = x_template(d), u₀ = u_template(d),
                        t_start = 0.0, t_end = 10.0, method = Tsit5(),
                        output_saveat = Float64[],
                        kwargs...) where {D<:Descriptor}

        function f_step!(ẋ, x, p, t)
            @unpack out, u, d = p
            f_output!(out, x, u, t, d)
            ẋ .= out.ẋ
        end

        function f_save(x, t, integrator)
            @unpack out, u, d = integrator.p
            f_output!(out, x, u, t, d)
            return deepcopy(out.y)
        end

        out₀ = init_output(x₀, u₀, t_start, d)

        params = (u = u₀, out = out₀, d = d)
        log = SavedValues(Float64, typeof(out₀.y))
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