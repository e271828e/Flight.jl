module System

using DifferentialEquations
using UnPack

export x_template, u_template, init_data, init_output, f_output!

abstract type Component end

x_template(::Type{<:Component}) = error("To be overridden by each Component subtype")
u_template(::Type{<:Component}) = error("To be overridden by each Component subtype")
x_template(::C) where {C<:Component} = x_template(C)
u_template(::C) where {C<:Component} = u_template(C)
# init_output(x, u, t, ::Component) = error("To be extended by each Component subtype")
# f_output!(out, x, u, t, ::Component) = error("To be extended by each Component subtype")
init_data(::Any, ::Any, ::Real, ::Component) = error("To be extended by each Component subtype")
init_output(::Any, ::Any, ::Real, ::Any, ::Component) = error("To be extended by each Component subtype")
f_output!(::Any, ::Any, ::Any, ::Any, ::Real, ::Any, ::Component) = error("To be extended by each Component subtype")

#also consider Hybrid{C} systems, which in addition to what a Continuous system
#has, should provide:
#1) an array xd of discrete states and an array ud of discrete inputs. these
#   should be also stored in the integrator parameters. ud, like u, is modified
#   externally
#2) a Discrete callback called on every integration step, which implements the
#   difference equation xd1 = f(xd0, x0, ud0, u0, t0). this callback should set
#   u_modified! to false
#3) an optional Iterative callback called periodically to handle numerical
#   errors (this is the way to implement quaternion renormalization)

struct Continuous{C}

    integrator::OrdinaryDiffEq.ODEIntegrator #just for annotation purposes
    log::SavedValues
    function Continuous(comp::C; x₀ = x_template(comp), u₀ = u_template(comp),
                        t_start = 0.0, t_end = 10.0, method = Tsit5(),
                        y_saveat = Float64[],
                        kwargs...) where {C<:Component}

        #pass the y cache for f_update! to have somewhere to write to, then
        #throw it away. what matters in this call is the update to the ẋ passed
        #by the integrator
        function f_step!(ẋ, x, p, t)
            @unpack y, u, data, comp = p
            f_output!(y, ẋ, x, u, t, data, comp)
        end

        #the dummy ẋ cache is passed for f_update! to have somewhere to write
        #to without cloberring the integrator's du, then it is thrown away. copy
        #and output the updated y
        function f_save(x, t, integrator)
            @unpack y, ẋ, u, data, comp = integrator.p
            f_output!(y, ẋ, x, u, t, data, comp)
            return deepcopy(y)
        end

        #initialize external data sources and output cache for this component
        d₀ = init_data(x₀, u₀, t_start, comp)
        (y₀, ẋ₀) = init_output(x₀, u₀, t_start, d₀, comp)

        #store
        params = (u = u₀, y = y₀, ẋ = ẋ₀, data = d₀, comp = comp)
        log = SavedValues(Float64, typeof(y₀))
        scb = SavingCallback(f_save, log, saveat = y_saveat) #ADD A FLAG TO DISABLE SAVING OPTIONALLY, IT REDUCES ALLOCATIONS

        problem = ODEProblem{true}(f_step!, x₀, (t_start, t_end), params)
        integrator = init(problem, method; callback = scb, save_everystep = false, kwargs...)
        new{C}(integrator, log)
    end
end

Base.getproperty(sys::Continuous, s::Symbol) = getproperty(sys, Val(s))

Base.getproperty(sys::Continuous, ::Val{:component}) = sys.integrator.p.comp
Base.getproperty(sys::Continuous, ::Val{:integrator}) = getfield(sys, :integrator)
Base.getproperty(sys::Continuous, ::Val{:log}) = getfield(sys, :log)

#forward everything else to the integrator...
Base.getproperty(sys::Continuous, ::Val{S}) where {S} = getproperty(getfield(sys, :integrator), S)

#...except for x and u (because DiffEqs calls the state u, instead of x)
Base.getproperty(sys::Continuous, ::Val{:x}) = sys.integrator.u #state vector
Base.getproperty(sys::Continuous, ::Val{:u}) = sys.integrator.p.u #input vector

Base.getproperty(sys::Continuous, ::Val{:ẋ}) = sys.integrator.p.ẋ #ẋ cache
Base.getproperty(sys::Continuous, ::Val{:y}) = sys.integrator.p.y #y cache
Base.getproperty(sys::Continuous, ::Val{:data}) = sys.integrator.p.data #external data cache

DifferentialEquations.step!(sys::Continuous, args...) = step!(sys.integrator, args...)
DifferentialEquations.reinit!(sys::Continuous, args...) = reinit!(sys.integrator, args...)

end