module System

using DifferentialEquations
using ComponentArrays
using UnPack

export AbstractComponent, ComponentGroup, ContinuousSystem
export x_template, u_template, y_template, d_template, init_output, f_output!

################## AbstractComponent Interface ###################

abstract type AbstractComponent end

x_template(::C) where {C<:AbstractComponent} = x_template(C)
u_template(::C) where {C<:AbstractComponent} = u_template(C)
y_template(::C) where {C<:AbstractComponent} = y_template(C)
d_template(::C) where {C<:AbstractComponent} = d_template(C)

x_template(::Type{<:AbstractComponent}) = error("To be overridden by each subtype")
u_template(::Type{<:AbstractComponent}) = error("To be overridden by each subtype")
y_template(::Type{<:AbstractComponent}) = error("To be overridden by each subtype")
d_template(::Type{<:AbstractComponent}) = error("To be overridden by each subtype")

f_output!(::Any, ::Any, ::Any, ::Any, ::Real, ::Any, ::AbstractComponent) = error("To be extended by each subtype")

function init_output(x::Any, u::Any, t::Real, data::Any, comp::AbstractComponent)
    ẋ = x_template(comp)
    y = y_template(comp)
    f_output!(y, ẋ, x, u, t, data, comp)
    return y, ẋ
end

#consider a HybridSystem{C}, which generalizes ContinuousSystem{C}, providing:
#1) an array xd of discrete states and an array ud of discrete inputs. these
#   should be also stored in the integrator parameters. ud, like u, is modified
#   externally
#2) a Discrete callback called on every integration step, which implements the
#   difference equation xd1 = f(xd0, x0, ud0, u0, t0). this callback should set
#   u_modified! to false
#3) an optional Iterative callback called periodically or after a certain number
#   of integration steps to handle numerical errors (this is the way to
#   implement quaternion renormalization)

############### Continuous System #####################

struct ContinuousSystem{C}

    integrator::OrdinaryDiffEq.ODEIntegrator #just for annotation purposes
    log::SavedValues

    function ContinuousSystem(comp::C;
        x₀ = x_template(comp), u₀ = u_template(comp),
        t_start = 0.0, t_end = 10.0, method = Tsit5(), y_saveat = Float64[],
        kwargs...) where {C<:AbstractComponent}

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

        data = d_template(comp)
        (y₀, ẋ₀) = init_output(x₀, u₀, t_start, data, comp)

        params = (u = u₀, y = y₀, ẋ = ẋ₀, data = data, comp = comp)
        log = SavedValues(Float64, typeof(y₀))
        scb = SavingCallback(f_save, log, saveat = y_saveat) #ADD A FLAG TO DISABLE SAVING OPTIONALLY, IT REDUCES ALLOCATIONS

        problem = ODEProblem{true}(f_step!, x₀, (t_start, t_end), params)
        integrator = init(problem, method; callback = scb, save_everystep = false, kwargs...)
        new{C}(integrator, log)
    end
end

Base.getproperty(sys::ContinuousSystem, s::Symbol) = getproperty(sys, Val(s))

Base.getproperty(sys::ContinuousSystem, ::Val{:component}) = sys.integrator.p.comp
Base.getproperty(sys::ContinuousSystem, ::Val{:integrator}) = getfield(sys, :integrator)
Base.getproperty(sys::ContinuousSystem, ::Val{:log}) = getfield(sys, :log)

#forward everything else to the integrator...
Base.getproperty(sys::ContinuousSystem, ::Val{S}) where {S} = getproperty(getfield(sys, :integrator), S)

#...except for x and u (because DiffEqs calls the state u, instead of x)
Base.getproperty(sys::ContinuousSystem, ::Val{:x}) = sys.integrator.u #state vector
Base.getproperty(sys::ContinuousSystem, ::Val{:u}) = sys.integrator.p.u #input vector

Base.getproperty(sys::ContinuousSystem, ::Val{:ẋ}) = sys.integrator.p.ẋ #ẋ cache
Base.getproperty(sys::ContinuousSystem, ::Val{:y}) = sys.integrator.p.y #y cache
Base.getproperty(sys::ContinuousSystem, ::Val{:data}) = sys.integrator.p.data #external data cache

DifferentialEquations.step!(sys::ContinuousSystem, args...) = step!(sys.integrator, args...)
DifferentialEquations.reinit!(sys::ContinuousSystem, args...) = reinit!(sys.integrator, args...)


#mostly useless
################## Generic Component Group ###################

struct ComponentGroup{T, L, C} <: AbstractComponent
    function ComponentGroup(nt::NamedTuple{L, NTuple{N, T}}) where {L, N, T <: AbstractComponent} #Dicts are not ordered, so they won't do
        new{T, L, values(nt)}()
    end
end

x_template(::ComponentGroup{T,L,C}) where {T,L,C} = ComponentVector(NamedTuple{L}(x_template.(C)))
u_template(::ComponentGroup{T,L,C}) where {T,L,C} = ComponentVector(NamedTuple{L}(u_template.(C)))
y_template(::ComponentGroup{T,L,C}) where {T,L,C} = NamedTuple{L}(y_template.(C))
d_template(::ComponentGroup{T,L,C}) where {T,L,C} = d_template(T)

#a single, shared data source is assumed, but this could be easily changed to:
# d_template(g::ComponentGroup{T,L,C}) where {T,L,C} = NamedTuple{L}(d_template.(C))

function f_output!(y::Any, ẋ::Any, x::Any, u::Any, t::Real, data::Any, ::ComponentGroup{T,L,C}) where{T,L,C}
    for (label, component) in zip(L, C)
        #could use map here on (y,ẋ, x, u) but it allocates and harms performance
        f_output!(getproperty(y,label), getproperty(ẋ, label), getproperty(x,label), getproperty(u,label), t, data, component)
    end
end

end