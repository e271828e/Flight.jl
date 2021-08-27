module System

using ComponentArrays

export x0, d0, u0, f_cont!, f_disc!
export AbstractComponent, AbstractSystem, DiscreteSystem, HybridSystem, AlgebraicSystem
export AbstractD, AbstractU, AbstractY
export ComponentGroup, CGroupD, CGroupU, CGroupY

# export plotlog


#pregunta: necesito realmente subsystems? que me aporta? por ejemplo, supongamos
#que me defino un tricycle landing gear. al final me tendre que hacer una
#f_cont! que considere la estructura especifica del state vector y elinput
#vector y reparta sus elementos a las distintas patas. si, pero en el momento en
#que quiera llamar a la f_cont! de cada una de esas patas, al final esa f_cont!
#acepta como argumento un sys. en algun sitio tendre que almacenar esos sys si
#quiero llamar a sus f_cont!s desde la f_cont! del parent system. por tanto si,
#si que aporta tener subsystems. de acuerdo. pero entonces eso significa que
#el parent system no debe almacenar en su field params directamente el
#AbstractComponent que lo genero, porque es de suponer que los distintos
#sub-fields de ese AbstractComponents habran sido repartidos por los diferentes
#subsystems. de hecho, ese es el caso por excelencia de ComponentGroup

no_extend_error(f::Function, ::Type{S}) where {S} = error("Function $f not implemented for subtype $S or incorrect call signature")

#anything around which we can build a System
abstract type AbstractComponent end #anything that can go in a HybridSystem

abstract type AbstractSystem{C<:AbstractComponent} end

#need the C type parameter for dispatch, the rest for type stability
struct HybridSystem{C, X, D, U, P, S} <: AbstractSystem{C}
    ẋ::X #continuous state vector derivative
    x::X #continuous state vector (to be used as a buffer for f_cont! evaluation)
    d::D #discrete state vector
    u::U #control inputs
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

struct DiscreteSystem{C, D, U, P, S} <: AbstractSystem{C}
    d::D #discrete state vector
    u::U #control inputs
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

struct AlgebraicSystem{C, U, P, S} <: AbstractSystem{C}
    u::U #control inputs
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

#output struct produced by AbstractSystem{C}
abstract type AbstractY{C<:AbstractComponent} end
#discrete state vector struct required by AbstractSystem{C}
abstract type AbstractD{C<:AbstractComponent} end
#input struct required by AbstractSystem{C}
abstract type AbstractU{C<:AbstractComponent} end

x0(::C) where {C<:AbstractComponent} = no_extend_error(x0, C)
d0(::C) where {C<:AbstractComponent} = nothing #systems are not required to have discrete states
u0(::C) where {C<:AbstractComponent} = nothing #sytems are not required to have control inputs

f_cont!(::S, args...) where {S<:AbstractSystem} = no_extend_error(f_cont!, S)
(f_disc!(::S, args...)::Bool) where {S<:AbstractSystem} = no_extend_error(f_disc!, S)

#this method is free to modify the system's discrete state, control inputs and
#continuous state. if it modifies the latter, it must return true, false
#otherwise. it is dangerous to provide a default fallback for f_disc!, because
#if the intended f_disc! implementation for the System has the wrong interface,
#the dispatch will revert to the fallback, which may not be obvious at all. it
#is safer to force each concrete System that does not require an actual f_disc!
#to implement a trivial f_disc! that returns false

######################### ComponentGroup ##############################

struct ComponentGroup{T<:AbstractComponent,L} <: AbstractComponent
    components::NamedTuple{L, M} where {L, M <: NTuple{N, T} where {N, T<:AbstractComponent}}
    function ComponentGroup(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, T}} where {N, T<:AbstractComponent}
        new{T,L}(nt)
    end
end

ComponentGroup(;kwargs...) = ComponentGroup((; kwargs...))

Base.length(::ComponentGroup{T,L}) where {T,L} = length(L)
Base.getindex(g::ComponentGroup, i) = getindex(getfield(g,:components), i)
Base.getproperty(g::ComponentGroup, i::Symbol) = getproperty(getfield(g,:components), i)
Base.keys(::ComponentGroup{T,L}) where {T,L} = L
Base.values(g::ComponentGroup) = values(getfield(g,:components))


struct CGroupU{U<:AbstractU,L} <: AbstractU{ComponentGroup}
    nt::NamedTuple{L, NTuple{N,U}} where {N}
    function CGroupU(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, U}} where {N, U<:AbstractU}
        new{U,L}(nt)
    end
end

struct CGroupD{D<:AbstractD,L} <: AbstractD{ComponentGroup}
    nt::NamedTuple{L, NTuple{N,D}} where {N}
    function CGroupD(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, D}} where {N, D<:AbstractD}
        new{D,L}(nt)
    end
end

struct CGroupY{Y<:AbstractY,L} <: AbstractY{ComponentGroup}
    nt::NamedTuple{L, NTuple{N,Y}} where {N}
    function CGroupY(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, Y}} where {N, Y<:AbstractY}
        new{Y,L}(nt)
    end
end

x0(g::ComponentGroup{T,L}) where {T,L} = ComponentVector(NamedTuple{L}(x0.(values(g))))
d0(g::ComponentGroup{T,L}) where {T,L} = CGroupD(NamedTuple{L}(d0.(values(g))))
u0(g::ComponentGroup{T,L}) where {T,L} = CGroupU(NamedTuple{L}(u0.(values(g))))

Base.getproperty(y::Union{CGroupY, CGroupD, CGroupU}, s::Symbol) = getproperty(getfield(y,:nt), s)
Base.getindex(y::Union{CGroupY, CGroupD, CGroupU}, s::Symbol) = getindex(getfield(y,:nt), s)

function HybridSystem(g::ComponentGroup{T,L},
    ẋ = x0(g), x = x0(g), d = d0(g), u = u0(g), t = Ref(0.0)) where {T,L}

    s_list = Vector{HybridSystem}()
    for label in L
        s_cmp = HybridSystem(map((λ)->getproperty(λ, label), (g, ẋ, x, d, u))..., t)
        push!(s_list, s_cmp)
    end

    params = nothing #everything is already stored in the subsystem's parameters
    subsystems = NamedTuple{L}(s_list)

    HybridSystem{map(typeof, (g, x, d, u, params, subsystems))...}(ẋ, x, d, u, t, params, subsystems)
end

@inline @generated function f_cont!(sys::HybridSystem{C}, args...) where {C<:ComponentGroup{T,L}} where {T,L}
    ex_tuple = Expr(:tuple) #construct a tuple around all the individual function calls
    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            f_cont!(sys.subsystems[$label], args...)
        end
        push!(ex_tuple.args, ex_ss)
    end
    #construct a NamedTuple from the labels L and the constructed tuple
    ex_nt = Expr(:call, Expr(:curly, NamedTuple, L), ex_tuple)
    ex_out = Expr(:call, CGroupY, ex_nt)
    return ex_out
end
#

@inline @generated function (f_disc!(sys::HybridSystem{C}, args...)::Bool) where {C<:ComponentGroup{T,L}} where {T,L}
    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))
    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            x_mod = x_mod || f_disc!(sys.subsystems[$label], args...)
        end
        push!(ex.args, ex_ss)
    end
    return ex
end

# #replace this with the appropriate overloads, Plot recipes, whatever
# plotlog(log, sys::AbstractSystem) = extend_error(S)




end