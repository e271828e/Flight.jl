module System

using ComponentArrays

export x0, d0, u0, f_cont!, f_disc!
export AbstractComponent, AbstractSystem, DiscreteSystem, HybridSystem, AlgebraicSystem
export AbstractD, AbstractU, AbstractY
export rplot

# export plotlog

no_extend_error(f::Function, ::Type{S}) where {S} = error("Function $f not implemented for subtype $S or incorrect call signature")

#anything around which we can build a System
abstract type AbstractComponent end #anything that can go in a HybridSystem

#output struct produced by AbstractSystem{C}
abstract type AbstractY{C<:AbstractComponent} end
#discrete state vector struct required by AbstractSystem{C}
abstract type AbstractD{C<:AbstractComponent} end
#input struct required by AbstractSystem{C}
abstract type AbstractU{C<:AbstractComponent} end

x0(::C) where {C<:AbstractComponent} = no_extend_error(x0, C)
d0(::C) where {C<:AbstractComponent} = nothing #systems are not required to have discrete states
u0(::C) where {C<:AbstractComponent} = nothing #sytems are not required to have control inputs

rplot(::AbstractVector{<:Real}, ::AbstractVector{Y}, args...) where {Y<:AbstractY} = no_extend_error(rplot, Y)

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


f_cont!(::S, args...) where {S<:AbstractSystem} = no_extend_error(f_cont!, S)
(f_disc!(::S, args...)::Bool) where {S<:AbstractSystem} = no_extend_error(f_disc!, S)

#this method is free to modify the system's discrete state, control inputs and
#continuous state. if it modifies the latter, it must return true, false
#otherwise. it is dangerous to provide a default fallback for f_disc!, because
#if the intended f_disc! implementation for the System has the wrong interface,
#the dispatch will revert to the fallback, which may not be obvious at all. it
#is safer to force each concrete System that does not require an actual f_disc!
#to implement a trivial f_disc! that returns false






end