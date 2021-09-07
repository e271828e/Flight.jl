module System

using ComponentArrays
import Flight.Plotting: plots

export get_x0, get_y0, get_u0, get_d0, f_cont!, f_disc!
export AbstractComponent, AbstractSystem, DiscreteSystem, HybridSystem, AlgebraicSystem
export AbstractD, AbstractU, AbstractY

no_extend_error(f::Function, ::Type{S}) where {S} = error(
    "Function $f not implemented for type $S or incorrect call signature")
no_extend_warning(f::Function, ::Type{S}) where {S} = println(
    "Warning: Function $f not implemented for type $S or incorrect call signature")

#anything around which we can build a System
abstract type AbstractComponent end #anything that can go in a HybridSystem

#output struct produced by AbstractSystem{C}
abstract type AbstractY{C<:AbstractComponent} end
#discrete state vector struct required by AbstractSystem{C}
abstract type AbstractD{C<:AbstractComponent} end
#input struct required by AbstractSystem{C}
abstract type AbstractU{C<:AbstractComponent} end

get_x0(::C) where {C<:AbstractComponent} = no_extend_error(get_x0, C)
get_y0(::C) where {C<:AbstractComponent} = nothing
get_u0(::C) where {C<:AbstractComponent} = nothing #sytems are not required to have control inputs
get_d0(::C) where {C<:AbstractComponent} = nothing #systems are not required to have discrete states

abstract type AbstractSystem{C<:AbstractComponent} end

#need the C type parameter for dispatch, the rest for type stability
mutable struct HybridSystem{C, X <: Union{Nothing, AbstractVector{<:Real}},
                    Y, U, D, P, S} <: AbstractSystem{C}
    ẋ::X #continuous state vector derivative
    x::X #continuous state vector (to be used as a buffer for f_cont! evaluation)
    y::Y #output state
    u::U #control inputs
    d::D #discrete state
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

mutable struct DiscreteSystem{C, Y, U, D, P, S} <: AbstractSystem{C}
    y::Y #control inputs
    u::U #control inputs
    d::D #discrete state vector
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

f_cont!(::S, args...) where {S<:AbstractSystem} = no_extend_error(f_cont!, S)
(f_disc!(::S, args...)::Bool) where {S<:AbstractSystem} = no_extend_error(f_disc!, S)

#f_disc! is free to modify a Hybrid system's discrete state, control inputs and
#continuous state. if it modifies the latter, it must return true, false
#otherwise. it is dangerous to provide a default fallback for f_disc!, because
#if the intended f_disc! implementation for the System has the wrong interface,
#the dispatch will revert to the fallback, which may not be obvious at all. it
#is safer to force each concrete System that does not require an actual f_disc!
#to implement a trivial f_disc! that returns false

function plots(::AbstractVector{<:Real}, ::AbstractVector{T}) where {T<:AbstractY}
    no_extend_warning(plots, T) #nothing to plot by default, warn about it
end




end