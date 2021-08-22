module System

export AbstractComponent, HybridSystem, X, U, f_cont!, f_disc!
# export plotlog

abstract type AbstractComponent end #anything that can go in a HybridSystem

no_extend_error(f::Function, ::Type{S}) where {S} = error("Function $f not implemented for subtype $S or incorrect call signature")

X(::C) where {C<:AbstractComponent} = no_extend_error(X, C)
D(::C) where {C<:AbstractComponent} = nothing #systems are not required to have discrete states
U(::C) where {C<:AbstractComponent} = nothing #sytems are not required to have control inputs

#need the C type parameter for dispatch, the rest for type stability
struct HybridSystem{C<:AbstractComponent, X, D, U, P, S}
    ẋ::X #continuous state vector derivative
    x::X #continuous state vector (to be used as a buffer for f_cont! evaluation)
    d::D #discrete state vector
    u::U #control inputs
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

f_cont!(::S, args...) where {S<:HybridSystem} = no_extend_error(f_cont!, S)

#this method is free to modify the system's discrete state, control inputs and
#continuous state. if it modifies the latter, it must return true, false
#otherwise. it is dangerous to provide a default fallback for f_disc!, because
#if the intended f_disc! implementation for the System has the wrong interface,
#the dispatch will revert to the fallback, which may not be obvious at all. it
#is safer to force each concrete System that does not require an actual f_disc!
#to implement a trivial f_disc! that returns false
(f_disc!(::S, args...)::Bool) where {S<:HybridSystem} = no_extend_error(f_cont!, S)


# #replace this with the appropriate overloads, Plot recipes, whatever
# plotlog(log, sys::AbstractSystem) = extend_error(S)

#the AbstractSystem requires defining at least a continuous state vector, which
#should be a ComponentVector returned by method X. to main reasons:
#1) DifferentialEquations seems to fail when the state vector is an array of
#   Union{Float64, Missing}. on the other hand, since the input and output
#   vectors are parameters from the integrator's perspective, it does not handle
#   them explicitly and they are no problem.
#2) If a continuous dynamical system has no continuous states, it is not a
#   continuous dynamical system, so it should not be defined at such. if it is
#   part of a parent system which does have states, it can be handled within the
#   parent System.



end