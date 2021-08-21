module System

export AbstractComponent, ContinuousSystem, X, U, Y, f_cont!, f_disc!
# export plotlog

abstract type AbstractComponent end #anything that can go in a ContinuousSystem

no_extend_error(f::Function, ::Type{S}) where {S} = error("Function $f not implemented for subtype $S or incorrect call signature")

X(::C) where {C<:AbstractComponent} = no_extend_error(X, C)
Y(::C) where {C<:AbstractComponent} = no_extend_error(Y, C)
U(::C) where {C<:AbstractComponent} = no_extend_error(U, C)

#need the C type parameter for dispatch, the rest for type stability
struct ContinuousSystem{C<:AbstractComponent, X, Y, U, P, S}
    ẋ::X
    x::X
    y::Y
    u::U
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

f_cont!(::S, args...) where {S<:ContinuousSystem} = no_extend_error(f_cont!, C)
(f_disc!(::S, args...)::Bool) where {S<:ContinuousSystem} = no_extend_error(f_cont!, C)

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

#an output vector is also required, because:
#1) Dealing with Union{Float64,Missing} is slower for some methods used by
#   Kinematics and Dynamics (RQuat() for example)
#2) A continuous dynamical system should have some observable output. if nothing
#   else, the state itself
#this should return true if x was modified, false otherwise.

# it is dangerous to provide a default fallback for f_disc!, because if the
#intended f_disc! implementation for the System has the wrong interface, the
#dispatch will revert to the fallback, which may not be obvious at all. it is
#safer to force each concrete System that does not require an actual f_disc! to
#implement a trivial f_disc! that returns false

# #replace this with the appropriate overloads, Plot recipes, whatever
# plotlog(log, sys::AbstractSystem) = extend_error(S)

end