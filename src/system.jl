module System

# using DiffEqCallbacks: SavedValues
export AbstractSystem, X, U, Y, D, f_cont!, f_disc!, plotlog

################## AbstractSystem Interface ###################

abstract type AbstractSystem end

extend_error(::Type{S}) where {S} = error("Method not implemented for subtype $S or incorrect call signature")

X(::S) where {S<:AbstractSystem} = extend_error(S)
U(::S) where {S<:AbstractSystem} = extend_error(S)
Y(::S) where {S<:AbstractSystem} = extend_error(S)

f_cont!(y, ẋ, x, u, t, s::S, args...) where {S<:AbstractSystem} = extend_error(S)
(f_disc!(x, u, t, s::S, args...)::Bool) where {S<:AbstractSystem} = extend_error(S)
#this should return true if x was modified, false otherwise.

# it is dangerous to provide a default fallback for f_disc!, because if the
#intended f_disc! implementation for the System has the wrong interface, the
#dispatch will revert to the fallback, which may not be obvious at all. it is
#safer to force each concrete System that does not require an actual f_disc! to
#implement a trivial f_disc! that returns false

#replace this with the appropriate overloads, Plot recipes, whatever
plotlog(log, sys::AbstractSystem) = extend_error(S)

end