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
(f_disc!(x, u, t, s::S, args...)::Bool) where {S<:AbstractSystem} = false #return true if u is modified

#replace this with the appropriate overloads, Plot recipes, whatever
plotlog(log, sys::AbstractSystem) = extend_error(S)

end