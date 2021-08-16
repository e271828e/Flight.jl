module System

# using DiffEqCallbacks: SavedValues
export AbstractSystem, X, U, Y, D, f_cont!, f_disc!, plotlog

################## AbstractSystem Interface ###################

abstract type AbstractSystem end

extend_error(::Type{S}) where {S} = error("Method not implemented for subtype $S or incorrect call signature")

X(::AbstractSystem) = extend_error(S)
U(::AbstractSystem) = extend_error(S)
Y(::AbstractSystem) = extend_error(S)
D(::AbstractSystem) = extend_error(S)

f_cont!(y::Any, ẋ::Any, x::Any, u::Any, t::Any, data::Any, s::AbstractSystem) = extend_error(S)
(f_disc!(x::Any, u::Any, t::Any, data::Any, s::AbstractSystem)::Bool) = false #return true if u is modified

#replace this with the appropriate overloads, Plot recipes, whatever
plotlog(log, sys::AbstractSystem) = extend_error(S)

end