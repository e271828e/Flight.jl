module System

# using DiffEqCallbacks: SavedValues
export AbstractSystem, X, U, Y, D, f_cont!, f_disc!, plotlog

################## AbstractSystem Interface ###################

abstract type AbstractSystem end

extend_error(::Type{S}) where {S} = error("Method not implemented for subtype $S or incorrect call signature")

X(::S) where {S<:AbstractSystem} = extend_error(S)
U(::S) where {S<:AbstractSystem} = extend_error(S)
Y(::S) where {S<:AbstractSystem} = extend_error(S)
D(::S) where {S<:AbstractSystem} = extend_error(S)

f_cont!(y::Any, ẋ::Any, x::Any, u::Any, t::Any, data::Any, s::S) where {S<:AbstractSystem} = extend_error(S)
(f_disc!(x::Any, u::Any, t::Any, data::Any, s::S)::Bool) where {S<:AbstractSystem} = false #return true if u is modified

#replace this with the appropriate overloads, Plot recipes, whatever
plotlog(log, sys::AbstractSystem) = extend_error(S)

end