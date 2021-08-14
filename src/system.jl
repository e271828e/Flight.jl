module System

export AbstractSystem, X, U, Y, D, f_output!

################## AbstractSystem Interface ###################

abstract type AbstractSystem end

extend_error(::Type{C}) where {C} = error("Method not implemented for subtype $C or incorrect call signature")

X(::C) where {C<:AbstractSystem} = extend_error(C)
U(::C) where {C<:AbstractSystem} = extend_error(C)
Y(::C) where {C<:AbstractSystem} = extend_error(C)
D(::C) where {C<:AbstractSystem} = extend_error(C)

f_output!(::Any, ::Any, ::Any, ::Any, ::Real, ::Any, ::C) where {C<:AbstractSystem} = extend_error(C)

end