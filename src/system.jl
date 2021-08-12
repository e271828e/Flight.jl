module System

export AbstractSystem
export X, U, D, typeof_y, output_init, f_output!

################## AbstractSystem Interface ###################

abstract type AbstractSystem end

extend_error(::Type{C}) where {C} = error("Method not implemented for subtype $C or incorrect call signature")

X(::C) where {C<:AbstractSystem} = extend_error(C)
U(::C) where {C<:AbstractSystem} = extend_error(C)
D(::C) where {C<:AbstractSystem} = extend_error(C)
typeof_y(::C) where {C<:AbstractSystem} = extend_error(C)

# X(::Type{C}) where {C<:AbstractSystem} = extend_error(C)
# U(::Type{C}) where {C<:AbstractSystem} = extend_error(C)
# D(::Type{C}) where {C<:AbstractSystem} = extend_error(C)
# typeof_y(::Type{C}) where {C<:AbstractSystem} = extend_error(C)

f_output!(::Any, ::Any, ::Any, ::Real, ::Any, ::C) where {C<:AbstractSystem} = extend_error(C)

#########################################################

function output_init(x::Any, u::Any, t::Real, data::Any, sys::AbstractSystem)
    ẋ = X(sys)
    y = f_output!(ẋ, x, u, t, data, sys)
    return y, ẋ
end

end