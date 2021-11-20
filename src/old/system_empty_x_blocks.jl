module System

using ComponentArrays
import Flight.Plotting: plots

export assemble_x0, init_x0, init_y0, init_u0, init_d0, f_cont!, f_disc!
export AbstractComponent, AbstractSystem, System

no_extend_error(f::Function, ::Type{S}) where {S} = error(
    "Function $f not implemented for type $S or incorrect call signature")
no_extend_warning(f::Function, ::Type{S}) where {S} = println(
    "Warning: Function $f not implemented for type $S or incorrect call signature")

#anything around which we can build a System
abstract type AbstractComponent end #anything that can go in a System

#init_x0 must return either nothing, an AbstractVector{<:Real} or a NT with a
#Union{AbstractVector{<:Real},Nothing}
init_x0(::C) where {C<:AbstractComponent} = nothing
init_y0(::C) where {C<:AbstractComponent} = nothing
init_u0(::C) where {C<:AbstractComponent} = nothing #sytems are not required to have control inputs
init_d0(::C) where {C<:AbstractComponent} = nothing #systems are not required to have discrete states

abstract type AbstractSystem{C<:AbstractComponent} end

#need the C type parameter for dispatch, the rest for type stability
mutable struct System{C, X <: Union{Nothing, AbstractVector{<:Real}},
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

f_cont!(::S, args...) where {S<:AbstractSystem} = no_extend_error(f_cont!, S)
(f_disc!(::S, args...)::Bool) where {S<:AbstractSystem} = no_extend_error(f_disc!, S)

#f_disc! is free to modify a Hybrid system's discrete state, control inputs and
#continuous state. if it modifies the latter, it must return true, false
#otherwise. it is dangerous to provide a default fallback for f_disc!, because
#if the intended f_disc! implementation for the System has the wrong interface,
#the dispatch will revert to the fallback, which may not be obvious at all. it
#is safer to force each concrete System that does not require an actual f_disc!
#to implement a trivial f_disc! that returns false

#when the System constructor for a certain Component is passed no
#parameters for dx and x, it calls the init_x0 method for that Component, which
#may return

assemble_x0(c::AbstractComponent) = assemble_x0(init_x0(c))

assemble_x0(::Nothing) = (s = nothing, ss = nothing)

assemble_x0(x::AbstractVector{<:Real}) = (s = x, ss = nothing)

function assemble_x0(x::ComponentVector)

    x_ss = Dict{Symbol,Any}()

    for id in keys(x)
        x_ss[id] = view(x, id)
    end

    return (s = x, ss = NamedTuple{Tuple(keys(x_ss))}(values(x_ss)))

end

# function assemble_x0(nt::NamedTuple{L,X}
#     ) where {L, X <: NTuple{N, Union{AbstractVector{<:Real}, Nothing}} where {N}}
function assemble_x0(nt::NamedTuple)
    #assembles the System's state vector from the non-empty subsystems' state
    #vectors. returns the System's state vector and the views to be assigned to
    #the subsystems
    x_ss = Dict{Symbol,Any}(pairs(nt))

    #para cada entry del nt, llamamos a assemble_x0, que nos devolvera un x_s y
    #un x_ss

    x_s_blocks = filter(x_ss) do item
        isa(item.second, AbstractVector{<:Real})
    end
    x_s = (!isempty(x_s_blocks) ? ComponentVector{Float64}(x_s_blocks) : nothing)

    for id in keys(x_ss)
        x_ss[id] = (id in keys(x_s_blocks) ? view(x_s, id) : x_ss[id])
    end

    return (s = x_s, ss = NamedTuple{Tuple(keys(x_ss))}(values(x_ss)))

end

function plots(::AbstractVector{<:Real}, ::AbstractVector{T}) where {T}
    no_extend_warning(plots, T) #nothing to plot by default, warn about it
end




end