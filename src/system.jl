module System

using ComponentArrays
import Flight.Plotting: plots

export get_x0, get_y0, get_u0, get_d0, f_cont!, f_disc!
export AbstractComponent, AbstractSystem, HybridSystem

no_extend_error(f::Function, ::Type{S}) where {S} = error(
    "Function $f not implemented for type $S or incorrect call signature")
no_extend_warning(f::Function, ::Type{S}) where {S} = println(
    "Warning: Function $f not implemented for type $S or incorrect call signature")

#anything around which we can build a System
abstract type AbstractComponent end #anything that can go in a HybridSystem

#every AbstractComponent's get_x0 must return an AbstractVector{<:Real}, even if
#its inherently discrete and its f_cont! does nothing. this is ugly but ensures
#composability of HybridSystems without the hassle of dealing automatically with
#empty state vector blocks, which is magnified by the need to assign views from
#the root System's state vector to each children in its hierarchy
get_x0(::AbstractComponent) = [0.0]
get_y0(::AbstractComponent) = nothing
get_u0(::AbstractComponent) = nothing #sytems are not required to have control inputs
get_d0(::AbstractComponent) = nothing #systems are not required to have discrete states

abstract type AbstractSystem{C<:AbstractComponent} end

#need the C type parameter for dispatch, the rest for type stability
mutable struct HybridSystem{C, X <: AbstractVector{<:Float64},
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

function HybridSystem(c::AbstractComponent, ẋ = get_x0(c), x = get_x0(c),
                    y = get_y0(c), u = get_u0(c), d = get_d0(c), t = Ref(0.0))
    params = c #assign the component descriptor itself as a System parameter
    subsystems = nothing
    HybridSystem{map(typeof, (c, x, y, u, d, params, subsystems))...}(
                                    ẋ, x, y, u, d, t, params, subsystems)
end

#f_disc! is free to modify a Hybrid system's discrete state, control inputs and
#continuous state. if it modifies the latter, it must return true, false
#otherwise. no fallbacks are provided for safety reasons: if the intended
#f_cont! or f_disc! implementations for the System have the wrong interface, the
#dispatch will silently revert to the fallback, which does nothing and may not
#be obvious at all.
f_cont!(::S, args...) where {S<:AbstractSystem} = no_extend_error(f_cont!, S)
(f_disc!(::S, args...)::Bool) where {S<:AbstractSystem} = no_extend_error(f_disc!, S)

function plots(::AbstractVector{<:Real}, ::AbstractVector{T}) where {T}
    no_extend_warning(plots, T) #nothing to plot by default, warn about it
end

"""
#alternative: abstract type MaybeContinuousDynamics end struct
HasContinuousDynamics end struct HasNoContinuousDynamics end
MaybeContinuousDynamics(::Type{<:Any}) = NoContinuousDynamics
MaybeContinuousDynamics(::Type{<:EThruster}) = ContinuousDynamics
# MaybeContinuousDynamics(::Type{<:SimpleISA}) = NoContinuousDynamics

#this is the first method that should be overridden by any HybridSystem that
#declares having continuous dynamics. function f_cont!(s::HybridSystem{C},
args...) where {C<:AbstractComponent} f_cont!(s, MaybeContinuousDynamics(C),
args...) end function f_cont!(::S, ::ContinuousDynamics, args...) where
{S<:HybridSystem} println("Warning: S declared itself as having
ContinuousDynamics but does not implement f_cont!") end f_cont!(::HybridSystem,
::NoContinuousDynamics, args...) = nothing #its OK #this assumes that when some
type calls f_cont! with variable arguments, but the #second one is NOT of type
MaybeContinuousDynamics, the compiler will dispatch #to the first method first.
this method should have been overridden

#we could do the same with MaybeDiscreteDynamics


#what purpose does it serve? it avoids an unwanted and inadverted dispatch to
#the fallback f_cont! or f_disc due to wrong interface in its own
#implementation. obviously, the component is assumed to have actually
#implemented f_cont! and f_disc! and declared itself as having Continous or
#DiscreteDynamics. a component that neither declares ContinuousDynamics nor
#defines f_cont! will receive no warning

#we could go another step and do the same with get_x0, which will also dispatch
#on MaybeContinuousDynamics, and get_d0, on MaybeDiscreteDynamics. then
#MaybeInput, MaybeOutput for get_u0 and get_d0

#however, all this does is force systems with continuous dyanmics to define the
#a MaybeContinuousDynamics method in addition to f_cont!. on exchange, systems
without continuous dynamics are not forced to implement f_cont!. not much to be
gained

"""



end