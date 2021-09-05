module Airframe #this module is really needed, do NOT merge it into aircraft

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack

using Flight.Attitude
using Flight.System
using Flight.Dynamics

import Flight.System: HybridSystem, get_x0, get_d0, get_u0, f_cont!, f_disc!
import Flight.Dynamics: get_wr_b, get_hr_b

export ACGroup, ACGroupD, ACGroupU, ACGroupY
export AbstractAirframeComponent


abstract type AbstractAirframeComponent <: AbstractComponent end

######################### AirframeComponentGroup #############################

#we must keep N as a type parameter, because it's left open in the components
#type declaration!
struct ACGroup{T<:AbstractAirframeComponent,N,L} <: AbstractAirframeComponent
    components::NamedTuple{L, M} where {L, M <: NTuple{N, T}}
    function ACGroup(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, T}} where {N, T<:AbstractAirframeComponent}
        new{T,N,L}(nt)
    end
end

ACGroup(;kwargs...) = ACGroup((; kwargs...))

Base.length(::ACGroup{T,N,L}) where {T,N,L} = N
Base.getindex(g::ACGroup, i) = getindex(getfield(g,:components), i)
Base.getproperty(g::ACGroup, i::Symbol) = getproperty(getfield(g,:components), i)
Base.keys(::ACGroup{T,N,L}) where {T,N,L} = L
Base.values(g::ACGroup) = values(getfield(g,:components))


struct ACGroupU{U<:AbstractU,N,L} <: AbstractU{ACGroup}
    nt::NamedTuple{L, NTuple{N,U}}
    function ACGroupU(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, U}} where {N, U}
        new{U,N,L}(nt)
    end
end

struct ACGroupD{D<:AbstractD,N,L} <: AbstractD{ACGroup}
    nt::NamedTuple{L, NTuple{N,D}}
    function ACGroupD(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, D}} where {N, D}
        new{D,N,L}(nt)
    end
end

struct ACGroupY{Y<:AbstractY,N,L} <: AbstractY{ACGroup}
    nt::NamedTuple{L, NTuple{N,Y}}
    function ACGroupY(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, Y}} where {N, Y}
        new{Y,N,L}(nt)
    end
end

get_x0(g::ACGroup{T,N,L}) where {T,N,L} = ComponentVector(NamedTuple{L}(get_x0.(values(g))))
get_d0(g::ACGroup{T,N,L}) where {T,N,L} = ACGroupD(NamedTuple{L}(get_d0.(values(g))))
get_u0(g::ACGroup{T,N,L}) where {T,N,L} = ACGroupU(NamedTuple{L}(get_u0.(values(g))))

Base.getproperty(y::Union{ACGroupY, ACGroupD, ACGroupU}, s::Symbol) = getproperty(getfield(y,:nt), s)
Base.getindex(y::Union{ACGroupY, ACGroupD, ACGroupU}, s::Symbol) = getindex(getfield(y,:nt), s)

function HybridSystem(g::ACGroup{T,N,L},
    ẋ = get_x0(g), x = get_x0(g), d = get_d0(g), u = get_u0(g), t = Ref(0.0)) where {T,N,L}

    s_list = Vector{HybridSystem}()
    for label in L
        s_cmp = HybridSystem(map((λ)->getproperty(λ, label), (g, ẋ, x, d, u))..., t)
        push!(s_list, s_cmp)
    end

    params = nothing #everything is already stored in the subsystem's parameters
    subsystems = NamedTuple{L}(s_list)

    HybridSystem{map(typeof, (g, x, d, u, params, subsystems))...}(ẋ, x, d, u, t, params, subsystems)
end

@inline @generated function f_cont!(sys::HybridSystem{C}, args...) where {C<:ACGroup{T,N,L}} where {T,N,L}
    ex_tuple = Expr(:tuple) #construct a tuple around all the individual function calls
    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            f_cont!(sys.subsystems[$label], args...)
        end
        push!(ex_tuple.args, ex_ss)
    end
    #construct a NamedTuple from the labels L and the constructed tuple
    ex_nt = Expr(:call, Expr(:curly, NamedTuple, L), ex_tuple)
    ex_out = Expr(:call, ACGroupY, ex_nt)
    return ex_out
end
#

@inline @generated function (f_disc!(sys::HybridSystem{C}, args...)::Bool) where {C<:ACGroup{T,N,L}} where {T,N,L}
    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))
    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            x_mod = x_mod || f_disc!(sys.subsystems[$label], args...)
        end
        push!(ex.args, ex_ss)
    end
    return ex
end

@inline @generated function get_wr_b(y::ACGroupY{Y,N,L}) where {Y<:AbstractY{<:AbstractAirframeComponent},N,L}

    ex = Expr(:block)
    push!(ex.args, :(wr = Wrench())) #allocate a zero wrench

    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            wr += get_wr_b(y[$label])
        end
        push!(ex.args, ex_ss)
    end
    return ex
end

@inline @generated function get_hr_b(y::ACGroupY{Y,N,L}) where {Y<:AbstractY{<:AbstractAirframeComponent},N,L}

    ex = Expr(:block)
    push!(ex.args, :(h = SVector(0., 0., 0.))) #allocate

    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            h += get_hr_b(y[$label])
        end
        push!(ex.args, ex_ss)
    end
    return ex
end



end