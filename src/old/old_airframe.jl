module Airframe #this module is really needed, do NOT merge it into aircraft

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.Attitude
using Flight.Dynamics
using Flight.ModelingTools
import Flight.ModelingTools: System, get_x0, get_y0, get_u0, get_d0, f_cont!, f_disc!

using Flight.Plotting
import Flight.Plotting: plots

export AbstractAirframeComponent, NullAirframeComponent, AbstractAirframeNode, AirframeGroup
export get_wr_b, get_hr_b


abstract type AbstractAirframeComponent <: AbstractComponent end

function get_wr_b(::T) where {T<:System{<:AbstractAirframeComponent}}
    error("Method get_wr_b not implemented for type $T or incorrect call signature")
end
function get_hr_b(::T) where {T<:System{<:AbstractAirframeComponent}}
    error("Method hr_b not implemented for type $T or incorrect call signature")
end

######################### NullAirframeComponent ############################

struct NullAirframeComponent <: AbstractAirframeComponent end

@inline get_wr_b(::System{NullAirframeComponent}) = Wrench()
@inline get_hr_b(::System{NullAirframeComponent}) = zeros(SVector{3})

@inline f_cont!(::System{NullAirframeComponent}, args...) = nothing
@inline (f_disc!(::System{NullAirframeComponent}, args...)::Bool) = false

######################### AirframeGroup #############################

#AirframeGroup must only group identical components, otherwise the arguments to
#f_cont! and f_disc! would generally be different for each one

#must keep N as a type parameter, because it's left open in the components
#type declaration
struct AirframeGroup{T<:AbstractAirframeComponent,N,L} <: AbstractAirframeComponent
    components::NamedTuple{L, M} where {L, M <: NTuple{N, T}}
    function AirframeGroup(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, T}} where {N, T<:AbstractAirframeComponent}
        new{T,N,L}(nt)
    end
end

AirframeGroup(;kwargs...) = AirframeGroup((; kwargs...))

Base.length(::AirframeGroup{T,N,L}) where {T,N,L} = N
Base.getindex(g::AirframeGroup, i) = getindex(getfield(g,:components), i)
Base.getproperty(g::AirframeGroup, i::Symbol) = getproperty(getfield(g,:components), i)
Base.keys(::AirframeGroup{T,N,L}) where {T,N,L} = L
Base.values(g::AirframeGroup) = values(getfield(g,:components))

get_x0(g::AirframeGroup{T,N,L}) where {T,N,L} = NamedTuple{L}(get_x0.(values(g))) |> ComponentVector
get_u0(g::AirframeGroup{T,N,L}) where {T,N,L} = NamedTuple{L}(get_u0.(values(g)))
get_y0(g::AirframeGroup{T,N,L}) where {T,N,L} = NamedTuple{L}(get_y0.(values(g)))
get_d0(g::AirframeGroup{T,N,L}) where {T,N,L} = NamedTuple{L}(get_d0.(values(g)))

function System(g::AirframeGroup{T,N,L},
                    ẋ = get_x0(g), x = get_x0(g), y = get_y0(g), u = get_u0(g),
                    d = get_d0(g), t = Ref(0.0)) where {T,N,L}

    ss_list = Vector{System}()
    for label in L
        s_cmp = System(map((λ)->getproperty(λ, label), (g, ẋ, x, y, u, d))..., t)
        push!(ss_list, s_cmp)
    end

    params = nothing #everything is already stored in the subsystem's parameters
    subsystems = NamedTuple{L}(ss_list)

    System{map(typeof, (g, x, y, u, d, params, subsystems))...}(ẋ, x, y, u, d, t, params, subsystems)
end

@inline @generated function f_cont!(sys::System{C}, args...
    ) where {C<:AirframeGroup{T,N,L}} where {T <: AbstractAirframeComponent,N,L}

    ex_main = Expr(:block)

    #call f_cont! on each subsystem
    ex_calls = Expr(:block)
    for label in L
        push!(ex_calls.args,
            :(f_cont!(sys.subsystems[$(QuoteNode(label))], args...)))
    end

    #retrieve the y from each subsystem and build a tuple with them
    ex_tuple = Expr(:tuple)
    for label in L
        push!(ex_tuple.args,
            :(sys.subsystems[$(QuoteNode(label))].y))
    end

    #build a NamedTuple from the subsystem's labels and the constructed tuple
    ex_y = Expr(:call, Expr(:curly, NamedTuple, L), ex_tuple)

    #assign the result to the parent system's y
    ex_assign_y = Expr(:(=), :(sys.y), ex_y)

    #pack everything into the main block expression
    push!(ex_main.args, ex_calls)
    push!(ex_main.args, ex_assign_y)
    push!(ex_main.args, :(return nothing))

    return ex_main

end


@inline @generated function (f_disc!(sys::System{C}, args...)::Bool
    ) where {C<:AirframeGroup{T,N,L}} where {T <:AbstractAirframeComponent,N,L}

    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))
    for label in L
        #we need all f_disc! calls executed, so | must be used instead of ||
        push!(ex.args,
            :(x_mod = x_mod | f_disc!(sys.subsystems[$(QuoteNode(label))], args...)))
    end
    return ex

end

@inline @generated function get_wr_b(sys::System{C}
    ) where {C<:AirframeGroup{T,N,L}} where {T <: AbstractAirframeComponent,N,L}

    ex = Expr(:block)
    push!(ex.args, :(wr = Wrench())) #allocate a zero wrench
    for label in L
        push!(ex.args,
            :(wr += get_wr_b(sys.subsystems[$(QuoteNode(label))])))
    end
    return ex

end

@inline @generated function get_hr_b(sys::System{C}
    ) where {C<:AirframeGroup{T,N,L}} where {T <: AbstractAirframeComponent,N,L}

    ex = Expr(:block)
    push!(ex.args, :(h = SVector(0., 0., 0.))) #allocate
    for label in L
        push!(ex.args,
            :(h += get_hr_b(sys.subsystems[$(QuoteNode(label))])))
    end
    return ex

end

end #module