module Airframe

using ComponentArrays
using StaticArrays

using Flight.Dynamics
using Flight.System
import Flight.System: HybridSystem, X, D, U, f_cont!, f_disc!

export AbstractAirframeComponent, AirframeComponentGroup
export get_wr_Ob_b, get_h_Gc_b


abstract type AbstractAirframeComponent <: AbstractComponent end

function get_wr_Ob_b(::HybridSystem{C}, args...) where {C<:AbstractAirframeComponent}
    error("Method get_wr_Ob_b not implemented for HybridSystem{$C} or incorrect call signature")
end

function get_h_Gc_b(::HybridSystem{C}, args...) where {C<:AbstractAirframeComponent }
    error("Method get_h_Gc_b not implemented for HybridSystem{$C} or incorrect call signature")
end

######################### AirframeComponentGroup ##############################

struct AirframeComponentGroup{T<:AbstractAirframeComponent,N,L} <: AbstractAirframeComponent
    components::NamedTuple{L, M} where {L, M <: NTuple{N, T} where {N, T<:AbstractAirframeComponent}}
    function AirframeComponentGroup(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, T}} where {N, T<:AbstractAirframeComponent}
        new{T,N,L}(nt)
    end
end
# (G::Type{<:AbstractComponentGroup})(;kwargs...) = G((; kwargs...))
AirframeComponentGroup(;kwargs...) = AirframeComponentGroup((; kwargs...))
Base.length(::AirframeComponentGroup{T,N,L}) where {T,N,L} = N
Base.getindex(g::AirframeComponentGroup, i) = getindex(getfield(g,:components), i)
Base.getproperty(g::AirframeComponentGroup, i::Symbol) = getproperty(getfield(g,:components), i)
Base.keys(::AirframeComponentGroup{T,N,L}) where {T,N,L} = L
Base.values(g::AirframeComponentGroup) = values(getfield(g,:components))

X(g::AirframeComponentGroup{T,N,L}) where {T,N,L} = ComponentVector(NamedTuple{L}(X.(values(g))))
D(g::AirframeComponentGroup{T,N,L}) where {T,N,L} = NamedTuple{L}(D.(values(g)))
U(g::AirframeComponentGroup{T,N,L}) where {T,N,L} = NamedTuple{L}(U.(values(g)))

function HybridSystem(g::AirframeComponentGroup{T,N,L},
    ẋ = X(g), x = X(g), d = D(g), u = U(g), t = Ref(0.0)) where {T,N,L}
    #having L allows us to know the length of g and therefore the number of
    #expressions we need to generate
    s_list = Vector{HybridSystem}()
    for label in L
        s_cmp = HybridSystem(map((λ)->getproperty(λ, label), (g, ẋ, x, d, u))..., t)
        push!(s_list, s_cmp)
    end
    params = nothing #everything is already stored in the subsystem's parameters
    subsystems = NamedTuple{L}(s_list)
    HybridSystem{map(typeof, (g, x, d, u, params, subsystems))...}(ẋ, x, d, u, t, params, subsystems)
end

@inline @generated function f_cont!(sys::HybridSystem{C}, args...) where {C<:AirframeComponentGroup{T,N,L}} where {T,N,L}
    ex_tuple = Expr(:tuple) #construct a tuple around all the individual function calls
    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            f_cont!(sys.subsystems[$label], args...)
        end
        push!(ex_tuple.args, ex_ss)
    end
    # Core.println(ex_tuple)
    #construct a NamedTuple from the labels L and the constructed tuple
    ex_all = Expr(:call, Expr(:curly, NamedTuple, L), ex_tuple)
    return ex_all
end

@inline @generated function (f_disc!(sys::HybridSystem{C}, args...)::Bool) where {C<:AirframeComponentGroup{T,N,L}} where {T,N,L}
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

#sys::HybridSystem{NewComponentGroup} would not work here due to the non-covariance of
#the type system
@inline @generated function get_wr_Ob_b(y::NamedTuple{L}) where {L}

    ex = Expr(:block)
    push!(ex.args, :(wr = Wrench())) #allocate a zero wrench

    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            wr += get_wr_Ob_b(y[$label])
        end
        push!(ex.args, ex_ss)
    end
    return ex
end

@inline @generated function get_h_Gc_b(y::NamedTuple{L}) where {L}

    ex = Expr(:block)
    push!(ex.args, :(h = SVector(0., 0., 0.))) #allocate

    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            h += get_h_Gc_b(y[$label])
        end
        push!(ex.args, ex_ss)
    end
    return ex
end

end