module Airframe

using ComponentArrays

using Flight.System
import Flight.System: ContinuousSystem, X, Y, U, f_cont!, f_disc!

export AbstractAirframeComponent, AirframeComponentGroup
export get_wr_Ob_b, get_h_Gc_b


abstract type AbstractAirframeComponent <: AbstractComponent end

function get_wr_Ob_b(::ContinuousSystem{C}, args...) where {C<:AbstractAirframeComponent}
    error("Method get_wr_Ob_b not implemented for ContinuousSystem{$C} or incorrect call signature")
end

function get_h_Gc_b(::ContinuousSystem{C}, args...) where {C<:AbstractAirframeComponent }
    error("Method get_h_Gc_b not implemented for ContinuousSystem{$C} or incorrect call signature")
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
Y(g::AirframeComponentGroup{T,N,L}) where {T,N,L} = ComponentVector(NamedTuple{L}(Y.(values(g))))
U(g::AirframeComponentGroup{T,N,L}) where {T,N,L} = NamedTuple{L}(U.(values(g)))

function ContinuousSystem(g::AirframeComponentGroup{T,N,L},
    ẋ = X(g), x = X(g), y = Y(g), u = U(g), t = Ref(0.0)) where {T,N,L}
    #having L allows us to know the length of g and therefore the number of
    #expressions we need to generate
    s_list = Vector{ContinuousSystem}()
    for label in L
        s_cmp = ContinuousSystem(map((λ)->getproperty(λ, label), (g, ẋ, x, y, u))..., t)
        push!(s_list, s_cmp)
    end
    params = nothing #everything is already stored in the subsystem's parameters
    subsystems = NamedTuple{L}(s_list)
    ContinuousSystem{map(typeof, (g, x, y, u, params, subsystems))...}(ẋ, x, y, u, t, params, subsystems)
end

@inline @generated function f_cont!(sys::ContinuousSystem{C}, args...) where {C<:AirframeComponentGroup{T,N,L}} where {T,N,L}
    ex = Expr(:block)
    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            f_cont!(sys.subsystems[$label], args...)
        end
        push!(ex.args, ex_ss)
    end
    return ex
end

@inline @generated function (f_disc!(sys::ContinuousSystem{C}, args...)::Bool) where {C<:AirframeComponentGroup{T,N,L}} where {T,N,L}
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

#sys::ContinuousSystem{NewComponentGroup} would not work here due to the non-covariance of
#the type system
@inline @generated function get_wr_Ob_b(sys::ContinuousSystem{C}, args...) where {C<:AirframeComponentGroup{T,N,L}} where {T,N,L}

    ex = Expr(:block)
    push!(ex.args, :(wr = Wrench())) #allocate a zero wrench

    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            wr .+= get_wr_Ob_b(sys.subsystems[$label], args...)
        end
        push!(ex.args, ex_ss)
    end
    return ex
end

@inline @generated function get_h_Gc_b(sys::ContinuousSystem{C}, args...) where {C<:AirframeComponentGroup{T,N,L}} where {T,N,L}
    ex = Expr(:block)
    push!(ex.args, :(h = SVector(0., 0., 0.))) #allocate

    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            h += get_h_Gc_b(sys.subsystems[$label], args...)
        end
        push!(ex.args, ex_ss)
    end
    return ex
end

end