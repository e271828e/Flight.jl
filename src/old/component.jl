module Component

using LinearAlgebra
using StaticArrays
using ComponentArrays

using Flight.Attitude
using Flight.Dynamics
using Flight.System
import Flight.System: X, Y, U, f_cont!, f_disc!

export AbstractComponent, AbstractComponentGroup
export get_wr_Ob_b, get_h_Gc_b


################ AbstractComponent Interface ################

abstract type AbstractComponent <: AbstractSystem end

get_wr_Ob_b(y, comp::AbstractComponent, args...) = error("Method get_wr_Ob_b not implemented for subtype $comp or incorrect call signature")
get_h_Gc_b(y, comp::AbstractComponent, args...) = error("Method get_h_Gc_b not implemented for subtype $comp or incorrect call signature")


################ ComponentGroup ###########################

abstract type AbstractComponentGroup{C} <: AbstractComponent end

(G::Type{<:AbstractComponentGroup})(;kwargs...) = G((; kwargs...))
Base.getindex(::AbstractComponentGroup{C}, i::Integer) where {C} = getindex(C, i)
Base.getproperty(::AbstractComponentGroup{C}, i::Symbol) where {C} = getproperty(C, i)
labels(::AbstractComponentGroup{C}) where {C} = keys(C)
components(::AbstractComponentGroup{C}) where {C} = values(C)

X(::AbstractComponentGroup{C}) where {C} = ComponentVector(NamedTuple{keys(C)}(X.(values(C))))
U(::AbstractComponentGroup{C}) where {C} = ComponentVector(NamedTuple{keys(C)}(U.(values(C))))
Y(::AbstractComponentGroup{C}) where {C} = ComponentVector(NamedTuple{keys(C)}(Y.(values(C))))

@inline @generated function f_cont!(y, ẋ, x, u, t, ::AbstractComponentGroup{C}, args...) where {C}
    ex = Expr(:block)
    for label in keys(C)
        label = QuoteNode(label)
        ex_comp = quote
            f_cont!(view(y, $label), view(ẋ, $label), view(x, $label), view(u,$label), t, C[$label], args...)
        end
        push!(ex.args, ex_comp)
    end
    return ex
end

@inline @generated function (f_disc!(x, u, t, ::AbstractComponentGroup{C}, args...)::Bool) where {C}
    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))
    for label in keys(C)
        label = QuoteNode(label)
        ex_comp = quote
            x_mod = x_mod || f_disc!(view(x, $label), view(u,$label), t, C[$label], args...)
        end
        push!(ex.args, ex_comp)
    end
    return ex
end

@inline @generated function get_wr_Ob_b(y, ::AbstractComponentGroup{C}, args...) where {C}

    ex = Expr(:block)
    push!(ex.args, :(wr = Wrench())) #allocate a zero wrench

    for label in keys(C)
        label = QuoteNode(label)
        ex_comp = quote
            #extract and perform in-place broadcasted addition of each
            #component's wrench
            wr .+= get_wr_Ob_b(view(y,$label), C[$label], args...)
        end
        push!(ex.args, ex_comp)
    end
    return ex
end

@inline @generated function get_h_Gc_b(y, ::AbstractComponentGroup{C}, args...) where {C}

    ex = Expr(:block)
    push!(ex.args, :(h = SVector(0., 0., 0.))) #allocate

    for label in keys(C)
        label = QuoteNode(label)
        ex_comp = quote
            h += get_h_Gc_b(view(y,$label), C[$label], args...)
        end
        push!(ex.args, ex_comp)
    end
    return ex
end

end