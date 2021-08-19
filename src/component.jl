module Component

using LinearAlgebra
using StaticArrays
using ComponentArrays

using Flight.Attitude
using Flight.System
import Flight.System: X, Y, U, f_cont!, f_disc!

export WrenchCV, Wrench, Frame, AbstractComponent, AbstractComponentGroup
export get_wr_Ob_b, get_h_Gc_b


################# Wrench ########################

const WrenchAxes = getaxes(ComponentVector(F = zeros(3), M = zeros(3)))
const WrenchCV{D} = ComponentVector{Float64, D, typeof(WrenchAxes)} where {D <: AbstractVector{Float64}}
const Wrench(v::AbstractVector{Float64}) = (@assert length(v) == 6; ComponentVector(v, WrenchAxes))
function Wrench(; F = SVector(0.0,0,0), M = SVector(0.0,0,0))
    wr = ComponentVector{Float64}(undef, WrenchAxes)
    wr.F = F; wr.M = M
    return wr
end


####################### Frame ###############

"""
#Specifies a local Frame fc(Oc, Ɛc) relative to the airframe reference
frame fb(Ob, Ɛb) by:
#a) the position vector of the local frame origin Oc relative to the reference
#frame origin Ob, projected in the reference frame axes
# b) the attitude of the local frame axes relative to the reference
#frame axes, given by rotation quaternion q_bc
"""
Base.@kwdef struct Frame
    r_ObOc_b::SVector{3,Float64} = zeros(SVector{3})
    q_bc::RQuat = RQuat()
end

"""
Translate a Wrench specified on a local Frame fc(Oc, εc) to the
airframe reference frame fb(Ob, εb) given the relative Frame
specification f_bc
"""

function Base.:*(f_bc::Frame, wr_Oc_c::WrenchCV)

    F_Oc_c = wr_Oc_c.F
    M_Oc_c = wr_Oc_c.M

    #project on the reference axes
    F_Oc_b = f_bc.q_bc * F_Oc_c
    M_Oc_b = f_bc.q_bc * M_Oc_c

    #translate them to airframe origin
    F_Ob_b = F_Oc_b
    M_Ob_b = M_Oc_b + f_bc.r_ObOc_b × F_Oc_b
    Wrench(F = F_Ob_b, M = M_Ob_b) #wr_Ob_b

end

################ AbstractComponent Interface ################

abstract type AbstractComponent <: AbstractSystem end

get_wr_Ob_b(y, comp::AbstractComponent, args...) = error("Method not implemented for subtype $comp or incorrect call signature")
get_h_Gc_b(y, comp::AbstractComponent, args...) = error("Method not implemented for subtype $comp or incorrect call signature")


################ ComponentGroup ###########################

abstract type AbstractComponentGroup{C} <: AbstractComponent end

(G::Type{<:AbstractComponentGroup})(;kwargs...) = G((; kwargs...))
Base.getindex(::AbstractComponentGroup{C}, i::Integer) where {C} = getindex(C, i)
Base.getproperty(::AbstractComponentGroup{C}, i::Symbol) where {C} = getproperty(C, i)
labels(::AbstractComponentGroup{C}) where {C} = keys(C)
components(::AbstractComponentGroup{C}) where {C} = values(C)

#ONLY CREATE BLOCKS FOR THOSE THAT RETURN != nothing
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