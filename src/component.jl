module Component

using LinearAlgebra
using StaticArrays
using ComponentArrays

using Flight.Attitude
using Flight.System
import Flight.System: X, Y, U, D, f_cont!, f_disc!

export WrenchCV, Wrench, Frame, AbstractComponent, ComponentGroup
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

get_wr_Ob_b(::Any, comp::AbstractComponent) = error("Method not implemented for subtype $comp or incorrect call signature")
get_h_Gc_b(::Any, comp::AbstractComponent) = error("Method not implemented for subtype $comp or incorrect call signature")


################ ComponentGroup ###########################

struct ComponentGroup{T <: AbstractComponent, C} <: AbstractComponent
    function ComponentGroup(nt::NamedTuple{L, NTuple{N, T}}) where {L, N, T <: AbstractComponent} #Dicts are not ordered, so they won't do
        new{T,nt}()
    end
end
ComponentGroup(;kwargs...) = ComponentGroup((;kwargs...))
Base.getindex(::ComponentGroup{T, C}, i::Integer) where {T, C} = getindex(C, i)
Base.getproperty(::ComponentGroup{T, C}, i::Symbol) where {T, C} = getproperty(C, i)
labels(::ComponentGroup{T, C}) where {T, C} = keys(C)
components(::ComponentGroup{T, C}) where {T, C} = values(C)

X(::ComponentGroup{T, C}) where {T, C} = ComponentVector(NamedTuple{keys(C)}(X.(values(C))))
U(::ComponentGroup{T, C}) where {T, C} = ComponentVector(NamedTuple{keys(C)}(U.(values(C))))
Y(::ComponentGroup{T, C}) where {T, C} = ComponentVector(NamedTuple{keys(C)}(Y.(values(C))))
D(::ComponentGroup{T, C}) where {T, C} = D(C[1]) #assume all components use the same external data sources
# D(::ComponentGroup{T, C}) where {T, C} = NamedTuple{L}(D.(C))

@inline @generated function f_cont!(y::Any, ẋ::Any, x::Any, u::Any, t::Any,
                                    data::Any, ::ComponentGroup{T,C}) where {T,C}
    ex = Expr(:block)
    for label in keys(C)
        label = QuoteNode(label)
        ex_comp = quote
            y_cmp = @view y[$label]; ẋ_cmp = @view ẋ[$label]
            x_cmp = @view x[$label]; u_cmp = @view u[$label]
            f_cont!(y_cmp, ẋ_cmp, x_cmp, u_cmp, t, data, C[$label])
        end
        push!(ex.args, ex_comp)
    end
    return ex
end

@inline @generated function (f_disc!(x::Any, u::Any, t::Any,
                                    data::Any, ::ComponentGroup{T,C})::Bool) where {T,C}
    ex = Expr(:block)
    push!(ex.args, :(modified_x = false))
    for label in keys(C)
        label = QuoteNode(label)
        ex_comp = quote
            x_cmp = @view x[$label]; u_cmp = @view u[$label]
            modified_x_cmp = f_disc!(x_cmp, u_cmp, t, data, C[$label])
            modified_x = modified_x || modified_x_cmp
        end
        push!(ex.args, ex_comp)
    end
    return ex
end

@inline @generated function get_wr_Ob_b(y::Any, ::ComponentGroup{T,C}) where {T,C}

    ex = Expr(:block)
    push!(ex.args, :(wr = Wrench())) #allocate a zero wrench

    for label in keys(C)
        label = QuoteNode(label)
        ex_comp = quote
            #extract and perform in-place broadcasted addition of each
            #component's wrench
            wr .+= get_wr_Ob_b(view(y,$label), C[$label])
        end
        push!(ex.args, ex_comp)
    end
    return ex
end

@inline @generated function get_h_Gc_b(y::Any, ::ComponentGroup{T,C}) where {T,C}

    ex = Expr(:block)
    push!(ex.args, :(h = SVector(0., 0., 0.))) #allocate

    for label in keys(C)
        label = QuoteNode(label)
        ex_comp = quote
            h += get_h_Gc_b(view(y,$label), C[$label])
        end
        push!(ex.args, ex_comp)
    end
    return ex
end

end