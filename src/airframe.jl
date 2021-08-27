module Airframe #this module is really needed, do NOT merge it into aircraft

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack

using Flight.Attitude
using Flight.System

import Flight.System: HybridSystem, x0, d0, u0, f_cont!, f_disc!

export AbstractAirframeComponent, FrameSpec, Wrench
export ACGroup, ACGroupD, ACGroupU, ACGroupY
export translate, get_wr_Ob_b, get_h_Gc_b


############################# FrameSpec ###############################

"""
Specifies a reference frame `fc(Oc, Ɛc)` relative to another `fb(Ob, Ɛb)`

Frame `fc(Oc, Ɛc)` is specified by:
- The position vector from fb's origin Ob to fc's origin Oc, projected on fb's
  axes εb (`r_ObOc_b`)
- The rotation quaternion from fb's axes εb to fc's axes εc (`q_bc`)
"""
Base.@kwdef struct FrameSpec
    r_ObOc_b::SVector{3,Float64} = zeros(SVector{3})
    q_bc::RQuat = RQuat()
end


############################## Wrench #################################

"""
Force and torque combination defined on a concrete reference frame

A `Wrench` is defined on reference frame fc(Oc,εc) when it is applied at its
origin Oc and projected on its axes εc
"""
Base.@kwdef struct Wrench
    F::SVector{3,Float64} = zeros(SVector{3})
    M::SVector{3,Float64} = zeros(SVector{3})
 end


 """
    Base.:+(wr1::Wrench, wr2::Wrench)

Add two compatible `Wrench` instances.

`Wrench` addition should only be performed between compatible `Wrench`
instances, i.e., those defined in the same reference frame
"""
Base.:+(wr1::Wrench, wr2::Wrench) = Wrench(F = wr1.F + wr2.F, M = wr1.M + wr2.M)


"""
    translate(f_bc::FrameSpec, wr_Oc_c::Wrench)

Translate a Wrench from one reference frame to another.

If `f_bc` is a `FrameSpec` specifying frame fc(Oc, εc) relative to fb(Ob, εb),
and `wr_Oc_c` is a `Wrench` defined on fc, then `wr_Ob_b = translate(f_bc,
wr_Oc_c)` is the equivalent `Wrench` defined on fb.

An alternative function call notation is provided:
`f_bc(wr_Oc_c) == translate(f_bc, wr_Oc_c)`
"""
function translate(f_bc::FrameSpec, wr_Oc_c::Wrench)

    @unpack q_bc, r_ObOc_b = f_bc
    F_Oc_c = wr_Oc_c.F
    M_Oc_c = wr_Oc_c.M

    #project onto airframe axes
    F_Oc_b = q_bc(F_Oc_c)
    M_Oc_b = q_bc(M_Oc_c)

    #translate to airframe origin
    F_Ob_b = F_Oc_b
    M_Ob_b = M_Oc_b + r_ObOc_b × F_Oc_b

    return Wrench(F = F_Ob_b, M = M_Ob_b) #wr_Ob_b

end

(f_bc::FrameSpec)(wr_Oc_c::Wrench) = translate(f_bc, wr_Oc_c)



abstract type AbstractAirframeComponent <: AbstractComponent end

#every airframe component must output a struct that implements these methods for
#retrieving wr_Ob_b and h_Gc_b
function get_wr_Ob_b(::T) where {T<:AbstractY{<:AbstractAirframeComponent}}
    error("Method get_wr_Ob_b not implemented for type $T or incorrect call signature")
end
function get_h_Gc_b(::T) where {T<:AbstractY{<:AbstractAirframeComponent}}
    error("Method h_Gc_b not implemented for type $T or incorrect call signature")
end


######################### AircraftComponentGroup ##############################

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

x0(g::ACGroup{T,N,L}) where {T,N,L} = ComponentVector(NamedTuple{L}(x0.(values(g))))
d0(g::ACGroup{T,N,L}) where {T,N,L} = ACGroupD(NamedTuple{L}(d0.(values(g))))
u0(g::ACGroup{T,N,L}) where {T,N,L} = ACGroupU(NamedTuple{L}(u0.(values(g))))

Base.getproperty(y::Union{ACGroupY, ACGroupD, ACGroupU}, s::Symbol) = getproperty(getfield(y,:nt), s)
Base.getindex(y::Union{ACGroupY, ACGroupD, ACGroupU}, s::Symbol) = getindex(getfield(y,:nt), s)

function HybridSystem(g::ACGroup{T,N,L},
    ẋ = x0(g), x = x0(g), d = d0(g), u = u0(g), t = Ref(0.0)) where {T,N,L}

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

@inline @generated function get_wr_Ob_b(y::ACGroupY{Y,N,L}) where {Y<:AbstractY{<:AbstractAirframeComponent},N,L}

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

@inline @generated function get_h_Gc_b(y::ACGroupY{Y,N,L}) where {Y<:AbstractY{<:AbstractAirframeComponent},N,L}

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