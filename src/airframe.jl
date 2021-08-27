module Airframe #this module is really needed, do NOT merge it into aircraft

using LinearAlgebra
using StaticArrays
using UnPack

using Flight.Attitude
using Flight.System

export AbstractAirframeComponent, FrameSpec, Wrench
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

@inline @generated function get_wr_Ob_b(y::CGroupY{Y,L}) where {L, Y<:AbstractY{<:AbstractAirframeComponent}}

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

@inline @generated function get_h_Gc_b(y::CGroupY{Y,L}) where {L, Y<:AbstractY{<:AbstractAirframeComponent}}

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