module Components

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.Modeling
using Flight.Plotting
using Flight.Attitude
using Flight.Dynamics

import Flight.Plotting: plots

export HasMass, HasNoMass, get_mp_b
export HasWrench, HasNoWrench, get_wr_b
export HasAngularMomentum, HasNoAngularMomentum, get_hr_b

############################ MassTrait #################################

abstract type MassTrait end
struct HasMass <: MassTrait end
struct HasNoMass <: MassTrait end

MassTrait(::S) where {S<:System} = error(
    "Please extend Components.MassTrait for $S")

"""
Notes:
- When get_mp_b is called on a System, the returned MassProperties
  instance must be expressed in the System's parent reference frame.
- At the root of the component hierarchy will typically be the airframe. The
  airframe is its own parent. Therefore, the MassProperties returned by the
  airframe must be expressed in its own reference frame: total airframe mass,
  position vector from the airframe origin to the airframe center of mass
  expressed in airframe axes, and inertia tensor of the airframe with respect to
  its origin, expressed in airframe axes. and these are the properties expected
  by the dynamics equations.
- Aircraft dynamics and kinematics are formulated on the airframe origin Ob
  instead of the aircraft's center of mass G. This allows for any of the
  aircraft's mass properties to change, either gradually (for example, due to
  fuel consumption) or suddenly (due to a payload release), without having to
  worry about discontinuities in the kinematic state vector.
"""

get_mp_b(sys::System) = get_mp_b(MassTrait(sys), sys)

get_mp_b(::HasMass, sys::System) = error(
    "$(typeof(sys)) has the HasMass trait, but no get_mp_b constructor is defined for it")

get_mp_b(::HasNoMass, sys::System) = MassProperties()

#default implementation for a SystemGroup with the HasMass trait, tries
#to compute the aggregate mass properties for all the subsystems
@inline @generated function get_mp_b(::HasMass, sys::System{D}) where {D<:SystemGroupDescriptor}

    # Core.print("Generated function called")
    ex = Expr(:block)
    push!(ex.args, :(p = MassProperties()))
    for label in fieldnames(D)
        push!(ex.args,
            :(p += get_mp_b(sys.subsystems[$(QuoteNode(label))])))
    end
    return ex

end
###################### AngularMomentumTrait ##########################

abstract type AngularMomentumTrait end
struct HasAngularMomentum <: AngularMomentumTrait end
struct HasNoAngularMomentum <: AngularMomentumTrait end

#prevents the trait system from failing silently when wrongly extended
AngularMomentumTrait(::S) where {S<:System} = error(
    "Please extend Components.AngularMomentumTrait for $S")

get_hr_b(sys::System) = get_hr_b(AngularMomentumTrait(sys), sys)

get_hr_b(::HasAngularMomentum, sys::System) = error(
    "$(typeof(sys)) has angular momentum, but no method get_hr_b was defined for it")

get_hr_b(::HasNoAngularMomentum, sys::System) = zeros(SVector{3})

#default implementation for a SystemGroup with the HasAngularMomentum trait, tries
#to sum the angular momentum from its individual components. override as required
@inline @generated function get_hr_b(::HasAngularMomentum, sys::System{D}) where {D<:SystemGroupDescriptor}

    # Core.print("Generated function called")
    ex = Expr(:block)
    push!(ex.args, :(h = SVector(0., 0., 0.))) #allocate
    for label in fieldnames(D)
        push!(ex.args,
            :(h += get_hr_b(sys.subsystems[$(QuoteNode(label))])))
    end
    return ex

end

###################### WrenchTrait ##########################

abstract type WrenchTrait end
struct HasWrench <: WrenchTrait end
struct HasNoWrench <: WrenchTrait end

#prevents the trait system from failing silently when wrongly extended
WrenchTrait(::S) where {S<:System} = error(
    "Please extend Components.WrenchTrait for $S")

get_wr_b(sys::System) = get_wr_b(WrenchTrait(sys), sys)

get_wr_b(::HasWrench, sys::System) = error(
    "$(typeof(sys)) is a Wrench source, but no method get_wr_b was defined for it")

get_wr_b(::HasNoWrench, sys::System) = Wrench()

#default implementation for a SystemGroup with the HasWrench trait, tries
#to sum all the Wrenches from its individual components. override as required
@inline @generated function get_wr_b(::HasWrench, sys::System{D}) where {D<:SystemGroupDescriptor}

    # Core.print("Generated function called")
    ex = Expr(:block)
    push!(ex.args, :(wr = Wrench())) #allocate a zero wrench
    for label in fieldnames(D)
        push!(ex.args,
            :(wr += get_wr_b(sys.subsystems[$(QuoteNode(label))])))
    end
    return ex

end


end #module