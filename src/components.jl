module Components

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.Modeling
using Flight.Plotting
import Flight.Plotting: plots

using Flight.Attitude

using Flight.Dynamics

export HasMass, HasNoMass, get_mass_properties
export HasWrench, HasNoWrench, get_wr_b
export HasAngularMomentum, HasNoAngularMomentum, get_hr_b

############################ MassTrait #################################

abstract type MassTrait end
struct HasMass <: MassTrait end
struct HasNoMass <: MassTrait end

MassTrait(::S) where {S<:System} = error(
    "Please extend Components.MassTrait for $S")

get_mass_properties(sys::System) = get_mass_properties(MassTrait(sys), sys)

get_mass_properties(::HasMass, sys::System) = error(
    "$(typeof(sys)) has the HasMass trait, but no get_mass_properties method is defined for it")

get_mass_properties(::HasNoMass, sys::System) = error(
    "$(typeof(sys)) has the HasNoMass trait")

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