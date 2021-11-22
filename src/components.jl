module Components

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.ModelingTools
import Flight.ModelingTools: System, init_x0, init_y0, init_u0, init_d0, f_cont!, f_disc!

using Flight.Plotting
import Flight.Plotting: plots

using Flight.Attitude

using Flight.Dynamics
import Flight.Dynamics: MassProperties

export MassTrait, HasMass, HasNoMass, get_mass_properties
export ExternalWrenchTrait, ExternalWrench, NoExternalWrench, get_wr_b
export AngularMomentumTrait, AngularMomentum, NoAngularMomentum, get_hr_b

############################ MassTrait #################################

abstract type MassTrait end
struct HasMass <: MassTrait end
struct HasNoMass <: MassTrait end

MassTrait(::Type{<:System}) = HasNoMass()

MassProperties(sys::S) where {S<:System} = MassProperties(MassTrait(S), sys)

MassProperties(::HasMass, sys::System) = error(
    "$(typeof(sys)) has the HasMass trait, but no MassProperties constructor is defined for it")

MassProperties(::HasNoMass, sys::System) = error(
    "$(typeof(sys)) has the HasNoMass trait, it cannot return a MassProperties instance")

#alternative implementation
# MassTrait(::System) = HasNoMass()
# MassProperties(sys::System) = MassProperties(MassTrait(sys), sys)

###################### AngularMomentumTrait ##########################

abstract type AngularMomentumTrait end
struct AngularMomentum <: AngularMomentumTrait end
struct NoAngularMomentum <: AngularMomentumTrait end

AngularMomentumTrait(::Type{<:System}) = NoAngularMomentum()

get_hr_b(sys::S) where {S<:System} = get_hr_b(AngularMomentumTrait(S), sys)

get_hr_b(::AngularMomentum, sys::System) = error(
    "$(typeof(sys)) has the AngularMomentum trait, but no method get_hr_b was defined for it")

get_hr_b(::NoAngularMomentum, ::System) = zeros(SVector{3})

#default implementation for a SystemGroup with the AngularMomentum trait, tries
#to sum the angular momentum from its individual components. override as required
@inline @generated function get_hr_b(sys::System{D}) where {D<:SystemGroupDescriptor}

    # Core.print("Generated function called")
    ex = Expr(:block)
    push!(ex.args, :(h = SVector(0., 0., 0.))) #allocate
    for label in fieldnames(D)
        push!(ex.args,
            :(h += get_hr_b(sys.subsystems[$(QuoteNode(label))])))
    end
    return ex

end

###################### ExternalWrenchTrait ##########################

abstract type ExternalWrenchTrait end
struct ExternalWrench <: ExternalWrenchTrait end
struct NoExternalWrench <: ExternalWrenchTrait end

ExternalWrenchTrait(::Type{<:System}) = NoExternalWrench()

get_wr_b(sys::S) where {S<:System} = get_wr_b(ExternalWrenchTrait(S), sys)

get_wr_b(::ExternalWrench, sys::System) = error(
    "$(typeof(sys)) has the ExternalWrench trait, but no method get_wr_b was defined for it")

get_wr_b(::NoExternalWrench, ::System) = Wrench()

#default implementation for a SystemGroup with the ExternalWrench trait, tries
#to sum all the Wrenches from its individual components. override as required
@inline @generated function get_wr_b(::ExternalWrench, sys::System{D}) where {D<:SystemGroupDescriptor}

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