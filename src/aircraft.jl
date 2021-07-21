module Aircraft

# using StaticArrays: SVector, SMatrix
using LinearAlgebra
using UnPack

using Flight.LBV
# using Flight.WGS84
# using Flight.Attitude
using Flight.Kinematics


#1) gathers all wrenches and additional angular momentum from all aircraft
#components, adds them up
#2) from aircrat components, it computes MassData
#3) from XAircraft.kin computes PVData

# and then calls the dynamics method x_vel_dot


#maybe Aircraft has a mass object as a field whose purpose is to output mass
#data

#see Python for the contents of an Aircraft. But it does not need to store
#"Airframe". by the very type of the XKin defined by the aircraft, dispatch
#takes care of the appropriate method of x_vel_dot that will be called

abstract type MassModel end
abstract type LdgGroup end
abstract type PwpGroup end
abstract type SrfGroup end

struct ConstantMassModel <: MassModel
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end

##### TestAircraft #####

#see if we can use type parameters to define the specific mass, ldg, pwp and srf
#groups we are using in a given aircraft. this allows us to have different ldg,
#pwp, srf models, etc, without making the aircraft abstract, and therefore reuse
#the same x_dot function, which assumes the existence of these fields. this
#requires LdgGroup, etc. to be isbits types. they can only hold numbers,
#StaticVectors, etc. no dicts, no lists, no generic arrays, etc. it seems doable!

@define_node XTestAircraft (kin = XKinWGS84, )

struct TestAircraft
    mass::ConstantMassModel
    # ldg_group::LdgGroup
    # pwp_group::PwpGroup
    # srf_group::SrfGroup
end

#TerrainModel does not belong to the Aircraft itself. it must be defined
#separately, and passed as an argument
#the same goes for Atmosphere
#there must be a level above the Aircraft, which will typically be the
#simulation, that defines the terrain and atmospheric models, and holds all the
#aircraft participating in the simulation. this may be a block based simulation
#or a custom made one. but it must exist in some form

#MassModel is an abstract type that, given some inputs (typically state of the
#fuel system and external payloads), returns a MassData struct. but for now, we
#can simply define a ConstantMassModel

end