module Aircraft

# using StaticArrays: SVector, SMatrix
using LinearAlgebra
using UnPack

using Flight.LBV
# using Flight.WGS84
# using Flight.Attitude
using Flight.Kinematics

export XTestAircraft
export dt

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
#StaticVectors, etc. no dicts, no lists, no generic arrays, etc. it seems
#doable! however, come to think of it, there is no actual need to store the
#different models as type parameters! well, because each

#the principle is that, regardless of any model's complexity, anything that can
#change during operation will be stored separately in its corresponding state
#and input vectors. and these will not be stored in the structs themselves. not
#only that,

#on the other hand, it seems desirable to store references to each state and
#input vector sub block as part of its owner subsystem struct.
#so what do we do?

#option 1)
#lets assume we want to compute x_ldg_dot. in our aircraft we have a TricycleLdg
#<: LdgGroup. for this TricycleLdg we will have a dt!(x_dot::XTricycleLdg,
#x::XTricycleLdg, u::UTricycleLdg, ::PVDataWGS84, ::TricycleLdg, ::TerrainModel) which updates
#x_dot::XTricycleLdg and also returns a TricycleLdgData. the reason to pass
#x_dot is that it will be a subblock view of a previously allocated overall
#aircraft state vector. this avoids having to allocate each block and sub-block
#down the hierarchy to finally assemble them in an upstream process. this way,
#the flow is only top to bottom.

#option 2) the alternative is the Python way: when the aircraft is created, we
#store the whole state vector LBV in a struct field x::XTestAircraft. we also
#define x::XTricycleLdg within aircraft.ldg::TricycleLdg. we define an
#assign_state_vector method that recursively assigns state vector sub blocks
#down the hierarchy. for example, when we call
#assign_state_vector(x::XTestAircraft, a::TestAircraft)), this method calls
#assign_state_vector(x.ldg, a.ldg), which is defined by TricycleLdg with the
#appropriate signature. this method in turn calls assign_state_vector(x.nlg,
#l.nlg), etc. we do the same for xdot, x and u. once this process is finished,
#we no longer need to pass the xdot, x and u vectors along the dot methods,
#because their corresponding sub blocks will be available as views in each of
#the subsystems down the hierarchy. this cleans up the dt method signature. now
#it would look like dt!(sys::TricycleLdg, ::PVDataWGS84, ::TerrainModel) but
#in exchange requires defining those recursive initialization methods for each
#subsystem. also, systems are no longer of isbits types. the advantage of having
#isbits types was that we could use them as type parameters, which in turn
#allows Aircraft not to be an abstract type. does this really matter that much?
#the only disadvantage of having an abstract type is the difficulty in reusing
#the main dt(::Aircraft) function, which is similar for every aircraft and has

#option 1a) this is a variation of 1) in which we define a mutable struct
#TricycleLdg with fields: x::XTricycleLdg{Float64}, u::UTricycleLdg{Float64},
#p::TricycleLdgParams, just to simplify the method argument list for dt!, which
#would look like dt!(x_dot::XTricycleLdg, sys::TricycleLdg, ::PVDataWGS84,
#::TerrainModel). middle of nowhere.

#Let's try first option 1. See if it does not get too ugly.

#now, it could be argued that if we stored references to

@define_node XTestAircraft (kin = XKinWGS84, )

#do we put the aircraft state and input vectors INSIDE the struct? Note that
#this will cause it to lose its itsbitstype condition

struct TestAircraft
    mass::ConstantMassModel
    # ldg_group::LdgGroup
    # pwp_group::PwpGroup
    # srf_group::SrfGroup
end

function dt(x::XTestAircraft, t::Real)
    x_kin = x.kin

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