module C172FBWBase

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using NLopt
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore
using Flight.FlightPhysics
using Flight.FlightComponents

using ..AircraftBase
using ..C172FBW

export Cessna172FBWBase


################################################################################
############################### Cessna172FBWBase ###############################

#Cessna172FBW with NoAvionics
const Cessna172FBWBase{K, T} = C172FBW.Aircraft{K, T, NoAvionics} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain}

function Cessna172FBWBase(kinematics = LTF(), terrain = HorizontalTerrain())
    C172FBW.Aircraft(kinematics, terrain, NoAvionics())
end

############################ Joystick Mappings #################################

#redirect input assignments directly to the actuation system
function IODevices.assign!(sys::System{<:Cessna172FBWBase}, joystick::Joystick,
                           mapping::InputMapping)
    IODevices.assign!(sys.physics.platform.act, joystick, mapping)
end

end #module