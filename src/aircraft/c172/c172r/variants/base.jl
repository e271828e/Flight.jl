module C172RBase

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using NLopt
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore
using Flight.FlightPhysics
using Flight.FlightComponents

using ..AircraftBase
using ..C172R

export Cessna172RBase


################################################################################
############################### Cessna172RBase #####################################

#Cessna172R with NoAvionics
const Cessna172RBase{K, T} = C172R.Aircraft{K, T, NoAvionics} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain}

function Cessna172RBase(kinematics = LTF(), terrain = HorizontalTerrain())
    C172R.Aircraft(kinematics, terrain, NoAvionics())
end

############################ Joystick Mappings #################################

#redirect input assignments directly to the actuation system
function IODevices.assign!(sys::System{<:Cessna172RBase}, joystick::Joystick,
                           mapping::InputMapping)
    IODevices.assign!(sys.physics.platform.act, joystick, mapping)
end

end #module