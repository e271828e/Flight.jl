module C172MCSBase

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using NLopt
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore.Systems
using Flight.FlightCore.IODevices
using Flight.FlightCore.Joysticks

using Flight.FlightPhysics.Attitude
using Flight.FlightPhysics.Geodesy
using Flight.FlightPhysics.Kinematics
using Flight.FlightPhysics.Environment

using Flight.FlightComponents.Control
using Flight.FlightComponents.Piston
using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.World

using ..C172MCS

export Cessna172MCSBase

################################################################################
############################### Cessna172RBase #####################################

#Cessna172MCS with NoAvionics
const Cessna172MCSBase{K} = C172MCS.Template{K, NoAvionics} where {K}
Cessna172MCSBase(kinematics = LTF()) = C172MCS.Template(kinematics, NoAvionics())

############################ Joystick Mappings #################################

#redirect input assignments directly to the actuation system
function IODevices.assign!(sys::System{<:Cessna172MCSBase}, joystick::Joystick,
                           mapping::InputMapping)
    IODevices.assign!(sys.airframe.act, joystick, mapping)
end




end #module