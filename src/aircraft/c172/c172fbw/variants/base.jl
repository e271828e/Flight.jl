module C172FBWBase

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

using ...C172
using ..C172FBW

export Cessna172FBWBase


################################################################################
############################### Cessna172RBase #####################################

#Cessna172FBW with NoAvionics
const Cessna172FBWBase{K} = C172FBW.Template{K, NoAvionics} where {K}
Cessna172FBWBase(kinematics = LTF()) = C172FBW.Template(kinematics, NoAvionics())

################################ Tools #########################################

function Aircraft.trim!(ac::System{<:Cessna172FBWBase}, args...; kwargs...)
    trim!(ac.physics, args...; kwargs...)
end

############################ Joystick Mappings #################################

#redirect input assignments directly to the actuation system
function IODevices.assign!(sys::System{<:Cessna172FBWBase}, joystick::Joystick,
                           mapping::InputMapping)
    IODevices.assign!(sys.airframe.act, joystick, mapping)
end




end #module