module C172R

using UnPack

using Flight.FlightCore
using Flight.FlightPhysics

using ..Aircraft

include("airframe.jl"); using .C172RAirframe
include("avionics.jl"); using .C172RAvionics

export Cessna172R

################################################################################
############################### Cessna172R #####################################

#Cessna172RTemplate requires a subtype of C172R.Airframe, allows installing any
#supported avionics and allows using different kinematic descriptions
const Cessna172RTemplate{K, F, V} = AircraftTemplate{K, F, V} where {K, F <: Airframe, V}

function Cessna172RTemplate(kinematics = LTF(), avionics = DirectControls())
    AircraftTemplate( kinematics, Airframe(), avionics)
end

#the default Cessna172R installs the C172R.DirectControls avionics, which
#provides a basic reversible direct control system
const Cessna172R{K, F} = Cessna172RTemplate{K, F, DirectControls} where {K, F <: Airframe}
Cessna172R(kinematics = LTF()) = Cessna172RTemplate(kinematics, DirectControls())

struct Cessna172RU{F, V}
    airframe::F
    avionics::V
end

Systems.init(::SystemU, ac::Cessna172R) = Cessna172RU(init_u(ac.airframe), init_u(ac.avionics))


################################################################################
############################ Joystick Mappings #################################

function IODevices.assign!(u::Cessna172RU, joystick::Joystick{XBoxControllerID},
                           mapping::InputMapping)
    IODevices.assign!(u.avionics, joystick, mapping)
end


################################################################################
############################### Analysis Tools #################################

include("tools/trim.jl"); using .Trim
include("tools/linear.jl"); using .Linear


end #module