module C172R

using UnPack

using Flight.Engine.Systems
using Flight.Engine.IODevices
using Flight.Engine.Joysticks

using Flight.Physics.Attitude
using Flight.Physics.Kinematics

using Flight.Components.Aircraft

include("c172r/airframe.jl"); using .C172RAirframe
include("c172r/avionics.jl"); using .C172RAvionics

export Cessna172R

################################################################################
############################### Cessna172R #####################################

#Cessna172RBase requires a subtype of C172R.Airframe, but allows installing any
#avionics and using different kinematic descriptions
const Cessna172RBase{K, F, V} = AircraftTemplate{K, F, V} where {K, F <: Airframe, V}

function Cessna172RBase(kinematics = LTF(), avionics = ReversibleControls())
    AircraftTemplate( kinematics, Airframe(), avionics)
end

#the default Cessna172R installing the C172R.ReversibleControls avionics, which
#provides a basic reversible direct control system
const Cessna172R{K, F} = Cessna172RBase{K, F, ReversibleControls} where {K, F}
Cessna172R(kinematics = LTF()) = Cessna172RBase(kinematics, ReversibleControls())

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

include("c172r/tools/trim.jl"); using .Trim
include("c172r/tools/linear.jl"); using .Linear




end #module