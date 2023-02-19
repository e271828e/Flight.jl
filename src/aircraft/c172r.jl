module C172R

using UnPack

export Cessna172R

include("c172r/airframe.jl");
include("c172r/avionics.jl");


################################################################################
############################### Cessna172R #####################################

#Cessna172RBase requires a subtype of C172R.Airframe, but allows installing any
#avionics and using different kinematic descriptions
const Cessna172RBase{K, F, V} = AircraftTemplate{K, F, V} where {K, F <: Airframe, V}

function Cessna172RBase(kinematics = LTF(), avionics = MechanicalControls())
    AircraftTemplate( kinematics, Airframe(), avionics)
end

#the default Cessna172R installs the C172R.MechanicalControls avionics (which
#provides only a basic reversible control system)
const Cessna172R{K, F} = Cessna172RBase{K, F, MechanicalControls} where {K, F}
Cessna172R(kinematics = LTF()) = Cessna172RBase(kinematics, MechanicalControls())

#By default, the AircraftTemplate's U type will be automatically determined by
#the System's generic initializers as a NamedTuple{(:airframe, :avionics),
#Tuple{typeof(airframe)}, typeof(avionics)}, where typeof(airframe) and
#typeof(avionics) will typically be themselves NamedTuples. this makes
#dispatching on our specific Cessna172R <: AircraftTemplate impractical. a
#workaround is to define:

# const Cessna172RTemplate = Cessna172R()
# const Cessna172RX = typeof(init_x(Cessna172RTemplate))
# const Cessna172RU = typeof(init_u(Cessna172RTemplate))
# const Cessna172RY = typeof(init_y(Cessna172RTemplate))
# const Cessna172RS = typeof(init_s(Cessna172RTemplate))

#...which works, but has the drawback of slowing down module imports. An
#alternative is to override the default AircraftTemplate SystemU initializer
#with a custom, dispatch-friendly type, with the same fields:

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