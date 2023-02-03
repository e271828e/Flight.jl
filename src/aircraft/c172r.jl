module C172R

# using LinearAlgebra
using UnPack
# using Reexport
# using Interpolations
# using HDF5
# using Printf
# using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic

export Cessna172R

include("c172r/airframe.jl");
include("c172r/avionics.jl");

############################# Update Methods ###################################

#to be extended by any avionics installed on the aircraft
function map_controls!(airframe::System{<:Airframe}, avionics::System{<:AbstractAvionics})
    println("Please extend C172R.map_controls! for $(typeof(airframe)), $(typeof(avionics))")
    MethodError(map_controls!, (airframe, avionics))
end

function Systems.f_ode!(airframe::System{<:Airframe}, avionics::System{<:AbstractAvionics},
                kin::KinematicData, air::AirData, trn::System{<:AbstractTerrain})

    @unpack aero, pwp, ldg, fuel, pld = airframe

    map_controls!(airframe, avionics)
    f_ode!(aero, pwp, air, kin, trn)
    f_ode!(ldg, kin, trn) #update landing gear continuous state & outputs
    f_ode!(pwp, air, kin) #update powerplant continuous state & outputs
    f_ode!(fuel, pwp) #update fuel system

    update_y!(airframe)

end

function Systems.f_step!(airframe::System{<:Airframe}, avionics::System{<:AbstractAvionics},
                         ::KinematicSystem)
    @unpack aero, pwp, fuel, ldg, fuel, pld = airframe

    x_mod = false
    map_controls!(airframe, avionics)
    x_mod |= f_step!(aero)
    x_mod |= f_step!(ldg)
    x_mod |= f_step!(pwp, fuel)
    return x_mod

end

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