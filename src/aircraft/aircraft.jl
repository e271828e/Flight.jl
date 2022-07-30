module Aircraft

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.Systems
using Flight.Attitude
using Flight.Geodesy, Flight.Terrain, Flight.Atmosphere
using Flight.Kinematics, Flight.RigidBody
using Flight.Input, Flight.Output

import Flight.Systems: init, f_cont!, f_disc!
import Flight.RigidBody: MassTrait, WrenchTrait, AngularMomentumTrait, get_mp_b
import Flight.Output: update!

export AircraftBase, AbstractAirframe, AbstractAerodynamics, AbstractAvionics


###############################################################################
############################## Airframe #######################################

abstract type AbstractAirframe <: SystemDescriptor end
MassTrait(::System{<:AbstractAirframe}) = HasMass()
WrenchTrait(::System{<:AbstractAirframe}) = GetsExternalWrench()
AngularMomentumTrait(::System{<:AbstractAirframe}) = HasAngularMomentum()

########################## EmptyAirframe ###########################

Base.@kwdef struct EmptyAirframe <: AbstractAirframe
    mass_distribution::RigidBodyDistribution = RigidBodyDistribution(1, SA[1.0 0 0; 0 1.0 0; 0 0 1.0])
end

WrenchTrait(::System{EmptyAirframe}) = GetsNoExternalWrench()
AngularMomentumTrait(::System{EmptyAirframe}) = HasNoAngularMomentum()

get_mp_b(sys::System{EmptyAirframe}) = MassProperties(sys.params.mass_distribution)

@inline f_cont!(::System{EmptyAirframe}, args...) = nothing
@inline (f_disc!(::System{EmptyAirframe}, args...)::Bool) = false

####################### AbstractAerodynamics ##########################

abstract type AbstractAerodynamics <: SystemDescriptor end

MassTrait(::System{<:AbstractAerodynamics}) = HasNoMass()
WrenchTrait(::System{<:AbstractAerodynamics}) = GetsExternalWrench()
AngularMomentumTrait(::System{<:AbstractAerodynamics}) = HasNoAngularMomentum()


###############################################################################
######################### AbstractAvionics ####################################

abstract type AbstractAvionics <: SystemDescriptor end

struct NoAvionics <: AbstractAvionics end

@inline f_cont!(::System{NoAvionics}, args...) = nothing
@inline (f_disc!(::System{NoAvionics}, args...)::Bool) = false

###############################################################################
############################## AircraftBase ###################################

struct AircraftBase{K <: AbstractKinematics,
                    F <: AbstractAirframe,
                    A <: AbstractAvionics} <: SystemDescriptor
    kinematics::K
    airframe::F
    avionics::A
end

function AircraftBase(kinematics::K = LTF(),
                      airframe::F = EmptyAirframe(),
                      avionics::A = NoAvionics()) where {I,K,F,A}
    AircraftBase{K,F,A}(kinematics, airframe, avionics)
end

#override the default SystemDescriptor implementation, because we need to
#add some stuff besides subsystem outputs
init(::SystemY, ac::AircraftBase) = (
    airframe = init_y(ac.airframe),
    avionics = init_y(ac.avionics),
    kinematics = init_y(ac.kinematics),
    rigidbody = RigidBodyData(),
    airflow = AirflowData(),
    )

function init!(ac::System{<:AircraftBase}, ic::KinematicInit)
    Kinematics.init!(ac.x.kinematics, ic)
end

function f_cont!(sys::System{<:AircraftBase}, atm::System{<:SimpleAtmosphere}, trn::AbstractTerrain)

    @unpack ẋ, x, subsystems = sys
    @unpack kinematics, airframe, avionics = subsystems

    #update kinematics
    f_cont!(kinematics)
    kin_data = kinematics.y.common
    air_data = AirflowData(kin_data, atm)

    #update avionics and airframe components
    f_cont!(avionics, airframe, kin_data, air_data, trn)
    f_cont!(airframe, avionics, kin_data, air_data, trn)

    mp_b = get_mp_b(airframe)
    wr_b = get_wr_b(airframe)
    hr_b = get_hr_b(airframe)

    #update velocity derivatives
    rb_data = f_dyn!(kinematics.ẋ.vel, kin_data, mp_b, wr_b, hr_b)

    sys.y = (airframe = airframe.y, avionics = avionics.y, kinematics = kinematics.y,
            rigidbody = rb_data, airflow = air_data,)

    return nothing

end

function f_disc!(sys::System{<:AircraftBase})
    @unpack kinematics, airframe, avionics = sys

    x_mod = false
    x_mod = x_mod || f_disc!(kinematics, 1e-8)
    x_mod = x_mod || f_disc!(airframe, avionics, kinematics)
    x_mod = x_mod || f_disc!(avionics, airframe, kinematics)

    return x_mod
end

function update!(xp::XPInterface, kin::KinematicData, aircraft::Integer = 0)

    ll = LatLon(kin.n_e)
    e_nb = REuler(kin.q_nb)

    lat = rad2deg(ll.ϕ)
    lon = rad2deg(ll.λ)
    h = kin.h_o

    psi = rad2deg(e_nb.ψ)
    theta = rad2deg(e_nb.θ)
    phi = rad2deg(e_nb.φ)

    Output.set_position!(xp; lat, lon, h, psi, theta, phi, aircraft)

end


end #module