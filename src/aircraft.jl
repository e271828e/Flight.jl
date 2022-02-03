module Aircraft

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.Modeling
using Flight.Plotting
using Flight.Terrain
using Flight.Attitude
using Flight.Geodesy
using Flight.Atmosphere
using Flight.Airdata
using Flight.Kinematics
using Flight.Dynamics
using Flight.Input
using Flight.Output

import Flight.Modeling: init_x, init_y, init_u, init_d,f_cont!, f_disc!
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_mp_b
import Flight.Plotting: plots
import Flight.Output: update!

export AircraftBase, AbstractAirframe, EmptyAirframe


abstract type AbstractAircraftID end

###############################################################################
############################## Airframe #######################################

abstract type AbstractAirframe <: SystemGroupDescriptor end
MassTrait(::System{<:AbstractAirframe}) = HasMass()

########################## EmptyAirframe ###########################

Base.@kwdef struct EmptyAirframe <: AbstractAirframe
    mass_distribution::RigidBody = RigidBody(1, SA[1.0 0 0; 0 1.0 0; 0 0 1.0])
end

WrenchTrait(::System{EmptyAirframe}) = HasNoWrench()
AngularMomentumTrait(::System{EmptyAirframe}) = HasNoAngularMomentum()

get_mp_b(sys::System{EmptyAirframe}) = MassProperties(sys.params.mass_distribution)

@inline f_cont!(::System{EmptyAirframe}, args...) = nothing
@inline (f_disc!(::System{EmptyAirframe}, args...)::Bool) = false


###############################################################################
############################## AircraftBase ###################################

struct GenericID <: AbstractAircraftID end

struct AircraftBase{I <: AbstractAircraftID,
                    K <: AbstractKinematics,
                    F <: AbstractAirframe,
                    V <: SystemDescriptor} <: SystemGroupDescriptor

    kinematics::K
    airframe::F
    avionics::V
end

function AircraftBase(     ::I = GenericID();
                        kin::K = KinLTF(),
                        afm::F = EmptyAirframe(),
                        avs::V = NullSystemDescriptor()) where {I,K,F,V}
    AircraftBase{I,K,F,V}(kin, afm, avs)
end

#override the default SystemGroupDescriptor implementation, because we need to
#add some stuff besides subsystem outputs
init_y(::Type{T}) where {T<:AircraftBase{I,K,F,V}} where {I,K,F,V} = (
    kinematics = init_y(K),
    airframe = init_y(F),
    avionics = init_y(V),
    dynamics = DynData(),
    air = AirData(),
    )

function init!(ac::System{T}, init::KinInit) where {T<:AircraftBase{I,K}} where {I,K}
    ac.x.kinematics .= init_x(K, init)
end

function f_cont!(sys::System{<:AircraftBase}, trn::AbstractTerrain, atm::AtmosphericSystem)

    @unpack ẋ, x, subsystems = sys
    @unpack kinematics, airframe, avionics = subsystems

    #update kinematics
    f_cont!(kinematics)
    kin_data = kinematics.y
    air_data = AirData(kin_data, atm)

    #update avionics and airframe components
    f_cont!(avionics, airframe, kin_data, air_data, trn)
    f_cont!(airframe, avionics, kin_data, air_data, trn)

    mp_b = get_mp_b(airframe)
    wr_b = get_wr_b(airframe)
    hr_b = get_hr_b(airframe)

    #update velocity derivatives
    dyn_data = f_dyn!(kinematics.ẋ.vel, kinematics.y, mp_b, wr_b, hr_b)

    sys.y = (kinematics = kinematics.y, airframe = airframe.y, avionics = avionics.y,
            dynamics = dyn_data, air = air_data,)

    return nothing

end

function f_disc!(sys::System{<:AircraftBase})
    @unpack kinematics, airframe, avionics = sys.subsystems

    x_mod = f_disc!(kinematics, 1e-8) |
            f_disc!(airframe, avionics) |
            f_disc!(avionics, airframe)

    return x_mod
end

function update!(xp::XPInterface, ac::System{<:AircraftBase}, ac_number::Integer = 0)
    update!(xp, ac.y.kinematics.pos, ac_number)
end

function update!(xp::XPInterface, pos::PosData, aircraft::Integer = 0)

    llh = Geographic(pos.ϕ_λ, pos.h_o)
    euler = REuler(pos.q_nb)

    lat = rad2deg(llh.l2d.ϕ)
    lon = rad2deg(llh.l2d.λ)
    alt = llh.alt

    psi = rad2deg(euler.ψ)
    theta = rad2deg(euler.θ)
    phi = rad2deg(euler.φ)

    Output.set_position!(xp; lat, lon, alt, psi, theta, phi, aircraft)

end

#no custom plots function required, all outputs are namedtuples

end #module