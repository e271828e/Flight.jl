module Aircraft

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.Modeling
using Flight.Plotting
using Flight.Terrain
using Flight.Atmosphere
using Flight.Airdata
using Flight.Kinematics
using Flight.Dynamics
# using Flight.Components

import Flight.Modeling: System, init_x, init_y, init_u, init_d,f_cont!, f_disc!
import Flight.Dynamics: MassTrait, get_mp_b
import Flight.Plotting: plots
import Flight.Kinematics: init!

export AircraftBase

######################### AircraftBase ########################

abstract type AbstractAircraftID end

struct GenericID <: AbstractAircraftID end

######################## EmptyAirframe ########################

Base.@kwdef struct EmptyAirframe <: SystemDescriptor
    mass_distribution::RigidBody = RigidBody(1, SA[1.0 0 0; 0 1.0 0; 0 0 1.0])
end

get_mp_b(sys::System{EmptyAirframe}) = MassProperties(sys.params.mass_distribution)




Base.@kwdef struct AircraftBase{I <: AbstractAircraftID,
                    K <: AbstractKinematics,
                    F <: SystemDescriptor,
                    C <: SystemDescriptor} <: SystemGroupDescriptor

    # identifier::I = GenericID()
    kinematics::K = KinLTF()
    airframe::F = EmptyAirframe()
    controls::C = NullSystemDescriptor()
end

init_x(::Type{T}) where {T<:AircraftBase{I,K,F,C}} where {I,K,F,C} =
    ComponentVector(
    kin = init_x(K),
    afm = init_x(F),
    ctl = init_x(C),
    )

init_u(::Type{T}) where {T<:AircraftBase{I,K,F,C}} where {I,K,F,C} = (
    afm = init_u(F),
    ctl = init_u(C),
    )

init_y(::Type{T}) where {T<:AircraftBase{I,K,F,C}} where {I,K,F,C} = (
    kin = KinData(),
    dyn = DynData(),
    air = AirData(),
    afm = init_y(F),
    ctl = init_y(C),
    )

init_d(::Type{T}) where {T<:AircraftBase{I,K,F,C}} where {I,K,F,C} = (
    afm = init_d(F),
    ctl = init_d(C),
    )

const AircraftBaseSys{I,K,F,C} = System{AircraftBase{I,K,F,C}} where {I,K,F,C}

function System(ac::T, ẋ = init_x(T), x = init_x(T), y = init_y(T),
                u = init_u(T), d = init_d(T), t = Ref(0.0)) where {T<:AircraftBase}

    params = ()
    subsystems = (
        afm = System(ac.airframe, ẋ.afm, x.afm, y.afm, u.afm, d.afm, t),
        ctl = System(ac.controls, ẋ.ctl, x.ctl, y.ctl, u.ctl, d.ctl, t),)

    println(x)
    System{map(typeof, (ac, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)
end

init!(ac::System{<:AircraftBase}, init::KinInit) = init!(ac.x.kin, init)

function f_cont!(sys::AircraftBaseSys, trn::AbstractTerrain, atm::AtmosphericSystem)

    @unpack ẋ, x, subsystems = sys
    @unpack afm, ctl = subsystems

    #update kinematics
    kin = f_kin!(ẋ.kin.pos, x.kin)

    air = AirData(kin, atm)

    #update controls
    f_cont!(ctl, afm, kin, air, trn)
    #update airframe components
    f_cont!(afm, ctl, kin, air, trn)

    mp_b = get_mp_b(afm)
    wr_b = get_wr_b(afm)
    hr_b = get_hr_b(afm)

    # update dynamics
    dyn = f_dyn!(ẋ.kin.vel, kin, mp_b, wr_b, hr_b)

    sys.y = (kin = kin, dyn = dyn, air = air, afm = afm.y, ctl = ctl.y)
    return nothing

end

function f_disc!(sys::AircraftBaseSys)
    @unpack afm, ctl = sys.subsystems

    x_mod = renormalize!(sys.x.kin, 1e-8) |
            f_disc!(afm, ctl) |
            f_disc!(ctl, afm)

    return x_mod
end

assign_joystick_inputs!(args...) = throw(MethodError(assign_joystick_inputs!, args))

#no custom plots function required, since all outputs are namedtuples

end #module