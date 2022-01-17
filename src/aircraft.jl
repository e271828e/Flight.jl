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
using Flight.Components

import Flight.Modeling: init_x, init_y, init_u, init_d,f_cont!, f_disc!
import Flight.Plotting: plots
import Flight.Kinematics: init!

######################### AircraftBase ########################

abstract type AbstractAircraftID end

struct GenericID <: AbstractAircraftID end

Base.@kwdef struct AircraftBase{I <: AbstractAircraftID,
                    K <: AbstractKinematics,
                    F <: SystemDescriptor,
                    C <: SystemDescriptor} <: SystemDescriptor

    identifier::I = GenericID()
    kinematics::K = KinLTF()
    airframe::F = RigidBody()
    controls::C = NullSystemDescriptor()
end

init_x(ac::AircraftBase) = ComponentVector(
    kin = init_x(ac.kinematics),
    afr = init_x(ac.airframe),
    ctl = init_x(ac.controls),
    )

init_u(ac::AircraftBase) = init_u(ac.controls)

init_y(ac::AircraftBase) = (
    kin = KinData(),
    dyn = DynData(),
    air = AirData(),
    afr = init_y(ac.airframe),
    ctl = init_y(ac.controls),
    )

init_d(ac::AircraftBase) = (
    afr = init_d(ac.airframe),
    ctl = init_d(ac.controls),
    )

const AircraftBaseSys{I,K,F,C} = System{AircraftBase{I,K,F,C}} where {I,K,F,C}

function System(ac::AircraftBase, ẋ = init_x(ac), x = init_x(ac),
                    y = init_y(ac), u = init_u(ac), d = init_d(ac), t = Ref(0.0))

    params = ()
    subsystems = (
        afr = System(ac.airframe, ẋ.afr, x.afr, y.afr, init_u(ac.airframe), d.afr, t),
        ctl = System(ac.controls, ẋ.ctl, x.ctl, y.ctl, u, d.ctl, t),)

    System{map(typeof, (ac, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)
end

init!(ac::System{<:AircraftBase}, init::KinInit) = init!(ac.x.kin, init)

function f_cont!(sys::AircraftBaseSys, trn::AbstractTerrain, atm::AtmosphericSystem)

    @unpack ẋ, x, subsystems = sys
    @unpack afr, ctl = subsystems

    #update kinematics
    kin = f_kin!(ẋ.kin.pos, x.kin)

    air = AirData(kin, atm)

    #update controls
    f_cont!(ctl, afr, kin, air, trn)
    #update airframe components
    f_cont!(afr, ctl, kin, air, trn)

    mass = get_mass_properties(afr)
    wr_ext_b = get_wr_b(afr)
    hr_b = get_hr_b(afr)

    # update dynamics
    dyn = f_dyn!(ẋ.kin.vel, kin, mass, wr_ext_b, hr_b)

    sys.y = (kin = kin, dyn = dyn, air = air, afr = afr.y, ctl = ctl.y)
    return nothing

end

function f_disc!(sys::AircraftBaseSys)
    @unpack afr, ctl = sys.subsystems

    x_mod = renormalize!(sys.x.kin, 1e-8) |
            f_disc!(afr, ctl) |
            f_disc!(ctl, afr)

    return x_mod
end

assign_joystick_inputs!(args...) = throw(MethodError(assign_joystick_inputs!, args))

#no custom plots function required, since all outputs are namedtuples

end #module