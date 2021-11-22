module Aircraft

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.ModelingTools
import Flight.ModelingTools: System, init_x0, init_y0, init_u0, init_d0,f_cont!, f_disc!

using Flight.Terrain
using Flight.Atmosphere
using Flight.Airdata

using Flight.Kinematics
using Flight.Dynamics
import Flight.Dynamics: MassProperties

using Flight.Components

using Flight.Plotting
import Flight.Plotting: plots

export AbstractAircraftID
export AircraftBase, AircraftAirframe

######################### AircraftBase ########################

abstract type AbstractAircraftID end

struct AircraftBaseID <: AbstractAircraftID end


struct AircraftBase{I <: AbstractAircraftID,
                    K <: AbstractKinematics,
                    C <: SystemDescriptor,
                    T <: SystemDescriptor} <: SystemDescriptor

    identifier::I
    kin::K
    cmp::C
    ctl::T
end

init_x0(ac::AircraftBase) = ComponentVector(
    kin = init_x0(ac.kin),
    cmp = init_x0(ac.cmp),
    ctl = init_x0(ac.ctl),
    )

init_u0(ac::AircraftBase) = init_u0(ac.ctl)

init_y0(ac::AircraftBase) = (
    kin = KinData(),
    dyn = DynData(),
    air = AirData(),
    cmp = init_y0(ac.cmp),
    ctl = init_y0(ac.ctl),
    )

init_d0(ac::AircraftBase) = (
    cmp = init_d0(ac.cmp),
    ctl = init_d0(ac.ctl),
    )

const AircraftBaseSys{I,K,C,T} = System{AircraftBase{I,K,C,T}} where {I,K,C,T}

function System(ac::AircraftBase, ẋ = init_x0(ac), x = init_x0(ac),
                    y = init_y0(ac), u = init_u0(ac), d = init_d0(ac), t = Ref(0.0))

    params = ()
    subsystems = (
        cmp = System(ac.cmp, ẋ.cmp, x.cmp, y.cmp, init_u0(ac.cmp), d.cmp, t),
        ctl = System(ac.ctl, ẋ.ctl, x.ctl, y.ctl, u, d.ctl, t),)

    System{map(typeof, (ac, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)
end


function f_cont!(sys::AircraftBaseSys, trn::AbstractTerrain, atm::AtmosphericSystem)

    @unpack ẋ, x, subsystems = sys
    @unpack cmp, ctl = subsystems

    #update kinematics
    kin = f_kin!(ẋ.kin.pos, x.kin)

    air = AirData(kin, atm)

    # #update controls
    f_cont!(ctl, cmp, kin, air, trn)
    #update components
    f_cont!(cmp, ctl, kin, air, trn)

    mass = MassProperties(cmp)
    wr_ext_b = get_wr_b(cmp)
    hr_b = get_hr_b(cmp)

    # update dynamics
    dyn = f_dyn!(ẋ.kin.vel, kin, mass, wr_ext_b, hr_b)

    sys.y = (kin = kin, dyn = dyn, air = air, cmp = cmp.y, ctl = ctl.y)
    return nothing

end

function f_disc!(sys::AircraftBaseSys)
    @unpack cmp, ctl = sys.subsystems

    x_mod = renormalize!(sys.x.kin, 1e-8) |
            f_disc!(cmp, ctl) |
            f_disc!(ctl, cmp)

    return x_mod
end

#no custom plots function required, since all outputs are namedtuples

end #module