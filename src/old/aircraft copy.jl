module Aircraft

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.Modeling

using Flight.Terrain
using Flight.Atmosphere

using Flight.StateMachine
using Flight.Airframe
using Flight.Mass
using Flight.Airdata
# using Flight.Propulsion
# using Flight.Aerodynamics
# using Flight.LandingGear

using Flight.Kinematics
using Flight.Dynamics

import Flight.Modeling: System, init_x0, init_y0, init_u0, init_d0,f_cont!, f_disc!
import Flight.Dynamics: MassData

using Flight.Plotting
import Flight.Plotting: plots

export AircraftBase, AircraftBaseD, AircraftBaseY
export AbstractMass, ConstantMass, TunableMass
export AbstractAircraftID

######################### AircraftBase ########################

abstract type AbstractAircraftID end

struct AircraftBaseID <: AbstractAircraftID end

Base.@kwdef struct AircraftBase{
                    Id <: AbstractAircraftID,
                    Kin <: AbstractKinematics,
                    #SMac <: AbstractStateMachine,
                    Mass <: AbstractMass,
                    Aero <: AbstractAirframeComponent,
                    Pwp <: AbstractAirframeComponent,
                    Ldg <: AbstractAirframeComponent,
                    Srf <: AbstractAirframeComponent} <: SystemDescriptor
    id::Id = AircraftBaseID()
    kin::Kin = KinLTF()
    mass::Mass = TunableMass()
    aero::Aero = NullAirframeComponent()
    pwp::Pwp = NullAirframeComponent()
    ldg::Ldg = NullAirframeComponent()
    srf::Srf = NullAirframeComponent()
end

struct AircraftBaseY{MassY, AeroY, PwpY, LdgY, SrfY}
    dyn::DynData
    kin::KinData
    air::AirData
    mass::MassY
    aero::AeroY
    pwp::PwpY
    ldg::LdgY
    srf::SrfY
end

struct AircraftBaseD{MassD, AeroD, PwpD, LdgD, SrfD}
    mass::MassD
    aero::AeroD
    pwp::PwpD
    ldg::LdgD
    srf::SrfD
end

#unlike X, D and Y, which must hold all their aircraft subsystems' counterparts,
#U is arbitrarily defined as a design choice

struct NoAircraftU end

init_u0(::AircraftBase) = NoAircraftU()

init_x0(ac::AircraftBase) = ComponentVector(
    kin = init_x0(ac.kin),
    mass = init_x0(ac.mass),
    aero = init_x0(ac.aero),
    pwp = init_x0(ac.pwp),
    ldg = init_x0(ac.ldg),
    srf = init_x0(ac.srf)
    )

init_y0(ac::AircraftBase) = AircraftBaseY(
    DynData(),
    KinData(),
    AirData(),
    init_y0(ac.mass),
    init_y0(ac.aero),
    init_y0(ac.pwp),
    init_y0(ac.ldg),
    init_y0(ac.srf)
    )

init_d0(ac::AircraftBase) = AircraftBaseD(
    init_d0(ac.mass),
    init_d0(ac.aero),
    init_d0(ac.pwp),
    init_d0(ac.ldg),
    init_d0(ac.srf)
    )


const AircraftBaseSys{I,K,M,A,P,L,S} = System{AircraftBase{I,K,M,A,P,L,S}} where {I,K,M,A,P,L,S}


function System(ac::AircraftBase, ẋ = init_x0(ac), x = init_x0(ac),
                    y = init_y0(ac), u = init_u0(ac), d = init_d0(ac), t = Ref(0.0))

    #each subsystem allocates its own u, then we can decide how the aircraft's u
    #should map onto it via assign_control_inputs!
    # smac = System(ac.smac, ẋ.smac, x.smac, y.smac, init_u0(ac.smac), d.smac, t)
    mass = System(ac.mass, ẋ.mass, x.mass, y.mass, init_u0(ac.mass), d.mass, t)
    aero = System(ac.aero, ẋ.aero, x.aero, y.aero, init_u0(ac.aero), d.aero, t)
    pwp = System(ac.pwp, ẋ.pwp, x.pwp, y.pwp, init_u0(ac.pwp), d.pwp, t)
    ldg = System(ac.ldg, ẋ.ldg, x.ldg, y.ldg, init_u0(ac.ldg), d.ldg, t)
    srf = System(ac.srf, ẋ.srf, x.srf, y.srf, init_u0(ac.srf), d.srf, t)

    params = ()

    subsystems = (mass = mass, aero = aero, pwp = pwp, ldg = ldg, srf = srf)

    System{map(typeof, (ac, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)
end


#for an aircraft with no specific mapping, the user is expected to act upon
#the subsystems inputs directly, so this function has nothing to do
assign_control_inputs!(::AircraftBaseSys) = nothing

#overriding this function for specific AircraftBase type parameter
#combinations allows us to customize how the aircraft's control inputs map
#into its subsystems's control inputs.

#=
# example:
const my_pwp = AirframeGroup(left = EThruster(), right = EThruster())
const MyPwp = typeof(my_pwp)
init_u0(::AircraftBase{MyAircraftID}) = ComponentVector(throttle = 0.0)
function assign_control_inputs!(::AircraftBase{MyAircraftID})
    ac.subsystems.pwp.left.throttle = ac.u.throttle
    ac.subsystems.pwp.right.throttle = ac.u.throttle
end
=#


function f_cont!(ac_sys::AircraftBaseSys, trn::AbstractTerrain, atm::AtmosphericSystem)

    @unpack ẋ, x, u, params, subsystems = ac_sys
    @unpack mass, aero, pwp, ldg, srf = subsystems

    kin_data = f_kin!(ẋ.kin.pos, x.kin)

    air_data = AirData(kin_data, atm)

    #before updating the subsystems, assign their control inputs from the
    #aircraft's control input vector
    assign_control_inputs!(ac_sys)

    #f_cont!(smac)
    f_cont!(mass)
    f_cont!(srf, air_data) #update surface actuator continuous state & outputs
    f_cont!(ldg, kin_data, trn) #update landing gear continuous state & outputs
    f_cont!(pwp, kin_data, air_data) #update powerplant continuous state & outputs
    f_cont!(aero, air_data, kin_data, srf, trn) #requires previous srf update

    mass_data = MassData(mass)

    #initialize external Wrench and additional angular momentum
    wr_ext_b = get_wr_b(pwp) + get_wr_b(ldg) + get_wr_b(aero)
    # wr_ext_b = get_wr_b(pwp) + get_wr_b(ldg)
    hr_b = get_hr_b(pwp) + get_hr_b(ldg) + get_hr_b(aero)

    # update dynamics
    dyn_data = f_dyn!(ẋ.kin.vel, kin_data, mass_data, wr_ext_b, hr_b)

    ac_sys.y = AircraftBaseY(dyn_data, kin_data, air_data, mass.y, aero.y, pwp.y, ldg.y, srf.y)

    return nothing
end


function f_disc!(ac_sys::AircraftBaseSys)
    @unpack mass, aero, pwp, ldg, srf = ac_sys.subsystems

            # f_disc!(smac) |
    x_mod = renormalize!(ac_sys.x.kin, 1e-8) |
            f_disc!(mass) |
            f_disc!(aero) |
            f_disc!(pwp) |
            f_disc!(ldg) |
            f_disc!(srf)

    return x_mod
end

function plots(t::AbstractVector{<:Real}, data::AbstractVector{<:AircraftBaseY};
    mode::Symbol, save_path::Union{String, Nothing}, kwargs...)

    sa = StructArray(data)

    ss_path = (dyn = "dynamics", kin = "kinematics", air = "airdata", mass = "mass",
        aero = "aerodynamics", pwp = "powerplant", ldg = "landinggear", srf = "surfaces")

    for (label, path) ∈ zip(keys(ss_path), values(ss_path))
        println("Generating plots for $label")
        plots(t, getproperty(sa, label);
                mode, save_path = mkpath(joinpath(save_path, path)), kwargs...)
    end

end



end