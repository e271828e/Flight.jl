module Aircraft

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.ModelingTools

using Flight.Terrain
using Flight.Atmosphere

using Flight.StateMachine
using Flight.Airframe
using Flight.Airdata
using Flight.Propulsion
using Flight.Aerodynamics
# using Flight.LandingGear

using Flight.Kinematics
using Flight.Dynamics

import Flight.ModelingTools: System, get_x0, get_y0, get_u0, get_d0,f_cont!, f_disc!
import Flight.Dynamics: MassData

using Flight.Plotting
import Flight.Plotting: plots

export TestAircraft, TestAircraftD, TestAircraftY
export NoControlMapping, ConstantMass

abstract type AbstractMass <: AbstractComponent end


#given the fuel system outputs and the installed payloads, a MassSystem must
#return a MassData instance. open questions:

#1) should this information be passed as inputs along with the MassSystem
#   instance to the MassData constructor? or should they be defined as part of
#   the MassSystem's inputs and assigned explicitly?
function MassData(::T) where {T<:System{<:AbstractMass}}
    error("MassData constructor not implemented for $T")
end

############################ ConstantMass ######################################

struct ConstantMass <: AbstractMass end

Base.@kwdef mutable struct ConstantMassU
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end

get_u0(::ConstantMass) = ConstantMassU()

f_cont!(::System{ConstantMass}) = nothing
f_disc!(::System{ConstantMass}) = false

function MassData(sys::System{ConstantMass})
    @unpack m, J_Ob_b, r_ObG_b = sys.u
    MassData(m, J_Ob_b, r_ObG_b)
end

########################## ControlMappings #######################

abstract type AbstractControlMapping end
#in a NoControlMapping control setup there are no aircraft controls! we act upon
#the subsystem's controls directly! this allows testing multiple subsystem
#configurations without the hassle of having to define a dedicated
#assign_control_inputs! method for each one
struct NoControlMapping <: AbstractControlMapping end


######################### TestAircraft ########################

Base.@kwdef struct TestAircraft{
                    Kin <: AbstractKinematics,
                    Ctrl <: AbstractControlMapping,
                    #SMac <: AbstractStateMachine,
                    Mass <: AbstractMass,
                    Aero <: AbstractAirframeComponent,
                    Pwp <: AbstractAirframeComponent,
                    Ldg <: AbstractAirframeComponent,
                    Srf <: AbstractAirframeComponent} <: AbstractComponent
    kin::Kin = KinLTF()
    ctrl::Ctrl = NoControlMapping()
    mass::Mass = ConstantMass()
    aero::Aero = NullAirframeComponent()
    pwp::Pwp = NullAirframeComponent()
    ldg::Ldg = NullAirframeComponent()
    srf::Srf = NullAirframeComponent()
end

struct TestAircraftY{MassY, AeroY, PwpY, LdgY, SrfY}
    dyn::DynData
    kin::KinData
    air::AirData
    #smac::SMacY
    mass::MassY
    aero::AeroY
    pwp::PwpY
    ldg::LdgY
    srf::SrfY
end

struct TestAircraftD{MassD, AeroD, PwpD, LdgD, SrfD}
    #smac::SMacD
    mass::MassD
    aero::AeroD
    pwp::PwpD
    ldg::LdgD
    srf::SrfD
end

#in contrast with X, D and Y, which must hold all their aircraft subsystems'
#counterparts, U is arbitrarily defined by us as a design choice, and will
#obviously be partially determined by the selected ControlMapping

#=
# example:
const my_pwp = AirframeGroup(left = EThruster(), right = EThruster())
const MyPwp = typeof(my_pwp)
struct MyMapping <: AbstractControlMapping
get_u0(::TestAircraft{MyMapping, Mass, MyPwp, Ldg}) where {Mass,Ldg} = ComponentVector(throttle = 0.0)
#and now we define: assign_control_method! dispatching on these type
function assign_control_inputs!(::TestAircraft{MyMapping, Mass, MyPwp, Ldg}) where {Mass,Ldg}
    ac.subsystems.pwp.left.throttle = ac.u.throttle
    ac.subsystems.pwp.right.throttle = ac.u.throttle
end
=#
struct NoAircraftU end

get_u0(::TestAircraft{Kin,NoControlMapping,Mass,Aero,Pwp,Ldg,Srf} where {Kin,Mass,Aero,Pwp,Ldg,Srf}) = NoAircraftU()

get_x0(ac::TestAircraft) = ComponentVector(
    kin = get_x0(ac.kin),
    # smac = get_x0(ac.smac),
    mass = get_x0(ac.mass),
    aero = get_x0(ac.aero),
    pwp = get_x0(ac.pwp),
    ldg = get_x0(ac.ldg),
    srf = get_x0(ac.srf)
    )

get_y0(ac::TestAircraft) = TestAircraftY(
    DynData(),
    KinData(),
    AirData(),
    # get_y0(ac.smac),
    get_y0(ac.mass),
    get_y0(ac.aero),
    get_y0(ac.pwp),
    get_y0(ac.ldg),
    get_y0(ac.srf)
    )

get_d0(ac::TestAircraft) = TestAircraftD(
    # get_d0(ac.smac),
    get_d0(ac.mass),
    get_d0(ac.aero),
    get_d0(ac.pwp),
    get_d0(ac.ldg),
    get_d0(ac.srf)
    )


const TestAircraftSys{K,C,M,A,P,L,S} = System{TestAircraft{K,C,M,A,P,L,S}} where {K,C,M,A,P,L,S}


function System(ac::TestAircraft, ẋ = get_x0(ac), x = get_x0(ac),
                    y = get_y0(ac), u = get_u0(ac), d = get_d0(ac), t = Ref(0.0))

    #each subsystem allocates its own u, then we can decide how the aircraft's u
    #should map onto it via assign_control_inputs!
    # smac = System(ac.smac, ẋ.smac, x.smac, y.smac, get_u0(ac.smac), d.smac, t)
    mass = System(ac.mass, ẋ.mass, x.mass, y.mass, get_u0(ac.mass), d.mass, t)
    aero = System(ac.aero, ẋ.aero, x.aero, y.aero, get_u0(ac.aero), d.aero, t)
    pwp = System(ac.pwp, ẋ.pwp, x.pwp, y.pwp, get_u0(ac.pwp), d.pwp, t)
    ldg = System(ac.ldg, ẋ.ldg, x.ldg, y.ldg, get_u0(ac.ldg), d.ldg, t)
    srf = System(ac.srf, ẋ.srf, x.srf, y.srf, get_u0(ac.srf), d.srf, t)

    params = ()

    subsystems = (mass = mass, aero = aero, pwp = pwp, ldg = ldg, srf = srf)

    System{map(typeof, (ac, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)
end


#for an aircraft with no specific mapping, the user is expected to act upon
#the subsystems inputs directly, so this function has nothing to do
#overriding this function for specific TestAircraft type parameter
#combinations allows us to customize how the aircraft's control inputs map
#into its subsystems's control inputs.
assign_control_inputs!(::TestAircraftSys{K,NoControlMapping,M,A,P,L,S} where {K,M,A,P,L,S}) = nothing


function f_cont!(ac_sys::TestAircraftSys, trn::AbstractTerrain, atm::AtmosphericSystem)

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

    ac_sys.y = TestAircraftY(dyn_data, kin_data, air_data, mass.y, aero.y, pwp.y, ldg.y, srf.y)

    return nothing
end


function f_disc!(ac_sys::TestAircraftSys)
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

function plots(t::AbstractVector{<:Real}, data::AbstractVector{<:TestAircraftY};
    mode::Symbol, save_path::Union{String, Nothing}, kwargs...)

    sa = StructArray(data)

    ss_path = (dyn = "dynamics", kin = "kinematics", air = "airdata", mass = "mass",
        aero = "aerodynamics", pwp = "powerplant", ldg = "landinggear", srf = "surfaces")

    for (label, path) ∈ zip(keys(ss_path), values(ss_path))
        plots(t, getproperty(sa, label);
                mode, save_path = mkpath(joinpath(save_path, path)), kwargs...)
    end

end



end