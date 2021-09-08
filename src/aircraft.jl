module Aircraft

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.Terrain
using Flight.Atmosphere

using Flight.System
using Flight.StateMachine
using Flight.Airframe
using Flight.Airdata
using Flight.Propulsion
# using Flight.LandingGear

using Flight.Kinematics
using Flight.Dynamics

import Flight.System: HybridSystem, get_x0, get_y0, get_u0, get_d0,f_cont!, f_disc!

using Flight.Plotting
import Flight.Plotting: plots

export TestAircraft, TestAircraftD, TestAircraftY
export NoMapping, ConstantMassModel

abstract type AbstractMassModel end

#given some inputs (typically state of the fuel system and external payloads),
#an AbstractMassModel returns a MassData struct (defined in the Dynamics
#module). for now, we can simply define a ConstantMassModel

Base.@kwdef struct ConstantMassModel <: AbstractMassModel
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end

get_mass_data(model::ConstantMassModel) = MassData(model.m, model.J_Ob_b, model.r_ObG_b)

abstract type AbstractControlMapping end
#in a NoMapping control setup there are no aircraft controls! we act upon the
#subsystem's controls directly! this allows testing multiple subsystem
#configurations without the hassle of having to define a dedicated
#assign_control_inputs! method for each one
struct NoMapping <: AbstractControlMapping end


struct TestAircraft{Ctl <: AbstractControlMapping,
                    StM <: AbstractStateMachine,
                    Mass <: AbstractMassModel,
                    Pwp <: AbstractAirframeComponent} <: AbstractComponent
    # Ldg} <: AbstractComponent
    ctl::Ctl
    stm::StM
    mass::Mass
    pwp::Pwp
    # ldg::Ldg
end

function TestAircraft()

    ctl = NoMapping()
    stm = NoStateMachine()
    mass = ConstantMassModel(m = 1, J_Ob_b = 1*Matrix{Float64}(I,3,3))
    pwp = ACGroup((
        left = EThruster(motor = ElectricMotor(α = CW)),
        right = EThruster(motor = ElectricMotor(α = CCW))))
    # ldg = ACGroup((
    #     lmain = LandingGearLeg(),
    #     rmain = LandingGearLeg(),
    #     nlg = LandingGearLeg()))

    # TestAircraft(mass, pwp, ldg)
    TestAircraft(ctl, stm, mass, pwp)
end

struct TestAircraftY{StmY, PwpY}
    dyn::DynY
    kin::KinY
    air::AirY
    stm::StmY
    pwp::PwpY
end

struct TestAircraftD{StmD, PwpD}
    stm::StmD
    pwp::PwpD
end

#in contrast with X, D and Y, which must hold all their aircraft subsystems'
#counterparts, U is arbitrarily defined by us as a design choice, and will
#obviously be partially determined by the selected ControlMapping
struct EmptyAircraftU end


function get_x0(ac::TestAircraft)
    ComponentVector(
        kin = get_x0(KinInit()),
        stm = get_x0(ac.stm),
        pwp = get_x0(ac.pwp)
        )
end
get_y0(ac::TestAircraft) = TestAircraftY(DynY(), KinY(), AirY(), get_y0(ac.stm), get_y0(ac.pwp))
get_d0(ac::TestAircraft) = TestAircraftD(get_d0(ac.stm), get_d0(ac.pwp))
get_u0(::TestAircraft{NoMapping,Mass,Pwp} where {Mass,Pwp}) = EmptyAircraftU()


#=
# example:
const my_pwp = ACGroup(left = EThruster(), right = EThruster())
const MyPwp = typeof(my_pwp)
struct MyMapping <: AbstractControlMapping
get_u0(::TestAircraft{MyMapping, Mass, MyPwp, Ldg}) where {Mass,Ldg} = ComponentVector(throttle = 0.0)
#and now we define: assign_control_method! dispatching on these type
function assign_control_inputs!(::TestAircraft{MyMapping, Mass, MyPwp, Ldg}) where {Mass,Ldg}
    ac.subsystems.pwp.left.throttle = ac.u.throttle
    ac.subsystems.pwp.right.throttle = ac.u.throttle
end
=#

const TestAircraftSys{C,S,M,P} = HybridSystem{TestAircraft{C,S,M,P}} where {C,S,M,P}


function HybridSystem(ac::TestAircraft, ẋ = get_x0(ac), x = get_x0(ac),
                    y = get_y0(ac), u = get_u0(ac), d = get_d0(ac), t = Ref(0.0))

    #each subsystem allocate its own u, then we can decide how the aircraft's u
    #should map onto it via assign_control_inputs!
    stm = HybridSystem(ac.stm, ẋ.stm, x.stm, y.stm, get_u0(ac.stm), d.stm, t)
    pwp = HybridSystem(ac.pwp, ẋ.pwp, x.pwp, y.pwp, get_u0(ac.pwp), d.pwp, t)
    # ldg = HybridSystem(ac.ldg, ẋ.ldg, x.ldg, d.ldg, get_u0(ac.ldg), t)
    params = (mass = ac.mass,)
    # subsystems = (pwp = pwp, ldg = ldg)
    subsystems = (stm = stm, pwp = pwp,)
    HybridSystem{map(typeof, (ac, x, y, u, d, params, subsystems))...}(ẋ, x, y, u, d, t, params, subsystems)
end

#for an aircraft with no specific mapping, the user is expected to act upon
#the subsystems inputs directly, so this function has nothing to do
#overriding this function for specific TestAircraft type parameter
#combinations allows us to customize how the aircraft's control inputs map
#into its subsystems's control inputs.
assign_control_inputs!(::TestAircraftSys{NoMapping,S,M,P} where {S,M,P}) = nothing


function f_cont!(ac_sys::TestAircraftSys{C,S,M,P} where {C,S,M,P},
                trn::AbstractTerrainModel,
                atm::AbstractAtmosphericModel)


    @unpack ẋ, x, u, params, subsystems = ac_sys
    @unpack stm, pwp = subsystems
    @unpack mass = params

    y_kin = f_kin!(ẋ.kin.pos, x.kin)

    #before updating the subsystems, assign their control inputs from the
    #aircraft's control input vector
    assign_control_inputs!(ac_sys)

    y_air = AirY() #replacej
    # y.air .= get_air_data(). #call air data system here to update air data, passing also as
    # argument data.atmospheric_model

    #get aerodynamics Wrench
    # y_aero = get_wr_b(Aero, y.air, y.srf, y.ldg, trn)

    y_stm = stm.y

    #we don't need to update or extract each subsystem's ẋ, x or y, because
    #they are either views or mutable fields
    f_cont!(pwp, y_air)
    y_pwp = pwp.y
    # #update landing gear
    # f_cont!(y.ldg, ẋ.ldg, x.ldg, u.ldg, t, Ldg, trn)

    #initialize external Wrench and additional angular momentum
    wr_ext_b = Wrench()
    hr_b = SVector(0.,0.,0.)

    #add powerplant contributions
    wr_ext_b += get_wr_b(pwp)
    hr_b += get_hr_b(pwp)

    # #add landing gear contributions
    # wr_ext_b .+= get_wr_b(y.ldg, Ldg)
    # hr_b += get_hr_b(y.ldg, Ldg)

    #mass data depends on the state of the systems, we need updated y to compute
    #it, so it should go after the systems
    mass_data = get_mass_data(mass)

    #update dynamics
    y_dyn = f_dyn!(ẋ.kin.vel, wr_ext_b, hr_b, mass_data, y_kin)

    ac_sys.y = TestAircraftY(y_dyn, y_kin, y_air, y_stm, y_pwp)

    return nothing
end


function f_disc!(ac_sys::TestAircraftSys{C,M,P} where {C,M,P})
    x_mod = renormalize!(ac_sys.x.kin, 1e-8)
    # println(x_mod)
    return x_mod
end

function plots(t::AbstractVector{<:Real}, data::AbstractVector{<:TestAircraftY};
    mode::Symbol, save_path::Union{String, Nothing}, kwargs...)

    sa = StructArray(data)

    plots(t, sa.dyn; mode, save_path = mkpath(joinpath(save_path, "dynamics")), kwargs...)
    plots(t, sa.kin; mode, save_path = mkpath(joinpath(save_path, "kinematics")), kwargs...)
    plots(t, sa.pwp; mode, save_path = mkpath(joinpath(save_path, "powerplant")), kwargs...)
end



end