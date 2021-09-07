module Aircraft

using StaticArrays: SVector, SMatrix
using LinearAlgebra
using ComponentArrays
using RecursiveArrayTools
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


#there are no aircraft's "own" continuous states (if there were, we could always
#create a subsystem to accomodate them). however, we cannot say the same of
#discrete states. there may be some discrete logic (a state machine) in the
#aircraft which has no continuous states. and if it has no continuous states, we
#cannot store it within any subsystem (all systems are required to have at least
#continuous states). we reserve a state machine field for this. the state
#machine only has discrete states of its own, but every time its f_disc! is
#called, it can (and should) accept the aircraft's continuous state (which it
#should not modify) and its input vector (which it shouldn't modify either,
#because if it is an internally modifiable input vector, it is not actually an
#input vector). the state machine, like

struct TestAircraftY{StmY, PwpY}
    kin::KinY
    acc::AccY
    air::AirY
    stm::StmY
    pwp::PwpY
end


struct TestAircraftD{StmD, PwpD}
    stm::StmD
    pwp::PwpD
end

#the crucial different wrt X, D and Y is that U is arbitrarily defined by
#ourselves, it is not determined by the aircraft subsystems (although obviously
#it should be defined taking them into account)

#here we should check which subsystems are hybrid (stateful), and add only those
#as x0 blocks
function get_x0(ac::TestAircraft)
    (
        kin = System.assemble_x0(Kin())[1],
        stm = System.assemble_x0(ac.stm)[1],
        pwp = System.assemble_x0(ac.pwp)[1]
        )
end

get_y0(ac::TestAircraft) = TestAircraftY(KinY(), AccY(), AirY(), get_y0(ac.stm), get_y0(ac.pwp))
get_d0(ac::TestAircraft) = TestAircraftD(get_d0(ac.stm), get_d0(ac.pwp))


#replace this with EmptyU <: AbstractU{EmptyComponent}
struct EmptyAircraftU end
get_u0(::TestAircraft{NoMapping,Mass,Pwp} where {Mass,Pwp}) = EmptyAircraftU()

#in a NoMapping aircraft, there are no aircraft controls! we act upon the
#subsystem's controls directly! this allows testing multiple subsystem
#configurations without the hassle of having to define a dedicated
#assign_control_inputs! method for each one

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
    (ẋ_s, ẋ_ss) = System.assemble_x0(ẋ)
    (x_s, x_ss) = System.assemble_x0(x)

    stm = HybridSystem(ac.stm, ẋ_ss.stm, x_ss.stm, y.stm, get_u0(ac.stm), d.stm, t)
    pwp = HybridSystem(ac.pwp, ẋ_ss.pwp, x_ss.pwp, y.pwp, get_u0(ac.pwp), d.pwp, t)
    # ldg = HybridSystem(ac.ldg, ẋ.ldg, x.ldg, d.ldg, get_u0(ac.ldg), t)
    params = (mass = ac.mass,)
    # subsystems = (pwp = pwp, ldg = ldg)
    subsystems = (stm = stm, pwp = pwp,)
    HybridSystem{map(typeof, (ac, x_s, y, u, d, params, subsystems))...}(ẋ_s, x_s, y, u, d, t, params, subsystems)
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
    y_acc = f_dyn!(ẋ.kin.vel, wr_ext_b, hr_b, mass_data, y_kin)

    ac_sys.y = TestAircraftY(y_kin, y_acc, y_air, y_stm, y_pwp)

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
    kin_data = sa.kin
    acc_data = sa.acc
    pwp_data = sa.pwp

    #put kinematics and acceleration outputs all in a single airframe folder
    plots(t, kin_data; mode, save_path, kwargs...)
    # plots(t, acc_data; mode, save_path, kwargs...)
    plots(t, pwp_data; mode, save_path, kwargs...)
end

######### Example: extracting y fields for plotting

# Base.@kwdef struct Output{Y, YD} #this is the type we pass to SavedValues
# y::Y = ComponentVector(a = fill(1.0, 2), b = fill(2.0, 3))
# yd::YD = ComponentVector(m = fill(4,3), n = fill(-2, 2))
# end
# log = collect(Output() for i in 1:5) #create some copies of it
# sa = StructArray(log) #now all the y fields of log lie in a contiguous array
# y = sa.y
# y_voa = VectorOfArray(y) #now we have a vector of Y's that indexes like a matrix
# y_mat = convert(Array, y_voa) #and a matrix of y's whose rows still preserve axis metadata
# function plotlog(log, aircraft::ParametricAircraft)

#############

#     y = log.y

#     #this could all be delegated to kin, which can return a handle to each plot
#     #it produces
#     #who saves the plots? how is the folder hierarchy generated?
#     kin = y[:kin, :]
#     pos = kin[:pos, :]
#     Δx = pos[:Δx, :]
#     Δy = pos[:Δy, :]
#     h = pos[:h, :]

#     #the only thing we ask of the log type is that it has fields :t and :saveval
#     #now we would construct a NamedTuple to delegate
#     return((log.t, h))

# end


end