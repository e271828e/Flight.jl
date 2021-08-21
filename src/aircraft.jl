module Aircraft

using StaticArrays: SVector, SMatrix
using LinearAlgebra
using ComponentArrays
using RecursiveArrayTools
using UnPack

using Flight.System
using Flight.Kinematics
using Flight.Dynamics
using Flight.Airdata
using Flight.Airframe
using Flight.Propulsion
# using Flight.LandingGear
using Flight.Terrain
using Flight.Atmosphere
import Flight.System: ContinuousSystem, X, Y, U, f_cont!, f_disc!
# import Flight.System: plotlog

export TestAircraft

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
struct NoControlMapping <: AbstractControlMapping end

struct TestAircraft{
    Ctl <: AbstractControlMapping,
    Mass <: AbstractMassModel,
    Pwp <: AbstractAirframeComponent} <: AbstractComponent
    # Ldg <: AbstractAirframeComponent} <: AbstractComponent
    # pwp::Pwp <: AbstractAirframeComponent
    ctl::Ctl
    mass::Mass
    pwp::Pwp
    # ldg::Ldg <: AbstractAirframeComponent
end

function TestAircraft()

    ctl = NoControlMapping()
    mass = ConstantMassModel(m = 1, J_Ob_b = 1*Matrix{Float64}(I,3,3))
    pwp = AirframeComponentGroup((
        left = EThruster(motor = ElectricMotor(α = CW)),
        right = EThruster(motor = ElectricMotor(α = CCW))))
    # ldg = AirframeComponentGroup((
    #     lmain = LandingGearLeg(),
    #     rmain = LandingGearLeg(),
    #     nlg = LandingGearLeg()))

    # TestAircraft(mass, pwp, ldg)
    TestAircraft(ctl, mass, pwp)
end

# X(ac::TestAircraft) = ComponentVector(kin = X(Kin()), pwp = X(ac.pwp), ldg = X(ac.ldg))
# Y(ac::TestAircraft) = ComponentVector(kin = Y(Kin()), acc = Y(Acc()), air = Y(AirData()), pwp = Y(ac.pwp), ldg = Y(ac.ldg))
# U(ac::TestAircraft) = (pwp = U(ac.pwp), ldg = U(ac.ldg))
X(ac::TestAircraft) = ComponentVector(kin = X(Kin()), pwp = X(ac.pwp))
Y(ac::TestAircraft) = ComponentVector(kin = Y(Kin()), acc = Y(Acc()), air = Y(AirData()), pwp = Y(ac.pwp))
U(::TestAircraft{NoControlMapping,Mass,Pwp} where {Mass,Pwp}) = nothing
#in a NoMapping aircraft, there are no aircraft controls! we act upon the
#subsystem's controls directly! this allows testing multiple subsystem
#configurations without the hassle of having to define a dedicated
#assign_control_inputs! method for each one

#for a direct mapping, we simply assemble the aircraft's u from each of its
#component's u's
# U(ac::TestAircraft{DirectMapping,Mass,Pwp} where {Mass,Pwp}) = (pwp = U(ac.pwp),)


#for an aircraft with no specific mapping, the user is expected to act upon
#the subsystems inputs directly, so this function has nothing to do
#overriding this function for specific TestAircraft type parameter
#combinations allows us to customize how the aircraft's control inputs map
#into its subsystems's control inputs.
assign_control_inputs!(::ContinuousSystem{TestAircraft{NoControlMapping,M,P}} where {M,P}) = nothing

#=
# example:
const my_pwp = AirframeComponentGroup(left = EThruster(), right = EThruster())
const MyPwp = typeof(my_pwp)
struct MyMapping <: AbstractControlMapping
U(::TestAircraft{MyMapping, Mass, MyPwp, Ldg}) where {Mass,Ldg} = ComponentVector(throttle = 0.0)
#and now we define: assign_control_method! dispatching on these type
function assign_control_inputs!(::TestAircraft{MyMapping, Mass, MyPwp, Ldg}) where {Mass,Ldg}
    ac.subsystems.pwp.left.throttle = ac.u.throttle
    ac.subsystems.pwp.right.throttle = ac.u.throttle
end
=#
function ContinuousSystem(ac::TestAircraft, ẋ = X(ac), x = X(ac), y = Y(ac), u = U(ac), t = Ref(0.0))
    #each subsystem allocate its own u, then we can decide how the aircraft's u
    #should map onto it via assign_control_inputs!
    pwp = ContinuousSystem(ac.pwp, ẋ.pwp, x.pwp, y.pwp, U(ac.pwp), t)
    # ldg = ContinuousSystem(ac.ldg, ẋ.ldg, x.ldg, y.ldg, u.ldg, t)
    params = (mass = ac.mass,)
    # subsystems = (pwp = pwp, ldg = ldg)
    subsystems = (pwp = pwp,)
    ContinuousSystem{map(typeof, (ac, x, y, u, params, subsystems))...}(ẋ, x, y, u, t, params, subsystems)
end


function f_cont!(ac_sy::ContinuousSystem{TestAircraft{C,M,P}} where {C,M,P},
                trn::AbstractTerrainModel,
                atm::AbstractAtmosphericModel)


    #update kinematics
    @unpack ẋ, x, y, u, params, subsystems = ac_sys
    @unpack pwp = subsystems
    @unpack mass = params

    f_kin!(y.kin, ẋ.kin.pos, x.kin)

    #before updating the subsystems, assign their control inputs from the
    #aircraft's control input vector
    assign_control_inputs!(ac_sys)


    # y.air .= get_air_data(). #call air data system here to update air data, passing also as
    # argument data.atmospheric_model

    #get aerodynamics Wrench
    # y_aero = get_wr_Ob_b(Aero, y.air, y.srf, y.ldg, trn)

    #we don't need to update or extract each subsystem's ẋ, x or y, because
    #they are either views or mutable fields
    f_cont!(pwp, y.air)
    # #update landing gear
    # f_cont!(y.ldg, ẋ.ldg, x.ldg, u.ldg, t, Ldg, trn)

    #initialize external Wrench and additional angular momentum
    wr_ext_Ob_b = Wrench()
    h_rot_b = SVector(0.,0.,0.)

    #add powerplant contributions
    wr_ext_Ob_b .+= get_wr_Ob_b(pwp)
    h_rot_b += get_h_Gc_b(pwp)

    # #add landing gear contributions
    # wr_ext_Ob_b .+= get_wr_Ob_b(y.ldg, Ldg)
    # h_rot_b += get_h_Gc_b(y.ldg, Ldg)

    #mass data depends on the state of the systems, we need updated y to compute
    #it, so it should go after the systems
    mass_data = get_mass_data(mass)

    #update dynamics
    f_dyn!(y.acc, ẋ.kin.vel, wr_ext_Ob_b, h_rot_b, mass_data, y.kin)

    return nothing
end

# degraded(nrm) = (abs(nrm - 1.0) > 1e-10)

# function f_disc!(x, u, t, aircraft::ParametricAircraft)

#     norm_q_lb = norm(x.kin.pos.q_lb)
#     norm_q_el = norm(x.kin.pos.q_el)
#     if degraded(norm_q_lb) || degraded(norm_q_el)
#         # println("Renormalized")
#         x.kin.pos.q_lb ./= norm_q_lb
#         x.kin.pos.q_el ./= norm_q_el
#         return true #x modified
#     else
#         return false #x not modified
#     end

# end

# function plotlog(log, aircraft::ParametricAircraft)

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