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
import Flight.System: HybridSystem, X, D, U, f_cont!, f_disc!
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

#there are no aircraft's "own" continuous states (if there were, we could always
#create a subsystem to accomodate them). however, we cannot say the same of
#discrete states. there may be some discrete logic in the aircraft which has no
#continuous states. and if it has no continuous states, we cannot store it
#within any subsystem (all systems are required to have at least continuous
#states). for these discrete states we can simply allocate a field in the NT
#returned by D
X(ac::TestAircraft) = ComponentVector(kin = X(Kin()), pwp = X(ac.pwp))
D(ac::TestAircraft) = (pwp = D(ac.pwp), )
U(::TestAircraft{NoControlMapping,Mass,Pwp} where {Mass,Pwp}) = nothing

#in a NoMapping aircraft, there are no aircraft controls! we act upon the
#subsystem's controls directly! this allows testing multiple subsystem
#configurations without the hassle of having to define a dedicated
#assign_control_inputs! method for each one

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

const TestAircraftSys{C,M,P} = HybridSystem{TestAircraft{C,M,P}} where {C,M,P}

function HybridSystem(ac::TestAircraft, ẋ = X(ac), x = X(ac), d = D(ac), u = U(ac), t = Ref(0.0))
    #each subsystem allocate its own u, then we can decide how the aircraft's u
    #should map onto it via assign_control_inputs!
    pwp = HybridSystem(ac.pwp, ẋ.pwp, x.pwp, d.pwp, U(ac.pwp), t)
    # ldg = HybridSystem(ac.ldg, ẋ.ldg, x.ldg, d.ldg, U(ac.ldg), t)
    params = (mass = ac.mass,)
    # subsystems = (pwp = pwp, ldg = ldg)
    subsystems = (pwp = pwp,)
    HybridSystem{map(typeof, (ac, x, d, u, params, subsystems))...}(ẋ, x, d, u, t, params, subsystems)
end

#for an aircraft with no specific mapping, the user is expected to act upon
#the subsystems inputs directly, so this function has nothing to do
#overriding this function for specific TestAircraft type parameter
#combinations allows us to customize how the aircraft's control inputs map
#into its subsystems's control inputs.
assign_control_inputs!(::TestAircraftSys{NoControlMapping,M,P} where {M,P}) = nothing


function f_cont!(ac_sys::TestAircraftSys{C,M,P} where {C,M,P},
                trn::AbstractTerrainModel,
                atm::AbstractAtmosphericModel)


    @unpack ẋ, x, u, params, subsystems = ac_sys
    @unpack pwp = subsystems
    @unpack mass = params

    y_kin = f_kin!(ẋ.kin.pos, x.kin)

    #before updating the subsystems, assign their control inputs from the
    #aircraft's control input vector
    assign_control_inputs!(ac_sys)

    y_air = AirDataY() #replacej
    # y.air .= get_air_data(). #call air data system here to update air data, passing also as
    # argument data.atmospheric_model

    #get aerodynamics Wrench
    # y_aero = get_wr_Ob_b(Aero, y.air, y.srf, y.ldg, trn)

    #we don't need to update or extract each subsystem's ẋ, x or y, because
    #they are either views or mutable fields
    y_pwp = f_cont!(pwp, y_air)
    # #update landing gear
    # f_cont!(y.ldg, ẋ.ldg, x.ldg, u.ldg, t, Ldg, trn)

    #initialize external Wrench and additional angular momentum
    wr_ext_Ob_b = Wrench()
    h_rot_b = SVector(0.,0.,0.)

    #add powerplant contributions
    wr_ext_Ob_b += get_wr_Ob_b(y_pwp)
    h_rot_b += get_h_Gc_b(y_pwp)

    # #add landing gear contributions
    # wr_ext_Ob_b .+= get_wr_Ob_b(y.ldg, Ldg)
    # h_rot_b += get_h_Gc_b(y.ldg, Ldg)

    #mass data depends on the state of the systems, we need updated y to compute
    #it, so it should go after the systems
    mass_data = get_mass_data(mass)

    #update dynamics
    y_acc = f_dyn!(ẋ.kin.vel, wr_ext_Ob_b, h_rot_b, mass_data, y_kin)

# Y(ac::TestAircraft) = ComponentVector(kin = Y(Kin()), acc = Y(Acc()), air =
# Y(AirData()), pwp = Y(ac.pwp))
    return (kin = y_kin, acc = y_acc, air = y_air, pwp = y_pwp)
end


function f_disc!(ac_sys::TestAircraftSys{C,M,P} where {C,M,P})
    x_mod = renormalize!(ac_sys.x.kin, 1e-8)
    # println(x_mod)
    return x_mod
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