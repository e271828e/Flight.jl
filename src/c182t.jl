module C182T

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack
using Unitful

using Flight.Modeling
using Flight.Plotting
using Flight.Misc
using Flight.Attitude
using Flight.Terrain
using Flight.Airdata
using Flight.Kinematics
using Flight.Dynamics
using Flight.Aerodynamics: AbstractAerodynamics
using Flight.Propulsion: EThruster, ElectricMotor, SimpleProp, CW, CCW
using Flight.LandingGear
using Flight.Aircraft: AircraftBase, AbstractAircraftID, AbstractAirframe
using Flight.Input: XBoxController, get_axis_value, is_released

import Flight.Modeling: init_x, init_y, init_u, init_d, f_cont!, f_disc!
import Flight.Plotting: plots
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_wr_b, get_mp_b
import Flight.Input: assign!

export C182TDescriptor

struct ID <: AbstractAircraftID end

############################## Controls #################################

struct Controls <: SystemDescriptor end

Base.@kwdef mutable struct ControlsU
    throttle::Bounded{Float64, 0, 1} = 0.0
    yoke_Δx::Bounded{Float64, -1, 1} = 0.0 #ailerons (+ bank right)
    yoke_x0::Bounded{Float64, -1, 1} = 0.0 #ailerons (+ bank right)
    yoke_Δy::Bounded{Float64, -1, 1} = 0.0 #elevator (+ pitch up)
    yoke_y0::Bounded{Float64, -1, 1} = 0.0 #elevator (+ pitch up)
    pedals::Bounded{Float64, -1, 1} = 0.0 #rudder and nose wheel (+ yaw right)
    brake_left::Bounded{Float64, 0, 1} = 0.0 #[0, 1]
    brake_right::Bounded{Float64, 0, 1} = 0.0 #[0, 1]
    flaps::Bounded{Float64, 0, 1} = 0.0 #[0, 1]
end

#const is essential when declaring type aliases!
Base.@kwdef struct ControlsY
    throttle::Float64
    yoke_Δx::Float64
    yoke_x0::Float64
    yoke_Δy::Float64
    yoke_y0::Float64
    pedals::Float64
    brake_left::Float64
    brake_right::Float64
    flaps::Float64
end

init_u(::Type{Controls}) = ControlsU()
init_y(::Type{Controls}) = ControlsY(zeros(SVector{9})...)


################################################################################
############################ Airframe ############################

############################## OEW #################################

#dummy subsystem to add OEW mass properties to the airframe; could instead
#customize the get_mp_b method for System{Airframe} but this way we can fall
#back on the default System{SystemGroupDescriptor} implementation for get_mp_b,
#which requires all subsystems to define their own get_mp_b methods
struct OEW <: SystemDescriptor end

const mp_b_OEW = let
    #define the empty airframe as a RigidBody
    OEW_G = RigidBody(894, SA[1285 0 0; 0 1825 0; 0 0 2667])
    #define the transform from the origin of the airframe reference frame (Ob)
    #to the empty airframe's center of mass (G)
    t_Ob_G = FrameTransform(r = SVector{3}(0.056, 0, 0.582))
    #compute the empty airframe mass properties at Ob
    MassProperties(OEW_G, t_Ob_G)
end

MassTrait(::System{OEW}) = HasMass()
WrenchTrait(::System{OEW}) = HasNoWrench()
AngularMomentumTrait(::System{OEW}) = HasNoAngularMomentum()

get_mp_b(::System{OEW}) = mp_b_OEW

##################################### Pwp ######################################

struct Pwp <: SystemGroupDescriptor
    left::EThruster
    right::EThruster
end

MassTrait(::System{Pwp}) = HasNoMass()
WrenchTrait(::System{Pwp}) = HasWrench()
AngularMomentumTrait(::System{Pwp}) = HasAngularMomentum()

function Pwp()

    prop = SimpleProp(kF = 2e-3, J = 0.25)

    left = EThruster(propeller = prop, motor = ElectricMotor(α = CW))
    right = EThruster(propeller = prop, motor = ElectricMotor(α = CCW))

    Pwp(left, right)

end

##### Landing Gear ####

#more flexible, but requires adding a type parameter for the ldg field in the
#Airframe. if we simply declare ldg::Ldg, Ldg is a UnionAll and the System
#constructor (appropriately) fails. also, requires System{<:Ldg} instead of
#System{Ldg} in method signatures
# struct Ldg{L <: LandingGearUnit, R <: LandingGearUnit,
#     N <: LandingGearUnit} <: SystemGroupDescriptor
#     left::L
#     right::R
#     nose::N
# end

#alternative, less flexible (implies type redefinition if the type parameters of
#the field types change, and Revise will complain)
struct Ldg <: SystemGroupDescriptor
    left::LandingGearUnit{Strut{SimpleDamper}, NoSteering, DirectBraking}
    right::LandingGearUnit{Strut{SimpleDamper}, NoSteering, DirectBraking}
    nose::LandingGearUnit{Strut{SimpleDamper}, DirectSteering, NoBraking}
end

MassTrait(::System{Ldg}) = HasNoMass()
WrenchTrait(::System{Ldg}) = HasWrench()
AngularMomentumTrait(::System{Ldg}) = HasNoAngularMomentum()

function Ldg()

    mlg_damper = SimpleDamper(k_s = 25000, k_d_ext = 2000, k_d_cmp = 2000)
    nlg_damper = SimpleDamper(k_s = 25000, k_d_ext = 2000, k_d_cmp = 2000)

    left = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-0.381, -1.092, 1.902], q = RQuat() ),
            l_0 = 0.0,
            damper = mlg_damper),
        braking = DirectBraking())

    right = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-0.381, 1.092, 1.902], q = RQuat() ),
            l_0 = 0.0,
            damper = mlg_damper),
        braking = DirectBraking())

    nose = LandingGearUnit(
        strut = Strut(
            # t_bs = FrameTransform(r = [1.27, 0, 2.004] , q = RQuat()),
            t_bs = FrameTransform(r = [1.27, 0, 1.9] , q = RQuat()),
            l_0 = 0.0,
            damper = nlg_damper),
        steering = DirectSteering())

    Ldg(left, right, nose)

end

##### Payload ####

struct Payload <: SystemDescriptor
    pilot::MassProperties #mp_b
    copilot::MassProperties #mp_b
    baggage::MassProperties

    function Payload( ; pilot::AbstractMassDistribution = PointMass(75),
                        copilot::AbstractMassDistribution = PointMass(75),
                        baggage::AbstractMassDistribution = PointMass(50))

        pilot_slot = FrameTransform(r = SVector{3}(0.183, -0.356, 0.899))
        copilot_slot = FrameTransform(r = SVector{3}(0.183, 0.356, 0.899))
        baggage_slot = FrameTransform(r = SVector{3}(-1.316, 0, 0.899))

        return new( MassProperties(pilot, pilot_slot),
                    MassProperties(copilot, copilot_slot),
                    MassProperties(baggage, baggage_slot))

    end
end

Base.@kwdef mutable struct PayloadD
    pilot::Bool = true
    copilot::Bool = true
    baggage::Bool = true
end

init_d(::Type{Payload}) = PayloadD()

MassTrait(::System{Payload}) = HasMass()
WrenchTrait(::System{Payload}) = HasNoWrench()
AngularMomentumTrait(::System{Payload}) = HasNoAngularMomentum()

function get_mp_b(sys::System{Payload})
    mp_b = MassProperties()
    sys.d.pilot ? mp_b += sys.params.pilot : nothing
    sys.d.copilot ? mp_b += sys.params.copilot : nothing
    sys.d.baggage ? mp_b += sys.params.baggage : nothing
    return mp_b
end


##### Fuel ####

struct Fuel <: SystemDescriptor end

init_x(::Type{Fuel}) = ComponentVector(m_left = 0.0, m_right = 0.0) #fuel tank contents

MassTrait(::System{Fuel}) = HasMass()
WrenchTrait(::System{Fuel}) = HasNoWrench()
AngularMomentumTrait(::System{Fuel}) = HasNoAngularMomentum()

function get_mp_b(sys::System{Fuel})
    mp_b = MassProperties()

    frame_left = FrameTransform(r = SVector{3}(0.325, -2.845, 0))
    frame_right = FrameTransform(r = SVector{3}(0.325, 2.845, 0))

    m_left = PointMass(sys.x.m_left)
    m_right = PointMass(sys.x.m_right)

    mp_b += MassProperties(m_left, frame_left)
    mp_b += MassProperties(m_right, frame_right)

    return mp_b
end

f_cont!(::System{Fuel}, ::System{Pwp}) = nothing

#here we would define f_cont! to account for fuel consumption. also could add
#discrete states to define which tank(s) we're drawing from

##### Aerodynamics ####


Base.@kwdef struct Aero <: AbstractAerodynamics
    S::Float64 = 16.165 #wing area
    b::Float64 = 10.912 #wingspan
    c::Float64 = 1.494 #mean aerodynamic chord
    δe_max::Float64 = 30 |> deg2rad #maximum elevator deflection (rad)
    δa_max::Float64 = 40 |> deg2rad #maximum (combined) aileron deflection (rad)
    δr_max::Float64 = 30 |> deg2rad #maximum rudder deflection (rad)
    δf_max::Float64 = 30 |> deg2rad #maximum flap deflection (rad)
    τ::Float64 = 0.1 #time constant for airflow angle filtering
end

Base.@kwdef mutable struct AeroU
    e::Bounded{Float64, -1, 1} = 0.0 #elevator control input (+ pitch down)
    a::Bounded{Float64, -1, 1} = 0.0 #aileron control input (+ roll left)
    r::Bounded{Float64, -1, 1} = 0.0 #rudder control input (+ yaw left)
    f::Bounded{Float64, 0, 1} = 0.0 # flap control input (+ flap down)
end

# const aero_derivatives =
#this should be a namedtuple with fields C_L, C_D, C_Y. Each one in turn should
#be a NamedTuple, and it should contain C_L.α, C_L.q, each of which can be a
#constant or an interpolator

#then within get_aero_coeffs we unpack everything and evaluate each interpolator
#with its appropriate independent variables.


Base.@kwdef struct AeroC
    C_L::Float64 = 0.0
    C_D::Float64 = 0.0
    C_Y::Float64 = 0.0
    C_l::Float64 = 0.0
    C_m::Float64 = 0.0
    C_n::Float64 = 0.0
end

Base.@kwdef struct AeroY
    α::Float64 = 0.0 #preprocessed AoA
    α_filt::Float64 = 0.0 #filtered AoA
    α_filt_dot::Float64 = 0.0 #filtered AoA derivative
    β::Float64 = 0.0 #preprocessed AoS
    β_filt::Float64 = 0.0 #filtered AoS
    β_filt_dot::Float64 = 0.0 #filtered AoS derivative
    e::Float64 = 0.0 #normalized elevator control input
    a::Float64 = 0.0 #normalized aileron control input
    r::Float64 = 0.0 #normalized rudder control input
    f::Float64 = 0.0 #normalized flap control input
    coefs::AeroC = AeroC() #aerodynamic coefficients
    wr_b::Wrench = Wrench() #aerodynamic Wrench, airframe
end

init_x(::Type{Aero}) = ComponentVector(α_filt = 0.0, β_filt = 0.0) #filtered airflow angles
init_y(::Type{Aero}) = AeroY()
init_u(::Type{Aero}) = AeroU()


################################ Airframe ######################################

Base.@kwdef struct Airframe <: AbstractAirframe
    oew::OEW = OEW()
    aero::Aero = Aero()
    pwp::Pwp = Pwp()
    ldg::Ldg = Ldg()
    fuel::Fuel = Fuel()
    pld::Payload = Payload()
end

MassTrait(::System{Airframe}) = HasMass()
WrenchTrait(::System{Airframe}) = HasWrench()
AngularMomentumTrait(::System{Airframe}) = HasAngularMomentum()


################################################################################
####################### Update functions #######################################


####################### REMOVE THIS ########################################
####################### REMOVE THIS ########################################
####################### REMOVE THIS ########################################
####################### REMOVE THIS ########################################
####################### REMOVE THIS ########################################
####################### REMOVE THIS ########################################

f_cont!(sys::System{Aero}, pwp::System{Pwp},
    air::AirData, kinematics::KinData, terrain::AbstractTerrain) = nothing

"""
function f_cont!(sys::System{Aero}, pwp::System{Pwp},
    air::AirData, kinematics::KinData, terrain::AbstractTerrain)

    #for this aircraft we have chosen the aerodynamics frame fa as the main
    #aircraft reference frame fb. all the AirData variables from the AirData
    #module are computed in frame fb, they are directly applicable here. if this
    #wasn't the case, either because the origin Oa is not coincident with Ob and
    #we care about the velocity lever arm or because the axes εa and εb are
    #different, we would need to compute v_wOa_a from v_wOb_b and then recompute
    #all the derived air data variables

    # v_wOb_b = v_eOb_b - v_ew_b
    # v_eOa_b = v_eOb_b + ω_eb_b × r_ObOa_b
    # v_wOa_b = v_eOa_b - v_ew_b = v_eOb_b + ω_eb_b × r_ObOa_b - v_ew_b
    # v_wOa_b = v_wOb_b + ω_eb_b × r_ObOa_b
    # v_wOa_a = q_ba'(v_wOa_b)

    #for near-zero TAS, the airflow angles are likely to chatter between 0, -π
    #and π. this introduces noise in airflow angle derivatives and general
    #unpleasantness. to fix it we fade them in from zero to a safe threshold. as
    #for V, it's only used for non-dimensionalization. we just need to avoid
    #dividing by zero. for TAS < TAS_min, dynamic pressure will be close to zero
    #and therefore forces and moments will vanish anyway.

    @unpack ẋ, x, u, params = sys
    @unpack α_filt, β_filt = x
    @unpack S, b, c, δe_max, δa_max, δr_max, δf_max, τ = params
    @unpack e, a, r, f = u
    @unpack TAS, q, α_b, β_b = air
    ω_lb_b = kinematics.vel.ω_lb_b

    #preprocess airflow angles and airspeed. looking at the lift/drag polar +/-1
    #seem like reasonable airflow angle limits
    TAS_min = 2.0
    χ_TAS = min(TAS/TAS_min, 1.0) #linear fade-in for airflow angles
    α = clamp(α_b * χ_TAS, -1.0, 1.0)
    β = clamp(β_b * χ_TAS, -1.0, 1.0)
    V = max(TAS, TAS_min)

    α_filt_dot = 1/τ * (α - α_filt)
    β_filt_dot = 1/τ * (β - β_filt)

    p_nd = ω_lb_b[1] * b / (2V) #non-dimensional roll rate
    q_nd = ω_lb_b[2] * c / (2V) #non-dimensional pitch rate
    r_nd = ω_lb_b[3] * b / (2V) #non-dimensional yaw rate

    α_dot_nd = α_filt_dot * c / (2V) #not used
    β_dot_nd = β_filt_dot * b / (2V)

    ẋ.α_filt = α_filt_dot
    ẋ.β_filt = β_filt_dot

    # T = get_wr_b(pwp).F[1]
    # C_T = T / (q * S) #thrust coefficient, not used
    δe = Float64(e) * δe_max
    δa = Float64(a) * δa_max
    δr = Float64(r) * δr_max
    δf = Float64(f) * δf_max

    aero_coeffs = get_aero_coeffs(; α, β, p_nd, q_nd, r_nd,
        δa, δr, δe, δf, α_dot_nd, β_dot_nd)

    q_as = get_stability_axes(α)
    F_aero_s = q * S * SVector{3,Float64}(-C_D, C_Y, -C_L)
    F_aero_a = q_as(F_aero_s)
    M_aero_a = q * S * SVector{3,Float64}(C_l * b, C_m * c, C_n * b)

    wr_b = wr_a = Wrench(F_aero_a, M_aero_a)

    sys.y = AeroY(; α, α_filt, α_filt_dot, β, β_filt, β_filt_dot,
        e, a, r, f, coeffs, wr_b)

end
"""

get_wr_b(sys::System{Aero}) = sys.y.wr_b


######################## Controls Update Functions ###########################

function f_cont!(ctl::System{Controls}, ::System{Airframe},
                ::KinData, ::AirData, ::AbstractTerrain)

    #here, controls do nothing but update their output state. for a more complex
    #aircraft a continuous state-space autopilot implementation could go here
    @unpack throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
            pedals, brake_left, brake_right, flaps = ctl.u

    return ControlsY(; throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
                            pedals, brake_left, brake_right, flaps)

end

f_disc!(::System{Controls}, ::System{Airframe}) = false

################################################################################
####################### Airframe Update Functions ##############################

#we can't fall back  on the default System Group implementation, because of the
#interactions between the different subsystems

function f_cont!(afm::System{Airframe}, ctl::System{Controls},
                kin::KinData, air::AirData, trn::AbstractTerrain)

    @unpack aero, pwp, ldg, fuel, pld = afm.subsystems

    assign_component_inputs!(afm, ctl)
    f_cont!(ldg, kin, trn) #update landing gear continuous state & outputs
    f_cont!(pwp, kin, air) #update powerplant continuous state & outputs
    f_cont!(fuel, pwp) #update fuel system
    f_cont!(aero, pwp, air, kin, trn)

    # afm.y = (aero = aero.y, pwp = pwp.y, ldg = ldg.y, fuel = fuel.y, pld = pld.y )
    afm.y = (aero = aero.y, pwp = pwp.y, ldg = ldg.y)

end

function assign_component_inputs!(afm::System{<:Airframe}, ctl::System{<:Controls})

    @unpack throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
            pedals, brake_left, brake_right, flaps = ctl.u
    @unpack aero, pwp, ldg = afm.subsystems

    #yoke_Δx is the offset with respect to the force-free position yoke_x0
    #yoke_Δy is the offset with respect to the force-free position yoke_y0

    pwp.u.left.throttle = throttle
    pwp.u.right.throttle = throttle
    ldg.u.nose.steering[] = pedals
    ldg.u.left.braking[] = brake_left
    ldg.u.right.braking[] = brake_right
    aero.u.e = -(yoke_y0 + yoke_Δy) #+yoke_Δy and +yoke_y0 are back and +δe is pitch down, need to invert it
    aero.u.a = -(yoke_x0 + yoke_Δx) #+yoke_Δx and +yoke_x0 are right and +δa is roll left, need to invert it
    aero.u.r = -pedals # +pedals is right and +δr is yaw left
    aero.u.f = flaps # +flaps is flaps down and +δf is flaps down

    return nothing
end

function f_disc!(afm::System{Airframe}, ::System{Controls})
    #fall back to the default SystemGroup implementation
    return f_disc!(afm)
end

f_disc!(::System{OEW}) = false
f_disc!(::System{Aero}) = false
f_disc!(::System{Fuel}) = false
f_disc!(::System{Payload}) = false
#Pwp and Ldg are SystemGroups, so we can rely on the fallback f_disc!

#get_mp_b, get_wr_b and get_hr_b use the fallback for SystemGroups, which in turn call
#get_mp_b, get_wr_b and get_hr_b on aero, pwp and ldg


################################################################################
############################# Input Interfaces ################################


####################### XBoxController Input Interface ########################

elevator_curve(x) = exp_axis_curve(x, strength = 0.5, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 0.5, deadzone = 0.05)
pedal_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 0, deadzone = 0.05)

function exp_axis_curve(x::Bounded{T}, args...; kwargs...) where {T}
    exp_axis_curve(T(x), args...; kwargs...)
end

function exp_axis_curve(x::Real; strength::Real = 0.0, deadzone::Real = 0.0)

    a = strength
    x0 = deadzone

    abs(x) <= 1 || throw(ArgumentError("Input to exponential curve must be within [-1, 1]"))
    (x0 >= 0 && x0 <= 1) || throw(ArgumentError("Exponential curve deadzone must be within [0, 1]"))

    if x > 0
        y = max(0, (x - x0)/(1 - x0)) * exp( a * (abs(x) -1) )
    else
        y = min(0, (x + x0)/(1 - x0)) * exp( a * (abs(x) -1) )
    end
end

function assign!(ac::System{<:AircraftBase{ID}}, joystick::XBoxController)

    u = ac.u.controls

    u.yoke_Δx = get_axis_value(joystick, :right_analog_x) |> aileron_curve
    u.yoke_Δy = get_axis_value(joystick, :right_analog_y) |> elevator_curve
    u.pedals = get_axis_value(joystick, :left_analog_x) |> pedal_curve
    u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

    u.yoke_x0 -= 0.01 * is_released(joystick, :dpad_left)
    u.yoke_x0 += 0.01 * is_released(joystick, :dpad_right)
    u.yoke_y0 -= 0.01 * is_released(joystick, :dpad_up)
    u.yoke_y0 += 0.01 * is_released(joystick, :dpad_down)

    u.throttle += 0.1 * is_released(joystick, :button_Y)
    u.throttle -= 0.1 * is_released(joystick, :button_A)

    # u.propeller_speed += 0.1 * is_released(joystick, :button_X) #rpms
    # u.propeller_speed -= 0.1 * is_released(joystick, :button_B)

    u.flaps += 0.5 * is_released(joystick, :right_bumper)
    u.flaps -= 0.5 * is_released(joystick, :left_bumper)

    # Y si quisiera landing gear up y down, podria usar option como
    #modifier

end

#Aircraft constructor override keyword inputs to customize

function C182TDescriptor(; id = ID(), kin = KinLTF(), afm = Airframe(), ctl = Controls())
    AircraftBase( id; kin, afm, ctl)
end


end #module