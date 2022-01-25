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
using Flight.Components
using Flight.Aerodynamics: AbstractAerodynamics
using Flight.Propulsion: EThruster, ElectricMotor, SimpleProp, CW, CCW
using Flight.LandingGear: LandingGearUnit, DirectSteering, DirectBraking, Strut, SimpleDamper
using Flight.Aircraft: AircraftBase, AbstractAircraftID
using Flight.Input: XBoxController, get_axis_value, is_released

import Flight.Modeling: init_x, init_y, init_u, init_d, f_cont!, f_disc!
import Flight.Plotting: plots
import Flight.Components: MassTrait, WrenchTrait, AngularMomentumTrait, get_wr_b, get_mp_b
import Flight.Aircraft: assign_joystick_inputs!

# export BeaverDescriptor

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



############################ Airframe ############################

##### Landing Gear ####

struct Ldg{L <: LandingGearUnit, R <: LandingGearUnit,
    N <: LandingGearUnit} <: SystemGroupDescriptor
    left::L
    right::R
    nose::N
end

WrenchTrait(::System{<:Ldg}) = HasWrench()
AngularMomentumTrait(::System{<:Ldg}) = HasNoAngularMomentum()

function Ldg()

    mlg_damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)
    nlg_damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)

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
            t_bs = FrameTransform(r = [1.27, 0, 2.004] , q = RQuat()),
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

#here we would define f_cont! to account for fuel consumption. also could add
#discrete states to define which tank(s) we're drawing from

##### Aerodynamics ####

Base.@kwdef struct Aero
    S::Float64 = 23.23 #wing area
    b::Float64 = 14.63 #wingspan
    c::Float64 = 1.5875 #mean aerodynamic chord
    t_ba::FrameTransform = FrameTransform() #airframe to aerodynamics frame
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
    coeffs::AeroC = AeroC() #aerodynamic coefficients
    wr_b::Wrench = Wrench() #aerodynamic Wrench, airframe
end

init_x(::Type{Aero}) = ComponentVector(α_filt = 0.0, β_filt = 0.0) #filtered airflow angles
init_y(::Type{Aero}) = AeroY()
init_u(::Type{Aero}) = AeroU()

MassTrait(::System{Aero}) = HasNoMass()
WrenchTrait(::System{Aero}) = HasWrench()
AngularMomentumTrait(::System{Aero}) = HasNoAngularMomentum()

get_wr_b(sys::System{Aero}) = sys.y.wr_b



####################### Update functions ###########################

function f_cont!(ctl::System{Controls}, ::System{Airframe},
                ::KinData, ::AirData, ::AbstractTerrain)

    #here, controls do nothing but update their output state. for a more complex
    #aircraft a continuous state-space autopilot implementation could go here
    @unpack throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
            pedals, brake_left, brake_right, flaps = ctl.u

    return ControlsY(; throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
                            pedals, brake_left, brake_right, flaps)

end

function f_cont!(sys::System{Aero}, pwp::System{Pwp},
    air::AirData, kinematics::KinData, terrain::AbstractTerrain)


    #for near-zero TAS, the airflow angles are likely to chatter between 0, -π
    #and π. this introduces noise in airflow angle derivatives and general
    #unpleasantness. to fix it we fade them in from zero to a safe threshold. as
    #for V, it's only used for non-dimensionalization. we just need to avoid
    #dividing by zero. for TAS < TAS_min, dynamic pressure will be close to zero
    #and therefore forces and moments will vanish anyway.

    @unpack ẋ, x, u, params = sys
    @unpack α_filt, β_filt = x
    @unpack S, b, c, t_ba, δe_max, δa_max, δr_max, δf_max, τ = params
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


    #OJO CON LAS ADIMENSIONALIZACIONS DE PQR y DE Cl, Cm, Cn, podrian no ser las
    #mismas que en Beaver

    aero_coeffs = get_aero_coeffs(; α, β, p_nd, q_nd, r_nd,
        δa, δr, δe, δf, α_dot_nd, β_dot_nd)

    q_as = get_stability_axes(α)
    F_aero_s = q * S * SVector{3,Float64}(-C_D, C_Y, -C_L)
    F_aero_a = q_as(F_aero_s)
    M_aero_a = q * S * SVector{3,Float64}(C_l * b, C_m * c, C_n * b)

    wr_a = Wrench(F_aero_a, M_aero_a)
    wr_b = t_ba(wr_a)

    sys.y = AeroY(; α, α_filt, α_filt_dot, β, β_filt, β_filt_dot,
        e, a, r, f, coeffs, wr_b)

end


f_disc!(::System{Controls}, ::System{Airframe}) = false
f_disc!(::System{Aero}) = false




c = RigidBody(894, SA[1285 0 0; 0 1825 0; 0 0 2667])
t_bc = FrameTransform(r = SVector{3}(0.056, 0, 0.582))
mp_b = MassProperties(c, t_bc)


MTOW = 1406

#aerodynamics (SI)
S = 16.165
b = 10.912
c = 1.494


end #module