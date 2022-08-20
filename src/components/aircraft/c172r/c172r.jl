module C172R

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack
using Reexport
using Interpolations
using HDF5
# using JLD

using Flight.Utils
using Flight.Input
using Flight.Systems
using Flight.Attitude
using Flight.Terrain
using Flight.Atmosphere
using Flight.Kinematics
using Flight.RigidBody
using Flight.LandingGear
using Flight.Propellers
using Flight.Piston
using Flight.Aircraft: AircraftBase, AbstractAirframe, AbstractAerodynamics, AbstractAvionics

import Flight.Systems: init, f_ode!, f_step!, f_disc!
import Flight.Kinematics: KinematicInit
import Flight.RigidBody: MassTrait, WrenchTrait, AngularMomentumTrait, get_wr_b, get_mp_b
import Flight.Piston: fuel_available
import Flight.Input: assign!

include("data/aero.jl")

export Cessna172R

############################## BasicAvionics #################################

struct BasicAvionics <: AbstractAvionics end

#elevator↑ (stick forward) -> e↑ -> δe↑ -> trailing edge down -> Cm↓ -> pitch down
#aileron↑ (stick right) -> a↑ -> δa↑ -> left trailing edge down, right up -> Cl↓ -> roll right
#pedals↑ (left pedal forward) -> r↓ -> δr↓ -> rudder trailing edge right -> Cn↑ -> yaw right
#pedals↑ (right pedal forward) -> nose wheel steering right -> yaw right
#flaps↑ -> δf↑ -> flap trailing edge down -> CL↑
Base.@kwdef mutable struct BasicAvionicsU
    throttle::Ranged{Float64, 0, 1} = 0.0
    aileron::Ranged{Float64, -1, 1} = 0.0
    Δ_aileron::Ranged{Float64, -1, 1} = 0.0 #incremental command, for input devices
    elevator::Ranged{Float64, -1, 1} = 0.0
    Δ_elevator::Ranged{Float64, -1, 1} = 0.0 #incremental command, for input devices
    pedals::Ranged{Float64, -1, 1} = 0.0
    Δ_pedals::Ranged{Float64, -1, 1} = 0.0 #incremental command, for input devices
    brake_left::Ranged{Float64, 0, 1} = 0.0
    brake_right::Ranged{Float64, 0, 1} = 0.0
    flaps::Ranged{Float64, 0, 1} = 0.0
    mixture::Ranged{Float64, 0, 1} = 0.5
    eng_start::Bool = false
    eng_stop::Bool = false
end

Base.@kwdef struct BasicAvionicsY
    throttle::Float64 = 0.0
    Δ_aileron::Float64 = 0.0
    aileron::Float64 = 0.0
    Δ_elevator::Float64 = 0.0
    elevator::Float64 = 0.0
    pedals::Float64 = 0.0
    Δ_pedals::Float64 = 0.0
    brake_left::Float64 = 0.0
    brake_right::Float64 = 0.0
    flaps::Float64 = 0.0
    mixture::Float64 = 0.5
    eng_start::Bool = false
    eng_stop::Bool = false
end

init(::SystemU, ::BasicAvionics) = BasicAvionicsU()
init(::SystemY, ::BasicAvionics) = BasicAvionicsY()


################################################################################
############################ Airframe Subsystems ################################

################################ Airframe ######################################

struct Structure <: Component end

# This component represents the airframe structure, together with any components
# rigidly attached to it, such as powerplant or landing gear, but not payload or
# fuel contents. Its mass corresponds roughly to the aircraft's Standard Empty
# Weight

#Structure mass properties computed in the vehicle reference frame b
const mp_b_str = let
    #define the structure as a RigidBodyDistribution
    str_G = RigidBodyDistribution(767.0, SA[820.0 0 0; 0 1164.0 0; 0 0 1702.0])
    #define the transform from the origin of the vehicle reference frame (Ob)
    #to the structure's center of mass (G)
    t_Ob_G = FrameTransform(r = SVector{3}(0.056, 0, 0.582))
    #compute the structure's mass properties at Ob
    MassProperties(str_G, t_Ob_G)
end

MassTrait(::System{Structure}) = HasMass()

#the structure itself receives no external actions. these are considered to act
#upon the vehicle's aerodynamics, power plant and landing gear. the same goes
#for rotational angular momentum.
WrenchTrait(::System{Structure}) = GetsNoExternalWrench()
AngularMomentumTrait(::System{Structure}) = HasNoAngularMomentum()

get_mp_b(::System{Structure}) = mp_b_str


############################ Aerodynamics ######################################

#the aircraft body reference frame fb is arbitrarily chosen to coincide with
#the aerodynamics frame fa, so the frame transform is trivial
const aero_lookup = generate_aero_lookup()
const f_ba = FrameTransform()

# if this weren't the case, and we cared not only about the rotation but also
#about the velocity lever arm, here's the rigorous way of computing v_wOa_a:
# v_wOb_b = v_eOb_b - v_ew_b
# v_eOa_b = v_eOb_b + ω_eb_b × r_ObOa_b
# v_wOa_b = v_eOa_b - v_ew_b = v_eOb_b + ω_eb_b × r_ObOa_b - v_ew_b
# v_wOa_b = v_wOb_b + ω_eb_b × r_ObOa_b
# v_wOa_a = q_ba'(v_wOa_b)

Base.@kwdef struct AeroCoeffs
    C_D::Float64 = 0.0
    C_Y::Float64 = 0.0
    C_L::Float64 = 0.0
    C_l::Float64 = 0.0
    C_m::Float64 = 0.0
    C_n::Float64 = 0.0
end

function get_aero_coeffs(; α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf, α_dot_nd, β_dot_nd, Δh_nd, stall)

    #set sensible bounds
    α = clamp(α, -0.1, 0.36) #0.36 is the highest value (post-stall) tabulated for C_L
    β = clamp(β, -0.2, 0.2)
    α_dot_nd = clamp(α_dot_nd, -0.04, 0.04)
    β_dot_nd = clamp(β_dot_nd, -0.2, 0.2)

    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = aero_lookup

    AeroCoeffs(
        C_D = C_D.z + C_D.ge(Δh_nd) * (C_D.α_δf(α,δf) + C_D.δf(δf)) + C_D.δe(δe) + C_D.β(β),
        C_Y = C_Y.δr * δr + C_Y.δa * δa + C_Y.β_δf(β,δf) + C_Y.p(α,δf) * p_nd + C_Y.r(α,δf) * r_nd,
        C_L = C_L.ge(Δh_nd) * (C_L.α(α,stall) + C_L.δf(δf)) + C_L.δe * δe + C_L.q * q_nd + C_L.α_dot * α_dot_nd,
        C_l = C_l.δa * δa + C_l.δr * δr + C_l.β * β + C_l.p * p_nd + C_l.r(α,δf) * r_nd,
        C_m = C_m.z + C_m.δe * δe + C_m.δf(δf) + C_m.α * α + C_m.q * q_nd + C_m.α_dot * α_dot_nd,
        C_n = C_n.δr * δr + C_n.δa * δa + C_n.β * β + C_n.p * p_nd + C_n.r * r_nd,
    )

end


Base.@kwdef struct Aero <: AbstractAerodynamics
    S::Float64 = 16.165 #wing area
    b::Float64 = 10.912 #wingspan
    c::Float64 = 1.494 #mean aerodynamic chord
    δe_range::NTuple{2,Float64} = deg2rad.((-28, 23)) #elevator deflection range (rad)
    δa_range::NTuple{2,Float64} = deg2rad.((-20, 20)) #aileron deflection range (rad)
    δr_range::NTuple{2,Float64} = deg2rad.((-16, 16)) #rudder deflection range (rad)
    δf_range::NTuple{2,Float64} = deg2rad.((0, 30)) #flap deflection range (rad)
    α_stall::NTuple{2,Float64} = (0.09, 0.36) #α values for stall hysteresis switching
    V_min::Float64 = 1.0 #lower airspeed threshold for non-dimensional angle rates
    τ::Float64 = 0.05 #time constant for filtered airflow angle derivatives
end

#e↑ -> δe↑ -> trailing edge down -> Cm↓ -> pitch down
#a↑ -> δa↑ -> left trailing edge down, right up -> Cl↓ -> roll right
#r↑ -> δr↑ -> rudder trailing edge left -> Cn↓ -> yaw left
#f↑ -> δf↑ -> flap trailing edge down -> CL↑

Base.@kwdef mutable struct AeroU
    e::Ranged{Float64, -1, 1} = 0.0
    a::Ranged{Float64, -1, 1} = 0.0
    r::Ranged{Float64, -1, 1} = 0.0
    f::Ranged{Float64, 0, 1} = 0.0
end

Base.@kwdef mutable struct AeroS #discrete state
    stall::Bool = false
end

Base.@kwdef struct AeroY
    e::Float64 = 0.0 #normalized elevator control input
    a::Float64 = 0.0 #normalized aileron control input
    r::Float64 = 0.0 #normalized rudder control input
    f::Float64 = 0.0 #normalized flap control input
    α::Float64 = 0.0 #clamped AoA, aerodynamic axes
    β::Float64 = 0.0 #clamped AoS, aerodynamic axes
    α_filt::Float64 = 0.0 #filtered AoA
    β_filt::Float64 = 0.0 #filtered AoS
    α_filt_dot::Float64 = 0.0 #filtered AoA derivative
    β_filt_dot::Float64 = 0.0 #filtered AoS derivative
    stall::Bool = false #stall state
    coeffs::AeroCoeffs = AeroCoeffs() #aerodynamic coefficients
    wr_b::Wrench = Wrench() #aerodynamic Wrench, vehicle frame
end

init(::SystemX, ::Aero) = init(SystemX(); α_filt = 0.0, β_filt = 0.0) #filtered airflow angles
init(::SystemY, ::Aero) = AeroY()
init(::SystemU, ::Aero) = AeroU()
init(::SystemS, ::Aero) = AeroS()


function f_ode!(sys::System{Aero}, ::System{<:Piston.Thruster},
    air::AirflowData, kinematics::KinematicData, terrain::System{<:AbstractTerrain})

    #for near-zero TAS, the airflow angles are likely to chatter between 0, -π
    #and π. this can cause lots of noise in airflow angle derivatives and
    #general unpleasantness. however, in this situation dynamic pressure will be
    #close to zero, so forces and moments will vanish anyway.

    @unpack ẋ, x, u, s, params = sys
    @unpack α_filt, β_filt = x
    @unpack e, a, r, f = u
    @unpack S, b, c, δe_range, δa_range, δr_range, δf_range, α_stall, V_min, τ = params
    @unpack TAS, q, v_wOb_b = air
    @unpack ω_lb_b, n_e, h_o = kinematics
    stall = s.stall

    v_wOb_a = f_ba.q'(v_wOb_b)
    α, β = get_airflow_angles(v_wOb_a)
    V = max(TAS, V_min) #avoid division by zero

    α_filt_dot = 1/τ * (α - α_filt)
    β_filt_dot = 1/τ * (β - β_filt)

    p_nd = ω_lb_b[1] * b / (2V) #non-dimensional roll rate
    q_nd = ω_lb_b[2] * c / (2V) #non-dimensional pitch rate
    r_nd = ω_lb_b[3] * b / (2V) #non-dimensional yaw rate

    α_dot_nd = α_filt_dot * c / (2V)
    β_dot_nd = β_filt_dot * b / (2V)

    δe = linear_scaling(e, δe_range)
    δa = linear_scaling(a, δa_range)
    δr = linear_scaling(r, δr_range)
    δf = linear_scaling(f, δf_range)

    #non-dimensional height above ground
    l2d_Oa = n_e #Oa = Ob
    h_Oa = h_o #orthometric
    h_trn_Oa = TerrainData(terrain, l2d_Oa).altitude #orthometric
    Δh_nd = (h_Oa - h_trn_Oa) / b

    # T = get_wr_b(pwp).F[1]
    # C_T = T / (q * S) #thrust coefficient, not used here

    coeffs = get_aero_coeffs(;
        α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf, α_dot_nd, β_dot_nd, Δh_nd, stall)

    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = coeffs

    q_as = get_stability_axes(α)
    F_aero_s = q * S * SVector{3,Float64}(-C_D, C_Y, -C_L)
    F_aero_a = q_as(F_aero_s)
    M_aero_a = q * S * SVector{3,Float64}(C_l * b, C_m * c, C_n * b)

    wr_b = wr_a = Wrench(F_aero_a, M_aero_a)

    ẋ.α_filt = α_filt_dot
    ẋ.β_filt = β_filt_dot

    sys.y = AeroY(; α, α_filt, α_filt_dot, β, β_filt, β_filt_dot,
        e, a, r, f, stall, coeffs, wr_b)

    return nothing

end

get_wr_b(sys::System{Aero}) = sys.y.wr_b

function f_step!(sys::System{Aero})
    #stall hysteresis
    α = sys.y.α
    α_stall = sys.params.α_stall
    if α > α_stall[2]
        sys.s.stall = true
    elseif α < α_stall[1]
        sys.s.stall = false
    end
    return false
end

# # splt_α = thplot(t, rad2deg.(α_b);
# #     title = "Angle of Attack", ylabel = L"$\alpha \ (deg)$",
# #     label = "", kwargs...)

# # splt_β = thplot(t, rad2deg.(β_b);
# #     title = "Angle of Sideslip", ylabel = L"$\beta \ (deg)$",
# #     label = "", kwargs...)

# # pd["05_α_β"] = plot(splt_α, splt_β;
# #     plot_title = "Airflow Angles [Airframe]",
# #     layout = (1,2),
# #     kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs



############################# Landing Gear ####################################

struct Ldg <: Component
    left::LandingGearUnit{NoSteering, DirectBraking, Strut{SimpleDamper}}
    right::LandingGearUnit{NoSteering, DirectBraking, Strut{SimpleDamper}}
    nose::LandingGearUnit{DirectSteering, NoBraking, Strut{SimpleDamper}}
end

MassTrait(::System{Ldg}) = HasNoMass()
WrenchTrait(::System{Ldg}) = GetsExternalWrench()
AngularMomentumTrait(::System{Ldg}) = HasNoAngularMomentum()

function Ldg()

    mlg_damper = SimpleDamper(k_s = 39404, #2700 lbf/ft
                              k_d_ext = 9340, #640 lbf/(ft/s)
                              k_d_cmp = 9340,
                              F_max = 50000)
    nlg_damper = SimpleDamper(k_s = 26269, #1800 lbf/ft
                              k_d_ext = 3503, #240 lbf/(ft/s)
                              k_d_cmp = 3503,
                              F_max = 50000)

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
            t_bs = FrameTransform(r = [1.27, 0, 1.9] , q = RQuat()),
            l_0 = 0.0,
            damper = nlg_damper),
        steering = DirectSteering())

    Ldg(left, right, nose)

end

############################## Payload ######################################

const pilot_slot = FrameTransform(r = SVector{3}(0.183, -0.356, 0.899))
const copilot_slot = FrameTransform(r = SVector{3}(0.183, 0.356, 0.899))
const psg_left_slot = FrameTransform(r = SVector{3}(-0.681, -0.356, 0.899))
const psg_right_slot = FrameTransform(r = SVector{3}(-0.681, 0.356, 0.899))
const baggage_slot = FrameTransform(r = SVector{3}(-1.316, 0, 0.899))

struct Payload <: Component
    pilot::MassProperties #mp_b
    copilot::MassProperties #mp_b
    psg_left::MassProperties #mp_b
    psg_right::MassProperties #mp_b
    baggage::MassProperties #mp_b

    function Payload( ; pilot::PointDistribution = PointDistribution(75),
                        copilot::PointDistribution = PointDistribution(75),
                        psg_left::PointDistribution = PointDistribution(75),
                        psg_right::PointDistribution = PointDistribution(75),
                        baggage::PointDistribution = PointDistribution(50))

        return new( MassProperties(pilot, pilot_slot),
                    MassProperties(copilot, copilot_slot),
                    MassProperties(psg_left, psg_left_slot),
                    MassProperties(psg_right, psg_right_slot),
                    MassProperties(baggage, baggage_slot)
                    )

    end
end

Base.@kwdef mutable struct PayloadU
    pilot::Bool = true
    copilot::Bool = true
    psg_left::Bool = false
    psg_right::Bool = false
    baggage::Bool = true
end

init(::SystemU, ::Payload) = PayloadU()

MassTrait(::System{Payload}) = HasMass()
WrenchTrait(::System{Payload}) = GetsNoExternalWrench()
AngularMomentumTrait(::System{Payload}) = HasNoAngularMomentum()

function get_mp_b(sys::System{Payload})
    mp_b = MassProperties()
    sys.u.pilot ? mp_b += sys.params.pilot : nothing
    sys.u.copilot ? mp_b += sys.params.copilot : nothing
    sys.u.psg_left ? mp_b += sys.params.psg_left : nothing
    sys.u.psg_right ? mp_b += sys.params.psg_right : nothing
    sys.u.baggage ? mp_b += sys.params.baggage : nothing
    return mp_b
end


############################## Fuel #########################################

#assumes fuel is drawn equally from both tanks, no need to model them
#individually for now
Base.@kwdef struct Fuel <: Piston.AbstractFuelSupply
    m_full::Float64 = 114.4 #maximum fuel mass (42 gal * 6 lb/gal * 0.454 kg/lb)
    m_empty::Float64 = 0.0 #residual fuel mass
end

Base.@kwdef struct FuelY
    m::Float64 = 0.0 #current fuel mass
end

#normalized fuel content (0: residual, 1: full)
init(::SystemX, ::Fuel) = [0.5] #cannot be a scalar, need an AbstractVector{<:Real}
init(::SystemY, ::Fuel) = FuelY()

function f_ode!(sys::System{Fuel}, pwp::System{<:Piston.Thruster})

    @unpack m_full, m_empty = sys.params #no need for subsystems
    m = m_empty + sys.x[1] * (m_full - m_empty) #current mass
    sys.ẋ .= -pwp.y.engine.ṁ / (m_full - m_empty)
    sys.y = FuelY(m)

end

fuel_available(sys::System{<:Fuel}) = (sys.y.m > 0)

function get_mp_b(fuel::System{Fuel})

    #in case x accidentally becomes negative (fuel is consumed beyond x=0 before
    #the engine dies)
    m_fuel = max(0.0, fuel.y.m)

    m_left = PointDistribution(0.5m_fuel)
    m_right = PointDistribution(0.5m_fuel)

    #fuel tanks reference frames
    frame_left = FrameTransform(r = SVector{3}(0.325, -2.845, 0))
    frame_right = FrameTransform(r = SVector{3}(0.325, 2.845, 0))

    mp_b = MassProperties()
    mp_b += MassProperties(m_left, frame_left)
    mp_b += MassProperties(m_right, frame_right)

    return mp_b
end

################################ Powerplant ####################################

Pwp() = Piston.Thruster(propeller = Propeller(t_bp = FrameTransform(r = [2.055, 0, 0.833])))

################################ Airframe ######################################

#P is introduced as a type parameter, because Piston.Thruster is itself a
#parametric type, and therefore not concrete
Base.@kwdef struct Airframe{P} <: AbstractAirframe
    str::Structure = Structure()
    aero::Aero = Aero()
    ldg::Ldg = Ldg()
    fuel::Fuel = Fuel()
    pld::Payload = Payload()
    pwp::P = Pwp()
end

################################################################################
####################### Update functions #######################################
################################################################################

################################################################################
######################## Avionics Update Functions ###########################

function f_ode!(avionics::System{BasicAvionics}, ::System{<:Airframe},
                ::KinematicData, ::AirflowData, ::System{<:AbstractTerrain})

    #here, avionics do nothing but update their output state. for a more complex
    #aircraft a continuous state-space autopilot implementation could go here
    @unpack throttle, Δ_aileron, aileron, Δ_elevator, elevator, pedals, Δ_pedals,
     brake_left, brake_right, flaps, mixture, eng_start, eng_stop = avionics.u

    return BasicAvionicsY(;
            throttle, Δ_aileron, aileron, Δ_elevator, elevator, pedals, Δ_pedals,
            brake_left, brake_right, flaps, mixture, eng_start, eng_stop)

end

#no digital components or state machines in BasicAvionics
@inline f_step!(::System{BasicAvionics}, ::System{<:Airframe}, ::System{<:AbstractKinematics}) = false
@inline f_disc!(::System{BasicAvionics}, ::System{<:Airframe}, ::System{<:AbstractKinematics}, Δt) = false


################################################################################
####################### Airframe Update Functions ##############################

function f_ode!(airframe::System{<:Airframe}, avionics::System{BasicAvionics},
                kin::KinematicData, air::AirflowData, trn::System{<:AbstractTerrain})

    @unpack aero, pwp, ldg, fuel, pld = airframe

    assign_component_inputs!(airframe, avionics)
    f_ode!(aero, pwp, air, kin, trn)
    f_ode!(ldg, kin, trn) #update landing gear continuous state & outputs
    f_ode!(pwp, air, kin) #update powerplant continuous state & outputs
    f_ode!(fuel, pwp) #update fuel system

    Systems.update_y!(airframe)

end

function f_step!(airframe::System{<:Airframe}, ::System{BasicAvionics}, ::System{<:AbstractKinematics})
    @unpack aero, pwp, fuel, ldg, fuel, pld = airframe

    x_mod = false
    x_mod = x_mod || f_step!(aero)
    x_mod = x_mod || f_step!(ldg)
    x_mod = x_mod || f_step!(pwp, fuel)
    return x_mod

end

#get_mp_b, get_wr_b and get_hr_b fall back to the @generated methods, which then
#recurse on Airframe subsystems

function assign_component_inputs!(airframe::System{<:Airframe}, avionics::System{BasicAvionics})

    @unpack throttle, Δ_aileron, aileron, Δ_elevator, elevator, pedals, Δ_pedals,
    brake_left, brake_right, flaps, mixture, eng_start, eng_stop = avionics.u
    @unpack aero, pwp, ldg = airframe

    pwp.u.engine.start = eng_start
    pwp.u.engine.stop = eng_stop
    pwp.u.engine.thr = throttle
    pwp.u.engine.mix = mixture
    ldg.u.nose.steering[] = (pedals + Δ_pedals) #pedals↑ (right pedal forward) -> nose wheel steering right
    ldg.u.left.braking[] = brake_left
    ldg.u.right.braking[] = brake_right
    aero.u.e = (elevator + Δ_elevator) #elevator↑ (stick forward) -> e↑
    aero.u.a = (aileron + Δ_aileron) #aileron↑ (stick right) -> a↑
    aero.u.r = -(pedals + Δ_pedals) #pedals↑ (left pedal forward) -> r↓
    aero.u.f = flaps #flaps↑ -> δf↑

    return nothing
end



################################################################################
############################# Input Interfaces ################################


####################### XBoxController Input Interface ########################

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
pedal_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function exp_axis_curve(x::Ranged{T}, args...; kwargs...) where {T}
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

function assign!(u::BasicAvionicsU, joystick::Joystick{XBoxController}, ::DefaultInputMapping)

    u.Δ_aileron = get_axis_value(joystick, :right_analog_x) |> aileron_curve
    u.Δ_elevator = -get_axis_value(joystick, :right_analog_y) |> elevator_curve
    u.Δ_pedals = get_axis_value(joystick, :left_analog_x) |> pedal_curve
    u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

    u.aileron -= 0.01 * was_released(joystick, :dpad_left)
    u.aileron += 0.01 * was_released(joystick, :dpad_right)
    u.elevator -= 0.01 * was_released(joystick, :dpad_down)
    u.elevator += 0.01 * was_released(joystick, :dpad_up)

    u.throttle += 0.1 * was_released(joystick, :button_Y)
    u.throttle -= 0.1 * was_released(joystick, :button_A)

    u.flaps += 0.3333 * was_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

end


################################################################################
############################### Cessna172R #####################################

const Cessna172R{K, F, V} = AircraftBase{K, F, V} where {K, F <: Airframe, V <: BasicAvionics}
Cessna172R(kinematics = LTF()) = AircraftBase( kinematics, Airframe(), BasicAvionics())

include("tools/trim.jl")
include("tools/linear.jl")
using .Trim
using .Linear

end #module