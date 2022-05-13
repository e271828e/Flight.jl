module C172R

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack
using Unitful
using Interpolations
using HDF5

using Flight.Systems
using Flight.Utils
using Flight.Attitude
using Flight.Terrain
using Flight.Air
using Flight.Kinematics
using Flight.Dynamics
using Flight.LandingGear
using Flight.Propellers
using Flight.Piston
using Flight.Aircraft: AircraftBase, AbstractAircraftID, AbstractVehicle, AbstractAerodynamics, AbstractAvionics
using Flight.Input: XBoxController, get_axis_value, is_released

import Flight.Systems: init, f_cont!, f_disc!
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_wr_b, get_mp_b
import Flight.Input: assign!
import Flight.Piston: fuel_available

export C172RAircraft

struct ID <: AbstractAircraftID end

############################## Avionics #################################

struct Avionics <: AbstractAvionics end

Base.@kwdef mutable struct AvionicsU
    throttle::Ranged{Float64, 0, 1} = 0.0
    yoke_Δx::Ranged{Float64, -1, 1} = 0.0 #ailerons (+ bank right)
    yoke_x0::Ranged{Float64, -1, 1} = 0.0 #ailerons (+ bank right)
    yoke_Δy::Ranged{Float64, -1, 1} = 0.0 #elevator (+ pitch up)
    yoke_y0::Ranged{Float64, -1, 1} = 0.0 #elevator (+ pitch up)
    pedals::Ranged{Float64, -1, 1} = 0.0 #rudder and nose wheel (+ yaw right)
    brake_left::Ranged{Float64, 0, 1} = 0.0 #[0, 1]
    brake_right::Ranged{Float64, 0, 1} = 0.0 #[0, 1]
    flaps::Ranged{Float64, 0, 1} = 0.0 #[0, 1]
end

Base.@kwdef struct AvionicsY
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

init(::Avionics, ::SystemU) = AvionicsU()
init(::Avionics, ::SystemY) = AvionicsY(zeros(SVector{9})...)


################################################################################
############################ Vehicle Subsystems ################################

################################ Airframe ######################################

struct Airframe <: SystemDescriptor end

#airframe mass properties. these are assumed to include any components attached
#rigidly to the airframe itself, such as power plant and landing gear, but not
#fuel contents or payload
const mp_b_afm = let
    #define the airframe as a RigidBody
    afm_G = RigidBody(767.0, SA[820.0 0 0; 0 1164.0 0; 0 0 1702.0])
    #define the transform from the origin of the vehicle reference frame (Ob)
    #to the airframe's center of mass (G)
    t_Ob_G = FrameTransform(r = SVector{3}(0.056, 0, 0.582))
    #compute the airframe's mass properties at Ob
    MassProperties(afm_G, t_Ob_G)
end

MassTrait(::System{Airframe}) = HasMass()

#the airframe itself receives no external actions. these are considered to act
#upon the vehicle's aerodynamics, power plant and landing gear. the same goes
#for rotational angular momentum.
WrenchTrait(::System{Airframe}) = GetsNoExternalWrench()
AngularMomentumTrait(::System{Airframe}) = HasNoAngularMomentum()

get_mp_b(::System{Airframe}) = mp_b_afm


############################ Aerodynamics ######################################

Base.@kwdef struct Aero <: AbstractAerodynamics
    S::Float64 = 16.165 #wing area
    b::Float64 = 10.912 #wingspan
    c::Float64 = 1.494 #mean aerodynamic chord
    δe_range::NTuple{2,Float64} = deg2rad.((-28, 23)) #elevator deflection range (rad)
    δa_range::NTuple{2,Float64} = deg2rad.((-20, 20)) #aileron deflection range (rad)
    δr_range::NTuple{2,Float64} = deg2rad.((-16, 16)) #rudder deflection range (rad)
    δf_range::NTuple{2,Float64} = deg2rad.((0, 30)) #flap deflection range (rad)
    α_bounds::NTuple{2,Float64} = (-0.1, 0.4) #α bounds for aerodynamic dataset input
    β_bounds::NTuple{2,Float64} = (-0.1, 0.4) #β bounds for aerodynamic dataset input
    α_stall::NTuple{2,Float64} = (0.09, 0.36) #α values for stall hysteresis switching
    V_min::Float64 = 1.0 #lower airspeed threshold for non-dimensional angle rates
    τ::Float64 = 0.1 #time constant for airflow angle filtering
end

Base.@kwdef mutable struct AeroU
    e::Ranged{Float64, -1, 1} = 0.0 #elevator control input (+ pitch down)
    a::Ranged{Float64, -1, 1} = 0.0 #aileron control input (+ roll right)
    r::Ranged{Float64, -1, 1} = 0.0 #rudder control input (+ yaw left)
    f::Ranged{Float64, 0, 1} = 0.0 # flap control input (+ flap down)
end

Base.@kwdef mutable struct AeroD #discrete state
    stall::Bool = false
end

Base.@kwdef struct AeroCoeffs
    C_D::Float64 = 0.0
    C_Y::Float64 = 0.0
    C_L::Float64 = 0.0
    C_l::Float64 = 0.0
    C_m::Float64 = 0.0
    C_n::Float64 = 0.0
end

Base.@kwdef struct AeroY
    e::Float64 = 0.0 #normalized elevator control input
    a::Float64 = 0.0 #normalized aileron control input
    r::Float64 = 0.0 #normalized rudder control input
    f::Float64 = 0.0 #normalized flap control input
    α::Float64 = 0.0 #clamped AoA
    β::Float64 = 0.0 #clamped AoS
    α_filt::Float64 = 0.0 #filtered AoA
    β_filt::Float64 = 0.0 #filtered AoS
    α_filt_dot::Float64 = 0.0 #filtered AoA derivative
    β_filt_dot::Float64 = 0.0 #filtered AoS derivative
    stall::Bool = false #stall state
    coeffs::AeroCoeffs = AeroCoeffs() #aerodynamic coefficients
    wr_b::Wrench = Wrench() #aerodynamic Wrench, vehicle frame
end

init(::Aero, ::SystemX) = ComponentVector(α_filt = 0.0, β_filt = 0.0) #filtered airflow angles
init(::Aero, ::SystemY) = AeroY()
init(::Aero, ::SystemU) = AeroU()
init(::Aero, ::SystemD) = AeroD()

#dataset taken from JSBSim's c172p
function load_aero_data()

    fname = "src/aircraft/c172r/c172r.h5"
    fid = h5open(fname, "r")

    gr_C_D = fid["C_D"]
    gr_C_Y = fid["C_Y"]
    gr_C_L = fid["C_L"]
    gr_C_l = fid["C_l"]
    gr_C_m = fid["C_m"]
    gr_C_n = fid["C_n"]

    C_D = (
        z = gr_C_D["zero"] |> read,
        β = LinearInterpolation(gr_C_D["β"]["β"] |> read, gr_C_D["β"]["data"] |> read, extrapolation_bc = Flat()),
        δe = LinearInterpolation(gr_C_D["δe"]["δe"] |> read, gr_C_D["δe"]["data"] |> read, extrapolation_bc = Flat()),
        δf = LinearInterpolation(gr_C_D["δf"]["δf"] |> read, gr_C_D["δf"]["data"] |> read, extrapolation_bc = Flat()),
        α_δf = LinearInterpolation((gr_C_D["α_δf"]["α"] |> read,  gr_C_D["α_δf"]["δf"] |> read), gr_C_D["α_δf"]["data"] |> read, extrapolation_bc = Flat()),
        ge = LinearInterpolation(gr_C_D["ge"]["Δh_nd"] |> read, gr_C_D["ge"]["data"] |> read, extrapolation_bc = Flat())
    )

    C_Y = (
        δr = gr_C_Y["δr"] |> read,
        δa = gr_C_Y["δa"] |> read,
        β_δf = LinearInterpolation((gr_C_Y["β_δf"]["β"] |> read,  gr_C_Y["β_δf"]["δf"] |> read), gr_C_Y["β_δf"]["data"] |> read, extrapolation_bc = Flat()),
        p = LinearInterpolation((gr_C_Y["p"]["α"] |> read,  gr_C_Y["p"]["δf"] |> read), gr_C_Y["p"]["data"] |> read, extrapolation_bc = Flat()),
        r = LinearInterpolation((gr_C_Y["r"]["α"] |> read,  gr_C_Y["r"]["δf"] |> read), gr_C_Y["r"]["data"] |> read, extrapolation_bc = Flat()),
    )

    C_L = (
        δe = gr_C_L["δe"] |> read,
        q = gr_C_L["q"] |> read,
        α_dot = gr_C_L["α_dot"] |> read,
        α = LinearInterpolation((gr_C_L["α"]["α"] |> read,  gr_C_L["α"]["stall"] |> read), gr_C_L["α"]["data"] |> read, extrapolation_bc = Flat()),
        δf = LinearInterpolation(gr_C_L["δf"]["δf"] |> read, gr_C_L["δf"]["data"] |> read, extrapolation_bc = Flat()),
        ge = LinearInterpolation(gr_C_L["ge"]["Δh_nd"] |> read, gr_C_L["ge"]["data"] |> read, extrapolation_bc = Flat())
    )

    C_l = (
        δa = gr_C_l["δa"] |> read,
        δr = gr_C_l["δr"] |> read,
        β = gr_C_l["β"] |> read,
        p = gr_C_l["p"] |> read,
        r = LinearInterpolation((gr_C_l["r"]["α"] |> read,  gr_C_l["r"]["δf"] |> read), gr_C_l["r"]["data"] |> read, extrapolation_bc = Flat()),
    )

    C_m = (
        z = gr_C_m["zero"] |> read,
        δe = gr_C_m["δe"] |> read,
        α = gr_C_m["α"] |> read,
        q = gr_C_m["q"] |> read,
        α_dot = gr_C_m["α_dot"] |> read,
        δf = LinearInterpolation(gr_C_m["δf"]["δf"] |> read, gr_C_m["δf"]["data"] |> read, extrapolation_bc = Flat()),
    )

    C_n = (
        δr = gr_C_n["δr"] |> read,
        δa = gr_C_n["δa"] |> read,
        β = gr_C_n["β"] |> read,
        p = gr_C_n["p"] |> read,
        r = gr_C_n["r"] |> read,
    )

    close(fid)

    return (C_D = C_D, C_Y = C_Y, C_L = C_L, C_l = C_l, C_m = C_m, C_n = C_n)

end

const aero_data = load_aero_data()

function get_aero_coeffs(; α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf, α_dot_nd, β_dot_nd, Δh_nd, stall)

    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = aero_data

    # @show α β p_nd q_nd r_nd δa δr δe δf α_dot_nd β_dot_nd Δh_nd stall

    AeroCoeffs(
        C_D = C_D.z + C_D.ge(Δh_nd) * (C_D.α_δf(α,δf) + C_D.δf(δf)) + C_D.δe(δe) + C_D.β(β),
        C_Y = C_Y.δr * δr + C_Y.δa * δa + C_Y.β_δf(β,δf) + C_Y.p(α,δf) * p_nd + C_Y.r(α,δf) * r_nd,
        C_L = C_L.ge(Δh_nd) * (C_L.α(α,stall) + C_L.δf(δf)) + C_L.δe * δe + C_L.q * q_nd + C_L.α_dot * α_dot_nd,
        C_l = C_l.δa * δa + C_l.δr * δr + C_l.β * β + C_l.p * p_nd + C_l.r(α,δf) * r_nd,
        C_m = C_m.z + C_m.δe * δe + C_m.δf(δf) + C_m.α * α + C_m.q * q_nd + C_m.α_dot * α_dot_nd,
        C_n = C_n.δr * δr + C_n.δa * δa + C_n.β * β + C_n.p * p_nd + C_n.r * r_nd,
    )

end

function f_cont!(sys::System{Aero}, ::System{<:Piston.Thruster},
    air::AirflowData, kinematics::KinData, terrain::AbstractTerrain)

    #for this aircraft we have chosen the aerodynamics frame fa as the main
    #aircraft reference frame fb. any airflow variables computed in frame fb are
    #directly usable here. if this wasn't the case, either because the origin Oa
    #is not coincident with Ob and we care about the velocity lever arm, or
    #because the axes εa and εb are different, we would need to compute v_wOa_a
    #from v_wOb_b and then recompute all the derived air data variables

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

    @unpack ẋ, x, u, d, params = sys
    @unpack α_filt, β_filt = x
    @unpack e, a, r, f = u
    @unpack S, b, c, δe_range, δa_range, δr_range, δf_range, α_bounds, β_bounds, α_stall, V_min, τ = params
    @unpack TAS, q, v_wOb_b = air
    @unpack ω_lb_b = kinematics.vel
    @unpack n_e, h_o = kinematics.pos
    stall = d.stall

    v_wOb_a = v_wOb_b
    α, β = get_airflow_angles(v_wOb_a)
    α = clamp(α, α_bounds[1], α_bounds[2])
    β = clamp(β, β_bounds[1], β_bounds[2])
    V = max(TAS, V_min)

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

    coeffs = get_aero_coeffs(; α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf,
                               α_dot_nd, β_dot_nd, Δh_nd, stall)

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

function f_disc!(sys::System{Aero})
    #stall hysteresis
    α = sys.y.α
    α_stall = sys.params.α_stall
    if α > α_stall[2]
        sys.d.stall = true
    elseif α < α_stall[1]
        sys.d.stall = false
    end
    return false
end


############################# Landing Gear ####################################

struct Ldg <: SystemGroupDescriptor
    left::LandingGearUnit{NoSteering, DirectBraking, Strut{SimpleDamper}}
    right::LandingGearUnit{NoSteering, DirectBraking, Strut{SimpleDamper}}
    nose::LandingGearUnit{DirectSteering, NoBraking, Strut{SimpleDamper}}
end

MassTrait(::System{Ldg}) = HasNoMass()
WrenchTrait(::System{Ldg}) = GetsExternalWrench()
AngularMomentumTrait(::System{Ldg}) = HasNoAngularMomentum()

function Ldg()

    mlg_damper = SimpleDamper(k_s = ustrip(u"N/m", 0.5*5400u"lbf/ft"),
                              k_d_ext = ustrip(u"N/(m/s)", 0.4*1600u"lbf/(ft/s)"),
                              k_d_cmp = ustrip(u"N/(m/s)", 0.4*1600u"lbf/(ft/s)"))
    nlg_damper = SimpleDamper(k_s = ustrip(u"N/m", 1800u"lbf/ft"),
                              k_d_ext = ustrip(u"N/(m/s)", 0.4*600u"lbf/(ft/s)"),
                              k_d_cmp = ustrip(u"N/(m/s)", 0.4*600u"lbf/(ft/s)"))

    left = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-0.381, -1.092, 1.902], q = RQuat() ),
            l_OsP = 0.0,
            damper = mlg_damper),
        braking = DirectBraking())

    right = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-0.381, 1.092, 1.902], q = RQuat() ),
            l_OsP = 0.0,
            damper = mlg_damper),
        braking = DirectBraking())

    nose = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [1.27, 0, 1.9] , q = RQuat()),
            l_OsP = 0.0,
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

struct Payload <: SystemDescriptor
    pilot::MassProperties #mp_b
    copilot::MassProperties #mp_b
    psg_left::MassProperties #mp_b
    psg_right::MassProperties #mp_b
    baggage::MassProperties #mp_b

    function Payload( ; pilot::AbstractMassDistribution = PointMass(75),
                        copilot::AbstractMassDistribution = PointMass(75),
                        psg_left::AbstractMassDistribution = PointMass(75),
                        psg_right::AbstractMassDistribution = PointMass(75),
                        baggage::AbstractMassDistribution = PointMass(50))

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

init(::Payload, ::SystemU) = PayloadU()

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
init(::Fuel, ::SystemX) = [0.5] #cannot be a scalar, need an AbstractVector{<:Real}
init(::Fuel, ::SystemY) = FuelY()

function f_cont!(sys::System{Fuel}, pwp::System{<:Piston.Thruster})

    @unpack m_full, m_empty = sys.params #no need for subsystems
    m = m_empty + sys.x[1] * (m_full - m_empty) #current mass
    sys.ẋ .= -pwp.y.engine.ṁ / (m_full - m_empty)
    sys.y = FuelY(m)

end

f_disc!(::System{Fuel}) = false

fuel_available(sys::System{<:Fuel}) = (sys.y.m > 0)

function get_mp_b(fuel::System{Fuel})

    #in case x accidentally becomes negative (fuel is consumed beyond x=0 before
    #the engine dies)
    m_fuel = max(0.0, fuel.y.m)

    m_left = PointMass(0.5m_fuel)
    m_right = PointMass(0.5m_fuel)

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

################################ Vehicle ######################################

#P is introduced as a type parameter, because Piston.Thruster is itself a
#parametric type, and therefore not concrete
Base.@kwdef struct Vehicle{P} <: AbstractVehicle
    afm::Airframe = Airframe()
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

function f_cont!(avionics::System{Avionics}, ::System{<:Vehicle},
                ::KinData, ::AirflowData, ::AbstractTerrain)

    #here, avionics do nothing but update their output state. for a more complex
    #aircraft a continuous state-space autopilot implementation could go here
    @unpack throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
            pedals, brake_left, brake_right, flaps = avionics.u

    return AvionicsY(; throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
                            pedals, brake_left, brake_right, flaps)

end

f_disc!(::System{Avionics}, ::System{<:Vehicle}) = false


################################################################################
####################### Vehicle Update Functions ##############################

function f_cont!(vehicle::System{<:Vehicle}, avionics::System{Avionics},
                kin::KinData, air::AirflowData, trn::AbstractTerrain)

    @unpack aero, pwp, ldg, fuel, pld = vehicle

    assign_component_inputs!(vehicle, avionics)
    f_cont!(aero, pwp, air, kin, trn)
    f_cont!(ldg, kin, trn) #update landing gear continuous state & outputs
    f_cont!(pwp, air, kin) #update powerplant continuous state & outputs
    f_cont!(fuel, pwp) #update fuel system

    Systems.assemble_y!(vehicle)

end

function f_disc!(vehicle::System{<:Vehicle}, ::System{Avionics})
    @unpack afm, aero, pwp, fuel, ldg, fuel, pld = vehicle

    x_mod = false
    x_mod = x_mod || f_disc!(aero)
    x_mod = x_mod || f_disc!(ldg)
    x_mod = x_mod || f_disc!(pwp, fuel)
    return x_mod

end

#get_mp_b, get_wr_b and get_hr_b fall back to the @generated methods, which then
#recurse on Vehicle subsystems

function assign_component_inputs!(vehicle::System{<:Vehicle}, avionics::System{<:Avionics})

    @unpack throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
            pedals, brake_left, brake_right, flaps = avionics.u
    @unpack aero, pwp, ldg = vehicle

    #yoke_Δx is the offset with respect to the force-free position yoke_x0
    #yoke_Δy is the offset with respect to the force-free position yoke_y0

    pwp.u.engine.thr = throttle
    ldg.u.nose.steering[] = pedals
    ldg.u.left.braking[] = brake_left
    ldg.u.right.braking[] = brake_right
    aero.u.e = -(yoke_y0 + yoke_Δy) #+yoke_Δy and +yoke_y0 are back and +δe is pitch down, need to invert it
    aero.u.a = (yoke_x0 + yoke_Δx) #+yoke_Δx and +yoke_x0 are right and +δa is roll right
    aero.u.r = -pedals # +pedals is right and +δr is yaw left
    aero.u.f = flaps # +flaps is flaps down and +δf is flaps down

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

function assign!(u::AvionicsU, joystick::XBoxController)

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

    u.flaps += 0.3333 * is_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * is_released(joystick, :left_bumper)

end

#Aircraft constructor override keyword inputs to customize

function C172RAircraft(; id = ID(), kinematics = KinLTF(), vehicle = Vehicle(), avionics = Avionics())
    AircraftBase( id; kinematics, vehicle, avionics)
end


end #module