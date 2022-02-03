module C182T

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack
using Unitful
using Interpolations
using HDF5

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

############################## Avionics #################################

struct Avionics <: SystemDescriptor end

Base.@kwdef mutable struct AvionicsU
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

init_u(::Type{Avionics}) = AvionicsU()
init_y(::Type{Avionics}) = AvionicsY(zeros(SVector{9})...)


################################################################################
############################ Type Definitions ##################################
################################################################################

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

    prop = SimpleProp(kF = 4e-3, J = 0.005)

    left = EThruster(frame = FrameTransform(r = [2.055, 0, 0.833]), propeller = prop, motor = ElectricMotor(α = CW))
    right = EThruster(frame = FrameTransform(r = [2.055, 0, 0.833]), propeller = prop, motor = ElectricMotor(α = CCW))

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

    mlg_damper = SimpleDamper(k_s = ustrip(u"N/m", 0.5*5400u"lbf/ft"),
                              k_d_ext = ustrip(u"N/(m/s)", 0.4*1600u"lbf/(ft/s)"),
                              k_d_cmp = ustrip(u"N/(m/s)", 0.4*1600u"lbf/(ft/s)"))
    nlg_damper = SimpleDamper(k_s = ustrip(u"N/m", 1800u"lbf/ft"),
                              k_d_ext = ustrip(u"N/(m/s)", 600u"lbf/(ft/s)"),
                              k_d_cmp = ustrip(u"N/(m/s)", 600u"lbf/(ft/s)"))

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

init_x(::Type{Fuel}) = ComponentVector(m_left = 50.0, m_right = 50.0) #fuel tank contents

MassTrait(::System{Fuel}) = HasMass()
WrenchTrait(::System{Fuel}) = HasNoWrench()
AngularMomentumTrait(::System{Fuel}) = HasNoAngularMomentum()

function get_mp_b(sys::System{Fuel})

    frame_left = FrameTransform(r = SVector{3}(0.325, -2.845, 0))
    frame_right = FrameTransform(r = SVector{3}(0.325, 2.845, 0))

    m_left = PointMass(sys.x.m_left)
    m_right = PointMass(sys.x.m_right)

    mp_b = MassProperties()
    mp_b += MassProperties(m_left, frame_left)
    mp_b += MassProperties(m_right, frame_right)

    return mp_b
end


##### Aerodynamics ####

#Cmdeltae es negativo. Por tanto, un deltae > 0 provoca pitch down. esto
#significa que deltae > 0 corresponde a stab trailing edge down. igual que en
#Beaver

#Cldeltaa > 0. es decir, un deltaa > 0 provoca roll right. Es decir, corresponde
#a left aileron trailing edge down. esto es al reves que en Beaver

#CYdeltar > 0. es decir, deltar > 0 provoca fuerza lateral a la derecha. esto
#significa trailing edge left. como en Beaver

#ojo: en JSBSim enchufan simplemente la deflexion del aleron izquierdo como
#delta_a en los coeficientes. pero esto tiene un problema: como el intervalo de
#delta_a es asimetrico, el avion va a tener distinto mando en roll a un lado y a
#otro. esto esta mal. deberia ser una suma o un promedio de deflexiones. y
#entonces el intervalo si seria simetrico. por tanto, voy a suponer que es un
#promedio (porque si no me saldria el doble de autoridad de la que tiene), y que
#el intervalo es simetrico

#aunque las tablas de interpolacion saturen en alpha y beta en sus extremos, hay
#contribuciones que son un producto de alpha o beta por una derivada de
#estabilidad sin mas, asi que conviene saturar alpha y beta aparte de esto

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
    e::Bounded{Float64, -1, 1} = 0.0 #elevator control input (+ pitch down)
    a::Bounded{Float64, -1, 1} = 0.0 #aileron control input (+ roll right)
    r::Bounded{Float64, -1, 1} = 0.0 #rudder control input (+ yaw left)
    f::Bounded{Float64, 0, 1} = 0.0 # flap control input (+ flap down)
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
    coeffs::AeroCoeffs = AeroCoeffs() #aerodynamic coefficients
    wr_b::Wrench = Wrench() #aerodynamic Wrench, airframe
end

init_x(::Type{Aero}) = ComponentVector(α_filt = 0.0, β_filt = 0.0) #filtered airflow angles
init_y(::Type{Aero}) = AeroY()
init_d(::Type{Aero}) = AeroD()
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
################################################################################

#aerodynamic dataset taken from JSBSim's c172p

function load_aero_data()

    fname = "src/c182t_aero.h5"
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

#scale a normalized control input to its corresponding surface deflection range
function linear_scaling(u::Bounded{T, UMin, UMax}, δ_range::NTuple{2,Real}) where {T, UMin, UMax}
    @assert UMin != UMax
    return δ_range[1] + (δ_range[2] - δ_range[1])/(UMax - UMin) * (T(u) - UMin)
end

function get_aero_coeffs(; α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf, α_dot_nd, β_dot_nd, Δh_nd, stall)

    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = aero_data

    # C_D.z |> println
    # C_D.ge(Δh_nd) |> println
    # C_D.α_δf(α,δf) |> println
    # C_D.δf(δf) |> println
    # C_D.δe(δe) |> println
    # C_D.β(β) |> println

    # @show α β p_nd q_nd r_nd δa δr δe δf α_dot_nd β_dot_nd Δh_nd stall

    AeroCoeffs(
        C_D = C_D.z + C_D.ge(Δh_nd) * (C_D.α_δf(α,δf) + C_D.δf(δf)) + C_D.δe(δe) + C_D.β(β),
        C_Y = C_Y.δr * δr + C_Y.δa * δa + C_Y.β_δf(β,δf) + C_Y.p(α,δf) * p_nd + C_Y.r(α,δf) * r_nd,
        C_L = C_L.ge(Δh_nd) * (C_L.α(α,stall) + C_L.δf(δf)) + C_L.δe * δe + C_L.q * q_nd + C_L.α_dot * α_dot_nd,
        C_l = C_l.δa * δa + C_l.δr * δr + C_l.β * β + C_l.p * p_nd + C_l.r(α,δf) * r_nd,
        C_m = C_m.z + C_m.δe * δe + C_m.δf(δf) + C_m.α * α + C_m.q * q_nd + C_m.α_dot * α_dot_nd,
        C_n = C_n.δr * δr + C_n.δa * δa + C_n.β * β + C_n.p * p_nd + C_n.r * r_nd,
        # C_n = C_n.δr * δr + C_n.δa * δa + C_n.p * p_nd + C_n.r * r_nd,
    )

end

# f_cont!(sys::System{Aero}, pwp::System{Pwp},
#     air::AirData, kinematics::KinData, terrain::AbstractTerrain) = nothing

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

    @unpack ẋ, x, u, d, params = sys
    @unpack α_filt, β_filt = x
    @unpack e, a, r, f = u
    @unpack S, b, c, δe_range, δa_range, δr_range, δf_range, α_bounds, β_bounds, α_stall, V_min, τ = params
    @unpack TAS, q, v_wOb_b = air
    @unpack ω_lb_b = kinematics.vel
    @unpack n_e, h_o = kinematics.pos

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
    h_trn_Oa = get_terrain_data(terrain, l2d_Oa).altitude #orthometric
    Δh_nd = (h_Oa - h_trn_Oa) / b

    #stall state
    stall = Float64(d.stall)

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
        e, a, r, f, coeffs, wr_b)

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

######################## Avionics Update Functions ###########################

function f_cont!(avs::System{Avionics}, ::System{Airframe},
                ::KinData, ::AirData, ::AbstractTerrain)

    #here, avionics do nothing but update their output state. for a more complex
    #aircraft a continuous state-space autopilot implementation could go here
    @unpack throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
            pedals, brake_left, brake_right, flaps = avs.u

    return AvionicsY(; throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
                            pedals, brake_left, brake_right, flaps)

end

f_disc!(::System{Avionics}, ::System{Airframe}) = false

################################################################################
####################### Airframe Update Functions ##############################

#here we could define f_cont! to account for fuel consumption
f_cont!(::System{Fuel}, ::System{Pwp}) = nothing


#we can't fall back  on the default System Group implementation, because of the
#interactions between the different subsystems

function f_cont!(afm::System{Airframe}, avs::System{Avionics},
                kin::KinData, air::AirData, trn::AbstractTerrain)

    @unpack aero, pwp, ldg, fuel, pld = afm.subsystems

    assign_component_inputs!(afm, avs)
    f_cont!(ldg, kin, trn) #update landing gear continuous state & outputs
    f_cont!(pwp, kin, air) #update powerplant continuous state & outputs
    f_cont!(fuel, pwp) #update fuel system
    f_cont!(aero, pwp, air, kin, trn)

    # afm.y = (aero = aero.y, pwp = pwp.y, ldg = ldg.y, fuel = fuel.y, pld = pld.y )
    afm.y = (aero = aero.y, pwp = pwp.y, ldg = ldg.y)

end

function assign_component_inputs!(afm::System{<:Airframe}, avs::System{<:Avionics})

    @unpack throttle, yoke_Δx, yoke_x0, yoke_Δy, yoke_y0,
            pedals, brake_left, brake_right, flaps = avs.u
    @unpack aero, pwp, ldg = afm.subsystems

    #yoke_Δx is the offset with respect to the force-free position yoke_x0
    #yoke_Δy is the offset with respect to the force-free position yoke_y0

    pwp.u.left.throttle = throttle
    pwp.u.right.throttle = throttle
    ldg.u.nose.steering[] = pedals
    ldg.u.left.braking[] = brake_left
    ldg.u.right.braking[] = brake_right
    aero.u.e = -(yoke_y0 + yoke_Δy) #+yoke_Δy and +yoke_y0 are back and +δe is pitch down, need to invert it
    aero.u.a = (yoke_x0 + yoke_Δx) #+yoke_Δx and +yoke_x0 are right and +δa is roll right
    aero.u.r = -pedals # +pedals is right and +δr is yaw left
    aero.u.f = flaps # +flaps is flaps down and +δf is flaps down

    return nothing
end

function f_disc!(afm::System{Airframe}, ::System{Avionics})
    #fall back to the default SystemGroup implementation
    return f_disc!(afm)
end

f_disc!(::System{OEW}) = false
f_disc!(::System{Fuel}) = false
f_disc!(::System{Payload}) = false
#Pwp and Ldg are SystemGroups, so we can rely on the fallback f_disc!

#for Fuel we could define discrete states to select which tank(s) we're drawing from

#get_mp_b, get_wr_b and get_hr_b use the fallback for SystemGroups, which in turn call
#get_mp_b, get_wr_b and get_hr_b on aero, pwp and ldg


################################################################################
############################# Input Interfaces ################################


####################### XBoxController Input Interface ########################

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
pedal_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

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

    u = ac.u.avionics

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

function C182TDescriptor(; id = ID(), kin = KinLTF(), afm = Airframe(), avs = Avionics())
    AircraftBase( id; kin, afm, avs)
end


    # splt_α = thplot(t, rad2deg.(α_b);
    #     title = "Angle of Attack", ylabel = L"$\alpha \ (deg)$",
    #     label = "", kwargs...)

    # splt_β = thplot(t, rad2deg.(β_b);
    #     title = "Angle of Sideslip", ylabel = L"$\beta \ (deg)$",
    #     label = "", kwargs...)

    # pd["05_α_β"] = plot(splt_α, splt_β;
    #     plot_title = "Airflow Angles [Airframe]",
    #     layout = (1,2),
    #     kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

end #module