module LandingGear

using LinearAlgebra, StaticArrays

using FlightCore
using FlightCore.GUI
using ..Types
using ..Attitude
using ..Geodesy
using ..Kinematics
using ..Dynamics
using ..Terrain
using ..Control: PIVector, PIVectorY

export LandingGearUnit, Strut, Contact, SimpleDamper
export NoSteering, NoBraking, DirectSteering, DirectBraking

#basis elements, for convenience
const e1 = SVector{3,Float64}(1,0,0)
const e3 = SVector{3,Float64}(0,0,1)


################################################################################
################################# Steering #####################################

abstract type AbstractSteering <: ModelDefinition end

################################ NoSteering ####################################

struct NoSteering <: AbstractSteering end

struct NoSteeringY end

Modeling.Y(::NoSteering) = NoSteeringY()

@no_updates NoSteering

get_steering_angle(::Model{NoSteering}, args...) = 0.0


############################### DirectSteering #################################

@kwdef struct DirectSteering <: AbstractSteering
    ψ_max::Float64 = π/6
end

@kwdef mutable struct DirectSteeringU
    engaged::Bool = true
    input::Ranged{Float64, -1., 1.} = Ranged(0.0, -1., 1.)
end

@kwdef struct DirectSteeringY
    engaged::Bool = true
    input::Float64 = 0.0
end

Modeling.U(::DirectSteering) = DirectSteeringU()
Modeling.Y(::DirectSteering) = DirectSteeringY()

@no_periodic DirectSteering
@no_step DirectSteering

function Modeling.f_ode!(mdl::Model{DirectSteering})
    (; engaged, input) = mdl.u
    mdl.y = DirectSteeringY(; engaged, input)
end

function get_steering_angle(mdl::Model{DirectSteering}, ψ_v::Real)
    (; engaged, input) = mdl.y
    ψ_sw = engaged ? Float64(input) * mdl.ψ_max : ψ_v
    return ψ_sw
end

function GUI.draw(mdl::Model{DirectSteering})

    (; engaged, input) = mdl.y
    TextFormatted("Engaged: $engaged")
    TextFormatted("Steering Input: $(Float64(input))")
end


################################################################################
################################# Braking ######################################

abstract type AbstractBraking <: ModelDefinition end

############################### NoBraking ######################################

struct NoBraking <: AbstractBraking end

struct NoBrakingY end

Modeling.Y(::NoBraking) = NoBrakingY()

@no_updates NoBraking

get_braking_factor(::Model{NoBraking}) = 0.0


############################# DirectBraking ####################################

@kwdef struct DirectBraking <: AbstractBraking
   η_br::Float64 = 1.0 #braking efficiency
end

@kwdef struct DirectBrakingY
   κ_br::Float64 = 0.0 #braking coefficient
end

Modeling.U(::DirectBraking) = Ref(Ranged(0.0, 0., 1.))
Modeling.Y(::DirectBraking) = DirectBrakingY()

@no_periodic DirectBraking
@no_step DirectBraking

function Modeling.f_ode!(mdl::Model{DirectBraking})
    mdl.y = DirectBrakingY(Float64(mdl.u[]) * mdl.η_br)
end

get_braking_factor(mdl::Model{DirectBraking}) = mdl.y.κ_br

function GUI.draw(mdl::Model{DirectBraking})

    TextFormatted("Braking Input: $(Float64(mdl.u[]))")
    TextFormatted("Braking Coefficient: $(mdl.y.κ_br)")

end

################################################################################
################################## Strut #######################################

################################### Damper #####################################

abstract type AbstractDamper end #not a Model!

@kwdef struct SimpleDamper <: AbstractDamper
    k_s::Float64 = 25000 #spring constant
    k_d_ext::Float64 = 1000 #extension damping coefficient
    k_d_cmp::Float64 = 1000 #compression damping coefficient
    F_max::Float64 = 50000 #maximum allowable damper force
end

#Force exerted by the damper along zs. The elongation ξ is positive along z_s.
# When the damper force is positive (weight on wheel), it pushes the piston rod
#assembly downwards along the positive z_s direction. When it is negative, it
#pulls the piston rod upwards along the negative z_s axis
function get_force(c::SimpleDamper, ξ::Real, ξ_dot::Real)
    k_d = (ξ_dot > 0 ? c.k_d_ext : c.k_d_cmp)
    F = -(c.k_s * ξ + k_d * ξ_dot)
    return F
end

########################### FrictionCoefficients ###############################

struct FrictionCoefficients
    μ_s::Float64 #static friction coefficient
    μ_d::Float64 #dynamic friction coefficient (μ_d < μ_s)
    v_s::Float64 #static friction upper velocity threshold
    v_d::Float64 #dynamic friction lower velocity threshold (v_d > v_s)

    function FrictionCoefficients(; μ_s = 0, μ_d = 0, v_s = 0, v_d = 1)
        @assert μ_s >= μ_d
        @assert v_s < v_d
        new(μ_s, μ_d, v_s, v_d)
    end
end

function get_μ(fr::FrictionCoefficients, v::Real)
    (; v_s, v_d, μ_s, μ_d) = fr
    κ_sd = clamp((abs(v) - v_s) / (v_d - v_s), 0, 1)
    return κ_sd * μ_d + (1 - κ_sd) * μ_s
end

abstract type RollingOrSkidding end
struct Rolling <: RollingOrSkidding end
struct Skidding <: RollingOrSkidding end

FrictionCoefficients(::Rolling, ::SurfaceType) =
    FrictionCoefficients(μ_s = 0.03, μ_d = 0.02, v_s = 0.005, v_d = 0.01)

function FrictionCoefficients(::Skidding, srf::SurfaceType)
    if srf == DryTarmac
        FrictionCoefficients( μ_s = 0.75, μ_d = 0.25, v_s = 0.005, v_d = 0.01)
    elseif srf == WetTarmac
        FrictionCoefficients(μ_s = 0.25, μ_d = 0.15, v_s = 0.005, v_d = 0.01)
    elseif srf == IcyTarmac
        FrictionCoefficients(μ_s = 0.075, μ_d = 0.025, v_s = 0.005, v_d = 0.01)
    else
        error("Unrecognized surface type")
    end
end


################################## Strut #######################################

@kwdef struct GroundCrash <: SimulationTermination
    msg::String = ""
end

Base.showerror(io::IO, e::GroundCrash) = print(io, "GroundCrash: ", e.msg)

@kwdef struct Strut{D<:AbstractDamper} <: ModelDefinition
    t_bs::FrameTransform = FrameTransform() #vehicle to strut frame transform
    l_0::Float64 = 0.0 #strut natural length from airframe attachment point to wheel endpoint
    damper::D = SimpleDamper()
    α_cs_max::Float64 = deg2rad(60) #maximum contact normal to strut angle
    ξ_dot_max::Float64 = 10.0 #maximum damper compression rate
end

@kwdef struct StrutY #defaults should be consistent with wow = 0
    wow::Bool = false #weight-on-wheel flag
    ξ::Float64 = 0.0 #damper elongation
    ξ_dot::Float64 = 0.0 #damper elongation rate
    F_dmp_zs::Float64 = 0.0 #axial damper force
    ψ_sw::Float64 = 0.0 #steering angle
    α_cs::Float64 = 0.0 #angle from contact frame z-axis (terrain normal) to strut axis
    t_sc::FrameTransform = FrameTransform() #strut to contact frame transform
    t_bc::FrameTransform = FrameTransform() #body to contact frame transform
    v_ec_xy::SVector{2,Float64} = zeros(SVector{2}) #contact point velocity
    surface::SurfaceType = DryTarmac #terrain surface type at the contact point
end

Modeling.Y(::Strut) = StrutY()

@no_periodic Strut

function Modeling.f_ode!(mdl::Model{<:Strut},
                        steering::Model{<:AbstractSteering},
                        terrain::Model{<:AbstractTerrain},
                        kin_data::KinData)

    (; t_bs, l_0, damper) = mdl.parameters
    (; q_eb, r_eb_e, v_eb_b, ω_eb_b) = kin_data

    q_bs = t_bs.q #body frame to strut frame rotation
    r_bs_b = t_bs.r #strut frame origin

    #strut line in ECEF coordinates
    q_es = q_eb ∘ q_bs
    ks_s = e3
    ks_e = q_es(ks_s) #strut axis
    r_bs_e = q_eb(r_bs_b)
    r_es_e = r_eb_e + r_bs_e #strut frame origin

    #do we have contact? cast a ray from the strut frame origin along the strut
    #axis, bounded by the strut's natural length
    hit = SurfaceIntersection(terrain, Geocentric(r_es_e), ks_e, l_0)
    wow = hit.valid

    if !wow #no contact
        mdl.y = StrutY(; wow) #everything else set to default
        return
    end

    (; P, kt_P_e, surface) = hit.data
    r_ec_e = Geocentric(P)[:] #contact frame origin
    kc_e = kt_P_e #contact frame z-axis (inward terrain normal)

    r_sc_e = r_ec_e - r_es_e
    l = r_sc_e ⋅ ks_e #distance from strut frame origin to contact point
    ξ = l - l_0 #non-positive, since l ≤ l_max = l_0

    cos_α = kc_e ⋅ ks_e #cosine of angle between contact z-axis and strut axis
    α_cs = acos(clamp(cos_α, -1, 1))

    r_sc_s = ks_s * l #contact frame position with respect to strut frame
    r_sc_b = q_bs(r_sc_s)
    r_bc_b = r_sc_b + r_bs_b #contact frame position with respect to body frame

    #contact frame origin velocity due to rigid body motion
    v_ec_b_body = v_eb_b + ω_eb_b × r_bc_b #body frame
    q_sb = q_bs'
    v_ec_s_body = q_sb(v_ec_b_body) #strut frame
    ψ_v = atan(v_ec_s_body[2], v_ec_s_body[1]) #azimuth

    #wheel frame axes
    ψ_sw = get_steering_angle(steering, ψ_v)
    q_sw = Rz(ψ_sw) #strut to wheel axes rotation
    q_ew = q_es ∘ q_sw #ECEF to wheel axes rotation

    #contact frame axes, ECEF coordinates
    iw_w = e1
    iw_e = q_ew(iw_w) #wheel x-axis
    iw_e_trn = iw_e - (iw_e ⋅ kc_e) * kc_e #projection of iw_e onto the terrain tangent plane
    ic_e = normalize(iw_e_trn) #contact x-axis
    jc_e = kc_e × ic_e #contact y-axis
    q_ec = RMatrix(SMatrix{3,3}([ic_e jc_e kc_e]), normalization = false)
    q_se = q_es'
    q_sc = q_se ∘ q_ec
    q_bc = q_bs ∘ q_sc

    #construct contact frame transforms
    t_sc = FrameTransform(r_sc_s, q_sc)
    t_bc = FrameTransform(r_bc_b, q_bc)

    #contact frame origin velocity due to rigid body motion, contact frame
    q_cb = q_bc'
    v_ec_c_body = q_cb(v_ec_b_body)

    #compute the damper elongation rate required to cancel the rigid body
    #contribution to the contact point velocity along the contact frame z axis.
    #the strut axis projected on the contact frame z-axis is ks_c[3] = cos_α
    ξ_dot = -v_ec_c_body[3] / cos_α

    #force exerted by the damper along the strut frame's z axis
    F_dmp_zs = get_force(damper, ξ, ξ_dot)

    #total contact point velocity, contact frame. its z-component is cancelled
    #out by the damper elongation rate, so only the in-plane components remain
    q_cs = q_sc'
    ks_c = q_cs(ks_s)
    v_ec_dmp_c = ks_c * ξ_dot #contact point velocity due to elongation rate
    v_ec_c = v_ec_c_body + v_ec_dmp_c
    v_ec_xy = v_ec_c[SVector(1,2)]

    mdl.y = StrutY(; wow, ξ, ξ_dot, F_dmp_zs, ψ_sw, α_cs, t_sc, t_bc, v_ec_xy, surface)

end

#sanity checks for crash detection
function Modeling.f_step!(mdl::Model{<:Strut})

    (; wow, α_cs, ξ_dot, F_dmp_zs) = mdl.y
    (; α_cs_max, ξ_dot_max, damper) = mdl.parameters

    wow || return nothing

    (α_cs > α_cs_max) && throw(GroundCrash(
        "Contact normal to strut angle α_cs = $(rad2deg(α_cs)) deg " *
        "at t = $(mdl.t[]) s"))

    (-ξ_dot > ξ_dot_max) && throw(GroundCrash(
        "Damper compression rate ξ_dot = $(-ξ_dot) m/s " *
        "at t = $(mdl.t[]) s"))

    (abs(F_dmp_zs) > damper.F_max) && throw(GroundCrash(
        "Damper force F_dmp_zs = $F_dmp_zs N " *
        "at t = $(mdl.t[]) s"))

    return nothing

end



#################################### GUI #######################################

function GUI.draw(mdl::Model{<:Strut}, window_label::String = "Strut")

    (; wow, ξ, ξ_dot, F_dmp_zs, ψ_sw, α_cs, v_ec_xy, surface) = mdl.y

        TextFormatted("Weight on Wheel: $wow")
        TextFormatted(@sprintf("Damper Elongation: %.7f m", ξ))
        TextFormatted(@sprintf("Damper Elongation Rate: %.7f m/s", ξ_dot))
        TextFormatted(@sprintf("Axial Damper Force: %.7f N", F_dmp_zs))
        TextFormatted(@sprintf("Wheel Steering Angle: %.7f deg", rad2deg(ψ_sw)))
        TextFormatted(@sprintf("Contact Normal to Strut Angle: %.7f deg", rad2deg(α_cs)))
        GUI.draw(v_ec_xy, "Contact Point Velocity (Oc / ECEF) [Contact]", "m/s")
        TextFormatted("Surface Type: $surface")

end

################################################################################
################################### Contact ####################################

@kwdef struct Contact <: ModelDefinition
    frc::PIVector{2} = PIVector{2}() #friction constraint compensator
end

@kwdef struct ContactY #defaults should be consistent with wow = 0
    μ_roll::Float64 = 0.0 #rolling friction coefficient
    μ_skid::Float64 = 0.0 #skidding friction coefficient
    κ_br::Float64 = 0.0 #braking factor
    ψ_cv::Float64 = 0.0 #tire slip angle
    μ_max::SVector{2,Float64} = zeros(SVector{2}) #maximum friction coefficient
    μ_eff::SVector{2,Float64} = zeros(SVector{2}) #effective friction coefficient
    f_c::SVector{3,Float64} = zeros(SVector{3}) #normalized contact force
    F_c::SVector{3,Float64} = zeros(SVector{3}) #contact force
    wr_b::Wrench = Wrench() #resulting Wrench on the vehicle frame
    frc::PIVectorY{2} = PIVectorY{2}() #contact friction regulator
end

Modeling.Y(::Contact) = ContactY()

@no_periodic Contact

function Modeling.f_init!(mdl::Model{Contact})
    #set up friction constraint compensator
    (; k_p, k_i, k_l, bound_lo, bound_hi) = mdl.frc.parameters
    k_p .= 5.0
    k_i .= 400.0
    k_l .= 0.2
    bound_lo .= -1
    bound_hi .= 1
end

function Modeling.f_ode!(mdl::Model{Contact},
                        strut::Model{<:Strut},
                        braking::Model{<:AbstractBraking})

    (; wow, F_dmp_zs, t_sc, t_bc, v_ec_xy, surface) = strut.y

    frc = mdl.frc
    frc.u.input .= -v_ec_xy #if !wow, v_ec_xy = [0,0]
    f_ode!(frc)

    if !wow
        mdl.y = ContactY(; frc = frc.y) #everything else default
        return
    end

    norm_v = norm(v_ec_xy)

    μ_roll = get_μ(FrictionCoefficients(Rolling(), surface), norm_v)
    μ_skid = get_μ(FrictionCoefficients(Skidding(), surface), norm_v)

    #longitudinal friction coefficient
    κ_br = get_braking_factor(braking)
    μ_x = μ_roll + (μ_skid - μ_roll) * κ_br

    #tire slip angle
    if norm_v < 1e-3 #prevents chattering in μ_y for near-zero contact velocity
        ψ_cv = π/2 #pure sideslip
    else
        ψ_cv = atan(v_ec_xy[2], v_ec_xy[1])
    end

    #lateral friction coefficient
    ψ_skid = deg2rad(10)
    ψ_abs = abs(ψ_cv)

    if ψ_abs < ψ_skid
        μ_y = μ_skid * ψ_abs / ψ_skid
    elseif ψ_abs > π - ψ_skid
        μ_y = μ_skid * (1 - (ψ_skid + ψ_abs - π)/ ψ_skid)
    else
        μ_y = μ_skid
    end

    μ_max = @SVector [μ_x, μ_y]
    μ_max *= min(1, μ_skid / norm(μ_max)) #scale μ_max so norm(μ_max) does not exceed μ_skid

    #scale μ_max with the feedback from the friction constraint compensator
    μ_eff = frc.y.output .* μ_max

    #normalized contact force projected on the contact frame
    f_c = SVector{3,Float64}(μ_eff[1], μ_eff[2], -1)
    f_s = t_sc.q(f_c) #project normalized force onto the strut frame
    # @assert f_s[3] < 0

    #the value of the ground's normal force must be such that its projection
    #along the strut cancels the damper's force
    N = -F_dmp_zs / f_s[3]
    N = max(0, N) #clamp negative values (could occur for large ξ_dot >0)
    F_c = f_c * N

    wr_c = Wrench(F = F_c)
    wr_b = t_bc(wr_c)

    mdl.y = ContactY(; μ_roll, μ_skid, κ_br, ψ_cv, μ_max, μ_eff, f_c, F_c, wr_b, frc = frc.y)

end

#here wow has its final value for the current integration step
function Modeling.f_step!(contact::Model{<:Contact}, strut::Model{<:Strut})

    !strut.y.wow && f_init!(contact.frc) #if !wow, reset friction regulator

end


################################# GUI ##########################################

function GUI.draw(mdl::Model{<:Contact}, window_label::String = "Contact")

    (; μ_roll, μ_skid, μ_max, κ_br, ψ_cv, μ_eff, f_c, F_c, wr_b) = mdl.y
    frc = mdl.frc

        TextFormatted(@sprintf("Rolling Friction Coefficient: %.7f", μ_roll))
        TextFormatted(@sprintf("Skidding Friction Coefficient: %.7f", μ_skid))
        TextFormatted(@sprintf("Braking Factor: %.7f", κ_br))
        TextFormatted(@sprintf("Tire Slip Angle: %.7f deg", rad2deg(ψ_cv)))
        GUI.draw(μ_max, "Maximum Friction Coefficient")
        GUI.draw(μ_eff, "Effective Friction Coefficient")
        GUI.draw(f_c, "Normalized Contact Force [Contact]")
        GUI.draw(F_c, "Contact Force [Contact]", "N")

        if TreeNode("Friction Regulator")
            GUI.draw(frc, window_label)
            TreePop()
        end

end


################################################################################
############################ LandingGearUnit ###################################

@kwdef struct LandingGearUnit{S <:AbstractSteering,
                                B <:AbstractBraking,
                                L <:Strut} <: ModelDefinition
    steering::S = NoSteering()
    braking::B = NoBraking()
    strut::L = Strut()
    contact::Contact = Contact()
end

@no_periodic LandingGearUnit

function Modeling.f_ode!(mdl::Model{<:LandingGearUnit},
                        terrain::Model{<:AbstractTerrain},
                        kin_data::KinData)

    (; strut, contact, steering, braking) = mdl

    f_ode!(steering)
    f_ode!(braking)
    f_ode!(strut, steering, terrain, kin_data)
    f_ode!(contact, strut, braking)

    f_output!(mdl)

end

function Modeling.f_step!(mdl::Model{<:LandingGearUnit})

    (; strut, contact, steering, braking) = mdl

    f_step!(steering)
    f_step!(braking)
    f_step!(strut)
    f_step!(contact, strut)

end

Dynamics.get_mp_b(::Model{<:LandingGearUnit}) = MassProperties()
Dynamics.get_hr_b(::Model{<:LandingGearUnit}) = zeros(SVector{3})
Dynamics.get_wr_b(mdl::Model{<:LandingGearUnit}) = mdl.y.contact.wr_b


################################################################################
################################# GUI ##########################################

function GUI.draw(mdl::Model{<:LandingGearUnit}, p_open::Ref{Bool} = Ref(true),
                 window_label::String = "Landing Gear Unit")

    BeginWindow(window_label, p_open)
        if TreeNode("Strut")
            GUI.draw(mdl.strut)
            TreePop()
        end
        if TreeNode("Steering")
            GUI.draw(mdl.steering)
            TreePop()
        end
        if TreeNode("Braking")
            GUI.draw(mdl.braking)
            TreePop()
        end
        if TreeNode("Contact")
            GUI.draw(mdl.contact)
            TreePop()
        end
    EndWindow()

end

end #module
