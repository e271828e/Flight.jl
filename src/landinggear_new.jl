module LandingGear

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack

using Flight.Modeling
using Flight.Plotting
using Flight.Attitude
using Flight.Geodesy
using Flight.Terrain
using Flight.Kinematics
using Flight.Dynamics
using Flight.Friction

import Flight.Modeling: init, f_cont!, f_disc!
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_wr_b
import Flight.Plotting: plots
import Flight.Friction: FrictionParameters

export LandingGearUnit, Strut, SimpleDamper, NoSteering, NoBraking, DirectSteering, DirectBraking

export DefaultFriction, FrictionRegulator

#basis elements, for convenience
const e1 = SVector{3,Float64}(1,0,0)
const e3 = SVector{3,Float64}(0,0,1)

########################### Steering #############################

abstract type AbstractSteering <: SystemDescriptor end

############## NoSteering ##############

struct NoSteering <: AbstractSteering end

get_steering_angle(::System{NoSteering}) = 0.0
f_cont!(::System{NoSteering}, args...) = nothing
f_disc!(::System{NoSteering}, args...) = false


############## DirectSteering ##############

Base.@kwdef struct DirectSteering <: AbstractSteering
    ψ_max::Float64 = π/6
end

Base.@kwdef struct DirectSteeringY
    ψ::Float64 = 0.0
end
#the contents of u must be mutable
init(::DirectSteering, ::SystemU) = Ref(Bounded(0.0, -1, 1))
init(::DirectSteering, ::SystemY) = DirectSteeringY(0.0) #steering angle

function f_cont!(sys::System{DirectSteering})
    sys.y = DirectSteeringY(Float64(sys.u[]) * sys.params.ψ_max)
end

f_disc!(::System{DirectSteering}, args...) = false

get_steering_angle(sys::System{DirectSteering}) = sys.y.ψ


########################### Braking #############################

abstract type AbstractBraking <: SystemDescriptor end

############## NoBraking ###############

struct NoBraking <: AbstractBraking end

get_braking_factor(::System{NoBraking}) = 0.0
f_cont!(::System{NoBraking}, args...) = nothing
f_disc!(::System{NoBraking}, args...) = false

########### DirectBraking #############

Base.@kwdef struct DirectBraking <: AbstractBraking
   η_br::Float64 = 1.0 #braking efficiency
end

Base.@kwdef struct DirectBrakingY
   κ_br::Float64 = 0.0 #braking coefficient
end

init(::DirectBraking, ::SystemU) = Ref(Bounded(0.0, 0, 1))
init(::DirectBraking, ::SystemY) = DirectBrakingY()

function f_cont!(sys::System{DirectBraking})
    sys.y = DirectBrakingY(Float64(sys.u[]) * sys.params.η_br)
end

f_disc!(::System{DirectBraking}, args...) = false

get_braking_factor(sys::System{DirectBraking}) = sys.y.κ_br


########################### Damper #############################

abstract type AbstractDamper end #not a System!
get_damper_force(args...) = throw(MethodError(get_damper_force, args))

Base.@kwdef struct SimpleDamper <: AbstractDamper
    l_0::Float64 = 0.0 #natural length
    k_s::Float64 = 25000 #spring constant
    k_d_ext::Float64 = 1000 #extension damping coefficient
    k_d_cmp::Float64 = 1000 #compression damping coefficient
    ξ_min::Float64 = -5 #compression below which the shock absorber is disabled
end

get_natural_length(c::SimpleDamper) = c.l_0

#Force exerted by the damper along zs. The deformation ξ is positive along z_s
#(elongation). The resulting force can be negative, meaning the damper pulls the
#piston rod assembly upwards along the negative z_s axis. This can happen when
#ξ_dot > 0 and sufficiently large
function get_damper_force(c::SimpleDamper, ξ::Real, ξ_dot::Real)
    k_d = (ξ_dot > 0 ? c.k_d_ext : c.k_d_cmp)
    F = -(c.k_s * ξ + k_d * ξ_dot)
    F = F * (ξ > c.ξ_min)
    return F
end


########################### FrictionParameters ################################

abstract type RollingOrSkidding end
struct Rolling <: RollingOrSkidding end
struct Skidding <: RollingOrSkidding end

FrictionParameters(::Rolling, ::SurfaceType) =
    FrictionParameters(μ_s = 0.03, μ_d = 0.02, v_s = 0.005, v_d = 0.01)

FrictionParameters(::Skidding, ::DryTarmac) =
    FrictionParameters( μ_s = 0.75, μ_d = 0.25, v_s = 0.005, v_d = 0.01)

FrictionParameters(::Skidding, ::WetTarmac) =
    FrictionParameters(μ_s = 0.25, μ_d = 0.15, v_s = 0.005, v_d = 0.01)

FrictionParameters(::Skidding, ::IcyTarmac) =
    FrictionParameters(μ_s = 0.075, μ_d = 0.025, v_s = 0.005, v_d = 0.01)



########################## Strut #########################

Base.@kwdef struct Strut{D<:AbstractDamper} <: SystemDescriptor
    t_bs::FrameTransform = FrameTransform() #airframe to strut frame transform
    l_OsP::Float64 = 0.0 #strut natural length from strut frame origin to wheel endpoint
    damper::D = SimpleDamper()
end

Base.@kwdef struct StrutY
    wow::Bool = false #weight-on-wheel flag
    ξ::Float64 = 0.0 #damper deformation
    ξ_dot::Float64 = 0.0 #damper deformation rate
    t_sc::FrameTransform = FrameTransform() #strut to contact frame transform
    t_bc::FrameTransform = FrameTransform() #body to contact frame transform
    F::Float64 = 0.0 #damper force on the ground along the strut z axis
    v_eOc_c::SVector{3,Float64} = zeros(SVector{3}) #contact point velocity
    srf::SurfaceType = Terrain.DryTarmac #surface type at the contact reference point
end

init(::Strut, ::SystemY) = StrutY()

function f_cont!(sys::System{<:Strut}, steering::System{<:AbstractSteering},
    terrain::AbstractTerrain, kin::KinData)

    @unpack t_bs, l_OsP, damper = sys.params
    @unpack n_e, h_e, q_eb, q_nb = kin.pos
    @unpack v_eOb_b, ω_eb_b = kin.vel

    #strut frame axes
    q_bs = t_bs.q
    r_ObOs_b = t_bs.r #strut frame origin

    #compute contact reference point P
    r_OsP_s = l_OsP * e3
    r_ObP_b = r_ObOs_b + q_bs(r_OsP_s)
    r_ObP_e = q_eb(r_ObP_b)
    r_OeOb_e = CartECEF(Geographic(n_e, h_e))
    r_OeP_e = r_OeOb_e + r_ObP_e
    P = Geographic(r_OeP_e)

    #get terrain data at the contact reference point (close enough to the actual
    #contact frame origin, which is still unknown)
    trn = get_terrain_data(terrain, P.l2d)
    h_trn = trn.altitude
    srf = trn.surface

    #compute damper deformation. if the projection of k_s onto z_n is close to
    #zero (strut nearly horizontal) or negative (strut upside down), set it to 0
    q_ns = q_nb ∘ q_bs
    k_s_zn = q_ns(e3)[3]
    Δh = AltOrth(P) - h_trn
    if k_s_zn > 1e-3
        ξ = min(0.0, Δh / k_s_zn)
    else
        ξ = 0.0 #0 instead of 0.0 causes type instability
    end
    wow = ξ < 0

    if !wow #we're done, assign only relevant outputs and return
        sys.y = StrutY(; wow, ξ, srf)
        return
    end

    #wheel frame axes
    ψ_sw = get_steering_angle(steering)
    q_sw = Rz(ψ_sw) #rotate strut axes to get wheel axes
    q_nw = q_ns ∘ q_sw #NED to contact axes rotation

    #contact frame axes
    k_c_n = trn.normal #NED components of contact z-axis
    i_w_n = q_nw(e1) #NED components of wheel x-axis
    i_w_n_trn = i_w_n - (i_w_n ⋅ k_c_n) * k_c_n #projection of i_w_n onto the terrain tangent plane
    i_c_n = normalize(i_w_n_trn) #NED components of contact x-axis
    j_c_n = k_c_n × i_c_n #NED components of contact y-axis
    R_nc = RMatrix(SMatrix{3,3}([i_c_n j_c_n k_c_n]), normalization = false)
    q_sc = q_ns' ∘ R_nc
    q_bc = q_bs ∘ q_sc

    #contact frame origin
    r_OsOc_s = e3 * (l_0 + ξ)
    r_OsOc_b = q_bs(r_OsOc_s)
    r_ObOc_b = r_OsOc_b + r_ObOs_b

    #construct contact frame transforms
    t_sc = FrameTransform(r_OsOc_s, q_sc)
    t_bc = FrameTransform(r_ObOc_b, q_bc)

    #contact frame origin velocity due to airframe motion
    v_eOc_afm_b = v_eOb_b + ω_eb_b × r_ObOc_b
    v_eOc_afm_c = q_bc'(v_eOc_afm_b)

    #compute the damper elongation rate required to cancel the airframe
    #contribution to the contact point velocity along the contact frame z axis
    q_sc = q_bs' ∘ q_bc
    k_s_c = q_sc'(e3)
    ξ_dot = -v_eOc_afm_c[3] / k_s_c[3]

    #compute contact point velocity
    v_eOc_dmp_c = k_s_c * ξ_dot
    v_eOc_c = v_eOc_afm_c + v_eOc_dmp_c

    F = get_damper_force(damper, ξ, ξ_dot)

    sys.y = StrutY(; wow, ξ, ξ_dot, t_sc, t_bc, F, v_eOc_c, srf)

end


f_disc!(::System{<:Strut}, args...) = false

########################## Contact ###############################


Base.@kwdef struct Contact <: SystemDescriptor
    friction::FrictionRegulator{2} = FrictionRegulator{2}(k_p = 5.0, k_i = 400.0)
end

Base.@kwdef struct ContactY
    friction::FrictionRegulatorY{2} = FrictionRegulatorY{2}()
    μ_roll::Float64 = 0.0 #rolling friction coefficient
    μ_skid::Float64 = 0.0 #skidding friction coefficient
    κ_br::Float64 = 0.0 #braking factor
    ψ_cv::Float64 = 0.0 #tire slip angle
    μ_max::SVector{2,Float64} = zeros(SVector{2}) #maximum friction coefficient
    μ_eff::SVector{2,Float64} = zeros(SVector{2}) #scaled friction coefficient
    f_c::SVector{3,Float64} = zeros(SVector{3}) #normalized contact force
    F_c::SVector{3,Float64} = zeros(SVector{3}) #contact force
    wr_b::Wrench = Wrench() #resulting airframe Wrench
end

#x should be initialized by the fallback methods
init(::Contact, ::SystemY) = ContactY()

function f_cont!(sys::System{Contact}, strut::System{<:Strut},
                braking::System{<:AbstractBraking})

    @unpack wow, t_sc, t_bc, F, v_eOc_c, srf = strut.y
    friction = sys.friction

    if !wow
        sys.y = ContactY(friction = friction.y) #rest of fields to default
        return
    end

    v = SVector{2,Float64}(v_eOc_c[1], v_eOc_c[2]) #contact point velocity

    μ_roll = get_μ(FrictionParameters(srf, Rolling()), norm(v))
    μ_skid = get_μ(FrictionParameters(srf, Skidding()), norm(v))

    #longitudinal friction coefficient
    κ_br = get_braking_factor(braking)
    μ_x = μ_roll + (μ_skid - μ_roll) * κ_br

    #tire slip angle
    if norm(v) < 1e-3 #prevents chattering in μ_y for near-zero contact velocity
        ψ_cv = π/2 #pure sideslip
    else
        ψ_cv = atan(v[2], v[1])
    end

    #lateral friction coefficient
    ψ_skid = deg2rad(10) #could be defined as a descriptor parameter
    ψ_abs = abs(ψ_cv)

    if ψ_abs < ψ_skid
        μ_y = μ_skid * ψ_abs / ψ_skid
    elseif ψ_abs > π - ψ_skid
        μ_y = μ_skid * (1 - (ψ_skid + ψ_abs - π)/ ψ_skid)
    else
        μ_y = μ_skid
    end

    #maximum friction coefficient vector
    μ_max = @SVector [μ_x, μ_y]

    #norm(μ_max) cannot exceed μ_skid; if it does, scale down its components
    μ_max *= min(1, μ_skid / norm(μ_max))

    #update friction regulator
    f_cont!(friction, v)

    #retrieve its scaling factor
    μ_eff = friction.y.α .* μ_max

    ############# HASTA AQUI
    ############# HASTA AQUI
    ############# HASTA AQUI
    ############# HASTA AQUI

    #normalized contact force projected on the contact frame
    f_c = SVector{3,Float64}(μ[1], μ[2], -1)
    f_s = t_sc.q(f_c) #project normalized force onto the strut frame
    @assert f_s[3] < 0

    N = -F / f_s[3]
    N = max(0, N) #cannot be negative (could happen with large xi_dot >0)
    F_c = f_c * N

    wr_c = Wrench(F = F_c)
    wr_b = t_bc(wr_c)

    sys.y = ContactY(; friction, μ_roll, μ_skid, κ_br, ψ_cv, μ_max, μ_eff, f_c, F_c, wr_b)

end

function f_disc!(contact::System{Contact}, strut::System{<:Strut})
    #if !wow, reset regulator and return x_mod flag
    return f_disc!(contact.friction, !strut.y.wow)
end


########################## LandingGearUnit #########################

Base.@kwdef struct LandingGearUnit{L<:Strut, S <: AbstractSteering,
                            B <: AbstractBraking} <: SystemDescriptor
    strut::L = Strut()
    contact::Contact = Contact()
    steering::S = NoSteering()
    braking::B = NoBraking()
end

#if we avoid the generic fallback for SystemGroup, we don't need to define
#traits for Steering, Braking, Contact or Strut
MassTrait(::System{<:LandingGearUnit}) = HasNoMass()
WrenchTrait(::System{<:LandingGearUnit}) = GetsExternalWrench()
AngularMomentumTrait(::System{<:LandingGearUnit}) = HasNoAngularMomentum()

get_wr_b(sys::System{<:LandingGearUnit}) = sys.y.contact.wr_b


function f_cont!(sys::System{<:LandingGearUnit}, kinematics::KinData,
                terrain::AbstractTerrain)

    @unpack strut, contact, steering, braking = sys.subsystems

    f_cont!(steering)
    f_cont!(braking)
    f_cont!(strut, steering, terrain, kinematics)
    f_cont!(contact, strut, braking)

    Modeling.assemble_y!(sys)

end

function f_disc!(sys::System{<:LandingGearUnit})

    return f_disc!(sys.contact, sys.strut)

end

#need to override the default AirframeNode implementation, since the underlying
#subcomponents are not AirframeComponents and only Contact contributes wr_b.
#could also define trivial get_wr_b and get_hr_b methods for steering, braking
#and strut, and let the default AirframeNode implementation do its thing
# get_wr_b(sys::System{<:LandingGearUnit}) = sys.y.contact.wr_b

end #module