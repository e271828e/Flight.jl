module LandingGear
#rename to ground, rename LandingGearUnit to LandingGear

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack

using Flight.ModelingTools
using Flight.Attitude
using Flight.Geodesy
using Flight.Kinematics
using Flight.Dynamics
using Flight.Airframe
using Flight.Terrain

#for method extension without module qualifier
import Flight.ModelingTools: System, get_x0, get_y0, get_u0, get_d0, f_cont!, f_disc!
import Flight.Airframe: get_wr_b, get_hr_b

export Strut, StrutY
export Contact, ContactY
export SimpleDamper, get_damper_force
export NoSteering, DirectSteering, DirectSteeringY, get_steering_angle, set_steering_input
export NoBraking, DirectBraking, DirectBrakingY, get_braking_coefficient, set_braking_input
export LandingGearUnit, LandingGearUnitY

#basis elements, for convenience
const e1 = SVector{3,Float64}(1,0,0)
const e3 = SVector{3,Float64}(0,0,1)

########################### Steering #############################

abstract type AbstractSteering <: AbstractComponent end

############## NoSteering ##############

struct NoSteering <: AbstractSteering end

set_steering_input(::System{NoSteering}, ::Real) = nothing
get_steering_angle(::System{NoSteering}) = 0.0
f_cont!(::System{NoSteering}, args...) = nothing


############## DirectSteering ##############

Base.@kwdef struct DirectSteering <: AbstractSteering
    ψ_max::Float64 = π/6
end

Base.@kwdef struct DirectSteeringY
    ψ::Float64 = 0.0
end
#we need to make the contents of u mutable
get_u0(::DirectSteering) = Ref(0.0) #steering input
get_y0(::DirectSteering) = DirectSteeringY(0.0) #steering angle

function set_steering_input(sys::System{DirectSteering}, u::Real)
    @assert abs(u) <= 1 "Steering input must be within [-1,1]"
    sys.u[] = u
end

function f_cont!(sys::System{DirectSteering})
    sys.y = DirectSteeringY(sys.u[] * sys.params.ψ_max)
end

get_steering_angle(sys::System{DirectSteering}) = sys.y.ψ

# struct ActuatedSteering{A <: AbstractActuator} <: AbstractSteering
#     limits::SVector{2,Float64} #more generally, transmission kinematics could go here
#     actuator::A #actuator model parameters go here
# end
# get_x0(steering::ActuatedSteering) = get_x0(steering.actuator)
# get_y0(steering::ActuatedSteering) = get_x0(steering.actuator)
# get_u0(steering::ActuatedSteering) = get_u0(steering.actuator) #typically, 1

########################### Braking #############################

abstract type AbstractBraking <: AbstractComponent end

############## NoBraking ###############

struct NoBraking <: AbstractBraking end

f_cont!(::System{NoBraking}, args...) = nothing
set_braking_input(::System{NoBraking}, ::Real) = nothing
get_braking_coefficient(::System{NoBraking}) = 0.0

########### DirectBraking #############

Base.@kwdef struct DirectBraking <: AbstractBraking
   η_br::Float64 = 1.0 #braking efficiency
end

Base.@kwdef struct DirectBrakingY
   α_br::Float64 = 0.0 #braking coefficient
end

get_u0(::DirectBraking) = Ref(0.0)
get_y0(::DirectBraking) = DirectBrakingY()

function set_braking_input(sys::System{DirectBraking}, u::Real)
    @assert (u <= 1 && u >= 0) "Braking input must be within [0,1]"
    sys.u[] = u
end

function f_cont!(sys::System{DirectBraking})
    sys.y = DirectBrakingY(sys.u[] * sys.params.η_br)
end

get_braking_coefficient(sys::System{DirectBraking}) = sys.y.α_br

########################### Damper #############################

abstract type AbstractDamper end #not a System!
get_damper_force(::A, args...) where {A<:AbstractDamper} = no_extend_error(get_force, A)

Base.@kwdef struct SimpleDamper <: AbstractDamper
    k_s::Float64 = 25000 #spring constant
    k_d_ext::Float64 = 1000 #extension damping coefficient
    k_d_cmp::Float64 = 1000 #compression damping coefficient
    ξ_min::Float64 = -1 #compression below which the shock absorber is disabled
end

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

########################## Strut #########################

Base.@kwdef struct Strut{D<:AbstractDamper} <: AbstractComponent
    t_bs::FrameTransform = FrameTransform()
    damper::D = SimpleDamper()
    l_0::Float64 = 1.0 #natural length
end

Base.@kwdef struct StrutY
    wow::Bool = false #weight-on-wheel flag
    ξ::Float64 = 0.0 #damper deformation
    ξ_dot::Float64 = 0.0 #damper deformation rate
    t_sc::FrameTransform = FrameTransform() #strut to contact frame transform
    t_bc::FrameTransform = FrameTransform() #body to contact frame transform
    F_dmp::Float64 = 0.0 #damper force on the ground along the strut z axis
    v_eOc_c::SVector{3,Float64} = zeros(SVector{3}) #contact point velocity
    srf::SurfaceType = Terrain.DryTarmac #surface type at the reference point
end

get_y0(::Strut) = StrutY()

function f_cont!(sys::System{<:Strut}, steering::System{<:AbstractSteering},
    terrain::AbstractTerrain, kin::KinData)

    @unpack t_bs, damper, l_0 = sys.params
    @unpack n_e, h_e, q_eb, q_nb = kin.pos
    @unpack v_eOb_b, ω_eb_b = kin.vel

    #strut frame axes
    q_bs = t_bs.q
    r_ObOs_b = t_bs.r #strut frame origin

    #compute contact reference point P
    r_OsP_s = l_0 * e3
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
    ξ = (k_s_zn > 1e-3 ? Δh / k_s_zn : 0.0) #0 instead of 0.0 causes type instability
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

    F_dmp = get_damper_force(damper, ξ, ξ_dot)

    sys.y = StrutY(; wow, ξ, ξ_dot, t_sc, t_bc, F_dmp, v_eOc_c, srf)

end

########################### Contact #############################

#static / dynamic friction coefficient set
struct StaticDynamic
    static::Float64
    dynamic::Float64
end

#static / dynamic friction interpolation from breakout coefficient
get_μ(f::StaticDynamic, α_bo::Real)::Float64 = α_bo * f.dynamic + (1-α_bo) * f.static

#pure rolling friction coefficient set
Base.@kwdef struct Rolling
    sd::StaticDynamic= StaticDynamic(0.03, 0.02)
end

get_μ(f::Rolling, ::SurfaceType)::StaticDynamic = f.sd

#pure skidding friction coefficient set
Base.@kwdef struct Skidding
    dry::StaticDynamic = StaticDynamic(0.75, 0.25)
    wet::StaticDynamic = StaticDynamic(0.25, 0.15)
    icy::StaticDynamic = StaticDynamic(0.075, 0.025)
end

function get_μ(f::Skidding, srf::SurfaceType)::StaticDynamic
    if srf === Terrain.DryTarmac
        return f.dry
    elseif srf === Terrain.WetTarmac
        return f.wet
    elseif srf === Terrain.IcyTarmac
        return f.icy
    else
        error("Invalid surface type")
    end
end

get_μ(f::Union{Rolling,Skidding}, srf::SurfaceType, α_bo::Real)::Float64 =
    get_μ(get_μ(f, srf), α_bo)

#this allocates. why??
# get_sd(f::Skidding, ::Val{Terrain.DryTarmac}) = f.dry
# get_sd(f::Skidding, ::Val{Terrain.WetTarmac}) = f.wet
# get_sd(f::Skidding, ::Val{Terrain.IcyTarmac}) = f.icy

Base.@kwdef struct Contact <: AbstractComponent
    rolling::Rolling = Rolling() #rolling friction coefficients
    skidding::Skidding = Skidding() #skidding friction coefficients
    v_bo::SVector{2,Float64} = [0.005, 0.01] #breakout velocity interval
    ψ_skid::Float64 = deg2rad(10) #skidding slip angle thresold
    k_p::Float64 = 5.0 #proportional gain for contact velocity regulator
    k_i::Float64 = 400.0 #integral gain for contact velocity regulator
end

Base.@kwdef struct ContactY
    v::SVector{2,Float64} = zeros(SVector{2}) #contact plane velocity
    s::SVector{2,Float64} = zeros(SVector{2}) #contact plane velocity integral
    α_p::SVector{2,Float64} = zeros(SVector{2}) #proportional μ scale factor
    α_i::SVector{2,Float64} = zeros(SVector{2}) #integral μ scale factor
    α_raw::SVector{2,Float64} = zeros(SVector{2}) #total scale factor, raw
    α::SVector{2,Float64} = zeros(SVector{2}) #total scale factor, clipped
    sat::SVector{2,Bool} = zeros(SVector{2,Bool}) #scale factor saturation flag
    ψ_cv::Float64 = 0.0 #tire slip angle
    α_bo::Float64 = 0.0 #breakout coefficient
    μ_max::SVector{2,Float64} = zeros(SVector{2}) #maximum friction coefficient
    μ::SVector{2,Float64} = zeros(SVector{2}) #scaled friction coefficient
    f_c::SVector{3,Float64} = zeros(SVector{3}) #non-dimensional contact force
    F_c::SVector{3,Float64} = zeros(SVector{3}) #contact force
    wr_b::Wrench = Wrench() #resulting airframe Wrench
end

get_x0(::Contact) = ComponentVector(x = 0.0, y = 0.0) #v regulator integrator states
get_y0(::Contact) = ContactY()

function f_cont!(sys::System{Contact}, strut::System{<:Strut},
                braking::System{<:AbstractBraking})


    @unpack rolling, skidding, v_bo, ψ_skid, k_p, k_i = sys.params
    @unpack wow, t_sc, t_bc, F_dmp, v_eOc_c, srf = strut.y

    if !wow
        sys.ẋ .= 0
        sys.y = ContactY()
        return
    end

    #breakout coefficient
    α_bo = clamp((norm(v_eOc_c) - v_bo[1]) / (v_bo[2] - v_bo[1]), 0, 1)

    μ_roll = get_μ(rolling, srf, α_bo)
    μ_skid = get_μ(skidding, srf, α_bo)

    #longitudinal friction coefficient
    α_br = get_braking_coefficient(braking)
    @assert (α_br >= 0 && α_br <= 1)
    μ_x = μ_roll + (μ_skid - μ_roll) * α_br

    #lateral friction coefficient
    ψ_cv = atan(v_eOc_c[2], v_eOc_c[1]) #tire slip angle
    ψ_abs = abs(ψ_cv)
    @assert (ψ_abs <= π)

    if α_bo < 1
        μ_y = μ_skid
    else
        if ψ_abs < ψ_skid
            μ_y = μ_skid * ψ_abs / ψ_skid
        elseif ψ_abs > π - ψ_skid
            μ_y = μ_skid * (1 - (ψ_skid + ψ_abs - π)/ ψ_skid)
        else
            μ_y = μ_skid
        end
    end

    #maximum friction coefficient vector
    μ_max = @SVector [μ_x, μ_y]

    #the magnitude of μ_max cannot exceed μ_skid; if it does, scale down its components
    μ_max *= min(1, μ_skid / norm(μ_max))

    v = SVector{2,Float64}(v_eOc_c[1], v_eOc_c[2]) #contact point velocity
    s = SVector{2,Float64}(sys.x) #velocity integrator state
    α_p = -k_p * v
    α_i = -k_i * s
    α_raw = α_p + α_i #raw μ scaling
    α = clamp.(α_raw, -1, 1) #clipped μ scaling
    μ = α .* μ_max

    #if not saturated, integrator accumulates
    sat = abs.(α_raw) .> abs.(α) #saturated?
    sys.ẋ .= v .* sat

    #non-dimensional contact force projected on the contact frame
    f_c = SVector{3,Float64}(μ[1], μ[2], -1)
    f_s = t_sc.q(f_c) #project non-dimensional force onto the strut frame
    @assert f_s[3] < 0

    N = -F_dmp / f_s[3]
    N = max(0, N) #cannot be negative (could happen with large xi_dot >0)
    F_c = f_c * N

    wr_c = Wrench(F = F_c)
    wr_b = t_bc(wr_c)

    sys.y = ContactY(; v, s, α_p, α_i, α_raw, α, sat, ψ_cv, α_bo, μ_max, μ, f_c, F_c, wr_b)

end

########################## LandingGearUnit #########################

Base.@kwdef struct LandingGearUnit{L<:Strut, S <: AbstractSteering,
                            B <: AbstractBraking} <: AbstractAirframeComponent
    strut::L = Strut()
    contact::Contact = Contact()
    steering::S = NoSteering()
    braking::B = NoBraking()
end

struct LandingGearUnitY{SteeringY, BrakingY}
    strut::StrutY
    contact::ContactY
    steering::SteeringY
    braking::BrakingY
end

get_y0(ldg::LandingGearUnit) = LandingGearUnitY(
    get_y0(ldg.strut),
    get_y0(ldg.contact),
    get_y0(ldg.steering),
    get_y0(ldg.braking),
    )

get_x0(ldg::LandingGearUnit) = ComponentVector(
    strut = get_x0(ldg.strut),
    contact = get_x0(ldg.contact),
    steering = get_x0(ldg.steering),
    braking = get_x0(ldg.braking),
    )

get_u0(ldg::LandingGearUnit) = (
    strut = get_u0(ldg.strut),
    contact = get_u0(ldg.contact),
    steering = get_u0(ldg.steering),
    braking = get_u0(ldg.braking),
    )

get_d0(ldg::LandingGearUnit) = (
    strut = get_d0(ldg.strut),
    contact = get_d0(ldg.contact),
    steering = get_d0(ldg.steering),
    braking = get_d0(ldg.braking),
    )


function System(ldg::LandingGearUnit, ẋ = get_x0(ldg), x = get_x0(ldg),
                    y = get_y0(ldg), u = get_u0(ldg), d = get_d0(ldg), t = Ref(0.0))

    ss_list = Vector{System}()
    ss_labels = (:strut, :contact, :steering, :braking)
    for label in ss_labels
        push!(ss_list, System(map((λ)->getproperty(λ, label), (ldg, ẋ, x, y, u, d))..., t))
    end

    params = nothing
    subsystems = NamedTuple{ss_labels}(ss_list)

    System{map(typeof, (ldg, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)
end

function f_cont!(sys::System{<:LandingGearUnit}, kinematics::KinData,
                terrain::AbstractTerrain)

    @unpack strut, contact, steering, braking = sys.subsystems

    #update steering and braking subsystems
    f_cont!(steering)
    f_cont!(braking)
    f_cont!(strut, steering, terrain, kinematics)
    f_cont!(contact, strut, braking)

    sys.y = LandingGearUnitY(strut.y, contact.y, steering.y, braking.y)

end

get_wr_b(sys::System{<:LandingGearUnit}) = sys.y.contact.wr_b
get_hr_b(::System{<:LandingGearUnit}) = zeros(SVector{3})

#TODO: add contact state reset for wow=false
f_disc!(::System{<:LandingGearUnit}) = false


end #module

#in System, define and extend f_branch!

# #individual Component
# f_branch!(y, dx, x, u, t, sys, args...) = f_branch!(Val(has_input(sys)), y, dx, x, u, t, args...)
# f_branch!(::Val{true}, y, dx, x, u, t, sys, args...) = f_cont!(y, dx, x, u, t, sys, args...)
# f_cont!(::HasInput, y, dx, x ,u, t, sys, args...) = f_cont!(y, dx, x, u, t, sys, args...)
# f_cont!(::HasNoInput, y, dx, x, u, t, sys, args...) = f_cont!(y, dx, x, t, sys, args...)

# #for a AirframeGroup
# f_cont!(MaybeInput(S), MaybeOutput(S), y, dx, x, u, t, sys, args...)
# f_cont!(::HasInput, ::HasOutput, y, dx, x ,u, t, sys, args...)
# #now, this method needs to consider the possibility for each component that it
# #may have or not Input or Output. so it must do
# for (label, component) in zip(keys(C), values(C))
#     if MaybeInput(typeof(component)) #need tocheck, because if it has no input, u[label] will not exist!
#         f_cont!(y_cmp, dx_cmp, x_cmp, u_cmp, t, cmp, args...)
#     else
#         f_cont!(y_cmp, dx_cmp, x_cmp, t, cmp, args...)
#     end
# end


#this method