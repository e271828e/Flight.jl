module LandingGear
#rename to ground, rename LandingGearLeg to LandingGear

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack

using Flight.ModelingTools
# using Flight.Airdata
using Flight.Kinematics
using Flight.Dynamics
using Flight.Airframe
using Flight.Terrain

#for method extension without module qualifier
import Flight.ModelingTools: System, get_x0, get_y0, get_u0, get_d0, f_cont!, f_disc!
import Flight.Airframe: get_wr_b, get_hr_b

export Contact, ContactY
export NoSteering, DirectSteering, get_steering_angle, set_steering_input
export NoBraking, DirectBraking, get_braking_coefficient, set_braking_input
export LandingGearLeg


########################### Contact #############################

struct WoW{x}
    WoW(x::Bool) = new{x}()
end

struct StaticDynamic
    static::Float64
    dynamic::Float64
end

Base.@kwdef struct Rolling
    sd::StaticDynamic= StaticDynamic(0.03, 0.02)
end

Base.@kwdef struct Skidding
    dry::StaticDynamic = StaticDynamic(0.75, 0.25)
    wet::StaticDynamic = StaticDynamic(0.25, 0.15)
    icy::StaticDynamic = StaticDynamic(0.075, 0.025)
end

get_μ(f::StaticDynamic, η_bo::Real)::Float64 = η_bo * f.dynamic + (1-η_bo) * f.static

get_μ(f::Rolling, ::SurfaceCondition)::StaticDynamic = f.sd

function get_μ(f::Skidding, cond::SurfaceCondition)::StaticDynamic
    if cond === Terrain.Dry
        return f.dry
    elseif cond === Terrain.Wet
        return f.wet
    elseif cond === Terrain.Icy
        return f.icy
    else
        error("Invalid surface condition")
    end
end

get_μ(f::Union{Rolling,Skidding}, cond::SurfaceCondition, η_bo::Real)::Float64 =
    get_μ(get_μ(f, cond), η_bo)

#this allocates. why??
# get_sd(f::Skidding, ::Val{Terrain.Dry}) = f.dry
# get_sd(f::Skidding, ::Val{Terrain.Wet}) = f.wet
# get_sd(f::Skidding, ::Val{Terrain.Icy}) = f.icy


Base.@kwdef struct Contact <: AbstractComponent
    rolling::Rolling = Rolling() #rolling friction coefficients
    skidding::Skidding = Skidding() #skidding friction coefficients
    v_bo::SVector{2,Float64} = [0.005, 0.01] #breakout velocity interval
    ψ_skid::Float64 = deg2rad(10) #skidding slip angle thresold
    k_p::Float64 = 5.0 #proportional gain for contact velocity regulator
    k_i::Float64 = 400.0 #integral gain for contact velocity regulator
end

Base.@kwdef struct ContactY
    v::SVector{2,Float64} = zeros(2)
    s::SVector{2,Float64} = zeros(2)
    α_p::SVector{2,Float64} = zeros(2)
    α_i::SVector{2,Float64} = zeros(2)
    α_raw::SVector{2,Float64} = zeros(2)
    α_sat::SVector{2,Bool} = zeros(2)
    α::SVector{2,Float64} = zeros(2)
    ψ_slip::Float64 = 0.0
    k_bo::Float64 = 0.0
    μ_max::SVector{2,Float64} = zeros(2)
    μ::SVector{2,Float64} = zeros(2)
end

get_x0(::Contact) = ComponentVector(x = 0.0, y = 0.0) #v regulator integrator states
get_y0(::Contact) = ContactY()


function f_cont!(sys::System{Contact}, ::WoW{true}, srf_cond::SurfaceCondition,
                 η_br::Real, v_eOc_c_in::AbstractVector{<:Real})

    @unpack rolling, skidding, v_bo, ψ_skid = sys.params

    v_eOc_c = SVector{3,Float64}(v_eOc_c_in)

    #breakout coefficient
    η_bo = clamp((norm(v_eOc_c) - v_bo[1]) / (v_bo[2] - v_bo[1]), 0, 1)

    μ_roll = get_μ(rolling, srf_cond, η_bo)
    μ_skid = get_μ(skidding, srf_cond, η_bo)

    #longitudinal friction coefficient
    @assert (η_br >= 0 && η_br <= 1)
    μ_x = μ_roll + (μ_skid - μ_roll) * η_br

    #lateral friction coefficient
    ψ_cv = atan(v_eOc_c[2], v_eOc_c[1]) #tire slip angle
    ψ_abs = abs(ψ_cv)
    @assert (ψ_abs <= π)

    if η_bo < 1
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

    # 2D friction coefficient vector
    μ = @SVector [μ_x, μ_y]

    # the magnitude of μ cannot exceed μ_skid; if it does, scale down its components
    μ *= min(1, μ_skid / norm(μ))

    # SEGUIR AQUI

end

function f_cont!(sys::System{Contact}, ::WoW{false})
    println("Tobeimplemented")
end

# ########################### Damper #############################

# abstract type AbstractDamper end #not a System!
# get_force(::A) where {A<:AbstractDamper} = no_extend_error(get_force, A)

# Base.@kwdef struct SimpleDamper <: AbstractDamper
#     l_0::Float64 = 1 #natural length
#     k_s::Float64 = 25000 #spring constant
#     k_d_ext::Float64 = 1000 #extension damping coefficient
#     k_d_cmp::Float64 = 1000 #compression damping coefficient
#     ξ_min::Float64 = -1 #compression below which the shock absorber is disabled
# end

# f_cont!(::System{SimpleDamper}, args...) = nothing

# ########################### Steering #############################

# abstract type AbstractSteering <: AbstractComponent end

# ############## NoSteering ##############

# struct NoSteering <: AbstractSteering end
# set_steering_input(::System{NoSteering}, ::Float64) = nothing
# get_steering_angle(::System{NoSteering}) = 0.0
# f_cont!(::System{NoSteering}, args...) = nothing


# ############## DirectSteering ##############

# Base.@kwdef struct DirectSteering <: AbstractSteering
#     ψ_max::Float64 = π/6
# end

# Base.@kwdef struct DirectSteeringY
#     ψ::Float64 = 0.0
# end
# #we need to make the contents of u mutable
# get_u0(::DirectSteering) = Ref(0.0) #steering input
# get_y0(::DirectSteering) = DirectSteeringY(0.0) #steering angle

# function f_cont!(sys::System{DirectSteering})
#     u = sys.u[]
#     abs(u) <= 1 ? sys.y = DirectSteeringY(u * sys.params.ψ_max) :
#         throw(ArgumentError("Steering input must be within [-1,1]"))
# end

# # struct ActuatedSteering{A <: AbstractActuator} <: AbstractSteering
# #     limits::SVector{2,Float64} #more generally, transmission kinematics could go here
# #     actuator::A #actuator model parameters go here
# # end
# # get_x0(steering::ActuatedSteering) = get_x0(steering.actuator)
# # get_y0(steering::ActuatedSteering) = get_x0(steering.actuator)
# # get_u0(steering::ActuatedSteering) = get_u0(steering.actuator) #typically, 1

# ########################### Braking #############################

# abstract type AbstractBraking <: AbstractComponent end

# ############## NoBraking ###############

# struct NoBraking <: AbstractBraking end
# f_cont!(::System{NoBraking}, args...) = nothing
# set_braking_input(::System{NoBraking}, ::Float64) = nothing
# get_braking_coefficient(::System{NoBraking}) = 0.0

# ########### DirectBraking #############

# Base.@kwdef struct DirectBraking <: AbstractBraking
#     η_brk::Float64 = 1.0 #braking efficiency
# end

# Base.@kwdef struct DirectBrakingY
#     k_brk::Float64 = 0.0 #braking coefficient
# end

# get_u0(::DirectBraking) = Ref(0.0)
# get_y0(::DirectBraking) = DirectBrakingY()

# function f_cont!(sys::System{DirectBraking})
#     u = sys.u[]
#     (u <= 1 && u >= 0) ? sys.y = DirectBrakingY(u * sys.params.η_brk) :
#         throw(ArgumentError("Braking Input must be within [0,1]"))
# end

# ########################## LandingGearLeg #########################

# Base.@kwdef struct LandingGearLeg{D<:AbstractDamper, S <: AbstractSteering,
#                             B <: AbstractBraking} <: AbstractAirframeComponent
#     frame::FrameSpec = FrameSpec()
#     contact::Contact = Contact()
#     damper::D = SimpleDamper()
#     steering::S = NoSteering()
#     braking::B = NoBraking()
# end

# struct LandingGearLegY{ContactY, SteeringY, BrakingY}
#     contact::ContactY
#     steering::SteeringY
#     braking::BrakingY
# end

# get_x0(ldg::LandingGearLeg) = ComponentVector(
#     contact = get_x0(ldg.contact),
#     steering = get_x0(ldg.steering),
#     braking = get_x0(ldg.braking)
#     )

# get_y0(ldg::LandingGearLeg) = LandingGearLegY(
#     get_y0(ldg.contact),
#     get_y0(ldg.steering),
#     get_y0(ldg.braking)
#     )

# get_u0(ldg::LandingGearLeg) = (
#     contact = get_u0(ldg.contact),
#     steering = get_u0(ldg.steering),
#     braking = get_u0(ldg.braking)
#     )

# get_d0(ldg::LandingGearLeg) = (
#     contact = get_d0(ldg.contact),
#     steering = get_d0(ldg.steering),
#     braking = get_d0(ldg.braking)
#     )


# function System(ldg::LandingGearLeg, ẋ = get_x0(ldg), x = get_x0(ldg),
#                     y = get_y0(ldg), u = get_u0(ldg), d = get_d0(ldg), t = Ref(0.0))

#     ss_list = Vector{System}()
#     ss_labels = (:contact, :steering, :braking)
#     for label in ss_labels
#         push!(ss_list, System(map((λ)->getproperty(λ, label), (ldg, ẋ, x, y, u, d))..., t))
#     end

#     params = (frame = ldg.frame, damper = ldg.damper)
#     subsystems = NamedTuple{ss_labels}(ss_list)

#     System{map(typeof, (ldg, x, y, u, d, params, subsystems))...}(
#                          ẋ, x, y, u, d, t, params, subsystems)
# end

# function f_cont!(sys::System{<:LandingGearLeg}, kin::KinData, trn::AbstractTerrain)

#     @unpack ẋ, x, u, params, subsystems = sys
#     @unpack frame, damper = params
#     @unpack contact, steering, braking = subsystems
#     @unpack n_e, h_e, q_eb, q_nb = kin.pos
#     @unpack v_eOb_b, ω_eb_b = kin.vel

#     #third basis element, for convenience
#     e3 = SVector{3,Float64}(0,0,1)

#     #update steering and braking subsystems
#     f_cont!(steering)
#     f_cont!(braking)

#     #strut frame axes
#     q_bs = frame.q_bc
#     q_ns = q_nb ∘ q_bs

#     #strut frame origin
#     r_ObOs_b = frame.r_ObOc_b
#     r_ObOs_e = q_eb(r_ObOs_b)
#     r_OeOb_e = CartECEF(Geographic(n_e, h_e))
#     Os = r_OeOb_e + r_ObOs_e

#     ###### WoW computation ######

#     #retrieve terrain data at the strut frame origin (close enough to the
#     #contact frame origin, which is still unknown)
#     terrain_data = get_terrain_data(trn, Os)
#     Δh = AltOrth(Os) - terrain_data.h

#     #project k_s onto z_n
#     k_s_zn = q_ns(e3)[3]

#     #if close to zero (strut nearly horizontal) or negative (strut upside down)
#     #we should set l=l_0 to yield WoW(false)
#     l = (k_s_zn > 1e-3 ? Δh / k_s_zn : damper.l_0)
#     ξ = l - damper.l_0
#     wow = WoW(ξ < 0)

#     if wow === WoW(false)
#         f_cont!(contact, WoW(false))
#         sys.y = LandingGearLegY(contact.y, steering.y, braking.y)
#         return
#     end

#     ###### contact frame construction ######

#     #contact frame axes
#     ψ_sw = get_steering_angle(steering)
#     q_sw = Rz(ψ_sw) #rotate strut axes to get wheel axes
#     q_nw = q_ns ∘ q_sw #NED to contact axes rotation
#     i_w_n = q_nw(e3) #NED components of wheel x-axis

#     k_t_n = terrain_data.k_t_n #terrain (inward pointing) normal
#     i_w_n_trn = i_w_n - (i_w_n ⋅ k_t_n) * k_t_n #projection of i_w_n onto the terrain tangent plane
#     i_c_n = normalize(i_w_n_trn) #NED components of contact x-axis
#     k_c_n = k_t_n #NED components of contact z-axis
#     j_c_n = k_c_n × i_c_n #NED components of contact y-axis
#     R_nc = RMatrix(SMatrix{3,3}([i_c_n j_c_n k_c_n]), normalization = false)
#     q_bc = q_nb' ∘ R_nc #airframe to contact axes rotation

#     #contact frame origin
#     r_OsOc_s = e3 * l
#     r_OsOc_b = q_bs(r_OsOc_s)
#     r_ObOc_b = r_OsOc_b + r_ObOs_b

#     bc = FrameSpec(r_ObOc_b, q_bc)

#     ###### contact point velocity computation ######

#     # contact frame origin velocity due to airframe motion
#     v_eOc_afm_b = v_eOb_b + ω_eb_b × r_ObOc_b
#     v_eOc_afm_c = q_bc'(v_eOc_af_b)

#     #compute the damper elongation rate required to cancel the airframe
#     #contribution to the contact point velocity along the contact frame z axis
#     q_cs = q_bc' ∘ q_bs
#     k_s_c = q_cs(e3)
#     ξ_dot = -v_eOc_afm_c[3] / k_s_c[3]
#     v_eOc_dmp_c = k_s_c * ξ_dot

#     v_eOc_c = v_eOc_afm_c + v_eOc_dmp_c

#     ###### contact force computation #####
#     η_br = get_braking_coefficient(braking)

#     #update contact model
#     f_cont!()

#     #we MUST NOT access contact.y as sys.y.contact, because the reference is not
#     #kept as it is with component vectors, and sys.y is only updated before
#     #returning from f_cont!
#     #in fact, we should only use this method:
#     get_normalized_contact_force()

#     #CONSIDERAR USAR RMatrix para todo. asi en vez de obtener k_s_n
#     #transformando 001 puedo coger directamente la 3a columna de R_ns

#     """
#     SET NORMALIZATION TO false!!!!!!!
#     """


#     #with the kinematics and the strut frame, determine WoW state

#     #the contact model gives us the normalized contact force in the contact
#     #frame


#     #now, from the kinematics of the contact point and the oleo parameters, we
#     #get F_shock_zs

#     sys.y = LandingGearLegY(contact.y, steering.y, braking.y)

# end

# #add f_disc! to reset the integrators for contact model
# f_disc!(::System{<:LandingGearLeg}) = false

# #CHANGE THIS
# #CHANGE THIS
# #CHANGE THIS
# #CHANGE THIS
# get_wr_b(::System{<:LandingGearLeg}) = Wrench()
# get_hr_b(::System{<:LandingGearLeg}) = zeros(SVector{3})


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