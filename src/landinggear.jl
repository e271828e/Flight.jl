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

import Flight.Modeling: init, f_cont!, f_disc!
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_wr_b
import Flight.Plotting: plots

export LandingGearUnit, Strut, SimpleDamper, NoSteering, NoBraking, DirectSteering, DirectBraking

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
#we need to make the contents of u mutable
init(::DirectSteering, ::SystemU) = Ref(0.0) #steering input
init(::DirectSteering, ::SystemY) = DirectSteeringY(0.0) #steering angle

function f_cont!(sys::System{DirectSteering})
    u = sys.u[]
    @assert abs(u) <= 1 "Steering input must be within [-1,1]"
    sys.y = DirectSteeringY(u * sys.params.ψ_max)
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
   η_br::Float64 = 0.0 #braking coefficient
end

init(::DirectBraking, ::SystemU) = Ref(0.0)
init(::DirectBraking, ::SystemY) = DirectBrakingY()

function f_cont!(sys::System{DirectBraking})
    u = sys.u[]
    @assert (u <= 1 && u >= 0) "Braking input must be within [0,1]"
    sys.y = DirectBrakingY(u * sys.params.η_br)
end

f_disc!(::System{DirectBraking}, args...) = false

get_braking_factor(sys::System{DirectBraking}) = sys.y.η_br


########################### Damper #############################

abstract type AbstractDamper end #not a System!
get_damper_force(args...) = throw(MethodError(get_damper_force, args))

Base.@kwdef struct SimpleDamper <: AbstractDamper
    k_s::Float64 = 25000 #spring constant
    k_d_ext::Float64 = 1000 #extension damping coefficient
    k_d_cmp::Float64 = 1000 #compression damping coefficient
    ξ_min::Float64 = -5 #compression below which the shock absorber is disabled
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

Base.@kwdef struct Strut{D<:AbstractDamper} <: SystemDescriptor
    t_bs::FrameTransform = FrameTransform()
    damper::D = SimpleDamper()
    l_0::Float64 = 0.0 #natural length
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

init(::Strut, ::SystemY) = StrutY()

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

    F_dmp = get_damper_force(damper, ξ, ξ_dot)

    sys.y = StrutY(; wow, ξ, ξ_dot, t_sc, t_bc, F_dmp, v_eOc_c, srf)

end

f_disc!(::System{<:Strut}, args...) = false

function plots(t, data::AbstractVector{<:StrutY}; mode, save_path, kwargs...)

    @unpack ξ, ξ_dot, F_dmp = StructArray(data)

    pd = Dict{String, Plots.Plot}()

    splt_ξ = thplot(t, ξ;
        title = "Elongation", ylabel = L"$\xi \ (m)$",
        label = "", kwargs...)

    splt_ξ_dot = thplot(t, ξ_dot;
        title = "Elongation Rate", ylabel = L"$\dot{\xi} \ (m/s)$",
        label = "", kwargs...)

    splt_F_dmp = thplot(t, F_dmp;
        title = "Force", ylabel = L"$F_{dmp} \ (N)$",
        label = "", kwargs...)

    pd["01_dmp"] = plot(splt_ξ, splt_ξ_dot, splt_F_dmp;
        plot_title = "Damper",
        layout = (1,3), link = :none,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    save_plots(pd; save_path)


end

########################### Contact #############################

#static / dynamic friction coefficient set
struct StaticDynamic
    static::Float64
    dynamic::Float64
end

#static / dynamic friction interpolation from breakout factor
get_μ(f::StaticDynamic, η_bo::Real)::Float64 = η_bo * f.dynamic + (1-η_bo) * f.static

#pure rolling friction coefficient set
Base.@kwdef struct Rolling
    sd::StaticDynamic= StaticDynamic(0.03, 0.02)
end

get_sd(f::Rolling, ::SurfaceType)::StaticDynamic = f.sd

#pure skidding friction coefficient set
Base.@kwdef struct Skidding
    dry::StaticDynamic = StaticDynamic(0.75, 0.25)
    wet::StaticDynamic = StaticDynamic(0.25, 0.15)
    icy::StaticDynamic = StaticDynamic(0.075, 0.025)
end

function get_sd(f::Skidding, srf::SurfaceType)::StaticDynamic
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

get_μ(f::Union{Rolling,Skidding}, srf::SurfaceType, η_bo::Real)::Float64 =
    get_μ(get_sd(f, srf), η_bo)

#why does this allocate??
# get_sd(f::Skidding, ::Val{Terrain.DryTarmac}) = f.dry
# get_sd(f::Skidding, ::Val{Terrain.WetTarmac}) = f.wet
# get_sd(f::Skidding, ::Val{Terrain.IcyTarmac}) = f.icy

#note: setting k_l > 0 asymptotically removes residual force fighting between
#contacts in equilibrium at the expense of some position creep under persistent
#external loads. k_l = 0.2 gives reasonable practical results

Base.@kwdef struct Contact <: SystemDescriptor
    rolling::Rolling = Rolling() #rolling friction coefficients
    skidding::Skidding = Skidding() #skidding friction coefficients
    ψ_skid::Float64 = deg2rad(10) #skidding slip angle thresold
    v_bo::SVector{2,Float64} = [0.005, 0.01] #breakout velocity interval
    k_p::Float64 = 5.0 #proportional gain for contact velocity regulator
    k_i::Float64 = 400.0 #integral gain for contact velocity regulator
    k_l::Float64 = 0.0 #integrator leak factor for contact velocity regulator
end

WrenchTrait(::System{<:Contact}) = GetsNoExternalWrench()

Base.@kwdef struct ContactY
    v::SVector{2,Float64} = zeros(SVector{2}) #contact plane velocity
    s::SVector{2,Float64} = zeros(SVector{2}) #contact plane velocity integral
    α_p::SVector{2,Float64} = zeros(SVector{2}) #proportional μ scale factor
    α_i::SVector{2,Float64} = zeros(SVector{2}) #integral μ scale factor
    α_raw::SVector{2,Float64} = zeros(SVector{2}) #total scale factor, raw
    α::SVector{2,Float64} = zeros(SVector{2}) #total scale factor, clipped
    sat::SVector{2,Bool} = zeros(SVector{2,Bool}) #scale factor saturation flag
    η_bo::Float64 = 0.0 #breakout factor
    μ_roll::Float64 = 0.0 #rolling friction coefficient
    μ_skid::Float64 = 0.0 #skidding friction coefficient
    η_br::Float64 = 0.0 #braking factor
    ψ_cv::Float64 = 0.0 #tire slip angle
    μ_max::SVector{2,Float64} = zeros(SVector{2}) #maximum friction coefficient
    μ::SVector{2,Float64} = zeros(SVector{2}) #scaled friction coefficient
    f_c::SVector{3,Float64} = zeros(SVector{3}) #normalized contact force
    F_c::SVector{3,Float64} = zeros(SVector{3}) #contact force
    wr_b::Wrench = Wrench() #resulting airframe Wrench
end

init(::Contact, ::SystemX) = ComponentVector(x = 0.0, y = 0.0) #v regulator integrator states
init(::Contact, ::SystemY) = ContactY()
init(::Contact, ::SystemD) = Ref(false) #contact active?

function f_cont!(sys::System{Contact}, strut::System{<:Strut},
                braking::System{<:AbstractBraking})


    @unpack rolling, skidding, v_bo, ψ_skid, k_p, k_i, k_l = sys.params
    @unpack wow, t_sc, t_bc, F_dmp, v_eOc_c, srf = strut.y

    if !wow #better than !sys.d[], which is set at f_disc! so has a one-step delay
        sys.ẋ .= 0
        sys.y = ContactY()
        return
    end

    v = SVector{2,Float64}(v_eOc_c[1], v_eOc_c[2]) #contact point velocity
    s = SVector{2,Float64}(sys.x) #velocity integrator state

    #breakout factor
    η_bo = clamp((norm(v) - v_bo[1]) / (v_bo[2] - v_bo[1]), 0, 1)

    μ_roll = get_μ(rolling, srf, η_bo)
    μ_skid = get_μ(skidding, srf, η_bo)

    #longitudinal friction coefficient
    η_br = get_braking_factor(braking)
    @assert (η_br >= 0 && η_br <= 1)
    μ_x = μ_roll + (μ_skid - μ_roll) * η_br

    #tire slip angle
    if η_bo < 1 #prevents chattering in μ_y for near-zero contact velocity
        ψ_cv = ψ_skid
    else
        ψ_cv = atan(v[2], v[1])
    end
    ψ_abs = abs(ψ_cv)
    @assert (ψ_abs <= π)

    #lateral friction coefficient
    if ψ_abs < ψ_skid
        μ_y = μ_skid * ψ_abs / ψ_skid
    elseif ψ_abs > π - ψ_skid
        μ_y = μ_skid * (1 - (ψ_skid + ψ_abs - π)/ ψ_skid)
    else
        μ_y = μ_skid
    end

    #maximum friction coefficient vector
    μ_max = @SVector [μ_x, μ_y]

    #the magnitude of μ_max cannot exceed μ_skid; if it does, scale down its components
    μ_max *= min(1, μ_skid / norm(μ_max))

    α_p = -k_p * v
    α_i = -k_i * s
    α_raw = α_p + α_i #raw μ scaling
    α = clamp.(α_raw, -1, 1) #clipped μ scaling
    #if not saturated, integrator accumulates
    sat = abs.(α_raw) .> abs.(α) #saturated?
    sys.ẋ .= (v - k_l * s) .* .!sat

    #normalized contact force projected on the contact frame
    μ = α .* μ_max
    f_c = SVector{3,Float64}(μ[1], μ[2], -1)
    f_s = t_sc.q(f_c) #project normalized force onto the strut frame
    @assert f_s[3] < 0

    N = -F_dmp / f_s[3]
    N = max(0, N) #cannot be negative (could happen with large xi_dot >0)
    F_c = f_c * N

    wr_c = Wrench(F = F_c)
    wr_b = t_bc(wr_c)

    sys.y = ContactY(; v, s, α_p, α_i, α_raw, α, sat, η_bo, μ_roll, μ_skid,
        η_br, ψ_cv, μ_max, μ, f_c, F_c, wr_b)

end

function f_disc!(sys::System{Contact}, strut::System{<:Strut})

    if strut.y.wow
        sys.d[] = true #declare contact as active
    else #no wow
        if sys.d[] #contact was active, declare it inactive and reset integrator states
            sys.d[] = false
            sys.x .= 0 #reset integrator states
            return true #continuous state was modified
        end
    end
    return false
end

function plots(t, data::AbstractVector{<:ContactY}; mode, save_path, kwargs...)

    @unpack v, s, α_p, α_i, α_raw, α, sat, η_bo, μ_roll, μ_skid,
        η_br, ψ_cv, μ_max, μ, f_c, F_c, wr_b = StructArray(data)

    extract_xy = (input) -> (input |> StructArray |> StructArrays.components)

    (v_x, v_y), (α_p_x, α_p_y), (α_i_x, α_i_y), (sat_x, sat_y),
        (α_raw_x, α_raw_y), (α_x, α_y), (μ_max_x, μ_max_y), (μ_x, μ_y) =
         map(extract_xy, (v, α_p, α_i, sat, α, α, μ_max, μ))

    pd = Dict{String, Plots.Plot}()

    splt_v_mag = thplot(t, norm.(v); title = "Contact Velocity Magnitude",
        ylabel = L"$\alpha_{bo}$", label = "", kwargs...)
    splt_η_bo = thplot(t, η_bo; title = "Breakout Factor",
        ylabel = L"$\alpha_{bo}$", label = "", kwargs...)
    splt_μ_roll = thplot(t, μ_roll; title = "Rolling Friction Coefficient",
        ylabel = L"$\mu_{roll}$", label = "", kwargs...)
    splt_μ_skid = thplot(t, μ_skid; title = "Skidding Friction Coefficient",
        ylabel = L"$\mu_{skid}$", label = "", kwargs...)

    pd["01_srf"] = plot(splt_v_mag, splt_η_bo, splt_μ_roll, splt_μ_skid;
        plot_title = "Rolling and Skidding Friction",
        layout = (2,2), link = :none,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    splt_v_x = thplot(t, v_x; title = "Contact Point Velocity",
        ylabel = L"$v^x \ (m/s)$", label = "", kwargs...)
    splt_α_p_x = thplot(t, α_p_x; title = "Proportional Term",
        ylabel = L"$\alpha_p^x$", label = "", kwargs...)
    splt_α_i_x = thplot(t, α_i_x; title = "Integral Term",
        ylabel = L"$\alpha_i^x$", label = "", kwargs...)
    splt_sat_x = thplot(t, sat_x; title = "Saturation",
        ylabel = L"$S^x$", label = "", kwargs...)
    splt_α_raw_x = thplot(t, α_raw_x; title = "Raw Output",
        ylabel = L"$\alpha_{raw}^{x}$", label = "", kwargs...)
    splt_α_x = thplot(t, α_x; title = "Clipped Output",
        ylabel = L"$\alpha^{x}$", label = "", kwargs...)

    pd["02_reg_x"] = plot(splt_v_x, splt_α_p_x, splt_α_i_x,
                        splt_sat_x, splt_α_raw_x, splt_α_x;
        plot_title = "Longitudinal Friction Regulator",
        layout = (2,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    splt_η_br = thplot(t, η_br; title = "Braking Coefficient",
        ylabel = L"$\alpha_{br}$", label = "", kwargs...)
    splt_μ_max_x = thplot(t, μ_max_x; title = "Maximum Friction Coefficient",
        ylabel = L"$\mu_{max}^{x}$", label = "", kwargs...)
    splt_μ_x = thplot(t, μ_x; title = "Effective Friction Coefficient",
        ylabel = L"$\mu^{x}$", label = "", kwargs...)

    pd["03_mu_x"] = plot(splt_η_br, splt_μ_max_x, splt_μ_x;
        plot_title = "Longitudinal Friction",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    splt_v_y = thplot(t, v_y; title = "Contact Point Velocity",
        ylabel = L"$v^y \ (m/s)$", label = "", kwargs...)
    splt_α_p_y = thplot(t, α_p_y; title = "Proportional Term",
        ylabel = L"$\alpha_p^y$", label = "", kwargs...)
    splt_α_i_y = thplot(t, α_i_y; title = "Integral Term",
        ylabel = L"$\alpha_i^y$", label = "", kwargs...)
    splt_sat_y = thplot(t, sat_y; title = "Saturation",
        ylabel = L"$S^y$", label = "", kwargs...)
    splt_α_raw_y = thplot(t, α_raw_y; title = "Raw Output",
        ylabel = L"$\alpha_{raw}^{y}$", label = "", kwargs...)
    splt_α_y = thplot(t, α_y; title = "Clipped Output",
        ylabel = L"$\alpha^{y}$", label = "", kwargs...)

    pd["04_reg_y"] = plot(splt_v_y, splt_α_p_y, splt_α_i_y,
                        splt_sat_y, splt_α_raw_y, splt_α_y;
        plot_title = "Lateral Contact Regulator",
        layout = (2,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    splt_ψ_cv = thplot(t, rad2deg.(ψ_cv); title = "Tire Slip Angle",
        ylabel = L"$\psi_{cv} \ (deg)$", label = "", kwargs...)
    splt_μ_max_y = thplot(t, μ_max_y; title = "Maximum Friction Coefficient",
        ylabel = L"$\mu_{max}^{y}$", label = "", kwargs...)
    splt_μ_y = thplot(t, μ_y; title = "Effective Friction Coefficient",
        ylabel = L"$\mu^{y}$", label = "", kwargs...)

    pd["05_mu_y"] = plot(splt_ψ_cv, splt_μ_max_y, splt_μ_y;
        plot_title = "Lateral Friction",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd["06_Fn_c"] = thplot(t, f_c;
        plot_title = "Normalized Contact Force",
        ylabel = [L"$f_{Oc \ (trn)}^{c}$" L"$f_{Oc \ (trn)}^{c}$" L"$f_{Oc \ (trn)}^{c}$"],
        th_split = :h, link = :none,
        kwargs...)

    pd["07_F_c"] = thplot(t, F_c;
        plot_title = "Contact Force",
        ylabel = [L"$F_{Oc \ (trn)}^{c} \ (N)$" L"$F_{Oc \ (trn)}^{c} \ (N)$" L"$F_{Oc \ (trn)}^{c} \ (N)$"],
        th_split = :h, link = :none,
        kwargs...)

    pd["08_wr_b"] = thplot(t, wr_b;
        plot_title = "Wrench [Airframe]",
        wr_source = "trn", wr_frame = "b",
        kwargs...)

    save_plots(pd; save_path)

end

########################## LandingGearUnit #########################

Base.@kwdef struct LandingGearUnit{L<:Strut, S <: AbstractSteering,
                            B <: AbstractBraking} <: NodeSystemDescriptor
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

    Modeling.update_y!(sys)

end

function f_disc!(sys::System{<:LandingGearUnit})

    @unpack strut, contact = sys.subsystems

    contact_modified = f_disc!(contact, strut)
    return contact_modified

end

#need to override the default AirframeNode implementation, since the underlying
#subcomponents are not AirframeComponents and only Contact contributes wr_b.
#could also define trivial get_wr_b and get_hr_b methods for steering, braking
#and strut, and let the default AirframeNode implementation do its thing
# get_wr_b(sys::System{<:LandingGearUnit}) = sys.y.contact.wr_b


end #module