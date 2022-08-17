module LandingGear

using StaticArrays, ComponentArrays, LinearAlgebra, UnPack, Plots

using Flight.Utils
using Flight.Systems
using Flight.Plotting

using Flight.Attitude
using Flight.Geodesy
using Flight.Terrain
using Flight.Kinematics
using Flight.RigidBody
using Flight.Common

import Flight.Systems: init, f_ode!, f_step!
import Flight.RigidBody: MassTrait, WrenchTrait, AngularMomentumTrait, get_wr_b
import Flight.Plotting: make_plots

export LandingGearUnit, Strut, SimpleDamper, NoSteering, NoBraking, DirectSteering, DirectBraking

#basis elements, for convenience
const e1 = SVector{3,Float64}(1,0,0)
const e3 = SVector{3,Float64}(0,0,1)


################################################################################
################################# Steering #####################################

abstract type AbstractSteering <: Component end

################################ NoSteering ####################################

struct NoSteering <: AbstractSteering end

get_steering_angle(::System{NoSteering}) = 0.0

############################### DirectSteering #################################

Base.@kwdef struct DirectSteering <: AbstractSteering
    ψ_max::Float64 = π/6
end

Base.@kwdef struct DirectSteeringY
    ψ::Float64 = 0.0
end
#the contents of u must be mutable
init(::SystemU, ::DirectSteering) = Ref(Ranged(0.0, -1, 1))
init(::SystemY, ::DirectSteering) = DirectSteeringY(0.0) #steering angle

function f_ode!(sys::System{DirectSteering})
    sys.y = DirectSteeringY(Float64(sys.u[]) * sys.params.ψ_max)
end

get_steering_angle(sys::System{DirectSteering}) = sys.y.ψ


################################################################################
################################# Braking ######################################

abstract type AbstractBraking <: Component end


############################### NoBraking ######################################

struct NoBraking <: AbstractBraking end

get_braking_factor(::System{NoBraking}) = 0.0


############################# DirectBraking ####################################

Base.@kwdef struct DirectBraking <: AbstractBraking
   η_br::Float64 = 1.0 #braking efficiency
end

Base.@kwdef struct DirectBrakingY
   κ_br::Float64 = 0.0 #braking coefficient
end

init(::SystemU, ::DirectBraking) = Ref(Ranged(0.0, 0, 1))
init(::SystemY, ::DirectBraking) = DirectBrakingY()

function f_ode!(sys::System{DirectBraking})
    sys.y = DirectBrakingY(Float64(sys.u[]) * sys.params.η_br)
end

get_braking_factor(sys::System{DirectBraking}) = sys.y.κ_br


################################################################################
################################## Strut #######################################

################################### Damper #####################################

abstract type AbstractDamper end #not a System!

get_force(args...) = throw(MethodError(get_force, args))

Base.@kwdef struct SimpleDamper <: AbstractDamper
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
    @assert abs(F) < c.F_max "Maximum allowable damper force exceeded"
    return F
end

################################## Strut #######################################

Base.@kwdef struct Strut{D<:AbstractDamper} <: Component
    t_bs::FrameTransform = FrameTransform() #vehicle to strut frame transform
    l_0::Float64 = 0.0 #strut natural length from airframe attachment point to wheel endpoint
    damper::D = SimpleDamper()
end

Base.@kwdef struct StrutY #defaults should be consistent with wow = 0
    wow::Bool = false #weight-on-wheel flag
    ξ::Float64 = 0.0 #damper elongation
    ξ_dot::Float64 = 0.0 #damper elongation rate
    t_sc::FrameTransform = FrameTransform() #strut to contact frame transform
    t_bc::FrameTransform = FrameTransform() #body to contact frame transform
    F::Float64 = 0.0 #damper force on the ground along the strut z axis
    v_eOc_c::SVector{2,Float64} = zeros(SVector{2}) #contact point velocity
    trn::TerrainData = TerrainData()
end

init(::SystemY, ::Strut) = StrutY()

function f_ode!(sys::System{<:Strut}, steering::System{<:AbstractSteering},
    terrain::System{<:AbstractTerrain}, kin::KinematicData)

    @unpack t_bs, l_0, damper = sys.params
    @unpack q_eb, q_nb, q_en, n_e, h_e, r_eOb_e, v_eOb_b, ω_eb_b = kin

    q_bs = t_bs.q #body frame to strut frame rotation
    r_ObOs_b = t_bs.r #strut frame origin

    #do we have ground contact?
    q_es = q_eb ∘ q_bs
    ks_e = q_es(e3)
    r_ObOs_e = q_eb(r_ObOs_b)
    r_OsOw0_e = l_0 * ks_e
    r_eOw0_e = r_eOb_e + r_ObOs_e + r_OsOw0_e
    Ow0 = r_eOw0_e |> Cartesian |> Geographic

    loc_Ot = NVector(Ow0)
    trn_data_Ot = TerrainData(terrain, loc_Ot)
    h_Ot = HEllip(trn_data_Ot)
    Ot = Geographic(loc_Ot, h_Ot)
    r_eOt_e = Cartesian(Ot)[:]

    r_eOs_e = r_eOb_e + r_ObOs_e
    r_OtOs_e = r_eOs_e - r_eOt_e

    ut_n = trn_data_Ot.normal
    ut_e = q_en(ut_n)
    l = -(ut_e ⋅ r_OtOs_e) / (ut_e ⋅ ks_e)

    Δl = l - l_0
    if Δl > 0 #no contact, assign non-default outputs and return
        sys.y = StrutY(; wow = false, trn = trn_data_Ot)
        return
    end
    ξ = Δl

    #wheel frame axes
    ψ_sw = get_steering_angle(steering)
    q_sw = Rz(ψ_sw) #rotate strut axes to get wheel axes
    q_ns = q_nb ∘ q_bs
    q_nw = q_ns ∘ q_sw #NED to contact axes rotation

    #contact frame axes
    kc_n = trn_data_Ot.normal #NED components of contact z-axis
    iw_n = q_nw(e1) #NED components of wheel x-axis
    iw_n_trn = iw_n - (iw_n ⋅ kc_n) * kc_n #projection of iw_n onto the terrain tangent plane
    ic_n = normalize(iw_n_trn) #NED components of contact x-axis
    jc_n = kc_n × ic_n #NED components of contact y-axis
    R_nc = RMatrix(SMatrix{3,3}([ic_n jc_n kc_n]), normalization = false)
    q_sc = q_ns' ∘ R_nc
    q_bc = q_bs ∘ q_sc

    #contact frame origin
    r_OsOc_s = e3 * (l_0 + ξ)
    r_OsOc_b = q_bs(r_OsOc_s)
    r_ObOc_b = r_OsOc_b + r_ObOs_b

    #construct contact frame transforms
    t_sc = FrameTransform(r_OsOc_s, q_sc)
    t_bc = FrameTransform(r_ObOc_b, q_bc)

    #contact frame origin velocity due to rigid body vehicle motion
    v_eOc_veh_b = v_eOb_b + ω_eb_b × r_ObOc_b
    v_eOc_veh_c = q_bc'(v_eOc_veh_b)

    #compute the damper elongation rate required to cancel the vehicle
    #contribution to the contact point velocity along the contact frame z axis
    q_sc = q_bs' ∘ q_bc
    ks_c = q_sc'(e3)
    ξ_dot = -v_eOc_veh_c[3] / ks_c[3]

    #compute contact point velocity
    v_eOc_dmp_c = ks_c * ξ_dot
    v_eOc_c_3D = v_eOc_veh_c + v_eOc_dmp_c
    v_eOc_c = v_eOc_c_3D[SVector(1,2)]
    @assert abs(v_eOc_c_3D[3]) < 1e-8

    F = get_force(damper, ξ, ξ_dot)

    sys.y = StrutY(; wow = true, ξ, ξ_dot, t_sc, t_bc, F, v_eOc_c, trn = trn_data_Ot)

end



################################################################################
############################### Contact ########################################

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
    @unpack v_s, v_d, μ_s, μ_d = fr
    κ_sd = clamp((norm(v) - v_s) / (v_d - v_s), 0, 1)
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


############################### Contact ########################################

Base.@kwdef struct Contact <: Component
    frc::PICompensator{2} = PICompensator{2}( #friction constraint compensator
        k_p = 5.0, k_i = 400.0, k_l = 0.2, bounds = (-1.0, 1.0))
end

Base.@kwdef struct ContactY
    frc::Common.PICompensatorY{2} = Common.PICompensatorY{2}()
    μ_roll::Float64 = 0.0 #rolling friction coefficient
    μ_skid::Float64 = 0.0 #skidding friction coefficient
    κ_br::Float64 = 0.0 #braking factor
    ψ_cv::Float64 = 0.0 #tire slip angle
    μ_max::SVector{2,Float64} = zeros(SVector{2}) #maximum friction coefficient
    μ_eff::SVector{2,Float64} = zeros(SVector{2}) #scaled friction coefficient
    f_c::SVector{3,Float64} = zeros(SVector{3}) #normalized contact force
    F_c::SVector{3,Float64} = zeros(SVector{3}) #contact force
    wr_b::Wrench = Wrench() #resulting Wrench on the vehicle frame
end

#x should be initialized by the default methods
init(::SystemY, ::Contact) = ContactY()

function f_ode!(sys::System{Contact}, strut::System{<:Strut},
                braking::System{<:AbstractBraking})

    @unpack wow, t_sc, t_bc, F, v_eOc_c, trn = strut.y

    frc = sys.frc
    frc.u.sat_enable .= true #this should always be enabled (could also be set in init_u)
    frc.u.input .= v_eOc_c #if !wow, v_eOc_c = [0,0]
    frc.u.reset .= !wow #if !wow, state will reset on the next call to f_step!

    #if !wow, call f_ode!(frc) once to let the reset input propagate to frc's
    #output. once it's updated, we have no need to call f_ode! until wow = true
    if !wow
        all(frc.y.reset) ? nothing : (f_ode!(frc); sys.y = ContactY(; frc = frc.y))
        return #we are done
    end

    f_ode!(frc) #update friction constraint compensator

    norm_v = norm(v_eOc_c)

    μ_roll = get_μ(FrictionCoefficients(Rolling(), trn.surface), norm_v)
    μ_skid = get_μ(FrictionCoefficients(Skidding(), trn.surface), norm_v)

    #longitudinal friction coefficient
    κ_br = get_braking_factor(braking)
    μ_x = μ_roll + (μ_skid - μ_roll) * κ_br

    #tire slip angle
    if norm_v < 1e-3 #prevents chattering in μ_y for near-zero contact velocity
        ψ_cv = π/2 #pure sideslip
    else
        ψ_cv = atan(v_eOc_c[2], v_eOc_c[1])
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

    μ_max = @SVector [μ_x, μ_y]
    μ_max *= min(1, μ_skid / norm(μ_max)) #scale μ_max so norm(μ_max) does not exceed μ_skid

    #scale μ_max with the feedback from the friction constraint compensator
    μ_eff = -frc.y.out .* μ_max

    #normalized contact force projected on the contact frame
    f_c = SVector{3,Float64}(μ_eff[1], μ_eff[2], -1)
    f_s = t_sc.q(f_c) #project normalized force onto the strut frame
    # @assert f_s[3] < 0

    #the value of the ground's normal force must be such that its projection
    #along the strut cancels the damper's force
    N = -F / f_s[3]
    N = max(0, N) #it must not be negative (might happen with large ξ_dot >0)
    F_c = f_c * N

    wr_c = Wrench(F = F_c)
    wr_b = t_bc(wr_c)

    sys.y = ContactY(; frc = frc.y,
                       μ_roll, μ_skid, κ_br, ψ_cv, μ_max, μ_eff, f_c, F_c, wr_b)

end


################################################################################
############################ LandingGearUnit ###################################

Base.@kwdef struct LandingGearUnit{S <: AbstractSteering,
                            B <: AbstractBraking, L<:Strut} <: Component
    steering::S = NoSteering()
    braking::B = NoBraking()
    strut::L = Strut()
    contact::Contact = Contact()
end

#if we avoid the generic fallback for SystemGroup, we don't need to define
#traits for Steering, Braking, Contact or Strut
MassTrait(::System{<:LandingGearUnit}) = HasNoMass()
WrenchTrait(::System{<:LandingGearUnit}) = GetsExternalWrench()
AngularMomentumTrait(::System{<:LandingGearUnit}) = HasNoAngularMomentum()

get_wr_b(sys::System{<:LandingGearUnit}) = sys.y.contact.wr_b

function f_ode!(sys::System{<:LandingGearUnit}, kinematics::KinematicData,
                terrain::System{<:AbstractTerrain})

    @unpack strut, contact, steering, braking = sys

    f_ode!(steering)
    f_ode!(braking)
    f_ode!(strut, steering, terrain, kinematics)
    f_ode!(contact, strut, braking)

    Systems.update_y!(sys)

end

function f_step!(sys::System{<:LandingGearUnit})

    x_mod = false

    x_mod = x_mod || f_step!(sys.steering)
    x_mod = x_mod || f_step!(sys.braking)
    x_mod = x_mod || f_step!(sys.strut)
    x_mod = x_mod || f_step!(sys.contact)

    return x_mod

end


################################################################################
############################ Plotting ##########################################

function make_plots(th::TimeHistory{<:StrutY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    subplot_ξ = plot(th.ξ;
        title = "Elongation", ylabel = L"$\xi \ (m)$",
        label = "", kwargs...)

    subplot_ξ_dot = plot(th.ξ_dot;
        title = "Elongation Rate", ylabel = L"$\dot{\xi} \ (m/s)$",
        label = "", kwargs...)

    subplot_F = plot(th.F;
        title = "Force", ylabel = L"$F \ (N)$",
        label = "", kwargs...)

    pd[:dmp] = plot(subplot_ξ, subplot_ξ_dot, subplot_F;
        plot_title = "Damper",
        layout = (1,3), link = :none,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end

function make_plots(th::TimeHistory{<:ContactY}; kwargs...)

    pd = OrderedDict{Symbol, Any}()

    pd[:frc] = make_plots(th.frc; kwargs...)

    (μ_max_x, μ_max_y) = Utils.get_scalar_components(th.μ_max)
    (μ_eff_x, μ_eff_y) = Utils.get_scalar_components(th.μ_eff)

    subplot_μ_roll = plot(th.μ_roll; title = "Rolling Friction Coefficient",
        ylabel = L"$\mu_{roll}$", label = "", kwargs...)
    subplot_μ_skid = plot(th.μ_skid; title = "Skidding Friction Coefficient",
        ylabel = L"$\mu_{skid}$", label = "", kwargs...)

    pd[:srf] = plot(subplot_μ_roll, subplot_μ_skid;
        plot_title = "Surface Friction",
        layout = (1,2), link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    subplot_κ_br = plot(th.κ_br; title = "Braking Coefficient",
        ylabel = L"$\alpha_{br}$", label = "", kwargs...)
    subplot_μ_max_x = plot(μ_max_x; title = "Maximum Friction Coefficient",
        ylabel = L"$\mu_{max}^{x}$", label = "", kwargs...)
    subplot_μ_eff_x = plot(μ_eff_x; title = "Effective Friction Coefficient",
        ylabel = L"$\mu^{x}$", label = "", kwargs...)

    pd[:μ_x] = plot(subplot_κ_br, subplot_μ_max_x, subplot_μ_eff_x;
        plot_title = "Longitudinal Friction",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    subplot_ψ_cv = plot(th._t, rad2deg.(th.ψ_cv._data); title = "Tire Slip Angle",
        ylabel = L"$\psi_{cv} \ (deg)$", label = "", kwargs...)
    subplot_μ_max_y = plot(μ_max_y; title = "Maximum Friction Coefficient",
        ylabel = L"$\mu_{max}^{y}$", label = "", kwargs...)
    subplot_μ_eff_y = plot(μ_eff_y; title = "Effective Friction Coefficient",
        ylabel = L"$\mu^{y}$", label = "", kwargs...)

    pd[:μ_y] = plot(subplot_ψ_cv, subplot_μ_max_y, subplot_μ_eff_y;
        plot_title = "Lateral Friction",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd[:f_c] = plot(th.f_c;
        plot_title = "Normalized Contact Force",
        ylabel = [L"$f_{Oc \ (trn)}^{c}$" L"$f_{Oc \ (trn)}^{c}$" L"$f_{Oc \ (trn)}^{c}$"],
        th_split = :h, link = :none,
        kwargs...)

    pd[:F_c] = plot(th.F_c;
        plot_title = "Contact Force",
        ylabel = [L"$F_{Oc \ (trn)}^{c} \ (N)$" L"$F_{Oc \ (trn)}^{c} \ (N)$" L"$F_{Oc \ (trn)}^{c} \ (N)$"],
        th_split = :h, link = :none,
        kwargs...)

    pd[:wr_b] = plot(th.wr_b;
        plot_title = "Wrench [Vehicle Axes]",
        wr_source = "trn", wr_frame = "b",
        kwargs...)

    return pd

end


end #module