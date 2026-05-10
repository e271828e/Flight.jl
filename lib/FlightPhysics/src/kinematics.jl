module Kinematics

using StaticArrays, StructArrays, ComponentArrays, LinearAlgebra

using FlightCore
using FlightCore.GUI

using ..Attitude
using ..Geodesy

export AbstractKinematicDescriptor, ECEF, WA, NED
export KinInit, KinData, KinSystem

const v_min_χγ = 0.1 #minimum speed for valid χ, γ


############################## Initialization ##################################
################################################################################

#user-friendly kinematic conditions for body frame initialization
struct Initializer
    q_nb::RQuat #attitude with respect to NED frame
    n_e::NVector #2D location, n-vector
    h_e::Altitude{Ellipsoidal} #ellipsoidal altitude
    ω_wb_b::SVector{3, Float64} #angular velocity with respect to local level frame, body coordinates
    v_eb_n::SVector{3, Float64} #Earth-relative velocity, NED coordinates
end

const KinInit = Initializer

function Initializer(;
    q_nb::Abstract3DRotation = RQuat(), location::Abstract2DLocation = LatLon(),
    h::Altitude = HEllip(), ω_wb_b::AbstractVector{<:Real} = zeros(SVector{3}),
    v_eb_n::AbstractVector{<:Real} = zeros(SVector{3}))

    n_e = NVector(location)
    h_e = HEllip(h, n_e)

    Initializer(q_nb, n_e, h_e, ω_wb_b, v_eb_n)
end


################################### KinData ####################################
################################################################################

struct KinData
    e_nb::REuler #body frame attitude with respect to NED frame, Euler angles
    q_nb::RQuat #body frame attitude with respect to NED frame, quaternion
    q_eb::RQuat #body frame attitude with respect to ECEF frame, quaternion
    q_en::RQuat #NED frame attitude with respect to ECEF frame, quaternion
    ϕ_λ::LatLon #2D location, latitude / longitude
    n_e::NVector #2D location, n-vector
    h_e::Altitude{Ellipsoidal} #ellipsoidal altitude
    h_o::Altitude{Orthometric} #orthometric altitude
    r_eb_e::SVector{3,Float64} #Cartesian ECEF position
    ω_wb_b::SVector{3,Float64} #angular velocity with respect to local level frame, body coordinates
    ω_eb_b::SVector{3,Float64} #angular velocity with respect to ECEF frame, body coordinates
    v_eb_b::SVector{3,Float64} #Earth-relative velocity, body coordinates
    v_eb_n::SVector{3,Float64} #Earth-relative velocity, NED coordinates
    v_gnd::Float64 #Earth-relative speed
    χ_gnd::Float64 #Earth-relative course angle
    γ_gnd::Float64 #Earth-relative flight path angle
end

function KinData(ic::KinInit = KinInit())

    (; q_nb, n_e, h_e, ω_wb_b, v_eb_n) = ic

    Ob = Geographic(n_e, h_e)
    e_nb = REuler(q_nb)
    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    ϕ_λ = LatLon(n_e)
    h_o = HOrth(h_e, n_e)
    r_eb_e = Cartesian(Ob)

    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_eb_b = ω_ew_b + ω_wb_b

    v_eb_b = q_nb'(v_eb_n)

    v_gnd = norm(v_eb_n)
    χ_gnd = v_gnd > v_min_χγ ? azimuth(v_eb_n) : 0.0
    γ_gnd = v_gnd > v_min_χγ ? inclination(v_eb_n) : 0.0

    KinData(   e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, r_eb_e,
                ω_wb_b, ω_eb_b, v_eb_b, v_eb_n, v_gnd, χ_gnd, γ_gnd)

end

Base.getproperty(data::KinData, s::Symbol) = getproperty(data, Val(s))

@generated function Base.getproperty(data::KinData, ::Val{S}) where {S}
    if S ∈ fieldnames(KinData)
        return :(getfield(data, $(QuoteNode(S))))
    elseif S === :psi || S === :ψ
        return :(getfield(data, :e_nb).ψ)
    elseif S === :theta || S === :θ
        return :(getfield(data, :e_nb).θ)
    elseif S === :phi || S === :φ
        return :(getfield(data, :e_nb).φ)
    elseif S === :lat || S === :ϕ
        return :(getfield(data, :ϕ_λ).ϕ)
    elseif S === :lon || S === :λ
        return :(getfield(data, :ϕ_λ).λ)
    else
        error("$(typeof(data)) has no property $S")
    end
end


function normalize_block!(x, ε)
    norm_x = norm(x)
    (abs(norm_x - 1.0) > ε) && (x ./= norm_x)
    return nothing
end

function Geodesy.gravity(kin_data::KinData)
    gravity(Geographic(kin_data.n_e, kin_data.h_e))
end

function Geodesy.G_n(kin_data::KinData)
    G_n(Geographic(kin_data.n_e, kin_data.h_e))
end


######################### AbstractKinematicDescriptor ##########################
################################################################################

abstract type AbstractKinematicDescriptor <: ModelDefinition end
@no_periodic AbstractKinematicDescriptor

const XVelTemplate = ComponentVector(ω_eb_b = zeros(3), v_eb_b = zeros(3))

Modeling.U(::AbstractKinematicDescriptor) = zero(XVelTemplate)
Modeling.Y(::AbstractKinematicDescriptor) = KinData()

const KinSystem = Model{<:AbstractKinematicDescriptor}

KinData(mdl::KinSystem) = mdl.y

########################### WA-based Kinematics #########################
##########################################################################

#fast, singularity-free (all-attitude, all-latitude) kinematic mechanization,
#appropriate for simulation. WA = wander-azimuth frame

struct WA <: AbstractKinematicDescriptor end

Modeling.X(::WA) = ComponentVector(
    q_wb = zeros(4), q_ew = zeros(4), h_e = 0.0)

function Modeling.f_init!(mdl::Model{WA}, ic::Initializer = Initializer())

    (; x, u) = mdl
    (; q_nb, n_e, h_e, ω_wb_b, v_eb_n) = ic

    Ob = Geographic(n_e, h_e)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_eb_b = ω_ew_b + ω_wb_b
    v_eb_b = q_nb'(v_eb_n)
    h_e = HEllip(Ob)

    q_wb = q_nb #arbitrarily initializes wander angle ψ_nw to 1

    u.ω_eb_b = ω_eb_b
    u.v_eb_b = v_eb_b

    x.q_wb = q_wb[:]
    x.q_ew = ltf(n_e)[:]
    x.h_e = h_e

    f_ode!(mdl) #update ẋ and y, not strictly necessary

end


function Modeling.f_ode!(mdl::Model{WA})

    (; ẋ, x, u) = mdl

    q_wb = RQuat(x.q_wb, normalization = false)
    q_ew = RQuat(x.q_ew, normalization = false)
    ω_eb_b = SVector{3}(u.ω_eb_b)
    v_eb_b = SVector{3}(u.v_eb_b)
    h_e = HEllip(x.h_e[1])

    ψ_nw = get_ψ_nw(q_ew)
    q_nw = Rz(ψ_nw)
    q_nb = q_nw ∘ q_wb
    q_eb = q_ew ∘ q_wb
    q_en = q_eb ∘ q_nb'
    e_nb = REuler(q_nb)

    n_e = NVector(q_ew)
    ϕ_λ = LatLon(n_e)
    h_o = HOrth(h_e, n_e)

    v_eb_n = q_nb(v_eb_b)
    Ob = Geographic(n_e, h_e)
    r_eb_e = Cartesian(Ob)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)

    ω_ew_w = q_nw'(ω_ew_n)
    ω_ew_b = q_wb'(ω_ew_w)
    ω_wb_b = ω_eb_b - ω_ew_b

    v_gnd = norm(v_eb_n)
    χ_gnd = v_gnd > v_min_χγ ? azimuth(v_eb_n) : 0.0
    γ_gnd = v_gnd > v_min_χγ ? inclination(v_eb_n) : 0.0

    ẋ.q_wb = Attitude.dt(q_wb, ω_wb_b)
    ẋ.q_ew = Attitude.dt(q_ew, ω_ew_w)

    ẋ.h_e = -v_eb_n[3]

    mdl.y = KinData( e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, r_eb_e,
                 ω_wb_b, ω_eb_b, v_eb_b, v_eb_n, v_gnd, χ_gnd, γ_gnd)

end


function Modeling.f_step!(mdl::Model{WA}, ε = 1e-8)
    normalize_block!(mdl.x.q_wb, ε)
    normalize_block!(mdl.x.q_ew, ε)
end


@inline function get_ω_ew_n(v_eb_n::AbstractVector{<:Real}, Ob::Geographic)

    (R_N, R_E) = radii(Ob)
    h_e = HEllip(Ob)

    return SVector{3}(
        v_eb_n[2] / (R_E + Float64(h_e)),
        -v_eb_n[1] / (R_N + Float64(h_e)),
        0.0)

end

########################## ECEF-based Kinematics #########################
##########################################################################

#fast, singularity-free (all-attitude, all-latitude) kinematic mechanization,
#appropriate for simulation.

struct ECEF <: AbstractKinematicDescriptor end

Modeling.X(::ECEF) = ComponentVector(
    q_eb = zeros(4), n_e = zeros(3), h_e = 0.0)

function Modeling.f_init!(mdl::Model{ECEF}, ic::Initializer = Initializer())

    (; x, u) = mdl
    (; q_nb, n_e, h_e, ω_wb_b, v_eb_n) = ic

    Ob = Geographic(n_e, h_e)

    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_eb_b = ω_ew_b + ω_wb_b
    v_eb_b = q_nb'(v_eb_n)

    u.ω_eb_b = ω_eb_b
    u.v_eb_b = v_eb_b

    x.q_eb = q_eb[:]
    x.n_e = n_e[:]
    x.h_e = h_e

    f_ode!(mdl) #update ẋ and y, not strictly necessary

end

function Modeling.f_ode!(mdl::Model{ECEF})

    (; ẋ, x, u) = mdl

    q_eb = RQuat(x.q_eb, normalization = false)
    n_e = NVector(x.n_e, normalization = false)
    ω_eb_b = SVector{3}(u.ω_eb_b)
    v_eb_b = SVector{3}(u.v_eb_b)
    h_e = HEllip(x.h_e[1])
    h_o = HOrth(h_e, n_e)
    ϕ_λ = LatLon(n_e)

    q_en = ltf(n_e)
    q_nb = q_en' ∘ q_eb
    e_nb = REuler(q_nb)

    Ob = Geographic(n_e, h_e)
    r_eb_e = Cartesian(Ob)
    v_eb_n = q_nb(v_eb_b)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_wb_b = ω_eb_b - ω_ew_b

    v_gnd = norm(v_eb_n)
    χ_gnd = v_gnd > v_min_χγ ? azimuth(v_eb_n) : 0.0
    γ_gnd = v_gnd > v_min_χγ ? inclination(v_eb_n) : 0.0

    ẋ.q_eb = Attitude.dt(q_eb, ω_eb_b)
    ẋ.n_e = q_en(ω_ew_n × SVector{3,Float64}(0,0,-1))
    ẋ.h_e = -v_eb_n[3]

    mdl.y = KinData(e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, r_eb_e,
                ω_wb_b, ω_eb_b, v_eb_b, v_eb_n, v_gnd, χ_gnd, γ_gnd)

end

function Modeling.f_step!(mdl::Model{ECEF}, ε = 1e-8)
    normalize_block!(mdl.x.q_eb, ε)
    normalize_block!(mdl.x.n_e, ε)
end


################################ NED Kinematics ################################
################################################################################

#non singularity-free kinematic mechanization. useful mostly for aircraft
#control analysis and design

struct NED <: AbstractKinematicDescriptor end

Modeling.X(::NED) = ComponentVector(ψ_nb = 0.0, θ_nb = 0.0, φ_nb = 0.0,
                                ϕ = 0.0, λ = 0.0, h_e = 0.0)


function Modeling.f_init!(mdl::Model{NED}, ic::Initializer = Initializer())

    (; x, u) = mdl
    (; q_nb, n_e, h_e, ω_wb_b, v_eb_n) = ic

    Ob = Geographic(n_e, h_e)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_eb_b = ω_ew_b + ω_wb_b
    v_eb_b = q_nb'(v_eb_n)

    e_nb = REuler(q_nb)
    ϕ_λ = LatLon(n_e)

    u.ω_eb_b = ω_eb_b
    u.v_eb_b = v_eb_b

    x.ψ_nb = e_nb.ψ
    x.θ_nb = e_nb.θ
    x.φ_nb = e_nb.φ
    x.ϕ = ϕ_λ.ϕ
    x.λ = ϕ_λ.λ
    x.h_e = h_e

    f_ode!(mdl) #update ẋ and y, not strictly necessary

end

function Modeling.f_ode!(mdl::Model{NED})

    (; ẋ, x, u) = mdl

    e_nb = REuler(x.ψ_nb, x.θ_nb, x.φ_nb)
    ϕ_λ = LatLon(x.ϕ, x.λ)
    h_e = HEllip(x.h_e[1])
    ω_eb_b = SVector{3}(u.ω_eb_b)
    v_eb_b = SVector{3}(u.v_eb_b)

    n_e = NVector(ϕ_λ)
    h_o = HOrth(h_e, n_e)

    q_nb = RQuat(e_nb)
    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    v_eb_n = q_nb(v_eb_b)
    Ob = Geographic(n_e, h_e)
    r_eb_e = Cartesian(Ob)

    ω_en_n = get_ω_en_n(v_eb_n, Ob)
    ω_en_b = q_nb'(ω_en_n)
    ω_nb_b = ω_eb_b - ω_en_b

    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_wb_b = ω_eb_b - ω_ew_b

    v_gnd = norm(v_eb_n)
    χ_gnd = azimuth(v_eb_n)
    γ_gnd = inclination(v_eb_n)

    ė_nb = Attitude.dt(e_nb, ω_nb_b)
    ϕ_λ_dot = Geodesy.dt(ϕ_λ, ω_en_n)

    ẋ.ψ_nb = ė_nb.ψ
    ẋ.θ_nb = ė_nb.θ
    ẋ.φ_nb = ė_nb.φ
    ẋ.ϕ = ϕ_λ_dot.ϕ
    ẋ.λ = ϕ_λ_dot.λ
    ẋ.h_e = -v_eb_n[3]

    mdl.y = KinData( e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, r_eb_e,
                 ω_wb_b, ω_eb_b, v_eb_b, v_eb_n, v_gnd, χ_gnd, γ_gnd)

end

Modeling.f_step!(::Model{NED}) = nothing


@inline function get_ω_en_n(v_eb_n::AbstractVector{<:Real}, Ob::Geographic)

    (R_N, R_E) = radii(Ob)
    h_e = HEllip(Ob)
    ϕ = LatLon(Ob).ϕ

    return SVector{3}(
        v_eb_n[2] / (R_E + Float64(h_e)),
        -v_eb_n[1] / (R_N + Float64(h_e)),
        -v_eb_n[2] * tan(ϕ) / (R_E + Float64(h_e))
        )
end


################################################################################
################################# GUI ##########################################


GUI.draw(dyn::Model{<:AbstractKinematicDescriptor}) = GUI.draw(KinData(dyn))

function GUI.draw(data::KinData, p_open::Ref{Bool} = Ref(true),
                    label::String = "Kinematic Data")

    (; e_nb, ϕ_λ, h_e, h_o, ω_wb_b, ω_eb_b, v_eb_b, v_eb_n) = data

    BeginWindow(label, p_open)

    if TreeNode("Angular Velocity (Body / WA) [Body]")

        TextFormatted(@sprintf("Yaw Rate: %.7f deg/s", rad2deg(ω_wb_b[1])))
        TextFormatted(@sprintf("Pitch Rate: %.7f deg/s", rad2deg(ω_wb_b[2])))
        TextFormatted(@sprintf("Roll Rate: %.7f deg/s", rad2deg(ω_wb_b[3])))

        TreePop()
    end

    GUI.draw(rad2deg.(ω_eb_b - ω_wb_b), "Transport Rate (WA / ECEF) [Body]", "deg/s")
    GUI.draw(v_eb_n, "Velocity (Ob / ECEF) [NED]", "m/s")
    GUI.draw(v_eb_b, "Velocity (Ob / ECEF) [Body]", "m/s")

    if TreeNode("Attitude (Body / NED)")

        (; ψ, θ, φ) = e_nb
        TextFormatted(@sprintf("Heading: %.7f deg", rad2deg(ψ)))
        TextFormatted(@sprintf("Inclination: %.7f deg", rad2deg(θ)))
        TextFormatted(@sprintf("Bank: %.7f deg", rad2deg(φ)))

        TreePop()
    end

    if TreeNode("Position (O / ECEF)")

        (; ϕ, λ) = ϕ_λ
        TextFormatted(@sprintf("Latitude: %.7f deg", rad2deg(ϕ)))
        TextFormatted(@sprintf("Longitude: %.7f deg", rad2deg(λ)))
        TextFormatted(@sprintf("Altitude (Ellipsoidal): %.7f m", Float64(h_e)))
        TextFormatted(@sprintf("Altitude (Orthometric): %.7f m", Float64(h_o)))

        TreePop()
    end

    EndWindow()

end


################################################################################
############################### XPlane12Control ######################################

function Network.XPlanePose(kin_data::KinData)

    (; ϕ_λ, e_nb, h_o) = kin_data

    ϕ = rad2deg(ϕ_λ.ϕ)
    λ = rad2deg(ϕ_λ.λ)
    h = Float64(h_o)

    ψ = rad2deg(e_nb.ψ)
    θ = rad2deg(e_nb.θ)
    φ = rad2deg(e_nb.φ)

    Network.XPlanePose(ϕ, λ, h, ψ, θ, φ)

end

end #module
