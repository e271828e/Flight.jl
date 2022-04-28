module Air

using LinearAlgebra
using StaticArrays
using UnPack

using Flight.Systems
using Flight.Attitude
using Flight.Geodesy
using Flight.Atmosphere
using Flight.Kinematics

using Flight.Atmosphere: T_std, p_std, ρ_std, γ

export get_airflow_angles, get_wind_axes, get_stability_axes
export AirData

#compute airflow angles at frame c from the c-frame aerodynamic velocity
@inline function get_airflow_angles(v_wOc_c::AbstractVector{<:Real})::Tuple{Float64, Float64}
    α = atan(v_wOc_c[3], v_wOc_c[1])
    β = atan(v_wOc_c[2], √(v_wOc_c[1]^2 + v_wOc_c[3]^2))
    return (α, β)
end

@inline function get_wind_axes(v_wOc_c::AbstractVector{<:Real})
    α, β = get_airflow_angles(v_wOc_c)
    get_wind_axes(α, β)
end

@inline function get_wind_axes(α::Real, β::Real)
    q_bw = Ry(-α) ∘ Rz(β)
    return q_bw
end

@inline function get_stability_axes(α::Real)
    q_bs = Ry(-α)
    return q_bs
end

struct AirData
    v_ew_n::SVector{3,Float64} #wind velocity, NED axes
    v_ew_b::SVector{3,Float64} #wind velocity, vehicle axes
    v_eOb_b::SVector{3,Float64} #vehicle velocity vector
    v_wOb_b::SVector{3,Float64} #vehicle aerodynamic velocity vector
    T::Float64 #static temperature
    p::Float64 #static pressure
    ρ::Float64 #density
    a::Float64 #speed of sound
    μ::Float64 #dynamic viscosity
    M::Float64 #Mach number
    Tt::Float64 #total temperature
    pt::Float64 #total pressure
    Δp::Float64 #impact pressure
    q::Float64 #dynamic pressure
    TAS::Float64 #true airspeed
    EAS::Float64 #equivalent airspeed
    CAS::Float64 #calibrated airspeed
end

AirData() = AirData(KinData(), AtmosphericData())

function AirData(kin_data::KinData, atm_sys::AtmosphericSystem)
    #the AtmosphericData constructor accepts any Geographic subtype, but it's
    #likely that ISA SL conditions and wind will be expressed in LatLon
    atm_data = AtmosphericData(atm_sys, Geographic(kin_data.pos.ϕ_λ, kin_data.pos.h_o))
    AirData(kin_data, atm_data)
end

function AirData(kin_data::KinData, atm_data::AtmosphericData)

    v_eOb_b = kin_data.vel.v_eOb_b
    v_ew_n = atm_data.wind.v_ew_n
    v_ew_b = kin_data.pos.q_nb'(v_ew_n)
    v_wOb_b = v_eOb_b - v_ew_b

    @unpack T, p, ρ, a, μ = atm_data.static
    TAS = norm(v_wOb_b)
    M = TAS / a
    Tt = T * (1 + (γ - 1)/2 * M^2)
    pt = p * (Tt/T)^(γ/(γ-1))
    Δp = pt - p
    q = 1/2 * ρ * TAS^2

    EAS = TAS * √(ρ / ρ_std)
    CAS = √(2γ/(γ-1) * p_std/ρ_std * ( (1 + q/p_std)^((γ-1)/γ) - 1) )

    AirData(v_ew_n, v_ew_b, v_eOb_b, v_wOb_b, T, p, ρ, a, μ, M, Tt, pt, Δp, q, TAS, EAS, CAS)

end

end #module