module Airdata

using LinearAlgebra
using StaticArrays
using UnPack

using Flight.Kinematics
using Flight.Geodesy
using Flight.Atmosphere
import Flight.Atmosphere: γ, p_std, ρ_std

export AirData

Base.@kwdef struct AirData
    v_wOb_b::SVector{3,Float64} = zeros(SVector{3}) #wind-relative velocity #v_wOb_b
    α_b::Float64 = 0.0 #airframe-axes AoA
    β_b::Float64 = 0.0 #airframe-axes AoS
    T::Float64 = 0.0 #static temperature
    p::Float64 = 0.0 #static pressure
    ρ::Float64 = 0.0 #density
    a::Float64 = 0.0 #speed of sound
    M::Float64 = 0.0 #Mach number
    Tt::Float64 = 0.0 #total temperature
    pt::Float64 = 0.0 #total pressure
    Δp::Float64 = 0.0 #impact pressure
    q::Float64 = 0.0 #dynamic pressure
    TAS::Float64 = 0.0 #true airspeed
    EAS::Float64 = 0.0 #equivalent airspeed
    CAS::Float64 = 0.0 #calibrated airspeed
end

function AirData(atm_sys::AtmosphericSystem, kin::KinY)
    #the AtmosphericData constructor accepts any Geographic subtype, but it's
    #likely that ISA SL conditions and wind will be expressed in LatLon
    atm_data = AtmosphericData(atm_sys, Geographic(kin.pos.ϕ_λ, kin.pos.h_o))
    AirData(atm_data, kin)
end

function AirData(atm::AtmosphericData, kin::KinY)

    v_eOb_b = kin.vel.v_eOb_b
    v_ew_b = kin.pos.q_nb'(atm.v_ew_n)
    v_wOb_b = v_eOb_b - v_ew_b
    (α_b, β_b) = get_airflow_angles(v_wOb_b)

    @unpack T, p, ρ, a = atm.isa_
    TAS = norm(v_wOb_b)
    M = TAS / a
    Tt = T * (1 + (γ - 1)/2 * M^2)
    pt = p * (Tt/T)^(γ/(γ-1))
    Δp = pt - p
    q = 1/2 * ρ * TAS^2

    EAS = TAS * √(ρ / ρ_std)
    CAS = √(2γ/(γ-1) * p_std/ρ_std * ( (1 + q/p_std)^((γ-1)/γ) - 1) )

    AirData(; v_wOb_b, α_b, β_b, T, p, ρ, a, M, Tt, pt, Δp, q, TAS, EAS, CAS)

end

#compute airflow angles at frame c from the c-frame aerodynamic velocity
function get_airflow_angles(v_wOc_c::AbstractVector{<:Real})
    (α_c = atan(v_wOc_c[3], v_wOc_c[1]),
    β_c = atan(v_wOc_c[2], √(v_wOc_c[1]^2 + v_wOc_c[3]^2)))
end

# having v_wOb_b in AirData, any AbstractComponent to which AirData is passed is
# free to transform v_wOb_b into its own frame. this may be particularly useful
# for an Aerodynamics component. if its frame is f(Oa, εa):

# v_wOb_b = v_eOb_b - v_ew_b
# v_eOa_b = v_eOb_b + ω_eb_b × r_ObOa_b
# v_wOa_b = v_eOa_b - v_ew_b = v_eOb_b + ω_eb_b × r_ObOa_b - v_ew_b
# v_wOa_b = v_wOb_b + ω_eb_b × r_ObOa_b
# v_wOa_a = q_ba'(v_wOa_b)

#generally, v_wOa ≈ v_wOb, so the whole conversion won't be necessary; at most,
#we may need to reproject v_wOb_b into v_wOb_a to comply with the axes used by
#the aerodynamics database

end