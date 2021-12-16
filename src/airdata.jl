module Airdata

using LinearAlgebra
using StaticArrays
using UnPack

using Flight.Plotting
using Flight.Attitude

using Flight.Geodesy
using Flight.Atmosphere
using Flight.Kinematics

using Flight.Atmosphere: T_std, p_std, ρ_std, γ

import Flight.Plotting: plots

export get_airflow_angles, get_wind_axes, get_stability_axes
export AirData

#compute airflow angles at frame c from the c-frame aerodynamic velocity
@inline function get_airflow_angles(v_wOc_c::AbstractVector{<:Real})::Tuple{Float64, Float64}
    #let the aerodynamics handle this
    # if norm(v_wOc_c) < 1
    #     α = β = 0.0
    # else
        α = atan(v_wOc_c[3], v_wOc_c[1])
        β = atan(v_wOc_c[2], √(v_wOc_c[1]^2 + v_wOc_c[3]^2))
    # end
    return (α, β)
end

@inline function get_wind_axes(v_wOc_c::AbstractVector{<:Real})
    α, β = get_airflow_angles(v_wOc_c)
    get_wind_axes(α, β)
end

@inline function get_wind_axes(α::Real, β::Real)
    q_cw = Ry(-α) ∘ Rz(β)
    return q_cw
end

@inline function get_stability_axes(α::Real)
    q_cs = Ry(-α)
    return q_cs
end

struct AirData
    v_ew_n::SVector{3,Float64} #wind velocity, NED axes
    v_ew_b::SVector{3,Float64} #wind-relative, airframe axes
    v_eOb_b::SVector{3,Float64} #airframe velocity vector
    v_wOb_b::SVector{3,Float64} #airframe aerodynamic velocity vector
    α_b::Float64 #airframe angle of attack
    β_b::Float64 #airframe angle of sideslip
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

function AirData()
    AirData(KinData(), AtmosphericData())
end

function AirData(kin::KinData, atm_sys::AtmosphericSystem)
    #the AtmosphericData constructor accepts any Geographic subtype, but it's
    #likely that ISA SL conditions and wind will be expressed in LatLon
    atm_data = AtmosphericData(atm_sys, Geographic(kin.pos.ϕ_λ, kin.pos.h_o))
    AirData(kin, atm_data)
end

function AirData(kin::KinData, atm_data::AtmosphericData)

    v_eOb_b = kin.vel.v_eOb_b
    v_ew_n = atm_data.wind.v_ew_n
    v_ew_b = kin.pos.q_nb'(v_ew_n)
    v_wOb_b = v_eOb_b - v_ew_b
    α_b, β_b = get_airflow_angles(v_wOb_b)

    @unpack T, p, ρ, a, μ = atm_data.isa_
    TAS = norm(v_wOb_b)
    M = TAS / a
    Tt = T * (1 + (γ - 1)/2 * M^2)
    pt = p * (Tt/T)^(γ/(γ-1))
    Δp = pt - p
    q = 1/2 * ρ * TAS^2

    EAS = TAS * √(ρ / ρ_std)
    CAS = √(2γ/(γ-1) * p_std/ρ_std * ( (1 + q/p_std)^((γ-1)/γ) - 1) )

    AirData(v_ew_n, v_ew_b, v_eOb_b, v_wOb_b, α_b, β_b,
            T, p, ρ, a, μ, M, Tt, pt, Δp, q, TAS, EAS, CAS)

end


function plots(t, data::AbstractVector{<:AirData}; mode, save_path, kwargs...)

    @unpack v_ew_n, v_ew_b, v_eOb_b, v_wOb_b, α_b, β_b, a, μ, ρ, TAS, EAS, CAS, M,
            p, pt, T, Tt, Δp, q = StructArray(data)

    pd = Dict{String, Plots.Plot}()

    pd["01_v_ew_n"] = thplot(t, v_ew_n;
        plot_title = "Velocity (Wind / ECEF) [NED]",
        label = ["North" "East" "Down"],
        ylabel = [L"$v_{ew}^{N} \ (m/s)$" L"$v_{ew}^{E} \ (m/s)$" L"$v_{ew}^{D} \ (m/s)$"],
        th_split = :h,
        kwargs...)

    pd["02_v_ew_b"] = thplot(t, v_ew_b;
        plot_title = "Velocity (Wind / ECEF) [Airframe]",
        ylabel = [L"$v_{ew}^{x_b} \ (m/s)$" L"$v_{ew}^{y_b} \ (m/s)$" L"$v_{ew}^{z_b} \ (m/s)$"],
        th_split = :h,
        kwargs...)

    pd["03_v_eOb_b"] = thplot(t, v_eOb_b;
        plot_title = "Velocity (Airframe / ECEF) [Airframe]",
        ylabel = [L"$v_{eb}^{x_b} \ (m/s)$" L"$v_{eb}^{y_b} \ (m/s)$" L"$v_{eb}^{z_b} \ (m/s)$"],
        th_split = :h,
        kwargs...)

    pd["04_v_wOb_b"] = thplot(t, v_wOb_b;
        plot_title = "Velocity (Airframe / Wind) [Airframe]",
        ylabel = [L"$v_{eb}^{x_b} \ (m/s)$" L"$v_{eb}^{y_b} \ (m/s)$" L"$v_{eb}^{z_b} \ (m/s)$"],
        th_split = :h,
        kwargs...)

    splt_α = thplot(t, rad2deg.(α_b);
        title = "Angle of Attack", ylabel = L"$\alpha \ (deg)$",
        label = "", kwargs...)

    splt_β = thplot(t, rad2deg.(β_b);
        title = "Angle of Sideslip", ylabel = L"$\beta \ (deg)$",
        label = "", kwargs...)

    pd["05_α_β"] = plot(splt_α, splt_β;
        plot_title = "Airflow Angles [Airframe]",
        layout = (1,2),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    splt_a = thplot(t, a;
        title = "Speed of Sound", ylabel = L"$a \ (m/s)$",
        label = "", kwargs...)

    splt_ρ = thplot(t, ρ;
        title = "Density", ylabel = L"$\rho \ (kg/m^3)$",
        label = "", kwargs...)

    splt_μ = thplot(t, μ;
        title = "Dynamic Viscosity", ylabel = L"$\mu \ (Pa \ s)$",
        label = "", kwargs...)

    pd["06_ρ_a"] = plot(splt_ρ, splt_a, splt_μ;
        plot_title = "Freestream Properties",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    splt_T = thplot(t, hcat(T, Tt);
        title = "Temperature",
        label = ["Static"  "Total"],
        ylabel = L"$T \ (K)$",
        th_split = :none, kwargs...)

    splt_p = thplot(t, hcat(p, pt)/1000;
        title = "Pressure",
        label = ["Static"  "Total"],
        ylabel = L"$p \ (kPa)$",
        th_split = :none, kwargs...)

    pd["07_T_p"] = plot(splt_T, splt_p;
        plot_title = "Freestream Properties",
        layout = (1,2),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    splt_airspeed = thplot(t, hcat(TAS,EAS,CAS);
        title = "Airspeed",
        label = ["True" "Equivalent" "Calibrated"],
        ylabel = L"$v \ (m/s)$",
        th_split = :none, kwargs...)

    splt_Mach = thplot(t, M;
        title = "Mach", ylabel = L"M",
        label = "", kwargs...)

    splt_q = thplot(t, q/1000;
        title = "Dynamic Pressure", ylabel = L"$q \ (kPa)$",
        label = "", kwargs...)

    l3 = @layout [a{0.5w} [b; c{0.5h}]]
    pd["08_airspeed_M_q"] = plot(splt_airspeed, splt_Mach, splt_q;
        layout = l3,
        plot_title = "Freestream Properties",
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    save_plots(pd; save_path)

end

end #module