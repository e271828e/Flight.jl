module Airdata

using LinearAlgebra
using StaticArrays
using UnPack

using Flight.Attitude
using Flight.Kinematics
using Flight.Geodesy
using Flight.Atmosphere
import Flight.Atmosphere: T_std, p_std, ρ_std, γ

using Flight.Plotting
import Flight.Plotting: plots

export AirData

struct AirData
    v_ew_n::SVector{3,Float64} #wind velocity, NED axes
    v_ew_b::SVector{3,Float64} #wind-relative, airframe axes
    v_eOb_b::SVector{3,Float64} #velocity vector, airframe axes
    v_wOb_b::SVector{3,Float64} #aerodynamic velocity vector, airframe axes
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

    @unpack T, p, ρ, a, μ = atm_data.isa_
    TAS = norm(v_wOb_b)
    M = TAS / a
    Tt = T * (1 + (γ - 1)/2 * M^2)
    pt = p * (Tt/T)^(γ/(γ-1))
    Δp = pt - p
    q = 1/2 * ρ * TAS^2

    EAS = TAS * √(ρ / ρ_std)
    CAS = √(2γ/(γ-1) * p_std/ρ_std * ( (1 + q/p_std)^((γ-1)/γ) - 1) )

    AirData(v_ew_n, v_ew_b, v_eOb_b, v_wOb_b,
            T, p, ρ, a, μ, M, Tt, pt, Δp, q, TAS, EAS, CAS)

end


function plots(t, data::AbstractVector{<:AirData}; mode, save_path, kwargs...)

    @unpack v_ew_n, v_ew_b, v_eOb_b, v_wOb_b, a, μ, ρ, TAS, EAS, CAS, M,
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
        plot_title = "Velocity (Airframe/ECEF) [Airframe]",
        ylabel = [L"$v_{eb}^{x_b} \ (m/s)$" L"$v_{eb}^{y_b} \ (m/s)$" L"$v_{eb}^{z_b} \ (m/s)$"],
        th_split = :h,
        kwargs...)

    pd["04_v_wOb_b"] = thplot(t, v_wOb_b;
        plot_title = "Velocity (Airframe / Wind) [Airframe]",
        ylabel = [L"$v_{eb}^{x_b} \ (m/s)$" L"$v_{eb}^{y_b} \ (m/s)$" L"$v_{eb}^{z_b} \ (m/s)$"],
        th_split = :h,
        kwargs...)

    splt_a = thplot(t, a;
        title = "Speed of Sound",
        label = "",
        ylabel = L"$a \ (m/s)$",
        kwargs...)

    splt_ρ = thplot(t, ρ;
        title = "Density",
        label = "",
        ylabel = L"$\rho \ (kg/m^3)$",
        kwargs...)

    splt_μ = thplot(t, μ;
        title = "Dynamic Viscosity",
        label = "",
        ylabel = L"$\mu \ (Pa \ s)$",
        kwargs...)

    pd["05_ρ_a"] = plot(splt_ρ, splt_a, splt_μ;
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

    pd["06_T_p"] = plot(splt_T, splt_p;
        plot_title = "Freestream Properties",
        layout = (1,2),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    splt_airspeed = thplot(t, hcat(TAS,EAS,CAS);
        title = "Airspeed",
        label = ["True" "Equivalent" "Calibrated"],
        ylabel = L"$v \ (m/s)$",
        th_split = :none, kwargs...)

    splt_Mach = thplot(t, M;
        title = "Mach",
        label = "",
        ylabel = L"M",
        kwargs...)

    splt_q = thplot(t, q/1000;
        title = "Dynamic Pressure",
        label = "",
        ylabel = L"$q \ (kPa)$",
        kwargs...)

    l3 = @layout [a{0.5w} [b; c{0.5h}]]
    pd["07_airspeed_M_q"] = plot(splt_airspeed, splt_Mach, splt_q;
        layout = l3,
        plot_title = "Freestream Properties",
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    save_plots(pd; save_path)

end

end #module