module Atmosphere

using StaticArrays, StructArrays, ComponentArrays, LinearAlgebra, UnPack

using Flight.FlightCore
using Flight.FlightCore.Utils: Ranged

using ..Attitude
using ..Geodesy
using ..Kinematics

export LocalAtmosphere, AtmData, ISAData, AirData
export p_std, T_std, g_std, ρ_std, ISA_layers
export get_velocity_vector, get_airflow_angles, get_wind_axes, get_stability_axes

### see ISO 2553

const R = 287.05287 #gas constant for dry ISA
const γ = 1.40 #heat capacity ratio for dry ISA
const βs = 1.458e-6 #Sutherland's empirical constant for dynamic viscosity
const S = 110.4 #Sutherland's empirical constant for dynamic viscosity

const T_std = 288.15
const p_std = 101325.0
const ρ_std = p_std / (R * T_std)
const g_std = 9.80665

@inline density(p,T) = p/(R*T)
@inline speed_of_sound(T) = √(γ*R*T)
@inline dynamic_viscosity(T) = (βs * T^1.5) / (T + S)
@inline SI2kts(v::Real) = 1.94384v


################################################################################
############################## LocalAtmosphere #################################

struct LocalAtmosphere <: SystemDefinition end

const T_sl_min = T_std - 50.0
const T_sl_max = T_std + 50.0
const p_sl_min = p_std - 10000.0
const p_sl_max = p_std + 10000.0

@kwdef mutable struct LocalAtmosphereU
    T_sl::Ranged{Float64, T_sl_min, T_sl_max} = T_std
    p_sl::Ranged{Float64, p_sl_min, p_sl_max} = p_std
    v_ew_n::MVector{3,Float64} = zeros(MVector{3})
end

@kwdef struct AtmData
    T_sl::Float64 = T_std
    p_sl::Float64 = p_std
    v_ew_n::SVector{3,Float64} = zeros(SVector{3})
end

Systems.U(::LocalAtmosphere) = LocalAtmosphereU()

function AtmData(sys::System{LocalAtmosphere})
    @unpack T_sl, p_sl, v_ew_n = sys.u
    AtmData(; T_sl, p_sl, v_ew_n = SVector{3,Float64}(v_ew_n))
end

function GUI.draw!(sys::System{<:LocalAtmosphere},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Local Atmosphere")

    u = sys.u

    CImGui.Begin(label, p_open)

    CImGui.PushItemWidth(-250)
    u.T_sl = GUI.safe_slider("Sea Level Temperature (K)", u.T_sl, "%.3f", true)
    u.p_sl = GUI.safe_slider("Sea Level Pressure (Pa)", u.p_sl, "%.3f", true)
    u.v_ew_n[1] = GUI.safe_slider("North Wind (m/s)", u.v_ew_n[1], -30, 30, "%.3f", true)
    u.v_ew_n[2] = GUI.safe_slider("East Wind (m/s)", u.v_ew_n[2], -30, 30, "%.3f", true)
    u.v_ew_n[3] = GUI.safe_slider("Down Wind (m/s)", u.v_ew_n[3], -30, 30, "%.3f", true)
    CImGui.PopItemWidth()

    CImGui.End()

    GUI.draw(sys, label)

end

################################################################################
############################### ISAData ########################################

struct ISAData
    p::Float64
    T::Float64
    ρ::Float64
    a::Float64
    μ::Float64
end

const ISA_layers = StructArray(
    β =      SVector{7,Float64}([-6.5e-3, 0, 1e-3, 2.8e-3, 0, -2.8e-3, -2e-3]),
    h_ceil = SVector{7,Float64}([11000, 20000, 32000, 47000, 51000, 71000, 84852]))

ISA_temperature_law(h::Real, T_b, h_b, β)::Float64 = T_b + β * (h - h_b)

function ISA_pressure_law(h::Real, g0, p_b, T_b, h_b, β)::Float64
    if β != 0.0
        p_b * (1 + β/T_b * (h - h_b)) ^ (-g0/(β*R))
    else
        p_b * exp(-g0/(R*T_b) * (h - h_b))
    end
end

#compute ISAData at a given geopotential altitude, using ISA_temperature_law and
#ISA_pressure_law to propagate the given sea level conditions upwards through
#the successive ISA_layers up to the requested altitude
function ISAData(h_geo::HGeop; T_sl::Real = T_std, p_sl::Real = p_std)

    h = Float64(h_geo)
    h_base = 0; T_base = T_sl; p_base = p_sl; g0 = g_std #g0 = sl.g

    for i in eachindex(ISA_layers)
        β, h_ceil = ISA_layers[i]
        if h < h_ceil
            T = ISA_temperature_law(h, T_base, h_base, β)
            p = ISA_pressure_law(h, g0, p_base, T_base, h_base, β)
            return ISAData(p, T, density(p, T), speed_of_sound(T), dynamic_viscosity(T) )
        end
        T_ceil = ISA_temperature_law(h_ceil, T_base, h_base, β)
        p_ceil = ISA_pressure_law(h_ceil, g0, p_base, T_base, h_base, β)
        h_base = h_ceil; T_base = T_ceil; p_base = p_ceil
    end

    throw(ArgumentError("Altitude out of bounds"))

end

ISAData(; kwargs...) = ISAData(HGeop(0); kwargs...)
ISAData(h::Real; kwargs...) = ISAData(HGeop(h); kwargs...)
ISAData(h_orth::HOrth; kwargs...) = ISAData(HGeop(h_orth); kwargs...)
ISAData(loc::Abstract3DLocation; kwargs...) = ISAData(HGeop(loc); kwargs...)

function ISAData(loc::Abstract3DLocation, atm::AtmData)
    ISAData(loc; T_sl = atm.T_sl, p_sl = atm.p_sl)
end

# #top-down / recursive implementation
# @inline function get_tp(h::Real, T0::Real = T0_std, p0::Real = p0_std, g0::Real = g0_std)

#     h == 0 && return (T0, p0)
#     (h_b, β) = layer_parameters(h)
#     (T_b, p_b) = get_tp(h_b, T0, p0, g0) #get pt at the layer base
#     T = ISA_temperature_law(h, T_b, h_b, β)
#     p = ISA_pressure_law(h, g0, p_b, T_b, h_b, β)
#     return (T, p)
# end


################################################################################
############################### AirData ########################################

const TAS_min_αβ = 0.1 #minimum TAS for valid α, β computation

struct AirData
    v_ew_n::SVector{3,Float64} #wind velocity, NED axes
    v_ew_b::SVector{3,Float64} #wind velocity, vehicle axes
    v_eOb_b::SVector{3,Float64} #vehicle velocity vector, vehicle axes
    v_wOb_b::SVector{3,Float64} #vehicle aerodynamic velocity, vehicle axes
    α_b::Float64 #vehicle frame AoA
    β_b::Float64 #vehicle frame AoS
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

AirData() = AirData(KinData(), AtmData())

TAS2EAS(TAS::Real; ρ::Real) = TAS * √(ρ / ρ_std)
EAS2TAS(TAS::Real; ρ::Real) = TAS * √(ρ_std / ρ)

function AirData(kin::KinData, atm::AtmData)

    @unpack h_o, v_eOb_b, q_nb = kin
    @unpack T_sl, p_sl, v_ew_n = atm

    v_ew_b = q_nb'(v_ew_n)
    v_wOb_b = v_eOb_b - v_ew_b
    α_b, β_b = get_airflow_angles(v_wOb_b)

    @unpack T, p, ρ, a, μ = ISAData(h_o; T_sl, p_sl)
    TAS = norm(v_wOb_b)
    M = TAS / a
    Tt = T * (1 + (γ - 1)/2 * M^2)
    pt = p * (Tt/T)^(γ/(γ-1))
    Δp = pt - p #impact pressure
    q = 1/2 * ρ * TAS^2 #true dynamic pressure

    EAS = TAS2EAS(TAS; ρ)
    CAS = √(2γ/(γ-1) * p_std/ρ_std * ( (1 + Δp/p_std)^((γ-1)/γ) - 1) )

    AirData(v_ew_n, v_ew_b, v_eOb_b, v_wOb_b, α_b, β_b, T, p, ρ, a, μ, M, Tt, pt, Δp, q, TAS, EAS, CAS)

end

#compute aerodynamic velocity vector from TAS and airflow angles
@inline function get_velocity_vector(TAS::Real, α::Real, β::Real)
    cos_β = cos(β)
    return TAS * SVector(cos(α) * cos_β, sin(β), sin(α) * cos_β)
end

#compute airflow angles at frame c from the c-frame aerodynamic velocity
@inline function get_airflow_angles(v_wOc_c::AbstractVector{<:Real})::Tuple{Float64, Float64}
    if norm(v_wOc_c) < TAS_min_αβ
        return (0.0, 0.0)
    else
        α = atan(v_wOc_c[3], v_wOc_c[1])
        β = atan(v_wOc_c[2], √(v_wOc_c[1]^2 + v_wOc_c[3]^2))
        return (α, β)
    end
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


################################## Plotting ####################################

function Plotting.make_plots(ts::TimeSeries{<:AirData}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    pd[:v_ew_n] = plot(ts.v_ew_n;
        plot_title = "Velocity (Wind / ECEF) [NED Axes]",
        label = ["North" "East" "Down"],
        ylabel = [L"$v_{ew}^{N} \ (m/s)$" L"$v_{ew}^{E} \ (m/s)$" L"$v_{ew}^{D} \ (m/s)$"],
        ts_split = :h,
        kwargs...)

    pd[:v_ew_b] = plot(ts.v_ew_b;
        plot_title = "Velocity (Wind / ECEF) [Vehicle Axes]",
        ylabel = [L"$v_{ew}^{x_b} \ (m/s)$" L"$v_{ew}^{y_b} \ (m/s)$" L"$v_{ew}^{z_b} \ (m/s)$"],
        ts_split = :h,
        kwargs...)

    pd[:v_eOb_b] = plot(ts.v_eOb_b;
        plot_title = "Velocity (Vehicle / ECEF) [Vehicle Axes]",
        ylabel = [L"$v_{eb}^{x_b} \ (m/s)$" L"$v_{eb}^{y_b} \ (m/s)$" L"$v_{eb}^{z_b} \ (m/s)$"],
        ts_split = :h,
        kwargs...)

    pd[:v_wOb_b] = plot(ts.v_wOb_b;
        plot_title = "Velocity (Vehicle / Wind) [Vehicle Axes]",
        ylabel = [L"$v_{eb}^{x_b} \ (m/s)$" L"$v_{eb}^{y_b} \ (m/s)$" L"$v_{eb}^{z_b} \ (m/s)$"],
        ts_split = :h,
        kwargs...)

        subplot_α = plot(TimeSeries(ts._t, rad2deg.(ts.α_b._data));
            title = "Angle of Attack", ylabel = L"$α_b \ (deg)$",
            label = "", kwargs...)

        subplot_β = plot(TimeSeries(ts._t, rad2deg.(ts.β_b._data));
            title = "Angle of Sideslip", ylabel = L"$β_b \ (deg)$",
            label = "", kwargs...)

    pd[:α_β] = plot(subplot_α, subplot_β;
        plot_title = "Airflow Angles [Vehicle Axes]",
        layout = (1,2),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

        subplot_a = plot(ts.a;
            title = "Speed of Sound", ylabel = L"$a \ (m/s)$",
            label = "", kwargs...)

        subplot_ρ = plot(ts.ρ;
            title = "Density", ylabel = L"$\rho \ (kg/m^3)$",
            label = "", kwargs...)

        subplot_μ = plot(ts.μ;
            title = "Dynamic Viscosity", ylabel = L"$\mu \ (Pa \ s)$",
            label = "", kwargs...)

    pd[:ρ_a] = plot(subplot_ρ, subplot_a, subplot_μ;
        plot_title = "Freestream Properties",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs


        subplot_T = plot(
            TimeSeries(ts._t, hcat(ts.T._data, ts.Tt._data)' |> collect);
            title = "Temperature",
            label = ["Static"  "Total"],
            ylabel = L"$T \ (K)$",
            ts_split = :none, kwargs...)

        subplot_p = plot(
            TimeSeries(ts._t, 1e-3*hcat(ts.p._data, ts.pt._data)' |> collect);
            title = "Pressure",
            label = ["Static"  "Total"],
            ylabel = L"$p \ (kPa)$",
            ts_split = :none, kwargs...)

    pd[:T_p] = plot(subplot_T, subplot_p;
        plot_title = "Freestream Properties",
        layout = (1,2),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

        subplot_airspeed = plot(
            TimeSeries(ts._t, hcat(ts.TAS._data, ts.EAS._data, ts.CAS._data)' |> collect);
            title = "Airspeed",
            label = ["True" "Equivalent" "Calibrated"],
            ylabel = L"$v \ (m/s)$",
            ts_split = :none, kwargs...)

        subplot_Mach = plot(ts.M;
            title = "Mach", ylabel = L"M",
            label = "", kwargs...)

        subplot_q = plot(ts._t, ts.q._data/1000;
            title = "Dynamic Pressure", ylabel = L"$q \ (kPa)$",
            label = "", kwargs...)

    l3 = @layout [a{0.5w} [b; c{0.5h}]]

    pd[:airspeed_M_q] = plot(
        subplot_airspeed, subplot_Mach, subplot_q;
        layout = l3,
        plot_title = "Freestream Properties",
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end

################################# GUI ##########################################

function GUI.draw(air::AirData, p_open::Ref{Bool} = Ref(true), label::String = "Air")

    @unpack v_ew_n, v_ew_b, v_wOb_b, α_b, β_b, T, p, ρ, a, μ, M, Tt, pt, Δp, q, TAS, EAS, CAS = air

    CImGui.Begin(label, p_open)

    GUI.draw(v_ew_n, "Velocity (Wind/ECEF) [NED]", "m/s")
    GUI.draw(v_ew_b, "Velocity (Wind/ECEF) [Body]", "m/s")
    GUI.draw(v_wOb_b, "Velocity (Body/Wind) [Body]", "m/s")
    CImGui.Text(@sprintf("AoA (Body/Wind): %.3f deg", rad2deg(α_b)))
    CImGui.Text(@sprintf("AoS (Body/Wind): %.3f deg", rad2deg(β_b)))

    CImGui.Text(@sprintf("Static Temperature: %.3f K", T))
    CImGui.Text(@sprintf("Total Temperature: %.3f K", Tt))
    CImGui.Text(@sprintf("Static Pressure: %.3f Pa", p))
    CImGui.Text(@sprintf("Total Pressure: %.3f Pa", pt))
    CImGui.Text(@sprintf("Impact Pressure: %.3f Pa", Δp))
    CImGui.Text(@sprintf("Dynamic Pressure: %.3f Pa", q))
    CImGui.Text(@sprintf("Density: %.3f kg/m3", ρ))
    CImGui.Text(@sprintf("Speed of Sound: %.3f m/s", a))
    CImGui.Text(@sprintf("Mach: %.3f", M))
    CImGui.Text(@sprintf("CAS: %.3f m/s", CAS))
    CImGui.Text(@sprintf("EAS: %.3f m/s", EAS))
    CImGui.Text(@sprintf("TAS: %.3f m/s", TAS))
    CImGui.Text(@sprintf("CAS: %.3f kts", SI2kts(CAS)))
    CImGui.Text(@sprintf("EAS: %.3f kts", SI2kts(EAS)))
    CImGui.Text(@sprintf("TAS: %.3f kts", SI2kts(TAS)))

    CImGui.End()

end


end