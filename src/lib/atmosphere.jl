module Atmosphere

using StaticArrays, StructArrays, ComponentArrays, LinearAlgebra, UnPack

using Flight.FlightCore

using ..Attitude
using ..Geodesy
using ..Kinematics

export p_std, T_std, g_std, ρ_std, ISA_layers
export ISAData, AtmosphericData, AirflowData

export SeaLevelStandard, TunableSeaLevel
export NoWind, TunableWind
export AbstractAtmosphere, SimpleAtmosphere

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
################################## ISA #########################################

@kwdef struct ISAData
    T::Float64 = T_std
    p::Float64 = p_std
end

########################## Sea Level Conditions ################################

#model the behavior of SL conditions as a function of time and 2D location
abstract type AbstractSeaLevelConditions <: SystemDefinition end

function ISAData(sys::System{<:AbstractSeaLevelConditions}, loc::Abstract2DLocation)
    MethodError(ISAData, (sys, loc,)) |> throw
end

############################## Sea Level Standard ##############################

struct SeaLevelStandard <: AbstractSeaLevelConditions end
@no_dynamics SeaLevelStandard

ISAData(::System{<:SeaLevelStandard}, ::Abstract2DLocation) = ISAData(T_std, p_std)


############################## Tunable Sea Level ################################

#simple model providing direct control over SL conditions
struct TunableSeaLevel <: AbstractSeaLevelConditions end
@no_dynamics TunableSeaLevel

const T_sl_min = T_std - 50.0
const T_sl_max = T_std + 50.0
const p_sl_min = p_std - 10000.0
const p_sl_max = p_std + 10000.0

@kwdef mutable struct TunableSeaLevelU
    T::Ranged{Float64, T_sl_min, T_sl_max} = T_std
    p::Ranged{Float64, p_sl_min, p_sl_max} = p_std
end

Systems.U(::TunableSeaLevel) = TunableSeaLevelU()

function ISAData(sys::System{<:TunableSeaLevel}, ::Abstract2DLocation)
    ISAData(sys.u.T, sys.u.p)
end

function GUI.draw!(sys::System{<:TunableSeaLevel},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Sea Level Conditions")

    u = sys.u
    CImGui.PushItemWidth(-50)
        # u.T = GUI.safe_slider("T (K)", u.T, "%.3f"; show_label = true)
        u.T = GUI.safe_slider("T (K)", Ref(Float64(u.T)), "%.3f"; show_label = true)[]
        u.p = GUI.safe_slider("p (Pa)", u.p, "%.3f"; show_label = true)
    CImGui.PopItemWidth()
end

############################## ISA Computation #################################

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

#compute T and p at a given geopotential altitude, using ISA_temperature_law and
#ISA_pressure_law to propagate the given sea level conditions upwards through
#the successive ISA_layers up to the requested altitude
function ISAData(h_geo::HGeop, sl::ISAData = ISAData())

    h = Float64(h_geo)
    h_base = 0; T_base = sl.T; p_base = sl.p; g0 = g_std #g0 = sl.g

    for i in eachindex(ISA_layers)
        β, h_ceil = ISA_layers[i]
        if h < h_ceil
            T = ISA_temperature_law(h, T_base, h_base, β)
            p = ISA_pressure_law(h, g0, p_base, T_base, h_base, β)
            return ISAData(; T, p)
        end
        T_ceil = ISA_temperature_law(h_ceil, T_base, h_base, β)
        p_ceil = ISA_pressure_law(h_ceil, g0, p_base, T_base, h_base, β)
        h_base = h_ceil; T_base = T_ceil; p_base = p_ceil
    end

    throw(ArgumentError("Altitude out of bounds"))

end

ISAData(h::Real, args...) = ISAData(HGeop(h), args...)
ISAData(h::Geodesy.AbstractAltitudeDatum, args...) = ISAData(HGeop(h), args...)


################################################################################
#################################### Wind ######################################

abstract type AbstractWind <: SystemDefinition end

function get_wind(sys::System{<:AbstractWind},
                  loc::Abstract3DLocation)::SVector{3,Float64}
    MethodError(get_wind, (sys, loc)) |> throw
end

############################## No Wind Model ###################################

struct NoWind <: AbstractWind end
@no_dynamics NoWind

get_wind(::System{NoWind}, ::Abstract3DLocation) = SVector{3,Float64}(0.0, 0.0, 0.0)

############################## Tunable Wind ####################################

struct TunableWind <: AbstractWind end
@no_dynamics TunableWind

get_wind(sys::System{TunableWind}, ::Abstract3DLocation) = SVector{3}(sys.u)

Systems.U(::TunableWind) = ComponentVector(N= 0.0, E = 0.0, D = 0.0)

function GUI.draw!(sys::System{<:TunableWind},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Wind")

    u = sys.u

    CImGui.PushItemWidth(-80)
        u.N = GUI.safe_slider("North (m/s)", u.N, -30, 30, "%.3f"; show_label = true)
        u.E = GUI.safe_slider("East (m/s)", u.E, -30, 30, "%.3f"; show_label = true)
        u.D = GUI.safe_slider("Down (m/s)", u.D, -30, 30, "%.3f"; show_label = true)
    CImGui.PopItemWidth()

    GUI.draw(sys, label)

end

################################################################################
############################## AtmosphericData #################################

@kwdef struct AtmosphericData #local atmospheric data
    T::Float64 = T_std
    p::Float64 = p_std
    ρ::Float64 = density(p_std, T_std)
    a::Float64 = speed_of_sound(T_std)
    μ::Float64 = dynamic_viscosity(T_std)
    v::SVector{3,Float64} = @SVector[0.0, 0.0, 0.0] #local wind velocity, NED axes
end

################################################################################
############################## AirflowData #####################################

struct AirflowData
    v_ew_n::SVector{3,Float64} #wind velocity, NED axes
    v_ew_b::SVector{3,Float64} #wind velocity, vehicle axes
    v_wb_b::SVector{3,Float64} #vehicle aerodynamic velocity, vehicle axes
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

TAS2EAS(TAS::Real; ρ::Real) = TAS * √(ρ / ρ_std)
EAS2TAS(TAS::Real; ρ::Real) = TAS * √(ρ_std / ρ)

function AirflowData(atm::AtmosphericData = AtmosphericData(),
                    kin::KinData = KinData())

    @unpack v_eb_b, q_nb = kin
    @unpack T, p, ρ, a, μ, v  = atm

    v_ew_n = v #ECEF-relative wind velocity, NED axes
    v_ew_b = q_nb'(v_ew_n)
    v_wb_b = v_eb_b - v_ew_b

    TAS = norm(v_wb_b)
    M = TAS / a
    Tt = T * (1 + (γ - 1)/2 * M^2)
    pt = p * (Tt/T)^(γ/(γ-1))
    Δp = pt - p #impact pressure
    q = 1/2 * ρ * TAS^2 #true dynamic pressure

    EAS = TAS2EAS(TAS; ρ)
    CAS = √(2γ/(γ-1) * p_std/ρ_std * ( (1 + Δp/p_std)^((γ-1)/γ) - 1) )

    AirflowData(v_ew_n, v_ew_b, v_wb_b, T, p, ρ, a, μ, M, Tt, pt, Δp, q, TAS, EAS, CAS)

end


################################################################################
############################ AbstractAtmosphere ################################

#can be static or dynamic, uniform or non-uniform
abstract type AbstractAtmosphere <: SystemDefinition end

function AtmosphericData(sys::System{<:AbstractAtmosphere}, loc::Abstract3DLocation)
    MethodError(AtmosphericData, (sys, loc,)) |> throw
end

################################################################################
############################## SimpleAtmosphere ################################

#Combines an ISA hydrostatic model with an arbitrary wind model

@kwdef struct SimpleAtmosphere{S <: AbstractSeaLevelConditions, W <: AbstractWind} <: AbstractAtmosphere
    sl::S = TunableSeaLevel()
    wind::W = TunableWind()
end

@ss_ode SimpleAtmosphere
@ss_disc SimpleAtmosphere
@no_step SimpleAtmosphere

function AtmosphericData(sys::System{<:SimpleAtmosphere},
                        loc::Abstract3DLocation)

    @unpack T, p = ISAData(HGeop(loc), ISAData(sys.sl, NVector(loc)))
    ρ = density(p, T)
    a = speed_of_sound(T)
    μ = dynamic_viscosity(T)
    v = get_wind(sys.wind, loc)
    AtmosphericData(; T, p, ρ, a, μ, v)
end

function AirflowData(sys::System{<:SimpleAtmosphere}, kin_data::KinData)
    @unpack n_e, h_o = kin_data
    AirflowData(AtmosphericData(sys, Geographic(n_e, h_o)), kin_data)
end

function GUI.draw!(sys::System{<:SimpleAtmosphere},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Simple Atmosphere")

    CImGui.Begin(label, p_open)
        if CImGui.CollapsingHeader("Sea-Level Conditions")
            GUI.draw!(sys.sl)
        end
        if CImGui.CollapsingHeader("Wind")
            GUI.draw!(sys.wind)
        end
    CImGui.End()

end


################################################################################
############################## Aerodynamics ####################################

#Aerodynamics frame (a): Reference frame, rigidly attached to the vehicle,
#wherein aerodynamic forces and moments are computed. It may or may not be the
#same as the vehicle frame b. The frame transform t_ba from b to a is constant
#and known a priori.

#Aerodynamic velocity vector (v_wa): ECEF-relative velocity of the aerodynamics
#reference frame v_ea minus the (local) ECEF-relative wind velocity vector v_ew

#Stability frame (s): Reference frame obtained from rotating the aerodynamics
#frame around its y-axis an angle -α such that the aerodynamic velocity
#vector lies in the resulting plane Os-xs-ys

#Wind frame (w): Reference frame obtained from rotating the stability frame
#around its z-axis an angle β such that the resulting xw axis is aligned with
#the aerodynamic velocity vector

const TAS_min_αβ = 0.1 #minimum TAS for valid α, β computation

#compute aerodynamic velocity vector from TAS and airflow angles
@inline function get_velocity_vector(TAS::Real, α::Real, β::Real)
    cos_β = cos(β)
    return TAS * SVector(cos(α) * cos_β, sin(β), sin(α) * cos_β)
end

#compute airflow angles at frame a from the a-frame aerodynamic velocity vector
@inline function get_airflow_angles(v_wa_a::AbstractVector{<:Real})::Tuple{Float64, Float64}
    if norm(v_wa_a) < TAS_min_αβ
        return (0.0, 0.0)
    else
        α = atan(v_wa_a[3], v_wa_a[1])
        β = atan(v_wa_a[2], √(v_wa_a[1]^2 + v_wa_a[3]^2))
        return (α, β)
    end
end

#compute the rotation from a-frame to w-frame given the a-frame components of
#the aerodynamic velocity vector
@inline function get_wind_axes(v_wa_a::AbstractVector{<:Real})
    α, β = get_airflow_angles(v_wa_a)
    get_wind_axes(α, β)
end

#compute the rotation from a-frame to w-frame from the a-frame airflow angles
@inline function get_wind_axes(α::Real, β::Real)
    q_aw = Ry(-α) ∘ Rz(β)
    return q_aw
end

#compute the rotation from a-frame to s-frame from the a-frame airflow angles
@inline function get_stability_axes(α::Real)
    q_as = Ry(-α)
    return q_as
end



################################## Plotting ####################################

function Plotting.make_plots(ts::TimeSeries{<:AirflowData}; kwargs...)

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

    pd[:v_wb_b] = plot(ts.v_wb_b;
        plot_title = "Velocity (Vehicle / Wind) [Vehicle Axes]",
        ylabel = [L"$v_{eb}^{x_b} \ (m/s)$" L"$v_{eb}^{y_b} \ (m/s)$" L"$v_{eb}^{z_b} \ (m/s)$"],
        ts_split = :h,
        kwargs...)

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

function GUI.draw(air::AirflowData, p_open::Ref{Bool} = Ref(true),
                label::String = "Airflow Data")

    @unpack v_ew_n, v_ew_b, v_wb_b, T, p, ρ, a, μ, M, Tt, pt, Δp, q, TAS, EAS, CAS = air

    CImGui.Begin(label, p_open)

    GUI.draw(v_ew_n, "Velocity (Wind/ECEF) [NED]", "m/s")
    GUI.draw(v_ew_b, "Velocity (Wind/ECEF) [Body]", "m/s")
    GUI.draw(v_wb_b, "Velocity (Body/Wind) [Body]", "m/s")

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