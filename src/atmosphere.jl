module Atmosphere

using StaticArrays, StructArrays, ComponentArrays

using Flight.Geodesy
using Flight.System

import Flight.System: HybridSystem, get_x0, get_y0, get_u0, get_d0, f_cont!, f_disc!

export SimpleISA, SimpleISASystem, get_ISA_data
export NoWind, get_wind_velocity
export ParametricAtmosphere

export AbstractAtmosphericModel, DummyAtmosphericModel

const R = 287.05287 #gas constant for dry air
const γ = 1.40 #heat capacity ratio for dry air
const T_std = 288.15
const p_std = 101325.0
const g_std = 9.80665

abstract type AbstractISA <: AbstractComponent end

Base.@kwdef struct ISAData
    p::Float64 = p_std
    T::Float64 = T_std
    ρ::Float64 = p_std / (R * T_std)
    a::Float64 = √(γ*R*T_std)
end

const ISA_layers = StructArray(
    β =      SVector{7,Float64}([-6.5e-3, 0, 1e-3, 2.8e-3, 0, -2.8e-3, -2e-3]),
    h_ceil = SVector{7,Float64}([11000, 20000, 32000, 47000, 51000, 71000, 80000]))

@inline ISA_temperature_law(h, T_b, h_b, β) = T_b + β * (h - h_b)

@inline function ISA_pressure_law(h, g0, p_b, T_b, h_b, β)
    if β != 0.0
        p_b * (1 + β/T_b * (h - h_b)) ^ (-g0/(β*R))
    else
       p_b * exp(-g0/(R*T_b) * (h - h_b))
    end
end

function get_SL_values(::T, ::Abstract2DLocation) where {T<:HybridSystem{<:AbstractISA}}
    error("get_SL_values not implemented for $T")
end

@inline function get_ISA_layer_parameters(h::Real)
    for i in 1:length(ISA_layers)
        h < ISA_layers.h_ceil[i] && return ISA_layers[i]
    end
    throw(ArgumentError("Altitude out of bounds"))
end

@inline function get_ISA_data(s::HybridSystem{<:AbstractISA}, p::Abstract3DPosition)

    p_nvg = Geographic{NVector, Geopotential}(p)
    (T_SL, p_SL, g_SL) = get_SL_values(s, p_nvg.loc)
    get_ISA_data(p_nvg.alt, T_SL, p_SL, g_SL)

end

@inline function get_ISA_data(h::AltGeop, T_sl = T_std, p_sl = p_std, g_sl = g_std)

    h_base = 0; T_base = T_sl; p_base = p_sl
    β, h_ceil = get_ISA_layer_parameters(h_base)

    while h > h_ceil
        T_ceil = ISA_temperature_law(h_ceil, T_base, h_base, β)
        p_ceil = ISA_pressure_law(h_ceil, g_sl, p_base, T_base, h_base, β)
        h_base = h_ceil; T_base = T_ceil; p_base = p_ceil
        β, h_ceil = get_ISA_layer_parameters(h_base)
    end
    T = ISA_temperature_law(h, T_base, h_base, β)
    p = ISA_pressure_law(h, g_sl, p_base, T_base, h_base, β)

    return ISAData(p = p, T = T, ρ = p / (R*T), a = √(γ*R*T) )

end

##################### ConstantUniformISA ###########################

#a HybridSystem{<:AbstractISA} may have an output type of its own, as any other
#HybridSystem. however, this output is NOT an ISAData instance. the ISA_data
#returned by a ISA System depends on the location specified in the query.
#instead, the output from an ISA System may hold quantities of interest related
#to its own internal state, if it has one due to it being a dynamic ISA
#implementation.

#AbstractISA subtypes: ConstantUniformISA (SimpleISA), ConstantFieldISA,
#DynamicUniformISA, DynamicFieldISA.

#the constancy of a SimpleISA System is not really such, it simply means that
#the System does not have a state and therefore cannot evolve on its own. but
#its input vector can still be used to change manually the SL conditions during
#simulation.

struct SimpleISA <: AbstractISA end #Constant, Uniform ISA

Base.@kwdef mutable struct USimpleISA
    T_sl::Float64 = T_std
    p_sl::Float64 = p_std
end

get_u0(::SimpleISA) = USimpleISA()
f_cont!(::HybridSystem{<:SimpleISA}, args...) = nothing
f_disc!(::HybridSystem{<:SimpleISA}, args...) = false

function get_SL_values(s::HybridSystem{<:SimpleISA}, ::Abstract2DLocation)
    return (T = s.u.T_sl, p = s.u.p_sl, g = g_std)
    #alternative using actual local SL gravity:
    # return (T = s.u.T_sl, p = s.u.p_sl, g = gravity(Geographic(loc, AltOrth(0.0))))
end



######################## AbstractWind ############################

abstract type AbstractWind <: AbstractComponent end

function get_wind_velocity(::T, ::Abstract2DLocation) where {T<:HybridSystem{<:AbstractWind}}
    error("get_wind_velocity not implemented for $T")
end

struct NoWind <: AbstractWind end
f_cont!(::HybridSystem{<:NoWind}, args...) = nothing
f_disc!(::HybridSystem{<:NoWind}, args...) = false

function get_wind_velocity(::HybridSystem{<:NoWind}, ::Abstract3DPosition)
    zeros(SVector{3,Float64})
end

#################### ParametricAtmosphere ############################

Base.@kwdef struct ParametricAtmosphere{I <: AbstractISA, W <: AbstractWind} <: AbstractComponent
    isa_::I = SimpleISA()
    wind::W = NoWind()
end

Base.@kwdef struct AtmosphericData
    isa_::ISAData = ISAData()
    v_ew_n::SVector{3,Float64} = zeros(3)
end

#this one's convenient for dispatching on plots() methods, but for U and D we
#can use simply NamedTuples
struct ParametricAtmosphereY{I, W}; isa_::I; wind::W; end

get_x0(atm::ParametricAtmosphere) = ComponentVector(isa_ = get_x0(atm.isa_), wind = get_x0(atm.wind))
get_y0(atm::ParametricAtmosphere) = ParametricAtmosphereY(get_y0(atm.isa_), get_y0(atm.wind))
get_u0(atm::ParametricAtmosphere) = (isa_ = get_u0(atm.isa_), wind = get_u0(atm.wind))
get_d0(atm::ParametricAtmosphere) = (isa_ = get_d0(atm.isa_), wind = get_d0(atm.wind))

function HybridSystem(atm::ParametricAtmosphere, ẋ = get_x0(atm), x = get_x0(atm),
                    y = get_y0(atm), u = get_u0(atm), d = get_d0(atm), t = Ref(0.0))

    isa_sys = HybridSystem(atm.isa_, ẋ.isa_, x.isa_, y.isa_, u.isa_, d.isa_, t)
    wind_sys = HybridSystem(atm.wind, ẋ.wind, x.wind, y.wind, u.wind, d.wind, t)
    params = nothing
    subsystems = (isa_ = isa_sys, wind = wind_sys,)
    HybridSystem{map(typeof, (atm, x, y, u, d, params, subsystems))...}(
                                ẋ, x, y, u, d, t, params, subsystems)
end

# function get_atmospheric_data(atm::AbstractAtmosphericModel, p::Abstract3DPosition)

#     # v_ew_n::SVector{3,Float64} #v_ew_n: wind velocity relative to the ECEF, NED axes
# end
#if we dont want to step the AtmosphericSystem in time because we know it is
#Constant, we don't create a World to compose them both. THE IMPORTANT THING IS
#THAT WHAT AIRCRAFT RECEIVES AS ATM MUST BE A SYSTEM, not a component. in case
#that we want to compose the aircraft with the environment in a World model, we
#simply pass it as a subsystem within f_cont!(::World). if the AtmosphericSystem
#is dynamic but still we don't want to compose it with the Aircraft, we can
#always put the atmospheric system in a separate model, step it separately, and
#pass a reference to its underlying AtmosphericSystem to the AircraftSystem that
#will be wrapped in its own Model


abstract type AbstractAtmosphericModel end

struct DummyAtmosphericModel <: AbstractAtmosphericModel end



# #top-down / recursive implementation
# @inline function get_tp(h::Real, T0::Real = T0_std, p0::Real = p0_std, g0::Real = g0_std)

#     h == 0 && return (T0, p0)
#     (h_b, β) = layer_parameters(h)
#     (T_b, p_b) = get_tp(h_b, T0, p0, g0) #get pt at the layer base
#     T = ISA_temperature_law(h, T_b, h_b, β)
#     p = ISA_pressure_law(h, g0, p_b, T_b, h_b, β)
#     return (T, p)
# end







end