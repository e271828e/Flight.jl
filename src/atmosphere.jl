module Atmosphere

using StaticArrays, StructArrays, ComponentArrays

using Flight.Geodesy
using Flight.System

import Flight.System: HybridSystem, get_x0, get_y0, get_u0, get_d0, f_cont!, f_disc!

export SimpleISA, get_ISA_data
export SimpleWind, get_wind_velocity
export AtmosphereCmp, AtmosphericData, AtmosphericSystem

const R = 287.05287 #gas constant for dry air
const γ = 1.40 #heat capacity ratio for dry air
const T_std = 288.15
const p_std = 101325.0
const ρ_std = p_std / (R * T_std)
const g_std = 9.80665

abstract type AbstractISA <: AbstractComponent end

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

@inline function get_ISA_layer_parameters(h::Real)
    for i in 1:length(ISA_layers)
        h < ISA_layers.h_ceil[i] && return ISA_layers[i]
    end
    throw(ArgumentError("Altitude out of bounds"))
end

Base.@kwdef struct SLConditions
    p::Float64 = p_std
    T::Float64 = T_std
    g::Float64 = g_std
end

Base.@kwdef struct ISAData
    p::Float64 = p_std
    T::Float64 = T_std
    ρ::Float64 = p_std / (R * T_std)
    a::Float64 = √(γ*R*T_std)
end

function SLConditions(::T, ::Abstract2DLocation) where {T<:HybridSystem{<:AbstractISA}}
    error("SLConditions constructor not implemented for $T")
end

@inline function ISAData(h::AltGeop; sl::SLConditions = SLConditions())

    h_base = 0; T_base = sl.T; p_base = sl.p; g_base = sl.g
    β, h_ceil = get_ISA_layer_parameters(h_base)

    while h > h_ceil
        T_ceil = ISA_temperature_law(h_ceil, T_base, h_base, β)
        p_ceil = ISA_pressure_law(h_ceil, g_base, p_base, T_base, h_base, β)
        h_base = h_ceil; T_base = T_ceil; p_base = p_ceil
        β, h_ceil = get_ISA_layer_parameters(h_base)
    end
    T = ISA_temperature_law(h, T_base, h_base, β)
    p = ISA_pressure_law(h, g_base, p_base, T_base, h_base, β)

    return ISAData(p, T, p / (R*T), √(γ*R*T) )

end

@inline function ISAData(sys::HybridSystem{<:AbstractISA}, p::Geographic)

    h_geop = Altitude{Geopotential}(p.alt, p.loc)
    sl = SLConditions(sys, p.loc)
    ISAData(h_geop; sl)

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

function SLConditions(s::HybridSystem{<:SimpleISA}, ::Abstract2DLocation)
    SLConditions(T = s.u.T_sl, p = s.u.p_sl, g = g_std)
    #alternative using actual local SL gravity:
    # return (T = s.u.T_sl, p = s.u.p_sl, g = gravity(Geographic(loc, AltOrth(0.0))))
end



######################## AbstractWind ############################

abstract type AbstractWind <: AbstractComponent end

Base.@kwdef struct WindData
    v_ew_n::SVector{3,Float64} = zeros(SVector{3})
end

function WindData(::T, ::Abstract2DLocation) where {T<:HybridSystem{<:AbstractWind}}
    error("WindData constructor not implemented for $T")
end

struct SimpleWind <: AbstractWind end

Base.@kwdef mutable struct USimpleWind
    v_ew_n::MVector{3,Float64} = zeros(3)
end

get_u0(::SimpleWind) = USimpleWind()
f_cont!(::HybridSystem{<:SimpleWind}, args...) = nothing
f_disc!(::HybridSystem{<:SimpleWind}, args...) = false

function WindData(wind::HybridSystem{<:SimpleWind}, ::Abstract3DPosition)
    wind.u.v_ew_n |> SVector{3,Float64} |> WindData
end

#################### AtmosphereCmp ############################

Base.@kwdef struct AtmosphereCmp{I <: AbstractISA, W <: AbstractWind} <: AbstractComponent
    isa_::I = SimpleISA()
    wind::W = SimpleWind()
end

#this one's convenient for dispatching on plots() methods, but for U and D we
#can use simply NamedTuples
struct AtmosphereCmpY{I, W}; isa_::I; wind::W; end

get_x0(atm::AtmosphereCmp) = ComponentVector(isa_ = get_x0(atm.isa_), wind = get_x0(atm.wind))
get_y0(atm::AtmosphereCmp) = AtmosphereCmpY(get_y0(atm.isa_), get_y0(atm.wind))
get_u0(atm::AtmosphereCmp) = (isa_ = get_u0(atm.isa_), wind = get_u0(atm.wind))
get_d0(atm::AtmosphereCmp) = (isa_ = get_d0(atm.isa_), wind = get_d0(atm.wind))

function HybridSystem(atm::AtmosphereCmp, ẋ = get_x0(atm), x = get_x0(atm),
                    y = get_y0(atm), u = get_u0(atm), d = get_d0(atm), t = Ref(0.0))

    isa_sys = HybridSystem(atm.isa_, ẋ.isa_, x.isa_, y.isa_, u.isa_, d.isa_, t)
    wind_sys = HybridSystem(atm.wind, ẋ.wind, x.wind, y.wind, u.wind, d.wind, t)
    params = nothing
    subsystems = (isa_ = isa_sys, wind = wind_sys,)
    HybridSystem{map(typeof, (atm, x, y, u, d, params, subsystems))...}(
                                ẋ, x, y, u, d, t, params, subsystems)
end

const AtmosphericSystem = HybridSystem{<:AtmosphereCmp}

struct AtmosphericData
    isa_::ISAData
    wind::WindData
end

function AtmosphericData(a::AtmosphericSystem, pos::Geographic)
    AtmosphericData(
        ISAData(a.subsystems.isa_, pos),
        WindData(a.subsystems.wind, pos))
end

#if we dont want to step the AtmosphericSystem in time because it is constant,
#we simply don't create a World to compose Aircraft and AtmosphericSystem. in
#case that we want to compose the aircraft with the environment in a World
#model, we simply pass it as a subsystem within f_cont!(::World). if the
#AtmosphericSystem is dynamic but still we don't want to compose it with the
#Aircraft, we can always put the atmospheric system in a separate model, step it
#separately, and pass a reference to its underlying AtmosphericSystem to the
#AircraftSystem that will be wrapped in its own Model




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