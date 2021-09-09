module Atmosphere

using StaticArrays, StructArrays

using Flight.Geodesy
using Flight.System

export ConstantUniformISA, get_ISA_data

export AbstractAtmosphericModel, DummyAtmosphericModel

const R = 287.05287 #gas constant for dry air
const γ = 1.40 #heat capacity ratio for dry air
const T0_std = 288.15
const p0_std = 101325.0
const g0_std = 9.80665

abstract type AbstractISA <: AbstractComponent end

const ISA_layers = StructArray(
    β =      SVector{7,Float64}([-6.5e-3, 0, 1e-3, 2.8e-3, 0, -2.8e-3, -2e-3]),
    h_ceil = SVector{7,Float64}([11000, 20000, 32000, 47000, 51000, 71000, 80000]))

get_T0(::AbstractISA, ::Abstract2DLocation) = error("get_T0 not implemented")
get_p0(::AbstractISA, ::Abstract2DLocation) = error("get_p0 not implemented")
get_g0(::AbstractISA, ::Abstract2DLocation) = error("get_g0 not implemented")

Base.@kwdef struct ISAData
    p::Float64 = p0_std
    T::Float64 = T0_std
    ρ::Float64 = p0_std / (R * T0_std)
    a::Float64 = √(γ*R*T0_std)
end

@inline function get_ISA_layer_parameters(h::Real)
    for i in 1:length(ISA_layers)
        h < ISA_layers.h_ceil[i] && return ISA_layers[i]
    end
    throw(ArgumentError("Altitude out of bounds"))
end

@inline ISA_temperature_law(h, T_b, h_b, β) = T_b + β * (h - h_b)

@inline function ISA_pressure_law(h, g0, p_b, T_b, h_b, β)
    if β != 0.0
        p_b * (1 + β/T_b * (h - h_b)) ^ (-g0/(β*R))
    else
       p_b * exp(-g0/(R*T_b) * (h - h_b))
    end
end

@inline function get_ISA_data(a::AbstractISA, p::Abstract3DPosition)

    p_nvg = Geographic{NVector, Geopotential}(p)

    get_ISA_data(p_nvg.alt,
                 get_T0(a, p_nvg.loc),
                 get_p0(a, p_nvg.loc),
                 get_g0(a, p_nvg.loc))

end

@inline function get_ISA_data(h::AltGeop, T0 = T0_std, p0 = p0_std, g0 = g0_std)

    h_base = 0; T_base = T0; p_base = p0
    β, h_ceil = get_ISA_layer_parameters(h_base)

    while h > h_ceil
        T_ceil = ISA_temperature_law(h_ceil, T_base, h_base, β)
        p_ceil = ISA_pressure_law(h_ceil, g0, p_base, T_base, h_base, β)
        h_base = h_ceil; T_base = T_ceil; p_base = p_ceil
        β, h_ceil = get_ISA_layer_parameters(h_base)
    end
    T = ISA_temperature_law(h, T_base, h_base, β)
    p = ISA_pressure_law(h, g0, p_base, T_base, h_base, β)

    return ISAData(p = p, T = T, ρ = p / (R*T), a = √(γ*R*T) )

end

Base.@kwdef struct ConstantUniformISA <: AbstractISA
    T0::Float64 = T0_std
    p0::Float64 = p0_std
end

get_T0(s::ConstantUniformISA, ::Abstract2DLocation) = s.T0
get_p0(s::ConstantUniformISA, ::Abstract2DLocation) = s.p0
get_g0(::ConstantUniformISA, ::Abstract2DLocation) = g0_std
#alternatively, use local MSL gravity
# get_p0(s::ConstantUniformISA, loc::Abstract2DLocation) = gravity(Geographic(loc,
# AltOrth(0.0)))

function HybridSystem(::ConstantUniformISA)
end

abstract type AbstractWind <: AbstractComponent end

struct ParametricAtmosphere{ISA <: AbstractISA, Wind <: AbstractWind} <: AbstractComponent
    _isa::AbstractISA
    _wind::AbstractWind
end
#NOW... we define a HybridSystem for ConstantUniformISA, which will be trivial then we
#could define ConstantUniformISA, ConstantFieldISA, DynamicUniformISA,
#DynamicFieldISA.

#then we define a ConstantUniformWind, ConstantFieldWind, DynamicUniformWind,
#DynamicFieldWind, also AbstractComponents

#an AtmosphericModel will be parameterized in these two, it doesn't need to be
#Abstract. its constructor will take its type parameters, construct the ISA and
#wind subsystems, and assign them as such. as parameters, the fields atm.isa and
#atm.wind

#WE NEED TO IMPLEMENT get_ISA_data, get_T0, etc for a System, instead of a
#component. internally, the System holds the time, so no need to pass it as an
#argument, only the location

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

function get_atmospheric_data(atm::AbstractAtmosphericModel, p::Abstract3DPosition)

end

# #top-down / recursive implementation
# @inline function layer_parameters(h::Real)
#     if h <= 11000
#         return (h_b = 0.0, β = -6.5e-3)
#     elseif h <= 20000
#         return (h_b = 11000.0, β = 0.0)
#     elseif h <= 32000
#         return (h_b = 20000.0, β = 1e-3)
#     elseif h <= 47000
#         return (h_b = 32000.0, β = 2.8e-3)
#     end
#     throw(ArgumentError("Altitude out of range"))
# end

# @inline function get_tp(h::Real, T0::Real = T0_std, p0::Real = p0_std, g0::Real = g0_std)

#     h == 0 && return (T0, p0)
#     (h_b, β) = layer_parameters(h)
#     (T_b, p_b) = get_tp(h_b, T0, p0, g0) #get pt at the layer base
#     T = ISA_temperature_law(h, T_b, h_b, β)
#     p = ISA_pressure_law(h, g0, p_b, T_b, h_b, β)
#     return (T, p)
# end


#if AtmosphericModel is a system, it will have its own output struct. but its
#own output struct is not necessarily the same as what it produces when queried
#at a specific location. for example, a more sophisticated AtmosphericModel may
#have states that capture the evolution in time of a wind velocity vector field.
#however, when queried at a specific time, it must return atmospheric quantities
#at a unique WGS84 location. the bottom line is that this output, which contains
#these atmospheric quantities, is NOT the same as the output of the
#AtmosphericModel as a dynamical system. we can call the first ones
#AtmosphericData, and the second one AtmosphericModelY
struct DummyAtmosphericModel <: AbstractAtmosphericModel end


    # v_ew_n::SVector{3,Float64} #v_ew_n: wind velocity relative to the ECEF, NED axes





end