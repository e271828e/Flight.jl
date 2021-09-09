module Atmosphere

using StaticArrays

using Flight.Geodesy
using Flight.System

export AbstractAtmosphericModel, DummyAtmosphericModel
export get_ISA_data

const R = 287.05287 #gas constant for dry air
const γ = 1.40 #heat capacity ratio for dry air
const T0_std = 288.15
const p0_std = 101325.0
const g0_std = 9.80665

abstract type AbstractISA <: AbstractComponent end

Base.@kwdef struct ISAData
    p::Float64 = p0_std
    T::Float64 = T0_std
    ρ::Float64 = p0_std / (R * T0_std)
    a::Float64 = √(γ*R*T0_std)
end

Base.@kwdef struct StaticISA <: AbstractISA
    T0::Float64 = T0_std
    p0::Float64 = p0_std
end

@inline temperature_law(h, T_b, h_b, β) = T_b + β * (h - h_b)

@inline function pressure_law(h, g0, p_b, T_b, h_b, β)
    if β != 0.0
        p_b * (1 + β/T_b * (h - h_b)) ^ (-g0/(β*R))
    else
       p_b * exp(-g0/(R*T_b) * (h - h_b))
    end
end

#bottom up / loop implementation
@inline function ISA_layer_parameters(h::Real)
    if h < 11000
        return (β = -6.5e-3, h_ceil = 11000.0)
    elseif h < 20000
        return (β = 0.0, h_ceil = 20000.0)
    elseif h < 32000
        return (β = 1e-3, h_ceil = 32000.0)
    elseif h < 47000
        return (β = 2.8e-3, h_ceil = 47000.0)
    elseif h < 51000
        return (β = 0.0, h_ceil = 51000.0)
    elseif h < 71000
        return (β = -2.8e-3, h_ceil = 71000.0)
    elseif h < 80000
        return (β = -2e-3, h_ceil = 80000.0)
    else
        throw(ArgumentError("Altitude out of range"))
    end
end

@inline function get_ISA_data(a::AbstractISA, p::Abstract3DPosition)

    p_nvo = Geographic{NVector, Orthometric}(p)
    n_e = p_nvo.loc
    h_orth = p_nvo.alt
    h_geop = AltGeop(h_orth)

    #need MSL gravity at the requested 2D location
    p_MSL = Geographic(n_e, AltOrth(0.0))
    g0 = gravity(p_MSL)
    #si procede, puedo sustituir esta g0 dependiente de la latitud por la g_std, que
    #es lo que hace el estandar ISA




end

@inline function get_ISA_data(h::Real, T0::Real = T0_std, p0::Real = p0_std, g0::Real = g0_std)

    h_base = 0; T_base = T0; p_base = p0
    β, h_ceil = ISA_layer_parameters(h_base)

    while h > h_ceil
        T_ceil = temperature_law(h_ceil, T_base, h_base, β)
        p_ceil = pressure_law(h_ceil, g0, p_base, T_base, h_base, β)
        h_base = h_ceil; T_base = T_ceil; p_base = p_ceil
        β, h_ceil = ISA_layer_parameters(h_base)
    end

    return ISAData(T = temperature_law(h, T_base, h_base, β),
                   p = pressure_law(h, g0, p_base, T_base, h_base, β),
                   ρ = p / (R*T),
                   a = √(γ*R*T) )
end

abstract type AbstractAtmosphericModel end

function get_atmospheric_data(atm::AbstractAtmosphericModel, p::Abstract3DPosition)

    #the geometric altitude in the context of the ISA model is orthometric
    h_orth = Altitude{Orthometric}(p)

    h_geop = get_h_geop(Ob)


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
#     T = temperature_law(h, T_b, h_b, β)
#     p = pressure_law(h, g0, p_b, T_b, h_b, β)
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