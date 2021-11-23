module Atmosphere

using StaticArrays, StructArrays, ComponentArrays

using Flight.Geodesy

import Flight.Modeling: SystemDescriptor, System, init_x0, init_y0, init_u0, init_d0, f_cont!, f_disc!

export TunableISA
export TunableWind
export AtmosphereDescriptor, AtmosphericData, AtmosphericSystem

const R = 287.05287 #gas constant for dry air
const γ = 1.40 #heat capacity ratio for dry air
const βs = 1.458e-6 #Sutherland's empirical constant for dynamic viscosity
const S = 110.4 #Sutherland's empirical constant for dynamic viscosity

const T_std = 288.15
const p_std = 101325.0
const ρ_std = p_std / (R * T_std)
const g_std = 9.80665

const ISA_layers = StructArray(
    β =      SVector{7,Float64}([-6.5e-3, 0, 1e-3, 2.8e-3, 0, -2.8e-3, -2e-3]),
    h_ceil = SVector{7,Float64}([11000, 20000, 32000, 47000, 51000, 71000, 84852]))

@inline density(p,T) = p/(R*T)
@inline speed_of_sound(T) = √(γ*R*T)
@inline dynamic_viscosity(T) = (βs * T^1.5) / (T + S)

@inline ISA_temperature_law(h::Real, T_b, h_b, β)::Float64 = T_b + β * (h - h_b)

@inline function ISA_pressure_law(h::Real, g0, p_b, T_b, h_b, β)::Float64
    if β != 0.0
        p_b * (1 + β/T_b * (h - h_b)) ^ (-g0/(β*R))
    else
        p_b * exp(-g0/(R*T_b) * (h - h_b))
    end
end

Base.@kwdef struct SLConditions
    p::Float64 = p_std
    T::Float64 = T_std
    g::Float64 = g_std
end

struct ISAData
    p::Float64
    T::Float64
    ρ::Float64
    a::Float64
    μ::Float64
end

ISAData() = ISAData(AltGeop(0))

@inline function ISAData(h_geo::AltGeop, sl::SLConditions = SLConditions())

    h = Float64(h_geo)
    h_base = 0; T_base = sl.T; p_base = sl.p; g_base = sl.g

    for i in 1:length(ISA_layers)
        β, h_ceil = ISA_layers[i]
        if h < h_ceil
            T = ISA_temperature_law(h, T_base, h_base, β)
            p = ISA_pressure_law(h, g_base, p_base, T_base, h_base, β)
            return ISAData(p, T, density(p, T), speed_of_sound(T), dynamic_viscosity(T) )
        end
        T_ceil = ISA_temperature_law(h_ceil, T_base, h_base, β)
        p_ceil = ISA_pressure_law(h_ceil, g_base, p_base, T_base, h_base, β)
        h_base = h_ceil; T_base = T_ceil; p_base = p_ceil
    end

    throw(ArgumentError("Altitude out of bounds"))

end

##################### AbstractISA ###########################

#a System{<:AbstractISA} may have an output type of its own, as any other
#System. however, this output generally will not be an ISAData instance. the
#ISAData returned by a ISA System depends on the location specified in the
#query. instead, the output from an ISA System may hold quantities of interest
#related to its own internal state, if it has one due to it being a dynamic ISA
#implementation.

#AbstractISA subtypes: ConstantUniformISA (TunableISA), ConstantFieldISA,
#DynamicUniformISA, DynamicFieldISA.

abstract type AbstractISA <: SystemDescriptor end

function SLConditions(::T, ::Abstract2DLocation) where {T<:System{<:AbstractISA}}
    error("SLConditions constructor not implemented for $T")
end

@inline function ISAData(sys::System{<:AbstractISA}, p::Geographic)

    h_geop = Altitude{Geopotential}(p.alt, p.l2d)
    sl = SLConditions(sys, p.l2d)
    ISAData(h_geop, sl)

end

##################### TunableISA ###########################

#a TunableISA System does not have a state and therefore cannot evolve on its
#own. but its input vector can still be used to manually tune the SL conditions
#during simulation.

struct TunableISA <: AbstractISA end #Constant, Uniform ISA

Base.@kwdef mutable struct UTunableISA #only allocates upon System instantiation
    T_sl::Float64 = T_std
    p_sl::Float64 = p_std
end

init_u0(::TunableISA) = UTunableISA()
f_cont!(::System{<:TunableISA}, args...) = nothing
f_disc!(::System{<:TunableISA}, args...) = false

function SLConditions(s::System{<:TunableISA}, ::Abstract2DLocation)
    SLConditions(T = s.u.T_sl, p = s.u.p_sl, g = g_std)
    #alternative using actual local SL gravity:
    # return (T = s.u.T_sl, p = s.u.p_sl, g = gravity(Geographic(l2d, AltOrth(0.0))))
end



######################## AbstractWind ############################

abstract type AbstractWind <: SystemDescriptor end

Base.@kwdef struct WindData
    v_ew_n::SVector{3,Float64} = zeros(SVector{3})
end

function WindData(::T, ::Abstract3DLocation) where {T<:System{<:AbstractWind}}
    error("WindData constructor not implemented for $T")
end


######################## TunableWind ############################

struct TunableWind <: AbstractWind end

Base.@kwdef mutable struct USimpleWind
    v_ew_n::MVector{3,Float64} = zeros(MVector{3}) #MVector allows changing single components
end

init_u0(::TunableWind) = USimpleWind()
f_cont!(::System{<:TunableWind}, args...) = nothing
f_disc!(::System{<:TunableWind}, args...) = false

function WindData(wind::System{<:TunableWind}, ::Abstract3DLocation)
    wind.u.v_ew_n |> SVector{3,Float64} |> WindData
end

#################### AtmosphereDescriptor ############################

Base.@kwdef struct AtmosphereDescriptor{I <: AbstractISA, W <: AbstractWind} <: SystemDescriptor
    isa_::I = TunableISA()
    wind::W = TunableWind()
end

#this one's convenient for dispatching on plots() methods, but for U and D we
#can use simply NamedTuples
struct AtmosphereCmpY{I, W}; isa_::I; wind::W; end

init_x0(atm::AtmosphereDescriptor) = ComponentVector(isa_ = init_x0(atm.isa_), wind = init_x0(atm.wind))
init_y0(atm::AtmosphereDescriptor) = AtmosphereCmpY(init_y0(atm.isa_), init_y0(atm.wind))
init_u0(atm::AtmosphereDescriptor) = (isa_ = init_u0(atm.isa_), wind = init_u0(atm.wind))
init_d0(atm::AtmosphereDescriptor) = (isa_ = init_d0(atm.isa_), wind = init_d0(atm.wind))

function System(atm::AtmosphereDescriptor, ẋ = init_x0(atm), x = init_x0(atm),
                    y = init_y0(atm), u = init_u0(atm), d = init_d0(atm), t = Ref(0.0))

    isa_sys = System(atm.isa_, ẋ.isa_, x.isa_, y.isa_, u.isa_, d.isa_, t)
    wind_sys = System(atm.wind, ẋ.wind, x.wind, y.wind, u.wind, d.wind, t)
    params = nothing
    subsystems = (isa_ = isa_sys, wind = wind_sys,)
    System{map(typeof, (atm, x, y, u, d, params, subsystems))...}(
                                ẋ, x, y, u, d, t, params, subsystems)
end

const AtmosphericSystem = System{<:AtmosphereDescriptor}

Base.@kwdef struct AtmosphericData
    isa_::ISAData = ISAData()
    wind::WindData = WindData()
end

function AtmosphericData(a::AtmosphericSystem, pos::Geographic)
    AtmosphericData(
        ISAData(a.subsystems.isa_, pos),
        WindData(a.subsystems.wind, pos))
end

function f_cont!(sys::AtmosphericSystem)
    f_cont!(sys.isa_)
    f_cont!(sys.wind)
end

#we can create a World <: SystemDescriptor to compose Aircraft and
#AtmosphericSystem, and step them together in a single Model. if the
#AtmosphericSystem is dynamic but still we don't want to compose it with the
#Aircraft, we can always put it in a separate model, step it separately, and
#simply pass a reference to its underlying AtmosphericSystem to the
#AircraftSystem, which will be wrapped in its own Model


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