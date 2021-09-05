module Atmosphere

using StaticArrays

using Flight.Geodesy

export AbstractAtmosphericModel, DummyAtmosphericModel

const R = 287.053 #gas constant for dry air

abstract type AbstractAtmosphericModel end

#if AtmosphericModel is a system, it will have its own output struct. but its
#own output struct is not necessarily the same as what it produces when queried
#at a specific location. for example, a more sophisticated AtmosphericModel may
#have states that capture the evolution in time of a wind velocity vector field.
#however, when queried at a specific time, it must return atmospheric quantities
#at a unique WGS84 location. the bottom line is that this output, which contains
#these atmospheric quantities, is NOT the same as the output of the
#AtmosphericModel as a dynamical system. we can call the first ones
#AtmosphericData, and the second one AtmosphericModelY

struct AtmosphericData
    p::Float64
    T::Float64
    ρ::Float64
    a::Float64
    v_ew_n::SVector{3,Float64} #v_ew_n: wind velocity relative to the ECEF, NED axes
end

get_atmospheric_data(::AbstractAtmosphericModel) = error("To be extended")


struct DummyAtmosphericModel <: AbstractAtmosphericModel end

#could probably separate the atmospheric model in fluidostatic model, wind
#model and possibly turbulence model. for now,

Base.@kwdef struct StaticISAModel <: AbstractAtmosphericModel
    p₀::Float64 = 101325
    T₀::Float64 = 288.15
end
# ρ₀ = p₀/R*T₀

function get_atmospheric_data(atm::StaticISAModel, Ob::WGS84Pos)

    #the geometric altitude in the context of the ISA model is orthometric: it
    #is measured with respect to MSL, that is, from the geoid's surface, rather
    #than the ellipsoid's. in contrast, the h provided by the Kinematics module
    #is relative to the WGS84 ellipsoid. so, the first step is to convert
    #altitude to orthometric. this method must be provided by the Geodesy module.
    #if properly implemented, it would have to make use of EGM96, EGM08 or EGM20
    #to evaluate the geoid's height N above the ellipsoid at the requested
    #location, then do: h_orth = h - N(ϕ, λ)
    h_orth = get_h_orth(Ob)

    #now, the fluidostatic equations used to compute pressure and density as
    #functions of altitude in the ISA model assume constant gravity. this is
    #only true if we use geopotential altitude. using a reasonable approximation
    #for the variation of gravity potential with altitude, we can compute it as:
    h_geop = get_h_geop(Ob)
    # h_geop = a * h_orth / (a + h_orth)

    #see Allerton

    #and with this we can now compute p, T, ρ for the different ISA layers

    #however...
    #ultimately, using ellipsoid or orthometric height in the ISA model is not
    #very important. if we use ellipsoid height, we could simply interpret p0
    #and T0 as being given at h=0 rather than h_orth=0 (MSL). in the end, these
    #are arbitrary standard values. in real world conditions, p0 and T0 will
    #have different values at different times.

    #also, at 20000 m (the start of the Stratosphere layer), the difference
    #between orthometric and geopotential altitude is 63m.

    #bottom line: for now, we can omit all this, and simply take the WGS84
    #ellipsoid altitude and plug it directly as orthometric-geopotential
    #altitude. improving this is NOT a priority

    #OTOH, computing geopotential altitude from geometric altitude is quite
    #easy, so we could do it. but always documenting its input as orthometric.
    #otherwise, it makes no sense conceptually. the simplification must lie in
    #the approximation h_orth ≈ h

    #NEED TO DOCUMENT ALL THIS

end

end