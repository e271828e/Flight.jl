module WGS84

using LinearAlgebra
using StaticArrays: SVector
using Flight.Attitude

export WGS84Pos, NVectorAlt, Cartesian

#WGS84 fundamental constants, SI units
const GM = 3.986005e+14 #Gravitational constant
const a = 6378137.0 #Equatorial radius
const f = 1/298.257223563 #Ellipsoid flattening
const ω_ie = 7.292115e-05 #Earth's angular velocity with respect to the ECI frame

#derived parameters
const b = a * (1 - f) #Polar semi-minor axis
const e² = 2f - f^2 #First eccentricity squared (^ operator calls _power_by_squaring, no need to write f*f)
const e = √e² #First eccentricity
const eʹ² = e² / (1 - e²) #Second eccentricity squared
const eʹ = √eʹ² #Second eccentricity

#convenience parameters
const a² = a^2
const b² = b^2
const q₀ = 0.5( (1 + 3/eʹ²) * atan(eʹ) - 3/eʹ ) #[Hof06]2-113
const q₀ʹ= 3(1 + 1/eʹ²) * (1 - 1/eʹ * atan(eʹ) ) - 1 #[Hof06]2-133 with u = b and [Hof06]2-138
const m = ω_ie^2 * a^2 * b / GM #[Hof06] 2-70

#normal gravity magnitudes
const γ_a = GM / (a * b) * (1 - m - m/6 * eʹ * q₀ʹ/q₀) #Normal gravity at the equator, [Hof06] 2-141
const γ_b = GM / a² * (1 + m/3 * eʹ * q₀ʹ/q₀) #Normal gravity at the poles, [Hof06] 2-142

#distance tolerance (unit in last place of a)
const ε_a = eps(a)

abstract type WGS84Pos end

########################### NVectorAlt #####################################


struct NVectorAlt <: WGS84Pos
    n_e::SVector{3,Float64}
    h::Float64
    function NVectorAlt(n_e::AbstractVector{T} where {T<:Real}, h::Real; normalization::Bool = true)
        n_e = SVector{3,Float64}(n_e) #normalization will be faster for an SVector
        @assert h>-200 "Altitude must be larger than -200 m"
        return normalization ? new(normalize(n_e), h) : new(n_e, h)
    end
end

NVectorAlt(r::WGS84Pos) = convert(NVectorAlt, r)
NVectorAlt(input::Tuple{Union{Nothing, AbstractVector{T} where T<:Real}, Real}) = NVectorAlt(input...)
NVectorAlt(n_e::AbstractVector{T} where {T<:Real}) = NVectorAlt(n_e, 0)
NVectorAlt(; n_e = [1, 0, 0], h = 0) = NVectorAlt(n_e, h)

Base.:(==)(p1::NVectorAlt, p2::NVectorAlt) = (p1.n_e == p2.n_e && p1.h == p2.h)
Base.:(≈)(p1::NVectorAlt, p2::NVectorAlt) = (p1.n_e ≈ p2.n_e && p1.h ≈ p2.h)

#any fallback methods go here

#require each Rotation subtype to implement conversions to and from RQuat
Base.convert(::Type{NVectorAlt}, p::P) where {P<:WGS84Pos} = error("Implement $P to NVectorAlt conversion")
Base.convert(::Type{P}, p::NVectorAlt) where {P<:WGS84Pos} = error("Implement NVectorAlt to $P conversion")
#this enables a two-step conversion via NVectorAlt as a fallback
Base.convert(::Type{P}, p::WGS84Pos) where {P<:WGS84Pos} = convert(P, NVectorAlt(p))
#trivial conversions
Base.convert(::Type{P}, p::P) where {P<:WGS84Pos} = p
Base.convert(::Type{NVectorAlt}, p::NVectorAlt) = p

########################### Cartesian ##############################

struct Cartesian <: WGS84Pos
    r_OeP_e::SVector{3,Float64}
end

Cartesian(r::WGS84Pos) = convert(Cartesian, r)
Cartesian() = convert(Cartesian, NVectorAlt())

Base.:(==)(p1::Cartesian, p2::Cartesian) = p1.r_OeP_e == p2.r_OeP_e
Base.:(≈)(p1::Cartesian, p2::Cartesian) = p1.r_OeP_e ≈ p2.r_OeP_e

function Base.convert(::Type{Cartesian}, p::NVectorAlt)

    n_e = p.n_e; h = p.h
    _, N = radii(p)

    return Cartesian([
        (N + h) * n_e[1],
        (N + h) * n_e[2],
        (N * (1 - e²) + h) * n_e[3] ])

end

function Base.convert(::Type{NVectorAlt}, r::Cartesian)

    #see Fukushima: Transformation_from_Cartesian_to_Geodetic_Coordinates_Accelerated_by_Halley's_Method

    x, y, z = r.r_OeP_e
    p = √(x^2 + y^2)

    c = a*e²
    ec² = 1 - e²
    ec = √ec²
    zc = ec * abs(z)

    s0 = abs(z)
    c0 = ec*p
    a0 = √(s0^2 + c0^2)
    a0³ = a0^3
    b0 = 1.5 * c * s0 * c0 * ((p * s0 - zc * c0) * a0 - c * s0 * c0)
    s1 = (zc * a0³ + c * s0^3) * a0³ - b0 * s0
    c1 = (p * a0³ - c * c0^3) * a0³ - b0 * c0

    cc = ec*c1
    s1² = s1^2
    cc² = cc^2
    h = (p*cc + s0*s1 - a*√(ec² * s1² + cc²)) / √(s1² + cc²)

    if s1 < cc #|ϕ < π/4|
        abs_tan_ϕ = s1 / cc
        cos_ϕ = 1 / √(1 + abs_tan_ϕ^2) #cos_ϕ >=0
        abs_sin_ϕ = abs_tan_ϕ * cos_ϕ
        sin_ϕ = abs_sin_ϕ * sign(z)
    else #|ϕ > π/4|
        abs_cot_ϕ = cc / s1
        abs_sin_ϕ = 1 / √(1 + abs_cot_ϕ^2)
        cos_ϕ = abs_cot_ϕ * abs_sin_ϕ
        sin_ϕ = abs_sin_ϕ * sign(z)
    end

    cos_λ = p > 0 ? (x / p) : 1
    sin_λ = p > 0 ? (y / p) : 0

    n_e = [cos_ϕ*cos_λ, cos_ϕ*sin_λ, sin_ϕ]

    return NVectorAlt(n_e, h)

end



function radii(p::NVectorAlt)

    f_den = √(1 - e² * p.n_e[3]^2)
    return ( M = a * (1 - e²) / f_den^3, N = a / f_den ) #(R_N, R_E)

end


end