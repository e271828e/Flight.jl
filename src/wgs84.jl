module WGS84

using LinearAlgebra
using StaticArrays: SVector
using ..Attitude

export WGS84Pos, NVectorAlt, Cartesian

#WGS84 fundamental constants, SI units
const GM = 3.986005e+14 #Gravitational constant
const a = 6378137 #Equatorial radius
const f = 1/298.257223563 #Ellipsoid flattening
const ω_ie = 7.292115e-05 #Earth's angular velocity with respect to the ECI frame

#derived parameters
const b = a * (1 - f) #Polar semi-minor axis
const e² = 2f - f^2 #First eccentricity squared (^ operator calls _power_by_squaring, no need to write f*f)
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

function Base.convert(::Type{NVectorAlt}, p::Cartesian)

    #a let block explicitly tells the compiler that the r_z and l_xy inside
    #these functions will never be changed, because they are new bindings which
    #live only within the block. but since the closure is not returned from this
    #function, so apparently no let block is required to stabilize the types

    ε_n = 1e-10
    ε_h = 1e-7
    max_iter = 10

    r_OeP_e = p.r_OeP_e
    r_x, r_y, r_z = r_OeP_e[1], r_OeP_e[2], r_OeP_e[3]
    l_xy = norm(r_OeP_e[1:2])

    #if abs(φ) closer to 0, use tan_φ as iteration variable
    function step_tan(tan_φ::Float64)::NTuple{3,Float64}
        cos_φ = 1 / √(1 + tan_φ^2) #positive within [-π/2, π/2]
        N = a / √(1 - e² * (1 - cos_φ^2))
        h = l_xy / cos_φ - N
        tan_φ = r_z / (l_xy * (1 - e² * N / (N + h)))
        return tan_φ, h, N
    end

    #if abs(φ) closer to π/2, use cot_φ as iteration variable
    function step_cot(cot_φ::Float64)::NTuple{3,Float64}
        sin_φ = sign(cot_φ) / √(1 + cot_φ^2) #within [-π/2, π/2] sign(sin_φ) = sign(cot_φ)
        N = a / √(1 - e² * sin_φ^2)
        h = r_z / sin_φ - (1 - e²) * N
        cot_φ = l_xy * (1 - e² * N / (N + h)) / r_z
        return cot_φ, h, N
    end

    h0 = 0.0
    if l_xy > abs(r_z)
        x0 = r_z / (l_xy * (1 - e²))
        step = step_tan
    else
        x0 = l_xy * (1 - e²) / r_z
        step = step_cot
    end

    for _ in 1:max_iter
        x, h, N = step(x0)
        if abs(x - x0) < ε_n && abs(h - h0) < ε_h
            n_e = [ r_x / (N + h),
                    r_y / (N + h),
                    r_z / ((1 - e²) * N + h)]
            return NVectorAlt(n_e, h)
        end
        x0 = x; h0 = h
    end
    error("Iteration terminated without convergence")

end



function radii(p::NVectorAlt)

    f_den = √(1 - e² * p.n_e[3]^2)
    return ( M = a * (1 - e²) / f_den^3, N = a / f_den ) #(R_N, R_E)

end


end