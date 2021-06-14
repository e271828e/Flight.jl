module WGS84

using Base: Real, Symbol
using LinearAlgebra
using StaticArrays: SVector
using Flight.Attitude

export NVector, WGS84Pos, Cartesian

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

########################### NVector ##################################

struct NVector <: AbstractVector{Float64}
    data::SVector{3,Float64}
    function NVector(data::AbstractVector{T} where {T<:Real}; normalization::Bool = true)
        data = SVector{3,Float64}(data) #normalization will be faster for an SVector
        return normalization ? new(normalize(data)) : new(data)
    end
end

function NVector(; ϕ::Real = 0., λ::Real = 0.)
    cos_ϕ = cos(ϕ)
    NVector([cos_ϕ * cos(λ), cos_ϕ * sin(λ), sin(ϕ)], normalization = false)
end

#extract NVector and Wander Angle from ECEF to Local Tangent Frame rotation
function NVector(r_el::RMatrix)
    n_e = NVector( -r_el[:,3], normalization = false)
    ψ_nl = atan( -r_el[3,2], r_el[3,1] )
    return (n_e = n_e, ψ_nl = ψ_nl)
end

function NVector(r_el::RQuat)
    #n_e is the third column of the R_el matrix. we don't need the complete
    #RQuat to RMatrix conversion
    q = r_el._quat
    dq12 = 2*q[1]*q[2]; dq13 = 2*q[1]*q[3]
    dq24 = 2*q[2]*q[4]; dq34 = 2*q[3]*q[4]

    n_e = NVector(-[dq24 + dq13, dq34 - dq12, 1 - 2*(q[2]^2 + q[3]^2)], normalization = false)
    ψ_nl = atan(-(dq34 + dq12), dq24 - dq13)
    return (n_e = n_e, ψ_nl = ψ_nl)
end

NVector(r_el::Rotation) = NVector(RQuat(r_el))

Base.:(==)(n1::NVector, n2::NVector) = n1.data == n2.data
Base.:(≈)(n1::NVector, n2::NVector) = n1.data ≈ n2.data

#AbstractArray
Base.size(x::NVector) = size(getfield(x,:data))
Base.getindex(x::NVector, i) = getindex(getfield(x,:data), i)

lat(n_e::NVector) = atan(n_e[3], √(n_e[1]^2 + n_e[2]^2))
lon(n_e::NVector) = atan(n_e[2], n_e[1])
ltf(n_e::NVector, ψ_nl::Real = 0.) = Rz( lon(n_e) ) ∘ Ry( -(lat(n_e) + 0.5π) ) ∘ Rz(ψ_nl)
ltf(nt::NamedTuple{(:n_e, :ψ_nl), Tuple{NVector, Float64}}) = ltf(nt[:n_e], nt[:ψ_nl])

function radii(n_e::NVector)
    f_den = √(1 - e² * n_e[3]^2)
    return ( M = a * (1 - e²) / f_den^3, N = a / f_den ) #(R_N, R_E)
end

#as long as s has an explicit value at compile time (which it will), the if will
#be optimized out by the compiler. no need to dispatch on Val(s) for efficiency.
Base.propertynames(::NVector, private::Bool = false) = (:lat, :lon, :ϕ, :λ, :ltf, :r_el)
function Base.getproperty(n_e::NVector, s::Symbol)
    if s == :lat || s == :ϕ
        return lat(n_e)
    elseif s == :lon || s == :λ
        return lon(n_e)
    elseif s == :ltf || s == :r_el
        return ltf(n_e)
    else
        error("NVector has no property $s")
    end
end

########################### WGS84Pos ##########################

const h_min = -20000

struct WGS84Pos
    n_e::NVector
    h::Float64
    function WGS84Pos(n_e::NVector, h::Real)
        h > h_min || throw(ArgumentError("Minimum altitude is $h_min m"))
        return new(n_e, h)
    end
end

WGS84Pos(; ϕ::Real = 0., λ::Real = 0., h::Real = 0.) = WGS84Pos(NVector(ϕ = ϕ, λ = λ), h)

function WGS84Pos(r::AbstractVector{<:Real})

    #NVector + Alt from ECEF Cartesian position vector. See Fukushima:
    #Transformation_from_Cartesian_to_Geodetic_Coordinates_Accelerated_by_Halley's_Method

    if length(r)!=3
        throw(ArgumentError("Constructor from Cartesian ECEF position expected a
        3-element vector, got length $(length(r))"))
    end

    x, y, z = r
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

    n_e = NVector([cos_ϕ*cos_λ, cos_ϕ*sin_λ, sin_ϕ])

    return WGS84Pos(n_e, h)

end

function cart_ECEF(p::WGS84Pos)

    n_e = p.n_e; h = p.h
    _, N = radii(p.n_e)

    return SVector{3, Float64}([
        (N + h) * n_e[1],
        (N + h) * n_e[2],
        (N * (1 - e²) + h) * n_e[3] ])

end

radii(p::WGS84Pos) = radii(p.n_e)

Base.:(==)(p1::WGS84Pos, p2::WGS84Pos) = (p1.n_e == p2.n_e && p1.h == p2.h)
Base.:(≈)(p1::WGS84Pos, p2::WGS84Pos) = (p1.n_e ≈ p2.n_e && p1.h ≈ p2.h)

Base.propertynames(p::WGS84Pos, private::Bool = false) = (:n_e, :h, :alt,
    propertynames(getfield(p, :n_e), private)...)

function Base.getproperty(p::WGS84Pos, s::Symbol)
    if s == :n_e
        return getfield(p, :n_e)
    elseif s == :h || s == :alt
        return getfield(p, :h)
    else #delegate to n_e
        return getproperty(getfield(p, :n_e), s)
    end
end


end