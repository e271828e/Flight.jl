module Geodesy

# using Base: Real, Symbol
using LinearAlgebra
using StaticArrays
using SHA
using UnPack
using Interpolations

using Flight.Attitude
using Flight.Plotting

export NVector, LatLon, Altitude, Ellipsoidal, Orthometric, Geopotential
export AltEllip, AltOrth, AltGeop
export Abstract3DPosition, Geographic, CartECEF
export ω_ie, gravity, g_l, G_l, ltf, radii, get_ψ_nl, get_geoid_offset

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



######################### Abstract2DLocation ###########################

abstract type Abstract2DLocation end

Base.convert(::Type{L}, loc::L) where {L<:Abstract2DLocation} = loc


#### NVector ####

Base.@kwdef struct NVector <: Abstract2DLocation
    data::SVector{3,Float64} = SVector{3}(1.0, 0.0, 0.0)
    function NVector(data::AbstractVector{T} where {T<:Real}; normalization::Bool = true)
        data = SVector{3,Float64}(data) #normalization will be faster for an SVector
        return normalization ? new(normalize(data)) : new(data)
    end
end

NVector(loc::Abstract2DLocation) = convert(NVector, loc)

#NVector from ECEF to Local Tangent Frame rotation
NVector(r_el::Rotation) = NVector(RMatrix(r_el))

NVector(r_el::RMatrix) = NVector( -r_el[:,3], normalization = false)

function NVector(r_el::RQuat)
    #n_e is simply the third column of the R_el rotation matrix. we don't need
    #the complete RQuat to RMatrix conversion
    q = r_el[:]
    dq12 = 2*q[1]*q[2]; dq13 = 2*q[1]*q[3]
    dq24 = 2*q[2]*q[4]; dq34 = 2*q[3]*q[4]
    NVector(-SVector{3}(dq24 + dq13, dq34 - dq12, 1 - 2*(q[2]^2 + q[3]^2)), normalization = false)
end

Base.:(==)(n1::NVector, n2::NVector) = n1.data == n2.data
Base.:(≈)(n1::NVector, n2::NVector) = n1.data ≈ n2.data
Base.:(-)(n::NVector) = NVector(-n.data)

#### AbstractArray interface
Base.size(::NVector) = (3,)
Base.length(::NVector) = 3
Base.getindex(n::NVector, i) = getindex(n.data, i)
#this allocates
# Base.iterate(n::NVector, state = 1) = (state > 3 ? nothing : (n.data[state], state + 1))


#### LatLon ####

Base.@kwdef struct LatLon <: Abstract2DLocation
    ϕ::Float64 = 0.0
    λ::Float64 = 0.0
    function LatLon(ϕ::Real, λ::Real)
        abs(ϕ) <= 0.5π || throw(ArgumentError("Latitude must be within [-π/2, π/2]"))
        abs(λ) <= π || throw(ArgumentError("Longitude must be within [-π, π]"))
        return new(ϕ, λ)
    end
end
LatLon(loc::Abstract2DLocation) = convert(LatLon, loc)

function Base.convert(::Type{NVector}, loc::LatLon)
    @unpack ϕ, λ = loc
    cos_ϕ = cos(ϕ)
    NVector(SVector{3,Float64}(cos_ϕ * cos(λ), cos_ϕ * sin(λ), sin(ϕ)), normalization = false)
end

function Base.convert(::Type{LatLon}, n_e::NVector)
    LatLon( ϕ = atan(n_e[3], √(n_e[1]^2 + n_e[2]^2)),
            λ = atan(n_e[2], n_e[1]))
end

#strict equality not implemented because comparison requires conversion
Base.:(≈)(loc1::LatLon, loc2::LatLon) = NVector(loc1) ≈ NVector(loc2)
Base.:(-)(loc::LatLon) = LatLon(-NVector(loc))


##### Generic Abstract2DLocation methods #####

#ellipsoid's radii of curvature
function radii(loc::Abstract2DLocation)
    n_e = NVector(loc)
    f_den = √(1 - e² * n_e[3]^2)
    return ( M = a * (1 - e²) / f_den^3, N = a / f_den ) #(R_N, R_E)
end

#local tangent frame from Abstract2DLocation
function ltf(loc::Abstract2DLocation, ψ_nl::Real = 0.)
    @unpack ϕ, λ = LatLon(loc)
    Rz(λ) ∘ Ry(-(ϕ + 0.5π)) ∘ Rz(ψ_nl)
end

#extract Wander Angle from ECEF to Local Tangent Frame rotation
get_ψ_nl(r_el::Rotation) = get_ψ_nl(RMatrix(r_el))
get_ψ_nl(r_el::RMatrix) = atan( -r_el[3,2], r_el[3,1] )
function get_ψ_nl(r_el::RQuat)
    #n_e is the third column of the R_el matrix. we don't need the complete
    #RQuat to RMatrix conversion
    q = r_el[:]
    dq12 = 2*q[1]*q[2]; dq13 = 2*q[1]*q[3]
    dq24 = 2*q[2]*q[4]; dq34 = 2*q[3]*q[4]
    atan(-(dq34 + dq12), dq24 - dq13)
end


#### AbstractAltitude ####

abstract type AbstractAltitudeDatum end
struct Ellipsoidal <: AbstractAltitudeDatum end
struct Orthometric <: AbstractAltitudeDatum end
struct Geopotential <: AbstractAltitudeDatum end

const h_min = -1000 #catches numerical catastrophes

function load_geoid_offset_interp(file_path = "src/ww15mgh_le.bin")
    #the target file stores a 721x1441 Matrix{Float32} in low-endian binary
    #format. the matrix holds the data points for the EGM96 geoid height offset
    #with respect to the WGS84 ellipsoid in 15 arc-minute resolution. latitude
    #goes from -π/2 to π/2, longitude from 0 to 2π
    tmp = Matrix{Float32}(undef, 721, 1441)
    open(file_path) do file
        if file |> sha2_256 |> bytes2hex != "9d190e021672769b508547021bcaebcc7d13558d66d215d019675a5f595f5cae"
            throw(ArgumentError("Wrong file hash"))
        else
            seekstart(file) #the hash check has moved the IOStream to EOF, reset position
            read!(file, tmp)
        end
    end
    #convert matrix elements to host's endianness and cast to Matrix{Float64}
    data = convert(Matrix{Float64}, ltoh.(tmp))
    ϕ_range = LinRange(-π/2, π/2, size(data, 1))
    λ_range = LinRange(0, 2π, size(data, 2))

    #return interpolator with extrapolation enabled to avoid machine precision
    #issues due to the boundaries being multiples of π
    LinearInterpolation((ϕ_range, λ_range), data, extrapolation_bc = Line())
    # CubicSplineInterpolation((ϕ_range, λ_range), data, extrapolation_bc = Line())
end

const geoid_offset_interp = load_geoid_offset_interp()

function get_geoid_offset(loc::Abstract2DLocation)
    #our longitude interval is [-π,π], but the table uses [0,2π], so we need to
    #correct for that
    latlon = LatLon(loc)
    ϕ = latlon.ϕ
    λ = mod(latlon.λ + 2π, 2π)
    geoid_offset_interp(ϕ, λ)
end

Base.@kwdef struct Altitude{D<:AbstractAltitudeDatum}
    _val::Float64 = 0.0
    #this constructor prevents the user from passing an Altitude with a
    #different datum, which would get implicitly converted to Float64 without
    #actually changing its value
    function Altitude{D}(h::Real) where {D}
        h >= h_min || throw(ArgumentError("Minimum altitude value is $h_min m, got $h"))
        return new{D}(h)
    end
end

const AltEllip = Altitude{Ellipsoidal}
const AltOrth = Altitude{Orthometric}
const AltGeop = Altitude{Geopotential}

Altitude{D}(h::Altitude{D}, args...) where {D} = Altitude{D}(h._val)

#Ellipsoidal and Orthometric altitudes are related by the geoid's offset at a
#given 2D location
Altitude{Ellipsoidal}(h::AltOrth, loc::Abstract2DLocation) = AltEllip(h._val + get_geoid_offset(loc))
Altitude{Orthometric}(h::AltEllip, loc::Abstract2DLocation) = AltOrth(h._val - get_geoid_offset(loc))

#Orthometric and Geopotential are directly related by the point-mass gravity
#approximation, 2D location not required
Altitude{Geopotential}(h::AltOrth) = AltGeop(h*a / (a+h))
Altitude{Orthometric}(h::AltGeop) = AltOrth(h*a / (a-h))

#still, for interface consistency we provide the two-argument method
Altitude{Geopotential}(h::AltOrth, ::Abstract2DLocation) = AltGeop(h)
Altitude{Orthometric}(h::AltGeop, ::Abstract2DLocation) = AltOrth(h)

#Geopotential and Ellipsoidal altitudes are related via Orthometric
Altitude{Geopotential}(h_ellip::AltEllip, loc::Abstract2DLocation) = AltGeop(AltOrth(h_ellip, loc))
Altitude{Ellipsoidal}(h_geop::AltGeop, loc::Abstract2DLocation) = AltEllip(AltOrth(h_geop), loc)

#operations between Altitude subtypes and Reals
Base.promote_rule(::Type{<:Altitude{D}}, ::Type{<:Real}) where {D} = Altitude{D}
Base.convert(::Type{<:Altitude{D}}, h::Real) where {D} = Altitude{D}(h)
Base.convert(::Type{T}, h::Altitude) where {T<:Real} = convert(T, h._val)
(::Type{<:T})(h::Altitude) where {T<:Real} = convert(T, h)

Base.:+(h::Altitude, Δh::Real) = h._val + Δh
Base.:-(h::Altitude, Δh::Real) = h._val - Δh
Base.:+(Δh::Real, h::Altitude) = Δh + h._val
Base.:-(Δh::Real, h::Altitude) = Δh - h._val

Base.:*(k::Real, h::Altitude) = k*h._val
Base.:*(h::Altitude, k::Real) = k*h
Base.:/(h::Altitude, k::Real) = h._val/k
Base.:-(h1::Altitude{D}, h2::Altitude{D}) where {D} = h1._val - h2._val

Base.:(==)(h1::Altitude{D}, h2::Altitude{D}) where {D} = h1._val == h2._val
Base.:≈(h1::Altitude{D}, h2::Altitude{D}) where {D} = h1._val ≈ h2._val
Base.:>(h1::Altitude{D}, h2::Altitude{D}) where {D} = h1._val > h2._val
Base.:<(h1::Altitude{D}, h2::Altitude{D}) where {D} = h1._val < h2._val

Base.:(==)(h1::Altitude{D}, h2::Real) where {D} = ==(promote(h1,h2)...)
Base.:≈(h1::Altitude{D}, h2::Real) where {D} = ≈(promote(h1,h2)...)
Base.:>(h1::Altitude{D}, h2::Real) where {D} = >(promote(h1,h2)...)
Base.:<(h1::Altitude{D}, h2::Real) where {D} = <(promote(h1,h2)...)

Base.:(==)(h1::Real, h2::Altitude{D}) where {D} = h2 == h1
Base.:≈(h1::Real, h2::Altitude{D}) where {D} = h2 ≈ h1
Base.:>(h1::Real, h2::Altitude{D}) where {D} = h2 < h1
Base.:<(h1::Real, h2::Altitude{D}) where {D} = h2 > h1


########################## Abstract3DPosition ##########################

abstract type Abstract3DPosition end

#avoid infinite recursion
Base.convert(::Type{P}, p::P) where {P<:Abstract3DPosition} = p


#### Geographic ####

#the default constructor generates a Geographic{NVector, EllipsoidalAlt} instance
Base.@kwdef struct Geographic{L <: Abstract2DLocation, H <: AbstractAltitudeDatum} <: Abstract3DPosition
    loc::L = NVector()
    alt::Altitude{H} = AltOrth()
end
Geographic(loc::Abstract2DLocation) = Geographic(loc, AltOrth())
Geographic(p::Abstract3DPosition) = Geographic{NVector,Ellipsoidal}(p)
Geographic{L,H}(p::Abstract3DPosition) where {L,H} = convert(Geographic{L,H}, p)

function Base.convert(::Type{Geographic{L,H}}, p::Geographic) where {L,H}
    Geographic(convert(L, p.loc), Altitude{H}(p))
end

Altitude{D}(p::Abstract3DPosition) where {D} = Altitude{D}(Geographic{NVector,D}(p))
Altitude{D}(p::Geographic) where {D} = Altitude{D}(p.alt, p.loc)

function Base.:(==)(p1::Geographic{NVector,H}, p2::Geographic{NVector,H}) where {H}
    return p1.alt == p2.alt && p1.loc == p2.loc
end

function Base.:(≈)(p1::Geographic{L,H}, p2::Geographic{L,H}) where {L,H}
    return p1.alt ≈ p2.alt && p1.loc ≈ p2.loc
end

Base.:(≈)(p1::Abstract3DPosition, p2::Abstract3DPosition) = CartECEF(p1) ≈ CartECEF(p2)

function Base.:(==)(p1::Abstract3DPosition, p2::Abstract3DPosition)
    throw(ArgumentError("Exact comparison between $(typeof(p1)) and $(typeof(p2)) not defined, use ≈ instead"))
end

Base.:(-)(p::T) where {T<:Abstract3DPosition} = convert(T, -CartECEF(p))



############################# CartECEF #############################

struct CartECEF <: Abstract3DPosition
    data::SVector{3,Float64}
end
CartECEF(p::Abstract3DPosition) = convert(CartECEF, p)
CartECEF() = CartECEF(Geographic())

Base.:(==)(r1::CartECEF, r2::CartECEF) = r1.data == r2.data
Base.:(≈)(r1::CartECEF, r2::CartECEF) = r1.data ≈ r2.data
Base.:(-)(r::CartECEF) = CartECEF(-r.data)

#### AbstractArray interface
Base.size(::CartECEF) = (3,)
Base.getindex(n::CartECEF, i) = getindex(n.data, i)

function Base.convert(::Type{Geographic{L,H}}, r::CartECEF) where {L,H}
    convert(Geographic{L,H}, convert(Geographic{NVector,Ellipsoidal}, r))
end

function Base.convert(::Type{Geographic{NVector,Ellipsoidal}}, r::CartECEF)

    #NVector + Alt from ECEF Cartesian position vector. See Fukushima:
    #Transformation_from_Cartesian_to_Geodetic_Coordinates_Accelerated_by_Halley's_Method

    x, y, z = r[:]
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

    return Geographic(
        NVector(SVector{3,Float64}(cos_ϕ*cos_λ, cos_ϕ*sin_λ, sin_ϕ)),
        Altitude{Ellipsoidal}(h))

end

function Base.convert(::Type{CartECEF}, p::Geographic)
    convert(CartECEF, convert(Geographic{NVector,Ellipsoidal}, p))
end

function Base.convert(::Type{CartECEF}, p::Geographic{NVector, Ellipsoidal})

    n_e = p.loc; h = Float64(p.alt)
    _, N = radii(n_e)

    return CartECEF(SVector{3, Float64}(
        (N + h) * n_e[1],
        (N + h) * n_e[2],
        (N * (1 - e²) + h) * n_e[3]))

end

##### Generic Abstract3DPosition methods ####

ltf(p::Abstract3DPosition, ψ_nl::Real = 0.0) = ltf(Geographic{NVector,Ellipsoidal}(p).loc, ψ_nl)
radii(p::Abstract3DPosition) = radii(Geographic{NVector,Ellipsoidal}(p).loc)

"""
    gravity(p::Abstract3DPosition)

Compute normal gravity.

Computation is based on Somigliana's formula for gravity at the ellipsoid
surface, with a second order altitude correction, accurate for small altitudes
above the WGS84 ellipsoid (h<<a). Direction is assumed normal to the WGS84
ellipsoid, a good enough approximation for most navigation applications. See
Hoffmann & Moritz.
"""
function gravity(p::Abstract3DPosition)

    p_nve = Geographic{NVector,Ellipsoidal}(p)
    n_e = p_nve.loc
    h = Float64(p_nve.alt)

    sin²ϕ = n_e[3]^2
    cos²ϕ = n_e[1]^2 + n_e[2]^2

    #gravity at the ellipsoid surface (Somigliana)
    γ_0 = (a * γ_a * cos²ϕ + b * γ_b * sin²ϕ) / √(a² * cos²ϕ + b² * sin²ϕ) #[Hof06] 2-146

    #altitude correction
    γ = γ_0 * (1 - 2/a * (1 + f + m - 2f * sin²ϕ) * h + 3/a² * h^2)

    return γ

end

"""
    g_l(p::Abstract3DPosition)

Compute gravity vector resolved in the local tangent frame.
"""
g_l(p::Abstract3DPosition) = SVector{3}(0, 0, gravity(p))

"""
    G_l(p::Abstract3DPosition)

Compute gravitational attraction resolved in the local tangent frame.
"""
function G_l(p::Abstract3DPosition)

    q_el = ltf(p)
    ω_ie_e = SVector{3, Float64}(0,0,ω_ie)
    r_eP_e = CartECEF(p)[:]
    G_l = g_l(p) + q_el'(ω_ie_e × (ω_ie_e × r_eP_e))
    return G_l

end

##################### Plotting ##############################

@recipe function plot_latlon(th::TimeHistory{<:AbstractVector{<:LatLon}})

    sa = StructArray(th.data)
    data = hcat(sa.ϕ/π, sa.λ/π)

    title --> ["Latitude" "Longitude"]
    label --> ["Latitude" "Longitude"]
    yguide --> [L"$\varphi \ (\pi \ rad)$" L"$\lambda \ (\pi \ rad)$"]
    th_split --> :v

    return TimeHistory(th.t, data)

end

@recipe function plot_altitude(th::TimeHistory{<:AbstractVector{<:Altitude{D}}}) where {D}

    title --> "Altitude ($String(D))"
    label --> "Altitude ($String(D))"
    yguide --> L"$h \ (m)$"

    return TimeHistory(th.t, StructArray(th.data)._val)

end

end #module