module Geodesy

# using Base: Real, Symbol
using LinearAlgebra
using StaticArrays
using SHA
using UnPack
using Interpolations
using HDF5
using Plots

using Flight.Utils
using Flight.Systems
using Flight.Plotting
using Flight.Attitude

import Flight.Plotting: make_plots

export Abstract2DLocation, NVector, LatLon
export Altitude, Ellipsoidal, Orthometric, Geopotential, AltE, AltO, AltG
export Abstract3DLocation, GeographicLocation, CartesianLocation
export ω_ie, gravity, g_n, G_n, ltf, radii, get_ψ_nl, get_geoid_offset

#WGS84 fundamental constants, SI units
const GM = 3.986005e+14 #Gravitational constant
const a = 6378137.0 #Equatorial radius
const f = 1/298.257223563 #Ellipsoid flattening
const ω_ie = 7.292115e-05 #Earth's angular velocity with respect to the ECI frame

#derived parameters
const b = a * (1 - f) #Polar semi-minor axis
const e² = 2f - f^2 #First eccentricity squared
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

Base.convert(::Type{L}, l2d::L) where {L<:Abstract2DLocation} = l2d


#### NVector ####

Base.@kwdef struct NVector <: Abstract2DLocation
    data::SVector{3,Float64} = SVector{3}(1.0, 0.0, 0.0)
    function NVector(data::AbstractVector{T} where {T<:Real}; normalization::Bool = true)
        data = SVector{3,Float64}(data) #normalization will be faster for an SVector
        return normalization ? new(normalize(data)) : new(data)
    end
end

NVector(l2d::Abstract2DLocation) = convert(NVector, l2d)

#NVector from ECEF to Local Tangent Frame rotation
NVector(r_el::Abstract3DRotation) = NVector(RMatrix(r_el))

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
LinearAlgebra.norm(n::NVector) = norm(getfield(n, :data)) #uses StaticArrays implementation
LinearAlgebra.normalize(n::NVector) = NVector(getfield(n, :data)) #let the constructor normalize


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
LatLon(l2d::Abstract2DLocation) = convert(LatLon, l2d)

function Base.convert(::Type{NVector}, l2d::LatLon)
    @unpack ϕ, λ = l2d
    cos_ϕ = cos(ϕ)
    NVector(SVector{3,Float64}(cos_ϕ * cos(λ), cos_ϕ * sin(λ), sin(ϕ)), normalization = false)
end

function Base.convert(::Type{LatLon}, n_e::NVector)
    LatLon( ϕ = atan(n_e[3], √(n_e[1]^2 + n_e[2]^2)),
            λ = atan(n_e[2], n_e[1]))
end

#strict equality not implemented because comparison requires conversion
Base.:(≈)(ll1::LatLon, ll2::LatLon) = NVector(ll1) ≈ NVector(ll2)
Base.:(-)(latlon::LatLon) = LatLon(-NVector(latlon))


##### Generic Abstract2DLocation methods #####

#ellipsoid's radii of curvature
function radii(l2d::Abstract2DLocation)
    n_e = NVector(l2d)
    f_den = √(1 - e² * n_e[3]^2)
    return ( M = a * (1 - e²) / f_den^3, N = a / f_den ) #(R_N, R_E)
end

#local tangent frame from Abstract2DLocation
function ltf(l2d::Abstract2DLocation, ψ_nl::Real = 0.)
    @unpack ϕ, λ = LatLon(l2d)
    Rz(λ) ∘ Ry(-(ϕ + 0.5π)) ∘ Rz(ψ_nl)
end

#extract Wander Angle from ECEF to Local Tangent Frame rotation
get_ψ_nl(r_el::Abstract3DRotation) = get_ψ_nl(RMatrix(r_el))
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
abstract type AbstractGeometricAltitudeDatum <: AbstractAltitudeDatum end
struct Ellipsoidal <: AbstractGeometricAltitudeDatum end
struct Orthometric <: AbstractGeometricAltitudeDatum end
struct Geopotential <: AbstractAltitudeDatum end

const h_min = -1000 #helps catch numerical catastrophes

#functionally equivalent to load_geoid_data_hdf5 but with the additional hash check
function load_geoid_data_bin(file_path = "src/ww15mgh_le.bin")
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

function load_geoid_data_hdf5(file_path = "src/common/ww15mgh_hdf5.h5")
    data = Matrix{Float32}(undef, 721, 1441)
    h5open(file_path) do file
        data .= file["geoid_height"] |> read
    end
    ϕ_range = LinRange(-π/2, π/2, size(data, 1))
    λ_range = LinRange(0, 2π, size(data, 2))
    LinearInterpolation((ϕ_range, λ_range), data, extrapolation_bc = Line())
end

const geoid_data = load_geoid_data_hdf5()

function get_geoid_offset(l2d::Abstract2DLocation)
    #our longitude interval is [-π,π], but the table uses [0,2π], so we need to
    #correct for that
    latlon = LatLon(l2d)
    ϕ = latlon.ϕ
    λ = mod(latlon.λ + 2π, 2π)
    geoid_data(ϕ, λ)
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

const AltE = Altitude{Ellipsoidal}
const AltO = Altitude{Orthometric}
const AltG = Altitude{Geopotential}

Altitude{D}(h::Altitude{D}, args...) where {D} = Altitude{D}(h._val)

#Ellipsoidal and Orthometric altitudes are related by the geoid's offset at a
#given 2D location
Altitude{Ellipsoidal}(h::AltO, l2d::Abstract2DLocation) = AltE(h._val + get_geoid_offset(l2d))
Altitude{Orthometric}(h::AltE, l2d::Abstract2DLocation) = AltO(h._val - get_geoid_offset(l2d))

#Orthometric and Geopotential are directly related by the point-mass gravity
#approximation, 2D location not required
Altitude{Geopotential}(h::AltO) = AltG(h._val*a / (a+h._val))
Altitude{Orthometric}(h::AltG) = AltO(h._val*a / (a-h._val))

#still, for interface consistency we provide the two-argument method
Altitude{Geopotential}(h::AltO, ::Abstract2DLocation) = AltG(h)
Altitude{Orthometric}(h::AltG, ::Abstract2DLocation) = AltO(h)

#Geopotential and Ellipsoidal altitudes are related via Orthometric
Altitude{Geopotential}(h_ellip::AltE, l2d::Abstract2DLocation) = AltG(AltO(h_ellip, l2d))
Altitude{Ellipsoidal}(h_geop::AltG, l2d::Abstract2DLocation) = AltE(AltO(h_geop), l2d)

#operations between Altitude subtypes and Reals
Base.promote_rule(::Type{<:Altitude{D}}, ::Type{<:Real}) where {D} = Altitude{D}
Base.convert(::Type{<:Altitude{D}}, h::Real) where {D} = Altitude{D}(h)
Base.convert(::Type{T}, h::Altitude) where {T<:Real} = convert(T, h._val)
(::Type{<:T})(h::Altitude) where {T<:Real} = convert(T, h)

Base.:+(h::T, Δh::Real) where {T<:Altitude} = T(h._val + Δh)
Base.:-(h::T, Δh::Real) where {T<:Altitude} = T(h._val - Δh)
Base.:+(Δh::Real, h::T) where {T<:Altitude} = T(Δh + h._val)
Base.:-(Δh::Real, h::T) where {T<:Altitude} = T(Δh - h._val)

Base.:*(k::Real, h::T) where {T<:Altitude} = T(k*h._val)
Base.:*(h::T, k::Real) where {T<:Altitude} = T(k*h)
Base.:/(h::T, k::Real) where {T<:Altitude} = T(h._val/k)
Base.:^(h::T, k::Real) where {T<:Altitude} = h._val^k

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


########################## Abstract3DLocation ##########################

abstract type Abstract3DLocation end

#avoid infinite recursion
Base.convert(::Type{P}, p::P) where {P<:Abstract3DLocation} = p

########################### GeographicLocation ###############################

#the default constructor generates a GeographicLocation{NVector, AltOrthometric} instance
Base.@kwdef struct GeographicLocation{L <: Abstract2DLocation, H <: AbstractAltitudeDatum} <: Abstract3DLocation
    l2d::L = NVector()
    alt::Altitude{H} = AltO()
end
GeographicLocation(l2d::Abstract2DLocation) = GeographicLocation(l2d, AltO())
GeographicLocation(l3d::Abstract3DLocation) = GeographicLocation{NVector,Ellipsoidal}(l3d)
GeographicLocation{L,H}(l3d::Abstract3DLocation) where {L,H} = convert(GeographicLocation{L,H}, l3d)

function Base.convert(::Type{GeographicLocation{L,H}}, geo::GeographicLocation) where {L,H}
    GeographicLocation(convert(L, geo.l2d), Altitude{H}(geo))
end

NVector(geo::GeographicLocation) = NVector(geo.l2d)
LatLon(geo::GeographicLocation) = LatLon(geo.l2d)
Altitude{D}(geo::GeographicLocation) where {D} = Altitude{D}(geo.alt, geo.l2d)

function Base.:(==)(geo1::GeographicLocation{NVector,H}, geo2::GeographicLocation{NVector,H}) where {H}
    return geo1.alt == geo2.alt && geo1.l2d == geo2.l2d
end

function Base.:(≈)(geo1::GeographicLocation{L,H}, geo2::GeographicLocation{L,H}) where {L,H}
    return geo1.alt ≈ geo2.alt && geo1.l2d ≈ geo2.l2d
end

Base.:(≈)(loc1::Abstract3DLocation, loc2::Abstract3DLocation) = CartesianLocation(loc1) ≈ CartesianLocation(loc2)

function Base.:(==)(loc1::Abstract3DLocation, loc2::Abstract3DLocation)
    throw(ArgumentError("Exact comparison between $(typeof(loc1)) and $(typeof(loc2)) not defined, use ≈ instead"))
end

Base.:(-)(l3d::T) where {T<:Abstract3DLocation} = convert(T, -CartesianLocation(l3d))


############################# CartesianLocation #############################

struct CartesianLocation <: Abstract3DLocation
    data::SVector{3,Float64}
end
CartesianLocation(l3d::Abstract3DLocation) = convert(CartesianLocation, l3d)
CartesianLocation() = CartesianLocation(GeographicLocation())

NVector(r::CartesianLocation) = GeographicLocation{NVector, Ellipsoidal}(r).l2d
LatLon(r::CartesianLocation) = GeographicLocation{LatLon, Ellipsoidal}(r).l2d
Altitude{D}(r::CartesianLocation) where {D} = Altitude{D}(GeographicLocation{NVector,D}(r))

Base.:(==)(r1::CartesianLocation, r2::CartesianLocation) = r1.data == r2.data
Base.:(≈)(r1::CartesianLocation, r2::CartesianLocation) = r1.data ≈ r2.data
Base.:(-)(r::CartesianLocation) = CartesianLocation(-r.data)
Base.:(+)(r1::CartesianLocation, r2::AbstractVector{<:Real}) = CartesianLocation(r1.data + SVector{3,Float64}(r2))
Base.:(+)(r1::AbstractVector{<:Real}, r2::CartesianLocation) = r2 + r1

#### AbstractArray interface
Base.size(::CartesianLocation) = (3,)
Base.getindex(n::CartesianLocation, i) = getindex(n.data, i)

function Base.convert(::Type{GeographicLocation{L,H}}, r::CartesianLocation) where {L,H}
    convert(GeographicLocation{L,H}, convert(GeographicLocation{NVector,Ellipsoidal}, r))
end

function Base.convert(::Type{GeographicLocation{NVector,Ellipsoidal}}, r::CartesianLocation)

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

    return GeographicLocation(
        NVector(SVector{3,Float64}(cos_ϕ*cos_λ, cos_ϕ*sin_λ, sin_ϕ)),
        Altitude{Ellipsoidal}(h))

end

function Base.convert(::Type{CartesianLocation}, geo::GeographicLocation)
    convert(CartesianLocation, convert(GeographicLocation{NVector,Ellipsoidal}, geo))
end

function Base.convert(::Type{CartesianLocation}, geo::GeographicLocation{NVector, Ellipsoidal})

    n_e = geo.l2d; h = Float64(geo.alt)
    _, N = radii(n_e)

    return CartesianLocation(SVector{3, Float64}(
        (N + h) * n_e[1],
        (N + h) * n_e[2],
        (N * (1 - e²) + h) * n_e[3]))

end

##### Generic Abstract3DLocation methods ####

#general conversion from 3D to 2D location
(::Type{L})(l3d::Abstract3DLocation) where {L<:Abstract2DLocation} = GeographicLocation{L,Ellipsoidal}(l3d).l2d

ltf(l3d::Abstract3DLocation, ψ_nl::Real = 0.0) = ltf(GeographicLocation{NVector,Ellipsoidal}(l3d).l2d, ψ_nl)
radii(l3d::Abstract3DLocation) = radii(GeographicLocation{NVector,Ellipsoidal}(l3d).l2d)

"""
    gravity(p::Abstract3DLocation)

Compute normal gravity.

Computation is based on Somigliana's formula for gravity at the ellipsoid
surface, with a second order altitude correction, accurate for small altitudes
above the WGS84 ellipsoid (h<<a). Direction is assumed normal to the WGS84
ellipsoid, a good enough approximation for most navigation applications. See
Hoffmann & Moritz.
"""
function gravity(l3d::Abstract3DLocation)

    p_nve = GeographicLocation{NVector,Ellipsoidal}(l3d)
    n_e = p_nve.l2d
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
    g_n(p::Abstract3DLocation)

Compute gravity vector resolved in the local tangent frame.
"""
g_n(l3d::Abstract3DLocation) = SVector{3}(0, 0, gravity(l3d))

"""
    G_n(p::Abstract3DLocation)

Compute gravitational attraction resolved in the local tangent frame.
"""
function G_n(l3d::Abstract3DLocation)

    q_en = ltf(l3d)
    ω_ie_e = SVector{3, Float64}(0,0,ω_ie)
    r_eP_e = CartesianLocation(l3d)[:]
    G_n = g_n(l3d) + q_en'(ω_ie_e × (ω_ie_e × r_eP_e))
    return G_n

end

############################### Plotting #######################################

@recipe function fp(th::TimeHistory{<:LatLon})

    title --> ["Latitude" "Longitude"]
    label --> ["Latitude" "Longitude"]
    yguide --> [L"$\varphi \ (\pi \ rad)$" L"$\lambda \ (\pi \ rad)$"]
    th_split --> :v

    data = hcat(th.ϕ._data, th.λ._data)'/π |> collect
    return TimeHistory(th._t, data)

end

@recipe function fp(th::TimeHistory{<:Altitude{D}}) where {D}

    title --> "Altitude ($(string(D)))"
    label --> "Altitude ($(string(D)))"
    yguide --> L"$h \ (m)$"

    return TimeHistory(th._t, Float64.(th._data))

end


end #module