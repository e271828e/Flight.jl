
"""
Position with respect to the WGS84 Ellipsoid
"""
# module Geodesy

using Base: Real, Symbol
using LinearAlgebra
using SHA
using StaticArrays: SVector
using UnPack
using Interpolations

using Flight.Attitude
using Flight.Plotting

export NVector, LatLon, EllipsoidalAlt, OrthometricAlt, Geographic, CartECEF
export ω_ie, gravity, g_l, G_l, ltf, radii, get_ψ_nl

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

abstract type AbstractAltitude end

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

Base.@kwdef struct EllipsoidalAlt <: AbstractAltitude
    _val::Float64 = 0.0
end

Base.@kwdef struct OrthometricAlt <: AbstractAltitude
    _val::Float64 = 0.0
end

Base.convert(::Type{H}, h::Real) where {H<:AbstractAltitude} = H(h)
Base.convert(::Type{T}, h::EllipsoidalAlt) where {T<:Real} = T(h._val)
Base.convert(::Type{T}, h::OrthometricAlt) where {T<:Real} = T(h._val)
(::Type{<:T})(h::EllipsoidalAlt) where {T<:Real} = T(h._val)
(::Type{<:T})(h::OrthometricAlt) where {T<:Real} = T(h._val)

Base.:+(h::H, Δh::Real) where {H<:AbstractAltitude} = H(h._val + Δh)
Base.:-(h::H, Δh::Real) where {H<:AbstractAltitude} = H(h._val - Δh)

Base.:+(h1::H, h2::H) where {H<:AbstractAltitude} = H(h1._val + h2._val)
Base.:-(h1::H, h2::H) where {H<:AbstractAltitude} = h1._val - h2._val #an altitude difference has no reference
Base.:(==)(h1::H, h2::H) where {H<:AbstractAltitude} = h1._val == h2._val
Base.:(≈)(h1::H, h2::H) where {H<:AbstractAltitude} = h1._val ≈ h2._val

function Base.:-(::H1, ::H2) where {H1<:AbstractAltitude, H2<:AbstractAltitude}
    throw(ArgumentError("Can only subtract altitudes with the same reference"))
end

########################## Abstract3DPosition ##########################

abstract type Abstract3DPosition end

#avoid infinite recursion
Base.convert(::Type{P}, p::P) where {P<:Abstract3DPosition} = p

#### Geographic ####

#the default constructor generates a Geographic{NVector, EllipsoidalAlt} instance
Base.@kwdef struct Geographic{L <: Abstract2DLocation, H <: AbstractAltitude} <: Abstract3DPosition
    loc::L = NVector()
    alt::H = EllipsoidalAlt()
    function Geographic(loc::L, alt::H) where {L,H}
        alt._val >= h_min || throw(ArgumentError("Minimum altitude is $h_min m, got $h"))
        return new{L,H}(loc, alt)
    end
end
Geographic{L,H}(p::Abstract3DPosition) where {L,H} = convert(Geographic{L,H}, p)

function Base.convert(::Type{Geographic{L,H}}, p::Geographic) where {L,H}
    Geographic(convert(L, p.loc), H(p))
end

EllipsoidalAlt(p::Geographic{L,EllipsoidalAlt}) where {L} = p.alt
EllipsoidalAlt(p::Geographic{L,OrthometricAlt}) where {L} = EllipsoidalAlt(p.alt._val + get_geoid_offset(p.loc))

OrthometricAlt(p::Geographic{L,OrthometricAlt}) where {L} = p.alt
OrthometricAlt(p::Geographic{L,EllipsoidalAlt}) where {L} = EllipsoidalAlt(p.alt._val + get_geoid_offset(p.loc))

function Base.:(==)(p1::Abstract3DPosition, p2::Abstract3DPosition)
    error("Exact comparison between $(typeof(p1)) and $(typeof(p2)) not defined, use ≈ instead")
end
Base.:(≈)(p1::Abstract3DPosition, p2::Abstract3DPosition) = CartECEF(p1) ≈ CartECEF(p2)
Base.:(-)(p::T) where {T<:Abstract3DPosition} = convert(T, -CartECEF(p))

ltf(p::Abstract3DPosition, ψ_nl::Real = 0.0) = ltf(Geographic{NVector,EllipsoidalAlt}(p).loc, ψ_nl)
radii(p::Abstract3DPosition) = radii(Geographic{NVector,EllipsoidalAlt}(p).loc)
gravity(p::Abstract3DPosition) = gravity(Geographic{NVector,EllipsoidalAlt}(p))
g_l(p::Abstract3DPosition) = g_l(Geographic{NVector,EllipsoidalAlt}(p))
G_l(p::Abstract3DPosition) = G_l(Geographic{NVector,EllipsoidalAlt}(p))

############################# CartECEF #############################

struct CartECEF <: Abstract3DPosition
    data::SVector{3,Float64}
end
CartECEF(p::Abstract3DPosition) = convert(CartECEF, p)

Base.:(==)(r1::CartECEF, r2::CartECEF) = r1.data == r2.data
Base.:(≈)(r1::CartECEF, r2::CartECEF) = r1.data ≈ r2.data
Base.:(-)(r::CartECEF) = CartECEF(-r.data)

#### AbstractArray interface
Base.size(::CartECEF) = (3,)
Base.getindex(n::CartECEF, i) = getindex(n.data, i)

function Base.convert(::Type{Geographic{L,H}}, r::CartECEF) where {L,H}
    convert(Geographic{L,H}, convert(Geographic{NVector,EllipsoidalAlt}, r))
end

function Base.convert(::Type{Geographic{NVector,EllipsoidalAlt}}, r::CartECEF)

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

    n_e = NVector(SVector{3,Float64}(cos_ϕ*cos_λ, cos_ϕ*sin_λ, sin_ϕ))

    return Geographic(n_e, h)

end

function Base.convert(::Type{CartECEF}, p::Geographic)
    convert(CartECEF, convert(Geographic{NVector,EllipsoidalAlt}, p))
end

function Base.convert(::Type{CartECEF}, p::Geographic{NVector, EllipsoidalAlt})

    n_e = p.loc; h = p.alt
    _, N = radii(n_e)

    return CartECEF(SVector{3, Float64}(
        (N + h) * n_e[1],
        (N + h) * n_e[2],
        (N * (1 - e²) + h) * n_e[3]))

end

##### Generic Abstract3DPosition methods ####

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

    p_nvhe = Geographic{NVector,EllipsoidalAlt}(p)
    n_e = p_nvhe.loc
    h = p_nvhe.alt._val

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



########################## Plotting #################################

#unless a more specialized method is defined, a TimeHistory{<:Abstract3DPosition} is
#converted to Geographic{} for plotting
@recipe function plot_geodesypos(th::TimeHistory{<:AbstractVector{<:Abstract3DPosition}})

    v_geo = Vector{Geographic{LatLon,EllipsoidalAlt}}(undef, length(th.data))
    for i in 1:length(v_geo)
        v_geo[i] = Geographic{LatLon,EllipsoidalAlt}(th.data[i])
    end
    geo_sa = StructArray(v_geo) #this is now a struct of two arrays: loc and alt
    latlon_sa = StructArray(geo_sa.loc) #this is now a struct of two arrays: ϕ and λ
    alt_sa = StructArray(geo_sa.alt) #this now has a struct of one array: _val
    data = hcat(latlon_sa.ϕ/π, latlon_sa.λ/π, alt_sa._val)

    #maybe convert to degrees
    label --> ["Latitude" "Longitude" "Altitude"]
    yguide --> [L"$\varphi \ (\pi \ rad)$" L"$\lambda \ (\pi \ rad)$" L"$h \ (m)$"]
    th_split --> :h
    link --> :none #when th_split link defaults to :y, but we need a different scale for h


    return TimeHistory(th.t, data)

end



# end