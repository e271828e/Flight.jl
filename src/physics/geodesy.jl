module Geodesy

using LinearAlgebra, StaticArrays, ComponentArrays, SHA, UnPack, Interpolations, HDF5

using Flight.FlightCore.Systems
using Flight.FlightCore.Plotting
using Flight.FlightCore.GUI

using ..Attitude

export Abstract2DLocation, NVector, LatLon
export Altitude, Ellipsoidal, Orthometric, Geopotential, HEllip, HOrth, HGeop
export Abstract3DPosition, Geographic, Cartesian
export ω_ie, gravity, g_n, G_n, ltf, radii, get_ψ_nl, get_geoid_height

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

Base.convert(::Type{L}, loc::L) where {L<:Abstract2DLocation} = loc


#### NVector ####

@kwdef struct NVector <: Abstract2DLocation
    data::SVector{3,Float64} = SVector{3}(1.0, 0.0, 0.0)
    function NVector(data::AbstractVector{T} where {T<:Real}; normalization::Bool = true)
        data = SVector{3,Float64}(data) #normalization will be faster for an SVector
        return normalization ? new(normalize(data)) : new(data)
    end
end

NVector(loc::Abstract2DLocation) = convert(NVector, loc)

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
Base.:(≈)(n1::NVector, n2::NVector; kwargs...) = ≈(n1.data, n2.data; kwargs...)
Base.:(-)(n::NVector) = NVector(-n.data)

#### AbstractArray interface
Base.size(::NVector) = (3,)
Base.length(::NVector) = 3
Base.getindex(n::NVector, i) = getindex(n.data, i)
#why does this allocate??
# Base.iterate(n::NVector, state = 1) = (state > 3 ? nothing : (n.data[state], state + 1))
LinearAlgebra.norm(n::NVector) = norm(getfield(n, :data)) #uses StaticArrays implementation
LinearAlgebra.normalize(n::NVector) = NVector(getfield(n, :data)) #let the constructor normalize

#### LatLon ####

@kwdef struct LatLon <: Abstract2DLocation
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
Base.:(≈)(ll1::LatLon, ll2::LatLon; kwargs...) = ≈(NVector(ll1), NVector(ll2); kwargs...)
Base.:(-)(latlon::LatLon) = LatLon(-NVector(latlon))

#time derivative of LatLon from NED frame transport rate (can be verified in
#[Groves])
function dt(ll::LatLon, ω_en_n::AbstractVector{<:Real})
    ϕ_dot = -ω_en_n[2]
    λ_dot = ω_en_n[1] / cos(ll.ϕ)
    # return SVector{2,Float64}(ϕ_dot, λ_dot)
    return ComponentVector(SVector(ϕ_dot, λ_dot), Axis(:ϕ, :λ))
end


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

#equivalent to egm_interp_from_hdf5 but with the additional hash check
function egm96_interp_from_bin()
    #the target file stores a 721x1441 Matrix{Float32} in low-endian binary
    #format. the matrix holds the data points for the EGM96 geoid height
    #measured from the WGS84 ellipsoid in 15 arc-minute resolution. latitude
    #goes from -π/2 to π/2, longitude from 0 to 2π
    file_path = joinpath(dirname(@__FILE__), "data", "ww15mgh_le.bin")
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
    ϕ_range = range(-π/2, π/2, size(data, 1))
    λ_range = range(0, 2π, size(data, 2))

    #return interpolator with extrapolation enabled to avoid machine precision
    #issues due to the boundaries being multiples of π
    linear_interpolation((ϕ_range, λ_range), data, extrapolation_bc = Line())
end

function egm96_interp_from_hdf5()

    data = Matrix{Float32}(undef, 721, 1441)
    file_path = joinpath(dirname(@__FILE__), "data", "ww15mgh_hdf5.h5")
    h5open(file_path) do file
        data .= file["geoid_height"] |> read
    end
    ϕ_range = range(-π/2, π/2, size(data, 1))
    λ_range = range(0, 2π, size(data, 2))
    linear_interpolation((ϕ_range, λ_range), data, extrapolation_bc = Line())
end

const egm96_interp = egm96_interp_from_hdf5() #geoid height interpolator

#need to pass geoid interpolator as an input, because since Julia 1.9 accessing
#the fields of a global const allocates, see:
#https://github.com/JuliaLang/julia/issues/49241
#https://github.com/JuliaLang/julia/issues/50317
function get_geoid_height(loc::Abstract2DLocation, geoid_height_interp = egm96_interp)
    #our longitude interval is [-π,π], but the table uses [0,2π], so we need to
    #correct for that
    latlon = LatLon(loc)
    ϕ = latlon.ϕ
    λ = mod(latlon.λ + 2π, 2π)
    geoid_height_interp(ϕ, λ)
end

@kwdef struct Altitude{D<:AbstractAltitudeDatum}
    _val::Float64 = 0.0
    #this constructor prevents the user from passing an Altitude with a
    #different datum, which would get implicitly converted to Float64 without
    #actually changing its value
    function Altitude{D}(h::Real) where {D}
        h >= h_min || throw(ArgumentError("Minimum altitude value is $h_min m, got $h"))
        return new{D}(h)
    end
end

const HEllip = Altitude{Ellipsoidal}
const HOrth = Altitude{Orthometric}
const HGeop = Altitude{Geopotential}

Altitude{D}(h::Altitude{D}, args...) where {D} = Altitude{D}(h._val)

#Ellipsoidal and Orthometric altitudes are related by the geoid's height at a
#given 2D location
Altitude{Ellipsoidal}(h::HOrth, loc::Abstract2DLocation) = HEllip(h._val + get_geoid_height(loc))
Altitude{Orthometric}(h::HEllip, loc::Abstract2DLocation) = HOrth(h._val - get_geoid_height(loc))

#Orthometric and Geopotential are directly related by the point-mass gravity
#approximation, 2D location not required
Altitude{Geopotential}(h::HOrth) = HGeop(h._val*a / (a+h._val))
Altitude{Orthometric}(h::HGeop) = HOrth(h._val*a / (a-h._val))

#still, for interface consistency we provide the two-argument method
Altitude{Geopotential}(h::HOrth, ::Abstract2DLocation) = HGeop(h)
Altitude{Orthometric}(h::HGeop, ::Abstract2DLocation) = HOrth(h)

#Geopotential and Ellipsoidal altitudes are related via Orthometric
Altitude{Geopotential}(h_ellip::HEllip, loc::Abstract2DLocation) = HGeop(HOrth(h_ellip, loc))
Altitude{Ellipsoidal}(h_geop::HGeop, loc::Abstract2DLocation) = HEllip(HOrth(h_geop), loc)

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
Base.:≈(h1::Altitude{D}, h2::Altitude{D}; kwargs...) where {D} = ≈(h1._val, h2._val; kwargs...)
Base.:>(h1::Altitude{D}, h2::Altitude{D}) where {D} = h1._val > h2._val
Base.:<(h1::Altitude{D}, h2::Altitude{D}) where {D} = h1._val < h2._val

Base.:(==)(h1::Altitude{D}, h2::Real) where {D} = ==(promote(h1,h2)...)
Base.:≈(h1::Altitude{D}, h2::Real; kwargs...) where {D} = ≈(promote(h1,h2)...; kwargs...)
Base.:>(h1::Altitude{D}, h2::Real) where {D} = >(promote(h1,h2)...)
Base.:<(h1::Altitude{D}, h2::Real) where {D} = <(promote(h1,h2)...)

Base.:(==)(h1::Real, h2::Altitude{D}) where {D} = h2 == h1
Base.:≈(h1::Real, h2::Altitude{D}; kwargs...) where {D} = ≈(h2, h1; kwargs...)
Base.:>(h1::Real, h2::Altitude{D}) where {D} = h2 < h1
Base.:<(h1::Real, h2::Altitude{D}) where {D} = h2 > h1


########################## Abstract3DPosition ##########################

abstract type Abstract3DPosition end

#avoid infinite recursion
Base.convert(::Type{P}, p::P) where {P<:Abstract3DPosition} = p

########################### Geographic ###############################

@kwdef struct Geographic{L <: Abstract2DLocation, H <: AbstractAltitudeDatum} <: Abstract3DPosition
    loc::L = NVector()
    h::Altitude{H} = HOrth()
end
Geographic(loc::Abstract2DLocation) = Geographic(loc, HOrth())
Geographic(pos::Abstract3DPosition) = Geographic{NVector,Ellipsoidal}(pos)
Geographic{L,H}(pos::Abstract3DPosition) where {L,H} = convert(Geographic{L,H}, pos)

function Base.convert(::Type{Geographic{L,H}}, geo::Geographic) where {L,H}
    Geographic(convert(L, geo.loc), Altitude{H}(geo))
end

NVector(geo::Geographic) = NVector(geo.loc)
LatLon(geo::Geographic) = LatLon(geo.loc)
Altitude{D}(geo::Geographic) where {D} = Altitude{D}(geo.h, geo.loc)

function Base.:(==)(geo1::Geographic{NVector,H}, geo2::Geographic{NVector,H}) where {H}
    return geo1.h == geo2.h && geo1.loc == geo2.loc
end

function Base.:(≈)(geo1::Geographic{L,H}, geo2::Geographic{L,H}; kwargs...) where {L,H}
    return ≈(geo1.h, geo2.h; kwargs...) && ≈(geo1.loc, geo2.loc; kwargs...)
end

function Base.:(≈)(loc1::Abstract3DPosition, loc2::Abstract3DPosition; kwargs...)
    ≈(Cartesian(loc1), Cartesian(loc2); kwargs...)
end

function Base.:(==)(loc1::Abstract3DPosition, loc2::Abstract3DPosition)
    throw(ArgumentError("Exact comparison between $(typeof(loc1)) and $(typeof(loc2)) not defined, use ≈ instead"))
end

Base.:(-)(pos::T) where {T<:Abstract3DPosition} = convert(T, -Cartesian(pos))


############################# Cartesian #############################

struct Cartesian <: Abstract3DPosition
    data::SVector{3,Float64}
end
Cartesian(pos::Abstract3DPosition) = convert(Cartesian, pos)
Cartesian() = Cartesian(Geographic())

NVector(r::Cartesian) = Geographic{NVector, Ellipsoidal}(r).loc
LatLon(r::Cartesian) = Geographic{LatLon, Ellipsoidal}(r).loc
Altitude{D}(r::Cartesian) where {D} = Altitude{D}(Geographic{NVector,D}(r))

Base.:(==)(r1::Cartesian, r2::Cartesian) = r1.data == r2.data
Base.:(≈)(r1::Cartesian, r2::Cartesian; kwargs...) = ≈(r1.data, r2.data; kwargs...)
Base.:(-)(r::Cartesian) = Cartesian(-r.data)
Base.:(+)(r1::Cartesian, r2::AbstractVector{<:Real}) = Cartesian(r1.data + SVector{3,Float64}(r2))
Base.:(+)(r1::AbstractVector{<:Real}, r2::Cartesian) = r2 + r1

#### AbstractArray interface
Base.size(::Cartesian) = (3,)
Base.getindex(n::Cartesian, i) = getindex(n.data, i)

function Base.convert(::Type{Geographic{L,H}}, r::Cartesian) where {L,H}
    convert(Geographic{L,H}, convert(Geographic{NVector,Ellipsoidal}, r))
end

function Base.convert(::Type{Geographic{NVector,Ellipsoidal}}, r::Cartesian)

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

function Base.convert(::Type{Cartesian}, geo::Geographic)
    convert(Cartesian, convert(Geographic{NVector,Ellipsoidal}, geo))
end

function Base.convert(::Type{Cartesian}, geo::Geographic{NVector, Ellipsoidal})

    n_e = geo.loc; h = Float64(geo.h)
    _, N = radii(n_e)

    return Cartesian(SVector{3, Float64}(
        (N + h) * n_e[1],
        (N + h) * n_e[2],
        (N * (1 - e²) + h) * n_e[3]))

end

Base.convert(::Type{SVector{3,Float64}}, r::Cartesian) = r.data

##### Generic Abstract3DPosition methods ####

#general conversion from 3D to 2D location
(::Type{L})(pos::Abstract3DPosition) where {L<:Abstract2DLocation} = Geographic{L,Ellipsoidal}(pos).loc

ltf(pos::Abstract3DPosition, ψ_nl::Real = 0.0) = ltf(Geographic{NVector,Ellipsoidal}(pos).loc, ψ_nl)
radii(pos::Abstract3DPosition) = radii(Geographic{NVector,Ellipsoidal}(pos).loc)

"""
    gravity(p::Abstract3DPosition)

Compute normal gravity.

Computation is based on Somigliana's formula for gravity at the ellipsoid
surface, with a second order altitude correction, accurate for small altitudes
above the WGS84 ellipsoid (h<<a). Direction is assumed normal to the WGS84
ellipsoid, a good enough approximation for most navigation applications. See
Hoffmann & Moritz.
"""
function gravity(pos::Abstract3DPosition)

    p_nve = Geographic{NVector,Ellipsoidal}(pos)
    n_e = p_nve.loc
    h = Float64(p_nve.h)

    sin²ϕ = n_e[3]^2
    cos²ϕ = n_e[1]^2 + n_e[2]^2

    #gravity at the ellipsoid surface (Somigliana)
    γ_0 = (a * γ_a * cos²ϕ + b * γ_b * sin²ϕ) / √(a² * cos²ϕ + b² * sin²ϕ) #[Hof06] 2-146

    #altitude correction
    γ = γ_0 * (1 - 2/a * (1 + f + m - 2f * sin²ϕ) * h + 3/a² * h^2)

    return γ

end

"""
    g_n(p::Abstract3DPosition)

Compute gravity vector resolved in the NED frame.
"""
g_n(pos::Abstract3DPosition) = SVector{3}(0, 0, gravity(pos))

"""
    G_n(p::Abstract3DPosition)

Compute gravitational attraction resolved in the NED frame.
"""
function G_n(pos::Abstract3DPosition)

    q_en = ltf(pos)
    ω_ie_e = SVector{3, Float64}(0,0,ω_ie)
    r_eP_e = Cartesian(pos)[:]
    G_n = g_n(pos) + q_en'(ω_ie_e × (ω_ie_e × r_eP_e))
    return G_n

end

############################### Plotting #######################################

#if no specific method available, convert to LatLon for plotting
@recipe function fp(th::TimeHistory{<:Abstract2DLocation})

    return TimeHistory(th._t, [LatLon(v) for v in th._data])

end

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