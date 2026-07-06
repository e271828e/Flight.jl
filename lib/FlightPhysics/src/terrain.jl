module Terrain

using LinearAlgebra, StaticArrays

using FlightCore
using FlightCore.GUI
using ..Geodesy
using ..Geodesy: AbstractGeometricAltitudeDatum

export AbstractTerrain, UniformTerrain
export TerrainData, SurfaceIntersection
export SurfaceType, DryTarmac, WetTarmac, IcyTarmac

@enum SurfaceType DryTarmac WetTarmac IcyTarmac

@kwdef struct TerrainData{T <: Abstract3DPosition}
    P::T = Geographic(NVector(), HOrth(0)) #position of surface point P
    kt_P_e::SVector{3,Float64} = -NVector()[:] #inward unit normal at P, ECEF coordinates
    surface::SurfaceType = DryTarmac
end

Geodesy.HOrth(data::TerrainData) = HOrth(data.P)

############################# AbstractTerrain ##################################

abstract type AbstractTerrain <: ModelDefinition end

function TerrainData(terrain::Model{<:AbstractTerrain}, location::Abstract2DLocation)
    throw(MethodError(TerrainData, (terrain, location)))
end

########################### SurfaceIntersection ################################

#result of intersecting a ray with the terrain surface. P denotes the
#intersection point. all other fields are meaningful only when valid == true,
#and default to values consistent with no intersection
@kwdef struct SurfaceIntersection{T <: Abstract3DPosition}
    valid::Bool = false #whether an intersection was found within bounds
    data::TerrainData{T} = TerrainData()
end

#intersect the ray with origin r_eO_e and unit direction u_e (both in ECEF
#coordinates) with the terrain surface. implementations must return the nearest
#intersection point P with ray parameter l ∈ [0, l_max], or an invalid
#SurfaceIntersection if there is none. kt_P_e must be unit-norm and
#inward-pointing, so that callers can adopt it directly as a frame axis
function SurfaceIntersection(terrain::Model{<:AbstractTerrain},
                            O::Abstract3DPosition,
                            u_e::AbstractVector{<:Real},
                            l_max::Real)
    throw(MethodError(SurfaceIntersection, (terrain, O, u_e, l_max)))
end

############################# UniformTerrain ################################

#terrain with uniform elevation and surface type. the elevation may be given as
#orthometric (constant altitude above the geoid) or ellipsoidal (constant
#altitude above the WGS84 ellipsoid; cheaper to query and analytically exact,
#since no geoid interpolation is involved)
@kwdef struct UniformTerrain{D <: AbstractGeometricAltitudeDatum} <: AbstractTerrain
    elevation::Altitude{D} = HOrth(0)
end

Modeling.U(::UniformTerrain) = Ref(DryTarmac)

@no_updates UniformTerrain

function TerrainData(terrain::Model{<:UniformTerrain}, loc::Abstract2DLocation)
    n_P_e = NVector(loc)
    TerrainData(; P = Geographic(n_P_e, terrain.elevation),
                kt_P_e = -n_P_e[:], surface = terrain.u[])
end

#the surface of constant orthometric altitude is approximated near the ray
#origin O by its local tangent plane, anchored at the terrain point P0 that
#shares O's 2D location. the plane's inward normal is the local (ellipsoidal)
#vertical at P0; the geoid slope is neglected, consistent with TerrainData.
#within this approximation the normal field is uniform over the ray span, so
#its value at P is that at P0. the approximation error grows with the
#intersection distance (~l/R_e rad in normal direction, ~l²/2R_e in surface
#height), so results are only meaningful for l ≪ R_e. should longer rays ever
#require it, the solution can be refined by re-anchoring the tangent plane at
#the candidate intersection's 2D location and intersecting again
function SurfaceIntersection(terrain::Model{<:UniformTerrain},
                            O::Abstract3DPosition,
                            u_e::AbstractVector{<:Real},
                            l_max::Real)

    r_eO_e = Geocentric(O)
    u_e = SVector{3,Float64}(u_e)

    #tangent plane anchor point P0 and inward surface normal
    n_P0_e = NVector(O)
    r_eP0_e = Geocentric(Geographic(n_P0_e, terrain.elevation))
    kt_P0_e = -n_P0_e[:] #inward surface normal is the reversed n-Vector

    #ray-plane intersection
    cos_α = kt_P0_e ⋅ u_e #cosine of angle between ray and terrain normal
    cos_α > 0 || return SurfaceIntersection() #parallel to or pointing away from the surface
    l = kt_P0_e ⋅ (r_eP0_e - r_eO_e) / cos_α
    0 <= l <= l_max || return SurfaceIntersection() #behind the origin or beyond l_max

    data = TerrainData(P = r_eO_e + l * u_e,  kt_P_e = kt_P0_e, surface = terrain.u[])
    SurfaceIntersection(; valid = true, data)

end

function GUI.draw!(mdl::Model{<:UniformTerrain},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Uniform Terrain")

    u = mdl.u
    BeginWindow(label, p_open)
        AlignTextToFramePadding(); TextFormatted("Surface Type"); SameLine()
        mode_button("Dry Tarmac", DryTarmac, DryTarmac, u[]; HSV_requested = HSV_gray)
        IsItemActive() && (u[] = DryTarmac); SameLine()
        mode_button("Wet Tarmac", WetTarmac, WetTarmac, u[]; HSV_requested = HSV_gray)
        IsItemActive() && (u[] = WetTarmac); SameLine()
        mode_button("Icy Tarmac", IcyTarmac, IcyTarmac, u[]; HSV_requested = HSV_gray)
        IsItemActive() && (u[] = IcyTarmac)
        datum = mdl.elevation isa HOrth ? "Orthometric" : "Ellipsoidal"
        TextFormatted("Elevation ($datum): $(Float64(mdl.elevation)) m")
    EndWindow()
end

end #module
