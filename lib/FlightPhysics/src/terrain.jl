module Terrain

using LinearAlgebra, StaticArrays

using FlightCore
using FlightCore.GUI
using ..Geodesy

export AbstractTerrain, HorizontalTerrain
export TerrainData, SurfaceIntersection
export SurfaceType, DryTarmac, WetTarmac, IcyTarmac

@enum SurfaceType DryTarmac WetTarmac IcyTarmac

@kwdef struct TerrainData
    elevation::Altitude{Orthometric} = HOrth(0)
    normal::SVector{3,Float64} = @SVector[0.0, 0.0, 1.0] #NED components, inward pointing
    surface::SurfaceType = DryTarmac
end

Geodesy.HOrth(data::TerrainData) = data.elevation

############################# AbstractTerrain ##################################

abstract type AbstractTerrain <: ModelDefinition end

function TerrainData(terrain::Model{<:AbstractTerrain}, location::Abstract2DLocation)
    throw(MethodError(TerrainData, (terrain, location)))
end

########################### SurfaceIntersection ################################

#result of intersecting a ray with the terrain surface. P denotes the
#intersection point. all other fields are meaningful only when valid == true,
#and default to values consistent with no intersection
@kwdef struct SurfaceIntersection
    valid::Bool = false #whether an intersection was found within bounds
    r_eP_e::SVector{3,Float64} = zeros(SVector{3}) #position of P, ECEF coordinates
    kt_P_e::SVector{3,Float64} = zeros(SVector{3}) #inward unit surface normal at P, ECEF coordinates
    surface::SurfaceType = DryTarmac #surface type at P
end

#intersect the ray with origin r_eO_e and unit direction u_e (both in ECEF
#coordinates) with the terrain surface. implementations must return the nearest
#intersection point P with ray parameter l ∈ [0, l_max], or an invalid
#SurfaceIntersection if there is none. kt_P_e must be unit-norm and
#inward-pointing, so that callers can adopt it directly as a frame axis
function SurfaceIntersection(terrain::Model{<:AbstractTerrain},
                            r_eO_e::AbstractVector{<:Real},
                            u_e::AbstractVector{<:Real},
                            l_max::Real)
    throw(MethodError(SurfaceIntersection, (terrain, r_eO_e, u_e, l_max)))
end

############################# HorizontalTerrain ################################

#flat terrain with constant orthometric elevation
@kwdef struct HorizontalTerrain <: AbstractTerrain
    elevation::Altitude{Orthometric} = HOrth(0)
end

Modeling.U(::HorizontalTerrain) = Ref(DryTarmac)

@no_updates HorizontalTerrain

function TerrainData(terrain::Model{<:HorizontalTerrain}, ::Abstract2DLocation)
    TerrainData(terrain)
end

function TerrainData(terrain::Model{<:HorizontalTerrain})
    TerrainData(terrain.elevation, SVector{3,Float64}(0,0,1), terrain.u[])
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
function SurfaceIntersection(trn::Model{<:HorizontalTerrain},
                            r_eO_e::AbstractVector{<:Real},
                            u_e::AbstractVector{<:Real},
                            l_max::Real)

    r_eO_e = Cartesian(r_eO_e)
    u_e = SVector{3,Float64}(u_e)

    #tangent plane anchor point P0 and inward surface normal
    n_P0_e = NVector(r_eO_e) #n-vector at O (and P0)
    r_eP0_e = Cartesian(Geographic(n_P0_e, trn.elevation))
    kt_P0_e = -SVector{3,Float64}(n_P0_e...) #inward surface normal is the reversed n-Vector

    #ray-plane intersection
    cos_α = kt_P0_e ⋅ u_e #cosine of angle between ray and terrain normal
    cos_α > 0 || return SurfaceIntersection() #parallel to or pointing away from the surface
    l = kt_P0_e ⋅ (r_eP0_e - r_eO_e) / cos_α
    0 <= l <= l_max || return SurfaceIntersection() #behind the origin or beyond l_max

    SurfaceIntersection(; valid = true, r_eP_e = r_eO_e + l * u_e,
                        kt_P_e = kt_P0_e, surface = trn.u[])

end

function GUI.draw!(mdl::Model{<:HorizontalTerrain},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Horizontal Terrain")

    u = mdl.u
    BeginWindow(label, p_open)
        AlignTextToFramePadding(); TextFormatted("Surface Type"); SameLine()
        mode_button("Dry Tarmac", DryTarmac, DryTarmac, u[]; HSV_requested = HSV_gray)
        IsItemActive() && (u[] = DryTarmac); SameLine()
        mode_button("Wet Tarmac", WetTarmac, WetTarmac, u[]; HSV_requested = HSV_gray)
        IsItemActive() && (u[] = WetTarmac); SameLine()
        mode_button("Icy Tarmac", IcyTarmac, IcyTarmac, u[]; HSV_requested = HSV_gray)
        IsItemActive() && (u[] = IcyTarmac)
        TextFormatted("Elevation (MSL): $(Float64(mdl.elevation)) m")
    EndWindow()
end

end #module
