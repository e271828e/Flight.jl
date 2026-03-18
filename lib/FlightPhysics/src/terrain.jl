module Terrain

using StaticArrays
using CImGui: IsItemActive, SameLine, Text, Begin, End, AlignTextToFramePadding

using FlightCore
using ..Geodesy

export AbstractTerrain, NoTerrain, HorizontalTerrain
export TerrainData, SurfaceType
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

function GUI.draw!(mdl::Model{<:HorizontalTerrain},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Horizontal Terrain")

    u = mdl.u
    Begin(label, p_open)
        AlignTextToFramePadding(); Text("Surface Type"); SameLine()
        mode_button("Dry Tarmac", DryTarmac, DryTarmac, u[]; HSV_requested = HSV_gray)
        IsItemActive() && (u[] = DryTarmac); SameLine()
        mode_button("Wet Tarmac", WetTarmac, WetTarmac, u[]; HSV_requested = HSV_gray)
        IsItemActive() && (u[] = WetTarmac); SameLine()
        mode_button("Icy Tarmac", IcyTarmac, IcyTarmac, u[]; HSV_requested = HSV_gray)
        IsItemActive() && (u[] = IcyTarmac)
        CImGui.Text("Elevation (MSL): $(Float64(mdl.elevation)) m")
    End()
end

end #module
