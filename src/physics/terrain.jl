module Terrain

using StaticArrays

using ..Geodesy

export AbstractTerrain, DummyTerrain, HorizontalTerrain
export TerrainData, SurfaceType
export SurfaceType, DryTarmac, WetTarmac, IcyTarmac

@enum SurfaceType DryTarmac WetTarmac IcyTarmac

struct TerrainData
    location::NVector
    altitude::Altitude{Orthometric}
    normal::SVector{3,Float64} #NED components, inward pointing
    surface::SurfaceType
end

function TerrainData(; location = NVector(),
                       altitude = HOrth(0),
                       normal = SVector{3,Float64}(0,0,1),
                       surface = DryTarmac)
    TerrainData(location, altitude, SVector{3,Float64}(normal), surface)
end

function Geodesy.Altitude{D}(data::TerrainData) where {D}
    Altitude{D}(data.altitude, data.location)
end

######################## AbstractTerrain ##########################

abstract type AbstractTerrain end
TerrainData(::AbstractTerrain, loc::Abstract2DLocation) = throw(MethodError(TerrainData, loc))

struct DummyTerrain <: AbstractTerrain end

@kwdef struct HorizontalTerrain <: AbstractTerrain
    altitude::Altitude{Orthometric} = HOrth(0)
    surface::SurfaceType = DryTarmac
end

function TerrainData(trn::HorizontalTerrain, loc::Abstract2DLocation = NVector())
    TerrainData(loc, trn.altitude, SVector{3,Float64}(0,0,1), trn.surface)
end

end #module
