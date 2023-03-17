module Terrain

using StaticArrays

using Flight.FlightCore.Systems

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

abstract type AbstractTerrain <: Component end
TerrainData(::System{<:AbstractTerrain}, args...) = throw(MethodError(TerrainData, args))

struct DummyTerrain <: AbstractTerrain end

struct HorizontalTerrain <: AbstractTerrain
    altitude::Altitude{Orthometric}
    surface::SurfaceType
end

HorizontalTerrain(; altitude = HOrth(0), surface = DryTarmac) =
    HorizontalTerrain(altitude, surface)

TerrainData(trn::System{<:HorizontalTerrain}, loc::Abstract2DLocation) =
    TerrainData(loc, trn.params.altitude, SVector{3,Float64}(0,0,1), trn.params.surface)

end #module
