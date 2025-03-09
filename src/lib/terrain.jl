module Terrain

using StaticArrays

using Flight.FlightCore
using ..Geodesy

export AbstractTerrain, NoTerrain, HorizontalTerrain
export TerrainData, SurfaceType
export SurfaceType, DryTarmac, WetTarmac, IcyTarmac

@enum SurfaceType DryTarmac WetTarmac IcyTarmac

@kwdef struct TerrainData
    altitude::Altitude{Orthometric} = HOrth(0)
    normal::SVector{3,Float64} = @SVector[0.0, 0.0, 1.0] #NED components, inward pointing
    surface::SurfaceType = DryTarmac
end

Geodesy.HOrth(data::TerrainData) = data.altitude

######################## AbstractTerrain ##########################

abstract type AbstractTerrain <: SystemDefinition end

function TerrainData(trn::System{<:AbstractTerrain}, loc::Abstract2DLocation)
    throw(MethodError(TerrainData, (trn, loc)))
end

# struct NoTerrain <: AbstractTerrain end

# function TerrainData(::NoTerrain, ::Abstract2DLocation = NVector())
#     TerrainData(Geodesy.h_min + 1, SVector{3,Float64}(0,0,1), DryTarmac)
# end

@kwdef struct HorizontalTerrain <: AbstractTerrain
    altitude::Altitude{Orthometric} = HOrth(0)
end

Systems.U(::HorizontalTerrain) = Ref(DryTarmac)

function TerrainData(trn::System{<:HorizontalTerrain}, ::Abstract2DLocation)
    TerrainData(trn.constants.altitude, SVector{3,Float64}(0,0,1), trn.u[])
end

end #module
