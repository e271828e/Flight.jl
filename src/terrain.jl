module Terrain

using StaticArrays

using Flight.Geodesy

export AbstractTerrain, DummyTerrain, HorizontalTerrain
export TerrainData, SurfaceType
export get_terrain_data

@enum SurfaceType DryTarmac WetTarmac IcyTarmac

struct TerrainData
    altitude::Altitude{Orthometric}
    normal::SVector{3,Float64} #NED components, inward pointing
    surface::SurfaceType
end

TerrainData(; altitude = AltOrth(0),
                normal = SVector{3,Float64}(0,0,1),
                surface = DryTarmac) =
    TerrainData(altitude, SVector{3,Float64}(normal), surface)


######################## AbstractTerrain ##########################

abstract type AbstractTerrain end
get_terrain_data(::AbstractTerrain, args...) = throw(MethodError(get_terrain_data, args))

struct DummyTerrain <: AbstractTerrain end

struct HorizontalTerrain <: AbstractTerrain
    altitude::Altitude{Orthometric}
    surface::SurfaceType
end

HorizontalTerrain(; altitude = AltOrth(0), surface = DryTarmac) =
    HorizontalTerrain(altitude, surface)

get_terrain_data(trn::HorizontalTerrain, ::Abstract2DLocation) =
    TerrainData(trn.altitude, SVector{3,Float64}(0,0,1), trn.surface)

end #module

#TerrainModel does not belong to the Aircraft itself. it must be defined
#separately, and passed as an external data source. the same goes for
#AtmosphereModel. there must be a level above the Aircraft, which will typically
#be the simulation scheduler, that defines both the environmental models and all
#the aircraft participating in the simulation. this may be a block based
#simulation engine or a custom made one.

#the AtmosphericModel contained in the Aircraft should behave as a gateway to
#the actual AtmosphericModel, through which the evolving atmospheric model can
#be queried by the Aircraft for the values at its current location. the
#TerrainModel can be handled similarly, because even if it is not evolving in
#time, it may be an arbitrarily complex terrain database which must serve all
#vehicles. initially, these Model "clients" will be simply references to the
#models themselves, because these will be constant, simple and shared with no
#one else

#an alternative approach, instead of querying the models in a client-server
#architecture, would be simply to provide access references to the environment
#model data and methods and store them in the external data source fields of the
#Aircraft System being simulated. then, the Aircraft can locally evaluate those
#methods for their particular position. in that case, the exchanged information
#is not the product of model evaluation at the Aircraft location, but the data
#required to perform those evaluations within the Aircraft model.

#both scenarios could be realized with an Observer model, with any vehicle in
#the simulation subscribing to the AtmosphericModel and the TerrainModel. each
#time these update, they notify all their subscribers so they can also be
#evolved by the scheduler. if this communication occurs through Channels, each
#Model (Atmospheric, Terrain, Vehicle...) can run on separate Tasks, and each of
#these in a different thread. Thread safety.
