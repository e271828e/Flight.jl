module Environment

using Reexport

using Flight.Engine.Systems
@reexport using ..Atmosphere
@reexport using ..Terrain

export AbstractEnvironment, SimpleEnvironment

abstract type AbstractEnvironment <: Component end

Base.@kwdef struct SimpleEnvironment{A <: AbstractAtmosphere, T <: AbstractTerrain} <: AbstractEnvironment
    atm::A = SimpleAtmosphere()
    trn::T = HorizontalTerrain()
end

#System methods auto-generated

end #module