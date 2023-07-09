module Environment

using Reexport

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI

@reexport using ..Atmosphere
@reexport using ..Terrain

export AbstractEnvironment, SimpleEnvironment

abstract type AbstractEnvironment <: Component end

get_atmosphere(::System{<:AbstractEnvironment}) = throw(MethodError(get_atmosphere, sys))
get_terrain(::System{<:AbstractEnvironment}) = throw(MethodError(get_terrain, sys))


################################################################################
############################# SimpleEnvironment #################################

Base.@kwdef struct SimpleEnvironment{A <: AbstractAtmosphere, T <: AbstractTerrain} <: AbstractEnvironment
    atm::A = SimpleAtmosphere()
    trn::T = HorizontalTerrain()
end

get_atmosphere(env::System{<:SimpleEnvironment}) = env.atm
get_terrain(env::System{<:SimpleEnvironment}) = env.trn

################################# GUI ##########################################

function GUI.draw!(sys::System{<:SimpleEnvironment}, label::String = "Environment")

    CImGui.Begin(label)

    show_atm = @cstatic check=false @c CImGui.Checkbox("Atmosphere", &check)
    show_trn = @cstatic check=false @c CImGui.Checkbox("Terrain", &check)

    show_atm && GUI.draw!(sys.atm, "Atmosphere")
    show_trn && GUI.draw(sys.trn, "Terrain")

    CImGui.End()

end

end #module