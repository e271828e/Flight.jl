module Environment

using Reexport

using Flight.Systems
@reexport using Flight.Atmosphere
@reexport using Flight.Terrain

import Flight.Systems: init, f_cont!, f_disc!

export AbstractEnvironment, SimpleEnvironment

abstract type AbstractEnvironment <: SystemDescriptor end

Base.@kwdef struct SimpleEnvironment{A <: AbstractAtmosphere, T <: AbstractTerrain} <: AbstractEnvironment
    atm::A = SimpleAtmosphere()
    trn::T = HorizontalTerrain()
end

function f_cont!(env::System{<:SimpleEnvironment})
    f_cont!(env.atm)
    f_cont!(env.trn)
end

function f_disc!(env::System{<:SimpleEnvironment})
    x_mod = false
    x_mod = x_mod || f_disc!(env.atm)
    x_mod = x_mod || f_disc!(env.trn)
end

end #module