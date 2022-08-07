module Environment

using Reexport

using Flight.Systems
@reexport using Flight.Atmosphere
@reexport using Flight.Terrain

import Flight.Systems: init, f_ode!, f_step!

export AbstractEnvironment, SimpleEnvironment

abstract type AbstractEnvironment <: SystemDescriptor end

Base.@kwdef struct SimpleEnvironment{A <: AbstractAtmosphere, T <: AbstractTerrain} <: AbstractEnvironment
    atm::A = SimpleAtmosphere()
    trn::T = HorizontalTerrain()
end

function f_ode!(env::System{<:SimpleEnvironment})
    f_ode!(env.atm)
    f_ode!(env.trn)
end

function f_step!(env::System{<:SimpleEnvironment})
    x_mod = false
    x_mod = x_mod || f_step!(env.atm)
    x_mod = x_mod || f_step!(env.trn)
end

end #module