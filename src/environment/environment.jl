module Environment

using Flight.Systems
using Flight.Atmosphere
using Flight.Terrain

import Flight.Systems: init, f_cont!, f_disc!

export Env

Base.@kwdef struct Env{A <: AbstractAtmosphere, T <: AbstractTerrain} <: SystemDescriptor
    atm::A = SimpleAtmosphere()
    trn::T = HorizontalTerrain()
end

function f_cont!(env::System{<:Env})
    f_cont!(env.atm)
    f_cont!(env.trn)
end

function f_disc!(env::System{<:Env})
    x_mod = false
    x_mod = x_mod || f_disc!(env.atm)
    x_mod = x_mod || f_disc!(env.trn)
end

end #module