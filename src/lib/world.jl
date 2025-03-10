module World

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightLib.AircraftBase: Aircraft

export AbstractWorld, SimpleWorld

################################################################################
################################### World ######################################

abstract type AbstractWorld <: SystemDefinition end

################################################################################
################################### World ######################################

@kwdef struct SimpleWorld{C <: Aircraft, A <: AbstractAtmosphere, T <: AbstractTerrain} <: AbstractWorld
    ac::C = Aircraft()
    atm::A = SimpleAtmosphere()
    trn::T = HorizontalTerrain()
end

function Systems.f_ode!(world::System{<:SimpleWorld})
    @unpack ac, atm, trn = world.subsystems
    f_ode!(atm)
    f_ode!(trn)
    f_ode!(ac, atm, trn)
    update_y!(world)
end

function Systems.f_disc!(::NoScheduling, world::System{<:SimpleWorld})
    @unpack ac, atm, trn = world.subsystems
    f_disc!(atm)
    f_disc!(trn)
    f_disc!(ac, atm, trn)
    update_y!(world)
end

function Systems.f_step!(world::System{<:SimpleWorld})
    @unpack ac, atm, trn = world.subsystems
    f_step!(atm)
    f_step!(trn)
    f_step!(ac, atm, trn)
end

function Systems.init!( world::System{<:SimpleWorld}, args...)
    @unpack ac, atm, trn = world.subsystems
    Systems.init!(atm)
    Systems.init!(trn)
    Systems.init!(ac, args...)
end

end
