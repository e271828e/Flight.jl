module World

using UnPack

using Flight.FlightCore
using Flight.FlightPhysics
using Flight.FlightAircraft

export SimpleWorld

abstract type AbstractWorld <: Component end

Base.@kwdef struct SimpleWorld{A <: AircraftTemplate, E <: AbstractEnvironment} <: AbstractWorld
    aircraft::A = Cessna172R()
    environment::E = SimpleEnvironment()
end

function Systems.f_ode!(sys::System{<:SimpleWorld})
    @unpack aircraft, environment = sys
    f_ode!(environment)
    f_ode!(aircraft, environment)
    return nothing
end

#f_disc! and f_step! revert to their fallback methods

end #module