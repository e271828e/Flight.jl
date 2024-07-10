module Navigation

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightPhysics
using Flight.FlightComponents

using ...AircraftBase
using ...C172
using ..C172RPA


################################################################################
################################## Navigator ###################################

@kwdef struct Navigator <: SystemDefinition end

@kwdef mutable struct NavigatorS <: SystemDefinition
    n_disc::Int = 0
end

@kwdef struct NavigatorY <: SystemDefinition
    n_disc::Int = 0
    Δt::Float64 = 0
end

Systems.S(::Navigator) = NavigatorS()
Systems.Y(::Navigator) = NavigatorY()

function Systems.reset!(sys::System{<:Navigator})
    sys.s.n_disc = 0
    foreach(ss -> Systems.reset!(ss), sys.subsystems) #reset subcontrollers
end

function Systems.f_disc!(nav::System{<:Navigator}, Δt::Real)
    nav.s.n_disc += 1
    nav.y = NavigatorY(; n_disc = nav.s.n_disc, Δt)
end


end #module