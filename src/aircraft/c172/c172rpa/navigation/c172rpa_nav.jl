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


end #module