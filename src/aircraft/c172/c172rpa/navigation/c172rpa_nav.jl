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

function Systems.f_disc!(::NoScheduling, nav::System{<:Navigator},
                        ::System{<:C172RPA.Vehicle})

end

function AircraftBase.trim!(nav::System{<:Navigator},
                            vehicle::System{<:C172RPA.Vehicle})

    @warn "Navigator trim not implemented"

end


end #module