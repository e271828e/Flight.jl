module Navigation

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightLib

using ...AircraftBase
using ...C172
using ..C172RPA


################################################################################
################################## Navigator ###################################

@kwdef struct Navigator <: SystemDefinition end

function Systems.f_disc!(::NoScheduling, nav::System{<:Navigator},
                        ::System{<:C172RPA.Vehicle})

end

function Systems.init!(nav::System{<:Navigator},
                            vehicle::System{<:C172RPA.Vehicle})

    @warn "Navigator init! not implemented"

end


end #module