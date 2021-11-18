module Mass

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.ModelingTools

using Flight.Airframe

using Flight.Dynamics

import Flight.ModelingTools: System, get_x0, get_y0, get_u0, get_d0,f_cont!, f_disc!
import Flight.Dynamics: MassData

using Flight.Plotting
import Flight.Plotting: plots

export AbstractMass, ConstantMass, TunableMass


############################ AbstractMass ######################################

abstract type AbstractMass <: AbstractComponent end

#given the fuel system outputs and the installed payloads, a MassSystem must
#return a MassData instance. open questions:

#1) should this information be passed as inputs along with the MassSystem
#   instance to the MassData constructor? or should they be defined as part of
#   the MassSystem's inputs and assigned explicitly?
function MassData(::T) where {T<:System{<:AbstractMass}}
    error("MassData constructor not implemented for $T")
end

############################ ConstantMass ######################################

Base.@kwdef struct ConstantMass <: AbstractMass
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end

f_cont!(::System{ConstantMass}) = nothing
f_disc!(::System{ConstantMass}) = false

function MassData(sys::System{ConstantMass})
    @unpack m, J_Ob_b, r_ObG_b = sys.params
    MassData(m, J_Ob_b, r_ObG_b)
end

############################ TunableMass ######################################

struct TunableMass <: AbstractMass end

Base.@kwdef mutable struct TunableMassU
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end

get_u0(::TunableMass) = TunableMassU()

f_cont!(::System{TunableMass}) = nothing
f_disc!(::System{TunableMass}) = false

function MassData(sys::System{TunableMass})
    @unpack m, J_Ob_b, r_ObG_b = sys.u
    MassData(m, J_Ob_b, r_ObG_b)
end


end #module