# module TestSysNew

using Flight

import Flight.Modeling: f_cont!, f_disc!, init_x, init_y, init_u, init_d
import Flight.ModelingTraits: ContinuousStateTrait, OutputTrait
# using Flight.ModelingNew: ContinuousStateTrait, OutputTrait
export StatefulGroup, OutputfulGroup, MixedGroup


############## NoSteering ##############

struct Stateful <: SysDescNew end

ContinuousStateTrait(::Type{Stateful}) = HasContinuousStates()
init_x(::Type{Stateful}) = 0.0

f_cont!(::SysNew{Stateful}, args...) = nothing
f_disc!(::SysNew{Stateful}, args...) = false

struct Outputful <: SysDescNew end

# ContinuousStateTrait(::Type{Outputful}) = HasContinuousStates()
OutputTrait(::Type{Outputful}) = HasOutputs()
init_y(::Type{Outputful}) = 0.0

function f_cont!(sys::SysNew{Outputful}, args...)
    sys.y += 1
end

f_disc!(::SysNew{Outputful}, args...) = true

#######################

Base.@kwdef struct StatefulGroup <: SysGroupDescNew
    stateful1::Stateful = Stateful()
    stateful2::Stateful = Stateful()
end

Base.@kwdef struct OutputfulGroup <: SysGroupDescNew
    outputful1::Outputful = Outputful()
    outputful2::Outputful = Outputful()
end

Base.@kwdef struct MixedGroup <: SysGroupDescNew
    stateful::Stateful = Stateful()
    outputful::Outputful = Outputful()
end

#what i need to do when composing systems is:

# end #module