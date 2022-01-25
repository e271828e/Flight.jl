# module TestSysNew

using Flight
using StaticArrays

import Flight.Modeling: f_cont!, f_disc!, init_x, init_y, init_u, init_d
export StatefulGroup, OutputfulGroup, MixedGroup


############## NoSteering ##############

struct Stateful <: SysDescNew end

init_x(::Type{Stateful}) = [0.0]
init_u(::Type{Stateful}) = Ref(true)

@inline function f_cont!(sys::SysNew{Stateful}, args...)
    sys.x .+= 1.5
end
@inline f_disc!(::SysNew{Stateful}, args...) = false

struct Outputful <: SysDescNew end

init_y(::Type{Outputful}) = 0.0
init_d(::Type{Outputful}) = [1,2,3]

@inline function f_cont!(sys::SysNew{Outputful}, args...)
    sys.y += 1
end

@inline f_disc!(::SysNew{Outputful}, args...) = true

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