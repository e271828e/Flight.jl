using Flight
using StaticArrays

import Flight.Modeling: f_cont!, f_disc!, init_x, init_y, init_u, init_d
export StatefulGroup, OutputfulGroup, MixedGroup


############## NoSteering ##############

struct Stateful <: SystemDescriptor end

init_x(::Type{Stateful}) = [0.0]
init_u(::Type{Stateful}) = Ref(true)

@inline function f_cont!(sys::System{Stateful}, args...)
    sys.x .+= 1.5
end
@inline f_disc!(::System{Stateful}, args...) = false

struct Outputful <: SystemDescriptor end

init_y(::Type{Outputful}) = 0.0
init_d(::Type{Outputful}) = [1,2,3]

@inline function f_cont!(sys::System{Outputful}, args...)
    sys.y += 1
end

@inline f_disc!(::System{Outputful}, args...) = true

#######################

Base.@kwdef struct StatefulGroup <: SystemGroupDescriptor
    stateful1::Stateful = Stateful()
    stateful2::Stateful = Stateful()
end

Base.@kwdef struct OutputfulGroup <: SystemGroupDescriptor
    outputful1::Outputful = Outputful()
    outputful2::Outputful = Outputful()
end

Base.@kwdef struct MixedGroup <: SystemGroupDescriptor
    stateful::Stateful = Stateful()
    outputful::Outputful = Outputful()
end
