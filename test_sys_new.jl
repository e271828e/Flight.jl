using Flight

import Flight.Modeling: f_cont!, f_disc!, init_x, init_y, init_u, init_d

############## NoSteering ##############

struct Stateful <: SysDescNew end

init_x(::Stateful) = [0.0]
init_y(::Stateful) = nothing #sytems are not required to have outputs

f_cont!(::SysNew{Stateful}, args...) = nothing
f_disc!(::SysNew{Stateful}, args...) = false

struct Outputful <: SysDescNew end

init_x(::Outputful) = nothing
init_y(::Outputful) = Ref(0.0)

f_cont!(::SysNew{Outputful}, args...) = nothing
f_disc!(::SysNew{Outputful}, args...) = false

Base.@kwdef struct TestGroup <: SysGroupDescNew
    stateful::Stateful = Stateful()
    outputful::Outputful = Outputful()
end

#what i need to do when composing systems is:
