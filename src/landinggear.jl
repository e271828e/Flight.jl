module LandingGear
#rename to ground, rename LandingGearLeg to LandingGear

using StaticArrays
using ComponentArrays

using Flight.ModelingTools
using Flight.Airdata
using Flight.Dynamics
using Flight.Airframe

#for method extension without module qualifier
import Flight.ModelingTools: get_x0, get_d0, get_y0, get_u0, f_cont!, f_disc!
import Flight.Airframe: get_wr_b, get_hr_b

export LandingGearLeg

"""
SET NORMALIZATION TO false!!!!!!!
"""

#first implement a LandingGearLeg with a basic shock absorber, direct steering
#and direct braking. then add the required generalizations to allow for
#NoSteering, NoBraking, etc.

#aqui faltan muchos parametros: los μ no deben ser Floats, sino permitir static,
#dynamic, falta breakout velocity interval
Base.@kwdef struct ContactModel
    μ_roll::Float64
    μ_skid::Float64
    k_p::Float64
    k_i::Float64
end
get_x0(::ContactModel) = ComponentVector(x = 0.0, y = 0.0) #v regulator integrator states
get_y0(::ContactModel) = ComponentVector(v = zeros(2), s = zeros(2), α_p = zeros(2), α_i = zeros(2),
                                    α_raw = zeros(2), α_sat = zeros(2), α = zeros(2))

#a shock absorber does not have to be a System, so strictly speaking we wouldn't
#need to derive from AbstractComponent
abstract type AbstractShockAbsorber <: AbstractComponent end
get_force(::A) where {A<:AbstractShockAbsorber} = no_extend_error(get_force, A)

Base.@kwdef struct StandardShockAbsorber <: AbstractShockAbsorber
    l_0::Float64 = 1 #natural length
    k_s::Float64 = 25000 #spring constant
    k_d_ext::Float64 = 1000 #extension damping coefficient
    k_d_cmp::Float64 = 1000 #compression damping coefficient
    ξ_min::Float64 = -1 #compression below which the shock absorber is disabled
end

abstract type AbstractSteering <: AbstractComponent end

struct NoSteering <: AbstractSteering end

struct DirectSteering <: AbstractSteering
    limits::SVector{2,Float64}
end
get_u0(::DirectSteering) = 0.0 #steering angle input
get_y0(::DirectSteering) = 0.0 #steering angle

# struct ActuatedSteering{A <: AbstractActuator} <: AbstractSteering
#     limits::SVector{2,Float64} #more generally, transmission kinematics could go here
#     actuator::A #actuator model parameters go here
# end
# get_x0(steering::ActuatedSteering) = get_x0(steering.actuator)
# get_y0(steering::ActuatedSteering) = get_x0(steering.actuator)
# get_u0(steering::ActuatedSteering) = get_u0(steering.actuator) #typically, 1

abstract type AbstractBraking <: AbstractComponent end

struct NoBraking <: AbstractBraking end

struct DirectBraking <: AbstractBraking
    efficiency::Float64
end
get_u0(::DirectBraking) = 0.0 #braking strength input
get_y0(::DirectBraking) = 0.0 #braking strength

Base.@kwdef struct LandingGearLeg{A<:AbstractShockAbsorber, S <: AbstractSteering, B <: AbstractBraking} <: AbstractAirframeComponent
    frame::FrameSpec
    contact::ContactModel #vRegulator goes here
    shock::A
    steering::S
    braking::B
end

function get_x0(ldg::LandingGearLeg{A,S,B}) where {A,S,B}
    x_blocks = Dict(:Symbol, Any)
    push!(x_blocks, :contact => get_x0(ldg.contact))
    push!(x_blocks, :shock => get_x0(ldg.shock))
    push!(x_blocks, :steering => get_x0(ldg.steering))
    push!(x_blocks, :braking => get_x0(ldg.braking))
    #remove null blocks, then convert to named tuple, then pass to
    #ComponentVector and instantiate
end



# get_steering_angle(x, u, ldg::LandingGearLeg{})
#if the steering is actuated, it will have its own state vector, and
#get_steering_angle will require it

#get_steering_angle(x, u, ldg::Ldg{A, NoSteering, C})
#get_steering_angle(x, u, ldg::Ldg{A, DirectSteering, C})
#get_steering_angle(x, u, ldg::Ldg{A, ActuatedSteering, C})

#methods that need to dispatch on the braking type parameter:
#get_braking_strength(x, u, ldg::Ldg{A,S,NoBraking} where {A,S} ) = 0
#get_braking_strength(x, u, ldg::Ldg{A,S,DirectBraking} where {A,S} ) = u.braking
#get_μ_max(x,u,ldg) calls braking_strength (from zero to one), but in the future
#could also #account for other braking model parameters


#################### AbstractSystem interface


end


#in System, define and extend f_branch!

# #individual Component
# f_branch!(y, dx, x, u, t, sys, args...) = f_branch!(Val(has_input(sys)), y, dx, x, u, t, args...)
# f_branch!(::Val{true}, y, dx, x, u, t, sys, args...) = f_cont!(y, dx, x, u, t, sys, args...)
# f_cont!(::HasInput, y, dx, x ,u, t, sys, args...) = f_cont!(y, dx, x, u, t, sys, args...)
# f_cont!(::HasNoInput, y, dx, x, u, t, sys, args...) = f_cont!(y, dx, x, t, sys, args...)

# #for a AirframeGroup
# f_cont!(MaybeInput(S), MaybeOutput(S), y, dx, x, u, t, sys, args...)
# f_cont!(::HasInput, ::HasOutput, y, dx, x ,u, t, sys, args...)
# #now, this method needs to consider the possibility for each component that it
# #may have or not Input or Output. so it must do
# for (label, component) in zip(keys(C), values(C))
#     if MaybeInput(typeof(component)) #need tocheck, because if it has no input, u[label] will not exist!
#         f_cont!(y_cmp, dx_cmp, x_cmp, u_cmp, t, cmp, args...)
#     else
#         f_cont!(y_cmp, dx_cmp, x_cmp, t, cmp, args...)
#     end
# end


#this method