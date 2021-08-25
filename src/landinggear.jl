module LandingGear
#rename to ground, rename LandingGearLeg to LandingGear

using StaticArrays
using ComponentArrays

using Flight.Airdata
using Flight.Dynamics
using Flight.Component
import Flight.Component: get_wr_Ob_b, get_h_Gc_b
import Flight.System: X, D, Y, U, f_cont!, f_disc!

export LandingGearLeg, LandingGearGroup

abstract type AbstractLandingGearLeg <: AbstractComponent end

Base.@kwdef struct LandingGearLeg <: AbstractLandingGearLeg
end

#=

#START BY REPLICATING PYTHON IMPLEMENTATION BUT WITHIN THIS FRAMEWORK
abstract type AbstractActuator end

Base.@kwdef struct LandingGearLeg{A<:AbstractShockAbsorber, S <: AbstractSteering, B <: AbstractBraking}
    frame::FrameSpec
    contact::ContactModel #vRegulator goes here
    shock::A
    steering::S
    braking::B
end

struct ContactModel
    μ_roll::Float64
    μ_skid::Float64
    k_p::Float64
    k_i::Float64
end
X(::ContactModel) = ComponentVector(x = 0.0, y = 0.0) #v regulator integrator states
Y(::ContactModel) = ComponentVector(v = zeros(2), s = zeros(2), α_p = zeros(2), α_i = zeros(2),
                                    α_raw = zeros(2), α_sat = zeros(2), α = zeros(2))
U(::ContactModel) = nothing

abstract type AbstractSteering end

struct NoSteering <: AbstractSteering end
X(::NoSteering) = nothing
Y(::NoSteering) = nothing
U(::NoSteering) = nothing

struct DirectSteering <: AbstractSteering
    limits::SVector{2,Float64}
end
X(::DirectSteering) = nothing
Y(::DirectSteering) = 0.0 #steering angle
U(::DirectSteering) = 0.0 #steering angle input

struct ActuatedSteering{A <: AbstractActuator} <: AbstractSteering
    limits::SVector{2,Float64} #more generally, transmission kinematics could go here
    actuator::A #actuator model parameters go here
end
X(steering::ActuatedSteering) = X(steering.actuator)
Y(steering::ActuatedSteering) = Y(steering.actuator)
U(steering::ActuatedSteering) = U(steering.actuator) #typically, 1

abstract type AbstractBraking end

struct NoBraking <: AbstractBraking end
X(::NoBraking) = nothing
Y(::NoBraking) = nothing
U(::NoBraking) = nothing

struct DirectBraking <: AbstractBraking
    efficiency::Float64
end
X(::DirectBraking) = nothing
Y(::DirectBraking) = 0.0 #braking strength
U(::DirectBraking) = 0.0 #braking strength input

function X(ldg::LandingGearLeg{A,S,B}) where {A,S,B}
    x_blocks = Dict(:Symbol, Any)
    push!(x_blocks, :contact => X(ldg.contact))
    push!(x_blocks, :shock => X(ldg.shock))
    push!(x_blocks, :steering => X(ldg.steering))
    push!(x_blocks, :braking => X(ldg.braking))
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
could #also #account for other braking model parameters

=#


#################### AbstractSystem interface
X(::LandingGearLeg) = ComponentVector(state = 0.0)
Y(::LandingGearLeg) = ComponentVector(output = 0.0) #both are valid, but not having missing elements in the overall Y vector is slightly faster for some operations
U(::LandingGearLeg) = missing
f_cont!(y, ẋ, x, u, t, ldg::LandingGearLeg, trn = nothing) = (ẋ.state = 0.001x.state)
f_disc!(x, u, t, ldg::LandingGearLeg, trn = nothing) = false

#################### AbstractComponent interface
get_wr_Ob_b(y, comp::LandingGearLeg) = Wrench()
get_h_Gc_b(y, comp::LandingGearLeg) = SVector(0.0, 0, 0)



struct LandingGearGroup{C} <: AbstractComponentGroup{C} end

function LandingGearGroup(nt::NamedTuple{L, T}  where {L, T<:NTuple{N,AbstractLandingGearLeg} where {N}})
    LandingGearGroup{nt}()
end

end


#in System, define and extend f_branch!

# #individual Component
# f_branch!(y, dx, x, u, t, sys, args...) = f_branch!(Val(has_input(sys)), y, dx, x, u, t, args...)
# f_branch!(::Val{true}, y, dx, x, u, t, sys, args...) = f_cont!(y, dx, x, u, t, sys, args...)
# f_cont!(::HasInput, y, dx, x ,u, t, sys, args...) = f_cont!(y, dx, x, u, t, sys, args...)
# f_cont!(::HasNoInput, y, dx, x, u, t, sys, args...) = f_cont!(y, dx, x, t, sys, args...)

# #for a ComponentGroup
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