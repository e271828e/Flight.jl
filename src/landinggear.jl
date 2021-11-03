module LandingGear
#rename to ground, rename LandingGearLeg to LandingGear

using StaticArrays
using ComponentArrays
using UnPack

using Flight.ModelingTools
# using Flight.Airdata
using Flight.Kinematics
using Flight.Dynamics
using Flight.Airframe
using Flight.Terrain

#for method extension without module qualifier
import Flight.ModelingTools: System, get_x0, get_y0, get_u0, get_d0,f_cont!, f_disc!
import Flight.Airframe: get_wr_b, get_hr_b

export LandingGearLeg, ContactModel
export NoSteering, DirectSteering, get_steering_angle, set_steering_input
export NoBraking, DirectBraking, get_braking_coefficient, set_braking_input

"""
SET NORMALIZATION TO false!!!!!!!
"""


########################### ContactModel #############################

abstract type AbstractContactModel <: AbstractComponent end

struct FrictionCoefficients
    static::Float64
    dynamic::Float64
end

SkiddingFriction = @NamedTuple begin
    dry::FrictionCoefficients
    wet::FrictionCoefficients
    icy::FrictionCoefficients
end

#wet: .25, .15, #icy: 0.075, 0.05

Base.@kwdef struct ContactModel <: AbstractContactModel
    μ_roll::FrictionCoefficients = FrictionCoefficients(0.03, 0.02)
    μ_skid::FrictionCoefficients = FrictionCoefficients(0.75, 0.5)
    v_bo::SVector{2,Float64} = [0.005, 0.01] #breakout velocity interval
    ψ_skid::Float64 = deg2rad(10)
    k_p::Float64 = 5.0
    k_i::Float64 = 400.0
end

Base.@kwdef struct ContactModelY
    v::SVector{2,Float64} = zeros(2)
    s::SVector{2,Float64} = zeros(2)
    α_p::SVector{2,Float64} = zeros(2)
    α_i::SVector{2,Float64} = zeros(2)
    α_raw::SVector{2,Float64} = zeros(2)
    α_sat::SVector{2,Bool} = zeros(2)
    α::SVector{2,Float64} = zeros(2)
    ψ_slip::Float64 = 0.0
    k_bo::Float64 = 0.0
    μ_max::SVector{2,Float64} = zeros(2)
    μ::SVector{2,Float64} = zeros(2)
end

get_x0(::ContactModel) = ComponentVector(x = 0.0, y = 0.0) #v regulator integrator states
get_y0(::ContactModel) = ContactModelY()
f_cont!(::System{ContactModel}, args...) = nothing


########################### Damper #############################

abstract type AbstractDamper end #not a System!
get_force(::A) where {A<:AbstractDamper} = no_extend_error(get_force, A)

Base.@kwdef struct SimpleDamper <: AbstractDamper
    l_0::Float64 = 1 #natural length
    k_s::Float64 = 25000 #spring constant
    k_d_ext::Float64 = 1000 #extension damping coefficient
    k_d_cmp::Float64 = 1000 #compression damping coefficient
    ξ_min::Float64 = -1 #compression below which the shock absorber is disabled
end

f_cont!(::System{SimpleDamper}, args...) = nothing

########################### Steering #############################

abstract type AbstractSteering <: AbstractComponent end

############## NoSteering ##############

struct NoSteering <: AbstractSteering end
set_steering_input(::System{NoSteering}, ::Float64) = nothing
get_steering_angle(::System{NoSteering}) = 0.0
f_cont!(::System{NoSteering}, args...) = nothing


############## DirectSteering ##############

Base.@kwdef struct DirectSteering <: AbstractSteering
    ψ_max::Float64 = π/6
end

Base.@kwdef struct DirectSteeringY
    ψ::Float64 = 0.0
end
#we need to make the contents of u mutable
get_u0(::DirectSteering) = Ref(0.0) #steering input
get_y0(::DirectSteering) = DirectSteeringY(0.0) #steering angle

function f_cont!(sys::System{DirectSteering})
    u = sys.u[]
    abs(u) <= 1 ? sys.y = DirectSteeringY(u * sys.params.ψ_max) :
        throw(ArgumentError("Steering input must be within [-1,1]"))
end

# struct ActuatedSteering{A <: AbstractActuator} <: AbstractSteering
#     limits::SVector{2,Float64} #more generally, transmission kinematics could go here
#     actuator::A #actuator model parameters go here
# end
# get_x0(steering::ActuatedSteering) = get_x0(steering.actuator)
# get_y0(steering::ActuatedSteering) = get_x0(steering.actuator)
# get_u0(steering::ActuatedSteering) = get_u0(steering.actuator) #typically, 1

########################### Braking #############################

abstract type AbstractBraking <: AbstractComponent end

############## NoBraking ###############

struct NoBraking <: AbstractBraking end
f_cont!(::System{NoBraking}, args...) = nothing
set_braking_input(::System{NoBraking}, ::Float64) = nothing
get_braking_coefficient(::System{NoBraking}) = 0.0

########### DirectBraking #############

Base.@kwdef struct DirectBraking <: AbstractBraking
    η_brk::Float64 = 1.0 #braking efficiency
end

Base.@kwdef struct DirectBrakingY
    k_brk::Float64 = 0.0 #braking coefficient
end

get_u0(::DirectBraking) = Ref(0.0)
get_y0(::DirectBraking) = DirectBrakingY()

function f_cont!(sys::System{DirectBraking})
    u = sys.u[]
    (u <= 1 && u >= 0) ? sys.y = DirectBrakingY(u * sys.params.η_brk) :
        throw(ArgumentError("Braking Input must be within [0,1]"))
end

########################## LandingGearLeg #########################

Base.@kwdef struct LandingGearLeg{D<:AbstractDamper, C<:AbstractContactModel,
        S <: AbstractSteering, B <: AbstractBraking} <: AbstractAirframeComponent
    frame::FrameSpec = FrameSpec()
    damper::D = SimpleDamper()
    contact::C = ContactModel()
    steering::S = NoSteering()
    braking::B = NoBraking()
end

struct LandingGearLegY{ContactY, SteeringY, BrakingY}
    contact::ContactY
    steering::SteeringY
    braking::BrakingY
end

#maybe put all subcomponents in a single field:
# Base.@kwdef struct LandingGearLegSubcomponents{C, A, S, B}
#     contact::C = ContactModel()
#     shock::A = ShockAbsorber()
#     steering::S = NoSteering()
#     braking::B = NoBraking()
# end

# const ldg_subsystems = (:contact, :shock, :steering, :braking)
# function get_u0(ldg::LandingGearLeg)
#     NamedTuple{ldg_ss_labels}(map(z -> get_u0(getproperty(ldg, z)), ldg_ss_labels))
# end

get_x0(ldg::LandingGearLeg) = ComponentVector(
    contact = get_x0(ldg.contact),
    steering = get_x0(ldg.steering),
    braking = get_x0(ldg.braking)
    )

#need the type definition for plot dispatch
get_y0(ldg::LandingGearLeg) = LandingGearLegY(
    get_y0(ldg.contact),
    get_y0(ldg.steering),
    get_y0(ldg.braking)
    )

get_u0(ldg::LandingGearLeg) = (
    contact = get_u0(ldg.contact),
    steering = get_u0(ldg.steering),
    braking = get_u0(ldg.braking)
    )

get_d0(ldg::LandingGearLeg) = (
    contact = get_d0(ldg.contact),
    steering = get_d0(ldg.steering),
    braking = get_d0(ldg.braking)
    )


function System(ldg::LandingGearLeg, ẋ = get_x0(ldg), x = get_x0(ldg),
                    y = get_y0(ldg), u = get_u0(ldg), d = get_d0(ldg), t = Ref(0.0))

    ss_list = Vector{System}()
    ss_labels = (:contact, :steering, :braking)
    for label in ss_labels
        push!(ss_list, System(map((λ)->getproperty(λ, label), (ldg, ẋ, x, y, u, d))..., t))
    end

    params = (frame = ldg.frame, damper = ldg.damper)
    subsystems = NamedTuple{ss_labels}(ss_list)

    System{map(typeof, (ldg, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)
end

function f_cont!(sys::System{<:LandingGearLeg}, kin::KinData, trn::AbstractTerrainModel)

    @unpack ẋ, x, u, params, subsystems = sys
    @unpack frame, damper = params
    @unpack contact, steering, braking = subsystems

    ##NEED TO DETERMINE THE RIGHT EXECUTION ORDER





    #update steering and braking subsystems, which we can assume independent
    #from contact status
    f_cont!(steering)
    f_cont!(braking)


    #with the kinematics and the strut frame, determine WoW state

    #the contact model gives us the normalized contact force in the contact
    #frame
    f_cont!(contact) #THIS SHOULD DEPEND ON WOW
    #add f_disc! to reset the integrators for contact model

    #now, from the kinematics of the contact point and the oleo parameters, we
    #get F_shock_zs

    sys.y = LandingGearLegY(contact.y, steering.y, braking.y)

end

f_disc!(::System{<:LandingGearLeg}) = false

#CHANGE THIS
#CHANGE THIS
#CHANGE THIS
#CHANGE THIS
get_wr_b(::System{<:LandingGearLeg}) = Wrench()
get_hr_b(::System{<:LandingGearLeg}) = zeros(SVector{3})


end #module

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