module FADS

using Random
using ComponentArrays
using StaticArrays
using UnPack
using Flight

export AirflowDynamics, AirflowProperties
export PressureSensor, PressureMeasurement, SensorArray

################################ AirflowDynamics ###############################

Base.@kwdef struct AirflowDynamics <: Component
    p::DoubleIntegrator{1} = DoubleIntegrator{1}(; k_u = 1.0, k_av = -0.3, k_ap = 0)
    # M::DoubleIntegrator
    α::DoubleIntegrator{1} = DoubleIntegrator{1}(; k_u = 1.0, k_av = -0.3, k_ap = -0.02)
    # β::DoubleIntegrator
end

# #we can use the fallback f_disc!, but we need to handle randn! manually
# function Random.randn!(rng::AbstractRNG, sys::System{<:AirflowDynamics})
#     map(s -> randn!(rng, s), values(sys.subsystems))
# end

############################### AirflowProperties ##############################

Base.@kwdef struct AirflowProperties
    p::Float64 = 0.0
    M::Float64 = 0.0
    α::Float64 = 0.0
    β::Float64 = 0.0
end

abstract type AbstractLocation end
struct LocationA <: AbstractLocation end
struct LocationB <: AbstractLocation end

function AirflowProperties(sys::System{<:AirflowDynamics})
    AirflowProperties(p = sys.y.p.p[1], α = sys.y.α.p[1])
end

function local_pressure(airflow::AirflowProperties, loc::Symbol)
    loc === :a && return local_pressure(airflow, LocationA())
    loc === :b && return local_pressure(airflow, LocationB())
    error("Invalid location")
end

function local_pressure(airflow::AirflowProperties, ::LocationA)
    @unpack p, M, α, β = airflow
    return p #just return the true p∞ for now
end

function local_pressure(airflow::AirflowProperties, ::LocationB)
    @unpack p, M, α, β = airflow
    return 2p
end

############################## PressureSensor ##################################

Base.@kwdef struct PressureSensor <: Component
    loc::Symbol
    noise::DiscreteGWN{1} = DiscreteGWN{1}(σ = [0.05])
end

Base.@kwdef struct PressureSensorY #drop the noise's output
    valid::Bool = true
    value::Float64 = 0.0
end

# Systems.init(::SystemU, cmp::PressureSensor) = (fail = Ref(false), noise = init_u(cmp.noise))
Systems.init(::SystemS, ::PressureSensor) = Ref(true) #healthy?
Systems.init(::SystemY, ::PressureSensor) = PressureSensorY()

function measure(sys::System{<:PressureSensor}, airflow::AirflowProperties, args...)
    @unpack s, params, noise = sys
    valid = s[]
    p_loc = local_pressure(airflow, params.loc)
    value = valid * p_loc + Stochastic.sample(noise, args...)[1]
    return PressureSensorY(; valid, value)
end

function Systems.f_disc!(sys::System{<:PressureSensor}, ::Float64,
                         airflow::AirflowProperties, args...)
    #call f_disc! on any stateful modeled errors here
    sys.y = measure(sys, airflow, args...)
    return true
end


################################# SensorArray ##################################

Base.@kwdef struct SensorArray <: Component
    a::PressureSensor = PressureSensor(loc = :a)
    b::PressureSensor = PressureSensor(loc = :b)
end

################################### Model ######################################

Base.@kwdef struct Model <: Component
    airflow::AirflowDynamics = AirflowDynamics()
    sensors::SensorArray = SensorArray()
end

#can we rely on the fallback f_disc!? we can for the method with AbstractRNG but
#for the method with provided u, we need to assign u to the overall Model input.
#each subsystm in the Model hierarchy will get its particular u block updated in
#this step. so it is ready for a f_disc! call without extra args. we can then
#call f_disc! without arguments, which will fall back to the generic method and
#propagate down the hierarchy.



################################# Estimator ####################################

Base.@kwdef struct Estimator <: Component end

################################### World ######################################

Base.@kwdef struct World <: Component
    model::Model = Model()
    estimator::Estimator = Estimator()
end

end #module