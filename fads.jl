module FADS

using Random
using ComponentArrays
using StaticArrays
using UnPack
using Flight
using Flight.Components.Discrete: OrnsteinUhlenbeck, SecondOrder

export AirflowDynamics, AirflowProperties
export SensorArray, PressureSensor, PressureMeasurement


################################ AirflowDynamics ###############################

Base.@kwdef struct AirflowDynamics <: Component
    p::SecondOrder{1} = SecondOrder{1}(; k_u = 1.0, k_av = -0.3, k_ap = 0)
    # M::SecondOrder
    α::SecondOrder{1} = SecondOrder{1}(; k_u = 1.0, k_av = -0.3, k_ap = -0.02)
    # β::SecondOrder
end


############################### AirflowProperties ##############################

Base.@kwdef struct AirflowProperties
    p::Float64 = 0.0
    M::Float64 = 0.0
    α::Float64 = 0.0
    β::Float64 = 0.0
end

function AirflowProperties(sys::System{<:AirflowDynamics})
    AirflowProperties(p = sys.y.p.p[1], α = sys.y.α.p)
end


############################## PressureSensor ##################################

Base.@kwdef struct PressureSensor <: Component
    loc::Symbol
    σ_y::Float64 = 0.05
end

Base.@kwdef struct PressureSensorY
    valid::Bool = true
    value::Float64 = 0.0
end

Systems.init(::SystemU, ::PressureSensor) = (fail = Ref(false), w = Ref(0.0)) #noise input
Systems.init(::SystemS, ::PressureSensor) = (fail = Ref(false),) #failed?
Systems.init(::SystemY, ::PressureSensor) = PressureSensorY()


#just update the fail state
function Systems.f_disc!(sys::System{<:PressureSensor}, ::Float64, airflow::AirflowProperties)
    sys.s.fail[] = sys.u.fail[]
    sys.y = measure(sys, airflow)
    return true
end

#randomize noise inputs
Random.randn!(rng::AbstractRNG, sys::System{<:PressureSensor}) = (sys.u.w[] = randn(rng))

function measure(sys::System{<:PressureSensor}, airflow::AirflowProperties)
    @unpack u, s, params = sys
    @unpack σ_y, loc = params

    healthy = !s.fail[]
    p_loc = local_pressure(airflow, loc)
    ỹ = p_loc * healthy + σ_y * u.w[]

    PressureSensorY(; valid = healthy, value = ỹ)

end

local_pressure(airflow::AirflowProperties, loc::Symbol) = local_pressure(airflow, Val(loc))

function local_pressure(airflow::AirflowProperties, ::Val{:a})
    @unpack p, M, α, β = airflow
    return p #just return the true p∞ for now
end

function local_pressure(airflow::AirflowProperties, ::Val{:b})
    @unpack p, M, α, β = airflow
    return 2p
end


################################# SensorArray ##################################

Base.@kwdef struct SensorSet <: Component
    a::PressureSensor = PressureSensor(loc = :a)
    b::PressureSensor = PressureSensor(loc = :b)
end

#we can use the fallback f_disc!, but randn! needs special handling
function Random.randn!(rng::AbstractRNG, sys::System{<:SensorSet})
    map(s -> randn!(rng, s), values(sys.subsystems))
end


################################# Estimator ####################################

Base.@kwdef struct Estimator <: Component end

################################### World ######################################

Base.@kwdef struct World <: Component
    airflow::AirflowDynamics = AirflowDynamics()
    sensors::SensorSet = SensorSet()
    estimator::Estimator = Estimator()
end

end #module