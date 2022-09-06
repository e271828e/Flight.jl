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
    # α::SecondOrder
    # β::SecondOrder
    p::SecondOrder{1} = SecondOrder{1}(; k_u = 1.0, k_av = -0.3, k_ap = 0)
    # M::SecondOrder
end


############################### AirflowProperties ##############################

Base.@kwdef struct AirflowProperties
    α::Float64 = 0.0
    β::Float64 = 0.0
    p::Float64 = 0.0
    M::Float64 = 0.0
end

function AirflowProperties(sys::System{<:AirflowDynamics})
    AirflowProperties(p = sys.y.p.p[1])
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
    @unpack α, β, p, M = airflow
    return p #just return the true p∞ for now
end

function local_pressure(airflow::AirflowProperties, ::Val{:b})
    @unpack α, β, p, M = airflow
    return 2p #just return the true p∞ for now
end


################################# SensorArray ##################################

Base.@kwdef struct SensorArray <: Component
    a::PressureSensor = PressureSensor(loc = :a)
    b::PressureSensor = PressureSensor(loc = :b)
end

#we can use the fallback f_disc!, but randn! needs special handling
function Random.randn!(rng::AbstractRNG, sys::System{<:SensorArray})
    map(s -> randn!(rng, s), values(sys.subsystems))
end



end #module