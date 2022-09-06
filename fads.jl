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
    p::GaussianDoubleIntegrator{1} = GaussianDoubleIntegrator{1}(; k_u = 1.0, k_av = -0.3, k_ap = 0)
    # M::GaussianDoubleIntegrator
    α::GaussianDoubleIntegrator{1} = GaussianDoubleIntegrator{1}(; k_u = 1.0, k_av = -0.3, k_ap = -0.02)
    # β::GaussianDoubleIntegrator
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

function AirflowProperties(sys::System{<:AirflowDynamics})
    AirflowProperties(p = sys.y.p.p[1], α = sys.y.α.p[1])
end


############################## PressureSensor ##################################

Base.@kwdef struct PressureSensor <: Component
    loc::Symbol
    noise::DiscreteGaussianWN{1} = DiscreteGaussianWN{1}(σ = [0.05])
end

Base.@kwdef struct PressureSensorY
    valid::Bool = true
    value::Float64 = 0.0
    noise::SVector{1,Float64} = zeros(SVector{1})
end

Systems.init(::SystemU, cmp::PressureSensor) = (fail = Ref(false), noise = init_u(cmp.noise))
Systems.init(::SystemS, ::PressureSensor) = (fail = Ref(false),) #failed?
Systems.init(::SystemY, ::PressureSensor) = PressureSensorY()

function Systems.f_disc!(sys::System{<:PressureSensor}, Δt::Float64, airflow::AirflowProperties, noise_args...)

    @unpack s, u, noise = sys

    #update noise y
    f_disc!(noise, Δt, noise_args...) #args can be a RNG or a N(0,1) sample

    s.fail[] = u.fail[]
    healthy = !s.fail[]
    p_loc = local_pressure(airflow, sys.params.loc)
    ỹ = p_loc * healthy + noise.y[1]

    sys.y = PressureSensorY(; valid = healthy, value = ỹ, noise = noise.y)
    return true
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

Base.@kwdef struct SensorArray <: Component
    a::PressureSensor = PressureSensor(loc = :a)
    b::PressureSensor = PressureSensor(loc = :b)
end


################################# Estimator ####################################

Base.@kwdef struct Estimator <: Component end

################################### World ######################################

Base.@kwdef struct World <: Component
    airflow::AirflowDynamics = AirflowDynamics()
    sensors::SensorArray = SensorArray()
    estimator::Estimator = Estimator()
end

function Systems.f_disc!(sys::System{<:World}, Δt::Real, rng::AbstractRNG)
    f_disc!(sys.airflow, Δt, rng)
    f_disc!(sys.sensors, Δt, rng)

end

end #module