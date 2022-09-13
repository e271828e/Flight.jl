module FADS

using Random
using ComponentArrays
using StaticArrays
using UnPack
using Flight
using Flight.Stochastic: StochasticProcess, σ²0

export AirflowDynamics, AirflowProperties
export PressureSensor, PressureMeasurement, SensorArray

################################ AirflowDynamics ###############################

Base.@kwdef struct AirflowDynamics <: StochasticProcess
    p::DoubleIntegrator{1} = DoubleIntegrator{1}(; k_u = 1.0, k_av = -0.3, k_ap = 0)
    # M::DoubleIntegrator
    α::DoubleIntegrator{1} = DoubleIntegrator{1}(; k_u = 1.0, k_av = -0.3, k_ap = -0.02)
    # β::DoubleIntegrator
end

function Stochastic.σ²0(sys::System{<:AirflowDynamics})
    ss = sys.subsystems
    (k => σ²0(v) for (k, v) in zip(keys(ss), values(ss))) |> OrderedDict |> ComponentVector
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

Base.@kwdef struct PressureSensor <: StochasticProcess
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

Base.@kwdef struct SensorArray <: StochasticProcess
    a::PressureSensor = PressureSensor(loc = :a)
    b::PressureSensor = PressureSensor(loc = :b)
end

########################### OverallModel ################################

Base.@kwdef struct Model <: StochasticProcess
    airflow::AirflowDynamics = AirflowDynamics()
    sensors::SensorArray = SensorArray()
end

function Stochastic.σ²0(sys::System{<:Model})
    #ComponentVector(airflow = σ²0(model.airflow), sensors = σ²(model.sensors))
    σ²0(sys.airflow) #only AirflowDynamics has states for now
end

function Systems.f_disc!(sys::System{<:Model}, Δt::Real,
                         u::Union{Real, AbstractVector{<:Real}})
    #assign the overall Model's noise input, then fall back to the generic,
    #recursive no-argument f_disc! method
    sys.u .= u
    f_disc!(sys, Δt)
end

################################################################################
############################## FADSFilter ######################################

Base.@kwdef struct FADSFilter{  D <: System{<:Model},
                                P <: SRUKF.StatePropagator,
                                M <: SRUKF.MeasurementProcessor} <: Component
    model::D #estimator Model
    sp::P #StatePropagator for the estimator Model
    mp_p::M #MeasurementProcessor for static pressure measurements
end

#sp and mp_p are sized during FADSFilter construction (not System!) from the
#model's x and u.

function Systems.init(::SystemS, fads::FADSFilter)

    @unpack model = fads
    LW = length(model.u) #noise length is given by Model input length

    #state SR-covariance for the filter model
    P_δx = diagm(σ²0(model))
    S_δx = cholesky(P_δx).L

    #normalized noise SR-covariance matrix for the filter model
    P_δw = SizedMatrix{LW, LW}(1.0I)
    S_δw = cholesky(P_δw).L

    #normalized noise SR-covariance matrix for static pressure measurements
    P_δv_p = SizedMatrix{1, 1}(1.0I)
    S_δv_p = cholesky(P_δv_p).L

    return (x̄ = x̄, S_δx = S_δx, S_δw = S_δw, S_δv_p = S_δv_p)

end

function Systems.f_disc!(sys::System{<:Filter}, Δt::Real)
end

function propagate!(filter::System{<:FADSFilter}, Δt::Real)

    @unpack model, sp_p = filter.params
    @unpack x̄, S_δx, S_δw = filter.s

    f! = let model = model, Δt = Δt
        function (x1, x0, w)
            #assign overall model state, including airflow and all sensor states
            model.x .= x0
            #call discrete dynamics with the normalized noise sample from the
            #unscented transform
            f_disc!(model, Δt, w)
            #retrieve the updated model state
            x1 .= model.x
        end
    end

    SRUKF.propagate!(sp_p, x̄, S_δx, S_δw_p, f!)

end

function correct!(filter::System{<:FADSFilter},
    measurements::NamedTuple{S, NTuple{N, PressureSensorY}}) where {S,N}

    @unpack model, mp_p = filter.params
    @unpack x̄, S_δx, S_δv_p = filter.s

    for (label, measurement) in measurements

        h! = let model = model, label = label
            function (y, x, v)
                #assign model state (includes airflow and any modeled sensor
                #states)
                model.x .= x
                #compute the airflow properties corresponding to the model state
                airflow = AirflowProperties(model.airflow)
                #retrieve the subsystem corresponding to the sensor that took
                #the measurement sample
                sensor = getproperty(model.sensors, label)
                #compute and assign the predicted measurement using the
                #single-sensor noise input from the unscented transform
                y .= measure(sensor, airflow, v).value
            end
        end

        measurement.valid || continue #if invalid, skip to the next measurement
        ỹ = measurement.value
        SRUKF.apply!(mp_p, x̄, S_δx, S_δv_p, ỹ, h!)

    end
end


################################### World ######################################

Base.@kwdef struct World{M, E} <: Component
    model::M = Model()
    estimator::E = Filter()
end

end #module