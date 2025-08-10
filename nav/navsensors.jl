module NavSensors

############################# IMPLEMENTATION ###################################

using Flight.FlightCore
using Flight.FlightLib

using StaticArrays, ComponentArrays
using LinearAlgebra, DataStructures
using UnPack

################################### IMU ########################################
################################################################################

struct IMUErrors <: ModelDefinition end
struct IMUErrorsY end

@kwdef struct IMUInputs
    q_eb::RQuat = RQuat() #ECEF to vehicle frame rotation
    r_eb_e::SVector{3, Float64} = zeros(SVector{3}) #ECEF position vector, ECEF frame
    ω_eb_b::SVector{3, Float64} = zeros(SVector{3}) #ECEF to vehicle frame angular velocity, vehicle frame
    a_ib_b::SVector{3, Float64} = zeros(SVector{3}) #inertial acceleration of vehicle frame
    α_ib_b::SVector{3, Float64} = zeros(SVector{3}) #ECEF to vehicle frame angular acceleration, vehicle frame
end

function IMUInputs(kin::KinData, dyn::Accelerations)
    @unpack q_eb, r_eb_e, ω_eb_b = kin
    @unpack a_ib_b, α_ib_b = dyn
    IMUInputs(; q_eb, r_eb_e, ω_eb_b, a_ib_b, α_ib_b)
end

@kwdef struct IMUSample
    ω̄_ic_c::SVector{3, Float64} = zeros(SVector{3}) #average angular velocity in case frame
    f̄_c_c::SVector{3, Float64} = zeros(SVector{3}) #average specific force in case frame
    ϑ_c::SVector{3, Float64} = zeros(SVector{3}) #raw angle increment
    ϑ_c_cc::SVector{3, Float64} = zeros(SVector{3}) #coning-corrected angle increment
    υ_c::SVector{3, Float64} = zeros(SVector{3}) #raw velocity increment
    υ_c_sc::SVector{3, Float64} = zeros(SVector{3}) #sculling-corrected velocity increment
end

@kwdef struct IMUOutputs
    sample::IMUSample = IMUSample()
    errors::IMUErrorsY = IMUErrorsY()
end

@kwdef struct IMU <: ModelDefinition
    t_bc::FrameTransform = FrameTransform() #vehicle frame to IMU case frame
    buffer::CircularBuffer{IMUSample} = CircularBuffer{IMUSample}(1)
    errors::IMUErrors = IMUErrors()
end

Modeling.U(::IMU) = Ref(IMUInputs())

Modeling.X(::IMU) = ComponentVector(
    ϑ_c = zeros(3), #raw angle increment
    q_c_cc = zeros(4), #coning-corrected inertial attitude increment
    υ_c = zeros(3), #raw velocity increment
    υ_c_sc = zeros(3), #sculling-corrected velocity increment
)

Modeling.Y(::IMU) = IMUOutputs()

function Modeling.f_ode!(mdl::Model{<:IMU})

    @unpack ẋ, x, u, constants = mdl
    @unpack q_eb, r_eb_e, ω_eb_b, a_ib_b, α_ib_b = u[]
    @unpack t_bc = constants

    q_bc = t_bc.q
    r_bc_b = t_bc.r

    q_c_cc = RQuat(x.q_c_cc, normalization = false)

    #compute geographic position of Oc
    r_bc_e = q_eb(r_bc_b)
    r_ec_e = r_eb_e + r_bc_e
    Oc = Geographic(Cartesian(r_ec_e))

    q_en = ltf(Oc) #local-level frame at Oc
    q_ec = q_eb * q_bc
    q_nc = q_en' ∘ q_ec

    G_c_n = G_n(Oc)
    G_c_c = q_nc'(G_c_n)

    ω_ie_e = SVector{3, Float64}(0, 0, ω_ie) #use WGS84 constant
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b
    ω_ic_b = ω_ib_b
    ω_ic_c = q_bc'(ω_ic_b)

    a_ic_b = a_ib_b + ω_ib_b × (ω_ib_b × r_bc_b) + α_ib_b × r_bc_b
    a_ic_c = q_bc'(a_ic_b)
    f_c_c = a_ic_c - G_c_c

    ẋ.ϑ_c = ω_ic_c
    ẋ.υ_c = f_c_c

    ẋ.q_c_cc = Attitude.dt(q_c_cc, ω_ic_c)
    ẋ.υ_c_sc = q_c_cc(f_c_c)

    #no updating mdl.y here, only in f_disc!

end

function Modeling.f_disc!(::NoScheduling, mdl::Model{<:IMU})

    @unpack ẋ, x, u, y, Δt, constants, submodels = mdl
    buffer = constants.buffer
    errors = submodels.errors

    ϑ_c = SVector{3}(x.ϑ_c)
    ϑ_c_cc = RVec(RQuat(x.q_c_cc, normalization = false))[:]
    υ_c = SVector{3}(x.υ_c)
    υ_c_sc = SVector{3}(x.υ_c_sc)

    #reset integrated quantities
    x.ϑ_c[:] .= 0
    x.q_c_cc[:] .= RQuat()[:] #identity quaternion
    x.υ_c[:] .= 0
    x.υ_c_sc[:] .= 0

    #TODO: apply error model to extracted integrated quantities

    #compute average angular velocity specific force over sampling interval
    ω̄_ic_c = ϑ_c / Δt
    f̄_c_c = υ_c / Δt

    #update model outputs
    sample = IMUSample(; ω̄_ic_c, f̄_c_c, ϑ_c, ϑ_c_cc, υ_c, υ_c_sc)
    push!(buffer, sample)

    f_disc!(errors)
    mdl.y = IMUOutputs(sample, errors.y)

end

############################# TESTS ###################################

using Test

function test_navsensors()
    # imu = IMU() |> Model
    mdl = IMUTestHarness() |> Model;
    return mdl

    #set up simulation with dt < Δt


end

end