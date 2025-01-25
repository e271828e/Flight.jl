module TestNavigationSensors

############################# IMPLEMENTATION ###################################

using Flight.FlightCore
using Flight.FlightPhysics

using StaticArrays, ComponentArrays
using LinearAlgebra, DataStructures
using UnPack

################################### IMU ########################################
################################################################################

struct IMUErrors <: SystemDefinition end
struct IMUErrorsY end

@kwdef struct IMUInputs
    q_eb::RQuat = RQuat() #ECEF to vehicle frame rotation
    q_nb::RQuat = RQuat() #NED to vehicle frame rotation
    r_eOb_e::SVector{3, Float64} = zeros(SVector{3}) #ECEF position vector, ECEF frame
    ω_eb_b::SVector{3, Float64} = zeros(SVector{3}) #ECEF to vehicle frame angular velocity, vehicle frame
    a_iOb_b::SVector{3, Float64} = zeros(SVector{3}) #inertial acceleration of Ob, vehicle frame
    α_ib_b::SVector{3, Float64} = zeros(SVector{3}) #ECEF to vehicle frame angular acceleration, vehicle frame
end

function IMUInputs(kin::KinData, dyn::DynDataOut)
    @unpack q_eb, q_nb, r_eOb_e, ω_eb_b = kin
    @unpack a_iOb_b, α_ib_b = dyn
    IMUInputs(; q_eb, q_nb, r_eOb_e, ω_eb_b, a_iOb_b, α_ib_b)
end

@kwdef struct IMUSample
    ω̄_ic_c::SVector{3, Float64} = zeros(SVector{3}) #average angular velocity in case frame
    f̄_Oc_c::SVector{3, Float64} = zeros(SVector{3}) #average specific force in case frame
    ϑ_c::SVector{3, Float64} = zeros(SVector{3}) #raw angle increment
    ϑ_c_cc::SVector{3, Float64} = zeros(SVector{3}) #coning-corrected angle increment
    υ_c::SVector{3, Float64} = zeros(SVector{3}) #raw velocity increment
    υ_c_sc::SVector{3, Float64} = zeros(SVector{3}) #sculling-corrected velocity increment
end

@kwdef struct IMUOutputs
    sample::IMUSample = IMUSample()
    errors::IMUErrorsY = IMUErrorsY()
end

@kwdef struct IMU <: SystemDefinition
    t_bc::FrameTransform = FrameTransform() #vehicle frame to IMU case frame
    buffer::CircularBuffer{IMUSample} = CircularBuffer{IMUSample}(1)
    errors::IMUErrors = IMUErrors()
end

Systems.U(::IMU) = Ref(IMUInputs())

Systems.X(::IMU) = ComponentVector(
    ϑ_c = zeros(3), #raw angle increment
    q_c_cc = zeros(4), #coning-corrected inertial attitude increment
    υ_c = zeros(3), #raw velocity increment
    υ_c_sc = zeros(3), #sculling-corrected velocity increment
)

Systems.Y(::IMU) = IMUOutputs()

function Systems.f_ode!(sys::System{<:IMU})

    @unpack ẋ, x, u, constants = sys
    @unpack q_eb, q_nb, r_eOb_e, ω_eb_b, a_iOb_b, α_ib_b = u[]
    @unpack t_bc = constants

    q_bc = t_bc.q
    r_ObOc_b = t_bc.r

    q_c_cc = RQuat(x.q_c_cc, normalization = false)

    #compute geographic position of Oc
    r_ObOc_e = q_eb(r_ObOc_b)
    r_eOc_e = r_eOb_e + r_ObOc_e
    Oc = Cartesian(r_eOc_e)

    q_nc = q_nb ∘ q_bc
    G_Oc_n = G_n(Oc)
    G_Oc_c = q_nc'(G_Oc_n)

    ω_ie_e = SVector{3, Float64}(0, 0, ω_ie) #use WGS84 constant
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b
    ω_ic_b = ω_ib_b
    ω_ic_c = q_bc'(ω_ic_b)

    a_iOc_b = a_iOb_b + ω_ib_b × (ω_ib_b × r_ObOc_b) + α_ib_b × r_ObOc_b
    a_iOc_c = q_bc'(a_iOc_b)
    f_Oc_c = a_iOc_c - G_Oc_c

    ẋ.ϑ_c = ω_ic_c
    ẋ.υ_c = f_Oc_c

    ẋ.q_c_cc = Attitude.dt(q_c_cc, ω_ic_c)
    ẋ.υ_c_sc = q_c_cc(f_Oc_c)

    #no updating sys.y here, only in f_disc!

end

function Systems.f_disc!(::NoScheduling, sys::System{<:IMU})

    @unpack ẋ, x, u, y, Δt, constants, subsystems = sys
    buffer = constants.buffer
    errors = subsystems.errors

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
    f̄_Oc_c = υ_c / Δt

    #update system outputs
    sample = IMUSample(; ω̄_ic_c, f̄_Oc_c, ϑ_c, ϑ_c_cc, υ_c, υ_c_sc)
    push!(buffer, sample)

    f_disc!(errors)
    sys.y = IMUOutputs(sample, errors.y)

end

############################# TESTS ###################################

using Test

@kwdef struct IMUTestHarness <: SystemDefinition
    imu::IMU = IMU()
    kin::LTF = LTF()
    dyn::RigidBodyDynamics = RigidBodyDynamics()
end


# struct IMUTestHarnessY{K}
#     imu::IMUOutputs
#     dynamics::DynDataOut
#     kinematics::K
# end

Systems.U(::IMUTestHarness) = Ref(Wrench())
# Systems.Y(imu::IMUTestHarness) = IMUTestHarnessY()

function Systems.f_ode!(sys::System{<:IMUTestHarness})

    @unpack ẋ, x, u, subsystems = sys
    @unpack imu, kin, dyn = subsystems

    rbd_Ob = RigidBodyDistribution(1.0, diagm(ones(3)))
    mp_Ob = MassProperties(rbd_Ob)
    wr_ext_Ob = u[]

    kin.u .= dyn.x
    f_ode!(kin) #update ẋ and y before extracting kinematics data

    kin_data = KinData(kin)
    rb_data = RigidBodyData(mp_Ob, wr_ext_Ob, zeros(SVector{3}))
    f_ode!(dyn, kin_data, rb_data)

    dyn_data = DynDataOut(dyn)

    imu.u[] = IMUInputs(kin_data, dyn_data)
    f_ode!(imu)

    sys.y = (imu = imu.y, kin = kin.y, dyn = dyn.y)

    return nothing

end

function Systems.f_disc!(::NoScheduling, sys::System{<:IMUTestHarness})

    @unpack imu, kin, dyn = sys.subsystems

    f_disc!(imu) #only imu has discrete dynamics
    sys.y = (imu = imu.y, kin = kin.y, dyn = dyn.y)

    return nothing

end

function test_navsensors()
    # imu = IMU() |> System
    sys = IMUTestHarness() |> System;
    return sys

    #set up simulation with dt < Δt


end

end