module NavigationSensors

using Flight.FlightCore

using StaticArrays

################################### IMU ########################################
################################################################################

struct IMUErrors <: SystemDefinition end

@kwdef struct IMU <: SystemDefinition
    t_bc::FrameTransform = FrameTransform() #vehicle frame to IMU frame
    errors::IMUErrors = IMUErrors()
end

@kwdef struct IMUInputs
    q_eb::RQuat = RQuat() #ECEF to vehicle frame rotation
    q_nb::RQuat = RQuat() #NED to vehicle frame rotation
    r_eOb_e::SVector{3, Float64} = zeros(SVector{3}) #vehicle frame ECEF position vector
    ω_ib_b::SVector{3, Float64} = zeros(SVector{3}) #vehicle frame inertial angular velocity
    a_iOb_b::SVector{3, Float64} = zeros(SVector{3}) #vehicle frame inertial acceleration
end

@kwdef struct IMUOutputs
    ω_ic_c::SVector{3, Float64} = zeros(SVector{3}) #angular velocity in IMU frame
    f_Oc_c::SVector{3, Float64} = zeros(SVector{3}) #specific force in IMU frame
    ϑ_c::SVector{3, Float64} = zeros(SVector{3}) #delta angle
    υ_iOc_c::SVector{3, Float64} = zeros(SVector{3}) #delta velocity
    β_c::SVector{3, Float64} = zeros(SVector{3}) #coning correction
    ζ_Oc_c::SVector{3, Float64} = zeros(SVector{3}) #sculling correction
end

Systems.U(::System{<:IMU}) = Ref(IMUInputs())

Systems.X(::System{<:IMU}) = ComponentVector(
    ϑ_c = zeros(3), #delta angle
    υ_iOc_c = zeros(3), #delta velocity
    β_c = zeros(3), #coning correction
    ζ_Oc_c = zeros(3), #sculling correction
)

Systems.Y(::System{<:IMU}) = IMUOutputs()

#on each call to Sensors f_ode!
function Systems.f_ode!(sys::System{<:IMU})

    @unpack ẋ, x, u, constants = sys
    @unpack q_eb, q_nb, r_eOb_e, ω_ib_b, a_iOb_b = u
    @unpack t_bc = constants

    q_bc = t_bc.q
    r_ObOc_b = t_bc.r

    q_ec = q_eb ∘ q_bc
    q_nc = q_nb ∘ q_bc


    r_eOc_e = r_eOb_e + q_eb(r_ObOc_b) #ECEF position of IMU frame origin Oc

    #TODO: implement vehicle to case frame transform
    ω_ic_b = ω_ib_b
    ω_ic_c = q_bc'(ω_ic_b) #inertial angular velocity of IMU frame
    a_iOc_c = a_iOb_b

    G_Oc_n = G_n(r_eOc_e) #gravitational acceleration at Oc in NED frame
    G_Oc_c = q_nc'(G_Oc_n) #gravitational acceleration at Oc in IMU frame
    f_Oc_c = a_iOc_c - G_Oc_c #specific force at Oc in IMU frame

    ϑ_c = SVector{3}(x.ϑ_c)
    υ_iOc_c = SVector{3}(x.υ_iOc_c)
    β_c = SVector{3}(x.β_c)
    ζ_Oc_c = SVector{3}(x.ζ_Oc_c)

    ẋ.ϑ_c = ω_ic_c
    ẋ.υ_iOc_c = f_Oc_c
    ẋ.β_c = 0.5 * ϑ_c × ω_ic_c
    ẋ.ζ_Oc_c = 0.5 * (ϑ_c × f_Oc_c + υ_iOc_c × ω_ic_c)







end


#this goes in IMURoot
    # buffer::CircularBuffer{IMUOutput} = CircularBuffer{IMUOutput}(1)


#on each IMUCore f_disc! call,


end