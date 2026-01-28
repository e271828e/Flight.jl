module Robot2D

using LinearAlgebra, StaticArrays, ComponentArrays

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightLib.Control.Discrete: LQR, PID, LQROutput, PIDOutput

const g = 9.80665 #m/s^2, standard gravity

################################################################################
################################## Vehicle2 #####################################

#body 1: main body, comprising the vehicle chassis and the DC motor's case and stator
#body 2: rolling body, comprising the wheels, axle, and the DC motor's rotor.

@kwdef struct Vehicle <: ModelDefinition
    L::Float64 = 0.15 #distance from main body's origin to its CoM (m)
    R::Float64 = 0.05 #wheel radius (m)
    m_1::Float64 = 1.0 #mass of main body (kg)
    m_2::Float64 = 0.1 #mass of rolling body (kg)
    J_1::Float64 = 1/12 * m_1 * (2L)^2 #main body's moment of inertia with respect to its CoM (kg*m^2)
    J_2::Float64 = 1/2 * m_2 * R^2 #rolling body's moment of inertia with respect to its CoM (kg*m^2)
    k_m::Float64 = 0.32 #motor's torque constant (N*m)
    b_m::Float64 = 0.0189 #motor's effective damping coefficient (N*m*s/rad)
    J_m::Float64 = 0.0014 #motor's effective moment of inertia (kg*m^2)
end

@kwdef struct VehicleY
    ω::Float64 = 0.0 #angular velocity of main body (rad/s)
    v::Float64 = 0.0 #horizontal velocity of the system's origin (m/s)
    θ::Float64 = 0.0 #angle of main body with respect to vertical (rad)
    η::Float64 = 0.0 #horizontal position of the system's origin (m)
    u_m::Float64 = 0.0 #motor control input (normalized)
    τ_m::Float64 = 0.0 #motor torque (N*m)
    ω_dot::Float64 = 0.0 #angular acceleration of main body (rad/s^2)
    v_dot::Float64 = 0.0 #linear acceleration of the system's origin (m/s^2)
end

Modeling.X(::Vehicle) = ComponentVector( ω = 0.0, v = 0.0, θ = 0.0, η = 0.0)
Modeling.U(::Vehicle) = Ref(Ranged(0.0, 0., 1.))
Modeling.Y(::Vehicle) = VehicleY()

function Modeling.f_ode!(mdl::Model{Vehicle})

    (; ẋ, x, u, parameters) = mdl
    (; ω, v, θ, η) = x
    (; L, R, m_1, m_2, J_1, J_2, k_m, b_m, J_m) = parameters

    u_m = Float64(u[])
    ω_m = v / R - ω
    τ_ss = k_m * u_m - b_m * ω_m #steady-state motor torque

    sθ = sin(θ)
    cθ = cos(θ)

    M_11 = m_1 * L^2 + J_1 + J_m
    M_22 = m_1 + m_2 + (J_2 + J_m) / R^2
    M_12 = m_1 * L * cθ - J_m / R

    M = @SMatrix[
        M_11    M_12
        M_12    M_22
    ]

    b = @SVector[
        -τ_ss + m_1 * L * g * sθ
        τ_ss / R + m_1 * L * ω^2 * sθ
    ]

    ω_dot, v_dot = M\b
    ω_m_dot = v_dot/R - ω_dot

    θ_dot = ω
    η_dot = v

    τ_m = τ_ss - J_m * ω_m_dot

    ẋ.ω = ω_dot
    ẋ.v = v_dot
    ẋ.θ = θ_dot
    ẋ.η = η_dot

    mdl.y = VehicleY(; ω, v, θ, η, u_m, τ_m, ω_dot, v_dot)

end

@no_step Vehicle
@no_periodic Vehicle


################################################################################
################################ Initialization ################################

@kwdef struct InitParameters <: FieldVector{5, Float64}
    u_m::Float64 = 0.0
    ω_dot::Float64 = 0.0
    ω::Float64 = 0.0
    θ::Float64 = 0.0
    η::Float64 = 0.0
end

function Modeling.init!( mdl::Model{Vehicle}, ip::InitParameters = InitParameters())

    (; u_m, ω_dot, ω, θ, η) = ip
    (; x, u, parameters) = mdl
    (; L, R, m_1, m_2, J_1, J_2, k_m, b_m, J_m) = parameters

    sθ = sin(θ)
    cθ = cos(θ)

    M_11 = m_1 * L^2 + J_1 + J_m
    M_22 = m_1 + m_2 + (J_2 + J_m) / R^2
    M_12 = m_1 * L * cθ - J_m / R

    A = @SMatrix[
        1       M_12
        -1/R    M_22
    ]

    b = @SVector[
        m_1 * L * g * sθ - M_11 * ω_dot
        m_1 * L * ω^2 * sθ - M_12 * ω_dot
    ]

    τ_ss, _ = A\b

    ω_m = (k_m * u_m - τ_ss) / b_m
    v = (ω + ω_m) * R

    (x.ω, x.v, x.θ, x.η)  = (ω, v, θ, η)
    u[] = u_m

    f_ode!(mdl)

end

################################################################################
################################ State Space ###################################

@kwdef struct XStateSpace <: FieldVector{4, Float64}
    ω::Float64 = 0.0
    v::Float64 = 0.0
    θ::Float64 = 0.0
    η::Float64 = 0.0
end

@kwdef struct UStateSpace <: FieldVector{1, Float64}
    m::Float64 = 0.0
end

@kwdef struct YStateSpace <: FieldVector{6, Float64}
    ω::Float64 = 0.0
    v::Float64 = 0.0
    θ::Float64 = 0.0
    η::Float64 = 0.0
    u_m::Float64 = 0.0
    τ_m::Float64 = 0.0
end

#get state space system's state vector derivative from Model
function get_ẋ_ss(mdl::Model{Vehicle})
    (; ω, v, θ, η) = mdl.ẋ
    XStateSpace(; ω, v, θ, η)
end

#get state space system's state vector from Model
function get_x_ss(mdl::Model{Vehicle})
    (; ω, v, θ, η) = mdl.x
    XStateSpace(; ω, v, θ, η)
end

#get state space system's input vector from Model
function get_u_ss(mdl::Model{Vehicle})
    UStateSpace(mdl.u[])
end

#get state space system's output vector from Model
function get_y_ss(mdl::Model{Vehicle})
    (; ω, v, θ, η, u_m, τ_m) = mdl.y
    YStateSpace(; ω, v, θ, η, u_m, τ_m)
end

#assign state space system's state vector to Model
function assign_x_ss!(mdl::Model{Vehicle}, x_ss::AbstractVector{<:Real})
    (; ω, v, θ, η) = XStateSpace(x_ss)
    mdl.x.ω, mdl.x.v, mdl.x.θ, mdl.x.η = ω, v, θ, η
end

#assign state space system's input vector to Model
function assign_u_ss!(mdl::Model{Vehicle}, u_ss::AbstractVector{<:Real})
    mdl.u[] = UStateSpace(u_ss[1])[1]
end

#build state space system's update function
function get_f_ss(mdl::Model{Vehicle})
    let mdl = mdl
        function (x::AbstractVector{<:Real}, u::AbstractVector{<:Real})
            assign_x_ss!(mdl, x)
            assign_u_ss!(mdl, u)
            f_ode!(mdl)
            get_ẋ_ss(mdl)
        end
    end
end

#build state space system's output function
function get_h_ss(mdl::Model{Vehicle})
    let mdl = mdl
        function (x::AbstractVector{<:Real}, u::AbstractVector{<:Real})
            assign_x_ss!(mdl, x)
            assign_u_ss!(mdl, u)
            f_ode!(mdl)
            get_y_ss(mdl)
        end
    end
end


################################################################################
################################ Linearization #################################

function Linearization.linearize(mdl::Model{Vehicle}, ip::InitParameters = InitParameters())

    init!(mdl, ip)

    #state space system's functions
    f = get_f_ss(mdl)
    h = get_h_ss(mdl)

    #linearization point
    x0 = get_x_ss(mdl)
    u0 = get_u_ss(mdl)

    linearize(f, h, x0, u0)

end


################################################################################
################################# Controller ###################################

@kwdef struct XController <: FieldVector{3, Float64}
    ω::Float64 = 0.0
    v::Float64 = 0.0
    θ::Float64 = 0.0
end

@kwdef struct ZController <: FieldVector{1, Float64}
    v::Float64 = 0.0
end

@kwdef struct UController <: FieldVector{1, Float64}
    m::Float64 = 0.0
end

function XController(vehicle::Model{<:Vehicle})
    (; ω, v, θ) = vehicle.x
    XController(; ω, v, θ)
end

ZController(vehicle::Model{<:Vehicle}) = ZController(vehicle.x.v)

################################################################################

@enum ControlMode mode_v = 1 mode_η = 2

@kwdef struct Controller{C} <: ModelDefinition
    v2m::C = LQR{3,1,1}()
    η2v::Float64 = 0.0
end

@kwdef mutable struct ControllerU
    mode::ControlMode = mode_v
    v_ref::Float64 = 0.0 #velocity reference
    η_ref::Float64 = 0.0 #position reference
end

@kwdef struct ControllerY{C}
    mode::ControlMode = mode_η
    v_ref::Float64 = 0.0 #velocity reference
    η_ref::Float64 = 0.0 #position reference
    v2m::C = LQROutput{3,1,1}()
end

Modeling.U(::Controller) = ControllerU()
Modeling.Y(::Controller) = ControllerY()

#implement init! to load gains and assign PID Output bounds. But careful, we
#also need to limit input to the LQR even when in velocity mode. So maybe we
#don't need the PID after all



@kwdef struct Robot{C} <: ModelDefinition
    vehicle::Vehicle = Vehicle()
    controller::C = Controller()
end


end #module
