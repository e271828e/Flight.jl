module Robot2D

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore
using Flight.FlightLib

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
Modeling.U(::Vehicle) = Ref{Float64}(0.0)
Modeling.Y(::Vehicle) = VehicleY()

#initialize vehicle at static vertical equilibrium
function Modeling.init!(mdl::Model{Vehicle})
    mdl.x .= 0
    mdl.u[] = 0
    f_ode!(mdl)
end

function Modeling.f_ode!(mdl::Model{Vehicle})

    @unpack ẋ, x, u, parameters = mdl
    @unpack L, R, m_1, m_2, J_1, J_2, k_m, b_m, J_m = parameters

    @unpack ω, v, θ, η = x
    u_m = u[]
    ω_m = v / R - ω
    τ_m_ss = k_m * u_m - b_m * ω_m #steady-state motor torque

    J = m_1 * L^2 + J_1 + J_m
    M = m_1 + m_2 + (J_2 + J_m) / R^2
    K = m_1 * L

    sθ = sin(θ)
    cθ = cos(θ)

    A = @SMatrix[
        J                   K * cθ - J_m / R
        K * cθ - J_m / R    M
    ]

    b = @SVector[
        -τ_m_ss + K * g * sθ
        1/R * τ_m_ss + K * ω^2 * sθ
    ]

    ω_dot, v_dot = A\b
    ω_m_dot = v_dot/R - ω_dot

    θ_dot = ω
    η_dot = v

    τ_m = τ_m_ss - J_m * ω_m_dot

    ẋ.ω = ω_dot
    ẋ.v = v_dot
    ẋ.θ = θ_dot
    ẋ.η = η_dot

    mdl.y = VehicleY(; ω, v, θ, η, u_m, τ_m, ω_dot, v_dot)

end

@no_step Vehicle
@no_periodic Vehicle

function test_vehicle()

    mdl = Vehicle() |> Model

    mdl.u[] = 0
    mdl.x.ω = 0.0
    mdl.x.v = 0
    mdl.x.θ = 0.0
    f_ode!(mdl)
    @show mdl.y

    mdl.u[] = 0.118125
    mdl.x.ω = 0.0
    mdl.x.v = 0.1
    mdl.x.θ = 0.0
    f_ode!(mdl)
    @show mdl.y

    mdl.u[] = 0.247714904
    mdl.x.ω = 0.0
    mdl.x.v = 0.1
    mdl.x.θ = 0.1
    f_ode!(mdl)
    @show mdl.y

    # mdl.u[] = 0.2
    # mdl.x.ω = 0.1
    # mdl.x.v = 0.1
    # mdl.x.θ = 0.05
    # f_ode!(mdl)
    # @show mdl.y

    nothing
end

#unknowns:
#omega1, omega2, theta, eta
#omega1_dot, omega2_dot, theta_dot, eta_dot
#F21x, F21z, Fi2x, Fi2z, tau_m
#u
#total: 14

#equations:
#6 dynamic equations
#1 motor model
#2 kinematic equations for theta_dot = omega1, eta_dot = omega2*R

#total: 9

#we need 5 constraints:
#1: η, which is decoupled and we can set arbitrarily
#2: ω2, which basically determines u
#345: ω1, ω1_dot, θ

#once we set ω1, ω1_dot and θ, ω2_dot is not free, because keeping the
#vehicle at a constant tilt away from the vertical requires linear acceleration

#how would we do this? set omega1, omega2, omega1_dot, theta as parameters.
#then, we should try to eliminate the four internal forces from the dynamic
#equations. this would leave us with 2 dynamic equations and 1 motor model,
#which contain theta, omega1, omega2, omega1dot, omega2dot, tau_m, u (7
#variables) theta, omega1, omega2, omega1dot are trim parameters omega2dot,
#tau_m and u are unknowns basically the cost function involves the squared error
#from these equations for some omega2dot, tau_m and u guesses with an extra step
#we can substitute the motor model for tau_m, and we have omega2dot and u as
#trim_state

#try Lagrangian dynamics derivation and compare results. it's more useful.

################################################################################
################################ State Space ###################################

#TODO: find the maximum steady-state velocity for abs(u) = 1
#we need to leave some margin for stabilization

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

function get_ẋ_ss(mdl::Model{Vehicle})
    @unpack ω, v, θ, η = mdl.ẋ
    XStateSpace(; ω, v, θ, η)
end

function get_x_ss(mdl::Model{Vehicle})
    @unpack ω, v, θ, η = mdl.x
    XStateSpace(; ω, v, θ, η)
end

function get_u_ss(mdl::Model{Vehicle})
    UStateSpace(mdl.u[])
end

function get_y_ss(mdl::Model{Vehicle})
    @unpack ω, v, θ, η, u_m, τ_m = mdl.y
    YStateSpace(; ω, v, θ, η, u_m, τ_m)
end

assign_x_ss!(mdl::Model{Vehicle}, x_ss::AbstractVector{<:Real}) = assign_x_ss!(mdl, XStateSpace(x_ss))

assign_u_ss!(mdl::Model{Vehicle}, u_ss::AbstractVector{<:Real}) = assign_u_ss!(mdl, UStateSpace(u_ss[1]))

function assign_x_ss!(mdl::Model{Vehicle}, x_ss::XStateSpace)
    @unpack ω, v, θ, η = x_ss
    @pack! mdl.x = ω, v, θ, η
end

assign_u_ss!(mdl::Model{Vehicle}, u_ss::UStateSpace) = (mdl.u[] = u_ss.m)


################################################################################
################################ Linearization #################################

function Linearization.linearize(mdl::Model{Vehicle})

    #define state space system's update function
    f = let mdl = mdl
        function (x, u)
            assign_x_ss!(mdl, x)
            assign_u_ss!(mdl, u)
            f_ode!(mdl)
            get_ẋ_ss(mdl)
        end
    end

    #define state space system's output function
    g = let mdl = mdl
        function (x, u)
            assign_x_ss!(mdl, x)
            assign_u_ss!(mdl, u)
            f_ode!(mdl)
            get_y_ss(mdl)
        end
    end

    #define linearization point
    x0 = get_x_ss(mdl)
    u0 = get_u_ss(mdl)

    linearize(f, g, x0, u0)

end



################################################################################
################################ Linearization #################################

@kwdef struct XController <: FieldVector{3, Float64}
    ω::Float64 = 0.0
    v::Float64 = 0.0
    θ::Float64 = 0.0
end

@kwdef struct UController <: FieldVector{1, Float64}
    m::Float64 = 0.0
end

@kwdef struct ZController <: FieldVector{1, Float64}
    v::Float64 = 0.0
end

struct Controller <: ModelDefinition end



#InvertedPendulum
@kwdef struct Robot{C} <: ModelDefinition
    vehicle::Vehicle = Vehicle()
    controller::C = Controller()
end

end