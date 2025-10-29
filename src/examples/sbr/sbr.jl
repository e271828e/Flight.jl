module SelfBalancingRobot

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore
using Flight.FlightLib

const g = 9.80665 #m/s^2, standard gravity

################################################################################
################################## Vehicle #####################################

#body 1: main body, comprising the vehicle chassis and the DC motor's case and stator
#body 2: rolling body, comprising the wheel, axle, and the DC motor's rotor.

@kwdef struct Vehicle <: ModelDefinition
    L::Float64 = 0.15 #distance from main body's origin to its CoM (m)
    R::Float64 = 0.05 #wheel radius (m)
    m1::Float64 = 0.5 #mass of main body (kg)
    m2::Float64 = 0.1 #mass of rolling body (kg)
    k_m::Float64 = 0.32 #motor's torque constant (N*m)
    b_m::Float64 = 0.0189 #motor's effective damping coefficient (N*m*s/rad)
    J_m::Float64 = 0.0014 #motor's effective moment of inertia (kg*m^2)
end

@kwdef struct VehicleY
    ω1::Float64 = 0.0 #angular velocity of main body (rad/s)
    ω2::Float64 = 0.0 #angular velocity of rolling body (rad/s)
    θ::Float64 = 0.0 #angle of main body with respect to vertical (rad)
    η::Float64 = 0.0 #horizontal position (m)
    v::Float64 = 0.0 #linear velocity (m/s)
    ω1_dot::Float64 = 0.0 #angular acceleration of main body (rad/s^2)
    ω2_dot::Float64 = 0.0 #angular acceleration of rolling body (rad/s^2)
    F_21x::Float64 = 0.0 #force exerted by body 2 on body 1 along inertial frame's x-axis (N)
    F_21z::Float64 = 0.0 #force exerted by body 2 on body 1 along inertial frame's z-axis (N)
    F_i2x::Float64 = 0.0 #force exerted by inertial frame on body 2 along inertial frame's x-axis (N)
    F_i2z::Float64 = 0.0 #force exerted by inertial frame on body 2 along inertial frame's z-axis (N)
    τ_m::Float64 = 0.0 #motor torque exerted by body 1 on body 2 (N*m)
end

Modeling.X(::Vehicle) = ComponentVector( ω1 = 0.0, ω2 = 0.0, θ = 0.0, η = 0.0)
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
    @unpack L, R, k_m, b_m, J_m, m1, m2 = parameters

    @unpack ω1, ω2, θ, η = x
    u_m = u[]

    #approximate main body's moment of inertia with respect to its CoM as
    #that of a thin rod of length 2L
    J1 = 1/12 * m1 * (2L)^2

    #approximate rolling body's moment of inertia with respect to its CoM as the
    #motor's effective moment of inertia plus the contribution from the body's
    #mass, assuming a solid disk distribution
    J2 = J_m + 1/2 * m2 * R^2

    A = @SMatrix[
        m1*L*cos(θ)  m1*R    -1         0           0     0     0
        m1*L*sin(θ)  0       0          1           0     0     0
        J1           0       L*cos(θ)   -L*sin(θ)   0     0     1
        0            m2*R    1          0           -1    0     0
        0            0       0          -1          0     1     0
        0            J2      0          0           R     0     -1
        0            0       0          0           0     0     1
    ]

    b = @SVector[
        m1*L*ω1^2*sin(θ),
        m1*g - m1*L*ω1^2*cos(θ),
        0,
        0,
        m2*g,
        0,
        k_m * u_m - b_m * (ω2 - ω1)
    ]

    ω1_dot, ω2_dot, F_21x, F_21z, F_i2x, F_i2z, τ_m = A\b
    v = ω2 * R

    θ_dot = ω1
    η_dot = v

    ẋ.ω1 = ω1_dot
    ẋ.ω2 = ω2_dot
    ẋ.θ = θ_dot
    ẋ.η = η_dot

    mdl.y = VehicleY(; ω1, ω2, θ, η, v, ω1_dot, ω2_dot, F_21x, F_21z, F_i2x, F_i2z, τ_m)

end

@no_step Vehicle
@no_periodic Vehicle


################################################################################
################################ State Space ###################################

#TODO: find the maximum steady-state velocity for abs(u) = 1
#TODO: keep τ_m in YStateSpace to check for saturation during controller design,
#we need to leave some margin for stabilization

@kwdef struct XStateSpace <: FieldVector{4, Float64}
    ω1::Float64 = 0.0
    ω2::Float64 = 0.0
    θ::Float64 = 0.0
    η::Float64 = 0.0
end

@kwdef struct UStateSpace <: FieldVector{1, Float64}
    motor::Float64 = 0.0
end

@kwdef struct YStateSpace <: FieldVector{6, Float64}
    ω1::Float64 = 0.0
    ω2::Float64 = 0.0
    θ::Float64 = 0.0
    η::Float64 = 0.0
    v::Float64 = 0.0
    τ_m::Float64 = 0.0
end

function get_ẋ_ss(mdl::Model{Vehicle})
    @unpack ω1, ω2, θ, η = mdl.ẋ
    XStateSpace(; ω1, ω2, θ, η)
end

function get_x_ss(mdl::Model{Vehicle})
    @unpack ω1, ω2, θ, η = mdl.x
    XStateSpace(; ω1, ω2, θ, η)
end

function get_u_ss(mdl::Model{Vehicle})
    UStateSpace(mdl.u[])
end

function get_y_ss(mdl::Model{Vehicle})
    @unpack ω1, ω2, θ, η, v, τ_m = mdl.y
    YStateSpace(; ω1, ω2, θ, η, v, τ_m)
end

assign_x_ss!(mdl::Model{Vehicle}, x::AbstractVector{<:Real}) = assign_x_ss!(mdl, XStateSpace(x))

function assign_x_ss!(mdl::Model{Vehicle}, x::XStateSpace)
    @unpack ω1, ω2, θ, η = x
    @pack! mdl.x = ω1, ω2, θ, η
end

assign_u_ss!(mdl::Model{Vehicle}, u::AbstractVector{<:Real}) = assign_u_ss!(mdl, UStateSpace(u[1]))

assign_u_ss!(mdl::Model{Vehicle}, u::UStateSpace) = (mdl.u[] = u.motor)


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


function test_vehicle()

    mdl = Vehicle() |> Model

    mdl.u[] = 0
    mdl.x.ω1 = 0.0
    mdl.x.ω2 = 0
    mdl.x.θ = 0.0
    f_ode!(mdl)
    @show mdl.ẋ
    @show mdl.y

    mdl.u[] = 0.118125
    mdl.x.ω1 = 0.0
    mdl.x.ω2 = 2.0
    mdl.x.θ = 0.0
    f_ode!(mdl)
    @show mdl.ẋ
    @show mdl.y

    mdl.u[] = 1.0
    mdl.x.ω1 = 0.0
    mdl.x.ω2 = 17.0
    mdl.x.θ = 0.0
    f_ode!(mdl)
    @show mdl.ẋ
    @show mdl.y

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

struct Controller <: ModelDefinition end

#InvertedPendulum
@kwdef struct Robot{C} <: ModelDefinition
    vehicle::Vehicle = Vehicle()
    controller::C = Controller()
end



end #module