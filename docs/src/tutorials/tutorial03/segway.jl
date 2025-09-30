using StaticArrays

function num_solve()

    g = 9.81
    m1 = 0.5
    m2 = 0.1
    L = 0.15
    R = 0.05
    k_motor = 0.32
    b_motor = 0.0189

    #approximate main body moment of inertia wrt its CoM as that of a thin rod
    #of length 2L
    @show J1 = 1/12 * m1 * (2L)^2

    @show J_motor = 0.0014
    @show J_wheel = 1/2 * m2 * R^2
    @show J2 = J_motor + J_wheel

    u = 0.118125
    ω1 = 0.0
    ω2 = 2.0
    θ = 0

    #with a nonzero theta, we need to accelerate to keep ω̇1=0
    u = 0.12839
    ω1 = 0.0
    ω2 = 2.0
    θ = 0.01

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
        k_motor * u - b_motor * (ω2 - ω1)
    ]

    ω̇1, ω̇2, F_21x, F_21z, F_i2x, F_i2z, τ_m = A\b

    #kinematics
    θ̇1 = ω1
    η̇ = ω2 * R

    return ω̇1, ω̇2, θ̇1, η̇, F_21x, F_21z, F_i2x, F_i2z, τ_m

end