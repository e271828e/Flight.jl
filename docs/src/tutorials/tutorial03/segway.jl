using StaticArrays

function num_solve()

    g = 9.81
    m1 = 0.5
    m2 = 0.1
    L = 0.3
    R = 0.05
    k_motor = 0.32
    b_motor = 0.0189
    J_motor = 0.0014
    @show J_wheel = 1/2 * m2 * R^2
    J1 = 1/12 * m1 * L^2
    J2 = J_motor + J_wheel

    #quiza mejor poner a capon los momentos de inercia y sustituir L/2 por la
    #distancia al cg del rod assembly. Y no lo llamemos rod. Llamemoslo main
    #body. Basta con decir que tiene simetria y cual es su MoI y su L al cg.
    #Rehacer eqs sustituyendo L/2 por una LCoM

    u = 0.118125
    ω1 = 0.0
    ω2 = 2.0
    θ = 0

    #with a nonzero theta, we need to accelerate to keep ω̇1=0
    # u = 0.12839
    # ω1 = 0.0
    # ω2 = 2.0
    # θ = 0.01

    A = @SMatrix[
        m1*L/2*cos(θ)  m1*R    -1            0              0     0     0
        m1*L/2*sin(θ)  0       0             1              0     0     0
        J1             0       L/2*cos(θ)   -L/2*sin(θ)     0     0     1
        0              m2*R    1             0              -1    0     0
        0              0       0             -1             0     1     0
        0              J2      0             0              R     0     -1
        0              0       0             0              0     0     1
    ]

    b = @SVector[
        m1*L/2*ω1^2*sin(θ),
        m1*g - m1*L/2*ω1^2*cos(θ),
        0,
        0,
        m2*g,
        0,
        k_motor * u - b_motor * (ω2 - ω1)
    ]

    ω̇1, ω̇2, F_21x, F_21z, F_i2x, F_i2z, τ_m = A\b

    #kinematics
    θ̇1 = ω1
    ṗ = ω2 * R

    return ω̇1, ω̇2, θ̇1, ṗ, F_21x, F_21z, F_i2x, F_i2z, τ_m

end