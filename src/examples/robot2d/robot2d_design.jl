module Robot2DControlDesign

using Flight.FlightCore
using Flight.FlightLib

using Flight.FlightLib.Linearization: delete_vars
using Flight.FlightLib.Control.Discrete: LQRParams
using ControlSystems, RobustAndOptimalControl, Plots, UnPack, ComponentArrays, LinearAlgebra

include("robot2d.jl")
using .Robot2D
using .Robot2D: Vehicle

function design_velocity_tracker(vehicle::Vehicle = Vehicle())

    mdl = Model(vehicle)
    lss = linearize(mdl)
    lss_red = delete_vars(lss, :η)
    P = named_ss(lss)

    P_v, params_v2m = let lss = lss_red

        x_trim = lss.x0
        n_x = length(x_trim)
        x_labels = collect(keys(x_trim))
        @assert tuple(x_labels...) === propertynames(Robot2D.XController())

        u_trim = lss.u0
        n_u = length(u_trim)
        u_labels = collect(keys(u_trim))
        @assert tuple(u_labels...) === propertynames(Robot2D.UController())

        z_labels = [:v, ]
        z_trim = lss.y0[z_labels]
        n_z = length(z_labels)
        @assert tuple(z_labels...) === propertynames(Robot2D.ZController())

        A = lss.A
        B = lss.B
        C = lss.C[z_labels, :]
        D = lss.D[z_labels, :]

        C_int = C[z_labels, :]
        D_int = D[z_labels, :]
        n_int, _ = size(C_int)

        A_aug = [A zeros(n_x, n_int); C_int zeros(n_int, n_int)]
        B_aug = [B; D_int]
        C_aug = [C zeros(n_z, n_int)]
        D_aug = D

        P_aug = ss(A_aug, B_aug, C_aug, D_aug)

        x_aug_labels = push!(copy(x_labels), :ξ_v)
        Q_diag = ComponentVector(zeros(length(x_aug_labels)), Axis(x_aug_labels))
        R_diag = ComponentVector(zeros(length(u_labels)), Axis(u_labels))

        Q_diag.ω = 1e-3
        Q_diag.v = 1e-2
        Q_diag.θ = 0
        Q_diag.ξ_v = 5e-2

        R_diag.m = 1e-1

        Q = diagm(Q_diag)
        R = diagm(R_diag)

        #compute gain matrix
        K_aug = lqr(P_aug, Q, R)

        L = [A B; C D]
        M = inv(L)
        M_12 = M[1:n_x, n_x+1:end]
        M_22 = M[n_x+1:end, n_x+1:end]

        #extract system state and integrator blocks from the feedback matrix
        K_x = K_aug[:, 1:n_x]
        K_ξ = K_aug[:, n_x+1:end]

        K_fbk = K_x
        K_fwd = M_22 + K_x * M_12
        K_int = K_ξ

        K_fbk_ss = named_ss(ss(K_fbk), u = x_labels, y = :m_fbk)
        K_fwd_ss = named_ss(ss(K_fwd), u = :v_ref, y = :m_fwd)
        K_int_ss = named_ss(ss(K_int), u = :v_err, y = :m_int_in)

        int_ss = named_ss(ss(tf(1, [1,0])) .* I(1),
                            x = :m_ξ,
                            u = :m_int_in,
                            y = :m_int_out);

        v_sum = sumblock("v_err = v - v_ref")
        m_sum = sumblock("m = m_fwd - m_fbk - m_int_out")

        connections = vcat(
            Pair.(x_labels, x_labels),
            :v_err => :v_err,
            :m_int_in=> :m_int_in,
            :m_int_out => :m_int_out,
            :m_fwd => :m_fwd,
            :m_fbk => :m_fbk,
            :m => :m,
            )

            #connect back to full plant
        P_v = connect([P, int_ss, K_fwd_ss, K_fbk_ss, K_int_ss,
                        v_sum, m_sum], connections;
                        w1 = :v_ref, z1 = P.y, unique = false)

        params = LQRParams(;
            K_fbk = Matrix(K_fbk), K_fwd = Matrix(K_fwd), K_int = Matrix(K_int),
            x_trim = Vector(x_trim), u_trim = Vector(u_trim), z_trim = Vector(z_trim))

        (P_v, params)

    end

    P_η, params_η2v = let

    end

    error("Add η PID output bounds and save together with LQR and PID design parameters")

end

end #module
