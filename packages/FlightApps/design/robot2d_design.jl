module Robot2DDesign

using FlightCore, FlightPhysics, FlightApps
using FlightPhysics.Control: LQRDataPoint
using FlightApps.Robot2D: Vehicle

using LinearAlgebra
using ControlSystems, RobustAndOptimalControl, ComponentArrays, HDF5

function design_v2m_lqr(mdl::Model{<:Vehicle} = Model(Vehicle()); save = true,
    file::String = joinpath(dirname(@__DIR__), normpath("src/robot2d/robot2d.h5")))

    lss = linearize(mdl)
    lss = delete_vars(lss, :η) #build reduced design model

    x_trim = lss.x0
    n_x = length(x_trim)
    x_labels = collect(keys(x_trim))
    @assert tuple(x_labels...) === propertynames(Robot2D.XController())

    u_labels = [:m, ]
    u_trim = lss.u0[u_labels]

    z_labels = [:v, ]
    z_trim = lss.y0[z_labels]

    A = lss.A
    B = lss.B
    C = lss.C[z_labels, :]
    D = lss.D[z_labels, :]

    C_int = C[z_labels, :]
    D_int = D[z_labels, :]
    n_int = size(C_int, 1)

    A_aug = [A zeros(size(A, 1), n_int); C_int zeros(n_int, n_int)]
    B_aug = [B; D_int]
    C_aug = [C zeros(size(C, 1), n_int)]
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

    #compute augmented feedback matrix
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

    point = LQRDataPoint(; K_fbk, K_fwd, K_int, x_trim, u_trim, z_trim)

    if save
        h5open(file, "w") do fid
            for name in propertynames(point)
                fid[string(name)] = Array(getproperty(point, name))
            end
        end
    end

    return point

end


end #module
