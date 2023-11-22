module MCSLookups

using Flight
using Flight.FlightCore

using Flight.FlightPhysics
using Flight.FlightComponents

using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172FBWMCS

using UnPack
using ControlSystems
using RobustAndOptimalControl
using StructArrays
using LinearAlgebra


function design_θv_p()

    ac = Cessna172FBWBase(NED()) |> System #linearization requires NED kinematics

    design_point = C172.TrimParameters()

    #complete longitudinal model
    P_lss_lon = Control.Continuous.LinearizedSS(ac, design_point; model = :lon);
    P_nss_lon = named_ss(P_lss_lon);

    x_labels_lon = keys(P_lss_lon.x0) |> collect
    y_labels_lon = keys(P_lss_lon.y0) |> collect
    u_labels_lon = keys(P_lss_lon.u0) |> collect
    z_labels_lon = [:θ, :EAS] #selected command variables

    x_labels = deleteat!(x_labels_lon, findfirst(isequal(:h), x_labels_lon))
    y_labels = deleteat!(y_labels_lon, findfirst(isequal(:h), y_labels_lon))
    u_labels = u_labels_lon
    z_labels = z_labels_lon

    n_x = length(x_labels)
    n_u = length(u_labels)
    n_z = length(z_labels)

    #design model
    P_lss = Control.Continuous.submodel(P_lss_lon; x = x_labels, u = u_labels, y = y_labels)
    P_nss = named_ss(P_lss);

    #before proceeding, ensure consistency in component selection and ordering
    #between our design model and MCS avionics implementation for state, control and
    #command vectors
    @assert tuple(x_labels...) === propertynames(C172FBWMCS.XLon())
    @assert tuple(u_labels...) === propertynames(C172FBWMCS.ULon())
    @assert tuple(z_labels...) === propertynames(C172FBWMCS.ZθEAS())

    #some useful signal labels
    u_labels_fbk = Symbol.(string.(u_labels) .* "_fbk")
    u_labels_fwd = Symbol.(string.(u_labels) .* "_fwd")
    z_labels_dmd = Symbol.(string.(z_labels) .* "_dmd")

    x_trim = P_lss.x0
    u_trim = P_lss.u0
    z_trim = P_lss.y0[z_labels]

    v_norm = norm([x_trim.v_x, x_trim.v_z])

    #weight matrices
    Q = C172FBWMCS.XLon(q = 1, θ = 5, v_x = 1/v_norm, v_z = 1/v_norm, α_filt = 0, ω_eng = 0, thr_v = 0, thr_p = 0, ele_v = 0, ele_p = 0) |> diagm
    R = C172FBWMCS.ULon(throttle_cmd = 1, elevator_cmd = 0.1) |> diagm

    #feedback gain matrix
    C_fbk = lqr(P_nss, Q, R)
    C_fbk_ss = named_ss(ss(C_fbk); u = x_labels, y = u_labels_fbk)

    elevator_sum = sumblock("elevator_cmd = elevator_cmd_fwd - elevator_cmd_fbk")
    throttle_sum = sumblock("throttle_cmd = throttle_cmd_fwd- throttle_cmd_fbk")
    connections = vcat(Pair.(x_labels, x_labels), Pair.(u_labels, u_labels), Pair.(C_fbk_ss.y, C_fbk_ss.y))

    P_nss_fbk = connect([elevator_sum, throttle_sum, P_nss, C_fbk_ss], connections; w1 = u_labels_fwd, z1 = vcat(y_labels, u_labels))

    ############################## Compute feedforward #########################
    F = P_lss.A
    G = P_lss.B
    H_x = P_lss.C[z_labels, :]
    H_u = zeros(n_z, n_u)

    A = [F G; H_x H_u]
    B = inv(A)
    B_12 = B[1:n_x, n_x+1:end]
    B_22 = B[n_x+1:end, n_x+1:end]
    C_fwd = B_22 + C_fbk * B_12
    C_fwd_ss = named_ss(ss(C_fwd), u = z_labels_dmd, y = u_labels_fwd)

    connections = Pair.(C_fwd_ss.y, C_fwd_ss.y)
    P_nss_fbk_θv = connect([C_fwd_ss, P_nss_fbk], connections; w1 = z_labels_dmd, z1 = P_nss_fbk.y)

    return C_fbk

end


end #module