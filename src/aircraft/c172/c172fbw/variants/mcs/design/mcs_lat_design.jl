module MCSLatDesign

using Flight
using Flight.FlightCore

using Flight.FlightPhysics
using Flight.FlightComponents
using Flight.FlightComponents.Control.Continuous: LinearizedSS
using Flight.FlightComponents.Control.Discrete: Integrator, IntegratorOutput, PID, PIDOutput, LQRTracker, LQRTrackerOutput

using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172FBWMCS

using UnPack
using ControlSystems
using RobustAndOptimalControl
using StructArrays
using ComponentArrays
using LinearAlgebra


function design_model(design_point::C172.TrimParameters = C172.TrimParameters())

    ac = Cessna172FBWBase(NED()) |> System #linearization requires NED kinematics

    P_lss_lat = Control.Continuous.LinearizedSS(ac, design_point; model = :lat);

    x_labels_lat = keys(P_lss_lat.x0) |> collect
    y_labels_lat = keys(P_lss_lat.y0) |> collect
    u_labels_lat = keys(P_lss_lat.u0) |> collect

    #remove heading from state vector to avoid quasi-static equilibrium
    x_labels = deleteat!(x_labels_lat, findfirst(isequal(:ψ), x_labels_lat))
    y_labels = deleteat!(y_labels_lat, findfirst(isequal(:ψ), y_labels_lat))
    u_labels = u_labels_lat

    #ensure consistency in component selection and ordering between our design model
    #and MCS avionics implementation for state and control vectors
    @assert tuple(x_labels...) === propertynames(C172FBWMCS.XLat())
    @assert tuple(u_labels...) === propertynames(C172FBWMCS.ULat())

    #extract design model
    P_lss = Control.Continuous.submodel(P_lss_lat; x = x_labels, u = u_labels, y = y_labels)

    return P_lss

end


function design_φβ(P_lss::LinearizedSS = design_model())

    P_nss = named_ss(P_lss)

    @show x_labels = P_nss.x
    @show u_labels = P_nss.u
    @show z_labels = [:φ, :β]
    @assert tuple(z_labels...) === propertynames(C172FBWMCS.ZLatPhiBeta())

    x_trim = P_lss.x0
    u_trim = P_lss.u0
    z_trim = P_lss.y0[z_labels]

    n_x = length(P_lss.x0)
    n_u = length(P_lss.u0)
    n_z = length(z_labels)

    F = P_lss.A
    G = P_lss.B
    Hx = P_lss.C[z_labels, :]
    Hu = P_lss.D[z_labels, :]

    Hx_int = Hx[:β, :]'
    Hu_int = Hu[:β, :]'
    n_int, _ = size(Hx_int)

    F_aug = [F zeros(n_x, n_int); Hx_int zeros(n_int, n_int)]
    G_aug = [G; Hu_int]
    Hx_aug = [Hx zeros(n_z, n_int)]
    Hu_aug = Hu

    P_aug = ss(F_aug, G_aug, Hx_aug, Hu_aug)

    @unpack v_x, v_y = P_lss.x0
    v_norm = norm([v_x, v_y])

    #weight matrices
    Q = ComponentVector(p = 0, r = 0.1, φ = 0.15, v_x = 0/v_norm, v_y = 0.1/v_norm, β_filt = 0, ail_v = 0, ail_p = 0, rud_v = 0, rud_p = 0, ξ_β = 0.001) |> diagm
    R = C172FBWMCS.ULat(aileron_cmd = 0.1, rudder_cmd = 0.05) |> diagm

    #compute gain matrix
    C_aug = lqr(P_aug, Q, R)

    #extract system state and integrator blocks from the feedback matrix
    C_x = C_aug[:, 1:n_x]
    C_ξ = C_aug[:, n_x+1:end]

    #construct
    A = [F G; Hx Hu]
    B = inv(A)
    B_12 = B[1:n_x, n_x+1:end]
    B_22 = B[n_x+1:end, n_x+1:end]

    C_fbk = C_x
    C_fwd = B_22 + C_x * B_12
    C_int = ComponentMatrix(zeros(n_u, n_z), Axis(u_labels), Axis(z_labels))
    C_int[:, :β] .= C_ξ

    params_φ_β = Control.Discrete.LQRTrackerParams(; C_fbk, C_fwd, C_int, x_trim, u_trim, z_trim)

    #some useful signal labels
    u_labels_fbk = Symbol.(string.(u_labels) .* "_fbk")
    u_labels_fwd = Symbol.(string.(u_labels) .* "_fwd")
    u_labels_sum = Symbol.(string.(u_labels) .* "_sum")
    u_labels_int_u = Symbol.(string.(u_labels) .* "_int_u")
    u_labels_int = Symbol.(string.(u_labels) .* "_int")
    u_labels_ξ = Symbol.(string.(u_labels) .* "_ξ")

    z_labels_sp = Symbol.(string.(z_labels) .* "_sp")
    z_labels_sp1 = Symbol.(string.(z_labels) .* "_sp1")
    z_labels_sp2 = Symbol.(string.(z_labels) .* "_sp2")
    z_labels_err = Symbol.(string.(z_labels) .* "_err")
    z_labels_sum = Symbol.(string.(z_labels) .* "_sum")
    z_labels_sp_fwd = Symbol.(string.(z_labels) .* "_sp_fwd")
    z_labels_sp_sum = Symbol.(string.(z_labels) .* "_sp_sum")


    C_fbk_ss = named_ss(ss(C_fbk), u = x_labels, y = u_labels_fbk)

    C_fwd_ss = named_ss(ss(C_fwd), u = z_labels_sp_fwd, y = u_labels_fwd)

    C_int_ss = named_ss(ss(C_int), u = z_labels_err, y = u_labels_int_u)

    int_ss = named_ss(ss(tf(1, [1,0])) .* I(2),
                        x = u_labels_ξ,
                        u = u_labels_int_u,
                        y = u_labels_int);

    P_nss = named_ss(P_lss);

    φ_err_sum = sumblock("φ_err = φ_sum - φ_sp_sum")
    β_err_sum = sumblock("β_err = β_sum - β_sp_sum")

    aileron_cmd_sum = sumblock("aileron_cmd_sum = aileron_cmd_fwd - aileron_cmd_fbk - aileron_cmd_int")
    rudder_cmd_sum = sumblock("rudder_cmd_sum = rudder_cmd_fwd - rudder_cmd_fbk - rudder_cmd_int")

    φ_sp_splitter = splitter(:φ_sp, 2)
    β_sp_splitter = splitter(:β_sp, 2)

    connections = vcat(
        Pair.(x_labels, x_labels),
        Pair.(z_labels, z_labels_sum),
        Pair.(z_labels_sp1, z_labels_sp_sum),
        Pair.(z_labels_sp2, z_labels_sp_fwd),
        Pair.(z_labels_err, z_labels_err),
        Pair.(u_labels_sum, u_labels),
        Pair.(u_labels_fwd, u_labels_fwd),
        Pair.(u_labels_fbk, u_labels_fbk),
        Pair.(u_labels_int, u_labels_int),
        Pair.(u_labels_int_u, u_labels_int_u),
        )

    P_nss_φβ = connect([P_nss, int_ss, C_fwd_ss, C_fbk_ss, C_int_ss,
                        φ_err_sum, β_err_sum,
                        aileron_cmd_sum, rudder_cmd_sum,
                        φ_sp_splitter, β_sp_splitter], connections;
                        w1 = z_labels_sp, z1 = vcat(y_labels))


    P_φ2p = P_nss_φβ[:p, :φ_sp];


end


end #module