module MCSDesign

using Flight
using Flight.FlightCore

using Flight.FlightPhysics
using Flight.FlightComponents
using Flight.FlightComponents.Control.Continuous: LinearizedSS
using Flight.FlightComponents.Control.Discrete: PIDParams, LQRTrackerParams
using Flight.FlightComponents.Control.PIDOpt: Settings, Metrics, optimize_PID, build_PID, check_results

using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172MCS

using HDF5
using Logging
using UnPack
using ControlSystems
using RobustAndOptimalControl
using StaticArrays
using StructArrays
using ComponentArrays
using LinearAlgebra
using Interpolations


function generate_lookups(
    EAS_range::AbstractRange{Float64} = range(25, 55, length = 7),
    h_range::AbstractRange{Float64} = range(50, 3050, length = 4);
    channel::Symbol = :lat,
    folder::String = dirname(@__DIR__)) #save to parent folder by default

    if channel === :lon
        f_opt = design_lon
    elseif channel === :lat
        f_opt = design_lat
    else
        error("Valid values for channel keyword: :lon, :lat")
    end

    results = map(Iterators.product(EAS_range, h_range)) do (EAS, h)

        println("Designing $channel controllers for EAS = $EAS, h = $h")

        #all other design point parameters at default
        flaps = EAS < 30 ? 1.0 : 0.0
        design_point = C172.TrimParameters(; Ob = Geographic(LatLon(), HEllip(h)), EAS, flaps)

        results = f_opt(design_point)
        return results

    end |> StructArray |> StructArrays.components

    filenames = joinpath.(dirname(@__DIR__), "data", string.(keys(results)) .* "_lookup.h5")

    bounds = ((EAS_range[1], EAS_range[end]), (h_range[1], h_range[end]))

    foreach(values(results), filenames) do data, fname
        C172MCS.save_lookup(data, bounds, joinpath(folder, fname))
    end

    return results

end

function design_lon(design_point::C172.TrimParameters = C172.TrimParameters())

    ac = Cessna172FBWBase(NED()) |> System #linearization requires NED kinematics

    lss_lon = Control.Continuous.LinearizedSS(ac, design_point; model = :lon);

    x_labels_lon = keys(lss_lon.x0) |> collect
    y_labels_lon = keys(lss_lon.y0) |> collect
    u_labels_lon = keys(lss_lon.u0) |> collect

    x_labels = copy(x_labels_lon)
    y_labels = copy(y_labels_lon)
    u_labels = copy(u_labels_lon)

    x_labels = deleteat!(x_labels, findfirst(isequal(:h), x_labels))
    y_labels = deleteat!(y_labels, findfirst(isequal(:h), y_labels))
    u_labels = u_labels_lon

    #ensure consistency in component selection and ordering between our design model
    #and MCS avionics implementation for state and control vectors
    @assert tuple(x_labels...) === propertynames(C172MCS.XLon())
    @assert tuple(u_labels...) === propertynames(C172MCS.ULon())

    lss_red = Control.Continuous.submodel(lss_lon; x = x_labels, u = u_labels, y = y_labels)
    x_trim = lss_red.x0
    u_trim = lss_red.u0

    n_x = length(lss_red.x0)
    n_u = length(lss_red.u0)

    P_lon = named_ss(lss_lon)
    P_red = named_ss(lss_red);

    ############################ thr+ele SAS ###################################

    P_te, params_te2te = let

        z_labels = [:throttle_cmd, :elevator_cmd]

        @assert tuple(z_labels...) === propertynames(C172MCS.ZLonThrEle())
        z_trim = lss_red.y0[z_labels]
        n_z = length(z_labels)

        F = lss_red.A
        G = lss_red.B
        Hx = lss_red.C[z_labels, :]
        Hu = lss_red.D[z_labels, :]

        @unpack v_x, v_z = lss_red.x0
        v_norm = norm([v_x, v_z])

        #weight matrices
        Q = C172MCS.XLon(q = 1, θ = 30, v_x = 0.01/v_norm, v_z = 1/v_norm, α_filt = 0,
            ω_eng = 0, thr_v = 0.0, thr_p = 0, ele_v = 0, ele_p = 0) |> diagm
        R = C172MCS.ULon(throttle_cmd = 2, elevator_cmd = 2) |> diagm

        #feedback gain matrix
        C_fbk = lqr(P_red, Q, R)

        #forward gain matrix
        A = [F G; Hx Hu]
        B = inv(A)
        B_12 = B[1:n_x, n_x+1:end]
        B_22 = B[n_x+1:end, n_x+1:end]
        C_fwd = B_22 + C_fbk * B_12

        #no integral control
        C_int = zeros(n_u, n_z)

        #useful signal labels for connections
        u_labels_fbk = Symbol.(string.(u_labels) .* "_fbk")
        u_labels_fwd = Symbol.(string.(u_labels) .* "_fwd")
        u_labels_sum = Symbol.(string.(u_labels) .* "_sum")

        #some useful signal labels
        u_labels_fbk = Symbol.(string.(u_labels) .* "_fbk")
        u_labels_fwd = Symbol.(string.(u_labels) .* "_fwd")
        u_labels_sum = Symbol.(string.(u_labels) .* "_sum")
        z_labels_sp = Symbol.(string.(z_labels) .* "_sp")

        C_fbk_ss = named_ss(ss(C_fbk), u = x_labels, y = u_labels_fbk)
        C_fwd_ss = named_ss(ss(C_fwd), u = z_labels_sp, y = u_labels_fwd)

        #summing junctions
        throttle_sum = sumblock("throttle_cmd_sum = throttle_cmd_fwd- throttle_cmd_fbk")
        elevator_sum = sumblock("elevator_cmd_sum = elevator_cmd_fwd - elevator_cmd_fbk")

        connections = vcat(
            Pair.(x_labels, x_labels),
            Pair.(u_labels_fbk, u_labels_fbk),
            Pair.(u_labels_sum, u_labels),
            Pair.(u_labels_fwd, u_labels_fwd)
            )

        P_te = connect([P_lon, throttle_sum, elevator_sum, C_fbk_ss, C_fwd_ss],
            connections; w1 = z_labels_sp, z1 = P_lon.y)

        params_te = LQRTrackerParams(; #export everything as plain arrays
            C_fbk = Matrix(C_fbk), C_fwd = Matrix(C_fwd), C_int = Matrix(C_int),
            x_trim = Vector(x_trim), u_trim = Vector(u_trim), z_trim = Vector(z_trim))

        (P_te, params_te)

    end

    P_tq, params_q2e = let

        P_e2q = P_te[:q, :elevator_cmd_sp]

        q2e_int = tf(1, [1, 0]) |> ss
        P_q2e_opt = series(q2e_int, ss(P_e2q))

        t_sim_q2e = 10
        lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDParams(; k_p = 10.0, k_i = 35.0, k_d = 1.5, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_q2e, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 2, ∫e = 10, ef = 2, ∫u = 0.1, up = 0.00)
        params_0 = PIDParams(; k_p = 1.5, k_i = 10, k_d = 0.2, τ_f = 0.01)

        q2e_results = optimize_PID(P_q2e_opt; params_0, settings, weights, global_search = false)

        params_q2e = q2e_results.params
        if !check_results(q2e_results, Metrics(; Ms = 1.3, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
            println("Warning: Checks failed for pitch rate control PID, design point $(design_point)")
            println(q2e_results.metrics)
        end

        q2e_pid = build_PID(q2e_results.params)
        C_q2e = named_ss(series(q2e_int, q2e_pid), :C_q2e; u = :q_err, y = :elevator_cmd_sp);

        q2e_sum = sumblock("q_err = q_sp - q")
        P_tq = connect([P_te, q2e_sum, C_q2e],
            [:q_err=>:q_err, :q=>:q, :elevator_cmd_sp=>:elevator_cmd_sp],
            w1 = [:throttle_cmd_sp, :q_sp], z1 = P_te.y)

        (P_tq, params_q2e)

    end

    P_vq, params_v2t = let

        P_t2v = P_tq[:EAS, :throttle_cmd]

        t_sim_v2t = 10
        lower_bounds = PIDParams(; k_p = 0.01, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDParams(; k_p = 1.0, k_i = 2.0, k_d = 0.0, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_v2t, maxeval = 5000, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 5.0, ∫e = 10.0, ef = 1.0, ∫u = 0.0, up = 0.0)
        params_0 = PIDParams(; k_p = 0.2, k_i = 0.1, k_d = 0.0, τ_f = 0.01)

        v2t_results = optimize_PID(P_t2v; params_0, settings, weights, global_search = false)

        params_v2t = v2t_results.params
        if !check_results(v2t_results, Metrics(; Ms = 1.3, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
            println("Warning: Checks failed for pitch rate control PID, design point $(design_point)")
            println(v2t_results.metrics)
        end

        v2t_pid = build_PID(v2t_results.params)
        C_v2t = named_ss(ss(v2t_pid), :C_v2t; u = :EAS_err, y = :throttle_cmd_sp)

        v2t_sum = sumblock("EAS_err = EAS_sp - EAS")
        P_vq = connect([P_tq, v2t_sum, C_v2t],
            [:EAS_err=>:EAS_err, :EAS=>:EAS, :throttle_cmd_sp=>:throttle_cmd_sp],
            w1 = [:EAS_sp, :q_sp], z1 = P_tq.y)

        (P_vq, params_v2t)

    end

    P_vc, params_vc2te = let

        z_labels = [:EAS, :climb_rate]

        @assert tuple(z_labels...) === propertynames(C172MCS.ZLonEASClm())
        z_trim = lss_red.y0[z_labels]
        n_z = length(z_labels)

        F = lss_red.A
        G = lss_red.B
        Hx = lss_red.C[z_labels, :]
        Hu = lss_red.D[z_labels, :]

        Hx_int = Hx[z_labels, :]
        Hu_int = Hu[z_labels, :]
        n_int, _ = size(Hx_int)

        F_aug = [F zeros(n_x, n_int); Hx_int zeros(n_int, n_int)]
        G_aug = [G; Hu_int]
        Hx_aug = [Hx zeros(n_z, n_int)]
        Hu_aug = Hu

        P_aug = ss(F_aug, G_aug, Hx_aug, Hu_aug)

        @unpack v_x, v_z = lss_red.x0
        v_norm = norm([v_x, v_z])

        #weight matrices
        Q = ComponentVector(q = 1, θ = 100, v_x = 10/v_norm, v_z = 1/v_norm, α_filt = 1, ω_eng = 0,
            thr_v = 0.0, thr_p = 0, ele_v = 0, ele_p = 0, ξ_EAS = 0.005, ξ_climb_rate = 0.001) |> diagm
        R = C172MCS.ULon(throttle_cmd = 1, elevator_cmd = 1) |> diagm

        C_aug = lqr(P_aug, Q, R)

        #extract system state and integrator blocks from the feedback matrix
        C_x = C_aug[:, 1:n_x]
        C_ξ = C_aug[:, n_x+1:end]

        #construct feedforward matrix blocks
        A = [F G; Hx Hu]
        B = inv(A)
        B_12 = B[1:n_x, n_x+1:end]
        B_22 = B[n_x+1:end, n_x+1:end]

        C_fbk = C_x
        C_fwd = B_22 + C_x * B_12
        C_int = C_ξ

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

        EAS_err_sum = sumblock("EAS_err = EAS_sum - EAS_sp_sum")
        climb_rate_err_sum = sumblock("climb_rate_err = climb_rate_sum - climb_rate_sp_sum")

        throttle_cmd_sum = sumblock("throttle_cmd_sum = throttle_cmd_fwd - throttle_cmd_fbk - throttle_cmd_int")
        elevator_cmd_sum = sumblock("elevator_cmd_sum = elevator_cmd_fwd - elevator_cmd_fbk - elevator_cmd_int")

        EAS_sp_splitter = splitter(:EAS_sp, 2)
        climb_rate_sp_splitter = splitter(:climb_rate_sp, 2)

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

        #disable warning about connecting single output to multiple inputs (here,
        #φ goes both to state feedback and command variable error junction)
        P_vc = connect([P_lon, int_ss, C_fwd_ss, C_fbk_ss, C_int_ss,
                            EAS_err_sum, climb_rate_err_sum,
                            throttle_cmd_sum, elevator_cmd_sum,
                            EAS_sp_splitter, climb_rate_sp_splitter], connections;
                            w1 = z_labels_sp, z1 = P_lon.y)

        #convert everything to plain arrays
        params_vc2te = LQRTrackerParams(;
            C_fbk = Matrix(C_fbk), C_fwd = Matrix(C_fwd), C_int = Matrix(C_int),
            x_trim = Vector(x_trim), u_trim = Vector(u_trim), z_trim = Vector(z_trim))

        (P_vc, params_vc2te)
    end

    P_vt, params_vt2te = let

        z_labels = [:EAS, :throttle_cmd]

        @assert tuple(z_labels...) === propertynames(C172MCS.ZLonEASThr())
        z_trim = lss_red.y0[z_labels]
        n_z = length(z_labels)

        F = lss_red.A
        G = lss_red.B
        Hx = lss_red.C[z_labels, :]
        Hu = lss_red.D[z_labels, :]

        #define the blocks corresponding to the subset of the command variables for
        #which integral compensation is required. in this case, only one of them (EAS).
        #Since the resulting blocks have only one row, we get vectors, which we need to
        #transpose to get the desired row matrices back.
        Hx_int = Hx[:EAS, :]'
        Hu_int = Hu[:EAS, :]'
        n_int, _ = size(Hx_int)

        F_aug = [F zeros(n_x, n_int); Hx_int zeros(n_int, n_int)]
        G_aug = [G; Hu_int]
        Hx_aug = [Hx zeros(n_z, n_int)]
        Hu_aug = Hu

        P_aug = ss(F_aug, G_aug, Hx_aug, Hu_aug)

        @unpack v_x, v_z = lss_red.x0
        v_norm = norm([v_x, v_z])

        #weight matrices
        Q = ComponentVector(q = 1, θ = 100, v_x = 10/v_norm, v_z = 1/v_norm, α_filt = 1, ω_eng = 0, thr_v = 0.0, thr_p = 0, ele_v = 0, ele_p = 0, ξ_EAS = 0.005) |> diagm
        R = C172MCS.ULon(throttle_cmd = 1, elevator_cmd = 1) |> diagm

        #compute gain matrix
        C_aug = lqr(P_aug, Q, R)

        #extract system state and integrator blocks from the feedback matrix
        C_x = C_aug[:, 1:n_x]
        C_ξ = C_aug[:, n_x+1:end]

        #construct feedforward matrix blocks
        A = [F G; Hx Hu]
        B = inv(A)
        B_12 = B[1:n_x, n_x+1:end]
        B_22 = B[n_x+1:end, n_x+1:end]

        C_fbk = C_x
        C_fwd = B_22 + C_x * B_12
        C_int = ComponentMatrix(zeros(n_u, n_z), Axis(u_labels), Axis(z_labels))
        C_int[:, :EAS] .= C_ξ

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

        EAS_err_sum = sumblock("EAS_err = EAS_sum - EAS_sp_sum")
        throttle_cmd_err_sum = sumblock("throttle_cmd_err = throttle_cmd_sum - throttle_cmd_sp_sum")

        throttle_cmd_sum = sumblock("throttle_cmd_sum = throttle_cmd_fwd - throttle_cmd_fbk - throttle_cmd_int")
        elevator_cmd_sum = sumblock("elevator_cmd_sum = elevator_cmd_fwd - elevator_cmd_fbk - elevator_cmd_int")

        EAS_sp_splitter = splitter(:EAS_sp, 2)
        throttle_cmd_sp_splitter = splitter(:throttle_cmd_sp, 2)

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

        P_vt = connect([P_lon, int_ss, C_fwd_ss, C_fbk_ss, C_int_ss,
                        EAS_err_sum, throttle_cmd_err_sum,
                        throttle_cmd_sum, elevator_cmd_sum,
                        EAS_sp_splitter, throttle_cmd_sp_splitter], connections;
                        w1 = z_labels_sp, z1 = P_lon.y)

        #convert everything to plain arrays
        params_vt2te = LQRTrackerParams(;
            C_fbk = Matrix(C_fbk), C_fwd = Matrix(C_fwd), C_int = Matrix(C_int),
            x_trim = Vector(x_trim), u_trim = Vector(u_trim), z_trim = Vector(z_trim))

        (P_vt, params_vt2te)
    end

    return (te2te = params_te2te, q2e = params_q2e, v2t = params_v2t,
            vc2te = params_vc2te, vt2te = params_vt2te)

end


function design_lat(design_point::C172.TrimParameters = C172.TrimParameters())

    ac = Cessna172FBWBase(NED()) |> System #linearization requires NED kinematics

    lss_lat = Control.Continuous.LinearizedSS(ac, design_point; model = :lat);

    x_labels_lat = keys(lss_lat.x0) |> collect
    y_labels_lat = keys(lss_lat.y0) |> collect
    u_labels_lat = keys(lss_lat.u0) |> collect

    x_labels = copy(x_labels_lat)
    y_labels = copy(y_labels_lat)
    u_labels = copy(u_labels_lat)

    x_labels = deleteat!(x_labels, findfirst(isequal(:ψ), x_labels))
    y_labels = deleteat!(y_labels, findfirst(isequal(:ψ), y_labels))
    y_labels = deleteat!(y_labels, findfirst(isequal(:χ), y_labels))

    #ensure consistency in component selection and ordering between our design model
    #and MCS avionics implementation for state and control vectors
    @assert tuple(x_labels...) === propertynames(C172MCS.XLat())
    @assert tuple(u_labels...) === propertynames(C172MCS.ULat())

    #extract design model
    lss_red = Control.Continuous.submodel(lss_lat; x = x_labels, u = u_labels, y = y_labels)
    x_trim = lss_red.x0
    u_trim = lss_red.u0

    n_x = length(lss_red.x0)
    n_u = length(lss_red.u0)

    P_lat = named_ss(lss_lat)
    P_red = named_ss(lss_red);

    ############################### φ + β ######################################

    P_φβ, params_φβ2ar = let

        z_labels = [:φ, :β]

        @assert tuple(z_labels...) === propertynames(C172MCS.ZLatPhiBeta())
        z_trim = lss_red.y0[z_labels]
        n_z = length(z_labels)

        F = lss_red.A
        G = lss_red.B
        Hx = lss_red.C[z_labels, :]
        Hu = lss_red.D[z_labels, :]

        Hx_int = Hx[:β, :]'
        Hu_int = Hu[:β, :]'
        n_int, _ = size(Hx_int)

        F_aug = [F zeros(n_x, n_int); Hx_int zeros(n_int, n_int)]
        G_aug = [G; Hu_int]
        Hx_aug = [Hx zeros(n_z, n_int)]
        Hu_aug = Hu

        P_aug = ss(F_aug, G_aug, Hx_aug, Hu_aug)

        @unpack v_x, v_y = lss_red.x0
        v_norm = norm([v_x, v_y])

        #weight matrices
        Q = ComponentVector(p = 0, r = 0.1, φ = 0.15, v_x = 0/v_norm, v_y = 0.1/v_norm, β_filt = 0, ail_v = 0, ail_p = 0, rud_v = 0, rud_p = 0, ξ_β = 0.001) |> diagm
        R = C172MCS.ULat(aileron_cmd = 0.1, rudder_cmd = 0.05) |> diagm

        #compute gain matrix
        C_aug = lqr(P_aug, Q, R)

        #extract system state and integrator blocks from the feedback matrix
        C_x = C_aug[:, 1:n_x]
        C_ξ = C_aug[:, n_x+1:end]

        #construct feedforward matrix blocks
        A = [F G; Hx Hu]
        B = inv(A)
        B_12 = B[1:n_x, n_x+1:end]
        B_22 = B[n_x+1:end, n_x+1:end]

        C_fbk = C_x
        C_fwd = B_22 + C_x * B_12
        C_int = ComponentMatrix(zeros(n_u, n_z), Axis(u_labels), Axis(z_labels))
        C_int[:, :β] .= C_ξ

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

        #disable warning about connecting single output to multiple inputs (here,
        #φ goes both to state feedback and command variable error junction)
        Logging.disable_logging(Logging.Warn)
        P_φβ = connect([P_lat, int_ss, C_fwd_ss, C_fbk_ss, C_int_ss,
                            φ_err_sum, β_err_sum,
                            aileron_cmd_sum, rudder_cmd_sum,
                            φ_sp_splitter, β_sp_splitter], connections;
                            w1 = z_labels_sp, z1 = P_lat.y)
        Logging.disable_logging(Logging.LogLevel(typemin(Int32)))

        #convert everything to plain arrays
        params_φβar = LQRTrackerParams(;
            C_fbk = Matrix(C_fbk), C_fwd = Matrix(C_fwd), C_int = Matrix(C_int),
            x_trim = Vector(x_trim), u_trim = Vector(u_trim), z_trim = Vector(z_trim))

        (P_φβ, params_φβar)

    end


    ############################### p + β ######################################

    P_pβ, params_p2φ = let

        P_φ2p = P_φβ[:p, :φ_sp];

        p2φ_int = tf(1, [1, 0]) |> ss
        P_p2φ_opt = series(p2φ_int, ss(P_φ2p))

        t_sim_p2φ = 10
        lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDParams(; k_p = 15.0, k_i = 35.0, k_d = 1.5, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_p2φ, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 1, ∫e = 10, ef = 2, ∫u = 0.1, up = 0.00)
        params_0 = PIDParams(; k_p = 1.5, k_i = 3, k_d = 0.1, τ_f = 0.01)

        p2φ_results = optimize_PID(P_p2φ_opt; params_0, settings, weights, global_search = false)
        params_p2φ = p2φ_results.params
        if !check_results(p2φ_results, Metrics(; Ms = 1.3, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
            println("Warning: Checks failed for roll rate control PID, design point $(design_point)")
            println(p2φ_results.metrics)
        end

        p2φ_PID = build_PID(p2φ_results.params)
        C_p2φ = named_ss(series(p2φ_int, p2φ_PID), :C_p2φ; u = :p_err, y = :φ_sp)

        p2φ_sum = sumblock("p_err = p_sp - p")
        P_pβ = connect([P_φβ, p2φ_sum, C_p2φ], [:p_err=>:p_err, :p=>:p, :φ_sp=>:φ_sp], w1 = [:p_sp, :β_sp], z1 = P_φβ.y)

        (P_pβ, params_p2φ)

    end

    ############################### χ + β ######################################

    P_χβ, params_χ2φ = let

        P_φ2χ = P_φβ[:χ, :φ_sp];

        t_sim_χ2φ = 30
        lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.3, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDParams(; k_p = 10.0, k_i = 0.3, k_d = 0.0, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_χ2φ, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 3, ∫e = 10, ef = 1, ∫u = 0.00, up = 0.01)
        params_0 = PIDParams(; k_p = 3., k_i = 0.3, k_d = 0.0, τ_f = 0.01)

        χ2φ_results = optimize_PID(P_φ2χ; params_0, settings, weights, global_search = false)

        params_χ2φ = χ2φ_results.params
        if !check_results(χ2φ_results, Metrics(; Ms = 1.4, ∫e = 0.2, ef = 0.02, ∫u = Inf, up = Inf))
            println("Warning: Checks failed for course angle control PID, design point $(design_point)")
            println(χ2φ_results.metrics)
        end

        χ2φ_PID = build_PID(χ2φ_results.params)
        C_χ2φ = named_ss(χ2φ_PID, :C_χ2φ; u = :χ_err, y = :φ_sp);

        χ2φ_sum = sumblock("χ_err = χ_sp - χ")
        P_χβ = connect([P_φβ, χ2φ_sum, C_χ2φ], [:χ_err=>:χ_err, :χ=>:χ, :φ_sp=>:φ_sp], w1 = [:χ_sp, :β_sp], z1 = P_φβ.y)

        (P_χβ, params_χ2φ)

    end

    return (φβ2ar = params_φβ2ar, p2φ = params_p2φ, χ2φ = params_χ2φ)

end


end #module