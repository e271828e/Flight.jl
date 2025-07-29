module C172XControlDesign

using Flight
using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightLib.Control.Continuous: LinearizedSS, submodel
using Flight.FlightLib.Control.Discrete: PIDParams, LQRTrackerParams, save_lookup
using Flight.FlightLib.Control.PIDOpt: Settings, Metrics, optimize_PID, build_PID, check_results

using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172X
using Flight.FlightAircraft.C172X.C172XControl

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


#sweeps the envelope in EAS and altitude to generate each control channel's gain
#scheduling lookup tables

function generate_lookups(
    EAS_range::AbstractRange{Float64} = range(25, 55, length = 7),
    h_range::AbstractRange{Float64} = range(50, 3050, length = 4);
    channel::Symbol = :lat,
    global_search::Bool = false,
    folder::String = joinpath(dirname(@__DIR__), "data"))

    if channel === :lon
        f_opt = design_lon
    elseif channel === :lat
        f_opt = design_lat
    else
        @error("Valid values for channel keyword: :lon, :lat")
        return
    end

    results = map(Iterators.product(EAS_range, h_range)) do (EAS, h)

        @info("Designing $channel controllers for EAS = $EAS, h = $h")

        #all other design point parameters at default
        flaps = EAS < 30 ? 1.0 : 0.0
        design_point = C172.TrimParameters(; Ob = Geographic(LatLon(), HEllip(h)), EAS, flaps)

        results = f_opt(; design_point, global_search)
        return results

    end |> StructArray |> StructArrays.components

    filenames = joinpath.(folder, string.(keys(results)) .* "_lookup.h5")

    bounds = ((EAS_range[1], EAS_range[end]), (h_range[1], h_range[end]))

    foreach(values(results), filenames) do data, fname
        save_lookup(data, bounds, joinpath(folder, fname))
    end

    return results

end

#automates the design process from c172x_lon.ipynb, generating the parameters
#for all controllers at the specified design point

function design_lon(; design_point::C172.TrimParameters = C172.TrimParameters(),
                    global_search = false)

    ac = Cessna172Xv0(NED()) |> Model #linearization requires NED kinematics

    lss_full = LinearizedSS(ac, design_point; model = :lon);
    P_full = named_ss(lss_full)

    x_labels_full = keys(lss_full.x0) |> collect
    y_labels_full = keys(lss_full.y0) |> collect
    u_labels_full = keys(lss_full.u0) |> collect

    #reduced model without h
    x_labels_red = copy(x_labels_full)
    y_labels_red = copy(y_labels_full)
    u_labels_red = copy(u_labels_full)

    x_labels_red = deleteat!(x_labels_red, findfirst(isequal(:h), x_labels_red))
    y_labels_red = deleteat!(y_labels_red, findfirst(isequal(:h), y_labels_red))

    lss_red = submodel(lss_full; x = x_labels_red, u = u_labels_red, y = y_labels_red)
    P_red = submodel(P_full; x = x_labels_red, u = u_labels_red, y = y_labels_red)

    #pitch dynamics model
    x_labels_pit = [:q, :θ, :v_x, :v_z, :α_filt, :ele_p]
    y_labels_pit = vcat(x_labels_pit, [:f_x, :f_z, :α, :EAS, :TAS, :γ, :climb_rate, :elevator_cmd])
    u_labels_pit = [:elevator_cmd,]

    lss_pit = Control.Continuous.submodel(lss_full; x = x_labels_pit, u = u_labels_pit, y = y_labels_pit)
    P_pit = named_ss(lss_pit);

    # ensure consistency in component selection and ordering between our pitch
    # design model and avionics implementation for state and control vectors
    @assert tuple(x_labels_pit...) === propertynames(C172XControl.XPitch())
    # @assert tuple(u_labels_pit...) === propertynames(C172XControl.UPitch())

    x_trim = lss_pit.x0
    u_trim = lss_pit.u0

    n_x = length(x_labels_pit)
    n_u = length(u_labels_pit)

    ############################ thr+ele SAS ###################################

    P_red_te, params_e2e = let

        @unpack v_x, v_z = lss_pit.x0
        v_norm = norm([v_x, v_z])

        #weight matrices
        Q = C172XControl.XPitch(q = 1, θ = 20, v_x = 0.1/v_norm, v_z = 1/v_norm, α_filt = 0, ele_p = 0) |> diagm
        # R = C172XControl.UPitch(elevator_cmd = 2) |> diagm
        R = [2.0] |> diagm #elevator_cmd weight

        #feedback gain matrix
        C_fbk = lqr(P_pit, Q, R)

        z_labels = [:elevator_cmd,]
        # @assert tuple(z_labels...) === propertynames(C172XControl.ZPitch())
        z_trim = lss_pit.y0[z_labels]
        n_z = length(z_labels)

        F = lss_pit.A
        G = lss_pit.B
        Hx = lss_pit.C[z_labels, :]
        Hu = lss_pit.D[z_labels, :]

        #forward gain matrix
        A = [F G; Hx Hu]
        B = inv(A)
        B_12 = B[1:n_x, n_x+1:end]
        B_22 = B[n_x+1:end, n_x+1:end]
        C_fwd = B_22 + C_fbk * B_12

        #no integral control
        C_int = zeros(n_u, n_z)

        #useful signal labels
        u_labels_pit_fbk = Symbol.(string.(u_labels_pit) .* "_fbk")
        u_labels_pit_fwd = Symbol.(string.(u_labels_pit) .* "_fwd")
        u_labels_pit_sum = Symbol.(string.(u_labels_pit) .* "_sum")
        z_labels_ref = Symbol.(string.(z_labels) .* "_ref")

        C_fbk_ss = named_ss(ss(C_fbk), u = x_labels_pit, y = u_labels_pit_fbk)
        C_fwd_ss = named_ss(ss(C_fwd), u = z_labels_ref, y = u_labels_pit_fwd)

        #summing junctions
        # throttle_sum = sumblock("throttle_cmd_sum = throttle_cmd_fwd- throttle_cmd_fbk")
        elevator_sum = sumblock("elevator_cmd_sum = elevator_cmd_fwd - elevator_cmd_fbk")

        connections = vcat(
            Pair.(x_labels_pit, x_labels_pit),
            Pair.(u_labels_pit_fbk, u_labels_pit_fbk),
            Pair.(u_labels_pit_sum, u_labels_pit),
            Pair.(u_labels_pit_fwd, u_labels_pit_fwd)
            )

        #connect to the complete longitudinal dynamics (this is the one we'll keep)
        P_red_te = connect([P_red, elevator_sum, C_fbk_ss, C_fwd_ss],
            connections; w1 = [:throttle_cmd, :elevator_cmd_ref], z1 = P_red.y)

        #add dummy system to rename input throttle_cmd to throttle_cmd_ref
        D_thr = named_ss(ss(1), u = [:throttle_cmd_ref,], y = [:throttle_cmd_ref, ])
        P_red_te = connect([P_red_te, D_thr], [:throttle_cmd_ref => :throttle_cmd,];
            w1 = [:throttle_cmd_ref, :elevator_cmd_ref], z1 = P_red_te.y)

        params_e2e = LQRTrackerParams(; #export everything as plain arrays
            C_fbk = Matrix(C_fbk), C_fwd = Matrix(C_fwd), C_int = Matrix(C_int),
            x_trim = Vector(x_trim), u_trim = Vector(u_trim), z_trim = Vector(z_trim))

        (P_red_te, params_e2e)

    end

    P_red_tq, params_q2e = let

        P_e2q = P_red_te[:q, :elevator_cmd_ref]

        q2e_int = tf(1, [1, 0]) |> ss
        P_q2e_opt = series(q2e_int, ss(P_e2q))

        t_sim_q2e = 10
        lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDParams(; k_p = 10.0, k_i = 50.0, k_d = 2.0, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_q2e, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 1, ∫e = 15, ef = 2, ∫u = 0.1, up = 0.00)
        params_0 = PIDParams(; k_p = 2.0, k_i = 15, k_d = 0.4, τ_f = 0.01)

        q2e_results = optimize_PID(P_q2e_opt; params_0, settings, weights, global_search)

        params_q2e = q2e_results.params
        if !check_results(q2e_results, Metrics(; Ms = 1.3, ∫e = 0.2, ef = 0.02, ∫u = Inf, up = Inf))
            @warn("Checks failed for pitch rate PID, design point $(design_point), final metrics $(q2e_results.metrics)")
        end

        q2e_pid = build_PID(q2e_results.params)
        C_q2e = named_ss(series(q2e_int, q2e_pid), :C_q2e; u = :q_err, y = :elevator_cmd_ref);

        q2e_sum = sumblock("q_err = q_ref - q")
        P_red_tq = connect([P_red_te, q2e_sum, C_q2e],
            [:q_err=>:q_err, :q=>:q, :elevator_cmd_ref=>:elevator_cmd_ref],
            w1 = [:throttle_cmd_ref, :q_ref], z1 = P_red_te.y)

        (P_red_tq, params_q2e)

    end

    P_red_tθ = let

        k_p_θ2q = 1
        C_θ2q = named_ss(ss(k_p_θ2q), :C_θ2q; u = :θ_err, y = :q_ref);

        θ2q_sum = sumblock("θ_err = θ_ref - θ")
        P_red_tθ = connect([P_red_tq, θ2q_sum, C_θ2q], [:θ_err=>:θ_err, :θ=>:θ, :q_ref=>:q_ref],
                        w1 = [:throttle_cmd_ref, :θ_ref], z1 = P_red_tq.y);

        P_red_tθ

    end

    P_red_tv, params_v2θ = let

        P_θ2v = P_red_tθ[:EAS, :θ_ref]
        P_θ2v_opt = -P_θ2v

        t_sim_v2θ = 20
        lower_bounds = PIDParams(; k_p = 0.01, k_i = 0.000, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDParams(; k_p = 0.2, k_i = 0.05, k_d = 0.0, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_v2θ, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 2.0, ∫e = 5.0, ef = 1.0, ∫u = 0.0, up = 0.0)
        params_0 = PIDParams(; k_p = 0.05, k_i = 0.01, k_d = 0.0, τ_f = 0.01)

        v2θ_results = optimize_PID(P_θ2v_opt; params_0, settings, weights, global_search)

        params_v2θ = v2θ_results.params
        if !check_results(v2θ_results, Metrics(; Ms = 1.3, ∫e = 0.2, ef = 0.02, ∫u = Inf, up = Inf))
            @warn("Checks failed for EAS to θ PID, design point $(design_point), final metrics $(v2θ_results.metrics)")
        end

        v2θ_pid = build_PID(v2θ_results.params)
        C_v2θ = -v2θ_pid
        C_v2θ = named_ss(ss(C_v2θ), :C_v2θ; u = :EAS_err, y = :θ_ref)

        v2θ_sum = sumblock("EAS_err = EAS_ref - EAS")
        P_red_tv = connect([P_red_tθ, v2θ_sum, C_v2θ], [:EAS_err=>:EAS_err, :EAS=>:EAS, :θ_ref=>:θ_ref],
        w1 = [:throttle_cmd_ref, :EAS_ref], z1 = P_red_tθ.y)

        (P_red_tv, params_v2θ)

    end

    P_red_vθ, params_v2t = let

        P_t2v = P_red_tθ[:EAS, :throttle_cmd]

        t_sim_v2t = 10
        lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDParams(; k_p = 1.5, k_i = 0.5, k_d = 0.0, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_v2t, maxeval = 5000, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 2.0, ∫e = 5.0, ef = 1.0, ∫u = 0.0, up = 0.0)
        params_0 = PIDParams(; k_p = 0.2, k_i = 0.1, k_d = 0.0, τ_f = 0.01)

        v2t_results = optimize_PID(P_t2v; params_0, settings, weights, global_search)

        params_v2t = v2t_results.params
        if !check_results(v2t_results, Metrics(; Ms = 1.3, ∫e = 0.2, ef = 0.02, ∫u = Inf, up = Inf))
            @warn("Checks failed for EAS to throttle PID, design point $(design_point), final metrics $(v2t_results.metrics)")
        end

        v2t_pid = build_PID(v2t_results.params)
        C_v2t = named_ss(ss(v2t_pid), :C_v2t; u = :EAS_err, y = :throttle_cmd_ref)

        v2t_sum = sumblock("EAS_err = EAS_ref - EAS")
        P_red_vθ = connect([P_red_tθ, v2t_sum, C_v2t],
            [:EAS_err=>:EAS_err, :EAS=>:EAS, :throttle_cmd_ref=>:throttle_cmd_ref],
            w1 = [:EAS_ref, :θ_ref], z1 = P_red_tθ.y)

        (P_red_vθ, params_v2t)

    end

    P_red_vc, params_c2θ = let

        P_θ2c = P_red_vθ[:climb_rate, :θ_ref]

        t_sim_c2θ = 10
        lower_bounds = PIDParams(; k_p = 0.001, k_i = 0.001, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDParams(; k_p = 0.05, k_i = 0.03, k_d = 0.0, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_c2θ, maxeval = 5000, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 2.0, ∫e = 3.0, ef = 1.0, ∫u = 0.0, up = 0.1)
        params_0 = PIDParams(; k_p = 0.02, k_i = 0.01, k_d = 0.0, τ_f = 0.01)

        c2θ_results = optimize_PID(P_θ2c; params_0, settings, weights, global_search)

        params_c2θ = c2θ_results.params
        if !check_results(c2θ_results, Metrics(; Ms = 1.3, ∫e = 0.2, ef = 0.02, ∫u = Inf, up = Inf))
            @warn("Checks failed for climb rate to θ PID, design point $(design_point), final metrics $(c2θ_results.metrics)")
        end

        c2θ_PID = build_PID(c2θ_results.params)
        C_c2θ = named_ss(ss(c2θ_PID), :C_c2θ; u = :climb_rate_err, y = :θ_ref)

        c2θ_sum = sumblock("climb_rate_err = climb_rate_ref - climb_rate")
        P_red_vc = connect([P_red_vθ, c2θ_sum, C_c2θ],
            [:climb_rate_err=>:climb_rate_err, :climb_rate=>:climb_rate, :θ_ref=>:θ_ref],
            w1 = [:EAS_ref, :climb_rate_ref], z1 = P_red_vθ.y)

        (P_red_vc, params_c2θ)
    end

    return (e2e = params_e2e, q2e = params_q2e, v2θ = params_v2θ,
            v2t = params_v2t, c2θ = params_c2θ)

end

#automates the design process from c172x_lat.ipynb, generating the parameters
#for all controllers at the specified design point
function design_lat(; design_point::C172.TrimParameters = C172.TrimParameters(),
                    global_search::Bool = false)

    ac = Cessna172Xv0(NED()) |> Model #linearization requires NED kinematics

    lss_full = Control.Continuous.LinearizedSS(ac, design_point; model = :lat);
    P_full = named_ss(lss_full)

    x_labels_full = keys(lss_full.x0) |> collect
    y_labels_full = keys(lss_full.y0) |> collect
    u_labels = keys(lss_full.u0) |> collect

    x_labels_red = copy(x_labels_full)
    y_labels_red = copy(y_labels_full)

    x_labels_red = deleteat!(x_labels_red, findfirst(isequal(:ψ), x_labels_red))
    y_labels_red = deleteat!(y_labels_red, findfirst(isequal(:ψ), y_labels_red))
    y_labels_red = deleteat!(y_labels_red, findfirst(isequal(:χ), y_labels_red))

    #ensure consistency in component selection and ordering between our design model
    #and avionics implementation for state and control vectors
    @assert tuple(x_labels_red...) === propertynames(C172XControl.XLat())
    @assert tuple(u_labels...) === propertynames(C172XControl.ULat())

    #extract reduced design model
    lss_red = Control.Continuous.submodel(lss_full; x = x_labels_red, u = u_labels, y = y_labels_full)
    P_red = named_ss(lss_red);

    x_trim = lss_red.x0
    u_trim = lss_red.u0

    n_x = length(lss_red.x0)
    n_u = length(lss_red.u0)

    ################################ SAS #######################################

    P_ar, params_ar2ar = let

        #useful signal labels for connections
        z_labels = [:aileron_cmd, :rudder_cmd]
        z_trim = lss_red.y0[z_labels]
        n_z = length(z_labels)

        @unpack v_x, v_y = lss_red.x0
        v_norm = norm([v_x, v_y])

        #weight matrices
        Q = C172XControl.XLat(p = 0, r = 0.1, φ = 0.1, v_x = 0/v_norm, v_y = 0/v_norm, β_filt = 0, ail_p = 0, rud_p = 0) |> diagm
        R = C172XControl.ULat(aileron_cmd = 0.1, rudder_cmd = 0.01) |> diagm

        #feedback gain matrix
        C_fbk = lqr(P_red, Q, R)

        #passthrough feedforward
        C_fwd = Matrix{Float64}(I, n_z, n_z)

        #no integral control
        C_int = zeros(n_u, n_z)

        u_labels_fbk = Symbol.(string.(u_labels) .* "_fbk")
        u_labels_fwd = Symbol.(string.(u_labels) .* "_fwd")
        u_labels_sum = Symbol.(string.(u_labels) .* "_sum")
        z_labels_ref = Symbol.(string.(z_labels) .* "_ref")

        C_fbk_ss = named_ss(ss(C_fbk); u = x_labels_red, y = u_labels_fbk)
        C_fwd_ss = named_ss(ss(C_fwd), u = z_labels_ref, y = u_labels_fwd)

        #summing junctions
        aileron_sum = sumblock("aileron_cmd_sum = aileron_cmd_fwd- aileron_cmd_fbk")
        rudder_sum = sumblock("rudder_cmd_sum = rudder_cmd_fwd - rudder_cmd_fbk")

        connections_fbk = vcat(
            Pair.(x_labels_red, x_labels_red),
            Pair.(u_labels_fbk, u_labels_fbk),
            Pair.(u_labels_fwd, u_labels_fwd),
            Pair.(u_labels_sum, u_labels),
            )

        P_ar = connect([P_full, aileron_sum, rudder_sum, C_fbk_ss, C_fwd_ss],
                        connections_fbk; w1 = z_labels_ref, z1 = P_full.y)

        params_ar2ar = LQRTrackerParams(;
            C_fbk = Matrix(C_fbk), C_fwd = Matrix(C_fwd), C_int = Matrix(C_int),
            x_trim = Vector(x_trim), u_trim = Vector(u_trim), z_trim = Vector(z_trim))

        (P_ar, params_ar2ar)

    end

    ############################### φ + β ######################################

    P_φβ, params_φβ2ar = let

        z_labels = [:φ, :β]

        @assert tuple(z_labels...) === propertynames(C172XControl.ZLat())
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

        @unpack v_x, v_y = lss_red.x0
        v_norm = norm([v_x, v_y])

        #weight matrices
        Q = ComponentVector(p = 0, r = 0.1, φ = 0.25, v_x = 0/v_norm, v_y = 0.1/v_norm,
                            β_filt = 0, ail_p = 0, rud_p = 0, ξ_φ = 0.1, ξ_β = 0.001) |> diagm
        R = C172XControl.ULat(aileron_cmd = 0.05, rudder_cmd = 0.05) |> diagm

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
        C_int = C_ξ

        #some useful signal labels
        u_labels_fbk = Symbol.(string.(u_labels) .* "_fbk")
        u_labels_fwd = Symbol.(string.(u_labels) .* "_fwd")
        u_labels_sum = Symbol.(string.(u_labels) .* "_sum")
        u_labels_int_u = Symbol.(string.(u_labels) .* "_int_u")
        u_labels_int = Symbol.(string.(u_labels) .* "_int")
        u_labels_ξ = Symbol.(string.(u_labels) .* "_ξ")

        z_labels_ref = Symbol.(string.(z_labels) .* "_ref")
        z_labels_ref1 = Symbol.(string.(z_labels) .* "_ref1")
        z_labels_ref2 = Symbol.(string.(z_labels) .* "_ref2")
        z_labels_err = Symbol.(string.(z_labels) .* "_err")
        z_labels_sum = Symbol.(string.(z_labels) .* "_sum")
        z_labels_ref_fwd = Symbol.(string.(z_labels) .* "_ref_fwd")
        z_labels_ref_sum = Symbol.(string.(z_labels) .* "_ref_sum")

        C_fbk_ss = named_ss(ss(C_fbk), u = x_labels_red, y = u_labels_fbk)
        C_fwd_ss = named_ss(ss(C_fwd), u = z_labels_ref_fwd, y = u_labels_fwd)
        C_int_ss = named_ss(ss(C_int), u = z_labels_err, y = u_labels_int_u)

        int_ss = named_ss(ss(tf(1, [1,0])) .* I(2),
                            x = u_labels_ξ,
                            u = u_labels_int_u,
                            y = u_labels_int);

        φ_err_sum = sumblock("φ_err = φ_sum - φ_ref_sum")
        β_err_sum = sumblock("β_err = β_sum - β_ref_sum")

        aileron_cmd_sum = sumblock("aileron_cmd_sum = aileron_cmd_fwd - aileron_cmd_fbk - aileron_cmd_int")
        rudder_cmd_sum = sumblock("rudder_cmd_sum = rudder_cmd_fwd - rudder_cmd_fbk - rudder_cmd_int")

        φ_ref_splitter = splitter(:φ_ref, 2)
        β_ref_splitter = splitter(:β_ref, 2)

        connections = vcat(
            Pair.(x_labels_red, x_labels_red),
            Pair.(z_labels, z_labels_sum),
            Pair.(z_labels_ref1, z_labels_ref_sum),
            Pair.(z_labels_ref2, z_labels_ref_fwd),
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
        P_φβ = connect([P_full, int_ss, C_fwd_ss, C_fbk_ss, C_int_ss,
                            φ_err_sum, β_err_sum,
                            aileron_cmd_sum, rudder_cmd_sum,
                            φ_ref_splitter, β_ref_splitter], connections;
                            w1 = z_labels_ref, z1 = P_full.y)
        Logging.disable_logging(Logging.LogLevel(typemin(Int32)))

        #convert everything to plain arrays
        params_φβar = LQRTrackerParams(;
            C_fbk = Matrix(C_fbk), C_fwd = Matrix(C_fwd), C_int = Matrix(C_int),
            x_trim = Vector(x_trim), u_trim = Vector(u_trim), z_trim = Vector(z_trim))

        (P_φβ, params_φβar)

    end


    ############################### p + β ######################################

    P_pβ, params_p2φ = let

        P_φ2p = P_φβ[:p, :φ_ref];

        p2φ_int = tf(1, [1, 0]) |> ss
        P_p2φ_opt = series(p2φ_int, ss(P_φ2p))

        t_sim_p2φ = 10
        lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDParams(; k_p = 10.0, k_i = 35.0, k_d = 1.5, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_p2φ, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 0, ∫e = 2, ef = 2, ∫u = 1, up = 0.00)
        params_0 = PIDParams(; k_p = 1.5, k_i = 3, k_d = 0.1, τ_f = 0.01)

        p2φ_results = optimize_PID(P_p2φ_opt; params_0, settings, weights, global_search)
        params_p2φ = p2φ_results.params
        if !check_results(p2φ_results, Metrics(; Ms = Inf, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
            @warn("Checks failed for p to φ PID, design point $(design_point), final metrics $(p2φ_results.metrics)")
        end

        p2φ_PID = build_PID(p2φ_results.params)
        C_p2φ = named_ss(series(p2φ_int, p2φ_PID), :C_p2φ; u = :p_err, y = :φ_ref)

        p2φ_sum = sumblock("p_err = p_ref - p")
        P_pβ = connect([P_φβ, p2φ_sum, C_p2φ], [:p_err=>:p_err, :p=>:p, :φ_ref=>:φ_ref], w1 = [:p_ref, :β_ref], z1 = P_φβ.y)

        (P_pβ, params_p2φ)

    end

    ############################### χ + β ######################################

    P_χβ, params_χ2φ = let

        P_φ2χ = P_φβ[:χ, :φ_ref];

        t_sim_χ2φ = 30
        lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.3, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDParams(; k_p = 10.0, k_i = 0.3, k_d = 0.0, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_χ2φ, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 3, ∫e = 10, ef = 1, ∫u = 0.00, up = 0.01)
        params_0 = PIDParams(; k_p = 3., k_i = 0.3, k_d = 0.0, τ_f = 0.01)

        χ2φ_results = optimize_PID(P_φ2χ; params_0, settings, weights, global_search)

        params_χ2φ = χ2φ_results.params
        if !check_results(χ2φ_results, Metrics(; Ms = 1.4, ∫e = 0.2, ef = 0.02, ∫u = Inf, up = Inf))
            @warn("Checks failed for χ to φ PID, design point $(design_point), final metrics $(χ2φ_results.metrics)")
        end

        χ2φ_PID = build_PID(χ2φ_results.params)
        C_χ2φ = named_ss(χ2φ_PID, :C_χ2φ; u = :χ_err, y = :φ_ref);

        χ2φ_sum = sumblock("χ_err = χ_ref - χ")
        P_χβ = connect([P_φβ, χ2φ_sum, C_χ2φ], [:χ_err=>:χ_err, :χ=>:χ, :φ_ref=>:φ_ref], w1 = [:χ_ref, :β_ref], z1 = P_φβ.y)

        (P_χβ, params_χ2φ)

    end

    return (ar2ar = params_ar2ar, φβ2ar = params_φβ2ar, p2φ = params_p2φ, χ2φ = params_χ2φ)

end


end #module