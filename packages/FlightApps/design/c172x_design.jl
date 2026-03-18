module C172XDesign

using FlightCore, FlightPhysics, FlightApps

using FlightPhysics.Control: PIDData, LQRDataPoint, save_lookup_data
using FlightApps.C172
using FlightApps.C172X
using FlightApps.C172X.C172XControl

using HDF5
using Logging
using ControlSystems
using RobustAndOptimalControl
using StaticArrays
using StructArrays
using ComponentArrays
using LinearAlgebra
using Interpolations

include("pidopt.jl")
using .PIDOpt: Settings, Metrics, optimize_PID, build_PID, check_results

export get_design_model!, generate_lookups


function get_design_model!(aircraft::Model{<:C172X.Aircraft{NED}}, args...; kwargs...)
    get_design_model!(aircraft.vehicle, args...; kwargs...)
end

function get_design_model!(
            vehicle::Model{<:C172X.Vehicle{NED}},
            trim_params::C172.TrimParameters = C172.TrimParameters();
            model::Symbol = :full)

    lss = linearize(vehicle, trim_params)

    (; A, B, C, D, ẋ0, x0, y0) = lss

    #replace v_x, v_y, v_z, ω_eng with TAS, α, β, n_eng as state variables
    xp_labels = x_labels = keys(lss.x0) |> collect
    for (old, new) in zip((:v_x, :v_y, :v_z, :ω_eng), (:EAS, :α, :β, :n_eng))
        xp_labels[findfirst(isequal(old), xp_labels)] = new
    end

    #define new trim values
    ẋp0 = ComponentVector(copy(getdata(lss.ẋ0)), Axis(xp_labels))
    ẋp0[(:EAS, :α, :β, :n_eng)] .= 0 #guaranteed by trim constraints
    xp0 = lss.y0[xp_labels]
    up0 = lss.u0
    yp0 = lss.y0

    #compute new state-space matrices via similarity transformation
    T = lss.C[xp_labels, :]
    T_inv = ComponentMatrix(inv(T), Axis(x_labels), Axis(xp_labels))

    Ap = T * A * T_inv
    Bp = T * B
    Cp = C * T_inv
    Dp = D

    #rebuild the model
    lss = LinearizedSS(ẋ0 = ẋp0, x0 = xp0, u0 = up0, y0 = yp0,
                      A = Ap, B = Bp, C = Cp, D = Dp)

    if model === :full
        return lss

    elseif model === :lon
        x_labels = [:q, :θ, :EAS, :α, :h, :α_filt, :n_eng, :thr_p, :ele_p]
        u_labels = [:throttle_cmd, :elevator_cmd]
        y_labels = vcat(x_labels, [:f_x, :f_z, :TAS, :γ, :climb_rate, :throttle_cmd, :elevator_cmd])
        return subsystem(lss; x = x_labels, u = u_labels, y = y_labels)

    elseif model === :lat
        x_labels = [:p, :r, :ψ, :φ, :EAS, :β, :β_filt, :ail_p, :rud_p]
        u_labels = [:aileron_cmd, :rudder_cmd]
        y_labels = vcat(x_labels, [:f_y, :χ, :aileron_cmd, :rudder_cmd])
        return subsystem(lss; x = x_labels, u = u_labels, y = y_labels)

    else
        error("Valid model keyword values: :full, :lon, :lat")

    end

end

#sweep the envelope in EAS and altitude to generate each control channel's gain
#scheduling lookup tables

function generate_lookups(
    EAS_range::AbstractRange{Float64} = range(25, 55, length = 2), #7
    h_range::AbstractRange{Float64} = range(50, 3050, length = 2); #4
    channel::Symbol = :lon,
    global_search::Bool = false,
    folder::String = joinpath(dirname(@__DIR__), normpath("src/c172/c172x/control/data")))


        P_aug = ss(A_aug, B_aug, C_aug, D_aug)

        #weight matrices
        Q = ComponentVector(q = 20, θ = 100, EAS = 0.06, α = 0, h = 0.04, α_filt = 0,
                            n_eng = 0, thr_p = 0, ele_p = 0,
                            ξ_EAS = 0.005, ξ_h = 0.001) |> diagm
        R = C172XControl.ULon(throttle_cmd = 0.1, elevator_cmd = 0.05) |> diagm

        #compute gain matrix
        K_aug = lqr(P_aug, Q, R)

        L =[A B; C D]
        M = inv(L)
        M_12 = M[1:n_x, n_x+1:end]
        M_22 = M[n_x+1:end, n_x+1:end]

        #extract system state and integrator blocks from the feedback matrix
        K_x = K_aug[:, 1:n_x]
        K_ξ = K_aug[:, n_x+1:end]

        K_fbk = K_x
        K_fwd = M_22 + K_x * M_12
        K_int = K_ξ

        #some useful signal labels
        u_labels_fbk = Symbol.(string.(u_labels) .* "_fbk")
        u_labels_fwd = Symbol.(string.(u_labels) .* "_fwd")
        u_labels_sum = Symbol.(string.(u_labels) .* "_sum")
        u_labels_int_in = Symbol.(string.(u_labels) .* "_int_in")
        u_labels_int_out = Symbol.(string.(u_labels) .* "_int_out")
        u_labels_ξ = Symbol.(string.(u_labels) .* "_ξ")

        z_labels_ref = Symbol.(string.(z_labels) .* "_ref")
        z_labels_err = Symbol.(string.(z_labels) .* "_err")

        K_fbk_ss = named_ss(ss(K_fbk), u = x_labels, y = u_labels_fbk)
        K_fwd_ss = named_ss(ss(K_fwd), u = z_labels_ref, y = u_labels_fwd)
        K_int_ss = named_ss(ss(K_int), u = z_labels_err, y = u_labels_int_in)

        int_ss = named_ss(ss(tf(1, [1,0])) .* I(2),
                            x = u_labels_ξ,
                            u = u_labels_int_in,
                            y = u_labels_int_out);

        EAS_err_sum = sumblock("EAS_err = EAS - EAS_ref")
        h_err_sum = sumblock("h_err = h - h_ref")

        throttle_cmd_sum = sumblock("throttle_cmd_sum = throttle_cmd_fwd - throttle_cmd_fbk - throttle_cmd_int_out")
        elevator_cmd_sum = sumblock("elevator_cmd_sum = elevator_cmd_fwd - elevator_cmd_fbk - elevator_cmd_int_out")

        connections = vcat(
            Pair.(x_labels, x_labels),
            Pair.(z_labels_err, z_labels_err),
            Pair.(u_labels_fwd, u_labels_fwd),
            Pair.(u_labels_fbk, u_labels_fbk),
            Pair.(u_labels_int_out, u_labels_int_out),
            Pair.(u_labels_int_in, u_labels_int_in),
            Pair.(u_labels_sum, u_labels),
            )

        #disable warning about connecting single output to multiple inputs
        #(here, EAS and h go both to state feedback and command variable error
        #junction)
        Logging.disable_logging(Logging.Warn)
        P_vh = connect([P_lon, int_ss, K_fwd_ss, K_fbk_ss, K_int_ss,
                        EAS_err_sum, h_err_sum,
                        throttle_cmd_sum, elevator_cmd_sum], connections;
                        w1 = z_labels_ref, z1 = P_lon.y, unique = false)
        Logging.disable_logging(Logging.LogLevel(typemin(Int32)))

        data_vh2te = LQRDataPoint(; K_fbk, K_fwd, K_int, x_trim, u_trim, z_trim)

        (P_vh, data_vh2te)

    end

    return (te2te = data_te2te, q2e = data_q2e, v2t = data_v2t,
            c2θ = data_c2θ, tv2te = data_tv2te, vh2te = data_vh2te)

end

#automate the design process from c172x_lat.ipynb, generating the parameters
#for all controllers at the specified design point
function design_lat(; design_point::C172.TrimParameters = C172.TrimParameters(),
                    global_search::Bool = false)

    aircraft = Cessna172Xv0(NED()) |> Model #linearization requires NED kinematics

    #complete lateral model
    lss_lat = get_design_model!(aircraft, design_point; model = :lat);
    P_lat = named_ss(lss_lat);

    #reduced lateral model
    lss_red = delete_vars(lss_lat, (:ψ, :χ))
    P_red = named_ss(lss_red);

    ################################ SAS #######################################

    P_ar, data_ar2ar = let

        x_trim = lss_red.x0
        n_x = length(x_trim)
        x_labels = collect(keys(x_trim))
        #ensure consistency in component selection and ordering between
        #design model and avionics implementation for state and control vectors
        @assert tuple(x_labels...) === propertynames(C172XControl.XLatRed())

        u_trim = lss_red.u0
        n_u = length(u_trim)
        u_labels = collect(keys(u_trim))
        @assert tuple(u_labels...) === propertynames(C172XControl.ULatRed())

        z_labels = [:aileron_cmd, :rudder_cmd]
        z_trim = lss_red.y0[z_labels]
        n_z = length(z_trim)
        @assert tuple(z_labels...) === propertynames(C172XControl.Zar())

        #weight matrices
        Q = C172XControl.XLatRed(p = 0, r = 0.1, φ = 0.1, EAS = 0.0, β = 0, β_filt = 0, ail_p = 0, rud_p = 0) |> diagm
        R = C172XControl.ULatRed(aileron_cmd = 0.1, rudder_cmd = 0.01) |> diagm

        #feedback gain matrix
        K_fbk = lqr(P_red, Q, R)

        #passthrough feedforward
        K_fwd = Matrix{Float64}(I, n_z, n_z)

        #no integral control
        K_int = zeros(n_u, n_z)

        u_labels_fbk = Symbol.(string.(u_labels) .* "_fbk")
        u_labels_fwd = Symbol.(string.(u_labels) .* "_fwd")
        u_labels_sum = Symbol.(string.(u_labels) .* "_sum")
        z_labels_ref = Symbol.(string.(z_labels) .* "_ref")

        K_fbk_ss = named_ss(ss(K_fbk); u = x_labels, y = u_labels_fbk)
        K_fwd_ss = named_ss(ss(K_fwd), u = z_labels_ref, y = u_labels_fwd)

        #summing junctions
        aileron_sum = sumblock("aileron_cmd_sum = aileron_cmd_fwd- aileron_cmd_fbk")
        rudder_sum = sumblock("rudder_cmd_sum = rudder_cmd_fwd - rudder_cmd_fbk")

        connections_fbk = vcat(
            Pair.(x_labels, x_labels),
            Pair.(u_labels_fbk, u_labels_fbk),
            Pair.(u_labels_fwd, u_labels_fwd),
            Pair.(u_labels_sum, u_labels),
            )

        P_ar = connect([P_lat, aileron_sum, rudder_sum, K_fbk_ss, K_fwd_ss],
                        connections_fbk; w1 = z_labels_ref, z1 = P_lat.y)

        data_ar2ar = LQRDataPoint(; K_fbk, K_fwd, K_int, x_trim, u_trim, z_trim)

        (P_ar, data_ar2ar)

    end

    ############################### φ + β ######################################

    P_φβ, data_φβ2ar = let lss = lss_red

        x_trim = lss.x0
        n_x = length(x_trim)
        x_labels = collect(keys(x_trim))
        @assert tuple(x_labels...) === propertynames(C172XControl.XLatRed())

        u_trim = lss.u0
        n_u = length(u_trim)
        u_labels = collect(keys(u_trim))
        @assert tuple(u_labels...) === propertynames(C172XControl.ULatRed())

        z_labels = [:φ, :β]
        z_trim = lss.y0[z_labels]
        n_z = length(z_labels)
        @assert tuple(z_labels...) === propertynames(C172XControl.Zφβ())

        ################################ feedback ###################################

        #weight matrices
        Q = C172XControl.XLatRed(p = 0, r = 0.1, φ = 2, EAS = 0, β = 5, β_filt = 0, ail_p = 0, rud_p = 0) |> diagm
        R = C172XControl.ULatRed(aileron_cmd = 0.1, rudder_cmd = 0.03) |> diagm

        #feedback gain matrix
        P = named_ss(lss)
        K_fbk = lqr(P, Q, R)

        ################################ feedforward ###########################

        A = lss.A
        B = lss.B
        C = lss.C[z_labels, :]
        D = lss.D[z_labels, :]

        #useful signal labels for connections
        u_labels_fbk = Symbol.(string.(u_labels) .* "_fbk") #outputs from feedback block
        u_labels_fwd = Symbol.(string.(u_labels) .* "_fwd") #outputs from feedforward block
        u_labels_sum = Symbol.(string.(u_labels) .* "_sum") #outputs from summing junctions
        u_labels_ref = Symbol.(string.(u_labels) .* "_ref") #references, inputs to P
        z_labels_ref = Symbol.(string.(z_labels) .* "_ref")

        L =[A B; C D]
        M = inv(L)
        M_12 = M[1:n_x, n_x+1:end]
        M_22 = M[n_x+1:end, n_x+1:end]
        K_fwd = M_22 + K_fbk * M_12
        K_fwd_ss = named_ss(ss(K_fwd), u = z_labels_ref, y = u_labels_fwd)

        #no integral control
        K_int = zeros(n_u, n_z)

        K_fbk_ss = named_ss(ss(K_fbk); u = x_labels, y = u_labels_fbk)
        K_fwd_ss = named_ss(ss(K_fwd), u = z_labels_ref, y = u_labels_fwd)

        #summing junctions
        aileron_sum = sumblock("aileron_cmd_sum = aileron_cmd_fwd- aileron_cmd_fbk")
        rudder_sum = sumblock("rudder_cmd_sum = rudder_cmd_fwd - rudder_cmd_fbk")

        connections_fbk = vcat(
            Pair.(x_labels, x_labels),
            Pair.(u_labels_fbk, u_labels_fbk),
            Pair.(u_labels_fwd, u_labels_fwd),
            Pair.(u_labels_sum, u_labels),
            )

        P_φβ = connect([P_lat, aileron_sum, rudder_sum, K_fbk_ss, K_fwd_ss],
                        connections_fbk; external_inputs = z_labels_ref, external_outputs = P_lat.y);

        data_φβar = LQRDataPoint(; K_fbk, K_fwd, K_int, x_trim, u_trim, z_trim)

        (P_φβ, data_φβar)

    end


    ############################### p + β ######################################

    P_pβ, data_p2φ = let

        P_φ2p = P_φβ[:p, :φ_ref];

        p2φ_int = tf(1, [1, 0]) |> ss
        P_p2φ_opt = series(p2φ_int, ss(P_φ2p))

        t_sim_p2φ = 10
        lower_bounds = PIDData(; k_p = 0.1, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDData(; k_p = 10.0, k_i = 35.0, k_d = 1.5, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_p2φ, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 0, ∫e = 2, ef = 2, ∫u = 1, up = 0.00)
        data_0 = PIDData(; k_p = 1.5, k_i = 3, k_d = 0.1, τ_f = 0.01)

        p2φ_results = optimize_PID(P_p2φ_opt; data_0, settings, weights, global_search)
        data_p2φ = p2φ_results.data
        if !check_results(p2φ_results, Metrics(; Ms = Inf, ∫e = 0.1, ef = 0.04, ∫u = Inf, up = Inf))
            @warn("Checks failed for p to φ PID, design point $(design_point), final metrics $(p2φ_results.metrics)")
        end

        p2φ_PID = build_PID(p2φ_results.data)
        C_p2φ = named_ss(series(p2φ_int, p2φ_PID), :C_p2φ; u = :p_err, y = :φ_ref)

        p2φ_sum = sumblock("p_err = p_ref - p")
        P_pβ = connect([P_φβ, p2φ_sum, C_p2φ], [:p_err=>:p_err, :p=>:p, :φ_ref=>:φ_ref], w1 = [:p_ref, :β_ref], z1 = P_φβ.y)

        (P_pβ, data_p2φ)

    end

    ############################### χ + β ######################################

    P_χβ, data_χ2φ = let

        P_φ2χ = P_φβ[:χ, :φ_ref];

        t_sim_χ2φ = 30
        lower_bounds = PIDData(; k_p = 0.1, k_i = 0.4, k_d = 0.0, τ_f = 0.01)
        upper_bounds = PIDData(; k_p = 10.0, k_i = 0.4, k_d = 1.5, τ_f = 0.01)
        settings = Settings(; t_sim = t_sim_χ2φ, lower_bounds, upper_bounds)
        weights = Metrics(; Ms = 3, ∫e = 10, ef = 1, ∫u = 0.00, up = 0.01)
        data_0 = PIDData(; k_p = 3., k_i = 0.4, k_d = 0.0, τ_f = 0.01)

        χ2φ_results = optimize_PID(P_φ2χ; data_0, settings, weights, global_search)

        data_χ2φ = χ2φ_results.data
        if !check_results(χ2φ_results, Metrics(; Ms = 2, ∫e = 0.2, ef = 0.04, ∫u = Inf, up = Inf))
            @warn("Checks failed for χ to φ PID, design point $(design_point), final metrics $(χ2φ_results.metrics)")
        end

        χ2φ_PID = build_PID(χ2φ_results.data)
        C_χ2φ = named_ss(χ2φ_PID, :C_χ2φ; u = :χ_err, y = :φ_ref);

        χ2φ_sum = sumblock("χ_err = χ_ref - χ")
        P_χβ = connect([P_φβ, χ2φ_sum, C_χ2φ], [:χ_err=>:χ_err, :χ=>:χ, :φ_ref=>:φ_ref], w1 = [:χ_ref, :β_ref], z1 = P_φβ.y)

        (P_χβ, data_χ2φ)

    end

    return (ar2ar = data_ar2ar, φβ2ar = data_φβ2ar, p2φ = data_p2φ, χ2φ = data_χ2φ)

end


end #module
