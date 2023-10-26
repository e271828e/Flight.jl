module CASLookups

using Flight
using Flight.FlightCore.Systems
using Flight.FlightCore.Plotting

using Flight.FlightPhysics.Geodesy
using Flight.FlightPhysics.Kinematics
using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.Control
using Flight.FlightComponents.Aircraft
using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172FBWCAS

using UnPack
using ControlSystems
using RobustAndOptimalControl
using StructArrays


function generate_lookups(
    EAS_range::AbstractRange{Float64} = range(25, 55, length = 4),
    h_range::AbstractRange{Float64} = range(10, 3010, length = 4);
    channel::Symbol = :pitch,
    global_search::Bool = false,
    folder::String = dirname(@__DIR__)) #save to parent folder by default

    ac = Cessna172FBWBase(NED()) |> System

    mid_cg_pld = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

    if channel === :pitch
        f_opt = optimize_pitch
    elseif channel === :roll
        f_opt = optimize_roll
    else
        error("Valid values for channel keyword: :pitch, :roll")
    end

    results = map(Iterators.product(EAS_range, h_range)) do (EAS, h)

        println("Optimizing for EAS = $EAS, h = $h")

        flaps = EAS < 30 ? 1.0 : 0.0

        point = C172FBW.TrimParameters(; Ob = Geographic(LatLon(), HOrth(h)),
            EAS, flaps, γ_wOb_n = 0.0, x_fuel = 0.5, payload = mid_cg_pld)

        results = f_opt(ac; point, global_search)
        return results

    end |> StructArray |> StructArrays.components

    filenames = joinpath.(dirname(@__DIR__), "data", string.(keys(results)) .* "_lookup.h5")

    EAS_bounds = (EAS_range[1], EAS_range[end])
    h_bounds = (h_range[1], h_range[end])

    map(values(results), filenames) do results, fname

        data = PIDParams(StructArrays.components(StructArray(StructArray(results).params))...)
        lookup = C172FBWCAS.Lookup(data, EAS_bounds, h_bounds)
        C172FBWCAS.save_lookup(lookup, joinpath(folder, fname))
        return lookup
    end

end

#optimizes PID parameters in the pitch channel compensators
function optimize_pitch(ac::System{<:Cessna172FBWBase{NED}};
                    point::C172FBW.TrimParameters = C172FBW.TrimParameters(),
                    global_search::Bool = false)

    thr_ele_MIMO = named_ss(ac, point; model = :lon);

    ################################ q2e #######################################
    P_e2q = thr_ele_MIMO[:q, :elevator_cmd]
    q2e_int = tf(1, [1, 0]) |> ss
    P_q2e_opt = series(q2e_int, ss(P_e2q))

    t_sim_q2e = 10
    lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
    upper_bounds = PIDParams(; k_p = 50.0, k_i = 35.0, k_d = 1.5, τ_f = 0.01)
    settings = PIDOpt.Settings(; t_sim = t_sim_q2e, lower_bounds, upper_bounds)
    weights = PIDOpt.Metrics(; Ms = 2, ∫e = 10, ef = 2, ∫u = 0.1, up = 0.00)
    params_0 = PIDParams(; k_p = 3, k_i = 15, k_d = 0.5, τ_f = 0.01)

    q2e_results = PIDOpt.optimize_PID(P_q2e_opt; params_0, settings, weights, global_search)

    if !PIDOpt.check_results(q2e_results, PIDOpt.Metrics(; Ms = 1.5, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
        println("Warning: Checks failed for q design point $point")
        println(q2e_results.metrics)
    end

    q2e_PID = PIDOpt.build_PID(q2e_results.params)
    C_q2e = named_ss(series(q2e_int, q2e_PID), :C_q2e; u = :q_err, y = :elevator_cmd);
    q2e_sum = sumblock("q_err = q_dmd - q")
    thr_q_MIMO = connect([q2e_sum, C_q2e, thr_ele_MIMO], [:q_err=>:q_err, :q=>:q, :elevator_cmd=>:elevator_cmd], w1 = [:throttle_cmd, :q_dmd], z1 = thr_ele_MIMO.y)

    ################################ θ2q #######################################
    P_q2θ = thr_q_MIMO[:θ, :q_dmd]

    t_sim_θ2q = 10
    upper_bounds = PIDParams(; k_p = 50.0, k_i = 0.0, k_d = 5.0, τ_f = 0.01)
    settings = PIDOpt.Settings(; t_sim = t_sim_θ2q, maxeval = 5000, upper_bounds)
    weights = PIDOpt.Metrics(; Ms = 2.0, ∫e = 5.0, ef = 1.0, ∫u = 0.0, up = 0.1)
    params_0 = PIDParams(; k_p = 2.0, k_i = 0.0, k_d = 0.0, τ_f = 0.01)

    θ2q_results = PIDOpt.optimize_PID(P_q2θ; params_0, settings, weights, global_search)

    if !PIDOpt.check_results(θ2q_results, PIDOpt.Metrics(; Ms = 1.5, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
        println("Warning: Checks failed for θ2q design point $point")
        println(θ2q_results.metrics)
    end

    θ2q_PID = PIDOpt.build_PID(θ2q_results.params)
    C_θ2q = named_ss(θ2q_PID, :C_θ2q; u = :θ_err, y = :q_dmd)
    θ2q_sum = sumblock("θ_err = θ_dmd - θ")
    thr_θ_MIMO = connect([θ2q_sum, C_θ2q, thr_q_MIMO], [:θ_err=>:θ_err, :θ=>:θ, :q_dmd=>:q_dmd], w1 = [:throttle_cmd, :θ_dmd], z1 = thr_q_MIMO.y)

    ################################ v2t #######################################
    P_t2v = thr_θ_MIMO[:EAS, :throttle_cmd]

    t_sim_v2t = 10
    lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
    upper_bounds = PIDParams(; k_p = 4.0, k_i = 0.5, k_d = 0.0, τ_f = 0.01)
    settings = PIDOpt.Settings(; t_sim = t_sim_v2t, maxeval = 5000, lower_bounds, upper_bounds)
    weights = PIDOpt.Metrics(; Ms = 1.0, ∫e = 10.0, ef = 1.0, ∫u = 0.0, up = 0.0)
    params_0 = PIDParams(; k_p = 0.5, k_i = 0.1, k_d = 0.0, τ_f = 0.01)

    v2t_results = PIDOpt.optimize_PID(P_t2v; params_0, settings, weights, global_search = false)
    v2t_PID = PIDOpt.build_PID(v2t_results.params)

    v2t_results = PIDOpt.optimize_PID(P_t2v; params_0, settings, weights, global_search)

    if !PIDOpt.check_results(v2t_results, PIDOpt.Metrics(; Ms = 1.5, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
        println("Warning: Checks failed for v2t design point $point")
        println(v2t_results.metrics)
    end

    v2t_PID = PIDOpt.build_PID(v2t_results.params)
    C_v2t = named_ss(ss(v2t_PID), :C_v2t; u = :EAS_err, y = :throttle_cmd)
    v2t_sum = sumblock("EAS_err = EAS_dmd - EAS")
    v_θ_MIMO = connect([v2t_sum, C_v2t, thr_θ_MIMO], [:EAS_err=>:EAS_err, :EAS=>:EAS, :throttle_cmd=>:throttle_cmd], w1 = [:EAS_dmd, :θ_dmd], z1 = thr_θ_MIMO.y)

    ################################ c2θ #######################################

    P_θ2c = v_θ_MIMO[:c, :θ_dmd]
    t_sim_c2θ = 20
    lower_bounds = PIDParams(; k_p = 0.001, k_i = 0.001, k_d = 0.0, τ_f = 0.01)
    upper_bounds = PIDParams(; k_p = 0.03, k_i = 0.015, k_d = 0.0, τ_f = 0.01)
    initial_step = PIDParams(; k_p = 0.001, k_i = 0.001, k_d = 0.001, τ_f = 0.001)
    settings = PIDOpt.Settings(; t_sim = t_sim_c2θ, lower_bounds, upper_bounds, initial_step)
    weights = PIDOpt.Metrics(; Ms = 1.5, ∫e = 5.0, ef = 1.0, ∫u = 0.0, up = 0.1)
    params_0 = PIDParams(; k_p = 0.5, k_i = 0.1, k_d = 0.0, τ_f = 0.01)

    c2θ_results = PIDOpt.optimize_PID(P_θ2c; params_0, settings, weights, global_search = false)

    if !PIDOpt.check_results(c2θ_results, PIDOpt.Metrics(; Ms = 1.5, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
        println("Warning: Checks failed for c2θ design point $point")
        println(c2θ_results.metrics)
    end

    c2θ_PID = PIDOpt.build_PID(c2θ_results.params)
    C_c2θ = named_ss(ss(c2θ_PID), :C_c2θ; u = :c_err, y = :θ_dmd)
    c2θ_sum = sumblock("c_err = c_dmd - c")
    v_c_MIMO = connect([c2θ_sum, C_c2θ, v_θ_MIMO], [:c_err=>:c_err, :c=>:c, :θ_dmd =>:θ_dmd], w1 = [:EAS_dmd, :c_dmd], z1 = thr_θ_MIMO.y)


    ################################ v2θ #######################################
    P_θ2v = thr_θ_MIMO[:EAS, :θ_dmd]
    P_θ2v_opt = -P_θ2v #sign inversion on account of negative DC gain

    t_sim_v2θ = 20
    lower_bounds = PIDParams(; k_p = 0.01, k_i = 0.001, k_d = 0.0, τ_f = 0.01)
    upper_bounds = PIDParams(; k_p = 0.2, k_i = 0.05, k_d = 0.0, τ_f = 0.01)
    settings = PIDOpt.Settings(; t_sim = t_sim_v2θ, lower_bounds, upper_bounds)
    weights = PIDOpt.Metrics(; Ms = 2.0, ∫e = 5.0, ef = 1.0, ∫u = 0.0, up = 0.0)
    params_0 = PIDParams(; k_p = 0.05, k_i = 0.01, k_d = 0.0, τ_f = 0.01)

    v2θ_results = PIDOpt.optimize_PID(P_θ2v_opt; params_0, settings, weights, global_search)

    if !PIDOpt.check_results(v2θ_results, PIDOpt.Metrics(; Ms = 1.5, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
        println("Warning: Checks failed for v2θ design point $point")
        println(v2θ_results.metrics)
    end

    v2θ_PID = PIDOpt.build_PID(v2θ_results.params)

    return (q2e = q2e_results, θ2q = θ2q_results, v2t = v2t_results, c2θ = c2θ_results, v2θ = v2θ_results)

end

#optimizes PID parameters in the roll channel compensators
function optimize_roll(   ac::System{<:Cessna172FBWBase{NED}};
                        point::C172FBW.TrimParameters = C172FBW.TrimParameters(),
                        global_search::Bool = false)

    ail_rud_MIMO = named_ss(ac, point; model = :lat);

    P_a2p = ail_rud_MIMO[:p, :aileron_cmd]

    t_sim_p2a = 5
    upper_bounds = PIDParams(; k_p = 50.0, k_i = 20.0, k_d = 0.0, τ_f = 0.01)
    settings = PIDOpt.Settings(; t_sim = t_sim_p2a, upper_bounds)
    weights = PIDOpt.Metrics(; Ms = 1, ∫e = 10, ef = 1, ∫u = 0, up = 0.1)
    params_0 = PIDParams(; k_p = 1, k_i = 5, k_d = 0.01, τ_f = 0.01)

    p2a_results = PIDOpt.optimize_PID(P_a2p; params_0, settings, weights, global_search)

    if !PIDOpt.check_results(p2a_results, PIDOpt.Metrics(; Ms = 1.4, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
        println("Warning: Checks failed for p2a design point $point")
        println(p2a_results.metrics)
    end

    p2a_PID = PIDOpt.build_PID(p2a_results.params)
    C_p2a = named_ss(p2a_PID, :C_p2a; u = :p_err, y = :aileron_cmd);
    p2a_sum = sumblock("p_err = p_dmd - p")
    p_rud_MIMO = connect([p2a_sum, C_p2a, ail_rud_MIMO], [:p_err=>:p_err, :p=>:p, :aileron_cmd=>:aileron_cmd], w1 = [:p_dmd, :rudder_cmd], z1 = ail_rud_MIMO.y)

    P_p2φ = p_rud_MIMO[:φ, :p_dmd]

    t_sim_φ2p = 5
    lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
    upper_bounds = PIDParams(; k_p = 50.0, k_i = 10.0, k_d = 0.5, τ_f = 0.01)
    settings = PIDOpt.Settings(; t_sim = t_sim_φ2p, lower_bounds, upper_bounds)
    weights = PIDOpt.Metrics(; Ms = 2, ∫e = 10, ef = 1, ∫u = 0.00, up = 0.1)
    params_0 = PIDParams(; k_p = 2., k_i = 0., k_d = 0.0, τ_f = 0.01)

    φ2p_results = PIDOpt.optimize_PID(P_p2φ; params_0, settings, weights, global_search = false)

    if !PIDOpt.check_results(φ2p_results, PIDOpt.Metrics(; Ms = 1.4, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
        println("Warning: Checks failed for φ2p design point $point")
        println(φ2p_results.metrics)
    end

    φ2p_PID = PIDOpt.build_PID(φ2p_results.params)
    C_φ2p = named_ss(φ2p_PID, :φcmp; u = :φ_err, y = :p_dmd);
    φ2p_sum = sumblock("φ_err = φ_dmd - φ")
    φ_rud_MIMO = connect([φ2p_sum, C_φ2p, p_rud_MIMO], [:φ_err=>:φ_err, :φ=>:φ, :p_dmd=>:p_dmd], w1 = [:φ_dmd, :rudder_cmd], z1 = p_rud_MIMO.y)

    P_φ2χ = φ_rud_MIMO[:χ, :φ_dmd]

    t_sim_χ2φ = 30
    lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.3, k_d = 0.0, τ_f = 0.01)
    upper_bounds = PIDParams(; k_p = 10.0, k_i = 5.0, k_d = 0.5, τ_f = 0.01)
    settings = PIDOpt.Settings(; t_sim = t_sim_χ2φ, lower_bounds, upper_bounds)
    weights = PIDOpt.Metrics(; Ms = 2, ∫e = 10, ef = 1, ∫u = 0.00, up = 0.1)
    params_0 = PIDParams(; k_p = 3., k_i = 0., k_d = 0.0, τ_f = 0.01)

    χ2φ_results = PIDOpt.optimize_PID(P_φ2χ; params_0, settings, weights, global_search = false)

    if !PIDOpt.check_results(χ2φ_results, PIDOpt.Metrics(; Ms = 1.4, ∫e = 0.15, ef = 0.02, ∫u = Inf, up = Inf))
        println("Warning: Checks failed for χ2φ design point $point")
        println(χ2φ_results.metrics)
    end

    χ2φ_PID = PIDOpt.build_PID(χ2φ_results.params)
    C_χ2φ = named_ss(χ2φ_PID, :C_χ2φ; u = :χ_err, y = :p_dmd);

    return (p2a = p2a_results, φ2p = φ2p_results, χ2ϕ = χ2φ_results)

end



end #module