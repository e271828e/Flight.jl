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
    global_search::Bool = false,
    folder::String = dirname(@__DIR__)) #save to parent folder by default

    ac = Cessna172FBWBase(NED()) |> System

    mid_cg_pld = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

    results = map(Iterators.product(EAS_range, h_range)) do (EAS, h)

        println("Optimizing for EAS = $EAS, h = $h")

        flaps = EAS < 30 ? 1.0 : 0.0

        point = C172FBW.TrimParameters(; Ob = Geographic(LatLon(), HOrth(h)),
            EAS, flaps, γ_wOb_n = 0.0, x_fuel = 0.5, payload = mid_cg_pld)

        q_results = optimize_q(ac; point, global_search)
        p_results, β_results = optimize_pβ(ac; point, global_search)
        return q_results, p_results, β_results

    end

    q_results, p_results, β_results = StructArray(results) |> StructArrays.components

    EAS_bounds = (EAS_range[1], EAS_range[end])
    h_bounds = (h_range[1], h_range[end])

    map((q_results, p_results, β_results),
                ("q_lookup.h5", "p_lookup.h5", "β_lookup.h5")) do results, fname

        # PIDOpt.check_results.(results)
        data = PIDParams(StructArrays.components(StructArray(StructArray(results).params))...)
        lookup = C172FBWCAS.Lookup(data, EAS_bounds, h_bounds)
        C172FBWCAS.save_lookup(lookup, joinpath(folder, fname))
        return lookup

    end

end

#optimizes the PID in the pitch rate compensator for a single design point
function optimize_q(ac::System{<:Cessna172FBWBase{NED}};
                    point::C172FBW.TrimParameters = C172FBW.TrimParameters(),
                    global_search::Bool = false)

    thr_ele_MIMO = named_ss(ac, point; model = :lon);

    P_e2q = thr_ele_MIMO[:q, :elevator_cmd]
    q_int = tf(1, [1, 0]) |> ss
    P_q_opt = series(q_int, ss(P_e2q))

    settings = PIDOpt.Settings(; t_sim = 10,
        lower_bounds = PIDParams(; k_p = 0.0, k_i = 0.0, k_d = 0.0, τ_f = 0.01),
        upper_bounds = PIDParams(; k_p = 50.0, k_i = 100.0, k_d = 10.0, τ_f = 0.01))
    weights = PIDOpt.Metrics(; Ms = 1, ∫e = 10, ef = 2)
    params_0 = PIDParams(; k_p = 5, k_i = 30, k_d = 0.5, τ_f = 0.01)

    q_results = PIDOpt.optimize_PID(P_q_opt; params_0, settings, weights, global_search)

    return q_results

    # q_PID = PIDOpt.build_PID(q_results.params)
    # C_q2e = named_ss(series(q_int, q_PID), :qcmp; u = :q_err, y = :elevator_cmd);
    # qsum = sumblock("q_err = q_dmd - q")
    # thr_q_MIMO = connect([qsum, C_q2e, thr_ele_MIMO], [:q_err=>:q_err, :q=>:q, :elevator_cmd=>:elevator_cmd], w1 = [:throttle_cmd, :q_dmd], z1 = thr_ele_MIMO.y)

    # P_q2θ = thr_q_MIMO[:θ, :q_dmd]

    # settings = PIDOpt.Settings(; t_sim = 10,
    #     upper_bounds = PIDParams(; k_p = 50.0, k_i = 0.0, k_d = 5.0, τ_f = 0.01)) #disallow integral gain
    # weights = PIDOpt.Metrics(; Ms = 1, ∫e = 10, ef = 1)
    # params_0 = PIDParams(; k_p = 2.0, k_i = 0.0, k_d = 0.0, τ_f = 0.01)

    # θ_results = PIDOpt.optimize_PID(P_q2θ; params_0, settings, weights, global_search)

    # θ_PID = PIDOpt.build_PID(θ_results.params)
    # C_θ2q = named_ss(θ_PID, :θcmp; u = :θ_err, y = :q_dmd)

    # return (q_results, θ_results)

end

#optimizes the roll rate and sideslip PIDs for a single design point
function optimize_pβ(   ac::System{<:Cessna172FBWBase{NED}};
                        point::C172FBW.TrimParameters = C172FBW.TrimParameters(),
                        global_search::Bool = false)

    ail_rud_MIMO = named_ss(ac, point; model = :lat);

    P_a2p = ail_rud_MIMO[:p, :aileron_cmd]

    settings = PIDOpt.Settings(; t_sim = 5)
    weights = PIDOpt.Metrics(; Ms = 1, ∫e = 1, ef = 1)
    params_0 = PIDParams(; k_p = 1, k_i = 15, k_d = 0.025, τ_f = 0.01)

    p_results = PIDOpt.optimize_PID(P_a2p; params_0, settings, weights, global_search)

    p_PID = PIDOpt.build_PID(p_results.params)
    C_p2a = named_ss(p_PID, :pcmp; u = :p_err, y = :aileron_cmd);
    psum = sumblock("p_err = p_dmd - p")
    p_rud_MIMO = connect([psum, C_p2a, ail_rud_MIMO], [:p_err=>:p_err, :p=>:p, :aileron_cmd=>:aileron_cmd], w1 = [:p_dmd, :rudder_cmd], z1 = ail_rud_MIMO.y)

    P_r2β = p_rud_MIMO[:β, :rudder_cmd]
    P_β_opt = -P_r2β

    settings = PIDOpt.Settings(; t_sim = 5,
        upper_bounds = PIDParams(; k_p = 50.0, k_i = 50.0, k_d = 20.0, τ_f = 0.01))
    weights = PIDOpt.Metrics(; Ms = 1, ∫e = 10, ef = 2)
    params_0 = PIDParams(; k_p = 5.0, k_i = 20.0, k_d = 5.0, τ_f = 0.01)

    β_results = PIDOpt.optimize_PID(P_β_opt; params_0, settings, weights, global_search)
    return (p_results, β_results)

    # β_PID = PIDOpt.build_PID(β_results.params)
    # C_β2r = named_ss(-β_PID, :βcmp; u = :β_err, y = :rudder_cmd)
    # βsum = sumblock("β_err = β_dmd - β")
    # p_β_MIMO = connect([βsum, C_β2r, p_rud_MIMO], [:β_err=>:β_err, :β=>:β, :rudder_cmd=>:rudder_cmd], w1 = [:p_dmd, :β_dmd], z1 = p_rud_MIMO.y)

    # P_p2φ = p_β_MIMO[:φ, :p_dmd]
    # settings = PIDOpt.Settings(; t_sim = 5,
    #     upper_bounds = PIDParams(; k_p = 50.0, k_i = 0.0, k_d = 20.0, τ_f = 0.01)) #disallow integral gain
    # weights = PIDOpt.Metrics(; Ms = 1, ∫e = 10, ef = 2)
    # params_0 = PIDParams()

    # φ_results = PIDOpt.optimize_PID(P_p2φ; params_0, settings, weights, global_search)

    # return (p_results, β_results, φ_results)

end



end #module