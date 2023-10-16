module FBWLon

using Flight
using Flight.FlightCore.Systems
using Flight.FlightCore.Plotting

using Flight.FlightPhysics.Geodesy
using Flight.FlightPhysics.Kinematics
using Flight.FlightComponents.Aircraft
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172

using ControlSystems
using RobustAndOptimalControl

include("pid_design.jl"); using .PIDDesign

struct PitchRateLookupData
    coefs::PIDDesign.Params{Array{Float64,3}}
    TAS_bounds::NTuple{2, Float64}
    Mt_bounds::NTuple{2, Float64}
end

# function LookupData(TAS_range::AbstractRange{Float64} = range(30, 60, length = 2),
#                     h_range::AbstractRange{Float64} = range(10, 4010, length = 2))
function LookupData(TAS_range::AbstractRange{Float64} = range(60, 60, length = 1),
                    h_range::AbstractRange{Float64} = range(10, 4010, length = 2))

    ac = Cessna172FBWBase(NED()) |> System

    # fwd_cg_pld = C172.PayloadU(m_pilot = 100, m_copilot = 100, m_baggage = 0)
    # aft_cg_pld = C172.PayloadU(m_pilot = 50, m_copilot = 50, m_baggage = 100)
    mid_cg_pld = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

    data_size = (length(TAS_range), length(h_range))
    metrics = fill(PIDDesign.Metrics(zeros(3)), data_size)
    params = fill(PIDDesign.Params(zeros(4)), data_size)

    params_0 = PIDDesign.Params(; k_p = 4, k_i = 30, k_d = 0.43, τ_f = 0.01)

    for (i, (TAS, h)) in enumerate(Iterators.product(TAS_range, h_range))

        point = C172FBW.TrimParameters(; Ob = Geographic(LatLon(), HOrth(h)),
            TAS, γ_wOb_n = 0.0, x_fuel = 0.5, flaps = 0.0, payload = mid_cg_pld)

        println("Optimization for TAS = $TAS, h = $h")
        metrics[i], params[i] = optimize_point_q(ac; point, params_0, global_search = true)
        # exit_flag[i] = :TEST_FLAG
        # metrics[i] = PIDDesign.Metrics(zeros(3))
        # params[i] = PIDDesign.Params(; k_p = TAS, k_i = h)

        # return (exit_flag, metrics_opt, PID_params_opt)
        # return (metrics_opt, PID_params_opt)

    end

    return metrics, params

    # coefs = Coefficients(StructArrays.components(StructArray(data_points))...)

    # LookupData(coefs, J_bounds, Mt_bounds, Δβ_bounds)

    TAS_bounds = (TAS_range[1], TAS_range[end])
    h_bounds = (h_range[1], h_range[end])

end

#optimizes the PID in the pitch rate compensator for a single design point
function optimize_point_q(  ac::System{<:Cessna172FBWBase{NED}};
                            point::C172FBW.TrimParameters = C172FBW.TrimParameters(),
                            params_0::PIDDesign.Params = PIDDesign.Params(),
                            global_search::Bool = true)

    thr_ele_MIMO = named_ss(ac, point; model = :lon);
    P_e2q = thr_ele_MIMO[:q, :elevator_cmd]
    q_int = tf(1, [1, 0]) |> ss
    P_q_opt = series(q_int, ss(P_e2q))

    settings = PIDDesign.Settings(; t_sim = 10, maxeval = 5000)
    weights = PIDDesign.Metrics(; Ms = 1, ∫e = 10, ef = 2)

    exit_flag, cost, metrics_opt, q_PID_params_opt, q_PID_opt = PIDDesign.optimize_PID(
        P_q_opt; params_0, settings, weights, global_search)
    println(exit_flag)

    # return exit_flag, metrics_opt, q_PID_params_opt
    return metrics_opt, q_PID_params_opt

end

end #module