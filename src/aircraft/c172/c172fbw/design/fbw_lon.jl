module FBWLon

using Flight
using Flight.FlightCore.Systems
using Flight.FlightCore.Plotting

using Flight.FlightPhysics.Geodesy
using Flight.FlightPhysics.Kinematics
using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.Control
using Flight.FlightComponents.Aircraft
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172

using UnPack
using ControlSystems
using RobustAndOptimalControl
using StructArrays


struct LookupData
    params::PIDOpt.Params{Array{Float64,2}}
    EAS_bounds::NTuple{2, Float64}
    Mt_bounds::NTuple{2, Float64}
end

# struct PitchLookupData
#     q::LookupData
#     θ::LookupData
# end

function generate_q_data(EAS_range::AbstractRange{Float64} = range(25, 55, length = 4),
                        h_range::AbstractRange{Float64} = range(10, 3010, length = 4))

    ac = Cessna172FBWBase(NED()) |> System

    # fwd_cg_pld = C172.PayloadU(m_pilot = 100, m_copilot = 100, m_baggage = 0)
    # aft_cg_pld = C172.PayloadU(m_pilot = 50, m_copilot = 50, m_baggage = 100)
    mid_cg_pld = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

    #set reasonable initial condition for quicker convergence
    params_0 = PIDOpt.Params(; k_p = 5, k_i = 30, k_d = 0.5, τ_f = 0.01)

    results = map(Iterators.product(EAS_range, h_range)) do (EAS, h)

        println("Optimizing for EAS = $EAS, h = $h")
        flaps = EAS < 30 ? 1.0 : 0.0

        point = C172FBW.TrimParameters(; Ob = Geographic(LatLon(), HOrth(h)),
            EAS, flaps, γ_wOb_n = 0.0, x_fuel = 0.5, payload = mid_cg_pld)

        optimize_point(ac; point, params_0, global_search = false)

    end

    q_results, θ_results = StructArray(results) |> StructArrays.components

    PIDOpt.check_results.(q_results)
    # check_results.(θ_results)

    q_params = PIDOpt.Params(StructArrays.components(StructArray(StructArray(q_results).params))...)
    # θ_params = PIDOpt.Params(StructArrays.components(StructArray(StructArray(θ_results).params))...)

    EAS_bounds = (EAS_range[1], EAS_range[end])
    h_bounds = (h_range[1], h_range[end])

    q_data = LookupData(q_params, EAS_bounds, h_bounds)

    return q_data

end

#optimizes the PID in the pitch rate compensator for a single design point
function optimize_point(  ac::System{<:Cessna172FBWBase{NED}};
                            point::C172FBW.TrimParameters = C172FBW.TrimParameters(),
                            params_0::PIDOpt.Params = PIDOpt.Params(),
                            global_search::Bool = true)

    thr_ele_MIMO = named_ss(ac, point; model = :lon);
    P_e2q = thr_ele_MIMO[:q, :elevator_cmd]
    q_int = tf(1, [1, 0]) |> ss
    P_q_opt = series(q_int, ss(P_e2q))

    lower_bounds = PIDOpt.Params(; k_p = 0.0, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
    upper_bounds = PIDOpt.Params(; k_p = 50.0, k_i = 100.0, k_d = 10.0, τ_f = 0.01)
    settings = PIDOpt.Settings(; t_sim = 10, maxeval = 5000, lower_bounds, upper_bounds)
    weights = PIDOpt.Metrics(; Ms = 1, ∫e = 10, ef = 2)

    q_results = PIDOpt.optimize_PID(P_q_opt; params_0, settings, weights, global_search)
    return (q_results, :θ_results)

end


end #module