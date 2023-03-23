module Trim

using LinearAlgebra
using StaticArrays
using ComponentArrays
using BenchmarkTools
using UnPack
using NLopt

using Flight.FlightCore
using Flight.FlightPhysics

using Flight.FlightAircraft.Piston
using Flight.FlightAircraft.Aircraft

export XTrimTemplate, XTrim, TrimParameters

using ..C172R

############################### Trimming #######################################
################################################################################

#for trimming, the incremental control column inputs Δelevator and Δaileron are
#set to zero. when a simulation from a trimmed state, leaving the controller
#sticks at neutral positions corresponds to zero Δelevator and Δaileron, and
#therefore will retain the trim control inputs. it is the aileron and elevator
#input values that set the trim surface positions

#first 2 are aircraft-agnostic
const StateTemplate = ComponentVector(
    α_a = 0.0, #angle of attack, aerodynamic axes
    φ_nb = 0.0, #bank angle
    n_eng = 0.5, #normalized engine speed (ω/ω_rated)
    throttle = 0.5,
    aileron = 0.0,
    elevator = 0.0,
    rudder = 0.0, #rudder↑ -> aero.u.r↓ -> right yaw
)

const State{T, D} = ComponentVector{T, D, typeof(getaxes(StateTemplate))} where {T, D}

function State(; α_a = 0.0848, φ_nb = 0.0, n_eng = 0.75,
    throttle = 0.62, aileron = 0.015, elevator = -0.006, rudder = -0.03)

    x = copy(StateTemplate)
    @pack! x = α_a, φ_nb, n_eng, throttle, aileron, elevator, rudder
    return x

end

#first 5 are aircraft-agnostic
struct Parameters
    Ob::Geographic{NVector, Ellipsoidal} #position
    ψ_nb::Float64 #geographic heading
    TAS::Float64 #true airspeed
    γ_wOb_n::Float64 #wind-relative flight path angle
    ψ_lb_dot::Float64 #LTF-relative turn rate
    θ_lb_dot::Float64 #LTF-relative pitch rate
    β_a::Float64 #sideslip angle measured in the aerodynamic reference frame
    fuel::Float64 #fuel load, 0 to 1
    mixture::Float64 #engine mixture control, 0 to 1
    flaps::Float64 #flap setting, 0 to 1
end

function Parameters(;
    loc::Abstract2DLocation = LatLon(), h::Altitude = HOrth(1000),
    ψ_nb = 0.0, TAS = 40.0, γ_wOb_n = 0.0, ψ_lb_dot = 0.0, θ_lb_dot = 0.0,
    β_a = 0.0, fuel = 0.5, mixture = 0.5, flaps = 0.0)

    Ob = Geographic(loc, h)
    Parameters(Ob, ψ_nb, TAS, γ_wOb_n, ψ_lb_dot, θ_lb_dot, β_a, fuel, mixture, flaps)
end

#given the body-axes wind-relative velocity, the wind-relative flight path angle
#and the bank angle, the pitch angle is unambiguously determined
function θ_constraint(; v_wOb_b, γ_wOb_n, φ_nb)
    TAS = norm(v_wOb_b)
    a = v_wOb_b[1] / TAS
    b = (v_wOb_b[2] * sin(φ_nb) + v_wOb_b[3] * cos(φ_nb)) / TAS
    sγ = sin(γ_wOb_n)

    return atan((a*b + sγ*√(a^2 + b^2 - sγ^2))/(a^2 - sγ^2))
    # return asin((a*sγ + b*√(a^2 + b^2 - sγ^2))/(a^2 + b^2)) #equivalent

end

function Kinematics.Initializer(state::State, params::Parameters,
                                env::System{<:AbstractEnvironment})

    v_wOb_a = Atmosphere.get_velocity_vector(params.TAS, state.α_a, params.β_a)
    v_wOb_b = C172R.C172RAirframe.f_ba.q(v_wOb_a) #wind-relative aircraft velocity, body frame

    θ_nb = θ_constraint(; v_wOb_b, params.γ_wOb_n, state.φ_nb)
    e_nb = REuler(params.ψ_nb, θ_nb, state.φ_nb)
    q_nb = RQuat(e_nb)

    e_lb = e_nb #initialize LTF arbitrarily to NED
    ė_lb = SVector(params.ψ_lb_dot, params.θ_lb_dot, 0.0)
    ω_lb_b = Attitude.ω(e_lb, ė_lb)

    loc = NVector(params.Ob)
    h = HEllip(params.Ob)

    v_wOb_n = q_nb(v_wOb_b) #wind-relative aircraft velocity, NED frame
    v_ew_n = AtmosphericData(env.atm, params.Ob).wind.v_ew_n
    v_eOb_n = v_ew_n + v_wOb_n

    Kinematics.Initializer(; q_nb, loc, h, ω_lb_b, v_eOb_n, Δx = 0.0, Δy = 0.0)

end

#assigns trim state and parameters to the aircraft system, and then updates it
#by calling f_ode!
function assign!(ac::System{<:Cessna172R}, env::System{<:AbstractEnvironment},
    params::Parameters, state::State)

    init_kinematics!(ac, Kinematics.Initializer(state, params, env))

    ω_eng = state.n_eng * ac.airframe.pwp.engine.params.ω_rated

    ac.x.airframe.aero.α_filt = state.α_a #ensures zero state derivative
    ac.x.airframe.aero.β_filt = params.β_a #ensures zero state derivative
    ac.x.airframe.pwp.engine.ω = ω_eng
    ac.x.airframe.fuel .= params.fuel

    ac.u.avionics.throttle = state.throttle
    ac.u.avionics.aileron_trim = state.aileron
    ac.u.avionics.elevator_trim = state.elevator
    ac.u.avionics.rudder_trim = state.rudder
    ac.u.avionics.flaps = params.flaps
    ac.u.avionics.mixture = params.mixture

    #incremental control inputs to zero
    ac.u.avionics.elevator = 0
    ac.u.avionics.aileron = 0
    ac.u.avionics.rudder = 0

    #engine must be running, no way to trim otherwise
    ac.s.airframe.pwp.engine.state = Piston.eng_running
    @assert ac.x.airframe.pwp.engine.ω > ac.airframe.pwp.engine.params.ω_idle

    #engine idle compensator: as long as the engine remains at normal
    #operational speeds, well above its nominal idle speed, the idle controller
    #compensator's output will be saturated at its lower bound by proportional
    #error. its integrator will be disabled, its state will not change nor have
    #any effect on the engine. we can simply set it to zero
    ac.x.airframe.pwp.engine.idle .= 0.0

    #engine friction compensator: with the engine running at normal operational
    #speeds, the engine's friction constraint compensator will be saturated, so
    #its integrator will be disabled and its state will not change. furthermore,
    #with the engine running friction is ignored. we can simply set it to zero.
    ac.x.airframe.pwp.engine.frc .= 0.0

    f_ode!(ac, env)

    #check assumptions concerning airframe systems states & derivatives
    wow = reduce(|, SVector{3,Bool}(leg.strut.wow for leg in ac.y.airframe.ldg))
    @assert wow === false
    @assert ac.ẋ.airframe.pwp.engine.idle[1] .== 0
    @assert ac.ẋ.airframe.pwp.engine.frc[1] .== 0

end

function cost(ac::System{<:Cessna172R})

    v_nd_dot = SVector{3}(ac.ẋ.kinematics.vel.v_eOb_b) / norm(ac.y.kinematics.common.v_eOb_b)
    ω_dot = SVector{3}(ac.ẋ.kinematics.vel.ω_eb_b) #ω should already of order 1
    n_eng_dot = ac.ẋ.airframe.pwp.engine.ω / ac.airframe.pwp.engine.params.ω_rated

    sum(v_nd_dot.^2) + sum(ω_dot.^2) + n_eng_dot^2

end

function get_target_function(ac::System{<:Cessna172R},
    env::System{<:AbstractEnvironment}, params::Parameters = Parameters())

    let ac = ac, env = env, params = params
        function (x::State)
            assign!(ac, env, params, x)
            return cost(ac)
        end
    end

end

function trim!(ac::System{<:Cessna172R},
    env::System{<:AbstractEnvironment} = System(SimpleEnvironment()),
    params::Parameters = Parameters(),
    state::State = State()) #optional initial guess

    f_target = get_target_function(ac, env, params)

    #wrapper around f_target with the interface required by NLopt
    ax = getaxes(state)
    function f_opt(x::Vector{Float64}, ::Vector{Float64})
        s = ComponentVector(x, ax)
        return f_target(s)
    end

    n = length(state)
    x0 = zeros(n); lower_bounds = similar(x0); upper_bounds = similar(x0); initial_step = similar(x0)

    x0[:] .= state

    lower_bounds[:] .= State(
        α_a = -π/12,
        φ_nb = -π/3,
        n_eng = 0.4,
        throttle = 0,
        aileron = -1,
        elevator = -1,
        rudder = -1)

    upper_bounds[:] .= State(
        α_a = ac.airframe.aero.params.α_stall[2], #critical AoA is 0.28 < 0.36
        φ_nb = π/3,
        n_eng = 1.1,
        throttle = 1,
        aileron = 1,
        elevator = 1,
        rudder = 1)

    initial_step[:] .= 0.05 #safe value for all optimization variables

    # @show initial_cost = f_opt(x0, x0)
    # @btime $f_opt($x0, $x0)

    #any of these three algorithms works
    # opt = Opt(:LN_NELDERMEAD, length(x0))
    opt = Opt(:LN_BOBYQA, length(x0))
    # opt = Opt(:GN_CRS2_LM, length(x0))
    opt.min_objective = f_opt
    opt.maxeval = 100000
    opt.stopval = 1e-14
    opt.lower_bounds = lower_bounds
    opt.upper_bounds = upper_bounds
    opt.initial_step = initial_step

    # @btime optimize($opt, $x0)

    (minf, minx, exit_flag) = optimize(opt, x0)

    if exit_flag != :STOPVAL_REACHED
        println("Warning: Optimization did not converge")
    end
    state_opt = ComponentVector(minx, ax)
    assign!(ac, env, params, state_opt)
    return (exit_flag = exit_flag, result = state_opt)


end


end #module