module Trim

using LinearAlgebra
using StaticArrays
using ComponentArrays
using BenchmarkTools
using UnPack
# using Optim
using NLopt

using Flight.Systems
using Flight.Attitude
using Flight.Geodesy
using Flight.Kinematics
using Flight.Air
using Flight.Terrain
using Flight.Piston
using Flight.Aircraft

import Flight.Kinematics: KinematicInit
export XTrimTemplate, XTrim, TrimParameters

using ..C172R

############################### Trimming #######################################
################################################################################

#for trimming, we set the yoke incremental commands to zero. the reason is that,
#when we start an interactive simulation from a trimmed state, leaving the
#controller sticks at neutral corresponds to zero yoke incremental commands. in
#this case, it is the yoke_x and yoke_y values that set the surfaces where they
#need to be


# maybe split trim state into common and aircraft dependent, that way the first
# ones can be reused
const StateTemplate = ComponentVector(
    α_a = 0.0, #angle of attack, aerodynamic axes
    φ_nb = 0.0, #bank angle
    n_eng = 0.5, #normalized engine speed (ω/ω_rated)
    throttle = 0.5,
    yoke_x = 0.0,
    yoke_y = 0.0,
    pedals = 0.0,
)

const State{T, D} = ComponentVector{T, D, typeof(getaxes(StateTemplate))} where {T, D}

function State(; α_a = 0.0848, φ_nb = 0.0, n_eng = 0.75,
    throttle = 0.62, yoke_x = 0.015, yoke_y = -0.006, pedals = -0.03)

    x = copy(StateTemplate)
    @pack! x = α_a, φ_nb, n_eng, throttle, yoke_x, yoke_y, pedals
    return x

end

#the first 5 do not depend on the aircraft type
struct Parameters
    Ob::GeographicLocation{NVector, Ellipsoidal} #position
    ψ_nb::Float64 #geographic heading
    TAS::Float64 #true airspeed
    γ_wOb_n::Float64 #wind-relative flight path angle
    ψ_lb_dot::Float64 #LTF-relative turn rate
    β_a::Float64 #sideslip angle measured in the aerodynamic reference frame
    fuel::Float64 #fuel load, 0 to 1
    mixture::Float64 #engine mixture control, 0 to 1
    flaps::Float64 #flap setting, 0 to 1
end

function Parameters(;
    l2d::Abstract2DLocation = LatLon(), h::Altitude = AltO(1000),
    ψ_nb = 0.0, TAS = 40.0, γ_wOb_n = 0.0, ψ_lb_dot = 0.0, β_a = 0.0,
    fuel = 0.5, mixture = 0.5, flaps = 0.0)

    Ob = GeographicLocation(l2d, h)
    Parameters(Ob, ψ_nb, TAS, γ_wOb_n, ψ_lb_dot, β_a, fuel, mixture, flaps)
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

function KinematicInit(state::State, params::Parameters, atm::System{<:Atmosphere})

    v_wOb_a = Air.get_velocity_vector(params.TAS, state.α_a, params.β_a)
    v_wOb_b = C172R.f_ba.q(v_wOb_a) #wind-relative aircraft velocity, body frame

    θ_nb = θ_constraint(; v_wOb_b, params.γ_wOb_n, state.φ_nb)
    e_nb = REuler(params.ψ_nb, θ_nb, state.φ_nb)
    q_nb = RQuat(e_nb)

    e_lb = e_nb #initialize LTF arbitrarily to NED
    ė_lb = SVector(params.ψ_lb_dot, 0.0, 0.0)
    ω_lb_b = Attitude.ω(e_lb, ė_lb)

    l2d = NVector(params.Ob)
    h = AltE(params.Ob)

    v_wOb_n = q_nb(v_wOb_b) #wind-relative aircraft velocity, NED frame
    v_ew_n = AtmosphericData(atm, params.Ob).wind.v_ew_n
    v_eOb_n = v_ew_n + v_wOb_n

    KinematicInit(; q_nb, l2d, h, ω_lb_b, v_eOb_n)

end

#assigns trim state and parameters to the aircraft system, and then updates it
#by calling its continuous dynamics function
function assign!(ac::System{<:Cessna172R}, atm::System{<:Atmosphere},
    trn::AbstractTerrain, state::State, params::Parameters)

    Aircraft.init!(ac, KinematicInit(state, params, atm))

    ω_eng = state.n_eng * ac.airframe.pwp.engine.params.ω_rated

    ac.x.airframe.aero.α_filt = state.α_a #ensures zero state derivative
    ac.x.airframe.aero.β_filt = params.β_a #ensures zero state derivative
    ac.x.airframe.pwp.ω = ω_eng
    ac.x.airframe.fuel .= params.fuel

    ac.u.avionics.throttle = state.throttle
    ac.u.avionics.yoke_x = state.yoke_x
    ac.u.avionics.yoke_y = state.yoke_y
    ac.u.avionics.pedals = state.pedals
    ac.u.avionics.flaps = params.flaps
    ac.u.avionics.mixture = params.mixture

    #incremental yoke positions to zero
    ac.u.avionics.yoke_Δx = 0
    ac.u.avionics.yoke_Δy = 0

    #engine must be running, no way to trim otherwise
    ac.d.airframe.pwp.engine.state = Piston.eng_running
    @assert ac.x.airframe.pwp.ω > ac.airframe.pwp.engine.idle.params.ω_target

    #as long as the engine remains above the idle controller's target speed, the
    #idle controller's output will be saturated at 0 by proportional error, so
    #the integrator will be disabled and its state will not change. we just set
    #it to zero and forget about it
    ac.x.airframe.pwp.engine.idle .= 0.0

    #powerplant friction regulator: with the propeller spinning, the friction
    #regulator's output will be saturated at -1 due to proportional error, so
    #the integrator will be disabled and its state will not change. we just set
    #it to zero and forget about it
    ac.x.airframe.pwp.friction .= 0.0

    f_cont!(ac, atm, trn)

    #check assumptions concerning airframe systems states & derivatives
    wow = reduce(|, SVector{3,Bool}(leg.strut.wow for leg in ac.y.airframe.ldg))
    @assert wow === false
    @assert ac.ẋ.airframe.pwp.engine.idle[1] .== 0
    @assert ac.ẋ.airframe.pwp.friction[1] .== 0

end

function cost(ac::System{<:Cessna172R})


    # necesitamos introducir la barrier functions ANTES de que los control
    # inputs saturen a 1, porque entonces la cost function se hace estrictamente
    # infinita. por tanto, no podemos coger los valores de los controles de
    # dentro de la avionica para despues penalizarlos. hay que cogerlos
    # directamente del trim_state ANTES de que se asignen al sistema y se
    # saturen. o sea que o bien el interfaz de la cost function anade el state
    # como input, y se mantiene la funcion assign! independiente, o redefinimos
    # la funcion cost como cost! y absorbemos dentro assign!



    v_nd_dot = SVector{3}(ac.ẋ.kinematics.vel.v_eOb_b) / norm(ac.y.kinematics.common.v_eOb_b)
    ω_dot = SVector{3}(ac.ẋ.kinematics.vel.ω_eb_b) #ω should already of order 1
    n_eng_dot = ac.ẋ.airframe.pwp.ω / ac.airframe.pwp.engine.params.ω_rated

    sum(v_nd_dot.^2) + sum(ω_dot.^2) + n_eng_dot^2

end

function get_target_function(ac::System{<:Cessna172R},
    atm::System{<:Atmosphere}, trn::AbstractTerrain,
    params::Parameters = Parameters())

    f_target = let ac = ac, atm = atm, trn = trn, params = params

        function (x::State)
            assign!(ac, atm, trn, x, params)
            return cost(ac)
        end

    end

    return f_target

end

# function trim!(ac::System{<:Cessna172R},
#     atm::System{<:Atmosphere}, trn::AbstractTerrain,
#     x0::State = State(), params::Parameters = Parameters();
#     a = 0.02, b = 0.05, g_tol = 1e-14, iterations = 5000)

#     f = get_target_function(ac, atm, trn, params)
#     initial_simplex = Optim.AffineSimplexer(; a, b)
#     @show initial_cost = f(x0)
#     result = optimize(f, x0, NelderMead(; initial_simplex), Optim.Options(; g_tol, iterations))
#     @show final_cost = f(result.minimizer)
#     if final_cost > 10g_tol
#         println("Warning: Optimization did not converge")
#     end
#     @show
#     return result
# end

function trim_new!(; ac::System{<:Cessna172R} = System(Cessna172R()),
    atm::System{<:Atmosphere} = System(Atmosphere()),
    trn::AbstractTerrain = HorizontalTerrain(),
    state::State = State(), params::Parameters = Parameters())

    f_target = get_target_function(ac, atm, trn, params)

    #wrapper function with the interface required by NLopt
    ax = getaxes(state)
    function f_opt(x::Vector{Float64}, ::Vector{Float64})
        s = ComponentVector(x, ax)
        return f_target(s)
    end

    n = length(state)
    x0 = zeros(n); lower_bounds = similar(x0); upper_bounds = similar(x0)

    x0[:] .= state

    lower_bounds[:] .= State(
        α_a = -π/6,
        φ_nb = -π/3,
        n_eng = 0,
        throttle = 0,
        yoke_x = -1,
        yoke_y = -1,
        pedals = -1)

    upper_bounds[:] .= State(
        α_a = π/6,
        φ_nb = π/3,
        n_eng = 1.2,
        throttle = 1,
        yoke_x = 1,
        yoke_y = 1,
        pedals = 1)

    # @show initial_cost = f_opt(x0, x0)
    # @btime $f_opt($x0, $x0)

    opt  = Opt(:LN_NELDERMEAD, length(x0))
    opt.min_objective = f_opt
    opt.maxeval = 10000
    opt.stopval = 1e-12

    println("TUNE opt.initial_step FOR THE DIFFERENT trim states and retry COBYLA")
    opt.lower_bounds = lower_bounds
    opt.upper_bounds = upper_bounds

    # @show x0
    # @show opt.lower_bounds
    # @show opt.upper_bounds

    # @btime $f($x0, $x0)

    (minf,minx,ret) = optimize(opt, x0)
    @show ret
    @show minf
    @show numevals = opt.numevals # the number of function evaluations
    # println("got $minf at $minx after $numevals iterations (returned $ret)")

    @show final_state = ComponentVector(minx, ax)
    # @btime optimize($opt, $x0)


end





end #module