module Trim

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack

using Flight.Systems
using Flight.Attitude
using Flight.Geodesy
using Flight.Kinematics
using Flight.Air
using Flight.Terrain

import Flight.Kinematics: KinematicInit
export XTrimTemplate, XTrim, TrimParameters

using ..C172R

############################### Trimming #######################################
################################################################################

const TrimStateTemplate = ComponentVector(
    α = 0.08, θ_nb = 0.08, φ_nb = 0.0, ω_eng = 215.0, ω_lb_b = [0, 0, 0],
    throttle = 0.61, yoke_Δx = 0.01, yoke_Δy = -0.025, pedals = 0.0,
)

const TrimState{T, D} = ComponentVector{T, D, typeof(getaxes(TrimStateTemplate))} where {T, D}
TrimState() = similar(TrimStateTemplate)

#these do not require evaluation of are assigned to the aircraft states and inputs
struct TrimParameters
    Ob::GeographicLocation{NVector, Ellipsoidal}
    TAS::Float64
    β::Float64
    ψ_nb::Float64 #geographic heading
    yoke_x0::Float64 #force-free (aileron-trimmed) yoke position
    yoke_y0::Float64 #force-free (elevator-trimmed) yoke position
    fuel::Float64 #fuel load, 0 to 1
    mixture::Float64 #engine mixture control, 0 to 1
    flaps::Float64 #flap setting, 0 to 1
end

function TrimParameters(;
    l2d::Abstract2DLocation = LatLon(), h::Altitude = AltO(1000),
    TAS = 40.0, β = 0.0, ψ_nb = 0.0, yoke_x0 = 0.0, yoke_y0 = 0.0,
    fuel = 0.5, mixture = 0.5, flaps = 0.0)

    Ob = GeographicLocation(l2d, h)
    TrimParameters(Ob, TAS, β, ψ_nb, yoke_x0, yoke_y0, fuel, mixture, flaps)
end

Base.@kwdef struct TrimConstraints
    h_dot::Float64 = 0.0 #climb rate
    ψ_nb_dot::Float64 = 0.0 #turn rate
end

function KinematicInit(state::TrimState, params::TrimParameters, atm::Atmosphere)

    q_nb = REuler(params.ψ_nb, state.θ_nb, state.φ_nb) |> RQuat
    l2d = NVector(params.Ob)
    h = AltE(params.Ob)
    ω_lb_b = state.ω_lb_b

    v_wOb_a = Air.get_velocity_vector(params.TAS, state.α, params.β)
    v_wOb_b = C172R.f_ba.q(v_wOb_a) #wind-relative aircraft velocity, body frame
    v_wOb_n = q_nb(v_wOb_b) #wind-relative aircraft velocity, NED frame
    v_ew_n = AtmosphericData(atm, params.Ob).wind.v_ew_n
    v_eOb_n = v_ew_n + v_wOb_n

    KinematicInit(; q_nb, l2d, h, ω_lb_b, v_eOb_n)

end

#assigns trim state and parameters to the aircraft system, and then updates it
#by calling its continuous dynamics function
function assign!(ac::System{<:Cessna172R}, atm::System{<:Atmosphere},
    trn::AbstractTerrain, state::TrimState, params::TrimParameters)

    init!(ac, KinematicInit(state, params, atm))

    ac.x.airframe.aero.α_filt = state.α #ensures zero state derivative
    ac.x.airframe.aero.β_filt = params.β #ensures zero state derivative
    ac.x.airframe.pwp.ω = state.ω_eng
    ac.x.airframe.fuel .= params.fuel

    ac.u.avionics.throttle = state.throttle
    ac.u.avionics.yoke_Δx = state.yoke_Δx
    ac.u.avionics.yoke_Δy = state.yoke_Δy
    ac.u.avionics.pedals = state.pedals
    ac.u.avionics.flaps = params.flaps
    ac.u.avionics.yoke_x0 = params.yoke_x0
    ac.u.avionics.yoke_y0 = params.yoke_y0
    ac.u.avionics.mixture = params.mixture

    #engine must be running, no way to trim otherwise
    ac.d.airframe.pwp.engine.state = Piston.eng_running
    @assert ac.x.airframe.pwp.ω > ac.airframe.pwp.engine.idle.params.ω_idle

    #as long as the engine remains above the idle controller's target speed, the
    #idle controller's output will be saturated at 0 by proportional error, so
    #the integrator will be disabled and its state will not change. we just set
    #it to zero and forget about it
    ac.x.airframe.pwp.engine.idle = 0.0

    #powerplant friction regulator: with the propeller spinning, the friction
    #regulator's output will be saturated at -1 due to proportional error, so
    #the integrator will be disabled and its state will not change. we just set
    #it to zero and forget about it
    ac.x.airframe.pwp.friction .= 0.0

    f_cont!(ac, atm, trn)

    println("Aqui habria que comprobar varias cosas: que WoW es cero, que x_idle y x_friction son cero, y que alpha filt y beta filt son cero")

end

function cost(ac::System{<:Cessna172R{NED}}, constraints::TrimConstraints)

    #this function assumes the trim_state and parameters have already been
    #assigned and the System evaluated; the cost can be calculated entirely from
    #the aircraft system and the user-specified constraints

    return

end


function trim!(ac::System{<:Cessna172R{NED}}, atm::System{<:Atmosphere}, trn::AbstractTerrain,
    state0::TrimState = TrimState(), params::TrimParameters = TrimParameters(),
    constraints::TrimConstraints = TrimConstraints(); kwargs...)

    f_target = let ac = ac, params = params, constraints = constraints, atm = atm, trn = trn

        function (x::TrimState)
            assign!(ac, x, params, atm, trn)
            return cost(ac, constraints)
        end

    end

end


function test_target_eval()

    aircraft = System(Cessna172R())
    atmosphere = System(Atmosphere())
    terrain = HorizontalTerrain() #zero orthometric altitude


    f = get_target_function(aircraft, atmosphere, terrain, trim_parameters, trim_configuration)


    f(x_trim_0)
    # println(@ballocated($f($x_trim_0))) #no allocations! nice!

end


function trim!(ac, atmosphere, terrain, trim_parameters, trim_configuration)


    #is this let really needed? or being inside a function already freezes values
    func = let  x = aircraft.x, d = aircraft.d, u = aircraft.u, atmosphere = atmosphere, terrain = terrain,
                trim_parameters = trim_parameters, trim_configuration = trim_configuration

        #set the engine running, no way to trim otherwise
        d.airframe.pwp.engine.state = Piston.eng_running

        u.avionics.flaps = trim_configuration.flaps
        u.avionics.yoke_x0 = trim_configuration.yoke_x0
        u.avionics.yoke_y0 = trim_configuration.yoke_y0
        u.airframe.pwp.engine.mix = trim_configuration.mixture
        x.airframe.fuel .= trim_configuration.fuel

        x.airframe.pwp.engine.idle .= trim_configuration.x_pwp_idle
        x.airframe.pwp.friction .= trim_configuration.x_pwp_friction

        function compute_derivative(x_trim)


            x.airframe.aero.α_filt = x_trim.α #this ensures zero state derivative
            x.airframe.aero.β_filt = trim_parameters.β #this ensures zero state derivative
            x.airframe.pwp.ω = x_trim.ω_eng

            u.avionics.throttle = x_trim.throttle
            u.avionics.yoke_Δx = x_trim.yoke_Δx
            u.avionics.yoke_Δy = x_trim.yoke_Δy
            u.avionics.pedals = x_trim.pedals

            f_cont!(aircraft, atmosphere, terrain)

            println(aircraft.x)
            println()
            # println("cuando acabe con trim, create XPosMinimal, validate against others") #minimal attitude and pos repr
            println("ensure x_idle and x_friction have indeed zero derivatives")
            println("do not forget psi_dot in the cost function")

            return aircraft.ẋ

        end

    end

    return func

end

end #module