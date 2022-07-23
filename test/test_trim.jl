module TestTrim

using Test
using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack

using Flight


#TrimState, variables to optimize
function test_target_eval()

    aircraft = System(C172RAircraft())
    atmosphere = System(Atmosphere())
    terrain = HorizontalTerrain() #zero orthometric altitude

    #TrimConfiguration
    trim_configuration = (
        fuel = 0.5, #fuel load, 0 to 1
        mixture = 0.5, #engine mixture control, 0 to 1
        flaps = 0, #flap setting, 0 to 1
        yoke_x0 = 0.0, #force-free (aileron-trimmed) yoke position
        yoke_y0 = 0.0, #force-free (elevator-trimmed) yoke position
        #engine idle controller state: with the engine throttle above idle, the
        #controller output will be saturated at 0 due to proportional error, so the
        #integrator will be disabled and the state will not change
        x_pwp_idle = 0.0,
        #powerplant friction regulator: with the propeller running, regulator output
        #will be saturated at -1 due to proportional error, so the integrator will not
        #accumulate and the state will not change
        x_pwp_friction = 0.0,
    )

    #TrimParameters
    trim_parameters = (
        TAS = 40,
        β = 0.0,
        Ob = GeographicLocation(l2d = LatLon(), alt = AltO(500)),
        ψ_lb = 0.0, #LTF-relative heading (LTF initialized to NED)
        brake_left = 0., #left mlg brake, 0 to 1, irrelevant in flight
        brake_right = 0., #right mlg brake, 0 to 1, irrelevant in flight
    )

    #the derivatives of x_idle and x_thr_friction could be included in the
    #optimization cost to ensure they will remain zero

    f = get_target_function(aircraft, atmosphere, terrain, trim_parameters, trim_configuration)

    x_trim_0 = ComponentVector(
    ω_lb_b = [0, 0, 0],
    α = 0.08,
    θ_lb = 0.08,
    φ_lb = 0.0,
    ω_eng = 215.0,
    throttle = 0.61,
    yoke_Δx = 0.01, #relative to force-free position
    yoke_Δy = -0.025, #relative to force-free position
    pedals = 0.0,
    )

    f(x_trim_0)
    # println(@ballocated($f($x_trim_0))) #no allocations! nice!

end

function get_target_function(aircraft, atmosphere, terrain, trim_parameters, trim_configuration)

    #is this let really needed? or being inside a function already freezes values
    func = let  x = aircraft.x, d = aircraft.d, u = aircraft.u, atmosphere = atmosphere, terrain = terrain,
                trim_parameters = trim_parameters, trim_configuration = trim_configuration

        #set the engine running, no way to trim otherwise
        d.vehicle.pwp.engine.state = Piston.eng_running

        u.avionics.flaps = trim_configuration.flaps
        u.avionics.yoke_x0 = trim_configuration.yoke_x0
        u.avionics.yoke_y0 = trim_configuration.yoke_y0
        u.vehicle.pwp.engine.mix = trim_configuration.mixture
        x.vehicle.fuel .= trim_configuration.fuel

        x.vehicle.pwp.engine.idle .= trim_configuration.x_pwp_idle
        x.vehicle.pwp.friction .= trim_configuration.x_pwp_friction

        function compute_derivative(x_trim)

            error("Convertir primero a KinInit para independizar de la cinematica elegida. Alternativamente, imponer que la cinematica sea LTF al crear aircraft. Mejor lo primero, porque KinInit allocates, habria que hacer una version in place")
            ψ_nl = 0 #align NED to LTF, so q_lb = q_nb
            q_nb = q_lb = REuler(trim_parameters.ψ_lb, x_trim.θ_lb, x_trim.φ_lb) |> RQuat
            h_e = AltE(trim_parameters.Ob)
            x.kinematics.pos.q_lb .= q_lb[:]
            x.kinematics.pos.q_el .= ltf(trim_parameters.Ob, ψ_nl)[:]
            x.kinematics.pos.h_e = Float64(h_e)

            v_ew_n = AtmosphericData(atmosphere, trim_parameters.Ob).wind.v_ew_n
            v_ew_b = q_nb'(v_ew_n) #wind velocity, body axes
            v_wOb_a = Air.get_velocity_vector(trim_parameters.TAS, x_trim.α, trim_parameters.β)
            v_wOb_b = v_wOb_a #wind-relative aircraft velocity, body frame (aerodynamic frame = body frame)
            v_eOb_b = v_ew_b + v_wOb_b
            x.kinematics.vel.v_eOb_b = v_eOb_b

            v_eOb_n = q_nb(v_eOb_b)
            (R_N, R_E) = radii(trim_parameters.Ob)
            ω_el_n = SVector{3}(
                v_eOb_n[2] / (R_E + Float64(h_e)),
                -v_eOb_n[1] / (R_N + Float64(h_e)),
                0.0)
            ω_el_b = q_nb'(ω_el_n)
            x.kinematics.vel.ω_eb_b = ω_el_b + x_trim.ω_lb_b

            x.vehicle.aero.α_filt = x_trim.α #this ensures zero state derivative
            x.vehicle.aero.β_filt = trim_parameters.β #this ensures zero state derivative
            x.vehicle.pwp.ω = x_trim.ω_eng

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