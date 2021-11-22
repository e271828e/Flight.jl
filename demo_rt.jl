using Flight
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using BenchmarkTools
using Base.Iterators


plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12)

################################## AC2 ######################################

trn = DummyTerrain()
atm_sys = System(AtmosphereDescriptor());

ac = AircraftBase(
    kin = KinLTF(),
    mass = TunableMass(),
    aero = SimpleDrag(),
    pwp = AirframeGroup((
        left = EThruster(motor = ElectricMotor(α = CW)),
        right = EThruster(motor = ElectricMotor(α = CCW)))),
);

kin_init = KinInit(v_eOb_b = [50, 0, 0],
                    ω_lb_b = [0, 0, 0.2],
                    q_nb = REuler(ψ = 0, θ = 0, φ = 0),
                    Ob = Geographic(LatLon(ϕ = deg2rad(40.531818), λ = deg2rad(-3.574862)),
                                    AltOrth(811)));
atm_sys.u.wind.v_ew_n[1] = -5
ac_sys = System(ac);
init!(ac_sys.x.kin, kin_init)
ac_sys.subsystems.pwp.u.left.throttle = .1
ac_sys.subsystems.pwp.u.right.throttle = .1

ac_mdl = Model(ac_sys, (trn, atm_sys); dt = 0.01, adaptive = false, method = Heun(), y_saveat = 0.1);

xp = XPInterface()
disable_physics(xp)
set_position(xp, ac_mdl.sys.y.kin.pos)
output_div = 4

t_wall = time()
t_wall_0 = t_wall

# for i in take(ac_mdl.integrator, 5)
for i in ac_mdl.integrator

    #when using the integrator as an iterator, we don't need to step, it is done
    #automatically at the beginning of each iteration. therefore, the dt we need
    #for pacing the simulation is not the proposed dt for the next step, but the
    #one in the step the integrator has already taken.

    #retrieve the dt just taken by the integrator
    dt = ac_mdl.dt

    #compute the wall time corresponding to the newly updated simulation
    t_wall_next = t_wall + dt

    #busy wait until the wall time catches up
    while (time() < t_wall_next) end
    # println(time()-t_wall_next)

    t_wall = t_wall_next
    # println(t_wall - t_wall_0)

    if ac_mdl.success_iter % output_div == 0
        set_position(xp, ac_mdl.sys.y.kin.pos)
    end

    #when GLFW is implemented, build the model with t_end = Inf and use this to break
    # if !GLFW.WindowShouldClose(window)
    #     break
    # end

end



# plots(ac_mdl; save_path = joinpath("tmp", "plots_ac_rt"), plot_settings...)