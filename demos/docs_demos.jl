using Flight
using ControlSystems, RobustAndOptimalControl, Plots, LaTeXStrings

function nlsim_q()

    #set up nonlinear simulation
    world = Model(SimpleWorld(aircraft = Cessna172Sv0(NED())))
    sim = Simulation(world; t_end = 10)

    #initialize at default trim condition
    trim_params = C172.TrimParameters()
    init!(sim, trim_params)

    #advance simulation for one second from trim condition
    step!(sim, 1, true)
    #apply 10% elevator increment
    world.aircraft.vehicle.systems.act.u.elevator += 0.1
    #run to completion
    run!(sim)
    #extract pitch rate from simulation results
    ts = TimeSeries(sim);
    _, q_nonlinear, _ = get_components(ts.aircraft.vehicle.kinematics.ω_wb_b)

    #extract aircraft submodel and linearize it around the trim condition
    lss = LinearizedSS(world.aircraft, trim_params)
    #convert to NamedStateSpace
    nss = named_ss(lss)
    #extract elevator to pitch rate SISO system
    e2q = nss[:q, :elevator]
    #simulate a 0.1 step input at t=1
    y, t, _, _ = lsim(e2q, (x, t)->[0.1]*(t>=1), 0:0.01:10)
    q_linear = vec(y)

    plot(q_nonlinear; plot_title = "Pitch Rate", label = "Nonlinear", ylabel=L"$q \ (rad/s)$") |> display
    plot!(t, q_linear, label = "Linear")

end

# Elevator step response in nonlinear and linearized Cessna172S models
function nlsim_θ()

#1. Set up and run nonlinear simulation

        #instantiate the aircraft with NED kinematics, required for linearization
        world = SimpleWorld(aircraft = Cessna172Sv0(NED())) |> Model
        sim = Simulation(world; t_end = 10)

        #define trim conditions and initialize Simulation
        trim_params = C172.TrimParameters()
        init!(sim, trim_params)

        #advance 1 second from trim condition
        step!(sim, 1, true)

        #apply 10% elevator increment
        world.aircraft.vehicle.systems.act.u.elevator += 0.1

        #run to completion
        run!(sim)

        #extract pitch angle TimeSeries from simulation results
        ts = TimeSeries(sim)
        θ_nonlinear = ts.aircraft.vehicle.kinematics.e_nb.θ

    #2. Obtain linear SISO system

        #extract aircraft submodel and linearize it around the trim condition
        lss = LinearizedSS(world.aircraft, trim_params)

        #convert to NamedStateSpace
        nss = named_ss(lss)

        #extract elevator-to-pitch angle SISO system
        e2θ = nss[:θ, :elevator]

    #3. Compute linear response to elevator step input

        #simulate a 0.1 step input applied at t=1
        y, t, _, _ = lsim(e2θ, (x, t)->[0.1]*(t>=1), 0:0.01:10)

        #get perturbation Δθ around trim condition
        Δθ_linear = vec(y)

        #retrieve trim θ value from linearized aircraft model
        θ_trim = lss.y0[:θ]

        #compute total linear θ response
        θ_linear = θ_trim .+ Δθ_linear

    #4. Compare responses

        plot(θ_nonlinear; plot_title = "Pitch Angle", label = "Nonlinear", ylabel=L"$\theta \ (rad)$", size = (900, 600)) |> display
        plot!(t, θ_linear; label = "Linear")

    # savefig("step_response.png")

end

#Simulate an automated turning climb under constant wind conditions on the
#fly-by-wire Cessna172X model

function turning_climb()

    #1. Set up simulation

        #custom fly-by-wire Cessna172 variant with default environment
        world = SimpleWorld(; aircraft = Cessna172Xv2()) |> Model
        sim = Simulation(world; t_end = 600)

        #initialize using default trim conditions
        init!(sim, C172.TrimParameters())

    #3. Set wind conditions

        world.atmosphere.wind.u.N = 1.0 #1 m/s North
        world.atmosphere.wind.u.E = 0.5 #0.5 m/s East

    #3. Configure autopilot for turning climb

        #extract control laws submodel
        ctl = world.aircraft.avionics.ctl

        #set longitudinal control laws to track airspeed and climb rate
        ctl.lon.u.mode_req = C172XControl.ModeControlLon.EAS_clm

        #set lateral control laws to track bank angle and sideslip angle
        ctl.lat.u.mode_req = C172XControl.ModeControlLat.φ_β

        #set climb rate reference to 2.0 m/s, keep trim airspeed
        ctl.lon.u.clm_ref = 2.0

        #set bank angle reference to 30 degrees, keep trim sideslip angle
        ctl.lat.u.φ_ref = deg2rad(30)

    #4. Run Simulation and extract results

        run!(sim)
        ts = TimeSeries(sim)

    #5. Plot 3D trajectory

        kin_plots = make_plots(ts.aircraft.vehicle.kinematics; size = (900, 600))
        display(kin_plots[:Ob_t3d])

        # savefig(kin_plots[:Ob_t3d],"turning_climb_3d.png")

end


function elevator_doublet()

    world = SimpleWorld(aircraft = Cessna172Xv1()) |> Model
    sim = Simulation(world; dt = 0.02, t_end = 60)
    init!(sim, C172.TrimParameters())

    step!(sim, 5)
    world.aircraft.avionics.ctl.u.elevator_offset = 0.1
    step!(sim, 2)
    world.aircraft.avionics.ctl.u.elevator_offset = -0.1
    step!(sim, 2)
    world.aircraft.avionics.ctl.u.elevator_offset = 0
    Sim.run!(sim)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/elevator_doublet/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/elevator_doublet/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/elevator_doublet/dyn"); Plotting.defaults...)

    return sim

end

function elevator_doublet_callback()

    world = SimpleWorld(aircraft = Cessna172Xv1()) |> Model

    user_callback! = function(world::Model)
        t = world.t[]
        u_ctl = world.aircraft.avionics.lon.u
        if 5 <= t < 7
            u_ctl.elevator_offset = 0.1
        elseif 7 <= t < 9
            u_ctl.elevator_offset = -0.1
        else
            u_ctl.elevator_offset = 0
        end
    end

    sim = Simulation(world; dt = 0.02, t_end = 60, user_callback!)
    init!(sim, C172.TrimParameters())
    Sim.run!(sim)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/elevator_doublet/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/elevator_doublet/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/elevator_doublet/dyn"); Plotting.defaults...)

    return sim

end
