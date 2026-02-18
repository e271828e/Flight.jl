module C172Demos

using ControlSystems, RobustAndOptimalControl, Plots, LaTeXStrings, JSON3, Sockets

using FlightCore
using FlightLib
using FlightApps
using FlightApps.C172: is_on_gnd
using FlightApps.C172X.C172XControl: ModeControlLon, ModeControlLat
using FlightApps.C172X.C172XGuidance: ModeGuidance, Segment, SegmentGuidanceData

export interactive_simulation, crosswind_landing, traffic_pattern, json_loopback

#position and geographic heading for Salzburg airport (LOWS) runway 15
const loc_LOWS15 = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997))
const h_LOWS15 = HOrth(427.2)
const ψ_LOWS15 = deg2rad(157)


######################### Documentation Demos ##################################

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
    lss = linearize(world.aircraft, trim_params)
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
        lss = linearize(world.aircraft, trim_params)

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


############################## Misc Demos ######################################

function interactive_simulation(; aircraft::Cessna172 = Cessna172Xv1(),
                situation::Symbol = :ground,
                xp12_address = IPv4("127.0.0.1"),
                xp12_port = 49000,
                )

    #default atmospheric model
    atmosphere = SimpleAtmosphere()

    #horizontal terrain model with elevation matching LOWS runway 15
    terrain = HorizontalTerrain(h_LOWS15)

    #define world and build Model for simulation
    world = SimpleWorld(aircraft, atmosphere, terrain) |> Model

    if situation === :ground
        #initial condition specified through aircraft frame kinematics
        initializer = KinInit(;
            location = loc_LOWS15, #2D location
            h = h_LOWS15 + C172.Δh_to_gnd, #altitude
            q_nb = REuler(ψ_LOWS15, 0, 0), #attitude with respect to NED frame
            ω_wb_b = zeros(3), #angular velocity
            v_eb_n = zeros(3), #velocity
            ) |> C172.Init

    elseif situation === :air
        #initial condition specified by trim parameters
        # initializer = C172.TrimParameters()
        initializer = C172.TrimParameters(;
            Ob = Geographic(loc_LOWS15, h_LOWS15 + 500), #500 m above LOWS runway 15
            EAS = 50.0, #equivalent airspeed
            ψ_nb = ψ_LOWS15, #geographic heading
            γ_wb_n = 0.0, #wind-relative flight path angle
            ψ_wb_dot = 0.0, #turn rate
            flaps = 0.0, #flap setting
            fuel_load = 0.5, #available fuel fraction
        )
    else
        error("Unknown situation: $situation")
    end

    #create a Simulation with specified integration step size and stop time
    sim = Simulation(world; dt = 0.02, t_end = 1000)

    init!(sim, initializer)

    xp = XPlane12Control(address = xp12_address, port = xp12_port)
    Sim.attach!(sim, xp)
    for joystick in update_connected_joysticks()
        # isa(joystick, Joysticks.T16000M) && Sim.attach!(sim, joystick)
        Sim.attach!(sim, joystick)
    end

    Sim.run!(sim; gui = true)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/interactive_sim/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/interactive_sim/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/interactive_sim/dyn"); Plotting.defaults...)
    # save_plots(TimeSeries(sim).aircraft.vehicle.dynamics; Plotting.defaults...)

    return sim

end


################################################################################

function crosswind_landing(; gui::Bool = false,
                            xp12_address = IPv4("127.0.0.1"),
                            xp12_port = 49000,)

    final_leg = -Segment(Geographic(loc_LOWS15, h_LOWS15), χ = ψ_LOWS15 + π, s = 3e3, γ = deg2rad(3))

    initializer = C172.TrimParameters(;
        Ob = final_leg.p1,
        EAS = 30.0, #equivalent airspeed
        ψ_nb = ψ_LOWS15, #geographic heading
        γ_wb_n = -deg2rad(3), #wind-relative flight path angle
        ψ_wb_dot = 0.0, #turn rate
        flaps = 1.0, #flap setting
        fuel_load = 0.5, #available fuel fraction
    )

    user_callback! = let phase = Ref(:init)

        function(mdl::Model)

            (; aircraft, atmosphere) = mdl
            (; vehicle, avionics) = aircraft
            (; gdc, ctl) = avionics
            (; seg) = gdc
            (; act) = vehicle.systems

            atmosphere.wind.u.E = 6

            if phase[] === :init

                gdc.u.mode_req = ModeGuidance.segment
                seg.u.target = final_leg
                seg.u.hor_gdc_req = true
                seg.u.vrt_gdc_req = true
                ctl.lon.u.EAS_ref = 30
                vehicle.systems.act.flaps.u[] = 1.0

                # println("Entering final")
                phase[] = :final

            elseif phase[] === :final

                if vehicle.y.kinematics.h_e - final_leg.p2.h < 6
                    seg.u.vrt_gdc_req = false
                    ctl.lon.u.mode_req = ModeControlLon.EAS_clm
                    ctl.lon.u.clm_ref = -0.3

                    ctl.lat.u.mode_req = ModeControlLat.φ_β
                    ψ_current = vehicle.y.kinematics.e_nb.ψ
                    ψ_seg = seg.y.data.χ_12
                    ctl.lat.u.β_ref = Attitude.wrap_to_π(ψ_current - ψ_seg)
                    ctl.lat.u.φ_ref = 0

                    # println("Entering flare")
                    phase[] = :flare
                end

            elseif phase[] === :flare

                if is_on_gnd(vehicle)
                    ctl.lon.u.throttle_axis = 0
                    ctl.lat.u.rudder_axis = -0.04
                    vehicle.systems.act.flaps.u[] = 0.0
                    # println("Touchdown ")
                    phase[] = :ground
                end

            elseif phase[] === :ground

                ctl.lon.u.throttle_axis = 0
                act.brake_left.u[] = 1
                act.brake_right.u[] = 1

            end

        end #function

    end

    mdl = SimpleWorld(Cessna172Xv2(), SimpleAtmosphere(), HorizontalTerrain(h_LOWS15)) |> Model

    sim = Simulation(mdl; dt = 0.02, t_end = 1000, user_callback!)

    xp = XPlane12Control(address = xp12_address, port = xp12_port)
    Sim.attach!(sim, xp)

    init!(sim, initializer)
    Sim.run!(sim; gui)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/crosswind_landing/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/crosswind_landing/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/crosswind_landing/dyn"); Plotting.defaults...)

    return sim

end


################################################################################

function traffic_pattern(; gui::Bool = false,
                            xp12_address = IPv4("127.0.0.1"),
                            xp12_port = 49000,)

    #build a standard traffic pattern around runway 15
    p_LOWS15 = Geographic(loc_LOWS15, h_LOWS15)
    final_leg = -Segment(p_LOWS15, χ = ψ_LOWS15 + π, s = 3e3, γ = deg2rad(3))
    base_leg = -Segment(final_leg.p1, χ = ψ_LOWS15 - π/2, s = 1e3, γ = 0)
    downwind_leg = -Segment(base_leg.p1, χ = ψ_LOWS15, s = 6e3, γ = 0)
    crosswind_leg = -Segment(downwind_leg.p1, χ = ψ_LOWS15 + π/2, s = 1e3, γ = 0)
    departure_leg = Segment(p_LOWS15, crosswind_leg.p1)
    # @show SegmentGuidanceData(departure_leg, departure_leg.p1).γ_12 |> rad2deg

    capture_threshold = -200

    user_callback! = let phase = Ref(:standby)
        function(mdl::Model)

            (; aircraft) = mdl
            (; vehicle, avionics) = aircraft
            (; gdc, ctl) = avionics
            (; seg) = gdc
            (; act, pwp) = vehicle.systems

            t = mdl.t[]

            if phase[] === :standby

                t >= 5 && (phase[] = :startup)

            elseif phase[] === :startup

                pwp.engine.u.start = true
                if pwp.engine.y.state === Piston.EngineState.running
                    pwp.engine.u.start = false
                    phase[] = :takeoff
                end

            elseif phase[] === :takeoff

                #preselect modes and command references
                gdc.u.mode_req = ModeGuidance.segment
                seg.u.target = departure_leg
                seg.u.hor_gdc_req = true
                seg.u.vrt_gdc_req = true
                ctl.lon.u.EAS_ref = 35
                ctl.lon.u.throttle_axis = 1
                if !is_on_gnd(vehicle)
                    # println("Lift-off")
                    phase[] = :departure
                end

            elseif phase[] === :departure

                if avionics.y.gdc.seg.data.s_2b > capture_threshold
                    seg.u.target = crosswind_leg
                    # println("Entering crosswind")
                    phase[] = :crosswind
                end

            elseif phase[] === :crosswind

                if avionics.y.gdc.seg.data.s_2b > capture_threshold
                    seg.u.target = downwind_leg
                    # println("Entering downwind")
                    phase[] = :downwind
                end

            elseif phase[] === :downwind

                avionics.ctl.lon.u.EAS_ref = 50
                if avionics.y.gdc.seg.data.s_2b > capture_threshold
                    seg.u.target = base_leg
                    # println("Entering base")
                    phase[] = :base
                end

            elseif phase[] === :base

                avionics.ctl.lon.u.EAS_ref = 30
                vehicle.systems.act.flaps.u[] = 1.0
                if avionics.y.gdc.seg.data.s_2b > capture_threshold
                    seg.u.target = final_leg
                    # println("Entering final")
                    phase[] = :final
                end

            elseif phase[] === :final

                if vehicle.y.kinematics.h_e - final_leg.p2.h < 6
                    seg.u.vrt_gdc_req = false
                    ctl.lon.u.mode_req = ModeControlLon.EAS_clm
                    ctl.lon.u.clm_ref = -0.3

                    ctl.lat.u.mode_req = ModeControlLat.φ_β
                    ψ_current = vehicle.y.kinematics.e_nb.ψ
                    ψ_seg = seg.y.data.χ_12
                    ctl.lat.u.β_ref = Attitude.wrap_to_π(ψ_current - ψ_seg)
                    ctl.lat.u.φ_ref = 0

                    # println("Entering flare")
                    phase[] = :flare
                end

            elseif phase[] === :flare

                if is_on_gnd(vehicle)
                    ctl.lon.u.throttle_axis = 0
                    ctl.lat.u.rudder_axis = -0.04
                    vehicle.systems.act.flaps.u[] = 0.0
                    # println("Touchdown")
                    phase[] = :ground
                end

            elseif phase[] === :ground

                ctl.lon.u.throttle_axis = 0
                act.brake_left.u[] = 1
                act.brake_right.u[] = 1

            end
        end #function
    end

    initializer = KinInit(;
        location = loc_LOWS15, #2D location
        h = h_LOWS15 + C172.Δh_to_gnd, #altitude
        q_nb = REuler(ψ_LOWS15, 0, 0), #attitude with respect to NED frame
        ω_wb_b = zeros(3), #angular velocity
        v_eb_n = zeros(3), #velocity
        ) |> C172.Init

    mdl = SimpleWorld(Cessna172Xv2(), SimpleAtmosphere(), HorizontalTerrain(h_LOWS15)) |> Model

    sim = Simulation(mdl; dt = 0.02, t_end = 1000, user_callback!)

    xp = XPlane12Control(address = xp12_address, port = xp12_port)
    Sim.attach!(sim, xp)

    init!(sim, initializer)
    Sim.run!(sim; gui)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/traffic_pattern/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/traffic_pattern/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/traffic_pattern/dyn"); Plotting.defaults...)

    return sim

end


############################# UDP/JSON Loopback Test ###########################

struct JSONTestMapping <: IOMapping end

function IODevices.extract_output(mdl::Model{<:SimpleWorld}, ::JSONTestMapping)
    freq = 0.1
    φ_ref_max = π/6
    φ_ref = φ_ref_max * sin(2π*freq*mdl.t[])
    clm_ref = 0.0

    #these are all valid empty JSON entities. when passed to JSON3.write, they
    #yield respectively "\"\"", "[]" and "{}", all of length 2
    cmd = ""
    # cmd = []
    # cmd = Dict()

    if mdl.t[] > 5
        #enums will be automatically cast to Ints per the StructTypes methods
        #defined in C172X.C172XControl
        cmd = (
            lon = (
                mode_req = ModeControlLon.EAS_clm, #7 would also work
                clm_ref = clm_ref,
            ),
            lat = (
                mode_req = ModeControlLat.φ_β,
                φ_ref = φ_ref,
            )
        )
    end

    return JSON3.write(cmd)
end

function IODevices.assign_input!(world::Model{<:SimpleWorld}, ::JSONTestMapping, data::String)
    #caution: String(data) empties the original data::Vector{UInt8}, so
    #additional calls would return an empty string
    str = String(data)
    u = JSON3.read(str)

    if !isempty(u)
        JSON3.read!(JSON3.write(u.lon), world.aircraft.avionics.u.lon)
        JSON3.read!(JSON3.write(u.lat), world.aircraft.avionics.u.lat)
    end
end

function json_loopback(; gui::Bool = true, xp12 = false, save::Bool = true)


    h_trn = HOrth(427.2);
    world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

    sim = Simulation(world; t_end = 30)

    #on air, automatically trimmed by init!
    initializer = C172.TrimParameters(
        Ob = Geographic(LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)), HEllip(650)))

    #initialize simulated system
    init!(sim, initializer)

    xp12 && Sim.attach!(sim, XPlane12Control(; port = 49000))

    #the loopback interface must not share its port with XPlane12Control!
    Sim.attach!(sim, UDPInput(; port = 49017), JSONTestMapping())
    Sim.attach!(sim, UDPOutput(; port = 49017), JSONTestMapping())

    #trigger compilation of parsing methods before launching the simulation
    JSON3.read!(JSON3.write(world.aircraft.avionics.u.lon, allow_inf=true),
                            world.aircraft.avionics.u.lon; allow_inf=true)
    JSON3.read!(JSON3.write(world.aircraft.avionics.u.lat, allow_inf=true),
                            world.aircraft.avionics.u.lat; allow_inf=true)

    #set non-Inf pace for headless runs to allow the UDP interface to keep up
    pace = (gui ? 1.0 : 30)
    Sim.run!(sim; gui, pace)

    save && save_plots(TimeSeries(sim).aircraft.vehicle.kinematics,
                        normpath("tmp/plots/udp_loopback/kin");
                        Plotting.defaults...)

    return nothing

end

end #module