module C172Demos

using Sockets
using UnPack, JSON3

using Flight
using Flight.FlightAircraft.C172Y.C172YControl: ModeControlLon, ModeControlLat, is_on_gnd
using Flight.FlightAircraft.C172Y.C172YGuidance: ModeGuidance, Segment, SegmentGuidanceData

#position and geographic heading for Salzburg airport (LOWS) runway 15
const loc_LOWS15 = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997))
const h_LOWS15 = HOrth(427.2)
const ψ_LOWS15 = deg2rad(157)


############################# Tutorial 01 ######################################

function tutorial01(; aircraft::Cessna172 = Cessna172Xv1(),
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

    Sim.init!(sim, initializer)

    xp = XPlane12Control(address = xp12_address, port = xp12_port)
    Sim.attach!(sim, xp)
    for joystick in update_connected_joysticks()
        # isa(joystick, Joysticks.T16000M) && Sim.attach!(sim, joystick)
        Sim.attach!(sim, joystick)
    end

    Sim.run!(sim; gui = true)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/tutorial01/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/tutorial01/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/tutorial01/dyn"); Plotting.defaults...)
    # save_plots(TimeSeries(sim).aircraft.vehicle.dynamics; Plotting.defaults...)

    return sim

end


############################# Tutorial 02 ######################################

function tutorial02a()

    mdl = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain()) |> Model
    sim = Simulation(mdl; dt = 0.02, t_end = 60)
    Sim.init!(sim, C172.TrimParameters())

    Sim.step!(sim, 5)
    mdl.aircraft.avionics.ctl.u.elevator_offset = 0.1
    Sim.step!(sim, 2)
    mdl.aircraft.avionics.ctl.u.elevator_offset = -0.1
    Sim.step!(sim, 2)
    mdl.aircraft.avionics.ctl.u.elevator_offset = 0
    Sim.run!(sim)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/tutorial02/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/tutorial02/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/tutorial02/dyn"); Plotting.defaults...)

    return sim

end

function tutorial02b()

    mdl = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain()) |> Model

    user_callback! = function(mdl::Model)
        t = mdl.t[]
        u_ctl = mdl.aircraft.avionics.ctl.u
        if 5 <= t < 7
            u_ctl.elevator_offset = 0.1
        elseif 7 <= t < 9
            u_ctl.elevator_offset = -0.1
        else
            u_ctl.elevator_offset = 0
        end
    end

    sim = Simulation(mdl; dt = 0.02, t_end = 60, user_callback!)
    Sim.init!(sim, C172.TrimParameters())
    Sim.run!(sim)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/tutorial02/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/tutorial02/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/tutorial02/dyn"); Plotting.defaults...)

    return sim

end

#############################

function crosswind_landing(; gui::Bool = false)

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

    user_callback! = let phase = Ref(:final)

        function(mdl::Model)

            @unpack aircraft = mdl
            @unpack vehicle, avionics = aircraft
            @unpack act, pwp = vehicle.systems

            if phase[] === :final
                #preselect modes and command references
                avionics.gdc.u.mode_req = ModeGuidance.segment
                avionics.gdc.seg.u.target = final_leg
                avionics.gdc.seg.u.horizontal_req = true
                avionics.gdc.seg.u.vertical_req = true
                avionics.ctl.lon.u.EAS_ref = 30
                vehicle.systems.act.flaps.u[] = 1.0
                if vehicle.y.kinematics.h_e - final_leg.p2.h < 3
                    println("Entering flare")
                    phase[] = :flare
                end
            elseif phase[] === :flare
                avionics.gdc.u.mode_req = ModeGuidance.direct
                avionics.ctl.lon.u.mode_req = ModeControlLon.thr_EAS
                avionics.ctl.lon.u.throttle_axis = 0
                avionics.ctl.lon.u.EAS_ref = 25
                avionics.ctl.lat.u.mode_req = ModeControlLat.φ_β
                avionics.ctl.lat.u.φ_ref = 0
                avionics.ctl.lat.u.β_ref = 0
            end

        end #function

    end

    mdl = SimpleWorld(Cessna172Yv2(), SimpleAtmosphere(), HorizontalTerrain(h_LOWS15)) |> Model

    sim = Simulation(mdl; dt = 0.02, t_end = 1000, user_callback!)
    Sim.init!(sim, initializer)
    Sim.run!(sim; gui)

    save_plots(TimeSeries(sim).aircraft.vehicle.kinematics, normpath("tmp/plots/traffic_pattern/kin"); Plotting.defaults..., linewidth = 2,)
    save_plots(TimeSeries(sim).aircraft.vehicle.airflow, normpath("tmp/plots/traffic_pattern/air"); Plotting.defaults...)
    save_plots(TimeSeries(sim).aircraft.vehicle.dynamics, normpath("tmp/plots/traffic_pattern/dyn"); Plotting.defaults...)

    return sim

end




#############################

function traffic_pattern(; gui::Bool = false)

    #build a standard traffic pattern around runway 15
    p_LOWS15 = Geographic(loc_LOWS15, h_LOWS15)
    final_leg = -Segment(p_LOWS15, χ = ψ_LOWS15 + π, s = 4e3, γ = deg2rad(3))
    base_leg = -Segment(final_leg.p1, χ = ψ_LOWS15 - π/2, s = 2e3, γ = 0)
    downwind_leg = -Segment(base_leg.p1, χ = ψ_LOWS15, s = 8e3, γ = 0)
    crosswind_leg = -Segment(downwind_leg.p1, χ = ψ_LOWS15 + π/2, s = 2e3, γ = 0)
    departure_leg = Segment(p_LOWS15, crosswind_leg.p1)
    # @show SegmentGuidanceData(departure_leg, departure_leg.p1).γ_12 |> rad2deg

    capture_threshold = -300

    user_callback! = let phase = Ref(:standby)
        function(mdl::Model)
            @unpack aircraft = mdl
            @unpack vehicle, avionics = aircraft
            @unpack act, pwp = vehicle.systems
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
                avionics.gdc.u.mode_req = ModeGuidance.segment
                avionics.gdc.seg.u.target = departure_leg
                avionics.gdc.seg.u.horizontal_req = true
                avionics.gdc.seg.u.vertical_req = true
                avionics.ctl.lon.u.EAS_ref = 35
                avionics.ctl.lon.u.throttle_axis = 1
                if !is_on_gnd(vehicle)
                    println("Lift-off")
                    phase[] = :departure
                end
            elseif phase[] === :departure
                if avionics.y.gdc.seg.data.s_2b > capture_threshold
                    avionics.gdc.seg.u.target = crosswind_leg
                    println("Entering crosswind")
                    phase[] = :crosswind
                end
            elseif phase[] === :crosswind
                avionics.ctl.lon.u.EAS_ref = 50
                if avionics.y.gdc.seg.data.s_2b > capture_threshold
                    avionics.gdc.seg.u.target = downwind_leg
                    println("Entering downwind")
                    phase[] = :downwind
                end
            elseif phase[] === :downwind
                if avionics.y.gdc.seg.data.s_2b > capture_threshold
                    avionics.gdc.seg.u.target = base_leg
                    println("Entering base")
                    phase[] = :base
                end
            elseif phase[] === :base
                avionics.ctl.lon.u.EAS_ref = 30
                vehicle.systems.act.flaps.u[] = 1.0
                if avionics.y.gdc.seg.data.s_2b > capture_threshold
                    avionics.gdc.seg.u.target = final_leg
                    println("Entering final")
                    phase[] = :final
                end
            elseif phase[] === :final
                if vehicle.y.kinematics.h_e - final_leg.p2.h < 10
                    println("Entering flare")
                    phase[] = :flare
                end
            elseif phase[] === :flare
                avionics.gdc.u.mode_req = ModeGuidance.direct
                avionics.ctl.lon.u.mode_req = ModeControlLon.thr_EAS
                avionics.ctl.lon.u.throttle_axis = 0
                avionics.ctl.lon.u.EAS_ref = 25
                avionics.ctl.lat.u.mode_req = ModeControlLat.φ_β
                avionics.ctl.lat.u.φ_ref = 0
                avionics.ctl.lat.u.β_ref = 0
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

    mdl = SimpleWorld(Cessna172Yv2(), SimpleAtmosphere(), HorizontalTerrain(h_LOWS15)) |> Model

    sim = Simulation(mdl; dt = 0.02, t_end = 1000, user_callback!)
    Sim.init!(sim, initializer)
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
        #defined in C172Y.C172YControl
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
    world = SimpleWorld(Cessna172Yv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

    sim = Simulation(world; t_end = 30)

    #on air, automatically trimmed by reinit!
    initializer = C172.TrimParameters(
        Ob = Geographic(LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)), HEllip(650)))

    #initialize simulated system
    Sim.init!(sim, initializer)

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
                        normpath("tmp/plots/test_c172y1/test_json_loopback/kin");
                        Plotting.defaults...)

    return nothing

end



end #module