module TestC172

using Test
using UnPack
using BenchmarkTools
using Sockets
using OrdinaryDiffEq: Heun, RK4, Tsit5

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

export test_c172


function test_c172()
    @testset verbose = true "Cessna 172" begin
    end
end

function test_sim(; ac::Cessna172 = Cessna172Sv0(),
                loc::Abstract2DLocation = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)),
                h_trn::HOrth = HOrth(427.2),
                ψ::Real = deg2rad(157),
                situation::Symbol = :ground,
                interactive::Bool = false,
                t_end::Real = 1000,
                save::Bool = true)

    world = SimpleWorld(ac, SimpleAtmosphere(), HorizontalTerrain(altitude = h_trn)) |> System

    if situation === :ground
        initializer = KinInit(; loc, q_nb = REuler(ψ, 0, 0), h = h_trn + C172.Δh_to_gnd,
            ω_wb_b = zeros(3), v_eb_n = zeros(3));
    elseif situation === :air
        EAS = 50.0
        flaps = 0.0
        γ_wb_n = 0.0
        x_fuel = 0.5
        payload = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)

        initializer = C172.TrimParameters(; Ob = Geographic(loc, HEllip(650)), EAS, γ_wb_n, x_fuel, flaps, payload)
    else
        error("Unknown situation: $situation")
    end

    user_callback! = let
        function (world)
            t = world.t[]
        end
    end

    sim = Simulation(world; algorithm = RK4(), dt = 1/60, t_end, user_callback!)
    # sim = Simulation(world; algorithm = Tsit5(), dt = 1/60, t_end, user_callback!)
    # sim = Simulation(world; algorithm = Heun(), dt = 0.01, t_end, user_callback!)
    Sim.init!(sim, initializer)

    if interactive
        xp = XPlane12Output(address = IPv4("127.0.0.1"), port = 49000)
        for joystick in update_connected_joysticks()
            Sim.attach!(sim, joystick)
        end
        Sim.attach!(sim, xp)
        Sim.run_interactive!(sim)
    else
        Sim.run!(sim)
    end

    kin_plots = make_plots(TimeSeries(sim).ac.vehicle.kinematics; Plotting.defaults...)
    air_plots = make_plots(TimeSeries(sim).ac.vehicle.airflow; Plotting.defaults...)
    dyn_plots = make_plots(TimeSeries(sim).ac.vehicle.dynamics; Plotting.defaults...)
    save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172", "sim", "kin"))
    save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172", "sim", "air"))
    save && save_plots(dyn_plots, save_folder = joinpath("tmp", "test_c172", "sim", "dyn"))

end

function test_xp12()

    # https://developer.x-plane.com/article/using-the-correct-wing-datarefs/
    override_pose = "sim/operation/override/override_planepath[0]"
    override_surf = "sim/operation/override/override_control_surfaces[0]"
    override_prop = "sim/flightmodel2/engines/prop_disc/override[0]"
    override_nws = "sim/operation/override/override_wheel_steer[0]"
    rele_pos = "sim/flightmodel2/wing/elevator1_deg[8]"
    lele_pos = "sim/flightmodel2/wing/elevator1_deg[9]"
    lflap_pos = "sim/flightmodel2/wing/flap1_deg[0]"
    rflap_pos = "sim/flightmodel2/wing/flap1_deg[1]"
    rud_pos = "sim/flightmodel2/wing/rudder1_deg[10]"
    lail_pos = "sim/flightmodel2/wing/aileron1_deg[2]"
    rail_pos = "sim/flightmodel2/wing/aileron1_deg[3]"
    prop_is_disc = "sim/flightmodel2/engines/prop_is_disc[0]"
    prop_angle = "sim/flightmodel2/engines/prop_rotation_angle_deg[0]"
    nws_angle = "sim/flightmodel2/gear/tire_steer_actual_deg[0]"

    δe = 0.5
    δa = 0.5
    δr = 0.5
    δf = 0.5

    t = 0.02
    ω_prop = 100 #rad/s
    ϕ_prop = ω_prop * t

    ψ_sw = 0.5 #rad

    address = IPv4("127.0.0.1")
    # address = IPv4("192.168.1.2")
    port = 49000
    xpc = XPlane12Output(; address, port)
    IODevices.init!(xpc)

    h_trn = HOrth(427.2);
    kin_init = KinInit(
        loc = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997)),
        q_nb = REuler(deg2rad(157), 0, 0),
        h = h_trn + C172.Δh_to_gnd + 0.5);

    kin_data = KinData(kin_init)

    msg_tuple = (
        Network.xpmsg_set_dref(override_pose, 1),
        Network.xpmsg_set_dref(override_surf, 1),
        Network.xpmsg_set_dref(override_prop, 1),
        Network.xpmsg_set_dref(override_nws, 0),
        Network.xpmsg_set_dref(lele_pos, rad2deg(δe)),
        Network.xpmsg_set_dref(rele_pos, rad2deg(δe)),
        Network.xpmsg_set_dref(lail_pos, rad2deg(δa)),
        Network.xpmsg_set_dref(rail_pos, rad2deg(-δa)),
        Network.xpmsg_set_dref(lflap_pos, rad2deg(δf)),
        Network.xpmsg_set_dref(rflap_pos, rad2deg(δf)),
        Network.xpmsg_set_dref(rud_pos, rad2deg(δr)),
        Network.xpmsg_set_dref(prop_is_disc, 0),
        Network.xpmsg_set_dref(prop_angle, rad2deg(ϕ_prop)),
        Network.xpmsg_set_dref(nws_angle, rad2deg(ψ_sw)),
        Network.xpmsg_set_pose(XPlanePose(kin_data)) #UDP message
    )
    IODevices.handle_data!(xpc, msg_tuple)

    IODevices.shutdown!(xpc)

end


end