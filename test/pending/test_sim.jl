module TestSim

using Test
using UnPack
using BenchmarkTools

using Flight

export test_sim

function get_input_callback()

    inputs = init_joysticks() |> values |> collect
    # inputs = [XBoxController(),]

    let inputs = inputs

        function (u) #anonymous
            for input in inputs
                Input.update!(input)
                Input.assign!(u.avionics, input)
            end
        end

    end

end

function get_output_callback()

    # outputs = [XPInterface(host = IPv4("192.168.1.2"))] #Parsec
    outputs = [XPInterface(),]
    Output.init!.(outputs)

    let outputs = outputs

        function (y) #anonymous
            for output in outputs
                Output.update!(output, y.kinematics.pos)
            end
        end

    end

end

function test_sim_nrt()

    h_trn = AltOrth(608.55);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(AtmosphereDescriptor());
    ac = System(C172RAircraft());
    kin_init = KinInit(
        v_eOb_n = [30, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        Ob = Geographic(
            LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
            h_trn + 1.9 + 2200.5));

    Aircraft.init!(ac, kin_init)

    atm.u.wind.v_ew_n[1] = 0

    sim_callback! = let

        function (u, t, y, params)

            # ac.u.avionics.pedals = -.05
            # ac.u.avionics.yoke_Δx = 0.1
            ac.u.avionics.yoke_Δy = 0.0
            ac.u.avionics.brake_left = 0
            ac.u.avionics.brake_right = 0
            ac.u.avionics.throttle = 0

            ac.u.avionics.yoke_Δx = (t < 5 ? 0.1 : 0.0)

        end
    end

    sim = Simulation(ac; args_c = (trn, atm), t_end = 30, callback = sim_callback!)

    @show @ballocated Sim.step!($sim) #this takes many steps!

    # Sim.run!(sim)
    # plots = make_plots(sim)
    # save_plots(plots, save_folder = joinpath("tmp", "nrt_sim_test"))
    # save_plots(plots, save_folder = joinpath("tmp", "sim_test", Dates.format(now(), "yyyy_mm_dd_HHMMSS")))

    return sim

end

function test_sim_rt()

    h_trn = AltOrth(608.55);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(AtmosphereDescriptor());
    ac = System(C172RAircraft());
    kin_init = KinInit(
        v_eOb_n = [30, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        Ob = Geographic(
            LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
            h_trn + 1.9 + 2200.5));

    Aircraft.init!(ac, kin_init)

    sim_callback! = let input_callback! = get_input_callback(), output_callback = get_output_callback()

        function (u, t, y, params)
            input_callback!(u)
            output_callback(y)
        end
    end

    sim = Simulation(ac;
        args_c = (trn, atm),
        t_end = 10,
        callback = sim_callback!,
        realtime = true,
        )


    Sim.run!(sim)
    plots = make_plots(sim)
    save_plots(plots, save_folder = joinpath("tmp", "rt_sim_test"))
    # save_plots(plots, save_folder = joinpath("tmp", "sim_test", Dates.format(now(), "yyyy_mm_dd_HHMMSS")))

    return sim

end


end #module