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

    h_trn = AltO(608.55);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(Atmosphere());
    ac = System(C172RAircraft());
    kin_init = KinInit(
        v_eOb_n = [30, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        Ob = GeographicLocation(
            LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
            h_trn + 1.9 + 2200.5));

    Aircraft.init!(ac, kin_init)

    atm.u.wind.v_ew_n[1] = 0

    callback! = let

        function (u, t, y, params)

            ac.u.avionics.yoke_Δx = (t < 5 ? 0.2 : 0.0)
            ac.u.avionics.yoke_Δy = 0.0
            ac.u.avionics.pedals = 0.0
            ac.u.avionics.brake_left = 0
            ac.u.avionics.brake_right = 0
            ac.u.avionics.throttle = 0


        end
    end

    sim = Simulation(ac; args_c = (atm, trn), t_end = 150, sim_callback = callback!)

    # @show @ballocated SciMLBase.step!($sim) #this takes many steps!

    Sim.run!(sim)
    plots = make_plots(sim)
    save_plots(plots, save_folder = joinpath("tmp", "nrt_sim_test"))
    # save_plots(plots, save_folder = joinpath("tmp", "sim_test", Dates.format(now(), "yyyy_mm_dd_HHMMSS")))

    return sim

end

function test_sim_rt()

    h_trn = AltO(608.55);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(Atmosphere());
    ac = System(C172RAircraft());
    kin_init = KinInit(
        v_eOb_n = [30, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        Ob = GeographicLocation(
            LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
            h_trn + 1.9 + 2200.5));

    Aircraft.init!(ac, kin_init)

    callback! = let input_callback! = get_input_callback(), output_callback = get_output_callback()

        function (u, t, y, params)
            input_callback!(u)
            output_callback(y)
        end
    end

    sim = Simulation(ac;
        args_c = (atm, trn),
        t_end = 10,
        sim_callback = callback!,
        realtime = true,
        )


    Sim.run!(sim)
    plots = make_plots(sim)
    save_plots(plots, save_folder = joinpath("tmp", "rt_sim_test"))
    # save_plots(plots, save_folder = joinpath("tmp", "sim_test", Dates.format(now(), "yyyy_mm_dd_HHMMSS")))

    return sim

end


end #module