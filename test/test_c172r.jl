module TestC172R

using Test
using UnPack
using BenchmarkTools

using Flight

export test_c172r

function test_c172r()
    @testset verbose = true "Cessna172R" begin

        @testset verbose = true "Performance" begin test_system() end
        @testset verbose = true "Simulation" begin test_sim_nrt(save = false) end

    end
end

function test_system()

        trn = HorizontalTerrain();
        atm = System(Atmosphere());

        ac_LTF = System(Cessna172R(LTF()));
        ac_ECEF = System(Cessna172R(ECEF()));
        ac_NED = System(Cessna172R(NED()));

        #all three kinematics implementations must be supported, no allocations
        @test @ballocated(f_cont!($ac_LTF, $atm, $trn)) == 0
        @test @ballocated(f_disc!($ac_LTF)) == 0

        @test @ballocated(f_cont!($ac_ECEF, $atm, $trn)) == 0
        @test @ballocated(f_disc!($ac_ECEF)) == 0

        @test @ballocated(f_cont!($ac_NED, $atm, $trn)) == 0
        @test @ballocated(f_disc!($ac_NED)) == 0

    return nothing

end

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

function test_sim_nrt(; save::Bool = true)

    h_trn = AltO(608.55);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(Atmosphere());
    ac = System(Cessna172R());
    kin_init = KinematicInit(
        v_eOb_n = [30, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        l2d = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 2200.5);

    Aircraft.init!(ac, kin_init)
    ac.u.avionics.eng_start = true #engine start switch on

    atm.u.wind.v_ew_n[1] = 0

    callback! = let

        function (u, t, y, params)

            ac.u.avionics.throttle = 0.2
            ac.u.avionics.yoke_Δx = (t < 5 ? 0.2 : 0.0)
            ac.u.avionics.yoke_Δy = 0.0
            ac.u.avionics.pedals = 0.0
            ac.u.avionics.brake_left = 0
            ac.u.avionics.brake_right = 0

        end
    end

    sim = Simulation(ac; args_c = (atm, trn), t_end = 150, sim_callback = callback!)

    Sim.run!(sim)
    plots = make_plots(sim; Plotting.defaults...)
    save ? save_plots(plots, save_folder = joinpath("tmp", "nrt_sim_test")) : nothing

    return sim

end

function test_sim_rt(; save::Bool = true)

    h_trn = AltO(608.55);

    trn = HorizontalTerrain(altitude = h_trn);
    atm = System(Atmosphere());
    ac = System(Cessna172R());
    kin_init = KinematicInit(
        v_eOb_n = [30, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        l2d = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.9 + 100.5);

    Aircraft.init!(ac, kin_init)
    ac.u.avionics.eng_start = true #engine start switch on
    ac.u.avionics.throttle = 1

    callback! = let input_callback! = get_input_callback(), output_callback = get_output_callback()

        function (u, t, y, params)
            input_callback!(u)
            output_callback(y)
        end
    end

    sim = Simulation(ac;
        args_c = (atm, trn),
        t_end = 120,
        sim_callback = callback!,
        realtime = true,
        )


    Sim.run!(sim)
    plots = make_plots(sim; Plotting.defaults...)
    save ? save_plots(plots, save_folder = joinpath("tmp", "rt_sim_test")) : nothing

    return sim

end


end #module