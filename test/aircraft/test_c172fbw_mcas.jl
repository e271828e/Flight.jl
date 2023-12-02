module TestC172MCAS

using Test, UnPack, BenchmarkTools, Sockets

using Flight.FlightCore
using Flight.FlightCore.Sim
using Flight.FlightCore.Visualization

using Flight.FlightPhysics
using Flight.FlightComponents

using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172MCAS

export test_c172_mcas


function test_c172_mcas()
    @testset verbose = true "Cessna172 MCAS" begin

        test_system_methods()
        test_sim(save = false)

    end
end

function test_c172mcas_modes()

    trn = HorizontalTerrain()
    loc = NVector()
    trn_data = TerrainData(trn, loc)
    kin_init_gnd = KinematicInit( h = trn_data.altitude + 1.8);
    kin_init_air = KinematicInit( h = trn_data.altitude + 1000);

    ac = System(Cessna172MCS());

    #to exercise all airframe functionality, including landing gear, we
    #need to be on the ground with the engine running
    init_kinematics!(ac, kin_init_gnd)

    sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 600)

    #start on the ground, assign different navigation and control modes, and
    #make sure the system remains in lon_SAS_off, lat_SAS_off, and that surface
    #setpoints pass through directly to the actuators


end

function test_system_methods()

        @testset verbose = true "System Methods" begin

            trn = HorizontalTerrain()
            loc = NVector()
            trn_data = TerrainData(trn, loc)
            kin_init_gnd = KinematicInit( h = trn_data.altitude + 1.8);
            kin_init_air = KinematicInit( h = trn_data.altitude + 1000);

            ac = System(Cessna172MCAS());

            #to exercise all airframe functionality, including landing gear, we
            #need to be on the ground with the engine running
            init_kinematics!(ac, kin_init_gnd)
            ac.avionics.u.inceptors.eng_start = true #engine start switch on
            f_disc!(ac, 0.02, env) #run avionics for the engine start signal to propagate
            f_ode!(ac, env)
            f_step!(ac)
            f_ode!(ac, env)
            f_step!(ac)
            @test ac.y.physics.airframe.ldg.left.strut.wow == true
            @test ac.y.physics.airframe.pwp.engine.state === Piston.eng_starting

            @test @ballocated(f_ode!($ac)) == 0
            @test @ballocated(f_step!($ac)) == 0
            @test @ballocated(f_disc!($ac, 0.02)) == 0

            #now we put the aircraft in flight
            init_kinematics!(ac, kin_init_air)
            f_ode!(ac, env)
            @test ac.y.physics.airframe.ldg.left.strut.wow == false
            @test @ballocated(f_ode!($ac)) == 0
            @test @ballocated(f_step!($ac)) == 0

            #testing the different avionics modes for allocations is a bit more
            #involved
        end

    return nothing

end

function test_cas(; save::Bool = true)

    @testset verbose = true "Simulation" begin

        ac = Cessna172MCAS() |> System;
        design_condition = C172.TrimParameters()

        exit_flag, trim_state = trim!(ac, design_condition)

        @test exit_flag === true

        sys_io! = let

            function (ac)

                t = ac.t[]

                u_inceptors = ac.avionics.u.inceptors
                u_digital = ac.avionics.u.digital

                u_digital.lon_mode_sel = C172MCAS.lon_θ_EAS

                # u_digital.lon_mode_sel = C172MCAS.lon_q_EAS
                # if 5 < t < 15
                #     u_inceptors.pitch_input = 0.001
                # elseif 15 < t < 25
                #     u_inceptors.pitch_input = -0.001
                # else
                #     u_inceptors.pitch_input = 0
                # end

            end
        end

        sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 30, sys_io!, adaptive = false)
        # sim = Simulation(ac; dt = 0.01, Δt = 0.01, t_end = 60, adaptive = false)
        Sim.run!(sim, verbose = true)

        # plots = make_plots(sim; Plotting.defaults...)
        kin_plots = make_plots(TimeHistory(sim).physics.kinematics; Plotting.defaults...)
        air_plots = make_plots(TimeHistory(sim).physics.air; Plotting.defaults...)
        rb_plots = make_plots(TimeHistory(sim).physics.rigidbody; Plotting.defaults...)
        save && save_plots(kin_plots, save_folder = joinpath("tmp", "test_c172_mcas", "cas", "kin"))
        save && save_plots(air_plots, save_folder = joinpath("tmp", "test_c172_mcas", "cas", "air"))
        save && save_plots(rb_plots, save_folder = joinpath("tmp", "test_c172_mcas", "cas", "rigidbody"))

        return nothing

    end
end




end #module