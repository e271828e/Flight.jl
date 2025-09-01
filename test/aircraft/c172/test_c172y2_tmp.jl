module TestC172Yv2

using Test, UnPack, BenchmarkTools, Sockets, JSON3

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

#non-exported stuff
using Flight.FlightLib.Control.Discrete: load_pid_lookup, load_lqr_tracker_lookup
using Flight.FlightAircraft.C172Y.C172YControl: ModeControlLon, ModeControlLat
using Flight.FlightAircraft.C172Y.C172YControl: ModeGuidanceLon, ModeGuidanceLat
using Flight.FlightAircraft.C172Y.C172YControl: FlightPhase

export test_c172y2


function test_c172y2()
    @testset verbose = true "Cessna 172Yv2" begin

        test_control_modes()
        test_guidance_modes()

    end
end

y_kin(aircraft::Model{<:Cessna172Yv2}) = aircraft.y.vehicle.kinematics
y_air(aircraft::Model{<:Cessna172Yv2}) = aircraft.y.vehicle.airflow
y_aero(aircraft::Model{<:Cessna172Yv2}) = aircraft.y.vehicle.systems.aero

function test_guidance_modes()

    @testset verbose = true "Guidance Modes" begin

    h_trn = HOrth(0.0)
    world = SimpleWorld(Cessna172Yv2(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

    aircraft = world.aircraft
    ctl = aircraft.avionics.ctl

    init_air = C172.TrimParameters()
    dt = Δt = 0.01

    sim = Simulation(world; dt, Δt, t_end = 600)

    @testset verbose = true "ModeGuidanceLon.alt" begin

        Sim.init!(sim, init_air)
        y_kin_trim = y_kin(aircraft)

        ctl.u.mode_gdc_lon_req = ModeGuidanceLon.alt
        ctl.u.mode_ctl_lat_req = ModeControlLat.φ_β
        step!(sim, ctl.Δt, true)
        @test ctl.y.mode_gdc_lon === ModeGuidanceLon.alt
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm

        #when trim reference values are kept, the guidance mode must activate without
        #transients
        step!(sim, 1, true)
        @test all(isapprox.(y_kin(aircraft).ω_wb_b[2], y_kin_trim.ω_wb_b[2]; atol = 1e-5))
        @test all(isapprox.(y_kin(aircraft).v_eb_b[1], y_kin_trim.v_eb_b[1]; atol = 1e-2))

        #all tests while turning
        ctl.u.φ_ref = π/12

        ctl.u.h_target = y_kin_trim.h_e + 100
        step!(sim, 1, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.thr_EAS
        step!(sim, 60, true) #altitude is captured
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm
        @test isapprox.(y_kin(aircraft).h_e - HEllip(ctl.u.h_target), 0.0; atol = 1e-1)

        #reference changes within the current threshold do not prompt a mode change
        ctl.u.h_target = y_kin(aircraft).h_e - ctl.gdc_lon_alt.constants.h_thr / 2
        step!(sim, 1, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm
        step!(sim, 30, true) #altitude is captured
        @test isapprox.(y_kin(aircraft).h_e - HEllip(ctl.u.h_target), 0.0; atol = 1e-1)

        ctl.u.h_target = y_kin_trim.h_e - 100
        step!(sim, 1, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.thr_EAS
        step!(sim, 80, true) #altitude is captured
        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm
        @test isapprox.(y_kin(aircraft).h_e - HEllip(ctl.u.h_target), 0.0; atol = 1e-1)

        @test ctl.y.mode_ctl_lon === ModeControlLon.EAS_clm

        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

        ctl.u.h_target = y_kin_trim.h_e + 100
        step!(sim, 1, true)
        @test ctl.y.mode_ctl_lon === ModeControlLon.thr_EAS

        #test for allocations in the current control mode
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

        # return TimeSeries(sim)

    end

    end #testset

end


end #module