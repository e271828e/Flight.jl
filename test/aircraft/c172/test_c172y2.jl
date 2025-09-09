module TestC172Yv2

using Test, UnPack, BenchmarkTools, Sockets, JSON3, Logging

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

using Flight.FlightAircraft.C172Y.C172YControl: ModeControlLon, ModeControlLat,
        AltTrackingState, is_on_gnd
using Flight.FlightAircraft.C172Y.C172YGuidance: Segment, SegmentCoords, ModeGuidance

export test_c172y2

function test_c172y2(; alloc::Bool = true)

    @testset verbose = true "Cessna 172Yv2" begin

        @testset verbose = true "Segment" begin

            #default constructor should return a valid segment
            @test_nowarn Segment()
            #segments with zero horizontal length are invalid
            @test_throws AssertionError Segment(Geographic(), Geographic())
            @test_throws AssertionError Segment(Geographic(), Geographic(h = HEllip(100)))

            χ = π/3 #test segment's azimuth
            Δχ = π/4 #test point's segment-relative azimuth
            l = 1000 #test point's distance along its azimuth
            seg = Segment(Geographic(); χ, l = 10000, γ = deg2rad(5))
            p = Segment(Geographic(); χ = χ + Δχ, l, γ = 0).p2
            coords = SegmentCoords(seg, p)

            @test coords.l_1b ≈ l * cos(Δχ) atol = 1e-2
            @test coords.e_1b ≈ l * sin(Δχ) atol = 1e-2
            @test coords.h_1b ≈ coords.l_1b * tand(5) atol = 1e-2

            seg_inv = -seg
            @test seg_inv.p1 ≈ seg.p2
            @test seg_inv.p2 ≈ seg.p1

            if alloc
                @test @ballocated(Segment(Geographic(); χ = 0, l = 10000, γ = deg2rad(5))) == 0
                @test @ballocated(SegmentCoords($seg, $p)) == 0
            end

        end #testset

        @testset verbose = true "SegmentGuidance" begin

            h_trn = HOrth(0.0)
            world = SimpleWorld(Cessna172Yv2(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

            aircraft = world.aircraft
            gdc = aircraft.avionics
            init_air = C172.TrimParameters()
            init_gnd = C172.Init(KinInit( h = h_trn + C172.Δh_to_gnd))

            dt = Δt = 0.01
            sim = Simulation(world; dt, Δt, t_end = 100)

            Sim.init!(sim, init_gnd)
            @test is_on_gnd(aircraft.vehicle)

            #request segment guidance
            gdc.u.mode_req = ModeGuidance.segment

            #step for one controller sample period
            step!(sim, Δt, true)

            #the mode request should have been ignored due to wow = gnd
            @test gdc.y.mode === ModeGuidance.direct

            # if alloc
            #     @test @ballocated(f_ode!($world)) == 0
            #     @test @ballocated(f_step!($world)) == 0
            #     @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0
            # end

        end #testset

    end #testset


end

y_kin(aircraft::Model{<:Cessna172Yv2}) = aircraft.y.vehicle.kinematics
y_air(aircraft::Model{<:Cessna172Yv2}) = aircraft.y.vehicle.airflow
y_aero(aircraft::Model{<:Cessna172Yv2}) = aircraft.y.vehicle.systems.aero


end #module