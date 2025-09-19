module TestC172Xv2

using Test, UnPack, BenchmarkTools, Sockets, JSON3, Logging

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

using Flight.FlightAircraft.C172: is_on_gnd
using Flight.FlightAircraft.C172X.C172XControl: ModeControlLon, ModeControlLat
using Flight.FlightAircraft.C172X.C172XGuidance: Segment, SegmentGuidanceData, ModeGuidance

export test_c172x2


y_kin(aircraft::Model{<:Cessna172Xv2}) = aircraft.y.vehicle.kinematics
y_air(aircraft::Model{<:Cessna172Xv2}) = aircraft.y.vehicle.airflow
y_aero(aircraft::Model{<:Cessna172Xv2}) = aircraft.y.vehicle.systems.aero


function test_c172x2(; alloc::Bool = true)

    @testset verbose = true "Cessna 172Xv2" begin

        @testset verbose = true "Segment" begin

            #default constructor should return a valid segment
            @test_nowarn Segment()
            #segments with zero horizontal length are invalid
            @test_throws ArgumentError Segment(Geographic(), Geographic())
            @test_throws ArgumentError Segment(Geographic(), Geographic(h = HEllip(100)))

            χ = π/3 #test segment's azimuth
            Δχ = π/4 #test point's segment-relative azimuth
            s = 1e3 #test point's distance along its azimuth
            seg = Segment(Geographic(); χ, s = 1e4, γ = deg2rad(5))
            p = Segment(Geographic(); χ = χ + Δχ, s, γ = 0).p2
            data = SegmentGuidanceData(seg, p)

            @test data.s_1b ≈ s * cos(Δχ) atol = 1e-2
            @test data.e_sb ≈ s * sin(Δχ) atol = 1e-2
            @test data.h_s ≈ data.s_1b * tand(5) atol = 1e-2

            seg_inv = -seg
            @test seg_inv.p1 ≈ seg.p2
            @test seg_inv.p2 ≈ seg.p1

            if alloc
                @test @ballocated(Segment(Geographic(); χ = 0, s = 1e4, γ = deg2rad(5))) == 0
                @test @ballocated(SegmentGuidanceData($seg, $p)) == 0
            end

        end #testset

        @testset verbose = true "SegmentGuidance" begin

            h_trn = HOrth(0.0)
            world = SimpleWorld(Cessna172Xv2(), SimpleAtmosphere(), HorizontalTerrain(h_trn)) |> Model

            aircraft = world.aircraft
            @unpack ctl, gdc = aircraft.avionics
            e_thr = gdc.seg.constants.e_thr

            init_air = C172.TrimParameters(; ψ_nb = 0)
            init_gnd = C172.Init(KinInit( h = h_trn + C172.Δh_to_gnd))

            dt = Δt = 0.01
            sim = Simulation(world; dt, Δt, t_end = 100)

            ############################## Gnd #################################

            init!(sim, init_gnd)
            @test is_on_gnd(aircraft.vehicle)

            #request segment guidance mode and enable horizontal and vertical guidance
            gdc.u.mode_req = ModeGuidance.segment
            gdc.seg.u.hor_gdc_req = true
            gdc.seg.u.vrt_gdc_req = true

            #step for one guidance sample period
            step!(sim, gdc.Δt, true)

            #the mode request should have been ignored due to wow = gnd
            @test gdc.y.mode === ModeGuidance.direct

            ############################## Air #################################

            init!(sim, init_air)

            kin_data = y_kin(aircraft)
            @unpack n_e, h_e, χ_gnd = kin_data
            Ob = Geographic(n_e, h_e)
            χ_ac = χ_gnd
            Δh = 100

            #construct and assign a segment roughly perpendicular to the current
            #course, e_thr/2 to the right and 100m above the current altitude
            aux_seg = Segment(Ob; χ = χ_ac + π/2, s = e_thr/2, Δh)
            gdc.seg.u.target = Segment(aux_seg.p2; χ = 0, s = 1e4, γ = deg2rad(5))

            #step for one guidance sample period
            Sim.step!(sim, gdc.Δt, true)

            #the guidance mode request should now have been honored
            @test gdc.y.mode === ModeGuidance.segment

            #horizontal guidance request is always honored
            @test gdc.seg.y.hor_gdc === true
            @test ctl.lat.y.mode === ModeControlLat.χ_β
            @test ctl.lat.y.χ_ref === gdc.seg.y.χ_ref

            #since we are within the vertical guidance threshold e_thr, vertical
            #guidance must have engaged as requested
            @test gdc.seg.y.vrt_gdc === true
            @test (ctl.lon.y.mode === ModeControlLon.EAS_alt ||
                   ctl.lon.y.mode === ModeControlLon.thr_EAS)
            @test ctl.lon.y.h_ref === gdc.seg.y.h_ref
            @test ctl.lon.y.h_ref ≈ h_e + Δh atol = 1

            #intercept angle should be positive
            @test gdc.seg.y.Δχ > 0

            #test with both horizontal and vertical guidance enabled
            if alloc
                @test @ballocated(f_ode!($world)) == 0
                @test @ballocated(f_step!($world)) == 0
                @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0
            end

            #now place the segment to the left of the aircraft, e_thr/2 away
            aux_seg = Segment(Ob; χ = χ_ac - π/2, s = e_thr/2, γ = 0)
            gdc.seg.u.target = Segment(aux_seg.p2; χ = 0, s = 1e4, γ = deg2rad(5))

            Sim.step!(sim, gdc.Δt, true)

            #intercept angle should be negative
            @test gdc.seg.y.Δχ < 0

            aux_seg = Segment(Ob; χ = χ_ac + π/2, s = 2e_thr, γ = 0)
            gdc.seg.u.target = Segment(aux_seg.p2; χ = 0, s = 1e4, γ = deg2rad(5))

            Sim.step!(sim, gdc.Δt, true)

            #since we are outside the vertical guidance threshold e_thr,
            #vertical guidance must have disengaged
            @test gdc.seg.y.vrt_gdc === false

            #disable vertical guidance, longitudinal control mode should remain
            #unchanged
            lon_mode_prev = ctl.lon.y.mode
            gdc.seg.u.vrt_gdc_req = false
            Sim.step!(sim, gdc.Δt, true)
            @test gdc.seg.y.vrt_gdc === false
            @test ctl.lon.y.mode === lon_mode_prev

            #once vertical guidance is disabled, we can once again set the
            #longitudinal control modes
            ctl.lon.u.mode_req = ModeControlLon.sas
            Sim.step!(sim, gdc.Δt, true)
            @test ctl.lon.y.mode === ModeControlLon.sas

            #disable horizontal guidance, lateral control mode should remain
            #unchanged
            lat_mode_prev = ctl.lat.y.mode
            gdc.seg.u.hor_gdc_req = false
            Sim.step!(sim, gdc.Δt, true)
            @test gdc.seg.y.hor_gdc === false
            @test ctl.lat.y.mode === lat_mode_prev

            #once horizontal guidance is disabled, we can once again set the
            #lateral control modes
            ctl.lat.u.mode_req = ModeControlLat.sas
            Sim.step!(sim, gdc.Δt, true)
            @test ctl.lat.y.mode === ModeControlLat.sas

        end #testset

    end #testset

end #function

end #module