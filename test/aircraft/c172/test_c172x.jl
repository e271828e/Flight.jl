module TestC172X

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

export test_c172x

function test_c172x()
    @testset verbose = true "Cessna 172X" begin
        test_trimming()
        test_linearization()
        test_update_methods()
    end
end

function test_trimming()

    @testset verbose = true "Trimming" begin

        #test on vehicle
        vehicle = System(C172X.Vehicle())
        atm = SimpleAtmosphere() |> System
        trn = HorizontalTerrain() |> System
        trim_params = C172.TrimParameters()
        state = C172.TrimState()

        f_target = C172.get_f_target(vehicle, trim_params, atm, trn)

        @test @ballocated($f_target($state)) === 0
        success, _ = Systems.init!(vehicle, trim_params)
        @test success

    end #testset

end #function

function test_linearization()

    @testset verbose = true "Linearization" begin

        ss = Cessna172Xv0(NED()) |> System |> Control.Continuous.LinearizedSS

    end #testset

end #function


function test_update_methods()

    @testset verbose = true "Update Methods" begin

        atm = SimpleAtmosphere() |> System
        trn = HorizontalTerrain() |> System

        loc = NVector()
        trn_data = TerrainData(trn, loc)
        vehicle_init = KinInit( h = trn_data.altitude + 1.8) |> C172.Init

        ac = System(Cessna172Xv0());

        Systems.init!(ac, vehicle_init, atm, trn)

        #ensure we are on the ground for full landing gear code coverage
        @test ac.y.vehicle.components.ldg.left.strut.wow == true

        @test @ballocated(f_ode!($ac, $atm, $trn)) == 0
        @test @ballocated(f_disc!($ac, $atm, $trn)) == 0
        @test @ballocated(f_step!($ac, $atm, $trn)) == 0

    end

    return nothing

end

end #module