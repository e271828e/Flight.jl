module TestC172S

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

export test_c172s


function test_c172s()
    @testset verbose = true "Cessna 172S" begin

        test_trimming()
        test_linearization()
        test_update_methods()

    end
end

function test_trimming()

    @testset verbose = true "Trimming" begin

        vehicle = Model(C172S.Vehicle())
        atm = Model(SimpleAtmosphere())
        trn = Model(HorizontalTerrain())
        trim_params = C172.TrimParameters()
        state = C172.TrimState()

        f_target = C172.get_f_target(vehicle, trim_params, atm, trn)

        @test @ballocated($f_target($state)) === 0

        success, _ = Modeling.init!(vehicle, trim_params)

        @test success

    end #testset

end #function

function test_linearization()

    @testset verbose = true "Linearization" begin
        ss = Cessna172Sv0(NED()) |> Model |> Control.Continuous.LinearizedSS
    end #testset

end #function

function test_update_methods()

    @testset verbose = true "Update Methods" begin

        atm = SimpleAtmosphere() |> Model
        trn = HorizontalTerrain() |> Model

        loc = NVector()
        trn_data = TerrainData(trn, loc)
        vehicle_init = KinInit( h = trn_data.elevation + 1.8) |> C172.Init

        ac = Model(Cessna172Sv0());

        Modeling.init!(ac, vehicle_init, atm, trn)

        #ensure we are on the ground for full landing gear code coverage
        @test ac.y.vehicle.systems.ldg.left.strut.wow == true

        @test @ballocated(f_ode!($ac, $atm, $trn)) == 0
        @test @ballocated(f_disc!($ac, $atm, $trn)) == 0
        @test @ballocated(f_step!($ac, $atm, $trn)) == 0

    end

    return nothing

end


end #module