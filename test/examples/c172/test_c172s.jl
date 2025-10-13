module TestC172S

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightExamples

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
        atmosphere = Model(SimpleAtmosphere())
        terrain = Model(HorizontalTerrain())
        params = C172.TrimParameters()
        state = C172.TrimState()

        f_target = C172.get_f_target(vehicle, params, atmosphere, terrain)

        @test @ballocated($f_target($state)) === 0

        success, _ = Modeling.init!(vehicle, params)

        @test success

    end #testset

end #function

function test_linearization()

    @testset verbose = true "Linearization" begin
        @test_nowarn Cessna172Sv0(NED()) |> Model |> linearize!
    end #testset

end #function

function test_update_methods()

    @testset verbose = true "Update Methods" begin

        atmosphere = SimpleAtmosphere() |> Model
        terrain = HorizontalTerrain() |> Model

        location = NVector()
        trn_data = TerrainData(terrain, location)
        vehicle_init = KinInit( h = trn_data.elevation + 1.8) |> C172.Init

        aircraft = Model(Cessna172Sv0());

        Modeling.init!(aircraft, vehicle_init, atmosphere, terrain)

        #ensure we are on the ground for full landing gear code coverage
        @test aircraft.y.vehicle.systems.ldg.left.strut.wow == true

        @test @ballocated(f_ode!($aircraft, $atmosphere, $terrain)) == 0
        @test @ballocated(f_step!($aircraft, $atmosphere, $terrain)) == 0
        @test @ballocated(f_periodic!(NoScheduling(), $aircraft, $atmosphere, $terrain)) == 0

    end

    return nothing

end


end #module