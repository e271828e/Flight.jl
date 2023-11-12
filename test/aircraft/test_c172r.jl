module TestC172R

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight.FlightCore.Systems
using Flight.FlightCore.Sim

using Flight.FlightPhysics.Environment

using Flight.FlightComponents.Aircraft

using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172R

export test_c172r

function test_c172r()
    @testset verbose = true "Cessna172R" begin
        test_trimming()
    end
end

function test_trimming()

    @testset verbose = true "Trimming" begin

        physics = System(C172R.Physics())
        env = System(SimpleEnvironment())
        trim_params = C172.TrimParameters()
        state = C172.TrimState()

        f_target = C172.get_f_target(physics, trim_params, env)

        @test @ballocated($f_target($state)) === 0

        success, _ = trim!(physics, trim_params, env)

        @test success

    end #testset

end #function

function test_linearization()

    @testset verbose = true "Linearization" begin

    end #testset

end #function

end #module