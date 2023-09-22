module TestC172FBW

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight.FlightCore.Systems
using Flight.FlightCore.Sim

using Flight.FlightPhysics.Environment

using Flight.FlightComponents.Aircraft

using Flight.FlightAircraft.C172FBW

export test_c172fbw

function test_c172fbw()
    @testset verbose = true "Cessna172FBW" begin
        test_trimming()
    end
end

function test_trimming()

    @testset verbose = true "Trimming" begin

        ac = System(Cessna172FBWBase())
        env = System(SimpleEnvironment())
        trim_params = C172FBW.TrimParameters()
        state = C172FBW.TrimState()

        f_target = C172FBW.get_f_target(ac, env, trim_params)

        @test @ballocated($f_target($state)) === 0

        success, _ = trim!(ac; env, trim_params)

        @test success

    end #testset

end #function

end #module