module TestWorld

using Test
using BenchmarkTools

using Flight.FlightCore
using Flight.FlightLib

export test_world

function test_world()
    @testset verbose = true "Simple World" begin
        test_update_methods()
    end
end

function test_update_methods()

    @testset verbose = true "Update Methods" begin

        world = SimpleWorld() |> Model;

        @test @ballocated(f_ode!($world)) == 0
        @test @ballocated(f_step!($world)) == 0
        @test @ballocated(f_periodic!(NoScheduling(), $world)) == 0

    end

    return nothing

end #function

end #module