module TestPlotting

using Test
using Plots
using FlightCore.Plotting

export test_plotting

function test_plotting()
    @testset verbose = true "Plotting" begin

        #extension loads correctly
        @testset "Initialization" begin
            @test !isempty(Plotting.defaults)
        end

        #actual tests needed

    end
end

end # module
