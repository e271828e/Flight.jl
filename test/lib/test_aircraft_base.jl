module TestAircraftBase

using Test
using BenchmarkTools

using Flight.FlightLib

export test_aircraft_base

function test_aircraft_base()
    @testset verbose = true "Aircraft" begin
        test_θ_constraint()
    end
end

function test_θ_constraint()

    @testset verbose = true "θ Constraint" begin

        #precompute v_wb_b
        α_a = 0.15
        β_a = -0.11
        TAS = 100
        v_wOa_a = Atmosphere.get_velocity_vector(TAS, α_a, β_a)
        v_wb_b = v_wOa_a

        #set γ_wb_n and φ_nb arbitrarily and compute θ_nb
        γ_wb_n = -0.07 #set arbitrarily
        ψ_nb = 0.3 #doesn't really matter
        φ_nb = 0.7
        θ_nb = AircraftBase.θ_constraint(; v_wb_b, γ_wb_n, φ_nb)

        #then construct e_nb, transform v_wb_b to v_wb_n, recompute γ_wb_n
        #and check it matches the original value
        e_nb = REuler(ψ_nb, θ_nb, φ_nb)
        v_wb_n = e_nb(v_wb_b)
        γ_wb_n_test = inclination(v_wb_n)

        @test γ_wb_n_test ≈ γ_wb_n
        @test @ballocated(AircraftBase.θ_constraint(; v_wb_b = $v_wb_b, γ_wb_n = $γ_wb_n, φ_nb = $φ_nb)) === 0

    end


end #function

end #module