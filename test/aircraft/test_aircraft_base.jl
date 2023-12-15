module TestAircraft

using Test
using UnPack
using BenchmarkTools

using Flight.FlightPhysics
using Flight.FlightAircraft.AircraftBase

export test_aircraft_base

function test_aircraft_base()
    @testset verbose = true "Aircraft" begin
        test_θ_constraint()
    end
end

function test_θ_constraint()

    @testset verbose = true "θ Constraint" begin

        #precompute v_wOb_b
        α_a = 0.15
        β_a = -0.11
        TAS = 100
        v_wOa_a = Atmosphere.get_velocity_vector(TAS, α_a, β_a)
        v_wOb_b = v_wOa_a

        #set γ_wOb_n and φ_nb arbitrarily and compute θ_nb
        γ_wOb_n = -0.07 #set arbitrarily
        ψ_nb = 0.3 #inconsequential
        φ_nb = 0.7
        θ_nb = AircraftBase.θ_constraint(; v_wOb_b, γ_wOb_n, φ_nb)

        #then construct e_nb, transform v_wOb_b to v_wOb_n, recompute γ_wOb_n
        #and check it matches the original value
        e_nb = REuler(ψ_nb, θ_nb, φ_nb)
        v_wOb_n = e_nb(v_wOb_b)
        γ_wOb_n_test = inclination(v_wOb_n)

        @test γ_wOb_n_test ≈ γ_wOb_n
        @test @ballocated(AircraftBase.θ_constraint(; v_wOb_b = $v_wOb_b, γ_wOb_n = $γ_wOb_n, φ_nb = $φ_nb)) === 0

    end


end #function

end #module