using Test

using LinearAlgebra
using Flight


function test_θ_constraint()

    #precompute v_wOb_b
    α_a = 0.15
    β_a = -0.11
    TAS = 100
    v_wOa_a = Air.get_velocity_vector(TAS, α_a, β_a)
    v_wOb_b = v_wOa_a #in general we'd have v_wOb_b = q_ba(v_wOa_a)

    #set γ_wOb_n and φ_nb arbitrarily, compute θ_nb, construct e_nb, transform
    #v_wOb_b to v_wOb_n, compute γ_wOb_n and check it matches
    @show γ_wOb_n = -0.07 #set arbitrarily
    ψ_nb = 0.3 #does not matter
    φ_nb = 0.7
    @show θ_nb = θ_constraint(; v_wOb_b, γ_wOb_n, φ_nb)
    e_nb = REuler(ψ_nb, θ_nb, φ_nb)
    v_wOb_n = e_nb(v_wOb_b)
    @show γ_wOb_n_test = Attitude.inclination(v_wOb_n)
    @test γ_wOb_n_test ≈ γ_wOb_n
    @test @ballocated(θ_constraint(; v_wOb_b = $v_wOb_b, γ_wOb_n = $γ_wOb_n, φ_nb = $φ_nb)) === 0

end
