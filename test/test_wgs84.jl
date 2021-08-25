module TestWGS84

using Test
using LinearAlgebra
using Flight.WGS84, Flight.Rotations

export test_wgs84

function test_wgs84()
    @testset verbose = true "WGS84" begin
        @testset verbose = true "NVector" begin test_NVector() end
        @testset verbose = true "WGS84Pos" begin test_WGS84Pos() end
    end
end

function test_NVector()

    #inner constructor
    v = [1, 2, -3]
    n_e = NVector(v)
    @test n_e.data ≈ normalize(v) #normalizes by default

    #keyword constructor and latitude & longitude conversions
    @test NVector(; ϕ = 1.2).ϕ ≈ 1.2
    @test NVector(; ϕ = 1.2).λ == 0
    @test NVector(; λ = 0.2).ϕ == 0
    @test NVector(; λ = 0.2).λ ≈ 0.2
    @test NVector(; ϕ = .3, λ = -1.5).ϕ ≈ .3
    @test NVector(; ϕ = .3, λ = -1.5).λ ≈ -1.5

    #LTF constructors and LTF conversion
    r_en = WGS84.ltf(n_e, 0)
    @test r_en == n_e.ltf
    @test WGS84.ψ_nl(r_en) ≈ 0

    ψ_nl = π/3
    r_nl = Rz(ψ_nl)
    r_el = r_en ∘ r_nl
    @test WGS84.ψ_nl(r_el) ≈ ψ_nl
    @test WGS84.ψ_nl(RMatrix(r_el)) ≈ ψ_nl
    @test WGS84.ψ_nl(RAxAng(r_el)) ≈ ψ_nl

    #r_en and r_el should return the same NVector
    @test NVector(r_el) ≈ n_e
    @test NVector(RMatrix(r_el)) ≈ n_e
    @test NVector(RAxAng(r_el)) ≈ n_e

    #equality, unary minus and isapprox
    @test (-n_e).data == -(n_e.data)
    @test n_e == n_e
    @test n_e != -n_e
    @test n_e ≈ NVector(v + [0, 0, 1e-10])

    #torture test: from lat, lon to NVector, to r_en, then r_el, back to NVector
    ϕ_range = range(-π/2, π/2, length = 10)
    λ_range = range(-π, π, length = 10)
    ψ_nl_range = range(-π, π, length = 10)

    n_array = [NVector(; ϕ, λ) for (ϕ, λ, _) in Iterators.product(ϕ_range, λ_range, ψ_nl_range)]
    ψ_array = [ψ_nl for (_, _, ψ_nl) in Iterators.product(ϕ_range, λ_range, ψ_nl_range)]
    n_array_test = [NVector(NVector(; n0.ϕ, n0.λ).ltf ∘ Rz(ψ_nl)) for (n0, ψ_nl) in zip(n_array, ψ_array)]

    @test all(n_array .≈ n_array_test)

    @test WGS84.radii(NVector(ϕ = 0)).N ≈ WGS84.a
    @test WGS84.radii(NVector(ϕ = π/2)).M ≈ WGS84.radii(NVector(ϕ = π/2)).N

end

function test_WGS84Pos()

    v = [1, 2, -3]
    h = 1000
    n_e = NVector(v)

    #inner constructor
    p = WGS84Pos(n_e, h)
    @test p.n_e == n_e
    @test p.h == h

    #equality and isapprox
    @test p == p
    @test p != WGS84Pos(-n_e, h)
    @test p ≈ WGS84Pos(n_e, h+1e-10)
    @test !(p ≈ WGS84Pos(-n_e, h))

    # #outer constructors
    @test p ≈ WGS84Pos(; n_e.ϕ, n_e.λ, h)

    ϕ_range = range(-π/2, π/2, length = 10)
    λ_range = range(-π, π, length = 10)
    h_range = range(-1000, 10000, length = 10)

    p_array = [WGS84Pos(NVector(; ϕ, λ), h) for (ϕ, λ, h) in Iterators.product(ϕ_range, λ_range, h_range)]
    p_array_test = [WGS84Pos(rECEF(p)) for p in p_array]

    @test all(p_array .≈ p_array_test)

end



end