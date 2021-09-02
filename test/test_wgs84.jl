module TestWGS84

using Test
using LinearAlgebra
using Flight.WGS84, Flight.Attitude

export test_wgs84

function test_wgs84()
    @testset verbose = true "WGS84" begin
        @testset verbose = true "NVector" begin test_NVector() end
        @testset verbose = true "NVectorAlt" begin test_NVectorAlt() end
        @testset verbose = true "LatLonAlt" begin test_LatLonAlt() end
        @testset verbose = true "CartECEF" begin test_CartECEF() end
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

function test_NVectorAlt()

    v = [1, 2, -3]
    h = 1000
    n_e = NVector(v)

    #inner constructor
    p = NVectorAlt(n_e, h)
    @test p.n_e == n_e
    @test p.h == h

    #equality and approximate equality
    @test p == p
    @test p != NVectorAlt(-n_e, h)
    @test p ≈ NVectorAlt(n_e, h+1e-10)
    @test !(p ≈ NVectorAlt(-n_e, h))

    #inversion
    @test -p == NVectorAlt(-n_e, h)

end

function test_LatLonAlt()

    #keyword constructor
    llh = LatLonAlt(ϕ = π/6, λ = -π/4, h = 15092.0)
    nvh = NVectorAlt(llh)

    #approximate equality
    @test llh ≈ LatLonAlt(llh.ϕ+1e-12, llh.λ-1e-12, llh.h+1e-10)
    @test !(llh ≈ LatLonAlt())
    @test llh ≈ nvh #mixed types

    #inversion
    @test -llh ≈ -nvh

    #general functions
    @test g_l(llh) ≈ g_l(nvh)
    @test G_l(llh) ≈ G_l(nvh)
    @test ltf(llh) ≈ ltf(nvh)
    @test radii(llh) == radii(nvh)

    #conversion to and from NVectorAlt
    ϕ_range = range(-π/2, π/2, length = 10)
    λ_range = range(-π, π, length = 10)
    h_range = range(WGS84.h_min + 1, 10000, length = 10)

    llh_array = [LatLonAlt(ϕ, λ, h) for (ϕ, λ, h) in Iterators.product(ϕ_range, λ_range, h_range)]
    llh_array_test = [llh |> NVectorAlt |> LatLonAlt for llh in llh_array]
    @test all(llh_array .≈ llh_array_test)

end

function test_CartECEF()

    #keyword constructor
    llh = LatLonAlt(ϕ = π/6, λ = -π/4, h = 15092.0)
    cart = CartECEF(llh)
    nvh = NVectorAlt(cart)

    #approximate equality
    @test cart ≈ CartECEF(cart.data .+ 1e-9)
    @test !(cart ≈ LatLonAlt())
    @test cart ≈ nvh #mixed types
    @test cart ≈ llh #mixed types

    #inversion
    @test -llh ≈ -nvh
    @test -llh ≈ -cart

    #general functions
    @test g_l(cart) ≈ g_l(nvh)
    @test G_l(cart) ≈ G_l(nvh)
    @test ltf(cart) ≈ ltf(nvh)
    @test radii(cart) == radii(nvh)

    #conversion to and from NVectorAlt, LatLonAlt
    ϕ_range = range(-π/2, π/2, length = 10)
    λ_range = range(-π, π, length = 10)
    h_range = range(WGS84.h_min + 1, 10000, length = 10)

    cart_array = [CartECEF(LatLonAlt(ϕ, λ, h)) for (ϕ, λ, h) in Iterators.product(ϕ_range, λ_range, h_range)]
    cart_array_test = [cart |> NVectorAlt |> CartECEF |> LatLonAlt |> CartECEF for cart in cart_array]
    @test all(cart_array .≈ cart_array_test)

end

end