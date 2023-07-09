module TestGeodesy

using Test
using LinearAlgebra

using Flight.FlightPhysics.Attitude
using Flight.FlightPhysics.Geodesy

export test_geodesy

function test_geodesy()
    @testset verbose = true "Geodesy" begin
        @testset verbose = true "NVector" begin test_NVector() end
        @testset verbose = true "LatLon" begin test_LatLon() end
        @testset verbose = true "Altitude" begin test_Altitude() end
        @testset verbose = true "Geographic" begin test_GeographicLocation() end
        @testset verbose = true "Cartesian" begin test_CartesianLocation() end
    end
end

function test_NVector()

    @test NVector() isa NVector #no-argument constructor

    n_e = NVector([1,2,3]) #normalize by default
    @test norm(n_e[:]) ≈ 1
    @test !(norm(NVector([1,2,3], normalization = false)[:]) ≈ 1) #bypass normalization

    #test basic operators
    @test n_e == n_e
    @test n_e ≈ n_e
    @test n_e != -n_e
    @test !(n_e ≈ -n_e)

    q_el = RQuat([1,2,3,4])
    m_el = RMatrix(q_el)
    a_el = RAxAng(q_el)

    #construction from ECEF to LTF rotation
    @test NVector(q_el) ≈ NVector(m_el) ≈ NVector(a_el)

    #LTF constructors and LTF conversion
    r_en = ltf(n_e, 0)
    @test r_en == ltf(n_e)
    @test get_ψ_nl(r_en) ≈ 0

    ψ_nl = π/3
    r_nl = Rz(ψ_nl)
    r_el = r_en ∘ r_nl
    @test get_ψ_nl(r_el) ≈ ψ_nl
    @test get_ψ_nl(RMatrix(r_el)) ≈ ψ_nl
    @test get_ψ_nl(RAxAng(r_el)) ≈ ψ_nl

    #r_en and r_el should return the same NVector
    @test NVector(r_el) ≈ n_e
    @test NVector(RMatrix(r_el)) ≈ n_e
    @test NVector(RAxAng(r_el)) ≈ n_e

    #radii
    @test radii(n_e) isa NamedTuple

end

function test_LatLon()

    @test LatLon() isa LatLon
    @test LatLon(ϕ = π/4, λ = -π/3) isa LatLon
    @test_throws ArgumentError LatLon(ϕ = 3π/2)
    @test_throws ArgumentError LatLon(λ = -5π/2)

    #test basic operators
    latlon = LatLon(ϕ = π/3, λ = -π/6)
    @test latlon ≈ latlon
    @test !(latlon ≈ -latlon)

    #ltf
    @test ltf(latlon) ≈ ltf(NVector(latlon))
    ψ_nl = π/3
    r_el = ltf(latlon, 0) ∘ Rz(ψ_nl)
    @test get_ψ_nl(r_el) ≈ ψ_nl

    #conversion torture test: from LatLon to NVector, to r_en, then r_el, then
    #NVector, back to LatLon
    ϕ_range = range(-π/2, π/2, length = 10)
    λ_range = range(-π, π, length = 10)
    ψ_nl_range = range(-π, π, length = 10)

    ll_array = [LatLon(ϕ, λ) for (ϕ, λ, _) in Iterators.product(ϕ_range, λ_range, ψ_nl_range)]
    ψ_array = [ψ_nl for (_, _, ψ_nl) in Iterators.product(ϕ_range, λ_range, ψ_nl_range)]
    ll_array_test = [ll |> NVector |> ltf |> (r -> r ∘ Rz(ψ_nl)) |> NVector |> LatLon for (ll, ψ_nl) in zip(ll_array, ψ_array)]

    @test all(ll_array .≈ ll_array_test)

end

function test_Altitude()

    h_ellip = Altitude{Ellipsoidal}(1500)
    h_orth = Altitude{Orthometric}(1500)
    Δh = 300
    loc = NVector([3,1,-5])

    #cannot convert from one datum to another without a 2D location
    @test_throws MethodError Altitude{Ellipsoidal}(h_orth)
    @test_throws MethodError Altitude{Orthometric}(h_ellip)

    #however, with a 2D location we can convert back and forth
    @test Altitude{Orthometric}(Altitude{Ellipsoidal}(h_orth, loc), loc) == h_orth
    @test Altitude{Orthometric}(h_ellip, loc) isa Altitude{Orthometric}
    @test Altitude{Ellipsoidal}(h_orth, loc) isa Altitude{Ellipsoidal}
    @test Altitude{Geopotential}(h_ellip, loc) isa Altitude{Geopotential}
    @test Altitude{Geopotential}(h_orth) isa Altitude{Geopotential}

    #operations and conversions
    @test (h_ellip + Δh) isa HEllip
    @test (h_orth + Δh) isa HOrth
    @test Float64(h_ellip + Δh) == Float64(h_ellip) + Δh
    @test Float64(h_orth + Δh) == Float64(h_orth) + Δh

    @test h_ellip - h_ellip/2 isa Float64
    @test h_ellip - h_ellip/2 == h_ellip/2
    @test h_ellip + Float64(h_ellip/2) isa HEllip
    @test 0.7h_ellip isa HEllip
    @test h_ellip/2 isa HEllip

    @test h_ellip + Float64(h_ellip/2) == 3h_ellip/2
    @test h_ellip + Float64(h_ellip/2) ≈ 3h_ellip/2
    @test h_ellip >= h_ellip
    @test h_ellip <= h_ellip
    @test !(h_ellip < h_ellip)
    @test !(h_ellip > h_ellip)

    @test h_orth == 1500
    @test h_orth ≈ 1500
    @test h_orth >= 1500
    @test h_orth <= 1500
    @test !(h_orth < 1500)
    @test !(h_orth > 1500)

end

function test_GeographicLocation()

    p_nve = Geographic(h = HOrth(1500))
    p_llo = Geographic{LatLon,Orthometric}(p_nve)
    @test p_nve isa Geographic{NVector,Orthometric}
    @test Geographic(loc = LatLon(), h = HOrth()) isa Geographic{LatLon,Orthometric}

    #conversion
    @test Geographic{NVector,Ellipsoidal}(p_llo) isa Geographic{NVector,Ellipsoidal}

    @test_throws ArgumentError p_llo == p_llo
    @test p_nve == p_nve #strict equality only defined for Geographic{NVector}
    @test p_llo ≈ p_llo
    @test p_llo ≈ p_nve
    @test -(-p_llo) ≈ p_llo #requires Cartesian
    @test -(-p_nve) ≈ p_nve #requires Cartesian

    @test ltf(p_nve, π/3) == ltf(p_nve.loc, π/3)
    @test radii(p_nve) == radii(p_nve.loc)
    @test gravity(p_llo) ≈ gravity(p_nve)
    @test g_n(p_llo) ≈ g_n(p_nve)
    @test G_n(p_llo) ≈ G_n(p_nve)

end

function test_CartesianLocation()

    #construction from Geographic
    p_nvo = Geographic(NVector([3,1,-3]), HOrth(10000))
    p_lle = Geographic{LatLon,Ellipsoidal}(p_nvo)
    r = Cartesian(p_nvo)
    @test r ≈ Cartesian(p_lle)

    @test r == r
    @test -(-r) == r
    @test !(-r == r)
    @test r ≈ r

    #conversion torture test
    ftest(p) = p |> Cartesian |> Geographic{NVector,Ellipsoidal} |> Cartesian |>
                Geographic{NVector,Orthometric} |> Cartesian |>
                Geographic{LatLon, Ellipsoidal} |> Geographic{LatLon, Geopotential} |>
                Geographic{LatLon, Ellipsoidal}

    ϕ_range = range(-π/2, π/2, length = 10)
    λ_range = range(-π, π, length = 10)
    h_range = range(Geodesy.h_min + 200, 10000, length = 10)

    p_array = [Geographic(LatLon(ϕ, λ), HOrth(h)) for (ϕ, λ, h) in Iterators.product(ϕ_range, λ_range, h_range)]
    p_array_test = [ftest(p) for p in p_array]
    @test all(p_array .≈ p_array_test)

end

end #module