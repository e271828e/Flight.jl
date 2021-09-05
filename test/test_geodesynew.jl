using Test
using LinearAlgebra

# using Flight.Attitude
# using Flight.Quaternions


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

#
############## LatLon

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


############### START TESTS FOR A3DLocation

h_ellip = EllipsoidalAlt(1500)
h_orth = OrthometricAlt(1500)
Δh = 300

@test (h_ellip + Δh) isa EllipsoidalAlt
@test (h_orth + Δh) isa OrthometricAlt
@test Float64(h_ellip + Δh) == Float64(h_ellip) + Δh
@test Float64(h_orth + Δh) == Float64(h_orth) + Δh