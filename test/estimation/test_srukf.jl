module TestSRUKF
using ComponentArrays
using StaticArrays
using LinearAlgebra
using Test

using Flight.FlightPhysics

const δx_template = ComponentVector(ρ = zeros(3), δv = zeros(3));
const x_template = ComponentVector(q = zeros(4), v = zeros(3));

const δx_axes = getaxes(δx_template)[1]
const x_axes = getaxes(x_template)[1]


function g_x_plus!(x1::AbstractVector{<:AbstractFloat}, x0::AbstractVector{<:AbstractFloat}, δx::AbstractVector{<:AbstractFloat})
    x0_ca = ComponentVector(x0, x_axes)
    x1_ca = ComponentVector(x1, x_axes)
    δx_ca = ComponentVector(δx, δx_axes)

    q0 = RQuat(x0_ca.q)
    v0 = SVector{3,Float64}(x0_ca.v)

    ρ = RVec(δx_ca.ρ)
    δv = SVector{3,Float64}(δx_ca.δv)

    q1 = q0 ∘ ρ
    v1 = v0 + δv

    #mutates x1
    x1_ca.q = q1[:]
    x1_ca.v = v1
end

function g_x_minus!(δx::AbstractVector{<:AbstractFloat}, x1::AbstractVector{<:AbstractFloat}, x0::AbstractVector{<:AbstractFloat})
    x0_ca = ComponentVector(x0, x_axes)
    x1_ca = ComponentVector(x1, x_axes)
    δx_ca = ComponentVector(δx, δx_axes)

    q0 = RQuat(x0_ca.q)
    v0 = SVector{3,Float64}(x0_ca.v)

    q1 = RQuat(x1_ca.q)
    v1 = SVector{3,Float64}(x1_ca.v)

    #here, q1 = q0 ∘ q01, so q01 = q0' ∘ q1
    ρ = RVec(q0' ∘ q1)
    δv = v1 - v0

    #mutates δx
    δx_ca.ρ = ρ[:]
    δx_ca.δv = δv
end


function test00()

    # A = rand(6, 6)
    A = 2e-2*diagm(ones(6))
    P_δx_0 = A * A'
    P_δx_0 = ComponentMatrix(P_δx_0, δx_axes, δx_axes)
    S_δx_0 = cholesky(P_δx_0).L
    S_w_0 = 1e-2 * LowerTriangular(diagm(ones(6)))

    x_0 = copy(x_template)
    x_1 = copy(x_template)
    δx = copy(δx_template)

    x_0.q .= RQuat([1,2,-2.5,3])[:]
    x_0.v .= [4, 5, -3]

    δx.ρ = rand(3)
    δx.δv = rand(3)
    δx_copy = copy(δx)

    g_x_plus!(x_1, x_0, δx)
    g_x_minus!(δx, x_1, x_0)

    #first test the generalized addition and subtraction
    @test all(isapprox.(δx, δx_copy))

    #from here on, we are inside the SRUT
    N_δx = size(S_δx_0)[1]
    N_w = size(S_w_0)[1]
    α = 1e-3
    β = 2
    κ = 0
    S_δa = [S_δx_0 zeros(N_δx, N_w); zeros(N_w, N_δx) S_w_0]

    #defaults:
    # g_u_plus! = (u1, u0, δu) -> (u1 .= u0 .+ δu)
    # g_u_minus! = (δx, u1, u0) -> (δx .= u1 .- u0)

end

end