# The stepper seam (§10.2): loop ownership stops at one operation — *advance
# the continuous state from `t` by `h`* — delegated across this narrow
# interface so the integration method can be replaced without the loop
# changing. The contract has three clauses, each answered by dispatch on the
# stepper:
#
#   - **advance by arbitrary `h`** — `step!(m, sim, h)`: the loop lands on tick
#     boundaries and resumes from localized event times;
#   - **dense output on demand over the last completed step** — `dense!`, built
#     lazily on the pair `startpoint` retains, because only event localization
#     (§10.4) ever asks;
#   - **one-step methods only** — handlers reset state discontinuously, and a
#     one-step method restarts from a new state for free (D-017).
#
# The seam is never entered empty (§10.2): the framework short-circuits an
# empty state on its own side — `step!(sim, h)` in sim.jl — so no backend ever
# faces N = 0. Both first-cut backends are fixed-step, zero-allocation and
# generic in the scalar; the method is a deployment binding (`method = RK4`,
# the default), and nothing outside this file knows which one ran.

abstract type AbstractStepper end

"""
    startpoint(m)

The seam's retained pair `(xₙ, ẋₙ)`: the state and derivative at the start of
the last completed step, surviving in the backend's own scratch. This is the
data clause of the dense-output obligation — the frame loop's θ = 0 validation
reads the state half directly, and the generic `dense!` builds on both halves.
"""
function startpoint end

"""
    RK4(T, n)

The classical fourth-order Runge–Kutta method over a flat buffer of `n`
scalars of type `T` — the default backend (§10.2). Owns exactly its own
scratch: the segment-start state and the four stage derivatives.
"""
struct RK4{T} <: AbstractStepper
    x₀::Vector{T}
    k₁::Vector{T}
    k₂::Vector{T}
    k₃::Vector{T}
    k₄::Vector{T}
end
RK4(::Type{T}, n::Int) where {T} = RK4{T}(ntuple(_ -> zeros(T, n), 5)...)

function step!(m::RK4, sim, h)
    x, ẋ = sim.xbuf, sim.ẋbuf
    (; x₀, k₁, k₂, k₃, k₄) = m
    t₀ = sim.clock.t
    copyto!(x₀, x)

    evaluate!(sim); copyto!(k₁, ẋ)
    _advance!(x, x₀, k₁, h / 2); sim.clock.t = t₀ + h / 2
    evaluate!(sim); copyto!(k₂, ẋ)
    _advance!(x, x₀, k₂, h / 2)
    evaluate!(sim); copyto!(k₃, ẋ)
    _advance!(x, x₀, k₃, h); sim.clock.t = t₀ + h
    evaluate!(sim); copyto!(k₄, ẋ)

    @inbounds for i in eachindex(x)
        x[i] = x₀[i] + (h / 6) * (k₁[i] + 2k₂[i] + 2k₃[i] + k₄[i])
    end
    nothing
end

startpoint(m::RK4) = (m.x₀, m.k₁)

"""
    Heun(T, n)

Heun's second-order method (explicit trapezoidal) over a flat buffer of `n`
scalars of type `T` — the other first-cut backend (§10.2). Owns exactly its
own scratch: the segment-start state and the two stage derivatives.
"""
struct Heun{T} <: AbstractStepper
    x₀::Vector{T}
    k₁::Vector{T}
    k₂::Vector{T}
end
Heun(::Type{T}, n::Int) where {T} = Heun{T}(ntuple(_ -> zeros(T, n), 3)...)

function step!(m::Heun, sim, h)
    x, ẋ = sim.xbuf, sim.ẋbuf
    (; x₀, k₁, k₂) = m
    t₀ = sim.clock.t
    copyto!(x₀, x)

    evaluate!(sim); copyto!(k₁, ẋ)
    _advance!(x, x₀, k₁, h); sim.clock.t = t₀ + h
    evaluate!(sim); copyto!(k₂, ẋ)

    @inbounds for i in eachindex(x)
        x[i] = x₀[i] + (h / 2) * (k₁[i] + k₂[i])
    end
    nothing
end

startpoint(m::Heun) = (m.x₀, m.k₁)

@inline function _advance!(x, x₀, k, dt)
    @inbounds for i in eachindex(x)
        x[i] = x₀[i] + dt * k[i]
    end
end

"""
    dense!(m, x̂, x₁, ẋ₁, θ, h′)

Dense output over the last completed step (§10.2): write x̂(θ), θ ∈ [0, 1],
into `x̂`, from the seam's retained `startpoint` pair and the arrival pair
`(x₁, ẋ₁)` the caller holds — the caller pays ẋₙ₊₁ lazily, which is why it is
an argument rather than a retention. The generic method is the cubic Hermite
over those four, derivatives scaled by the segment width: uniform accuracy
O(h⁴), the standard pairing for a one-step method of order ≤ 4 — one order
below RK4's discrete solution, above Heun's — and why nothing more expensive
is worth running trials against. Both first-cut backends inherit it; a
higher-order backend would override with its own formula.
"""
function dense!(m::AbstractStepper, x̂, x₁, ẋ₁, θ::Float64, h′)
    (x₀, ẋ₀) = startpoint(m)
    θ² = θ * θ; θ³ = θ² * θ
    b₀ = 2θ³ - 3θ² + 1
    b₁ = 3θ² - 2θ³
    d₀ = (θ³ - 2θ² + θ) * h′
    d₁ = (θ³ - θ²) * h′
    @inbounds for i in eachindex(x̂)
        x̂[i] = b₀ * x₀[i] + d₀ * ẋ₀[i] + b₁ * x₁[i] + d₁ * ẋ₁[i]
    end
    nothing
end
