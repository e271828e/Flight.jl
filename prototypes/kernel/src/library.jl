# The coverage-driven component set for increment 2. Three types, chosen so
# that between them they exercise every continuous-tier shape the executor has
# to handle — and no more.

using StaticArrays

"""
Damped second-order plant. Carries state, publishes a **stage-1** port (`y`,
state-derived, no feedthrough — the port that lets a feedback loop close
legally) and a **stage-2** port (`power`, input-dependent), and defines `f`.
"""
struct Plant
    ω::Float64
    ζ::Float64
    q₀::SVector{2,Float64}
end

Plant(; ω = 2.0, ζ = 0.1, q₀ = SVector(0.0, 0.0)) = Plant(ω, ζ, q₀)

init_x(c::Plant) = (q = c.q₀,)
input_types(::Plant, ::Type{T}) where {T <: Real} = (u = T,)
output_types(::Plant, ::Type{T}) where {T <: Real} = (y = T, power = T)

h_x(::Plant, (; x)) = (y = x.q[1],)
h_xu(::Plant, (; x, u)) = (power = u.u * x.q[2],)

function f(c::Plant, (; x, u))
    q, ω, ζ = x.q, c.ω, c.ζ
    (q = SVector(q[2], -ω^2 * q[1] - 2ζ * ω * q[2] + u.u),)
end

"""
Proportional gain: **stateless**, stage 2 only. The three-level funnel of §5.2
in its smallest instance — a component that legitimately writes `h_xu` while
owning no state at all, so its bundle carries `u` and `t` and nothing else.
"""
struct Gain
    k::Float64
end

input_types(::Gain, ::Type{T}) where {T <: Real} = (e = T,)
output_types(::Gain, ::Type{T}) where {T <: Real} = (out = T,)

h_xu(c::Gain, (; u)) = (out = c.k * u.e,)

"""
Two-input summing junction (§6.2): aggregation is an explicit, ordered entry in
the schedule, never an implicit fan-in.
"""
struct Sum
    sa::Float64
    sb::Float64
end

Sum(; sa = 1.0, sb = -1.0) = Sum(sa, sb)

input_types(::Sum, ::Type{T}) where {T <: Real} = (a = T, b = T)
output_types(::Sum, ::Type{T}) where {T <: Real} = (e = T,)

h_xu(c::Sum, (; u)) = (e = c.sa * u.a + c.sb * u.b,)

# --- the discrete tier --------------------------------------------------------
# The plain declaration arities are the tier (D-166/D-167): these components
# declare the pinned world, and nothing about them walks with the activation.

"""
Discrete integrator: publishes its state from **stage 1** — the loop-breaking
port, exactly as on the continuous tier — and accumulates in `g`, which is
where `Δt` earns its place in the bundle.
"""
struct DiscreteIntegrator
    k::Float64
end

init_x(::DiscreteIntegrator) = (s = 0.0,)
input_types(::DiscreteIntegrator) = (e = Float64,)
output_types(::DiscreteIntegrator) = (u = Float64,)

h_x(::DiscreteIntegrator, (; x)) = (u = x.s,)
g(c::DiscreteIntegrator, (; x, u, Δt)) = (s = x.s + c.k * Δt * u.e,)

"""
Tick counter: `Int` state, `Int` and `Bool` ports. The second and third store
buffers of the bundle, which a continuous model never needs (D-162).
"""
struct TickCounter end

init_x(::TickCounter) = (n = 0,)
output_types(::TickCounter) = (n = Int, even = Bool)

h_x(::TickCounter, (; x)) = (n = x.n, even = iseven(x.n))
g(::TickCounter, (; x)) = (n = x.n + 1,)

"""
Two-channel exponential smoother, written in §7.3's blessed idiom: the in-place
math runs on the workspace, and what reaches the store is an isbits snapshot.
Nothing carries between calls — the scratch is garbage until written.
"""
struct Smoother
    α::Float64
end

init_x(::Smoother) = (v = SVector(0.0, 0.0),)
input_types(::Smoother) = (a = Float64, b = Float64)
output_types(::Smoother) = (v = SVector{2,Float64},)
workspace(::Smoother) = (tmp = Vector{Float64}(undef, 2),)

h_x(::Smoother, (; x)) = (v = x.v,)

function g(c::Smoother, (; x, u, ws))
    ws.tmp[1] = u.a
    ws.tmp[2] = u.b
    for i in 1:2
        ws.tmp[i] = c.α * ws.tmp[i] + (1 - c.α) * x.v[i]
    end
    (v = SVector{2,Float64}(ws.tmp[1], ws.tmp[2]),)
end

# --- the other two declarations the bundle law owes -------------------------

"""
A continuous workspace user: the same contract at the other arity, allocating
at the activation scalar so the scratch follows `Dual` when the activation does.
"""
struct WorkGain
    k::Float64
end

input_types(::WorkGain, ::Type{T}) where {T <: Real} = (in = T,)
output_types(::WorkGain, ::Type{T}) where {T <: Real} = (out = T,)
workspace(::WorkGain, ::Type{T}) where {T <: Real} = (tmp = Vector{T}(undef, 1),)

function h_xu(c::WorkGain, (; u, ws))
    ws.tmp[1] = c.k * u.in
    (out = ws.tmp[1],)
end

"""
A mode-carrying source: modes are continuous-only, and in this increment they
are read but never written — handlers are the event system, which is absent.
The branch returning a literal is the constant-branch idiom (D-166).
"""
struct ModedSource end

init_m(::ModedSource) = (phase = :idle,)
output_types(::ModedSource, ::Type{T}) where {T <: Real} = (out = T,)

h_x(::ModedSource, (; m)) = (out = m.phase === :idle ? 0.0 : 1.0,)

# --- the reference model ------------------------------------------------------

"""
    feedback_model(; k, feedback_port = :y)

Closed loop: `sum` differences a root-slot reference against the plant's
feedback port, `ctl` scales the error, the plant integrates it.

The default `feedback_port = :y` closes the loop through the plant's *stage-1*
port, which carries no input dependence — so the loop is legal and the schedule
is `sum → ctl → plant`. Passing `:power` closes it through a stage-2 port
instead, which is a genuine algebraic loop and must be rejected at build time
(§5.5).
"""
function feedback_model(; k = 4.0, ω = 2.0, ζ = 0.1, q₀ = SVector(0.0, 0.0),
                        feedback_port::Symbol = :y)
    ModelSpec([
        ComponentSpec(:plant, Plant(; ω, ζ, q₀), [:u => (:ctl, :out)]),
        ComponentSpec(:ctl, Gain(k), [:e => (:sum, :e)]),
        # `sum.a` is deliberately unwired: it is a root slot, its initial value
        # synthesized by `probe_value` and written by `set_slot!` (§9.3, §11.3).
        ComponentSpec(:sum, Sum(), [:b => (:plant, feedback_port)]),
    ])
end

"""
    sampled_loop(; k, kI, ω, ζ)

The sampled-data closed loop: a continuous plant under a discrete integrator's
zero-order-held command. `sum` differences the root-slot reference against the
plant's stage-1 port, the discrete `ctl` accumulates that error at its own tick
and holds `u` between ticks, and the plant integrates it.

The hold is not implemented anywhere: `ctl`'s entries are absent from the
interior sweep, so its cell simply cannot change between boundaries (§10.5).
"""
function sampled_loop(; kI = 3.0, ω = 2.0, ζ = 0.1)
    ModelSpec([
        ComponentSpec(:plant, Plant(; ω, ζ), [:u => (:ctl, :u)]),
        ComponentSpec(:ctl, DiscreteIntegrator(kI), [:e => (:sum, :e)]),
        ComponentSpec(:sum, Sum(), [:b => (:plant, :y)]),
    ])
end
