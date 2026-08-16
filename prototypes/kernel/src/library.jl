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
