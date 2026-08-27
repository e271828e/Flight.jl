# The coverage-driven component set for increment 2. Three types, chosen so
# that between them they exercise every continuous-tier shape the executor has
# to handle — and no more.

using StaticArrays

"""
Damped second-order plant. Carries state, publishes a **stage-1** port (`y`,
state-derived, no feedthrough — the port that lets a feedback loop close
legally) and a **stage-2** port (`power`, input-dependent), and defines `f`.
"""
struct Plant <: AbstractComponent
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
struct Gain <: AbstractComponent
    k::Float64
end

input_types(::Gain, ::Type{T}) where {T <: Real} = (e = T,)
output_types(::Gain, ::Type{T}) where {T <: Real} = (out = T,)

h_xu(c::Gain, (; u)) = (out = c.k * u.e,)

"""
Two-input summing junction (§6.2): aggregation is an explicit, ordered entry in
the schedule, never an implicit fan-in.
"""
struct Sum <: AbstractComponent
    sa::Float64
    sb::Float64
end

Sum(; sa = 1.0, sb = -1.0) = Sum(sa, sb)

input_types(::Sum, ::Type{T}) where {T <: Real} = (a = T, b = T)
output_types(::Sum, ::Type{T}) where {T <: Real} = (e = T,)

h_xu(c::Sum, (; u)) = (e = c.sa * u.a + c.sb * u.b,)

# --- the discrete tier --------------------------------------------------------
# The tier's own name family (D-195) — `init_s`, `h_s`/`h_su`, `g` — beside the
# plain declaration arities (D-166/D-167): these components declare the pinned
# world, and nothing about them walks with the activation.

"""
Discrete integrator: publishes its state from **stage 1** — the loop-breaking
port, exactly as on the continuous tier — and accumulates in `g`, which is
where `Δt` earns its place in the bundle.
"""
struct DiscreteIntegrator <: AbstractComponent
    k::Float64
end

init_s(::DiscreteIntegrator) = (acc = 0.0,)
input_types(::DiscreteIntegrator) = (e = Float64,)
output_types(::DiscreteIntegrator) = (u = Float64,)

h_s(::DiscreteIntegrator, (; s)) = (u = s.acc,)
g(c::DiscreteIntegrator, (; s, u, Δt)) = (acc = s.acc + c.k * Δt * u.e,)

"""
Tick counter: `Int` state, `Int` and `Bool` ports. The second and third store
buffers of the bundle, which a continuous model never needs (D-162).
"""
struct TickCounter <: AbstractComponent end

init_s(::TickCounter) = (n = 0,)
output_types(::TickCounter) = (n = Int, even = Bool)

h_s(::TickCounter, (; s)) = (n = s.n, even = iseven(s.n))
g(::TickCounter, (; s)) = (n = s.n + 1,)

"""
Two-channel exponential smoother, written in §7.3's blessed idiom: the in-place
math runs on the workspace, and what reaches the store is an isbits snapshot.
Nothing carries between calls — the scratch is garbage until written.
"""
struct Smoother <: AbstractComponent
    α::Float64
end

init_s(::Smoother) = (v = SVector(0.0, 0.0),)
input_types(::Smoother) = (a = Float64, b = Float64)
output_types(::Smoother) = (v = SVector{2,Float64},)
workspace(::Smoother) = (tmp = Vector{Float64}(undef, 2),)

h_s(::Smoother, (; s)) = (v = s.v,)

function g(c::Smoother, (; s, u, ws))
    ws.tmp[1] = u.a
    ws.tmp[2] = u.b
    for i in 1:2
        ws.tmp[i] = c.α * ws.tmp[i] + (1 - c.α) * s.v[i]
    end
    (v = SVector{2,Float64}(ws.tmp[1], ws.tmp[2]),)
end

# --- the other two declarations the bundle law owes -------------------------

"""
A continuous workspace user: the same contract at the other arity, allocating
at the activation scalar so the scratch follows `Dual` when the activation does.
"""
struct WorkGain <: AbstractComponent
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
A mode-carrying source: modes are continuous-only, read here and written by
nothing — a component may legitimately declare modes no event of its own
transitions (§8.2). The branch returning a literal is the constant-branch
idiom (D-166).
"""
struct ModedSource <: AbstractComponent end

init_m(::ModedSource) = (phase = :idle,)
output_types(::ModedSource, ::Type{T}) where {T <: Real} = (out = T,)

h_x(::ModedSource, (; m)) = (out = m.phase === :idle ? 0.0 : 1.0,)

# --- the event coverage set (§2.1, §10.6) -------------------------------------
# Guards and handlers are ordinary named functions referenced by `events` —
# nothing global-generic about them, which is why they carry component-prefixed
# names here rather than adding methods to a framework surface.

"""
Threshold trigger: mode-only state, a `Bool` guard on its input against a
level, and a *sticky* predicate — the input stays above the level after the
transition, so edge semantics must make it fire exactly once. `count` records
firings, which is what the tests read.
"""
struct Trigger <: AbstractComponent
    level::Float64
end

init_m(::Trigger) = (state = :armed, count = 0)
input_types(::Trigger, ::Type{T}) where {T <: Real} = (sig = T,)
output_types(::Trigger, ::Type{T}) where {T <: Real} = (on = Bool,)

h_x(::Trigger, (; m)) = (on = m.state === :fired,)

trigger_guard(c::Trigger, (; u)) = u.sig ≥ c.level
trigger_handler(::Trigger, (; m)) = (m = (state = :fired, count = m.count + 1),)
events(::Trigger) = (fire = Event(trigger_guard, trigger_handler),)

"""
Cascade follower: goes `:idle → :on` on the rising edge of its `Bool` input. A
chain of these under a `Trigger` is §10.6's logically-simultaneous cascade —
settled within one boundary, one iteration round per link, independent of `h`.
"""
struct Follower <: AbstractComponent end

init_m(::Follower) = (state = :idle,)
input_types(::Follower, ::Type{T}) where {T <: Real} = (go = Bool,)
output_types(::Follower, ::Type{T}) where {T <: Real} = (on = Bool,)

h_x(::Follower, (; m)) = (on = m.state === :on,)

follower_guard(::Follower, (; u)) = u.go
follower_handler(::Follower, (; m)) = (m = (state = :on,),)
events(::Follower) = (engage = Event(follower_guard, follower_handler),)

"""
Sawtooth: `q̇ = rate`, and a **sign-form** guard `q − 1` whose handler carries
the overshoot across the reset, `q ← q − 1`. The sign form declares the event
localized (§10.4, D-179) — and the carrying handler makes the *trajectory*
invariant to where the firing lands, `q(t) = rate·t − #wraps` either way, so
the boundary-detected recursion stays its exact reference. What localization
changes here is only the timing resolution `Bouncer` makes observable.
"""
struct Sawtooth <: AbstractComponent
    rate::Float64
end

init_x(::Sawtooth) = (q = 0.0,)
output_types(::Sawtooth, ::Type{T}) where {T <: Real} = (q = T,)

h_x(::Sawtooth, (; x)) = (q = x.q,)
f(c::Sawtooth, (; x)) = (q = c.rate,)

sawtooth_guard(::Sawtooth, (; x)) = x.q - 1.0
sawtooth_handler(::Sawtooth, (; x)) = (x = (q = x.q - 1.0,),)
events(::Sawtooth) = (wrap = Event(sawtooth_guard, sawtooth_handler),)

"""
Unit-circle rotor: `ċ = -ω s, ṡ = ω c` under a renormalizing `project` — the
state manifold in miniature. Boundary zero already projects (§14.5), so a
deliberately off-manifold `init_x` lands on the circle before the first step.
"""
struct Rotor <: AbstractComponent
    ω::Float64
    r₀::SVector{2,Float64}
end

Rotor(; ω = 1.0, r₀ = SVector(1.0, 0.0)) = Rotor(ω, r₀)

init_x(c::Rotor) = (r = c.r₀,)
output_types(::Rotor, ::Type{T}) where {T <: Real} = (c = T,)

h_x(::Rotor, (; x)) = (c = x.r[1],)
f(c::Rotor, (; x)) = (r = SVector(-c.ω * x.r[2], c.ω * x.r[1]),)
project(::Rotor, x) = (r = x.r / sqrt(x.r[1]^2 + x.r[2]^2),)

"""
Two events toggling one mode: each transition re-arms the other, so the pair
chatters without bound and each spends `firing_budget` at the first boundary —
the §10.6 degradation path, warned and bounded while the rest of the model
iterates untouched. `flips` counts every firing.
"""
struct Chatterer <: AbstractComponent end

init_m(::Chatterer) = (v = false, flips = 0)
output_types(::Chatterer, ::Type{T}) where {T <: Real} = (v = Bool,)

h_x(::Chatterer, (; m)) = (v = m.v,)

chatter_up(::Chatterer, (; m)) = !m.v
chatter_down(::Chatterer, (; m)) = m.v
chatter_flip(::Chatterer, (; m)) = (m = (v = !m.v, flips = m.flips + 1),)
events(::Chatterer) = (up = Event(chatter_up, chatter_flip),
                       down = Event(chatter_down, chatter_flip))

"""
Two events on one predicate: both edges rise in the same round, the first fires
by declaration order, and the second is eligible but blocked — its sample is
*not* overwritten (D-191), so the standing edge fires it in the next round.
"""
struct TwoShot <: AbstractComponent end

init_m(::TwoShot) = (a = false, b = false)
input_types(::TwoShot, ::Type{T}) where {T <: Real} = (sig = T,)
output_types(::TwoShot, ::Type{T}) where {T <: Real} = (a = Bool,)

h_x(::TwoShot, (; m)) = (a = m.a,)

twoshot_guard(::TwoShot, (; u)) = u.sig ≥ 1.0
twoshot_a(::TwoShot, (; m)) = (m = (; a = true),)
twoshot_b(::TwoShot, (; m)) = (m = (; b = true),)
events(::TwoShot) = (first = Event(twoshot_guard, twoshot_a),
                     second = Event(twoshot_guard, twoshot_b))

"""
The blocked-then-falsified variant: the second event's premise includes the
first's *not* having fired, so the round-1 block defers it into a round where
the premise is gone — it re-decides against the post-transition sweep and never
fires, where a within-round sequence would have fired it on the stale premise.
"""
struct Preempted <: AbstractComponent end

init_m(::Preempted) = (a = false, b = false)
input_types(::Preempted, ::Type{T}) where {T <: Real} = (sig = T,)
output_types(::Preempted, ::Type{T}) where {T <: Real} = (a = Bool,)

h_x(::Preempted, (; m)) = (a = m.a,)

preempted_guard_a(::Preempted, (; u)) = u.sig ≥ 1.0
preempted_guard_b(::Preempted, (; u, m)) = u.sig ≥ 1.0 && !m.a
preempted_a(::Preempted, (; m)) = (m = (; a = true),)
preempted_b(::Preempted, (; m)) = (m = (; b = true),)
events(::Preempted) = (first = Event(preempted_guard_a, preempted_a),
                       second = Event(preempted_guard_b, preempted_b))

# --- the localization coverage set (§10.4) --------------------------------------

"""
Crossing stamper: a **sign-form** guard on its input against a level, and a
handler that records the boundary time it fired at — the direct observable for
`t*`. A localized firing stamps within `localization_tol` of the true crossing
(exactly, when the trajectory feeding it is polynomial of degree ≤ 3, where the
Hermite interpolant is exact); an epoch-caused or degenerate edge stamps its
frame top exactly.
"""
struct Stamper <: AbstractComponent
    level::Float64
end

init_m(::Stamper) = (t_fired = -1.0, count = 0)
input_types(::Stamper, ::Type{T}) where {T <: Real} = (sig = T,)
output_types(::Stamper, ::Type{T}) where {T <: Real} = (armed = Bool,)

h_x(::Stamper, (; m)) = (armed = m.count == 0,)

stamper_guard(c::Stamper, (; u)) = u.sig - c.level
stamper_handler(::Stamper, (; m, t)) = (m = (t_fired = t, count = m.count + 1),)
events(::Stamper) = (cross = Event(stamper_guard, stamper_handler),)

"""
The gate idiom (§10.4): a mixed predicate in its blessed spelling,
`(gate) ? σ : -one(σ)` — the `Bool` factor rides the branch, the continuous
factor rides the value, and the guard's return type stays the nominal scalar,
so the event is localized. Trial evaluations vary only θ, with `u` fixed
through a localization, so the gate is constant over the bracket and σ
restricted to it is the continuous atom.
"""
struct GatedStamper <: AbstractComponent
    level::Float64
end

init_m(::GatedStamper) = (t_fired = -1.0, count = 0)
input_types(::GatedStamper, ::Type{T}) where {T <: Real} = (sig = T, gate = Bool)
output_types(::GatedStamper, ::Type{T}) where {T <: Real} = (armed = Bool,)

h_x(::GatedStamper, (; m)) = (armed = m.count == 0,)

function gated_stamper_guard(c::GatedStamper, (; u))
    σ = u.sig - c.level
    u.gate ? σ : -one(σ)
end
gated_stamper_handler(::GatedStamper, (; m, t)) = (m = (t_fired = t, count = m.count + 1),)
events(::GatedStamper) = (cross = Event(gated_stamper_guard, gated_stamper_handler),)

"""
Resetting ramp: `q̇ = rate` against a sign guard at `level`, and a handler that
*discards* the overshoot, `q ← 0` — unlike `Sawtooth`'s carrying handler, so
the trajectory itself depends on where the firing lands. Localized resets give
the exact period `level/rate`; boundary-resolution resets accumulate the
overshoot as phase error, which is the observable §10.4 buys.
"""
struct Bouncer <: AbstractComponent
    rate::Float64
    level::Float64
end

init_x(::Bouncer) = (q = 0.0,)
init_m(::Bouncer) = (count = 0,)
output_types(::Bouncer, ::Type{T}) where {T <: Real} = (q = T,)

h_x(::Bouncer, (; x)) = (q = x.q,)
f(c::Bouncer, (; x)) = (q = c.rate,)

bouncer_guard(c::Bouncer, (; x)) = x.q - c.level
bouncer_handler(::Bouncer, (; x, m)) = (x = (q = 0.0,), m = (count = m.count + 1,))
events(::Bouncer) = (reset = Event(bouncer_guard, bouncer_handler),)

"""
Relaxation chatterer: `q̇ = rate` against a sign guard at `level`, re-armed by
its own handler to `level − drop` — each remainder step re-crosses within the
same frame, so localizations pile up until `localization_budget` is spent and
the frame degrades to boundary granularity under a `ChatteringBudget` warning
(§10.4), while the run proceeds deterministically.
"""
struct Relaxer <: AbstractComponent
    rate::Float64
    level::Float64
    drop::Float64
end

init_x(::Relaxer) = (q = 0.0,)
init_m(::Relaxer) = (count = 0,)
output_types(::Relaxer, ::Type{T}) where {T <: Real} = (q = T,)

h_x(::Relaxer, (; x)) = (q = x.q,)
f(c::Relaxer, (; x)) = (q = c.rate,)

relaxer_guard(c::Relaxer, (; x)) = x.q - c.level
relaxer_handler(c::Relaxer, (; x, m)) =
    (x = (q = c.level - c.drop,), m = (count = m.count + 1,))
events(::Relaxer) = (pop = Event(relaxer_guard, relaxer_handler),)

# --- the termination coverage set (§13.5, §13.6) -------------------------------

"""
Overload monitor: §13.5's touchdown archetype — a **sign-form** guard on its
input against a level, a handler that latches `m.tripped`, and the sticky
`Bool` output face a `stop_on` policy names. The sign form declares the event
localized, so a run stopped on `tripped` ends at the crossing's `t*` boundary
with the crossing state as the terminal snapshot.
"""
struct Overload <: AbstractComponent
    level::Float64
end

init_m(::Overload) = (tripped = false,)
input_types(::Overload, ::Type{T}) where {T <: Real} = (sig = T,)
output_types(::Overload, ::Type{T}) where {T <: Real} = (tripped = Bool,)

h_x(::Overload, (; m)) = (tripped = m.tripped,)

overload_guard(c::Overload, (; u)) = u.sig - c.level
overload_handler(::Overload, (; m)) = (m = (tripped = true,),)
events(::Overload) = (trip = Event(overload_guard, overload_handler),)

"""
Exploder: the §13.6 specimen — `q̇ = 1` until its `arm` input goes true, then
its RHS throws `Exploded`. The throw escapes mid-integration, so the failing
frame has published nothing: what the abnormal tail leaves as final is the
last completed boundary's snapshot, which is the §13.6 discard-and-promote
made observable.
"""
struct Exploder <: AbstractComponent end

"Exploder's own exception type, so the tests assert the retained cause's identity."
struct Exploded <: Exception end

init_x(::Exploder) = (q = 0.0,)
input_types(::Exploder, ::Type{T}) where {T <: Real} = (arm = Bool,)
output_types(::Exploder, ::Type{T}) where {T <: Real} = (q = T,)

h_x(::Exploder, (; x)) = (q = x.q,)
f(::Exploder, (; x, u)) = u.arm ? throw(Exploded()) : (q = one(x.q),)

# --- the anonymous assembly (§8.5) --------------------------------------------

"""
`Group`: the on-the-fly assembly, one library type whose *values* are the ad-hoc
topologies. It needs no new rule — the container-children rule makes the
`children` field's elements children path-named `"children/key"`, and the three
declarations are ordinary functions of the instance, free to read its fields.

The type parameters carry the children's concrete types, so specialization is
unchanged; what is given up against a named type is dispatch, which exploratory
composition does not want.
"""
struct Group{C <: NamedTuple, W, I, O, R <: NamedTuple} <: AbstractComponent
    children::C      # component-typed elements → children by the container rule
    wires::W         # inert parameter data
    inputs::I
    outputs::O
    rates::R         # the ad-hoc rate scope; container elements are keyed `var"children/key"` (§8.7)
end

Group(children, wires, inputs, outputs) = Group(children, wires, inputs, outputs, (;))

child_connections(g::Group) = g.wires
input_connections(g::Group) = g.inputs
output_connections(g::Group) = g.outputs
sample_times(g::Group) = g.rates

# --- the reference models -----------------------------------------------------

"""
    feedback_model(; k, feedback_port = "y")

Closed loop as a `Group`: `sum` differences the `"ref"` input face against the
plant's feedback port, `ctl` scales the error, the plant integrates it. `"y"` is
its one output face.

The default `feedback_port = "y"` closes the loop through the plant's *stage-1*
port, which carries no input dependence — so the loop is legal and the schedule
is `sum → ctl → plant`. Passing `"power"` closes it through a stage-2 port
instead, which is a genuine algebraic loop and must be rejected at build time
(§5.5).
"""
function feedback_model(; k = 4.0, ω = 2.0, ζ = 0.1, q₀ = SVector(0.0, 0.0),
                        feedback_port::String = "y")
    Group((plant = Plant(; ω, ζ, q₀), ctl = Gain(k), sum = Sum()),
          ("children/ctl/out" => "children/plant/u",
           "children/sum/e" => "children/ctl/e",
           "children/plant/$feedback_port" => "children/sum/b"),
          # `sum.a` is claimed by no wire: the obligation is handed up to this
          # face, and at the root a face is a root input — its cell seeded by
          # `probe_value` for the build's own probes, and its initial value
          # authored by the init service's condition (§6.1, §11.3, §14.6).
          ("ref" => "children/sum/a",),
          ("children/plant/y" => "y",))
end

"""
    sampled_loop(; kI, ω, ζ)

The sampled-data closed loop: a continuous plant under a discrete integrator's
zero-order-held command. `sum` differences the `"ref"` face against the plant's
stage-1 port, the discrete `ctl` accumulates that error at its own tick and holds
`u` between ticks, and the plant integrates it.

Its two output faces are the mixed-tier pair: `"y"` sourced from the continuous
plant, `"cmd"` from the discrete controller (§8.6 — an assembly is tier-neutral,
and a face's type and tier are its internal endpoint's).

The hold is not implemented anywhere: `ctl`'s entries are absent from the
interior sweep, so its cell simply cannot change between boundaries (§10.5).
"""
function sampled_loop(; kI = 3.0, ω = 2.0, ζ = 0.1)
    Group((plant = Plant(; ω, ζ), ctl = DiscreteIntegrator(kI), sum = Sum()),
          ("children/ctl/u" => "children/plant/u",
           "children/sum/e" => "children/ctl/e",
           "children/plant/y" => "children/sum/b"),
          ("ref" => "children/sum/a",),
          ("children/plant/y" => "y", "children/ctl/u" => "cmd"))
end

# --- the named two-level assembly ---------------------------------------------
# Class by declaration shape (§8.5): `child_connections` and nothing else, on a
# plain struct whose component-typed fields are its children.

"""
The sampled loop as a *named* assembly, holding its three children in concretely
declared fields. An ancestor reaches none of them: every signal leaving this
level does so through a face declared here (§6.1, D-207), which is why the
plant's `power` port is re-exported below beside `y` and `cmd`.

`ctl_rate` is §10.5's exposed-multiplier idiom: a deployment preference surfaces
as a constructor parameter, the declaration stays the assembly's own
`sample_times`, and the component type stays rate-agnostic — `ctl` reads its
period from its bundle's `Δt` and nowhere else.
"""
struct SampledLoop <: AbstractComponent
    plant::Plant
    ctl::DiscreteIntegrator
    sum::Sum
    ctl_rate::Relative           # inert parameter data, read by `sample_times`
end

SampledLoop(; kI = 3.0, ω = 2.0, ζ = 0.1, ctl_rate = Relative(1)) =
    SampledLoop(Plant(; ω, ζ), DiscreteIntegrator(kI), Sum(), ctl_rate)

child_connections(::SampledLoop) =
    ("ctl/u" => "plant/u", "sum/e" => "ctl/e", "plant/y" => "sum/b")
input_connections(::SampledLoop) = ("ref" => "sum/a",)
output_connections(::SampledLoop) =
    ("plant/y" => "y", "ctl/u" => "cmd", "plant/power" => "power")
sample_times(l::SampledLoop) = (ctl = l.ctl_rate,)

"""
    Vehicle(; k, kI, ω, ζ)

Two levels: the vehicle scales its own `"ref"` face and hands the product to the
loop's input face, and re-exports three of the loop's signals — each through the
loop's own face, one level at a time, the plant's `power` having acquired a face
on the loop's boundary to be re-exported through (§6.1, D-207).

Its faces are the mixed-tier set: `"y"` derives from the continuous plant, `"cmd"`
from the discrete controller, and neither is declared anywhere — a face's type and
tier are its ultimate internal endpoint's (§8.6).
"""
struct Vehicle <: AbstractComponent
    loop::SampledLoop
    trim::Gain
end

Vehicle(; k = 1.0, kI = 3.0, ω = 2.0, ζ = 0.1) = Vehicle(SampledLoop(; kI, ω, ζ), Gain(k))

child_connections(::Vehicle) = ("trim/out" => "loop/ref",)
input_connections(::Vehicle) = ("ref" => "trim/e",)
output_connections(::Vehicle) =
    ("loop/y" => "y", "loop/cmd" => "cmd", "loop/power" => "power")

# --- the fragment-function idiom (§14.2) ----------------------------------------
# User-idiom material, not framework API: `condition` is an ordinary function
# dispatched on the component, shipped beside it, and nothing in `src/` outside
# this file knows the name exists. What it buys is locality — the caller says
# "start the plant at this displacement", and only `Plant` knows displacement and
# rate pack into `q` — composed by *pull* from the structure's owner, never by a
# schema routing sub-specs down the tree (D-064).

"""The plant's own vocabulary: displacement and rate, which it packs into `q`."""
condition(::Plant; y = 0.0, v = 0.0) = fragment(x = (q = SVector(y, v),))

"""The integrator's: the held command, which it accumulates in `acc`."""
condition(::DiscreteIntegrator; cmd = 0.0) = fragment(s = (acc = cmd,))

"""
Composition by pull (§14.2): the owner of the structure names its children and
scopes their fragments with `at`, and `combine` collects the siblings. It
authors no root input — `ref` is this assembly's *input face*, which is a root
input only when `SampledLoop` is itself the root, so the level that knows is
the one that owns the boundary.
"""
condition(l::SampledLoop; y = 0.0, v = 0.0, cmd = 0.0) =
    combine(at("plant", condition(l.plant; y, v)), at("ctl", condition(l.ctl; cmd)))

"""
The second level, and the one that owns the root boundary: it pulls the loop's
fragment under `at("loop", …)` — deep paths are compiled derivatives of this
nesting, never written by hand — and authors the root input its own contract
declares.
"""
condition(veh::Vehicle; ref = 0.0, kw...) =
    combine(at("loop", condition(veh.loop; kw...)), fragment(inputs = (ref = ref,)))

# --- the multi-rate coverage set (§10.5) ----------------------------------------

"""
Zero-order hold: a stateless discrete pass-through. Its one cell takes a fresh
sample at its own ticks and holds in between, so the tick pattern and the
deterministic aging of a stagger are directly observable through cell reads.
"""
struct ZOH <: AbstractComponent end

input_types(::ZOH) = (in = Float64,)
output_types(::ZOH) = (out = Float64,)
h_su(::ZOH, (; u)) = (out = u.in,)

"""Affine clock publisher: continuous, stateless, stage 1 — `out = c₀ + t`."""
struct Ramp <: AbstractComponent
    c₀::Float64
end

output_types(::Ramp, ::Type{T}) where {T <: Real} = (out = T,)
h_x(c::Ramp, (; t)) = (out = c.c₀ + t,)

"""
The rate scope of the spec's worked example (§9.2, §10.5): `inner` at the scope
base, `outer` staggered at `Relative(5, 2)`. `outer` samples a signal handed in
through the `g` face, so what it reads between that producer's ticks is the ZOH
aging the example computes.
"""
struct FCS <: AbstractComponent
    inner::ZOH
    outer::ZOH
end

child_connections(::FCS) = ()
input_connections(::FCS) = ("in" => "inner/in", "g" => "outer/in")
output_connections(::FCS) = ("inner/out" => "y_inner", "outer/out" => "y_outer")
sample_times(::FCS) = (inner = Relative(1), outer = Relative(5, 2))

"""
    MultiRate()

The §9.2 worked example: three discrete components under two scopes, `fcs` on
the root grid and `gnss` anchored at `Hz(50)`. Deployed at `Δt_base = 2 ms` the
compiled pairs are inner `(1, 0)`, outer `(5, 2)`, gnss `(10, 0)`, and one
hyperperiod is `lcm(Dᵢ) = 10` base ticks. The ramp source makes every sample
carry its own acquisition time, so the chart's dots — and the stagger's aging —
are readable off the cells.
"""
struct MultiRate <: AbstractComponent
    src::Ramp
    fcs::FCS
    gnss::ZOH
end

MultiRate(; c₀ = 1.0) = MultiRate(Ramp(c₀), FCS(ZOH(), ZOH()), ZOH())

child_connections(::MultiRate) =
    ("src/out" => "fcs/in", "src/out" => "gnss/in", "gnss/out" => "fcs/g")
output_connections(::MultiRate) =
    ("fcs/y_inner" => "inner", "fcs/y_outer" => "outer", "gnss/out" => "gnss")
sample_times(::MultiRate) = (fcs = Relative(1), gnss = Absolute(Hz(50)))

# --- the periphery's coverage set: devices and bindings (§11.3, §11.6) ----------

"""
    Pad(name)

A stub device: no hardware — the identity carrier for the roster. Mutable
deliberately, because identity is the instance (`===`, §11.3): an immutable
stub with equal fields would be egal to its twin, and two same-named `Pad`s
are meant to be two devices. Its loop body returns at once — §12.4(6)'s
voluntary exit, honest for a stub with nothing to wait on — so a rostered
`Pad` in a run is an identity holding a claim, its task already departed and
its staging cell writable from any task through the handle.
"""
mutable struct Pad <: AbstractDevice
    name::String
end
loop(::Pad, handle) = nothing

"""
    Panel(name)

A stub device declaring the calling-task affinity (§11.6): at most one holder
per roster, the calling task being a single-slot resource (§11.1, §11.3).
Its loop body is `Pad`'s immediate voluntary return, run inline on the
calling task by the wrapper.
"""
mutable struct Panel <: AbstractDevice
    name::String
end
needs_calling_task(::Panel) = true
loop(::Panel, handle) = nothing

"""
    Enumerated(faces...)

The returned claim source (§11.3): an input-side binding whose `claims` names
exactly the faces it was built with — the enumeration is the interface, and
the empty enumeration is the honest may-write-nothing degenerate, not a back
door (§11.6).
"""
struct Enumerated <: AbstractBinding
    faces::Vector{String}
end
Enumerated(faces::AbstractString...) = Enumerated(collect(String, faces))
is_input(::Enumerated) = true
claims(b::Enumerated) = b.faces

"""
    Greedy()

The computed claim source (§11.3): `is_greedy` declared, no `claims` of its
own — the framework computes the unclaimed complement at the attach point,
the shipped GUI binding's shape (§11.6, §11.7).
"""
struct Greedy <: AbstractBinding end
is_input(::Greedy) = true
is_greedy(::Greedy) = true

"""
    Readout(; label = selector, ...)

The output side's coverage binding (§11.2, §11.6): `reads` names exactly the
labeled selectors it was built with — the enumeration is the interface — and
its `map_output` is the identity, the wire datum being the gather's
NamedTuple itself. A snapshot-consuming telemetry peer differs only in what
its own `map_output` does with the same NamedTuple.
"""
struct Readout{R<:NamedTuple} <: AbstractBinding
    r::R
end
Readout(; sel...) = Readout(NamedTuple(sel))
is_output(::Readout) = true
reads(b::Readout) = b.r
map_output(nt, ::Readout) = nt
