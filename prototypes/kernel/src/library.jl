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
A mode-carrying source: modes are continuous-only, and in this increment they
are read but never written — handlers are the event system, which is absent.
The branch returning a literal is the constant-branch idiom (D-166).
"""
struct ModedSource <: AbstractComponent end

init_m(::ModedSource) = (phase = :idle,)
output_types(::ModedSource, ::Type{T}) where {T <: Real} = (out = T,)

h_x(::ModedSource, (; m)) = (out = m.phase === :idle ? 0.0 : 1.0,)

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
          # face, and at the root a face is a slot — its initial value
          # synthesized by `probe_value`, written by `set_slot!` (§6.1, §11.3).
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
declared fields — so an ancestor may deep-route into them (§6.1).

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
output_connections(::SampledLoop) = ("plant/y" => "y", "ctl/u" => "cmd")
sample_times(l::SampledLoop) = (ctl = l.ctl_rate,)

"""
    Vehicle(; k, kI, ω, ζ)

Two levels: the vehicle scales its own `"ref"` face and hands the product to the
loop's input face, and re-exports three of the loop's signals — two of them
through the loop's own faces, the third by deep-routing to `"loop/plant/power"`,
a grandchild's port reached through concretely declared fields.

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
    ("loop/y" => "y", "loop/cmd" => "cmd", "loop/plant/power" => "power")

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
