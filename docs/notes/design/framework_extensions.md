# Framework Extensions Charter

**Status: non-normative.** This document charters *future* design increments — worked
proposals, not decisions. Nothing here is part of the framework specification
(`framework_spec.md`, which wins on any conflict), nothing has decision-log rows, and
per the guarded-additions rule none of it is built until a demonstrated need arrives.
What this document preserves is the *analysis*: where the current design sits relative
to a general-purpose causal simulation engine, and which gaps are reachable through
seams the spec deliberately left open. For the sample-time gaps, the worked proposal
lives in its own companion, [`sample_time_proposal.md`](sample_time_proposal.md)
(adopted into the spec 2026-08-12, decision [D-185][d-185]–[D-187][d-187]).

Provenance: distilled from the gap-analysis discussion of 2026-08-08/09.

---

## 1. The gap survey: how far from a general-purpose causal engine?

The framework began as a rewrite of `FlightCore`, and its early decisions were shaped
by what porting `FlightPhysics` and `FlightApps` would need. The question this section
answers is how far the resulting design actually sits from a *general-purpose* causal
simulation engine — the reference point being Simulink's engine proper (its solver
suite, sample-time engine, conditional execution, algebraic-loop solver, zero-crossing
detection, signal engine and simulation control), not the ecosystem around it.

The short answer, developed below: the design is best described as **Simulink's
fixed-step, loop-free, harmonic-multirate causal subset, plus guarantees Simulink does
not offer** (deterministic bit-exact replay, a zero-allocation invariant, exact AD
Jacobians through the sweep, build-time error locality). The distance to full
generality is concentrated in two commitments — no algebraic-loop solving, no
implicit/DAE form — that are architectural, not incidental. Everything else is either
already covered, reachable through deliberate seams, or a library matter.

### 1.1 Already at parity or ahead

A surprising amount of the engine-level surface is covered, sometimes more rigorously
than the reference covers it:

- **Hybrid execution** — continuous flow, multi-rate discrete updates, zero-crossing
  events, ZOH between ticks — is the core formalism ([§2][s2]), with the ZOH delivered *by
  construction* through the interior/boundary sweep split ([§10.5][s10-5]) rather than by
  runtime gating.
- **Event detection** with cheap boundary detection plus opt-in localization
  ([§2.1][s2-1], [§10.4][s10-4]) matches Simulink's zero-crossing machinery in semantics; the explicit
  two-policy split is arguably cleaner than a global zero-crossing option.
- **Virtual hierarchy**: assemblies flatten for scheduling and survive for navigation
  ([§3.3][s3-3]) — exactly Simulink's virtual subsystem, and deliberately *only* that. [§10.5][s10-5]
  records why atomic subsystems are an artificial-loop factory we refuse to
  reproduce.
- **FSMs** fall out of the mode/event facets of the continuous primitive ([§3.1][s3-1]) —
  Stateflow's core semantics without a second formalism.
- **Trim and linearization** ([§14.7][s14-7]–[§14.10][s14-10]) are engine-integrated, with exact AD
  Jacobians through the eltype-generic sweep ([§7.2][s7-2]); perturbation-based `linmod` is
  strictly weaker.
- Guarantees the reference lacks: deterministic replay ([§2.2][s2-2], [§11.5][s11-5]), the
  zero-allocation invariant ([§7.5][s7-5]), and the error-locality discipline of [§13][s13].

### 1.2 Reachable through seams: the extension list

These gaps can be closed with limited, seam-respecting changes. Roughly in order of
increasing effort:

1. **Solver suite (adaptive, stiff, implicit).** The stepper seam ([§10.2][s10-2]) was shaped
   for exactly this: arbitrary `h`, dense output on demand, one-step methods only. An
   `OrdinaryDiffEq`-backed stepper as a package extension is already sanctioned prose
   in the spec. Adaptive stepping composes with the tick grid the same way Simulink
   handles sample-time hits: stretch steps between boundaries, land on them exactly.
   The multistep exclusion is a method-class restriction (event resets invalidate
   history), not a capability loss. Offline-mode only, which is fine.
2. **Sample-time phase offsets.** Recorded in [§10.5][s10-5] as "no phase offsets in the first
   cut (no demonstrated use)". The gating machinery extends naturally; worked out in
   [`sample_time_proposal.md`](sample_time_proposal.md).
3. **Non-trivial rate ratios.** Less of a gap than it looks: the harmonic grid
   already admits any *commensurable* rate set (30 Hz and 50 Hz coexist over a 150 Hz
   base as sibling divisors). What is awkward is the declaration ergonomics at the
   root — also worked out in [`sample_time_proposal.md`](sample_time_proposal.md).
4. **Continuous transport delay.** The engine reference supports delay blocks with
   interpolated history; we have only the discrete `UnitDelay` ([§5.5][s5-5], [§13.7][s13-7]). A
   continuous transport delay is a *library component* with a ring buffer in its
   discrete state `x` —
   it strains the immutable-state and allocation-policy corners but touches no kernel.
5. **Enabled-subsystem behavior** (hold or reset states while a control signal is
   low). The pattern already exists as mode logic plus reset events — the
   `PIVector(; reset = true)` shape from [§16][s16]. What is missing is only sugar for
   applying it to a whole assembly at once; additive if ever wanted.
6. **Refined output between boundaries.** The seam's dense output exists for
   localization; exposing it to logging is additive, and only matters once adaptive
   steps exist.
7. **Triggered (edge-driven) discrete updates** — the borderline item, *moderate*
   rather than limited. A triggered FSM is already covered (continuous component
   with no `x`, events at boundary granularity). Running a discrete component's `g`
   off a signal edge rather than a timer is closer than it first appears — the due
   set is already per-boundary runtime data ([§10.5][s10-5]), and the interior/boundary sweep
   split does not care *why* an entry is due — but the `Δt` story breaks: a
   triggered component has no period, and [§10.5][s10-5] makes `Δt` a schedule-derived bundle
   field that discretized laws consume. It would need elapsed-time-since-last-firing
   semantics and a declaration surface for trigger wiring. Composable with the
   architecture, but real spec work, not a patch.

### 1.3 Closed axes: redesign-level, and deliberately so

These would require reversing load-bearing commitments. Each was rejected with
recorded rationale; this list exists so a future reader knows the closure was a
choice, not an oversight.

- **Algebraic loop solving.** The reference runs a per-step fixed-point/Newton
  iteration on instantaneous cycles; we reject them at build ([§5.5][s5-5]). The divergence
  is structural: the entire [§5][s5] scheduling story — one topologically-ordered sweep
  over an immutable signal table, compiled statically, bounded per-step cost — rests
  on every cell being written exactly once per sweep. A loop solver needs *iterated*
  re-evaluation of a cut subgraph with mutable trial values and data-dependent
  convergence. It cannot be bolted on; it forks the execution model.
- **DAEs / implicit constraints.** Excluded in [§2.2][s2-2]. The sweep computes explicit
  `ẋ = f(...)`; implicit form changes the sweep contract itself (residual equations,
  algebraic variables — which are, definitionally, algebraic loops). The same fork as
  the previous item, larger.
- **Arbitrary/asynchronous sample times** via a time-ordered tick queue. Rejected in
  [§10.5][s10-5]; analyzed in depth in section 1.4 below, because it is the closed axis whose boundary
  is least obvious.
- **Variable-size signals / dynamic structure.** Fundamental conflict with the cell
  store, the flat state backing ([§7.1][s7-1]) and the compiled executor ([§9.7][s9-7]) — all sized
  at build — and with the zero-allocation invariant.
- **Shared mutable state** (Data Store Memory, Goto/From). Excluded on principle:
  no component reads another's state ([§3.2][s3-2]), no shared mutable model ([§11.1][s11-1]).
  Emulation via ports is always available, and the reference's own documentation
  treats data stores as a diagnostics-laden escape hatch.

### 1.4 The closed axis worth understanding: arbitrary sample times

Under the harmonic grid, every tick instant lies on a **static lattice** known at
build time: `t = (k·K + φ)·Δt_base`. Everything downstream exploits this: the due set
is a pure function of the frame index (counter modulo, [§10.5][s10-5]), frames have fixed
length for real-time pacing ([§10.7][s10-7]), and the trace needs no tick timestamps because
tick times are derivable. Ratios and offsets *enlarge* the lattice; they never change
its nature.

The tick-queue model abandons the lattice: tick instants become entries in a runtime
priority queue — pop the earliest, integrate exactly to it with a variable-length
step, service it, push the next. The next tick time can be computed at runtime, by
anything. The essential difference is **tick times as build-time constants versus
tick times as runtime data**. This buys three distinct things:

- **(a) Awkward period sets without a tiny forced step.** The harmonic grid's real
  cost is not "irrational ratios are forbidden" (irrelevant in practice) but "the
  base tick is the GCD of all periods, and the integrator inherits it": periods of
  0.1 s and 0.1001 s force `Δt_base = 100 µs` and hence `h ≤ 100 µs` — a 1000×
  stepping overhead paying for nothing. A tick queue adapts the step per interval and
  needs no common base.
- **(b) State-dependent tick instants.** The canonical case is crank-angle-domain
  engine control: a task that runs every 120° of crankshaft rotation has a period
  that varies continuously along the trajectory. No fixed lattice contains those
  instants, because they are not known until the state that determines them has been
  integrated.
- **(c) Asynchronous tasks**: independent clocks with no fixed phase relation
  (interrupt-driven boxes, drifting sensor clocks). Tick times are not merely
  runtime-computed but external to model time.

How much of this is genuinely lost? **(c)** the design answers by position rather
than omission: asynchrony is an I/O phenomenon and lives at the periphery ([§11][s11]) —
staged writes entering at boundaries, timestamped in the trace — keeping
model-internal time deterministic and replayable. **(b)** is mostly recoverable at
boundary resolution with machinery we have: a guard on the wrapped crank angle is a
zero-crossing event (localizable when timing precision matters), and firing a
discrete update off it is precisely item 7 of section 1.2 — a self-triggered component *is*
a variable-sample-time block, quantized to step boundaries unless localized. **(a)**
is the only unrecoverable item, and it is an efficiency corner case, not a semantic
one. This is why the axis stays closed: ratios + offsets exhaust everything a static
schedule can express, and of the queue's three exclusive uses, one is better solved
at the periphery, one is nearly reachable via events, and one is rare.

### 1.5 The bias, named

Re-examination showed the two decisions most suspected of `FlightCore` inheritance
were in fact re-derived on independent grounds: parent-relative rates ([§10.5][s10-5]
explicitly identifies `Subsampled`'s parent-relative multipliers as a call-tree
artifact, then re-chooses relative declaration because *ratios* are intrinsic to a
control architecture while absolute rates are deployment decisions) and virtual
assemblies (a departure from FlightCore's execution-through-composition, coinciding
with the reference's own virtual-subsystem default).

Where the `FlightPhysics`/`FlightApps` assumption genuinely biased the design is the
**authoring layer's single persona**: the *library author*, for whom everything is a
named type, structure is code, and rates are ratios baked into definitions. The
persona the assumption under-serves is the *model assembler* — someone composing
topologies and rate maps as data, interactively or programmatically. Serving the
second persona requires no change to the execution semantics and no reversal of any
recorded rejection: a `Group` library type (Addendum A), the rate-declaration
extensions of `sample_time_proposal.md`, and item 7 of section 1.2 would cover it. The load-bearing observation
throughout: **Stratum A flattens whatever was declared into the same `Build`
artifact** ([§9.1][s9-1]–[§9.2][s9-2]) — resolved wires, absolute divisors, schedule — so the
sweep, the cell store and the executor never see how the model was written down.
These are authoring-surface increments, not execution-model changes.

---

## 2. Sample-time extensions: phase offsets and absolute rates

**Superseded (2026-08-12).** The worked proposal for items 2 and 3 of section 1.2 —
phase offsets and absolute rate declarations — now lives in its own companion,
[`sample_time_proposal.md`](sample_time_proposal.md), which absorbed this section's
material (the root-declaration pain, the affine composition law, the register
asymmetry, the staggered-sensor-suite example) and extended it with the
`Relative`/`Absolute` declaration surface, mid-tree anchors, the full compilation
pipeline (the `(anchor, m, c)` intermediate, the derive-vs-declare rule for
`Δt_base`, the runtime gate), a nested-anchor worked example, and the grid
diagnostics. Everything there still stays on the static lattice: nothing from
section 1.4's closed axis is dragged in.

---

## 3. Parameter sensitivities: walking the parameter fields

Provenance: distilled from the parameter-AD discussion of 2026-08-12. Status as
everything here: analysis preserved, nothing built, guarded-additions rule applies.

### 3.1 The capability and its consumers

[§7.2][s7-2] scopes the eltype-generic surface deliberately: payload/value types are
**walked** (they follow the activation scalar), parameters are **pinned** ("stay
`Float64`; promotion handles mixing"), and [§14.10][s14-10] makes differentiation
participation a per-invocation *seeding* fact over `x`-taps and slots. The
consequence: a component parameter can only ever enter an evaluation as a
zero-partial constant. That is correct for every current service — in `∂ẋ/∂x` and
`∂ẋ/∂u` the parameter genuinely is a constant — but structurally incapable of being
the seeded direction itself. Forward-mode differentiation with respect to a quantity
means that quantity carries the unit partial, and a `k::Float64` field cannot hold a
`Dual`: the wall is at construction, before any math runs.

The extension admits parameter fields to the walked set by making them
type-parametric —

```julia
struct Gain{T <: Real}   # Float64 nominally; Dual under a seeded build
    k::T
end
```

— which enables exact `∂(anything on the continuous path)/∂k`: `ẋ` leaves, outputs,
trim residuals ([§14.7][s14-7]), and — with the section 3.4 caveat — trajectory values. Three
consumers this project would actually use: **parametric linearization**
(`∂A/∂(aero coefficient)`; `∂(trim point)/∂(mass, CG)` via the implicit function
theorem on the trim residual Jacobian), **gradient-based tuning** (FCS gains against
a closed-loop cost, exact gradients instead of one finite difference per gain), and
**parameter estimation** (fitting to flight data needs `∂(trajectory)/∂(parameters)`).

Not a Simulink-parity item: this is an axis where section 1.1 already places the
design ahead (exact AD Jacobians through the sweep); the extension widens that
surface rather than closing a gap.

### 3.2 Mechanism: the seed moves from a tap into the world value

Everything rides seams already in place; no new evaluation machinery. Mechanism,
verified end to end in plain Julia (2026-08-12):

```julia
struct Actuator{G,L}     # parent generic over children — the §8.5 shape
    gain::G
    lag::L
end
f(a::Actuator, x, u) = (a.gain.k * u - x) / a.lag.τ   # ẋ = (k·u − x)/τ

kd = Dual{T}(2.0, 1.0)                 # the parameter: unit partial
a  = Actuator(Gain(kd), Lag(0.5))      # the rest: plain Float64
xd = Dual{T}(1.0, 0.0)                 # state at the operating point
ud = Dual{T}(3.0, 0.0)                 # slot value
partials(f(a, xd, ud), 1)              # ∂ẋ/∂k = u/τ = 6.0, exact
```

- **Seeding.** Construct the target component around a unit-partial `Dual` and build
  the world from it; states and slots enter as zero-partial constants. This is
  [§14.10][s14-10]'s seeded pass verbatim, with the unit partial living in a field instead of
  an `x`-tap.
- **Mixed construction is semantically right, not merely convenient.** Unseeded
  components keep plain `Float64` parameters; promotion lifts each one into a
  zero-partial constant on contact with `Dual` traffic — the mathematically correct
  role for a parameter not being differentiated. [§7.2][s7-2]'s "promotion handles mixing"
  clause is already doing exactly this job.
- **Parents must be generic over child types.** `Actuator{G,L}` absorbs a
  `Gain{Dual}` child; a parent pinning `gain::Gain{Float64}` re-erects the wall one
  level up (loudly, at assembly construction). [§8.5][s8-5]'s type-based assembly shape
  already has this property — it is [§7.2][s7-2]'s "no `::SomeType{Float64}` annotations"
  rule transposed from method signatures to field declarations.
- **Parameters mix; storage does not.** Cell types and the `x`-buffer eltype are
  declarations *evaluated at the activation scalar* ([§7.2][s7-2], [§8.2][s8-2]), never inferred
  from traffic. A seeded world therefore runs under the matching `Dual` activation:
  homogeneous `Dual` storage, heterogeneous parameters. Under the nominal `Float64`
  activation the first `Dual` product to reach a `Float64` cell throws — the same
  loud-failure channel as [§9.4][s9-4]'s probe. Corollary: parameter duals carry the
  framework's canonical probe tag family, not private tags.
- **Nominal cost is zero.** At the `Float64` activation a parametric field
  specializes to the same concrete struct; `@kwdef` defaults pin the no-argument
  case — the pattern [§7.2][s7-2] already prescribes for the walked payload types. The
  pinned-parameters stance stays the right nominal default; this extension is
  additive.

### 3.3 Economics: seeding patterns are world types

`Actuator{Gain{Dual}, Lag{Float64}}` and `Actuator{Gain{Float64}, Lag{Dual}}` are
distinct types — one build, and one compilation of everything downstream, per choice
of target component. Fine for a one-off study. For scanning many parameters the
amortized shape is: construct **all** candidate parameters as `Dual`s — one world
type, one compilation — and select the differentiated subset by *value* (unit
vs. zero partials), chunked exactly as [§14.10][s14-10] chunks its tap directions. Which
parameters carry partials is a value fact, not a type fact.

### 3.4 Boundaries

- **Discrete-side parameters are excluded, and correctly so.** Under frozen-exact, a
  digital compensator coefficient's influence on the continuous state is temporal,
  not instantaneous; its point-wise partial is genuinely zero. Differentiating with
  respect to it means differentiating the sampled-data step map — [§14.10][s14-10]'s own
  "recorded, not built" extension. This section is that seam's parameter-flavored
  companion: if the sampled-data `Dual` activation is ever built, walked parameter
  fields on participating discrete components ride along.
- **Trajectory sensitivities through events need jump corrections.** Point-wise
  Jacobians are unconditionally clean, and `Dual` propagation through an integrator
  is exact along smooth flow (verified: the lag's `∂x/∂k` converges to the analytic
  steady-state sensitivity). But at an event crossing the crossing *time* depends on
  the parameter, and naive propagation misses that term; the standard sensitivity
  jump correction is research-grade machinery deliberately out of scope — consistent
  with [§10.4][s10-4]'s rejection of Newton/AD localization on C⁰ grounds.

### 3.5 What adoption would still need

The unresolved question is the user-facing surface, and per the guarded-additions
rule it stays unresolved until a consumer arrives: how a user *asks* for `∂y/∂k`.
Candidates, none worked out: a parameter-tap register extending [§14.10][s14-10]'s selector
family, with the framework performing the seeded construction internally; a
documented construct-with-`Dual`s convention (the capability with no service
wrapper); parameter-path overlays in [§14.1][s14-1]'s condition system. Also unexamined:
which library structs get the parametric spelling in the [§13.7][s13-7] migration — whether
[§7.2][s7-2]'s mechanical parametrization of the walked payload list simply extends to
parameter carriers, or participation is opted into per component.

---

## Addendum A. The `Group` component: on-the-fly assemblies

**Folded into the spec (2026-08-12, [D-184][d-184]).** `Group` is now normative: the
pattern, its sketch and the anonymous-functions framing live in [§8.5][s8-5],
and the component joins [§13.7][s13-7]'s starting inventory as its one
persona-admitted member. Nothing of the addendum's argument was lost in the
move; this heading remains so section 1.5's persona discussion keeps its
pointer.

<!-- citation link definitions — generated by tools/linkify.jl; do not edit -->
[d-184]: framework_decisions.md#d-184--fold-group-into-the-spec-as-an-ordinary-library-component
[d-185]: framework_decisions.md#d-185--adopt-the-phased-two-register-sample-time-declaration
[d-187]: framework_decisions.md#d-187--make-the-bound-schedule-a-named-artifact-with-exact-grid-diagnostics
[s10-2]: framework_spec.md#102-the-stepper-seam
[s10-4]: framework_spec.md#104-localization-mechanics
[s10-5]: framework_spec.md#105-multi-rate-tick-scheduling
[s10-7]: framework_spec.md#107-real-time-pacing
[s11]: framework_spec.md#11-runtime-periphery-the-data-plane
[s11-1]: framework_spec.md#111-no-shared-mutable-model-staged-writes-snapshot-reads
[s11-5]: framework_spec.md#115-inbound-the-input-trace
[s13]: framework_spec.md#13-error-discipline
[s13-7]: framework_spec.md#137-tooling-consequences-provenance-and-the-component-library
[s14-1]: framework_spec.md#141-conditions-are-path-addressed-overlays-on-the-declared-defaults
[s14-10]: framework_spec.md#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query
[s14-7]: framework_spec.md#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals
[s16]: framework_spec.md#16-open-axes
[s2]: framework_spec.md#2-formalism
[s2-1]: framework_spec.md#21-events-two-detection-policies
[s2-2]: framework_spec.md#22-exclusions-deliberate
[s3-1]: framework_spec.md#31-continuous-component-the-hybrid-primitive
[s3-2]: framework_spec.md#32-periodic-discrete-component
[s3-3]: framework_spec.md#33-assembly
[s5]: framework_spec.md#5-evaluation-order-and-feedthrough
[s5-5]: framework_spec.md#55-algebraic-loop-policy-reject-at-build-time
[s7-1]: framework_spec.md#71-continuous-state-structured-immutable-flat-backing
[s7-2]: framework_spec.md#72-numeric-genericity-eltype
[s7-5]: framework_spec.md#75-allocation-policy-a-scoped-invariant
[s8-2]: framework_spec.md#82-the-declaration-inventory
[s8-5]: framework_spec.md#85-assembly-declaration-type-based-class-by-declaration-shape
[s9-1]: framework_spec.md#91-three-strata
[s9-2]: framework_spec.md#92-the-build-artifact
[s9-4]: framework_spec.md#94-activations-executable-sets-laziness-caching
[s9-7]: framework_spec.md#97-the-compiled-executor
