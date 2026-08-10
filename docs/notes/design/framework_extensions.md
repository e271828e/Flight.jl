# Framework Extensions Charter

**Status: non-normative.** This document charters *future* design increments — worked
proposals, not decisions. Nothing here is part of the framework specification
(`framework_spec.md`, which wins on any conflict), nothing has decision-log rows, and
per the guarded-additions rule none of it is built until a demonstrated need arrives.
What this document preserves is the *analysis*: where the current design sits relative
to a general-purpose causal simulation engine, which gaps are reachable through seams
the spec deliberately left open, and — for the sample-time gaps — a worked proposal
detailed enough that a future increment can start from it rather than re-derive it.

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
  events, ZOH between ticks — is the core formalism ([§2](framework_spec.md#2-formalism)), with the ZOH delivered *by
  construction* through the interior/boundary sweep split ([§8.5](framework_spec.md#85-multi-rate-tick-scheduling)) rather than by
  runtime gating.
- **Event detection** with cheap boundary detection plus opt-in localization
  ([§2.1](framework_spec.md#21-events-two-detection-policies), [§8.4](framework_spec.md#84-localization-mechanics)) matches Simulink's zero-crossing machinery in semantics; the explicit
  two-policy split is arguably cleaner than a global zero-crossing option.
- **Virtual hierarchy**: assemblies flatten for scheduling and survive for navigation
  ([§3.3](framework_spec.md#33-assembly)) — exactly Simulink's virtual subsystem, and deliberately *only* that. [§8.5](framework_spec.md#85-multi-rate-tick-scheduling)
  records why atomic subsystems are an artificial-loop factory we refuse to
  reproduce.
- **FSMs** fall out of the mode/event facets of the continuous primitive ([§3.1](framework_spec.md#31-continuous-component-the-hybrid-primitive)) —
  Stateflow's core semantics without a second formalism.
- **Trim and linearization** ([§14.7](framework_spec.md#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)–[§14.10](framework_spec.md#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)) are engine-integrated, with exact AD
  Jacobians through the eltype-generic sweep ([§7.2](framework_spec.md#72-numeric-genericity-eltype)); perturbation-based `linmod` is
  strictly weaker.
- Guarantees the reference lacks: deterministic replay ([§2.2](framework_spec.md#22-exclusions-deliberate), [§9.5](framework_spec.md#95-inbound-the-input-trace)), the
  zero-allocation invariant ([§7.5](framework_spec.md#75-allocation-policy-a-scoped-invariant)), and the error-locality discipline of [§13](framework_spec.md#13-error-discipline).

### 1.2 Reachable through seams: the extension list

These gaps can be closed with limited, seam-respecting changes. Roughly in order of
increasing effort:

1. **Solver suite (adaptive, stiff, implicit).** The stepper seam ([§8.2](framework_spec.md#82-the-stepper-seam)) was shaped
   for exactly this: arbitrary `h`, dense output on demand, one-step methods only. An
   `OrdinaryDiffEq`-backed stepper as a package extension is already sanctioned prose
   in the spec. Adaptive stepping composes with the tick grid the same way Simulink
   handles sample-time hits: stretch steps between boundaries, land on them exactly.
   The multistep exclusion is a method-class restriction (event resets invalidate
   history), not a capability loss. Offline-mode only, which is fine.
2. **Sample-time phase offsets.** Recorded in [§8.5](framework_spec.md#85-multi-rate-tick-scheduling) as "no phase offsets in the first
   cut (no demonstrated use)". The gating machinery extends naturally; this is the
   subject of section 2 of this document.
3. **Non-trivial rate ratios.** Less of a gap than it looks: the harmonic grid
   already admits any *commensurable* rate set (30 Hz and 50 Hz coexist over a 150 Hz
   base as sibling divisors). What is awkward is the declaration ergonomics at the
   root — also the subject of section 2 of this document.
4. **Continuous transport delay.** The engine reference supports delay blocks with
   interpolated history; we have only the discrete `UnitDelay` ([§5.5](framework_spec.md#55-algebraic-loop-policy-reject-at-build-time), [§13.7](framework_spec.md#137-tooling-consequences-provenance-and-the-component-library)). A
   continuous transport delay is a *library component* with a ring buffer in `z` —
   it strains the immutable-`z` and allocation-policy corners but touches no kernel.
5. **Enabled-subsystem behavior** (hold or reset states while a control signal is
   low). The pattern already exists as mode logic plus reset events — the
   `PIVector(; reset = true)` shape from [§16](framework_spec.md#16-open-axes). What is missing is only sugar for
   applying it to a whole assembly at once; additive if ever wanted.
6. **Refined output between boundaries.** The seam's dense output exists for
   localization; exposing it to logging is additive, and only matters once adaptive
   steps exist.
7. **Triggered (edge-driven) discrete updates** — the borderline item, *moderate*
   rather than limited. A triggered FSM is already covered (continuous component
   with no `x`, events at boundary granularity). Running a discrete component's `g`
   off a signal edge rather than a timer is closer than it first appears — the due
   set is already per-boundary runtime data ([§8.5](framework_spec.md#85-multi-rate-tick-scheduling)), and the interior/boundary sweep
   split does not care *why* an entry is due — but the `Δt` story breaks: a
   triggered component has no period, and [§8.5](framework_spec.md#85-multi-rate-tick-scheduling) makes `Δt` a schedule-derived bundle
   field that discretized laws consume. It would need elapsed-time-since-last-firing
   semantics and a declaration surface for trigger wiring. Composable with the
   architecture, but real spec work, not a patch.

### 1.3 Closed axes: redesign-level, and deliberately so

These would require reversing load-bearing commitments. Each was rejected with
recorded rationale; this list exists so a future reader knows the closure was a
choice, not an oversight.

- **Algebraic loop solving.** The reference runs a per-step fixed-point/Newton
  iteration on instantaneous cycles; we reject them at build ([§5.5](framework_spec.md#55-algebraic-loop-policy-reject-at-build-time)). The divergence
  is structural: the entire [§5](framework_spec.md#5-evaluation-order-and-feedthrough) scheduling story — one topologically-ordered sweep
  over an immutable signal table, compiled statically, bounded per-step cost — rests
  on every cell being written exactly once per sweep. A loop solver needs *iterated*
  re-evaluation of a cut subgraph with mutable trial values and data-dependent
  convergence. It cannot be bolted on; it forks the execution model.
- **DAEs / implicit constraints.** Excluded in [§2.2](framework_spec.md#22-exclusions-deliberate). The sweep computes explicit
  `ẋ = f(...)`; implicit form changes the sweep contract itself (residual equations,
  algebraic variables — which are, definitionally, algebraic loops). The same fork as
  the previous item, larger.
- **Arbitrary/asynchronous sample times** via a time-ordered tick queue. Rejected in
  [§8.5](framework_spec.md#85-multi-rate-tick-scheduling); analyzed in depth in section 1.4 below, because it is the closed axis whose boundary
  is least obvious.
- **Variable-size signals / dynamic structure.** Fundamental conflict with the cell
  store, the flat state backing ([§7.1](framework_spec.md#71-continuous-state-structured-immutable-flat-backing)) and the compiled executor ([§12.7](framework_spec.md#127-the-compiled-executor)) — all sized
  at build — and with the zero-allocation invariant.
- **Shared mutable state** (Data Store Memory, Goto/From). Excluded on principle:
  no component reads another's state ([§3.2](framework_spec.md#32-periodic-discrete-component)), no shared mutable model ([§9.1](framework_spec.md#91-no-shared-mutable-model-staged-writes-snapshot-reads)).
  Emulation via ports is always available, and the reference's own documentation
  treats data stores as a diagnostics-laden escape hatch.

### 1.4 The closed axis worth understanding: arbitrary sample times

Under the harmonic grid, every tick instant lies on a **static lattice** known at
build time: `t = (k·K + φ)·Δt_base`. Everything downstream exploits this: the due set
is a pure function of the frame index (counter modulo, [§8.5](framework_spec.md#85-multi-rate-tick-scheduling)), frames have fixed
length for real-time pacing ([§8.7](framework_spec.md#87-real-time-pacing)), and the trace needs no tick timestamps because
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
than omission: asynchrony is an I/O phenomenon and lives at the periphery ([§9](framework_spec.md#9-runtime-periphery-the-data-plane)) —
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
were in fact re-derived on independent grounds: parent-relative rates ([§8.5](framework_spec.md#85-multi-rate-tick-scheduling)
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
extensions of section 2, and item 7 of section 1.2 would cover it. The load-bearing observation
throughout: **Stratum A flattens whatever was declared into the same `Build`
artifact** ([§12.1](framework_spec.md#121-three-strata)–[§12.2](framework_spec.md#122-the-build-artifact)) — resolved wires, absolute divisors, schedule — so the
sweep, the cell store and the executor never see how the model was written down.
These are authoring-surface increments, not execution-model changes.

---

## 2. Sample-time extensions: phase offsets and absolute rates

This section works out how items 2 and 3 of section 1.2 — phase offsets and convenient
non-trivial ratios — would be incorporated *within* the parent-relative sampling
paradigm, in enough detail that a future increment can start here. Everything below
stays on the static lattice: nothing from section 1.4's closed axis is dragged in.

### 2.1 Where the parent-relative paradigm actually hurts

Inside an authored assembly, relative integers are the right language: an inner/outer
cascade at 1:5 *is* a 1:5 design, and whoever writes the assembly owns that
arithmetic. The convention that makes `K ≥ 1` livable should be stated explicitly:
**a scope's base rate is its fastest member** (that member gets `K = 1`). With that
discipline, refactoring never fights the tree.

The degenerate case is the top. Suppose avionics must run at 50 Hz (the actuator bus
says so) and the nav filter at 30 Hz (the GPS says so). Today the user does the
arithmetic the framework should do:

```julia
# by hand: the common grid is 150 Hz, so...
rates(::Aircraft) = (avionics = 3, nav = 5)
Simulation(world; Δt_base = 1//150, h = 1//150)
```

Those `3` and `5` are not design ratios — nothing in the aircraft is "3:5"; they are
residue of a GCD computation. And they are *unstable* residue: change nav to 40 Hz
and the grid becomes 200 Hz, `avionics = 4`, `nav = 5`, `Δt_base = 1//200` — three
coordinated edits at two sites to express one changed fact. The pain is confined to
the boundary where relative structure meets absolute deployment facts, and the spec's
own doctrine already names that boundary: "absolute rates are deployment decisions
made once at the root."

### 2.2 Absolute rates at the root

The extension is not to abandon relative declaration but to let the root scope's rate
entries be **absolute, typed values**, and make the build do the numerology:

```julia
rates(::Aircraft) = (avionics = Hz(50), nav = Hz(30))   # or Period(1//50), etc.
Simulation(world; h = ...)   # Δt_base derived: gcd(1//50, 1//30) = 1//150
```

Mechanics:

- A plain integer in a `rates` entry keeps today's meaning (multiplier relative to
  the enclosing scope). A `Hz`/`Period` value declares an absolute rate for that
  subtree scope; the wrapper type distinguishes the two registers — no mode flag, no
  ambiguity.
- Periods are **exact rationals** (`Rational` seconds or integer Hz), never floats:
  GCD over floats is ill-defined, and `Rational` makes the derivation exact and the
  diagnostics printable.
- At deployment binding, `Δt_base` is derived as the GCD of all absolute
  declarations — or validated, if the user declared one; a declared `Δt_base` that
  fails to divide some declared period (or offset, section 2.3) is a `DeploymentInvalid`,
  collected per [§12.1](framework_spec.md#121-three-strata)'s style. Divisors for absolutely-declared scopes are computed
  there, and relative multipliers below them compose as they always did.
- The one structural shift is small and localized: divisor resolution for absolute
  entries moves from Stratum A to deployment binding, since it needs `Δt_base`.
  Everything downstream is indifferent — the executor's static split is
  continuous-vs-discrete membership, and divisors are values it gates by at runtime.

Every property the current design was chosen for survives. Library assemblies still
carry their ratios in their types and remain deployment-portable. A component author
still *cannot* know their absolute rate — absolute declarations exist only at the
root, which is deployment's territory anyway — so [§8.5](framework_spec.md#85-multi-rate-tick-scheduling)'s structural argument against
caching `Δt`-derived coefficients is untouched. And the top-level integers that were
never ratios stop existing.

If a deployment needs to pin a rate *deep* inside a library assembly ("the EKF must
be 30 Hz regardless of what encloses it"), the disciplined answer is the one the spec
uses for every other internal degree of freedom: the assembly exposes that `K` as a
constructor parameter. Absolute pinning from outside a subtree's contract would be
action at a distance.

### 2.3 Phase offsets and the composition law

Offsets come in the same two registers.

**Relative (interior) form.** An instantiation carries a pair instead of a bare
multiplier: `(K, φ)` with `0 ≤ φ < K`, meaning "fire on every K-th tick of my scope,
starting from its φ-th" — scope-relative, exactly like the multiplier, so a library
assembly can stagger its internal children without knowing any absolute rate.

The reason this is clean and not merely plausible is that the composition law works
out affinely. Write a scope's compiled schedule as an absolute divisor and phase
`(D_s, Φ_s)`: it ticks at base indices `n = Φ_s + m·D_s`. A child declared
`(K_c, φ_c)` inside it fires on scope ticks with `m ≡ φ_c (mod K_c)`, so its base
indices are

```
n = Φ_s + (φ_c + l·K_c)·D_s  =  (Φ_s + φ_c·D_s) + l·(K_c·D_s)
```

That is: `D_child = K_c·D_s` (the existing multiplicative rule, untouched) and
`Φ_child = Φ_s + φ_c·D_s`. Multipliers compose multiplicatively, phases affinely, and
both collapse at build into one `(D, Φ)` pair per discrete component — the same "all
scoping compiles away" shape [§8.5](framework_spec.md#85-multi-rate-tick-scheduling) already has. The boundary sweep's gating test
changes from `idx % D == 0` to `(idx − Φ) % D == 0`. One subtraction; the lattice
stays static.

**Absolute (root) form.** `Hz(f; offset = τ)` / `Period(T; offset = τ)` declares
tick instants `t = τ + k·T`, with `0 ≤ τ < T` and τ an exact rational time. Since
ticks may only land on base ticks, **absolute offsets join the GCD pool alongside
periods** when `Δt_base` is derived. This is the one subtlety with teeth: an offset
of `T/2` costs a 2× finer grid; an offset of `T/1000` costs 1000×. Offsets should be
simple fractions of their period.

**Spec amendments this would require**, all small but all honest:

- *Boundary zero.* [§8.5](framework_spec.md#85-multi-rate-tick-scheduling)'s "at boundary zero the due set is everything" refines to
  "everything with `Φ = 0`": an offset component's first tick is at `Φ·Δt_base`, and
  until then its cells hold the values the build probe populated from `init_z`
  ([§12.3](framework_spec.md#123-probing-and-input-synthesis)) — a coherent ZOH story, since those are exactly the values a tick at
  `t₀⁻` would have produced. [§14.5](framework_spec.md#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)'s boundary-zero transitions and the trim
  read-back need a sentence each acknowledging offset components are not due there.
- *`Δt` semantics.* Unchanged: the bundle's `Δt` is still `D·Δt_base`; the offset
  shifts firing instants, not the period.
- *Documentation next to [§8.5](framework_spec.md#85-multi-rate-tick-scheduling)'s simultaneous-tick rule.* The choice between
  coincident ticks (fresh same-instant reads via topological order) and staggered
  ticks (pipelined, deterministically aged reads) is a modeling decision with
  observable consequences — see section 2.4 — and should be documented where a user
  debugging "why does my controller see stale data" will find it.

### 2.4 Worked example: a staggered sensor suite

The sensor suite runs at a 0.002 s period with a 0.001 s offset; the control laws at
0.02 s, offset 0:

```julia
rates(::Vehicle) = (
    sensors = Period(1//500,  offset = 1//1000),   # 0.002 s, offset 0.001 s
    ctrl    = Period(1//50),                        # 0.02 s
)
```

**Grid derivation.** The GCD pool is {1/500, 1/50, 1/1000}, giving
`Δt_base = 1//1000`. Note that neither period alone needs a 1 kHz grid —
`gcd(1/500, 1/50) = 1/500` — it is the offset, half the sensor period, that refines
it. Deployment then needs `h ≤ 1 ms`; take `h = 1 ms`, `n = 1`.

**Compiled pairs**: `sensors → (D = 2, Φ = 1)` (due at odd base indices);
`ctrl → (D = 20, Φ = 0)` (due at multiples of 20). Timeline, with base index in
milliseconds:

```
t (ms):    0   1   2   3   4  ...  19  20  21  ...
ctrl:      C                           C
sensors:       S       S      ...  S       S
```

The two due sets are *disjoint at every boundary* — the designed effect of a
half-period stagger. Boundary by boundary:

- **t = 0**: due = {ctrl} (the amended boundary-zero rule). The sensors have never
  ticked; ctrl reads their build-probed `init_z` cells for the first millisecond.
- **Odd ms**: due = {sensors}; ctrl's cells hold.
- **Even ms, not multiples of 20**: due set empty — the boundary exists for events
  and frame structure only. These 1 kHz boundaries in a model whose fastest task is
  500 Hz are the overhead the offset-refined grid charges.
- **Multiples of 20 ms**: due = {ctrl}, whose `h_zu` reads sensor cells written at
  the previous odd boundary — data exactly 1 ms old, every time, deterministically.

**What the offset bought.** With `offset = 0`, sensors and ctrl coincide at every
20 ms boundary and [§8.5](framework_spec.md#85-multi-rate-tick-scheduling)'s simultaneous-tick rule gives ctrl same-instant sensor data
via topological order — the idealized synchronous-sampling picture. With the offset,
that becomes "pipelined, deterministically 1 ms stale" — which is precisely the knob:
a real ADC-plus-bus pipeline delivers data some fixed fraction of a sample period
old, and the offset expresses that latency *structurally*, in the schedule, rather
than through delay blocks. The second classic use is load shaping under real-time
pacing ([§8.7](framework_spec.md#87-real-time-pacing)): with the stagger, no single 1 ms frame executes both the sensor stack
and the control stack, so worst-case frame cost is the `max` of the two rather than
their sum.

### 2.5 Cascading into a scope: relative offsets at work

Continue the example into the `sensors` assembly:

```julia
rates(::Sensors) = (imu = 1, gnss = (20, 3))
```

The scope arrives compiled as `(D_s, Φ_s) = (2, 1)`: its grid is `t = 1, 3, 5, … ms`,
and that pair is what the children inherit — a child never sees seconds, only scope
ticks. Applying the composition law:

- **`imu`**: `(K=1, φ=0)` → `D = 2, Φ = 1`. Every scope tick: 1, 3, 5, … ms — the
  scope's declared rate *is* its fastest member, per the section 2.1 convention.
- **`gnss`**: `(K=20, φ=3)` → `D = 40, Φ = 1 + 3·2 = 7`. Ticks at 7, 47, 87, … ms —
  25 Hz, on the 3rd scope tick of each cycle of 20.

The `Sensors` author declared ratios and a tick-count stagger; they did not know —
could not know, preserving the [§8.5](framework_spec.md#85-multi-rate-tick-scheduling) property — that a scope tick would turn out to be
2 ms wide or displaced by 1 ms. Three properties of the stacking are worth making
explicit, because they are what distinguish the relative register from the absolute
one:

1. **Relative offsets never refine the base grid.** `φ` is an integer count of scope
   ticks — it selects among instants that already exist. Grid cost is confined
   entirely to the absolute register at the root. This asymmetry is the answer to
   "why not let everyone declare offsets in seconds."
2. **Relative offsets cannot leave the scope grid.** `gnss` can fire on the 3rd or
   19th scope tick of a cycle, never *between* scope ticks. Staggering off-grid
   requires expressing the offset where the grid is finer: at the root in absolute
   time (paying refinement there), or by declaring the scope base finer than its
   fastest member so unused slots exist.
3. **`K = 1` admits no stagger** (`0 ≤ φ < K` forces `φ = 0`): a component occupying
   its scope's full grid cannot be phase-shifted against it. Two same-rate siblings
   are staggered by giving the scope twice their rate and declaring `(2, 0)` and
   `(2, 1)` — the same "stagger costs resolution" rule one level down.

The merged schedule: at 7, 47, 87 … ms both `imu` and `gnss` are due — gnss's
instants are a subset of the scope grid, which imu fully occupies, so they always
coincide, and the simultaneous-tick rule orders them topologically within one sweep
(if imu feeds gnss, gnss reads same-boundary fresh imu output: a clean fusion cascade
with zero intra-suite latency). Meanwhile `ctrl`, reading at 20, 40, 60, 80 … ms,
sees gnss data aged alternately 13 ms and 33 ms (writes at 7, 47, …) — not an
artifact but the honest schedule of a multi-rate system, and *fully determined at
build* from the `(D, Φ)` table. That determinism suggests a diagnostic: the [§5.6](framework_spec.md#56-diagnostics-feedthrough-tracing)
feedthrough tracer could print per-wire data ages ("this connection delivers data
aged 13–33 ms"), turning the staleness question from a debugging session into a
table lookup.

---

## Addendum A. The `Group` component: on-the-fly assemblies

*(Parked here pending a proper home — likely the [§13.7](framework_spec.md#137-tooling-consequences-provenance-and-the-component-library) component library's starting
inventory.)*

The spec rejected the mutable builder (`Assembly()` + `add!`/`connect!`) in [§11.5](framework_spec.md#115-assembly-declaration-type-based-class-by-declaration-shape) on
grounds that survive any reframing: the dispatched-on type and the recipe defining
its structure drift apart, mutable state threads through declaration code, and it
does not even capture source locations. But the *immutable* version of "grouping
components by plain calls" needs no builder — it is already expressible inside the
current spec as a single library component, because [§11.5](framework_spec.md#115-assembly-declaration-type-based-class-by-declaration-shape)'s container-children rule
makes a `NamedTuple` field contribute its elements as path-named children, and
declarations are ordinary functions of the *instance*, free to read its fields:

```julia
struct Group{C <: NamedTuple, W, I, O} <: AbstractComponent
    children::C      # component-typed elements → children by the container rule
    wires::W         # inert parameter data
    inputs::I
    outputs::O
end
child_connections(g::Group)  = g.wires
input_connections(g::Group)  = g.inputs
output_connections(g::Group) = g.outputs

world = Group(
    (; plant = Plant(), ctrl = PID(kp = 2.0)),
    (("ctrl/u", "plant/u"), ("plant/y", "ctrl/y")),
    (;),
    (;),
)
```

One type, defined once; every ad-hoc topology is a *value* of it. The type parameters
still carry the children's concrete types, so Stratum C specialization and the
compiled executor work unchanged; wiring validation, did-you-mean errors and the
two-producer check all run at build against the instance exactly as for a named
assembly. What is given up relative to a named type is exactly what named types are
*for* — dispatching domain code on `::Cessna172X`, a reusable identity for the
topology — which exploratory work does not want anyway.

The framing that earns `Group` its place: the spec's bias was never the type-based
*semantics* (sound on general-purpose grounds too) but making *named* types the only
spelled-out route. A general-purpose engine ships `Group` in its standard library on
day one, the way Julia ships anonymous functions alongside named ones — serving the
model-assembler persona of section 1.5 with a library addition and perhaps a paragraph in
[§11.5](framework_spec.md#115-assembly-declaration-type-based-class-by-declaration-shape) acknowledging the pattern. No spec surgery.
