# A Modeling & Simulation Framework for Flight.jl — Design Document

**Status:** first checkpoint (v0.1). Axes 1–4 settled; axes 5–7 and the migration
outline pending (see [Open axes](#open-axes)).

---

## 1. Purpose and method

This document specifies the design of a modeling and simulation framework intended to
replace `FlightCore` as the substrate for `FlightPhysics` and `FlightApps`. The new
framework must match or surpass `FlightCore` in functionality, performance and
flexibility, while being more rigorous and explicit — reducing the learning curve and
the number of latent footguns for model authors.

Ground rules adopted for this design effort:

- **Independence.** This is one of (at least) two design solutions developed from
  scratch and in isolation from each other, to be compared and possibly merged later.
  Nothing here derives from, or has been checked against, the `core-redesign` branch.
- **Capability grounding, not interface grounding.** Requirements are derived from what
  `FlightPhysics` and `FlightApps` demonstrably *do* (code, unit tests, demos). Every
  `FlightCore` call site in the consumers is read as evidence of a capability the
  substrate must provide, never as a prescription for how it should be spelled.
- **No interface compatibility.** The new framework need not be source-compatible with
  the current consumers. A non-trivial migration of `FlightPhysics` and `FlightApps` is
  expected and accepted.
- **Guarded additions.** Whenever the design admits functionality beyond what the
  consumers demonstrate, it must be weighed against the fundamental strengths of
  Flight.jl: zero-allocation stepping, type stability, real-time interactive operation,
  live introspection (GUI), and compositional flexibility.

---

## 2. Formalism

The framework simulates **hybrid causal systems**, composed of:

- **Continuous dynamics**: `ẋ = f(x, m, u, t)` with algebraic outputs.
- **Multi-rate periodic discrete dynamics**: `z⁺ = h(z, u, t)` at declared rates, with
  outputs held zero-order between ticks.
- **Zero-crossing events**: guard functions with handlers, in two tiers (below).
- **Post-step manifold projection**: an optional per-component hook `x ← project(x)`
  applied after each accepted step (quaternion renormalization, DCM orthonormalization,
  any manifold-valued state). This is the cheap end of the projection-methods family
  from geometric integration.
- **External inputs**: injected asynchronously by the runtime (pilot controls, network),
  under rules to be fixed in the execution/periphery axes.

### 2.1 Events: two-tier detection

Both tiers share one declaration (guard function + handler); only the detection policy
differs:

- **Tier 1 (default, cheap):** guards are checked for sign changes at step boundaries
  only. No root-finding, no step rejection; the handler fires at the end of the step in
  which the crossing occurred. Cost: one guard evaluation per event per step. Fully
  compatible with fixed-step real-time execution.
- **Tier 2 (opt-in, per event):** localization of the crossing instant by root-finding,
  for events where timing precision genuinely matters. Available in offline/adaptive
  runs; degrades gracefully to Tier 1 in real-time fixed-step mode rather than blowing
  the frame budget.

This gives step-boundary logic *well-defined semantics*: the transition is defined by
the crossing; detection resolution is an execution-policy detail.

### 2.2 Exclusions (deliberate)

- **No DAEs / algebraic constraints.** Projection covers the actual need (state
  manifolds) at near-zero cost and zero solver complexity.
- **No SDEs / stochastic integrators.** Noise processes (Dryden/von Kármán turbulence,
  sensor noise) are modeled as ordinary RNG-driven discrete processes (shaping filters),
  which is both faithful to how they are specified and cheap. Consequence elevated to a
  framework guarantee: **deterministic replay** — RNG state lives in component discrete
  state (`z`), never in ambient globals, so same seed ⇒ bit-identical trajectory.
- **No unconditional per-step hook** (no `f_step!` equivalent). Every current use
  decomposes into projection (quaternion renorm) or Tier-1 events (engine phase
  transitions, stall hysteresis latch). Dropping the hook eliminates the footgun of
  model semantics that depend on the integrator's step size.

---

## 3. Component taxonomy

Three kinds, two of them leaves with crisp, closed semantics, one pure composition:

### 3.1 Continuous component (the hybrid primitive)

A classical hybrid automaton:

- **continuous state** `x` (isbits struct of real scalars — see §8),
- **mode variables** `m`: piecewise-constant values (enums, integers, flags) that
  parametrize the flow and change *only* through event handlers,
- **flow** `ẋ = f(x, m, u, t)`,
- **output groups** (see §5),
- **events**: guards (may read own state, modes, inputs, time) + handlers (update `m`,
  may reset own `x`),
- optional **projection**.

Any facet may be empty. In particular, a component with *no* continuous state — only
modes, events and mode-valued outputs — is an FSM. Both factorings of mode logic are
therefore supported with a single primitive:

- **internal modes** for tightly-coupled cases (stall hysteresis inside the aero
  component) — preserves cohesion and enables reset maps;
- **external FSM component** feeding modes through a port — maximal purity, independent
  testability, swappable supervision logic.

Rule of meaning: a supervisor *commanding* a mode change is an ordinary **input**; a
component *detecting* its own transition is an **event**. Two mechanisms, two meanings,
no overlap.

### 3.2 Periodic discrete component

- **discrete state** `z`: any immutable value (see §9),
- **update** `z⁺ = h(z, u, t)` at a declared rate,
- **output groups** (feedthrough applies at update instants: a proportional path is
  direct feedthrough; a state-only output is not).

`z` influences continuous dynamics only **through signals** (outputs held zero-order
between ticks); no component ever reads another's state.

### 3.3 Assembly

Pure composition: submodels + connections + exported ports. **No dynamics of its own.**
Hybridness emerges at the assembly level (an aircraft = continuous vehicle parts +
discrete avionics parts). Assemblies are flattened away for scheduling but retained as
the navigation/introspection hierarchy (GUI, logging, paths).

---

## 4. Ports and signals

### 4.1 Immutable value semantics

Ports exchange **immutable values** — typically isbits structs (floats, `SVector`s,
enums, nested immutables). The framework owns a **signal table**: one concretely-typed
slot per output port in the flattened model. A producer's output-group function returns
a named tuple of fresh values; the framework writes each into its slot; consumers read
slots.

Consequences:

- no aliasing, ever — nothing can be mutated under a consumer's feet;
- safe concurrent reads (GUI/logging threads) by construction;
- zero allocation for isbits payloads (named tuples of isbits are isbits);
- each slot has a definite freshness tied to its producer's position in the schedule
  (unlike a monolithic `y` struct, which can be half-fresh mid-sweep with no way to
  tell).

The signal requirement, stated precisely, is **immutability plus frozen references**:
signals may reference bulk data (see §7) provided that data is read-only for the
duration of the run. `isbits` is the common case, not the rule.

### 4.2 Consumers see ports, not groups

The port is the addressable unit. A component's outputs appear to consumers, GUI and
logs as one flat namespace (`dyn.vel`, `dyn.accel.f_c_c`, materializable lazily as a
view); which output group computes which port is a scheduling annotation, invisible
outside the component. Regrouping outputs is a non-breaking change.

---

## 5. Evaluation order and feedthrough

### 5.1 The scheduling problem

At every evaluation instant, all signals must be computed consistently: every consumer
reads values already produced at that instant. Build the directed graph of (a) wiring
edges and (b) intra-component feedthrough relations; if acyclic, a topological sort
yields a **static evaluation schedule**, computed once at build time. The hot loop runs
a flat list of `(component, group)` entries — zero runtime graph logic.

### 5.2 Output groups (structural feedthrough)

Feedthrough is declared **structurally, by function signature**:

- Each component partitions its outputs into **output groups**. Each group is a pure
  function receiving *only its declared inputs*; the **stage-1** group receives no
  inputs at all (`g_state(x, m)`), so "no feedthrough" is unfalsifiable — the function
  cannot read what it is not passed.
- **Default factoring is two-stage**: one state-only group, one direct group receiving
  all wired inputs. Finer factorings (more groups) are introduced only when needed
  (§5.3), and are invisible to consumers (§4.2).
- Group functions must be pure (read state/modes/declared inputs, return values, no
  side effects).

The schedule: all stage-1 groups first (any order), then direct groups in topological
order, then all derivative functions against the now-consistent signal table.
Derivative functions, guards, handlers and projections consume signals but impose no
ordering constraints — they run after the sweep.

**Shared expensive computations.** When a derivative function and an output group need
the same expensive result (the canonical case: Newton–Euler must be solved both for
`ẋ` and for the acceleration-like outputs), the framework offers two cache-free
resolutions: **derivative binding** (§9.4 — publish the derivative as an output port
and bind `ẋ` to it; the default idiom for cohesive components) and the
**computer/integrator split** (§9.4 — factor the math into a stateless component whose
outputs feed a trivial state-holding component; the Simulink-diagram idiom, useful when
the factoring earns reuse). Purity rules forbid the third classic resolution, mutable
caching, by design.

### 5.3 Artificial loops and the escape hatch

A component that bundles a no-feedthrough output with a feedthrough output in one
atomic evaluation unit can be **port-level acyclic yet unschedulable** (Simulink's
"artificial algebraic loop"). The canonical instance in this domain is rigid-body
dynamics: velocity out (pure state) + acceleration out (feedthrough from total force).
The two-stage default resolves it. If a component's direct group itself cross-couples
through a neighbor (rare), the author splits it into further groups; the build
diagnostic says so explicitly ("cycle through X.g_direct is artificial at port level —
split the group").

### 5.4 Algebraic loop policy: reject at build time

A genuine cycle in the instantaneous dependency graph is a **build error** with a
diagnostic naming the full path (`aero.F → dyn.a → aero.α̇ → aero.F`). The user breaks
it explicitly: insert dynamics (the α-filter idiom — already standard practice in the
domain and in the current C172 model), insert an explicit unit delay, or restructure.

Rejected alternatives:

- *Implicit unit delays* silently change the model's mathematics at a location the
  framework picked — the archetypal hidden footgun.
- *Numerical loop solving* (Simulink-style fixed-point/Newton per step) has
  data-dependent per-step cost and runtime convergence failures — hostile to real-time
  budgets — and conflicts with immutable signals.
- Implicit *algebraic balances* inside a component (e.g. a turbomachinery operating
  point) remain the component author's business: local, owned, bounded. Rejecting
  framework-level loops does not forbid such models.

---

## 6. Aggregation: reduce-ports

N-to-1 physical aggregation (total wrench, total mass properties, total internal
angular momentum — today's generated `get_wr_b`/`get_mp_b`/`get_hr_b` tree walks) is
expressed by **reduce-ports**:

- An input port may declare a reduction operation (`+`, or any declared
  **commutative-associative** op). Multiple connections to such a port are legal; the
  framework folds the producers' values before the consumer's group runs.
- Every contribution is an explicit wire (`connect(aero.wrench => dyn.wr_Σ_b)`); adding
  a contributor is one line; the scheduler sees one edge per producer.
- **Canonical fold order** (e.g. sorted by producer path): floating-point addition is
  not bit-associative, so the fold order is fixed to preserve bit-reproducibility
  across wiring reorderings.
- **Zero connections is a build error** unless the port explicitly opts into "may be
  unconnected" (identity element). A gear-less vehicle should say so, not silently sum
  nothing.

Rejected alternatives: explicit Σ-junction components (remain trivially expressible as
ordinary components, but positional-slot wiring and arity bookkeeping make them a poor
default); contribution buses (today's mechanism portified — dataflow invisible in the
graph, scoping rules, accidental contributions; reintroduces exactly the implicitness
being removed, at the highest-stakes signal in the model).

---

## 7. Environment access: function-valued signals

Atmosphere and terrain are **query-shaped**: consumers evaluate them at arguments of
their own choosing (each gear strut at its own contact point; airflow at the vehicle
pose). They are therefore carried by ordinary ports as **immutable query objects**
("field handles"):

- An environment component's output group emits a field value (`ISAField(T_sl, p_sl,
  wind)`, `TerrainField(...)`); consumers receive it through ordinary input ports and
  call query functions on it (`airdata(field, pos, vel)`, `ray_intersect(field, p, u)`)
  inside their own group functions.
- **Parametric models are isbits** (ISA, uniform wind, horizontal terrain). **Bulk-data
  models use the handle pattern**: an immutable struct combining isbits parameters with
  references to bulk data (heightmaps, wind grids, the geoid undulation grid) loaded at
  build time and frozen. Handles are rebuilt per evaluation allocation-free (immutable
  structs with existing references — never `Ref`s, whose mutable cell allocates).
- **No mutable caches inside field objects** (memoizing interpolators, lazy loaders):
  concurrent consumers and the GUI thread would race. Caches belong in the consumer's
  state, or the interpolant is restructured to be pure.
- Loggers treat field-handle signals specially (skip or summarize).

Pre-sampling — a component consuming the field and a pose and emitting plain data
(`Airflow` emitting `AirData` for the whole vehicle) — is an **idiom built on top**,
used where natural; not a separate mechanism. Resource injection (declare-and-resolve
service registries) was considered and rejected for the first cut: it is a second
composition mechanism, invisible in the graph — today's argument-threading, automated.

This replaces threading `atmosphere`/`terrain` as arguments through every update
signature, and dovetails with the terrain ray-query direction of the landing-gear
redesign.

---

## 8. Connections and hierarchy

- **Deep connection paths** are allowed, with one structural rule: both endpoints must
  resolve **within the assembly type being defined**. You may deep-route into structure
  you declared (`Cessna172` routing its single `trn` input to
  `systems.ldg.{left,right,nose}.trn_field` in one visible block — no per-level
  re-exports); you may only connect port-level into submodels held **generically**
  (`World` connects `terrain.field => aircraft.trn` and knows nothing more). This kills
  the re-export ceremony where it is ceremony, and preserves substitutability where the
  boundary is load-bearing.
- Paths are validated at build time; renames break loudly.
- Fan-out is free (one producer, many consumers).
- No auto-bubbling of unconnected inputs in the first cut (implicit interface growth).
- Unconnected output ports: build-time warning. Unconnected input ports: build error
  (no silent defaults).

---

## 9. State and data representation

### 9.1 Continuous state: structured immutable, flat backing

Each continuous component declares its state as an **isbits struct whose leaves are all
real scalars of a common eltype `T`** (`SVector`s, quaternions, `Ranged` values —
ultimately reals; `Int`s/enums/`Bool`s belong in modes). The framework:

- computes a **flat layout** at build time (compile-time offsets over one contiguous
  `Vector{T}` buffer it owns);
- **reconstructs** the typed immutable state value for a component at each evaluation
  (field loads at known offsets — register-level, zero cost);
- receives immutable results back: derivative functions return an `Ẋ`-typed value
  (scatter-stored into the flat `ẋ` buffer); event handlers and projection return a new
  `X` (written back).

**The buffer is authoritative; typed values are ephemeral reconstructions.** Nobody
outside the framework ever holds a mutable reference to state.

What this buys, against today's flat-`Vector` + `ComponentArrays`-views pattern:

- no aliased mutable views (who-writes-what by convention);
- derivative completeness is **structural** — the returned `Ẋ` has every field by
  construction; a forgotten `ẋ` entry is impossible rather than silently stale;
- the `SVector{3,Float64}(x.ω_eb_b)` conversion boilerplate disappears — state fields
  already are the types the math wants;
- the flat vector still exists: integrator compatibility (OrdinaryDiffEq or custom),
  trim solvers, HDF5 logging, and linearization all get their arrays;
- the hand-written per-aircraft state-space mapping layer
  (`get_x_ss`/`assign_x_ss!`/`get_u_ss`/...) is deleted, replaced by the framework's
  canonical layout.

### 9.2 Numeric genericity (eltype)

The state buffer, pack/unpack machinery, and the entire **continuous evaluation path**
are generic over `T <: Real`. One design property, four consumers:

1. exact Jacobians for **linearization** (ForwardDiff duals through the whole model,
   replacing finite differences),
2. derivatives for **trim** solvers,
3. the **feedthrough tracer** (§10),
4. a trivially checkable **CI invariant**: one evaluation sweep with `T = Dual` fails
   loudly (`MethodError`/`InexactError` at the offending line) on any Float64-pinning.

Scoping (what actually needs genericity — roughly half the type inventory):

- **Tier 1 — payload/value types constructed during evaluation (~25 structs):**
  the quaternion/attitude family (`Quaternion` becomes
  `Quaternion{N,T} <: AbstractVector{T}` — by invariance, `Float64` instances still
  match every existing `AbstractVector{Float64}` method, so existing behavior is
  untouched), `Wrench`, `FrameTransform`, `MassProperties`, `KinData`, `AirData`,
  geodesy value types, `TerrainData`, continuous output structs. Mechanical
  parametrization; constructors infer `T`, so call sites don't change; `@kwdef`
  defaults pin the no-argument case to `Float64`.
- **Tier 2 — parameters and definitions:** stay `Float64` (promotion handles mixing);
  no migration.
- **Tier 3 — the discrete side** (compensators, avionics): exempt; linearization and
  trim differentiate continuous dynamics only.

Lookups: **table data is a Tier-2 parameter; the query coordinate is Tier-1 traffic.**
Interpolations.jl evaluates generically over the coordinate (`itp(x::Dual)` works
through the `BSpline`/`scale`/`extrapolate` compositions in use). Caveats: `Linear()`
interpolants have kinked derivatives at knots (no regression vs. finite differences;
upgrade to `Cubic` where Jacobian quality near a lookup matters); a manual chain rule
via `Interpolations.gradient` is the escape hatch for anything exotic (and the pattern
for wrapping non-Julia black boxes).

Author-facing rules (three): no `::Float64` argument annotations in math (use `<:Real`
or nothing — the codebase already mostly does); no `Float64`-pinned intermediates
(`zero(SVector{3,T})`); **no `::SomeType{Float64}` return-type annotations** on the
continuous path (they force converts → `InexactError`; see `attitude.jl` `*` for the
live example pattern).

### 9.3 Discrete state, modes, and workspace

- **Storage:** `z` and `m` live in **typed cells** (same discipline as signal slots),
  overwritten by the framework when an update/handler returns a new value. They never
  touch the integrator buffer; no arithmetic is ever done on them.
- **Type freedom:** `z` may be *any immutable value* (frozen-reference rule; isbits not
  required). Enums, integers, nested structs, RNG state (four `UInt64`s of `Xoshiro` —
  required to live in `z` for deterministic replay).
- **Snapshots are free:** copying a cell copies a reference to immutable data —
  checkpoint/replay of the entire discrete side is "copy the cell values."
- **Workspace** (for heavy periodic algorithms, e.g. an n≈20 Kalman filter):
  component-declared mutable scratch, preallocated by the framework, handed to the
  update function, and **excluded from state semantics** — not snapshotted, not
  replayed, must carry no information between calls. Checkable: in debug mode the
  framework **NaN-poisons the workspace before every call**, so read-before-write of
  stale scratch detonates immediately.
- **Blessed idiom — zero-allocation ticks with immutable `z`:** do the in-place math
  (`mul!`, `cholesky!`, BLAS) on the workspace; at the end, snapshot into an isbits
  container (`z = KFState(SVector{20}(x̂_ws), SMatrix{20,20}(P_ws))`) and return it.
  Construction and storage of large `SArray`s are cheap and compile fine — the
  StaticArrays "codegen catastrophe" lives in its *operations* (unrolled matmuls),
  which are never called on snapshots. Discipline: snapshot values are for storage,
  logging, and element access only — never arithmetic. Optionally enforceable by an
  op-forbidding `ValueSnapshot{N,T}` wrapper (an `NTuple` with only
  `getindex`/iteration — structurally what `SArray` is, minus the methods). Practical
  ceiling: a few KB comfortable, tens of KB defensible, beyond that value semantics
  stop making sense.
- **Double-buffered mutable state**: documented as a possible future extension only.
  It reintroduces aliasing/publication questions (GUI reads vs. buffer flips) that the
  first cut deliberately avoids, and nothing in the plausible workload needs it.

### 9.4 Derivative binding

A continuous component may declare that some or all fields of its `ẋ` are **bound to
designated output ports of its own**, instead of being computed by a derivative
function:

```julia
ẋ_bindings(::NewtonEuler) = (ω_eb_b = :ω̇_eb_b, v_eb_b = :v̇_eb_b)   # sketch syntax
```

Semantics and rules:

- Bindings are applied by the framework **after the output sweep**, when the
  component's own outputs are guaranteed fresh (the same freshness argument that
  orders derivative functions after the sweep).
- **Mixing is allowed**: bind some state fields, compute the rest in `f`. A component
  whose fields are all bound needs no `f` at all.
- Build-time validation: every state field is covered by exactly one binding or by
  `f`; every referenced port exists and is type-compatible.

Motivation: derivatives and outputs often share expensive intermediates — rigid-body
dynamics must solve Newton–Euler both for `ẋ` and for the acceleration outputs
(`DynamicsData`). With separate `f` and output groups, the solve would run twice.
Binding publishes the derivative *as a port* (so accelerometer-style consumers can
wire to it) and reduces `f` to a framework-performed copy.

Prior art, for orientation: every causal framework meets this problem and resolves it
per its architecture — Simulink diagrams make integrators explicit blocks (derivatives
are ordinary wires into `1/s`); S-functions and FMUs use sanctioned **mutable caches**
(DWork vectors; FMI's lazy-evaluation caching); Modelica/MTK have `der(x) = expr`
natively with symbolic CSE. Derivative binding is the cache-free formulation that fits
this design's purity rules — Modelica's convenience brought to a causal component API.

The **computer/integrator split** remains fully expressible without any framework
support (a stateless component computing derivatives as outputs, wired into a trivial
state-holding component whose `f` copies them) and is the idiom of choice when the
factoring earns reuse — e.g. one Newton–Euler solver shared across vehicle variants,
or swappable kinematic descriptors against a common integrator shape. See
`sketch.jl` (split form) and `sketch_binding.jl` (merged form with bindings) for the
worked example; the merged form has half the components and wiring, and everything
derivable from pose alone migrates to stage 1, shortening the stage-2 chain.

### 9.5 Allocation policy: a scoped invariant

Not dogma — three reasons, only one about speed: (1) GC-pause jitter control for
real-time; (2) throughput for batch runs; (3) **the canary**: an unexpected allocation
is Julia's most reliable symptom of type instability, so a zero baseline makes
`@allocated == 0` a CI-testable invariant that catches inference regressions at the
offending commit.

- **Continuous hot path** (per-stage evaluation): exactly zero, CI-enforced.
- **Periodic ticks**: zero by idiom (workspace + snapshot pattern); documented
  tolerance for the rare exception.
- **Logging**: amortized-zero — snapshots are isbits records stored *inline* in a
  `Vector` (no boxing, no per-element garbage); `sizehint!` for the expected duration
  makes regrowth a non-event.
- **Tools where garbage is unavoidable**: arena allocation (Bumper.jl-style) for scoped
  temporaries; scheduled `GC.gc(false)` at frame boundaries to move collection out of
  the critical path. (Julia has no per-object freeing; these are the honest levers.)

---

## 10. Diagnostics: feedthrough tracing

Tracing is **diagnostic only, never load-bearing**: scheduling correctness comes
exclusively from structural output groups; tracing improves error messages and
verification. Triggered when the scheduler finds a cycle, to classify it (genuine →
"insert a state"; artificial → "split this group").

Two modes, degrading gracefully:

- **Global (value-blind) set-tracer** — a `Real` subtype carrying a set of input
  indices, unioned by every operation. **May-depend semantics**: sound
  over-approximation (saturated `clamp` still reports; `u²` at `u = 0` still reports).
  Exact, one evaluation. Requires the traced stage-2 code to be free of
  branches/lookups on *input-tainted* values.
- **Local (primal-carrying) set-tracer at sampled states** — the fallback whenever the
  global tracer hits an undecidable branch (piecewise friction, stall blending, any
  gridded lookup at an input-tainted coordinate). Reports the dependence pattern of the
  taken paths, sampled across randomized states. Strictly dominates Dual-based tracing
  (no derivative-zero blind spot; only untaken-branch misses).

Boundaries: only *inputs* are seeded, so branching on state/modes/parameters/time never
interferes. Stage-1 groups are never traced (nothing to seed). Derivatives, guards,
handlers, projections are outside tracing's jurisdiction entirely. Known tracer blind
spot: value-severing operations (dependence passed through a bare `Int` index, e.g.
nearest-neighbor lookup) — documented; linear/cubic interpolation is immune (dependence
flows through the fractional weights).

Both modes ride the same `T <: Real` genericity as `Dual`; Dual-cleanliness in CI
effectively guarantees traceability.

---

## 11. Case study: `Vehicle` today → this framework

The grounding exercise that validated §5. Current `Vehicle.f_ode!`
(`aircraftbase.jl:142-170`) is a hand-woven instance of the machinery specified here:

| Today (convention) | This design (checked structure) |
|---|---|
| `kinematics.u .= dynamics.x` — velocity extracted directly from the state vector because `f_ode!(dynamics)` can't run yet | `dyn`'s stage-1 output group, scheduled first by construction; the artificial loop in `VehicleDynamics` (velocity state-only, accelerations feedthrough) dissolves |
| Hand-ordered `f_ode!` body (kinematics → airdata → systems → route five `dynamics.u` assignments → dynamics last) | Build-time topological sort; wrong wiring = build error naming the cycle or dangling port |
| Velocity state duplicated (`dynamics.x` and `kinematics.u`) with manual sync, incl. `dynamics.x .= kinematics.u  #essential` in `f_init!` | One state, one owner; consumers wire to `dyn.vel` |
| `get_wr_b`/`get_mp_b`/`get_hr_b` generated tree-walk sums | Reduce-ports: one explicit wire per contributor, canonical fold |
| `f_step!` quaternion renorm + engine-phase/stall-latch checks | `@project` hook + Tier-1 events with defined semantics |
| `Aircraft.f_ode!` runs avionics before the vehicle → continuous avionics reads one-stage-stale `vehicle.y` (implicit delay) | Avionics scheduled inside the sweep after the stage-1 outputs it consumes — no delay; or declared periodic and samples post-step by stated semantics |
| `atmosphere`/`terrain` threaded as arguments through every signature | Field-handle signals through ordinary ports (§7) |

The migration cost surfaced by the same exercise: today's monolithic `KinData` splits
into `pose` (stage 1: `q_eb`, `r_eb_e`, `ϕ_λ_h`, ...) and `kin_vel` (stage 2:
`v_eb_n`, `v_gnd`, `χ`, `γ`, ...) because the parts genuinely have different
dependencies. The recurring trade, stated once: the framework asks authors to write
down structure previously kept in their heads, and pays them back by never letting it
silently rot.

The genuine algebraic loop in the domain — α̇-dependent aerodynamics — is already
broken in the current C172 model by a filter state, exactly the explicit break §5.4
prescribes. Evidence that reject-loops matches domain practice rather than fighting it.

---

## 12. Decision log (condensed)

| # | Decision | Rejected alternatives (why) |
|---|---|---|
| 1 | Hybrid causal formalism; two-tier events; projection; no DAE/SDE/per-step hook | DAEs (projection suffices), SDEs (shaping filters suffice), `f_step!` (step-size-dependent semantics) |
| 2 | Causal port-based paradigm | Acausal/MTK (fights interactivity, discrete logic, live introspection); hierarchical callables (rigor by convention — today's footguns); thin SciML library (nothing for GUI/logging/hierarchy to hang onto) |
| 3 | Taxonomy: hybrid continuous primitive + periodic discrete + assemblies; both mode factorings | Strict purity/no modes (loses reset maps; latch logic becomes wiring ceremony); uniform hybrid kind (intra-component ordering semantics murky) |
| 4 | Immutable value signals in a typed signal table | Shared mutable buffers (aliasing/staleness, concurrent-read hazards); mixed semantics (second thing to document/test) |
| 5 | Reject algebraic loops at build, explicit breaks | Implicit delays (silent math changes); per-step numerical solving (jitter, runtime failures) |
| 6 | Structural output groups, two-stage default | Per-output declarations + tracer verification (under-declaration = silent wrongness); traced multi-pass (hot-path cost, branch unsoundness); component-atomic conservative (false loops everywhere) |
| 7 | Reduce-ports with canonical fold | Σ-junctions as default (arity/positional ceremony); contribution buses (invisible dataflow) |
| 8 | Function-valued environment signals + handle pattern | Resource injection (second composition mechanism, invisible); pre-sampling as mechanism (dependency inversion at struts) |
| 9 | Deep paths within owned assembly types only | Unrestricted deep wiring (breaks substitutability across generic boundaries); strict one-level (re-export ceremony) |
| 10 | Structured immutable state over framework-owned flat `Vector{T}` | Mutable views (aliasing, silent missing-ẋ); fully structured/no flat vector (same machinery needed anyway, loses standard integrator interface) |
| 11 | Eltype genericity on the continuous path, three-tier scoping | Float64-only + finite differences (FD noise, keeps hand-written state-space layer, no tracer) |
| 12 | Set-propagation tracers (global + sampled-local), diagnostic-only | Dual-based tracing (derivative-zero blind spot); tracing as scheduling input (soundness) |
| 13 | Immutable `z` in cells + workspace + snapshot idiom | Mutable discrete state (aliasing, snapshot cost); double-buffering (deferred; publication races) |
| 14 | Scoped allocation invariant, CI-enforced on the hot path | Blanket dogma (fights logging reality); no policy (loses the type-instability canary) |
| 15 | Derivative binding (ẋ fields bound to own output ports, applied post-sweep); computer/integrator split as reuse idiom | Mutable caches between `f` and outputs (S-function/FMI style — hidden state, purity violation); accepting duplicate computation (Newton–Euler solved twice) |

---

## 13. Open axes

To be settled in subsequent sessions:

- **Time & execution (axis 5) — next up.** Integrator ownership (build on
  OrdinaryDiffEq vs. purpose-built fixed/adaptive stepping core); multi-rate tick
  scheduling against integrator steps; event-iteration semantics at step boundaries
  (single pass vs. iterate-to-quiescence, bounded — externalized FSMs make cascades
  more common); real-time pacing.
- **Runtime periphery (axis 6).** GUI, input devices, network I/O, logging; the
  concurrency model binding them to the sim loop; when the runtime may write
  externally-driven input slots (staging vs. step boundaries) to preserve determinism.
- **User-facing surface (axis 7).** Declaration layer (macros vs. plain
  constructors — all syntax in this document is illustrative sketch, not committed);
  build/verification tooling and error-message quality; initialization, trim and
  linearization workflows as first-class framework services.
- **Migration.** Outline for FlightPhysics/FlightApps (including the Tier-1
  parametrization pass and the `KinData`-style output splits); comparison criteria
  against FlightCore's demonstrated strengths (zero-alloc stepping, flexibility,
  interactive operation).
