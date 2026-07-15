# A Modeling & Simulation Framework for Flight.jl — Design Document

**Status:** second checkpoint (v0.2). Axes 1–5 settled; axes 6–7 and the migration
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
  for events where timing precision genuinely matters (mechanics in §11.4). Available
  in offline runs regardless of stepping mode; degrades gracefully to Tier 1 in
  real-time mode rather than blowing the frame budget.

This gives step-boundary logic *well-defined semantics*: the transition is defined by
the crossing; detection resolution is an execution-policy detail.

Guards may be **boolean predicates** (checked for becoming true at step boundaries) or
continuous sign-crossing functions. Tier 1 handles both; Tier 2 localization requires
the continuous form. This matters in practice: most transitions in FlightPhysics mix
input predicates with state thresholds (e.g. the piston engine's `starting → running`
fires on `ω > ω_idle && fuel_available`).

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

- **continuous state** `x` (isbits struct of real scalars — see §9),
- **mode variables** `m`: piecewise-constant values (enums, integers, flags) that
  parametrize the flow and change *only* through event handlers,
- **flow** `ẋ = f(x, m, u, t)`,
- **two output stages** (see §5.2),
- **events**: guards + handlers (update `m`, may reset own `x`); both read the fresh
  boundary signal table (§5.2),
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
- **two output stages** (feedthrough applies at update instants: a proportional path
  is direct feedthrough; a state-only output is not).

`z` influences continuous dynamics only **through signals** (outputs held zero-order
between ticks); no component ever reads another's state.

### 3.3 Assembly

Pure composition: submodels + connections + exported ports. **No dynamics of its own.**
Hybridness emerges at the assembly level (an aircraft = continuous vehicle parts +
discrete avionics parts). Assemblies are flattened away for scheduling but retained as
the navigation/introspection hierarchy (GUI, logging, paths) and as declaration-level
rate scopes (§11.5).

---

## 4. Ports and signals

### 4.1 Immutable value semantics

Ports exchange **immutable values** — typically isbits structs (floats, `SVector`s,
enums, nested immutables). The framework owns a **signal table**: one concretely-typed
slot per output port in the flattened model. A producer's output-stage function returns
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

### 4.2 Consumers see ports, not stages

The port is the addressable unit. A component's outputs appear to consumers, GUI and
logs as one flat namespace (`dyn.vel`, `dyn.f_c_c`, materializable lazily as a view);
which output stage computes which port is a scheduling annotation, invisible outside
the component. Moving an output between stages is a non-breaking change.

**Unlisted ports.** A port may be marked *unlisted* (surface syntax TBD, axis 7):
skipped by default in logging and GUI panels, and exempt from unconnected-output
warnings. Purely presentational — unlisted ports remain connectable and inspectable on
request. Intended for state published only to feed the component's own dynamics (§5.2),
e.g. RNG state.

---

## 5. Evaluation order and feedthrough

### 5.1 The scheduling problem

At every evaluation instant, all signals must be computed consistently: every consumer
reads values already produced at that instant. Build the directed graph of (a) wiring
edges and (b) intra-component feedthrough relations; if acyclic, a topological sort
yields a **static evaluation schedule**, computed once at build time. The hot loop runs
a flat list of `(component, stage)` entries — zero runtime graph logic.

### 5.2 Two-stage outputs and the state decoder (structural feedthrough)

Every component provides exactly **two output stages**, and feedthrough is declared
**structurally, by function signature** — there are no dependency annotations anywhere
in the design:

```julia
# continuous component                      # discrete component
y_s1 = g_s1(comp, x, m)                     y_s1 = g_s1(comp, z)
y_s2 = g_s2(comp, u, y_s1)                  y_s2 = g_s2(comp, u, y_s1)
ẋ    = f(comp, y, u, t)                     z⁺   = h(comp, y, u, t)

# event system (continuous side only) — readers of the same fresh table:
guard(comp, y, u, t)
(x⁺, m⁺) = handler(comp, y, u, t)           # constructs successors from y; may reset x

# the sole raw-state function:
x⁺ = project(comp, x)                        # manifold projection
```

- **`g_s1` is the state decoder** — the *only* place where raw state becomes
  information. It receives no inputs, so "no feedthrough" is unfalsifiable: the
  function cannot read what it is not passed. **Default `g_s1` is identity
  publication**: absent an explicit method, state fields *and mode fields* (`x`, `m`
  for continuous components; `z` for discrete) are published as ports named after
  them. The commitment this creates: `y_s1` must carry whatever the component's own
  downstream functions need — satisfied automatically by the default, and violated
  *loudly* (missing-field error on first evaluation) otherwise.
- **`g_s2` receives all wired inputs plus `y_s1`** — its own stage-1 results, so
  state-derived intermediates are computed once, not re-derived. It has no direct
  state access; conservatively, every stage-2 output is presumed dependent on every
  wired input.
- **`f` and `h` are pure signal consumers.** They run after the sweep, when the full
  signal table — including the component's own `y_s2` — is complete and fresh. `ẋ` and
  `z⁺` are computed from the same published values every other consumer sees:
  single-source-of-truth is mechanical, and the duplicated-math drift class (edit the
  torque law in `f`, forget the copy in `g`) cannot exist.
- **Guards and handlers are readers too.** At a step boundary, the order is
  *integrate → project → boundary sweep → guards*, so by guard/handler time `y` is a
  fresh decode of exactly the state being transformed; handlers construct `(x⁺, m⁺)`
  from it under the same commitment as `h`. **`project` is the sole raw-state
  function**, structurally: it runs *between* a state write and its decode (after
  integration, and after any handler `x`-reset), i.e. at the only positions in the
  schedule where no fresh `y` of the new state can exist yet.
- All of `g_s1`, `g_s2` must be pure (no side effects); state types make mutation
  impossible anyway (§9).

The schedule: all `g_s1` (any order), then `g_s2` in topological order, then all `f`
against the now-consistent signal table. Note the systemic consequence: *evaluating
the RHS means running the sweep* — there is no incremental `f`-only re-evaluation.
Implicit solvers, linearization and trim already work this way (seed `x`, run the
composite), so nothing is lost; axis 5 should restate it as a property of the
execution model.

**Step-boundary semantics.** At each boundary: integrate → project → boundary sweep →
evaluate **all guards once** against that sweep → for each fired event, in declaration
order: `handler → project → re-run the component's g_s1 and g_s2`. The per-event
re-decode keeps `y` fully fresh for any subsequent handler of the *same* component
(sequential composition, no lost updates) and leaves the signal table
post-transition-consistent for whatever else the boundary does (discrete ticks,
logging). Events are rare and the re-run touches one component, so the cost is noise.
Across components, all guards and handlers read the same boundary snapshot (plus their
own component's refreshed ports); one component's transition reaches others through
the next sweep. Newly-enabled guards fire within the *same* boundary: the
sweep → guards → handlers phase iterates to quiescence, with each event firing at
most once per boundary (settled in §11.6).

**Departure from the orthodox formalism, stated openly.** The textbook form is
`ẋ = f(x, u)`, `y = g(x, u)`; this design's is `y = g(x, u)`, `ẋ = f(y, u)`. The
composite map `x ↦ ẋ` is mathematically identical (linearization, trim and AD are
untouched), but the spelling is heterodox and every literate newcomer will pause at
it once. The teaching line: *"stage 1 publishes what you know from state alone — by
default, the state itself; stage 2 adds what needs inputs; your dynamics then read the
same table everyone else does."* The decision was grounded in a component-by-component
survey of FlightPhysics/FlightApps (§12.2): derivative/output overlap is the *norm*
in this domain (Newton–Euler, kinematics, piston engine, gear friction, every discrete
compensator), so the orthodox split would force either systematic duplicated math (with
its silent-drift bug class), systematic component atomization, or a cache mechanism —
all rejected. FlightCore's fused `f_ode!` already embodied the same economics; this
design keeps them while adding checked scheduling.

**Shared expensive computations** are thereby solved uniformly: compute once in
`g_s2`, publish, and let `f`/`h` (and external consumers, e.g. an accelerometer model
reading `f_c_c`) consume the ports. The **computer/integrator split** — a stateless
component computing derivatives as outputs, wired into a trivial state-holding
component — remains fully expressible without framework support and is the idiom of
choice when the factoring earns reuse (one Newton–Euler solver shared across vehicle
variants; swappable kinematic descriptors against a common integrator shape). Purity
rules forbid the third classic resolution, mutable caching, by design.

### 5.3 Artificial loops and the escape hatch

A component that bundles a no-feedthrough output with a feedthrough output in one
atomic evaluation unit can be **port-level acyclic yet unschedulable** (Simulink's
"artificial algebraic loop"). The canonical instance in this domain is rigid-body
dynamics: velocity out (pure state) + acceleration out (feedthrough from total force).
The two-stage split resolves it. In the rare case where a single component's stage-2
outputs cross-couple through a neighbor (port-level acyclic, stage-level cyclic), the
remedy is **splitting the component** — which documents real structure — and the build
diagnostic says so explicitly ("cycle through X.g_s2 is artificial at port level —
split the component"). One consequence of stage-2 conservatism worth recording: an
input consumed only by `f` (never by `g_s2`) still creates a scheduling edge if the
component has stage-2 outputs; in practice such components are integrator-shaped and
have none, and the remedy, if ever needed, is the same split.

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

- An environment component emits a field value (`ISAField(T_sl, p_sl,
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
  and passes it to the state readers/writers (`g_s1`, handlers, projection — field
  loads at known offsets, register-level, zero cost);
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

### 9.4 The fused-evaluation lineage (prior art and how we got here)

The §5.2 interfaces are the end point of a three-step simplification arc, recorded
here because each step replaced a mechanism with something smaller:

1. **N output groups → exactly two.** General output groups (each declaring its input
   subset) handled a cross-coupling case that never materialized in the domain; strict
   two-stage eliminates all dependency declarations, at the price of an occasional
   component split (§5.3).
2. **Derivative binding → own-output access.** An earlier revision had a declaration
   feature binding `ẋ` fields to output ports (`ẋ_bindings(::C) = (ω = :ω̇,)`). Passing
   the fresh signal table to `f`/`h` subsumes it: the "binding" is a one-line function
   body, no validation machinery, strictly more general (an `f` may combine published
   values with extra terms).
3. **Separate state arguments → the `g_s1` decoder.** With `y` in hand, passing `x` to
   `f` duplicated information (FlightPhysics culturally publishes state in `y` anyway);
   removing it produced the uniform continuous/discrete shapes and the
   single-computation-site guarantee.

Prior art, for orientation — every causal framework meets the shared-computation
problem and resolves it per its architecture: **Simulink diagrams** make integrators
explicit blocks (derivatives are ordinary wires into `1/s` — the computer/integrator
split is their native idiom); **S-functions and FMUs** use sanctioned *mutable caches*
(DWork vectors; FMI's lazy-evaluation caching) between their `mdlDerivatives`/
`mdlOutputs`-style callback pairs; **Modelica/MTK** write `der(x) = expr` natively
with symbolic CSE. The fused sweep + signal-consuming `f`/`h` is the cache-free
formulation that fits this design's purity rules — and it is also what FlightCore's
fused `f_ode!` did economically, minus the checked scheduling.

The **computer/integrator split** remains fully expressible without any framework
support (a stateless component computing derivatives as outputs, wired into a trivial
state-holding component) and is the idiom of choice when the factoring earns reuse —
e.g. one Newton–Euler solver shared across vehicle variants, or swappable kinematic
descriptors against a common integrator shape. See `sketch.jl` (split form) and
`sketch_decoder.jl` (merged form under the §5.2 interfaces) for the worked example;
the merged form has half the components and wiring, and everything derivable from
pose alone migrates to stage 1, shortening the stage-2 chain.

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
exclusively from the structural two-stage split; tracing improves error messages and
verification. Triggered when the scheduler finds a cycle, to classify it (genuine →
"insert a state"; artificial → "split this component").

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
interferes (under the §5.2 decoder, state-derived values arrive through `y_s1`, which
is never seeded — same conclusion). Stage-1 functions are never traced (nothing to
seed). Derivatives, guards, handlers, projections are outside tracing's jurisdiction
entirely. Known tracer blind
spot: value-severing operations (dependence passed through a bare `Int` index, e.g.
nearest-neighbor lookup) — documented; linear/cubic interpolation is immune (dependence
flows through the fractional weights).

Both modes ride the same `T <: Real` genericity as `Dual`; Dual-cleanliness in CI
effectively guarantees traceability.

---

## 11. Time and execution (axis 5)

### 11.1 Loop ownership: the framework owns the simulation loop

The simulation loop — the §5.2 boundary sequence, tick dispatch, event handling,
logging, input staging, pacing — is **framework code, unconditionally**. The step-
boundary contract is this design's central invariant; expressing it as an ordered
`CallbackSet` inside a third-party solver would put the framework's semantics back
into convention territory, enforced by callback-registration order in a foreign event
loop. That is the same "rigor by convention" failure mode the redesign exists to
eliminate.

The evidence that this is not hypothetical comes from FlightCore's own `sim.jl`: the
periodic callback is hand-rolled with explicit `add_tstop!` bookkeeping because a
DiffEqCallbacks release moved `PeriodicCallback` onto `task_local_storage`, breaking
the init-on-one-task/step-from-another pattern that interactive operation requires;
the RHS wrapper copies state in and out of the integrator per evaluation; stateless
models integrate a dummy `[0.0]`; log saving detours through `deepcopy` in a
`SavingCallback`. Each is a tax paid to express our semantics in someone else's loop —
and `run!` already drives the integrator interface manually anyway.

`OrdinaryDiffEq` is therefore **dropped as a dependency** of the new core.

### 11.2 The stepper seam

The one delegated operation is *advance the continuous state from `t` by `h`*, behind
a narrow internal interface (the **stepper seam**). Its contract:

- **advance by arbitrary `h`** — required anyway to land on tick boundaries and to
  resume from a localized event time;
- **dense output on demand over the last completed step** — required only by Tier-2
  localization (§11.4), constructed lazily;
- **one-step methods only** — event handlers reset state discontinuously, and a
  one-step method restarts from a new state for free; multistep methods would need
  history-rebuild machinery after every handler and are excluded.

First cut ships **in-house fixed-step RK4 and Heun** over the flat state buffer:
~a hundred lines, trivially zero-allocation (auditable for the §9.5 CI invariant),
trivially `T`-generic — though genericity is not even required of the stepper, since
linearization and the tracer drive the *sweep*, never the integrator. An
`OrdinaryDiffEq`-backed stepper can exist later as a package extension if an offline
study genuinely demands adaptive or stiff methods; per the guarded-additions rule it
is not built until then.

**Why fixed-step low-order suffices** (the domain argument, recorded because it is
load-bearing for the whole axis):

1. **The closed-loop tick cap.** Every application beyond bare propagation runs
   periodic avionics (50 Hz today), whose commands are zero-order-held signals:
   integrating past a tick with stale commands is wrong, so the integrator must land
   on every tick boundary regardless of method. Adaptive and high-order methods pay
   off exactly when steps can stretch; the execution model forbids the stretch by
   construction.
2. **A piecewise-smooth RHS starves high order.** Linearly-interpolated lookup tables
   (C¹-kinked at every knot), clamps, friction blends and mode branches deny
   high-order error estimators and implicit-solver Newton iterations the smoothness
   they assume. RK4 at 50 Hz already puts integration error orders of magnitude below
   the model uncertainty of a coefficient-table aircraft model.
3. **Stiffness has a remedy ladder.** The fastest continuous dynamics in the current
   codebase (actuator poles ~31 rad/s, gear damper decay, friction compensators) sit
   inside RK4's stability region at `h = 0.02` — the crosswind-landing demo is the
   empirical proof. If a future model exceeds it: first shrink `h` (the RHS costs
   microseconds; 500 Hz real-time is unremarkable), then subcycle the stepper against
   the tick grid, and only then reach for an implicit method through the adapter.
   If that day comes, §9.2 genericity supplies exact ForwardDiff Jacobians through
   the sweep for free.

### 11.3 Signal-table consistency is a boundary property

During a step, RK stages evaluate the sweep at internal stage states — the signal
table is transiently **integrator scratch**. The boundary sweep in the §5.2 sequence
is what restores consistency at each accepted boundary. The rule, binding for axis 6:
**external readers (GUI, logging, network output) observe the signal table only at
step boundaries.** Mid-step contents carry no meaning.

### 11.4 Tier-2 localization mechanics

Trigger: a Tier-2 guard changed sign across an accepted step `[tₙ, tₙ₊₁]`.

- **Interpolant (lazy).** Build the cubic Hermite continuous extension `x̂(θ)`,
  `θ = (t − tₙ)/h ∈ [0,1]`, from `(xₙ, ẋₙ, xₙ₊₁, ẋₙ₊₁)`; `ẋₙ` is the step's first
  stage, `ẋₙ₊₁` costs one sweep, paid only on trigger. Uniform accuracy O(h⁴) — one
  order below the discrete solution, the standard pairing, and the event time can
  only ever be as accurate as the interpolant, so nothing more expensive is worth
  probing.
- **Probes run the sweep.** Guards read `y`, so evaluating a guard at an interpolated
  state means writing `x̂(θ)` into the state buffer and running the sweep — the same
  rule as the RHS. One sweep per probe.
- **Root-finding: bracketed and derivative-free** (ITP or Brent; bisection is an
  acceptable fallback). The observed sign change *is* an unconditional convergence
  certificate. Newton/AD localization is rejected: guards are guaranteed C⁰, not C¹
  (clamps, table knots, saturated stretches where σ′ = 0), Newton discards the
  bracket for merely local guarantees, and its superlinear convergence saves a
  handful of microsecond probes per rare event. AD earns its keep in Jacobians, not
  in root-polishing a possibly-kinked bracketed scalar.
- **Post-event.** Handler at `t*` → project → re-decode (per §5.2) → **interpolant
  invalidated** (the handler made it a lie for `t > t*`) → resume integration from
  `t*` with the remainder step `h′ = tₙ₊₁ − t*` → re-check guards on the remainder,
  under a bounded per-step event budget with a chattering diagnostic.
- **Shared blind spot, documented:** an even number of crossings within one step
  defeats sign-change detection in both tiers; the mitigation is step size, not
  machinery.

### 11.5 Multi-rate tick scheduling

**Harmonic grid.** Every discrete component's period is an integer multiple of a base
tick period `Δt_base`, itself an integer multiple of the continuous step
(`Δt_base = n·h`, `n ≥ 1`). Ticks therefore land on step boundaries — the only place
anything discrete ever happens. Rejected: arbitrary periods via a time-ordered tick
queue, which forces variable-length steps and irregular real-time frames for a
generality nothing demonstrated wants.

**Discrete stages are gated to tick instants.** A discrete component's `g_s1`/`g_s2`
run only at its own ticks; its slots hold in between (ZOH, stated in sweep terms). The
alternative — re-running its stages at every boundary — would let outputs drift
between ticks as fresh continuous inputs flow in, silently un-sampling a sampled-data
controller. Implementation consequence, recorded: the boundary sweep is not one fixed
list — discrete stage entries are gated by tick counters, so different boundaries run
different subsets of the schedule.

**Simultaneous ticks are already well-defined** by settled machinery: all due
components run their `g` stages in topological order within the sweep, then all due
`h` updates run after it, in any order — each `h` reads the table and writes only its
own `z` cell. The FCS cascade's intra-tick ordering is a sweep property, not an
update-order property.

**Assemblies: virtual for execution, rate scopes for declaration.** There are no
atomic assemblies, and no opt-in variant. Execution atomicity coarsens the schedulable
unit, which is exactly how artificial algebraic loops are manufactured — §5.3 at
assembly scale, a hazard Simulink documents for its Atomic Subsystems — while the
thing it protects, non-interleaved execution, protects nothing here: the signal table
makes interleaving semantically invisible (consumers read slots whose freshness is
guaranteed by topological order, not contiguity). Its other Simulink roles — code
generation units, enabled/triggered execution — have no counterpart in the consumers;
conditional behavior ("on ground, force `direct`") is mode logic. FlightCore's
whole-tree atomicity (children's `f_periodic!` called from the parent's, hence
`Subsampled`'s parent-relative multipliers) was an artifact of call-tree execution,
not a design requirement; the signal table dissolves it.

**Rate declaration is relative, at composition.** A discrete component or sub-assembly
is instantiated with an integer multiplier `K ≥ 1` (default 1) relative to its
enclosing scope; multipliers compose multiplicatively down the tree; the root fixes
`Δt_base` in seconds. Rationale: in a layered control architecture the *ratios* are
intrinsic to the design and travel with the assembly type (inner loop at `K = 1`,
outer loops at `K = 5`, whatever the deployment), while absolute rates are deployment
decisions made once at the root. Absolute-first declaration was rejected: it welds
deployment rates into reusable definitions, and replicating relative structure with a
shared base-period variable does not compose across independently authored assemblies
(parameter-threading ceremony, cf. §7). At build, all scoping compiles away to **one
absolute tick divisor per discrete component**; the executor gates entries by counter
modulo. Recorded limitations: a child cannot tick faster than its scope (`K ≥ 1`) —
soft, since assembly multipliers are declaration sugar and factors can be refactored
onto siblings; and no phase offsets in the first cut (no demonstrated use).

**`Δt` has a single source of truth: the compiled schedule.** Each discrete
component's effective period is exposed read-only through the component handle
(`comp.Δt` — the same virtual-property move as FlightCore's `mdl.Δt`), available in
`g_s1`/`g_s2`/`h` and an error to touch on a continuous component. It must be readable
in the *stages*, not just `h`: per §12.2, the discretized laws that actually consume
`Δt` — a PID's backward-difference coefficients, a LeadLag's Tustin transform — run in
`g_s2`; `h` is a copy. Author rule: **never store `Δt`, or any `Δt`-derived
coefficient, as a component parameter** — recomputing derived coefficients per tick is
a few arithmetic ops, and a cached copy is a second thing for gain-scheduled
`assign!`-style updates to chase. Relative declaration structurally enforces the rule
for the period itself: under scoped multipliers a component author *cannot* know their
absolute rate — it does not exist until composition.

### 11.6 Event iteration at boundaries: to quiescence, once per event

Resolves the question deferred in §5.2. At each boundary, the event phase **iterates**:
rounds of *(re-run the boundary sweep → evaluate guards → fire newly-fired events in
declaration order, each with its per-event re-decode)* until a round fires nothing,
under the rule that **each declared event fires at most once per boundary**.

**Why iterate.** Under a single pass, a cascade of N logically-simultaneous
transitions (supervisor FSM → subordinate FSM → …) completes in N steps: latency N·h,
with h an execution parameter. That is model semantics depending on the integrator's
step size — the same footgun class §2.2 cited when killing `f_step!` — and §3.1's
blessing of externalized FSM components makes cross-component cascades the expected
idiom, not a corner case. Orthodoxy concurs: hybrid automata take sequences of
instantaneous transitions at one time point; Modelica iterates events to quiescence;
Stateflow runs charts to completion within a tick. (Tier-1 detection timing remains
h-dependent — that is the *resolution* at which a physical crossing is noticed; the
cascade delay would have been structure the framework inserts between transitions the
model declares simultaneous.)

**Why a full re-sweep per round.** The per-event re-decode refreshes only the
transitioning component's own ports; downstream stage-2 chains reading them stay
stale. A round therefore re-runs the whole gated schedule. Sweeps are microseconds
and rounds beyond the first require an actual cascade, so the cost is noise.

**Why once-per-event rather than a round cap.** Termination becomes structural —
rounds are bounded by the number of declared events, with no arbitrary K knob — and a
livelock (two FSMs toggling each other) resolves deterministically into "each fired
once, quiescence, re-fire next boundary," i.e. degrades to Tier-1 granularity instead
of burning a budget and erroring. The cost: an event legitimately re-enabled within
the same boundary waits one step — accepted at the same granularity Tier 1 accepts
physical re-crossings, and flagged by a diagnostic when it occurs.

**Ticks stay outside the iteration, after quiescence.** The two possible couplings
resolve asymmetrically:

- *Events → ticks: handled by machinery already in place.* Due discrete components'
  `g` stages are gated *into* the sweep (§11.5), so every iteration round refreshes
  them for free — same `z` (their `h` has not run), post-transition inputs. At
  quiescence, their published outputs reflect the settled boundary instant, which is
  exactly what "sampling at t" should mean for a logically-instantaneous cascade.
  Earlier rounds' tentative values are internal scratch, like RK stage evaluations;
  §11.3 extends naturally: external readers observe the table only after the boundary
  sequence *completes*.
- *Ticks → events: structurally impossible.* A tick's `g` stages contribute nothing
  guards have not already seen (they run inside the sweep, from current `z`); its `h`
  update writes `z⁺`, which is invisible until the component's *next* tick decodes it
  (`g_s1(comp, z)` is the sole reader of `z`) — the standard one-sample `z⁻¹` delay of
  sampled-data control, here enforced by construction. Nothing that happens after
  quiescence can flip a guard, so no combined event/tick fixed point exists to
  iterate.

The boundary macro-sequence, final form:

> integrate → project → **[sweep → guards → handlers]** iterated to quiescence
> (once per event) → all due `h` updates → logging / I/O staging.

The mixed case — a continuous component's handler and its discrete observers' ticks
landing on one boundary (engine `starting → running` under a 50 Hz FCS) — is decided
by the sequence: the transition fires in the iteration segment, the re-sweep re-runs
the FCS's stages against `running`-mode ports, and its `h` then updates from
post-transition values.

### 11.7 Real-time pacing

**The invariant: pacing is outside the semantics.** The pacer inserts waits between
completed boundaries and never reorders, skips or alters the boundary sequence. A
paced and an unpaced run with identical input traces produce bit-identical
trajectories — deterministic replay (§2.2) extends over pace. Interactive runs differ
only because their *inputs* differ.

**Wall-clock mapping: piecewise affine, re-anchored at every knee.** The map is
`τ(t) = τ_anchor + (t − t_anchor)/p`, with the anchor pair as its reference point. A
live pace change re-establishes the anchor at the current `(t, τ)` so the new slope
applies only forward; keeping the old anchor would retroactively reinterpret the
entire elapsed run at the new pace (deadlines minutes in the past or future after a
long session). Un-pause re-anchors for the same reason. Debt is cleared at re-anchor:
a deliberate user action is a natural sync point, and the counters record what was
forgiven.

**Deadline law: absolute schedule with bounded debt.** Boundary deadlines come from
the map; a frame exceeding its wall budget `h/p` leaves debt that subsequent frames
repay by running short or waitless — the long-run rate is exact and ms-scale hiccups
(GC, scheduler) are invisible. Debt beyond a threshold (a few frames' worth) is
forgiven by re-anchor plus warning, so long stalls (debugger, laptop sleep) do not
trigger catch-up bursts. Rejected: relative deadlines (next = last completion +
budget), under which every overrun permanently slips sim time against wall time.

**`p = ∞` is pacer-off, not a limit value.** FlightCore's arithmetic trick
(`τ_next = τ_last + dt/∞` collapses the wait predicate) does not survive debt
accounting: every deadline would sit perpetually in the past and the diagnostics
would faithfully report garbage. Unpaced mode is the explicit *absence* of deadlines —
no waits, no debt, no warnings; by the invariant, the same execution with the waits
deleted.

**Wait mechanism: hybrid sleep-then-spin, one knob.** Non-realtime OSes guarantee
only a lower bound on sleep: the thread becomes runnable no earlier than requested;
the wake-up is best-effort (timer granularity, scheduler load, macOS timer
coalescing), with no hard upper bound. Measured on the dev machine (idle, 2 ms
requests, 2026-07): Julia `sleep` overshoots ≈ 1.4 ms median — libuv's
millisecond-granularity timers; sub-ms requests are accepted and rounded up —
and `Libc.systemsleep` ≈ 0.5 ms; spikes under load are unbounded. The pacer therefore
sleeps toward `deadline − margin` and spins the remainder:

```julia
remaining = deadline - margin - τ()
remaining > 0 && sleep(remaining)   # coarse phase: cheap, lower-bound-only (runs at most once)
while τ() < deadline end            # spin phase: µs-precise, CPU cost bounded by margin
```

`margin` is a single constant calibrated to cover the primitive's granularity *plus*
typical overshoot (no second threshold — the resolution floor is absorbed into the
calibration, and a margin below the primitive's granularity defeats the spin phase's
purpose). It spans the whole design space:

- **`margin = 0` — pure sleep:** cheapest CPU; bursty boundary spacing, but the
  absolute schedule still delivers the exact *average* rate through debt repayment.
  The spin phase buys regularity, never rate correctness.
- **`margin` = a few ms — hybrid (default):** sleeps ~90% of a 20 ms budget, lands
  within µs of the deadline at a few percent of one core.
- **`margin = ∞` — pure busy-wait:** FlightCore's behavior, maximum boundary
  regularity at one pinned core; the "best attempt at real time" mode is the knob's
  endpoint, not a separate mechanism.

When the frame budget is at or below the margin (e.g. `h = 0.01` at `p = 5` → 2 ms),
the hybrid degenerates to pure spin per frame by construction. Rare wake-ups past the
deadline are overruns, absorbed as debt. Which primitive the coarse phase uses —
task-yielding `sleep` vs. thread-blocking `Libc.systemsleep` — is an axis-6
concurrency decision.

**Diagnostics.** Overrun count, current and peak debt, forgiven-debt events, wait
statistics — published as framework status for GUI and logs (today's `SimControl`
fields are the precedent).

**Forward pointers.** The wait interval is the natural staging slot for externally
injected inputs, applied at the next boundary; the staging rules — and the
concurrency model generally, which §11.3 constrains but does not decide — belong to
axis 6.

---

## 12. Case studies

### 12.1 `Vehicle` today → this framework

The grounding exercise that validated §5. Current `Vehicle.f_ode!`
(`aircraftbase.jl:142-170`) is a hand-woven instance of the machinery specified here:

| Today (convention) | This design (checked structure) |
|---|---|
| `kinematics.u .= dynamics.x` — velocity extracted directly from the state vector because `f_ode!(dynamics)` can't run yet | `dyn`'s stage-1 output, scheduled first by construction; the artificial loop in `VehicleDynamics` (velocity state-only, accelerations feedthrough) dissolves |
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

### 12.2 Torture tests for the §5.2 interfaces: `PistonEngine` and the FCS PID cascade

Two components were transliterated to validate the decoder interfaces before adoption.

**`PistonEngine`** (piston.jl:310-449) — mode enum with three flow regimes, four table
lookups, two embedded continuous PI compensators, boolean transitions, argument-threaded
`fuel_available`:

- The compensator paths (`idle`, `frc`) are pure functions of the engine's own state
  (`ω`), so their complete PI laws — outputs *and* state derivatives — evaluate in
  `g_s1`. (The alternative factoring, compensators as child components of an engine
  assembly, also schedules cleanly from the core's stage-1 ports.)
- `g_s2` runs the lookup chain and the mode branch once; `f` is a three-field copy
  (`ω̇`, `ẋ_idle`, `ẋ_frc`). Under the orthodox split, `f` would reproduce essentially
  the whole `f_ode!` body — four lookups and the mode branch — ×4 RK stages per step.
- `f_step!`'s transitions become Tier-1 events with mixed predicate/threshold guards
  (§2.1); `fuel_available` becomes an ordinary port (state-derived at the fuel system,
  hence stage-1 — no loop).
- Forced publications: none — everything `f` reads was already in `PistonEngineY`.

**`PID`** (control.jl:431-471) and the C172X FCS — the discrete side's representative:

- The current update entangles outputs and next state by construction (`y_i = x_i`:
  this tick's integral-path output *is* the updated integrator state). Under §5.2 the
  law runs once in `g_s2`, publishing paths, saturation and the updated states;
  `h` is a three-field copy. Under the orthodox split, `h(z, u, t)` would reproduce
  the entire law per compensator per tick.
- **Discovered latent delay.** The FCS chains anti-windup: outer compensators take
  `sat_ext` from the inner LQR's `sat_out` (c172x_ctl.jl:332,345,...). Wired naively,
  this is a *genuine* tick-domain algebraic loop
  (`outer.output → inner.input → inner.sat_out → outer.sat_ext → outer.int_halted →
  outer.y_i → outer.output`), which the build correctly rejects. Today's code escapes
  it only through hand-managed call order: the outer loops read the LQR's `sat_out`
  *before* the LQR updates, silently consuming the **previous tick's** value — a unit
  delay that exists nowhere in the code, only in statement ordering. Under this design
  the fix is one visible wire: connect `outer.sat_ext` to the inner compensator's
  stage-1 port for the previous saturation (`sat_out_0`, published by the identity
  decoder). The delay becomes an explicit property of the wiring. (The loop and its
  fix are formalism-independent; the framework's contribution is refusing to let the
  ambiguity through, and the decoder's is having the delayed value already on a port.)

Both components passed without blockers, with zero publications forced beyond current
practice — the empirical basis for §5.2's claim that derivative/output overlap is the
domain norm and the decoder matches the codebase's grain.

---

## 13. Decision log (condensed)

| # | Decision | Rejected alternatives (why) |
|---|---|---|
| 1 | Hybrid causal formalism; two-tier events; projection; no DAE/SDE/per-step hook | DAEs (projection suffices), SDEs (shaping filters suffice), `f_step!` (step-size-dependent semantics) |
| 2 | Causal port-based paradigm | Acausal/MTK (fights interactivity, discrete logic, live introspection); hierarchical callables (rigor by convention — today's footguns); thin SciML library (nothing for GUI/logging/hierarchy to hang onto) |
| 3 | Taxonomy: hybrid continuous primitive + periodic discrete + assemblies; both mode factorings | Strict purity/no modes (loses reset maps; latch logic becomes wiring ceremony); uniform hybrid kind (intra-component ordering semantics murky) |
| 4 | Immutable value signals in a typed signal table | Shared mutable buffers (aliasing/staleness, concurrent-read hazards); mixed semantics (second thing to document/test) |
| 5 | Reject algebraic loops at build, explicit breaks | Implicit delays (silent math changes); per-step numerical solving (jitter, runtime failures) |
| 6 | Strict two-stage structural feedthrough (state decoder `g_s1` + all-inputs `g_s2`); component split as the refinement | N output groups with declared input subsets (declaration/validation surface for a case that never materialized); per-output declarations + tracer verification (under-declaration = silent wrongness); traced multi-pass (hot-path cost, branch unsoundness); component-atomic conservative (false loops everywhere) |
| 7 | Reduce-ports with canonical fold | Σ-junctions as default (arity/positional ceremony); contribution buses (invisible dataflow) |
| 8 | Function-valued environment signals + handle pattern | Resource injection (second composition mechanism, invisible); pre-sampling as mechanism (dependency inversion at struts) |
| 9 | Deep paths within owned assembly types only | Unrestricted deep wiring (breaks substitutability across generic boundaries); strict one-level (re-export ceremony) |
| 10 | Structured immutable state over framework-owned flat `Vector{T}` | Mutable views (aliasing, silent missing-ẋ); fully structured/no flat vector (same machinery needed anyway, loses standard integrator interface) |
| 11 | Eltype genericity on the continuous path, three-tier scoping | Float64-only + finite differences (FD noise, keeps hand-written state-space layer, no tracer) |
| 12 | Set-propagation tracers (global + sampled-local), diagnostic-only | Dual-based tracing (derivative-zero blind spot); tracing as scheduling input (soundness) |
| 13 | Immutable `z` in cells + workspace + snapshot idiom | Mutable discrete state (aliasing, snapshot cost); double-buffering (deferred; publication races) |
| 14 | Scoped allocation invariant, CI-enforced on the hot path | Blanket dogma (fights logging reality); no policy (loses the type-instability canary) |
| 15 | Fused evaluation: `f`/`h` are pure signal consumers reading the fresh table (own `y` included); single computation site for derivatives and outputs | Mutable caches between `f` and outputs (S-function/FMI style — hidden state, purity violation); accepting duplicate computation (drift bug class: edited law in `f`, stale copy in `g`); derivative binding (superseded — a declaration feature subsumed by `y`-access); orthodox `f(x,u)` + atomization as standard idiom (2× components/wiring for the domain-normal overlap case) |
| 16 | Uniform component interfaces via the `g_s1` state decoder (default identity publication of state and modes); guards and handlers read the fresh boundary table, with per-event re-decode (`handler → project → g_s1 → g_s2`); `project` is the sole raw-state function (schedule-structural); unlisted-port convention for interface noise | Passing state alongside `y` (double-passing; two idioms in the wild); fully private discrete state (breaks uniformity; the codebase culturally publishes state anyway); one-handler-per-component-per-boundary restriction (superseded by the cheap per-event re-decode); exposing z *without* an unlisted convention (RNG/log noise) |
| 17 | Framework-owned simulation loop; stepper seam (advance by arbitrary `h` + on-demand dense output over the last step; one-step methods only); in-house fixed-step RK4/Heun as the sole first-cut backends; `OrdinaryDiffEq` dropped from dependency to possible future extension adapter | `OrdinaryDiffEq` as substrate with `CallbackSet` choreography (semantics by convention in a foreign event loop; demonstrated churn — the `task_local_storage` regression — in exactly the interactive multi-task usage); fused loop without the seam (loses the adaptive/stiff escape hatch for ~zero savings); multistep methods (history rebuild after every handler) |
| 18 | Tier-2 localization: lazy cubic Hermite dense output + bracketed derivative-free root-finding (ITP/Brent) on guard probes that run the sweep; post-event interpolant invalidation + remainder step + bounded event budget | Newton/AD localization (guards C⁰ not C¹ — kinks and σ′ = 0 stretches; discards the bracket certificate for local guarantees; negligible savings on rare microsecond probes); re-integration probes (4× cost; σ becomes trial-h-dependent); solver-matched high-order interpolants (only matter above order 4) |
| 19 | Harmonic tick grid on step boundaries; discrete stages gated to own tick instants (ZOH by construction); assemblies virtual for execution, rate scopes for declaration (integer multipliers `K ≥ 1` composing down the tree, compiled to absolute divisors); `comp.Δt` as single source of truth, no stored `Δt`-derived parameters | Atomic assemblies, incl. opt-in (coarsened schedulable unit → §5.3 artificial loops at assembly scale; interleaving protection meaningless under the signal table; FlightCore's whole-tree atomicity was a call-tree artifact); arbitrary tick periods via time queue (variable `h`, irregular frames, no demonstrated need); absolute-period declaration as default (welds deployment rates into reusable designs; base-period variables don't compose across independently authored assemblies); re-running discrete stages every boundary (un-samples sampled-data semantics); phase offsets (no demonstrated use); `Δt` via `h`-argument only (discretized laws live in `g_s2`) |
| 20 | Boundary event phase iterates to quiescence — rounds of full re-sweep → guards → handlers (declaration order, per-event re-decode) — each event firing at most once per boundary; due `h` updates run after quiescence, outside the iteration | Single pass per boundary (cascade latency N·h — step-size-dependent semantics, the §2.2 `f_step!` footgun class, made common by §3.1 externalized FSMs); bounded-rounds cap (arbitrary K knob; livelock burns the budget then errors instead of degrading to Tier-1 granularity); event/tick fixed-point iteration (structurally unnecessary — `z⁺` is invisible until the next tick decode) |
| 21 | Pacing outside the semantics (bit-identical paced/unpaced trajectories); piecewise-affine wall-clock map, anchor re-established at pace change and un-pause (debt cleared, counted); absolute deadlines with bounded debt + re-anchor on excess; `p = ∞` as explicit pacer-off; hybrid sleep-then-spin toward `deadline − margin`, with `margin` the single knob (0 = pure sleep, ∞ = pure busy-wait = FlightCore) | Relative deadlines (permanent sim-vs-wall slip); unbounded catch-up (burst after long stalls); keeping the anchor across pace changes (retroactively reinterprets elapsed history at the new pace); `p = ∞` as arithmetic limit (perpetual-overrun diagnostics under debt accounting); dedicated busy-wait mode flag (subsumed by `margin = ∞`); separate primitive-resolution threshold (absorbed into `margin` calibration) |

---

## 14. Open axes

To be settled in subsequent sessions:

- **Runtime periphery (axis 6) — next up.** GUI, input devices, network I/O, logging; the
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
