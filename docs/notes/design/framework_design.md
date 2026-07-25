# A Modeling & Simulation Framework for Flight.jl — Design Document

**Status:** seventh checkpoint (v0.7). Axes 1–6 settled; axis 7 parts 1 and 2 —
the component declaration layer (§13.1–§13.4) and the assembly declaration layer
(§13.5–§13.8) — settled, with periphery amendments to §4, §8 and §12 (write-side
granularity, root slots as exported faces, slot exclusivity, GUI liveness as a
derived property, stage-on-interaction, face-name/path notation split) grounded in
the completed interactive C172X case study (§14.4: inventory, surface walkthrough,
frame anatomies); build tooling, stopped-sim services and the migration outline
pending (see [Open axes](#open-axes)).

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

**Visibility.** Which ports exist at all is a declaration-layer decision: the output
contract *is* the public interface, and stage-function results outside it are private
intermediates — non-connectable, presentation-filtered, see §13.3. (Revision note: an
earlier version specified a presentational *unlisted* flag — skipped in logs and GUI
but still connectable. Retired in v0.5: it pretended privacy without enforcing it, and
its motivating case, RNG state feeding the component's own update, dissolved entirely
once the v0.5 prototypes let `h` read `z` directly.)

### 4.3 Table mechanics and port granularity

- **Scatter/gather is the whole protocol.** A stage function returns a named tuple;
  the framework scatters each field into that port's concretely-typed cell. Every
  reader — the next stage, `f`/`h`, guards, wired consumers, snapshot capture —
  gathers views from cells. The component's aggregate `y` is `merge(y_s1, y_s2)`
  *semantically* but virtual *physically*: reconstructed per call from cells (field
  loads, register-level, zero cost for isbits), never stored as an object. Name
  collisions across a component's stages are a build error.
- **Stage returns are named tuples of port values, period.** A custom struct is a
  first-class port *value* — one field of the returned tuple, one declared port, one
  cell (`pose = KinPose{T}`). Nested fields get no cells of their own; GUI and logs
  drill into them lazily (§4.2's view clause). Bare-struct returns are rejected:
  field-splatting would be ambiguous (one anonymous port vs. a splat), type-lossy
  under the merge (the struct type, and the methods that justify its existence,
  cannot be reassembled from a flat namespace), and reflection-hungry where the
  named-tuple form needs none.
- **Wiring is port-granular.** No sub-field connections: a consumer that wants less
  than a bundle asks the producer for a loose port, or takes the bundle and
  destructures. A field-projection connector is a guarded addition with an obvious
  shape, not built.
- **Granularity guideline** for authors: bundle what *shares a stage* (trivially
  enforced — each port has exactly one producing function) *and is consumed
  together*. Bundling across dependency footprints is the `KinData` mistake (§14.1:
  pose is stage 1, velocity-derived quantities are stage 2 — it must split). Fan-out
  is free, so publishing both a bundle and a hot loose field (`pose` *and* `q_eb`)
  is legitimate — one extra isbits cell.
- **Write-side corollary** (v0.6, from §14.4): **bundle what is written
  together.** The port is the atomic unit of the entire periphery — one cell, one
  root slot, one staged write, one device claim (§12.3), one trace address, one
  GUI liveness verdict (§12.5). Data written by different external writers, or at
  different cadences, must not share a port: pilot commands are scalar faces under
  a namespace prefix, and the convenient bundle is assembled *downstream*, inside
  the graph, by an ordinary component (single producer, consumed together — legal
  by the read-side rule). The two guidelines compose into one principle: a port's
  granularity is set by the finest-grained party owning either end — producers on
  the read side, external writers on the write side. Field-addressed staging (a
  lens into struct slots) stays a recorded guarded addition, unbuilt.

---

## 5. Evaluation order and feedthrough

### 5.1 The scheduling problem

At every evaluation instant, all signals must be computed consistently: every consumer
reads values already produced at that instant. Build the directed graph of (a) wiring
edges and (b) intra-component feedthrough relations; if acyclic, a topological sort
yields a **static evaluation schedule**, computed once at build time. The hot loop runs
a flat list of `(component, stage)` entries — zero runtime graph logic.

### 5.2 Two-stage outputs and structural feedthrough

Every component provides exactly **two output stages**, and feedthrough is declared
**structurally, by function signature** — there are no dependency annotations anywhere
in the design:

```julia
# continuous component                      # discrete component
y_s1     = g_s1(comp, x, m)                 y_s1 = g_s1(comp, z)
y_s2     = g_s2(comp, x, m, u, y_s1)        y_s2 = g_s2(comp, z, u, y_s1)
ẋ        = f(comp, x, m, y, u, t)           z⁺   = h(comp, z, y, u, t)

# event system (continuous side only) — same fresh table, same state views:
fired    = guard(comp, x, m, y, u, t)
(x⁺, m⁺) = handler(comp, x, m, y, u, t)     # constructs successors; may reset x
x⁺       = project(comp, x)                 # manifold projection
```

**The argument rule (stores and views).** A function's arguments are zero-copy views
of the distinct stores it reads: own state (`x`, `m` from the flat buffer and mode
cells; `z` from its cell on the discrete side), own published signals (`y`, gathered
from own table cells), inputs (`u`, gathered from foreign cells through the wiring's
name binding — composition information only the build knows), and the clock. The
signal table holds only *produced* signals, never transported ones: each datum has
exactly one home — buffer for `x`, cells for `z`/`m`, table for signals — and no
store mirrors another. Every argument earns its place as a view genuinely read, and
no further "simplification" of these signatures exists that does not introduce a
copy: eliminating `u` would mean republishing foreign cells under local names;
eliminating `x` was the pre-v0.5 identity transport, retired for exactly that reason
(§9.4, step 4).

- **`g_s1` is the no-feedthrough stage** — defined entirely by what it cannot see:
  it receives no inputs, so "no feedthrough" is unfalsifiable, and that structural
  guarantee is what its ports contribute to the schedule (they break would-be
  loops). It exists when the component has state-derived ports or shared
  state-derived intermediates for `g_s2`; otherwise it is simply absent. **A
  declared output that matches a state or mode field by name and type, and that no
  stage produces, is auto-published** by the framework from the state stores at
  stage-1 position (§13.3) — the useful residue of the retired identity default,
  now driven by the public contract instead of publishing everything.
- **`g_s2` receives all wired inputs plus `y_s1`** — its own stage-1 results, so
  shared intermediates are computed once, not re-derived — plus the state views;
  conservatively, every stage-2 output is presumed dependent on every wired input.
- **`f` and `h` run after the sweep**, when the full signal table — including the
  component's own `y_s2` — is complete and fresh. The fused idiom stands: compute
  each law once, in a stage; publish it; let `f`/`h` copy from `y`. The interfaces
  *reward* single-source-of-truth (nothing ever needs computing twice) rather than
  claiming to make duplication unwritable — a claim the pre-v0.5 prototypes could
  not honestly make either, since `f` always had `u` and the published state.
- **Guards and handlers read the same fresh world.** At a step boundary, the order
  is *integrate → project → boundary sweep → guards*, so by guard/handler time `y`
  is a fresh decode of exactly the state being transformed, and the state views are
  that state itself. Handlers construct `(x⁺, m⁺)` from raw state naturally — a
  reset map is `(; x..., ω = 0.0)`, no reassembly from published fields.
- **`project` runs between a state write and its decode** (after integration, and
  after any handler `x`-reset) — the only positions in the schedule where no fresh
  `y` of the new state can exist yet. Since v0.5 it is no longer *unique* in
  receiving raw state, but it remains unique in that schedule position.
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
`ẋ = f(x, u)`, `y = g(x, u)`; this design's `f` receives the orthodox arguments
*plus* the published table: `ẋ = f(x, m, y, u, t)`. The composite map `x ↦ ẋ` is
mathematically identical (linearization, trim and AD are untouched); the heterodox
element is only that derivatives may read outputs. The teaching line: *"stage 1
publishes what you know from state alone; stage 2 adds what needs inputs; your
dynamics read your own published results instead of recomputing them."* The decision
was grounded in a component-by-component survey of FlightPhysics/FlightApps (§14.2):
derivative/output overlap is the *norm* in this domain (Newton–Euler, kinematics,
piston engine, gear friction, every discrete compensator), so the orthodox split
would force either systematic duplicated math (with its silent-drift bug class),
systematic component atomization, or a cache mechanism — all rejected. FlightCore's
fused `f_ode!` already embodied the same economics; this design keeps them while
adding checked scheduling.

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

## 6. Aggregation: explicit summing junctions

N-to-1 physical aggregation (total wrench, total mass properties, total internal
angular momentum — today's generated `get_wr_b`/`get_mp_b`/`get_hr_b` tree walks) is
expressed by **ordinary junction components and explicit wires**. There is no
framework aggregation mechanism: no multi-connection ports, no declared fold ops, no
identity-element opt-outs. Every input port takes exactly one connection, everywhere.

```julia
struct SumJunction{V, N} end        #value type, arity; library-provided

inputs(::SumJunction{V, N}) where {V, N} =
    NamedTuple{ntuple(i -> Symbol(:in, i), N)}(ntuple(_ -> V, N))
g_s2(::SumJunction, u, y_s1) = (; Σ = +(u...))
```

- **Every mistake is loud** under the declaration layer: a forgotten contributor is
  an unconnected-input error naming `in4`; a double-wired slot violates
  single-connection; a stale arity surfaces as one or the other. The bookkeeping is
  ceremony, never silence.
- **The aggregate is a first-class signal.** `wr_sum.Σ` is an ordinary port:
  loggable, GUI-visible, fanned out to a second consumer (a loads monitor) for one
  wire. (Under consumer-side folding, the total was ephemeral gather scratch.)
- **Aggregation logic is arbitrary `g_s2` code** — mass-properties composition with
  its transport terms, weighted blends — not restricted to a declared
  commutative-associative binary op.
- **Fold order is author-visible**: the positional order of the junction's inputs.
  Reassigning contributors to different slots changes summation order, hence bits
  (float non-associativity) — deterministic per configuration and under author
  control, which is strictly more explicit than a framework-canonical order.
- For the handful of real sites, a **named site-specific junction**
  (`inputs(::VehicleWrenchSum) = (aero = …, ldg = …, pwp = …)`) documents the
  contributor set better than generated slots, at the price of hard-coding it into a
  type; the generic positional form remains the tool for configuration-variable
  sites. Both are plain components; the framework is not involved.

**The hierarchical aggregation idiom** (what replaces the tree walk). Only physical
contributors publish these ports — a strut publishes `wr_b`, avionics publishes
nothing — and each assembly that *owns* contributors aggregates them with an internal
junction and **exports the total** (§3.3: the junction is a component inside the
assembly; the assembly exports its `Σ` port). The §8 connection rules force this
shape: a generically-held submodel is opaque, so every generic boundary must export
its aggregate. `Ldg` sums its three struts and exports `wr_b`; the systems assembly
sums `aero + ldg + pwp`; the vehicle wires the systems totals into Newton–Euler.
Each recursion step of FlightCore's tree walk becomes one visible junction at the
level that owns the contributors — for the C172, about four junctions and fifteen
wires, written once, reading as a manifest of what weighs, what pushes and what
spins. Frame responsibility is unchanged: contributors publish in the common body
frame, applying their own mounting transforms at source. Do **not** bundle the three
quantities into one contribution struct: contributors are ragged (aero has wrench
but no mass, fuel the reverse, only `pwp` has angular momentum), and a bundle forces
zero-filled identity noise through every port — the "silently sum nothing" hazard in
a new coat. A `sum_ports!`-style helper (instantiate + wire + export in one call) is
guarded-addition sugar, added when migration shows the pattern repeated.

The ledger against FlightCore's tree walk, recorded: its zero-wiring convenience and
its worst failure mode were the same property — a contributor with a forgotten trait
method contributes *silently nothing* (a lighter vehicle, no diagnostic, ever).
Explicit wiring inverts every silence into a warning or error, and makes
per-contributor values and intermediate totals observable ports; the cost is that
adding a deep contributor edits one assembly level (its owner's wiring) instead of
zero. Rejected at every revision: **contribution buses** (today's mechanism
portified — dataflow invisible in the graph, scoping rules, accidental
contributions).

*Revision note (v0.5).* v0.1–v0.4 specified **reduce-ports**: consumer-declared
commutative-associative folds with multi-connection legality, canonical fold order
and identity-element opt-outs. Reversed during axis 7 (row 7 → row 37): the use-site
census never grew beyond the three Newton–Euler aggregations in one library
assembly; the mandatory typed declarations of §13 made junction mistakes loud,
dissolving the original "positional ceremony" objection; and `Reduce` had become the
declaration vocabulary's last remaining wrapper — killing it left `inputs()` as bare
types with zero framework vocabulary and retired the canonical-fold and
identity-element rules wholesale.

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
  `systems/ldg/{left,right,nose}/trn_field` in one visible block — no per-level
  re-exports); you may only connect port-level into submodels held **generically**
  (`World` connects `terrain/field => aircraft/trn` and knows nothing more). This kills
  the re-export ceremony where it is ceremony, and preserves substitutability where the
  boundary is load-bearing. Operationally: a path may traverse any chain of
  concretely-typed fields and stops at the first generically-held child, whose
  exported faces are the only things addressable beyond that point — a rule about
  the *declaration's* knowledge, not the build's (a deep path into a generic child
  is forbidden even where the concrete instantiation would resolve it, because it
  hard-codes one implementation and breaks on substitution).
- Paths are validated at build time; renames break loudly.
- Fan-out is free (one producer, many consumers). The converse is strict: every
  input port takes **exactly one** connection, no exceptions (aggregation is
  junctions, §6). The rule spans levels: an input fed both inside a sub-assembly
  and by an ancestor's deep route is a two-producers build error — deep routing
  cannot silently double-feed.
- No auto-bubbling of unconnected inputs (implicit interface growth — and worse:
  the forgotten-wire error, §13.4 walkthrough 2, would be silently *promoted to
  interface*, climbing level by level into a live root slot nobody feeds, instead
  of failing at build).
- Unconnected output ports: build-time warning. Unconnected input ports: build error
  (no silent defaults). **The check is a whole-tree property, not a
  per-declaration one**: within a single assembly declaration an unfed child input
  is simply *awaiting a claim from above* — a sibling wire, an ancestor's deep
  route, or an `exports` entry handing the obligation up one level (§13.6). The
  error fires at the root build for any input whose obligation chain never
  terminates. The one legitimate terminus fed by no component is the root
  assembly's own exported input face — a root slot (§12.3).

---

## 9. State and data representation

### 9.1 Continuous state: structured immutable, flat backing

Each continuous component declares its state as an **isbits struct whose leaves are all
real scalars of a common eltype `T`** (`SVector`s, quaternions, `Ranged` values —
ultimately reals; `Int`s/enums/`Bool`s belong in modes). The framework:

- computes a **flat layout** at build time (compile-time offsets over one contiguous
  `Vector{T}` buffer it owns);
- **reconstructs** the typed immutable state value for a component at each evaluation
  and passes it to every function receiving state views (§5.2 argument rule — field
  loads at known offsets, register-level, zero cost);
- receives immutable results back: derivative functions return an `Ẋ`-typed value
  (scatter-stored into the flat `ẋ` buffer); event handlers and projection return a new
  `X` (written back).

**The buffer is authoritative; typed values are ephemeral reconstructions.** Nobody
outside the framework ever holds a mutable reference to state. "Ephemeral" is
literal: an isbits view materializes in the caller's frame (registers or spilled
stack, the compiler's business) for exactly the duration of the call and has no
existence between calls — re-materializing is the same loads, value-identical
because the value is immutable and the buffer unchanged within a sweep. Whether the
executor rebuilds per call (natural for a data-driven schedule) or hoists once per
component per sweep (natural for an unrolled generated schedule) is codegen freedom,
semantically invisible. The complementary rule (§5.2): **one home per datum** —
buffer for `x`, cells for `z`/`m`, table for produced signals — and no store ever
mirrors another; in particular there are no state cells in the table beyond
contract-driven auto-published ports, which are interface, not transport.

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

The declaration layer makes this scoping visible per port: continuous components
declare outputs as *functions of the sweep scalar* (`outputs(::C, ::Type{T})`,
§13.2), so slot types per activity come from evaluating declarations at that `T` —
never from substitution machinery or from inference through user code.

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
4. **Decoder-exclusive state access → stores-and-views arguments (v0.5, axis 7).**
   Step 3's second half was reversed. Its justification — "state is published
   anyway" — quietly depended on default identity publication; once §13.3 made
   publication a deliberate interface act, the identity decode stood revealed as
   *transport*: copying the buffer into cells so a buffer view could be replaced by
   a slot view — ceremony of exactly the kind step 2 removed elsewhere, now with
   dead stores (post-v0.5, no own-function reads state through the table). The
   reductio that decided it: the same minimalist logic would pack `u` into `y_s2`
   and end at `f(comp, y, t)`, republishing foreign cells under local names. The
   fixed point is §5.2's argument rule — zero-copy views of the stores a function
   genuinely reads. What survives of step 3: the uniform shapes, the fused
   economics, and `g_s1` itself — no longer the sole state gate, but the
   no-feedthrough stage.

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
`sketch_decoder.jl` (merged form; both predate the v0.5 revisions — reduce-ports,
identity publication and the stateless prototypes appear in them) for the worked
example;
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
interferes (under the v0.5 prototypes stage-2 functions also receive state views, but
neither those nor `y_s1` are ever seeded — same conclusion). Stage-1 functions are never traced (nothing to
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
in the *stages*, not just `h`: per §14.2, the discretized laws that actually consume
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
axis 6 (§12).

---

## 12. Runtime periphery (axis 6)

GUI, input devices, network I/O and logging, and the concurrency model binding them
to the §11 loop. Settled in full: architecture and data exchange (§12.1–§12.6), loop
scheduling and thread policy (§12.7), the next-snapshot wait (§12.8), shutdown
(§12.9), and the mid-run mutation doctrine (§12.10).

### 12.1 No shared mutable model: staged writes, snapshot reads

FlightCore's periphery is one big lock: `SimControl` and the live `Model`, guarded by
`io_lock`, with one task per attached interface reading or mutating the model under it
(sim.jl). The lock does enforce §11.3's boundary-visibility rule — it is only ever
free between steps — but its costs are structural:

- **The loop's frame budget is hostage to its readers.** A slow GUI frame or a stalled
  `extract_output` holds the lock and the sim cannot step; blocking time is
  indistinguishable from overrun in any accounting.
- **Input timing is scheduler-determined and unrecorded.** Writes land between
  whichever boundaries the OS interleaving produced; there is no defined input trace,
  so §11.7's bit-identical replay is unachievable *in principle* for interactive runs.
- **It protects an idiom that no longer exists.** `assign_input!` and GUI widgets poke
  the live model; under the immutable signal table there is nothing to poke — the
  periphery needs a defined write path regardless.

The replacement has five planes:

1. **Staging (inbound):** devices submit pending input writes at any wall-clock
   moment, never touching live slots (§12.3).
2. **The drain:** exactly one point, at the top of each boundary sequence, where the
   loop takes the staged batches and applies them to the root input slots. Between
   drains the loop owns its data exclusively — no lock is held during stepping, ever.
3. **Publication (outbound):** at the end of each boundary sequence the loop publishes
   an immutable snapshot; readers observe it without coordinating with the loop
   (§12.2).
4. **Control:** pause/pace/stop on a separate few-word atomic surface (§12.6).
5. **Task topology:** one loop task, one task per attached device (the GUI's on the
   main thread — CImGui's constraint), start gate and shutdown protocol as today.

Two rules bind the implementation:

- **Every handoff is one atomic reference operation.** Both shared structures reduce
  their mutable surface to single words (release/acquire `@atomic` fields), and the
  GC is the reclamation mechanism: an object a reader still holds is reachable and
  therefore never recycled, which dissolves the reclamation problem (hazard pointers,
  epochs, RCU grace periods) that makes these patterns hard in non-GC languages.
  Deep immutability of the exchanged objects is what makes this sound.
- **No user code, no unbounded work, ever, inside a framework critical section.**
  FlightCore's pathologies are all "arbitrary code under a shared lock"; here
  mappings run on device tasks, rendering runs against snapshots, and the loop's
  frame contains framework and model code only. A stalled device produces stale
  snapshots and late staging; it cannot stall the loop. (Residual, pre-existing
  exposures: GC pauses and OS scheduling — the pacer's debt absorption is the
  mitigation.)

Consequence, recorded because it collapses an API axis: interactive and batch
simulation stop being different execution modes. A batch run is the same loop with
empty staging and no snapshot readers; a replayed interactive session is the same
loop with staging fed from a recording (§12.3).

### 12.2 Outbound: snapshot publication

The loop builds each snapshot — boundary-consistent signal table, `t`, framework
status (§11.7 diagnostics) — in private memory, then publishes it with a single
release-store to an `@atomic latest` reference; readers acquire-load and then work
with an immutable, coherent world for as long as they like. Wait-free in both
directions: a wedged reader cannot delay publication by a nanosecond; the loop cannot
tear a reader's view. Publication happens only after the boundary sequence completes
(§11.3 as extended by §11.6).

**Binding rule: nothing reachable from a published snapshot is ever written again.**
The table's immutable values (§4.1, §9) make the compiler enforce most of it; the
rule is what the soundness of lock-free reading rests on.

The captured table includes private intermediate cells (§13.3) — the copy is
mechanical, and they serve the author's own debug panels; presentation layers (log
export, GUI listings) filter to the public contract by default. **It also includes
the root slots** (v0.7, made explicit by §14.4): slots are source cells of the
table, not state stores, so they ride along — and this is load-bearing, not
incidental: the §12.5 peek's else-snapshot fallback is what an idle live widget
displays, and read-only mirrors of claimed slots (the axis sliders under joystick
claim) show the applied slot value from the snapshot. Slot values in the log are
derived data (recomputable from the trace), which is consistent — snapshots are
derived wholesale. The snapshot
deliberately does **not** carry the state stores (`x`, `m`, `z`): the state
trajectory is *derived* data — recomputable from the trace header plus the batches
(§12.3) by bit-identical replay — and per-boundary capture would systematically
record derived data, the same asymmetry the trace-default decision (row 29) refuses
in the other direction. "What was the private state at t = 37.2?" is answered by
replaying to 37.2 and inspecting the live stores; a state field wanted in logs or
GUI has the honest remedy of being declared public (one auto-published cell per
sweep). Post-run continuation reads the live stores directly; periodic full-state
checkpoints (warm restart without replay-from-zero) are a guarded addition shaped as
an opt-in log policy, and a dev-mode flag auto-publishing all state fields is a
possible diagnostic kin to the NaN-poisoned workspace.

Logging dissolves into publication: the log is a vector of retained snapshot
references — same objects, zero extra copies; the per-step `deepcopy` detour of the
`SavingCallback` disappears. Cost: one snapshot allocation per boundary, on the
framework side of the §9.5 scope (which already carved out logging); logged snapshots
are not garbage at all, unlogged ones die young. Rejected: preallocated snapshot
buffers (double/triple ring) — reuse reintroduces exactly the reader-liveness proof
the GC provides for free, to save an allocation profiling has not indicted.

### 12.3 Inbound: root input slots and per-device staging

**The write surface is root input slots** — and a root slot *is* the root
assembly's exported input face (§13.6): routed inward to consumers, produced by no
component, fed by the parent's wire at every non-root level — and at the root
there is no parent. No dedicated vocabulary survives (`add_input!` in the early
sketches is dead). Slots are sources to the build-time scheduler, constants within
a frame, and the *only* thing the periphery may write (the GUI reaches them
through §12.5's resolution; control commands are not writes, §12.6); devices,
mappings, the trace and the GUI write path address them by **face name** (§13.6):
structural slash paths never cross the periphery boundary — the periphery speaks
the root contract's names only.

**Slot exclusivity: one writer per slot at any time** (v0.6, from §14.4). A
device claims its slots at attach; claiming an already-claimed slot is an
attach-time error, and detaching releases the claims (a released slot's GUI
widgets re-enable, §12.5). This supersedes the cross-device conflict *policy* —
attachment-order precedence at drain — which resolved races the case study shows
nobody wants: every dual-writer field in the C172X demo is a joystick stream
shadowed by a GUI mirror, where simultaneous live writing is a bug. Per-device
cells, the CAS merge and the atomicswap drain all stay — they serve atomicity and
coalescing, not arbitration.

**Slot initial values are owned by the init/trim services** (v0.7, resolved by
§14.4). Input declarations are bare types (§13.2) and carry no defaults, but a
slot unfed by any device must hold a defined value from the first frame (today's
`U()` constructors provide these: `mixture = 0.5`). Export-entry defaults were
rejected: the trim service writes slot values it *solved for* (throttle,
elevator) — not declaration constants. `init!` establishes every slot and the
trace header captures the result; the concrete service spelling remains with §16.

**Staging: one atomic cell per attached device**, in attachment order fixed at build.
Each cell has a single writer — its own device task — and holds that device's latest
pending batch of slot writes:

- a *complete* writer (a joystick: full write-set every poll) overwrites its cell;
  coalescing between drains is ZOH-correct at boundary granularity;
- a *sparse* writer (the GUI: only what the user touched) accumulates via CAS merge
  on its own cell; the CAS can fail only because a drain intercepted the old batch,
  so the retry is bounded and the failure case is precisely correct — intercepted
  writes are already applied and must not be re-staged.

**Doctrine: staged values are levels, never deltas** (`press_count = 17`, never
`presses += 1`) — levels are idempotent and survive coalescing; button edges ride as
monotonic counters.

**The drain** (frame top): one `atomicswap(cell, nothing)` per device — an
indivisible take, no lost-write window — applied **in attachment order**, retained
as a deterministic application order (under slot exclusivity, cross-device writes
to one slot no longer arise; the order matters only for diagnostics). Which
*frame* a write lands in remains wall-clock reality; what the drain guarantees is
that the frame's outcome is a pure function of the drained batches.

**Mappings run on the device task**: today's `assign_input!(mdl, mapping, data)`
becomes pure `map_input(data, mapping) → batch`. User-extensible code thereby never
executes inside the loop's frame, and the trace consists of slot-level batches.

**Mappings are binding data, not shaping code** (v0.6, from §14.4). A mapping is
a declarative table — axis/button → slot path, plus per-axis conditioning
parameters (deadzone, expo strength) applied by a shared pure helper on the
device task. The boundary is set by the face contract: **a face's meaning is
writer-independent**, so faces carry *post-conditioning* semantics — a GUI slider
or a script writes the same command a curved stick delivers (running a mouse drag
through a deadzone would be absurd); this GUI-parity test is what places
conditioning upstream. Aircraft-semantic derivation (the C172X `q_ref = q_sf ·
axis` fan-out) must *not* ride along: it is FCS design and lives in-model — in
the avionics, or accepted as a small per-aircraft×device mapping entry (an
aircraft-design fork, §14.4). The trace records post-conditioning levels —
exactly what the model consumed, so replay is exact; raw-stick provenance (re-run
a session through *different* curves) is the known, accepted loss. Edge logic
follows the levels doctrine: devices stage monotonic press counters; accumulators
(trim offsets, flap detents) are model state, not mapping state (§14.4).

**The input trace** is the sequence of drained, device-tagged batches per frame. It
extends §11.7's determinism end-to-end: replaying a recorded interactive session —
staging fed from the recording, no devices or mappings present — reproduces the
trajectory bit-identically.

**The trace header captures the full initial state** `(x, m, z)` **plus the
initial root-slot values** at `init!` (v0.7 — an unfed `mixture = 0.5` never
appears in any batch, so replay is broken without them; the init/trim services own
slot initialization, §16, and the header capture extends naturally) — the one
full-state capture in a normal run, and the other half of what "given the
initial state and the trace, the log is recomputable" requires. Header plus batches
are the *primary* record; everything else, the state trajectory included, is
derived (§12.2).

**Trace recording is on by default** (cleared at `init!`, retrievable after the run,
plain kill switch for memory-constrained marathon sessions). The asymmetry that
decides the default: the trace is *primary* data and the log *derived* — given the
initial state and the trace, the log is recomputable (that is what bit-identical
replay means), while an untraced interactive session is unreproducible, permanently.
The cost supports it: the trace retains batches the devices already allocated, at
drain-rate × device-count — tens of MB per hour worst case, two orders of magnitude
below the snapshot log. No sampling, no rolling window (complexity without a
customer).

Rejected shapes (both torture-tested in §14.3): **per-slot atomic cells** — the
simplest (no merge machinery, and a per-slot layout cannot lose independent writes)
but same-slot conflicts resolve by hardware store order, i.e. sub-frame wall-clock
phase (run-to-run behavioral variance, §14.3), peeks are cross-device, the trace
loses provenance, and wide slot types hit Julia's atomic-width lock fallback; **a
shared lock-free batch stack** (CAS-push, swap-drain) — whole-batch atomicity and
the richest trace, but conflict order is still temporal (push order), and pending
memory is unbounded while paused, taxing every peek that must walk the chain.

### 12.4 Device taxonomy: one kind

FlightCore's input/output/GUI trichotomy is lock choreography, not modeling:
`get_data!` may block, so it runs outside the lock; `extract_output` must not block,
because it runs inside; the GUI breaks both rules at once and gets a third interface
(`render!` under lock, `sync = 0` so VSync cannot stall the sim, a manual framerate
sleep). With no lock, the protocol the taxonomy encoded has no referent.

**Every attached device receives the same handle**, carrying the two primitive
capabilities — read (latest snapshot; optionally wait-for-next-boundary, §12.8) and
stage —
plus control access (observe running, request shutdown; `should_abort` stays a
per-attachment flag). Input-only and output-only devices are degenerate uses, not
framework kinds; a bidirectional network peer is *one* device with one socket and one
lifecycle, not two framework devices sharing state. The GUI is an ordinary device —
the paradigm one, using every capability — with exactly two genuine peculiarities,
neither taxonomic: main-thread affinity (a launch concern) and read-modify-write
widgets (§12.5).

### 12.5 The GUI write path: port resolution, peek, active widgets

Panels remain per-component extensions in FlightCore's style — `GUI.draw!(ctx,
::LowPassFilter)`, discovered by walking the assembly — but widgets name **the
component's own ports**, never root slots. The build-time wiring answers, statically
and exactly: *is this port transitively driven by a root input slot, and which one?*
Every input port has exactly one source (§5), so the resolution is total:

- **root-driven → live widget**: peeks and stages the resolved slot through the
  GUI's own cell;
- **component-driven → read-only rendering**: displays the driven value from the
  snapshot, visually distinct, with the source as provenance ("driven by
  `avionics.throttle_cmd`").

This retires FlightCore's dead-slider convention (the `Cessna172Xv1` throttle: the
engine panel's slider visually live, silently overwritten by the avionics every
cycle — who commands what living in the user's head) and replaces it with checked
structure: **a widget is live exactly when the underlying input is yours to command
in this configuration.** User-commandability is a wiring decision made where
configurations are made; command-plus-manual-override is a mux component with a
root-wired select — explicit structure, not two writers racing (the same race as
§14.3's drag phase, retired by the same rule). The obligation this places on the
GUI: read-only rendering is first-class, not an error state — the author of
`input_slider!` cannot know at authoring time whether it will be live.

**Liveness is a derived property, and resolution is transitive** (v0.6). A widget
is live iff its port's feed chain — walked through wires and exports across *all*
levels, not just the local assembly — terminates in a root slot, *and* that slot
is currently unclaimed by a device (§12.3 exclusivity). There is no per-port
"GUI-controlled" marking anywhere: the export chain is the marking, written by
the one author entitled to write it (a component's ports become GUI-commandable
exactly when the assemblies above surface them). The switch between "driven by
its own panel" and "driven by an external provider" is therefore automatic — at
build time by wiring archetype (a scripted `World` wires a scenario component
into the same faces the interactive `World` exports to root), at run time by
device claim state. Rejected (v0.6): nominally-connected ports with a GUI
*override* channel — a second write path that breaks the
pure-function-of-drained-batches frame semantics, needs a parallel trace and
replay mechanism, and cannot resolve the producer conflict (either a dead widget
or a silently discarded wire); made safe — staged, frame-top, traced, exclusive —
it collapses into the root-slot mechanism it tried to bypass. The honest cost
stands: **unexported ports are unpokeable** — FlightCore's poke-any-`u` workflow
does not survive contract visibility (§13.3), deliberately.

**Peek rule:** a widget displays its **own pending write if any, else the snapshot
value**. Own-cell only: another device's pending write is invisible by design (its
applied value arrives via the snapshot one frame later); cross-device peek would
re-couple devices for sub-perceptual benefit. While paused, staged edits display
indefinitely and apply at the un-pause drain. Fan-out is consistent for free:
widgets on ports resolving to the same slot peek the same pending value.

**Staging contract: widgets stage on interaction events only** (v0.7, superseding
the active-widget stage-every-pass contract). Value widgets (sliders, drags) stage
the new absolute level on edit; edge widgets (buttons) stage on activation, as a
level computed from the peek — a flaps button peeks the current counter `k` and
stages `k+1`. The levels doctrine makes this safe by construction: repeated staging
of the same level within a drain window is idempotent (no repeat-increment hazard),
and multi-click within one window counts correctly through the own-pending-first
peek (`k` → stage `k+1`; second click peeks pending `k+1` → stages `k+2`). Held
buttons do not re-stage — after the drain applies and the snapshot catches up,
re-staging from the peek would auto-repeat at frame rate; the activation edge is
the intent.

The superseded contract's motivation died with slot exclusivity (§12.3): stage-
every-pass existed to win every drain against a streaming device sharing the slot
for the grab's duration (§14.3's drag phase) — but a slot the GUI can write is now
by definition unclaimed, so once staged and drained a value simply *stays*; there
is nothing to reassert against. Nor is it worth keeping as insurance: if an
anomalous writer ever touched an unclaimed slot through a framework defect, the
slot visibly fighting is a diagnosable anomaly — continuous re-staging would mask
the invariant violation at render rate. Side benefit: staging traffic (and trace
noise) drops from render-rate-while-grabbed to actual edits. **Claim-transition
policy:** if a device claims a slot mid-interaction (attach during a drag — a
deliberate act concurrent with a held grab, vanishingly rare), the widget flips
read-only at the next render and the drain discards stale GUI entries to
newly-claimed slots with a warning.

### 12.6 Control plane

Pause, un-pause, pace changes, stop: a few scalar fields on a separate atomic
surface, consulted by the loop at frame top and inside its wait and pause states.
**Not staging, structurally:** staged writes apply at drains, and a paused loop
drains nothing — un-pause via staging would deadlock by construction. Riding outside
the drain/trace path is safe for determinism precisely because §11.7 put pacing
outside the semantics: control changes *when* frames execute, never what they
compute (stop merely truncates the trajectory). While paused the loop blocks on a
condition (notified on un-pause and stop), not a spin.

### 12.7 Loop scheduling: wait primitive, yields, thread budget

**Coarse phase = task-yielding `sleep`; no `systemsleep` variant.** The precision
argument for `Libc.systemsleep` (≈0.5 ms vs ≈1.4 ms median overshoot, §11.7) is
worth ~1.5 ms of `margin` — a few percent of one core at 50 Hz. Against it: `sleep`
releases the loop's thread, making the pacer's wait slot the natural scheduling
window for co-resident device tasks (the design already spends that slot twice:
§11.7's staging slot, §12.3's drain source); `systemsleep` turns the slot into dead
time and makes the periphery correct only when every device task has its own
thread — resurrecting FlightCore's hard `nthreads` check as a correctness
precondition. The failure asymmetry decides: `sleep`'s worst case is a late
wake-up → an overrun → absorbed as debt and *diagnosed*; `systemsleep`'s is starved
device tasks → silent functional degradation. And §11.7 committed to `margin` as
the single knob; a primitive selector is a second knob hiding inside the first (its
two settings differ by a margin recalibration). A `systemsleep` variant for
dedicated-thread hard-RT deployments is a guarded addition.

**Yield rule: with devices attached, every frame yields at least once** —
implicitly via the coarse-phase `sleep` when it runs, via an explicit `yield()`
otherwise (unpaced runs; pure-spin frames with budget ≤ margin). Zero semantic
cost: pacing, hence yielding, is outside the semantics (§11.7). The spin phase
itself never yields — that would trade its µs precision for scheduler noise.
Consequence: the loop occupies a thread for at most one frame before the scheduler
can run anyone else, so the thread-monopolist precondition for Julia's
cooperative-scheduler freeze (a never-yielding task holds its thread forever) is
structurally absent from framework tasks.

**Thread budget: a documented sizing rule and a startup warning, not a hard
error.** FlightCore's `nthreads` error was honest for its architecture: under-lock
blocking plus a busy-wait pacer wedges *totally* when tasks share threads — the GUI
queues on the lock behind the stall, the window freezes, no escape and no message.
Here the freeze cannot reproduce: the loop yields every frame, nothing couples a
stall to anyone else (the GUI waits on nothing, ever — snapshot acquire-load, own
staging cell, atomic control), and the GUI runs on the *calling* task, so it cannot
fail to be scheduled: under any starvation the window keeps rendering and the stop
button keeps working. Undersized sessions degrade to laggy inputs and stale
snapshots — visible, recoverable states. `run!` warns when `Threads.nthreads()` is
tight for the attached population, naming the `julia -t` remedy; sizing guidance:
one thread for the loop, the main thread for the GUI, headroom for compute-heavy or
blocking-ccall devices (libuv-backed I/O yields; raw blocking ccalls pin their
thread for the duration). No pinning, no sticky tasks.

**Liveness heartbeat.** Since starvation is survivable it must be diagnosable: the
published framework status includes per-device liveness (last-staged / last-read
wall time, task state) next to the pacer diagnostics. A starved, blocked or crashed
device task shows in the GUI as a stale heartbeat with a name on it, not as
mysteriously frozen physics.

### 12.8 The next-snapshot wait

Rate-matched output devices (telemetry, disk streaming) act once per boundary
without polling: a monotonic frame counter published with the snapshot, plus one
`Threads.Condition`. The loop's publication is `lock; counter += 1; notify;
unlock` — nanoseconds of framework-only code, never blocked by waiters (a waiter
parked in `wait` has released the lock as part of parking). The device-side
`wait_next_snapshot(handle)` blocks until `counter > last_seen && running` under
the canonical predicate-loop idiom, which handles waiters at different paces,
frames skipped while transmitting, and shutdown (§12.9 wakes all waiters; each
predicate routes its owner out) with no per-frame reset. An `Event` latch is the
wrong primitive here: recurring signals require un-latching, and the reset has no
correct placement under asynchronously arriving waiters (the `io_start` reset
comments in FlightCore's sim.jl document the once-per-run version of exactly this
race). Conditions carry no facts, only "look again"; the facts (counter, running)
live in state each waiter tests privately.

**Semantics: newest-wins, no queues.** A slow consumer skips frames and always
receives the current world. This mirrors the inbound side: coalescing to the
newest batch (in) and to the newest snapshot (out) are the same ZOH decision; no
backpressure exists in either direction, and the loop never waits on anyone.
Rejected: per-consumer every-boundary queues (unbounded under a slow consumer —
the batch stack's pause pathology again; complete history *is* the log). The GUI
does not use the wait (VSync-paced, it reads `latest` each render).

### 12.9 Shutdown protocol

1. **Initiation:** `t_end` reached, or a control-plane stop (GUI, device handle,
   code). The loop always completes the current boundary sequence — never stops
   mid-frame — publishes the final snapshot, then sets the sticky stopped status.
   Publishing first guarantees output devices can flush the true final state.
2. **Wake all framework waits** (next-snapshot, pause): waiters observe the
   stopped status and unwind — a stop while paused therefore works.
3. **Unblock device-specific blocking calls** via an `unblock!(device)` hook,
   default no-op; a network input's override closes its own socket, raising in
   the blocked task (caught by the framework device loop, treated as shutdown).
   This demotes FlightCore's EOT convention from load-bearing shutdown mechanism
   to an optional wire-protocol courtesy between remote peers.
4. **Device loops exit:** `while running(handle)` with all blocking points
   interruptible per (2)–(3); `finally shutdown!(device)` guaranteed.
5. **Join with a timeout:** a device task exceeding it is reported *by name*
   (§12.7 heartbeat) and abandoned with a warning rather than hanging `run!`.
6. **Device-initiated paths:** `should_close` (window ✕, peer EOT) exits the
   device's own loop; with `should_abort` set it also requests a sim stop,
   otherwise the sim continues with the device absent (its cell stops filling —
   the loop is structurally indifferent). A crashing device task is caught by the
   framework wrapper and follows the same path, logged with the device's name.
7. **Loop-side failure** runs (1)–(5) from the catch path (the
   `SimulationTermination` machinery is the precedent), so devices unwind cleanly
   regardless of who died.

### 12.10 Scripts and the mid-run mutation doctrine

What the consumers demonstrably mutate mid-run, surveyed: FlightCore's
`user_callback!` has exactly two archetypes — the timetable script
(c172_demos.jl:290: `elevator_offset` as a function of `t`) and the synthetic
pilot (c172_demos.jl:423, 525: a phase FSM reading `y` and writing mode requests,
references, flaps, wind). Both write only `u` fields; no demo, test or GUI path
pokes `x`/`z`/`s` mid-run, and `init!`/trim appear only between construction and
`run!` (c172_demos.jl:303).

**Sim-time scripts are model behavior: scenario components.** Both archetypes are
clocked by *sim time* (`t`, the trajectory). Mapping them to devices fails
outright unpaced — the demos run at `pace = Inf`, where frames take microseconds
and a wall-clock task's staging lands at scheduler-determined sim times,
differently every run. The clock is the criterion: **sim-time scripts →
source/supervisor components** (periodic discrete, `K = 1` for today's
`dt = 0.02` callbacks), executed synchronously in the loop, deterministic paced or
unpaced, replayed by recomputation with no trace; **wall-clock interactions →
devices**, traced and replayed from the trace. The component mapping is strictly
richer than the callback it replaces: the `Ref(:init)` phase closure becomes
honest `z` (visible in snapshots, logs and plots); inputs arrive same-boundary
fresh by topological order (the callback ran post-step, one boundary staler); and
in a scenario configuration the script drives the avionics' input ports, so §12.5
renders the corresponding GUI widgets read-only with provenance — today's
demo-vs-GUI dead-slider fight, resolved by the port-resolution rule.

**`user_callback!` is eliminated.** It is the periphery's `f_step!`: arbitrary
unrecorded mutation, ordered by convention, invisible to replay (§2.2). Its
historical justification was FlightCore's composition cost — a supervisor required
a full `System` declaration against a ten-line closure; this framework prices a
component at roughly the closure's weight, removing the pressure. Its call sites
migrate to scenario components, not devices.

**Manual event triggering needs no mechanism:** a root input slot plus a Tier-1
guard reading it (levels doctrine: latched commands or counters), already
expressible in settled machinery — the demos' engine start/stop buttons are
`u`-writes today.

**Mid-run re-initialization is not built, because it is not demonstrated.**
Initialization and trim are stopped-sim workflows (axis 7's first-class services),
where no concurrency perimeter exists — no loop, no devices, plain single-task
code. The guarded-addition shape is on record should demand appear: a traced,
boundary-executed intervention command applied through project → sweep → publish,
so no consumer ever observes un-decoded state.

**The doctrine, final form:** while a simulation runs, the periphery stages
root-input writes and issues control commands — nothing else, structurally.
Anything that wants to poke the model mid-run is an *input* in disguise (wire a
slot and guard), *model behavior* in disguise (add a scenario component), or a
*wall-clock interaction* (attach a device).

---

## 13. The declaration layer (axis 7): components and assemblies

How an author spells a component: where the structural facts live, what the build
takes as authoritative, and what is checked against what. §13.1–§13.4 settle the
component side (v0.5); §13.5–§13.8 settle the assembly side (v0.6); the build
pipeline itself and the stopped-sim services remain open (§16). Concrete syntax
below is near-final in shape but still illustrative in spelling. The sketches
(`sketch.jl`, `sketch_decoder.jl`, `sketch_io.jl`) predate this section — they
still show reduce-ports, identity publication, the stateless prototypes and the
builder-style assembly — and will be refreshed after the §14.4 walkthroughs.

### 13.1 Position: a declarative trait layer — plain Julia, no macros

Stage functions are ordinary multiple-dispatch methods (the `GUI.draw!` precedent).
Structural facts are declared through a small set of well-known functions returning
plain values, defined alongside the methods. No macro DSL: generated code is opaque
to the debugging, tooling and comprehension workflows the charter protects (§1), and
a macro can only ever *lower to* a layer like this one — so a convenience macro
remains addable a posteriori as pure sugar (the `@kwdef` precedent), while never
becoming load-bearing. Redundancy between declarations and function bodies is
accepted deliberately, under one non-negotiable condition: **every inconsistency
fails loudly**, at build time where possible, at first execution otherwise.

**The schema-authority principle.** Declarations *define* the model's structure;
evaluation *checks* conformance against them — never the reverse. The build probes
user functions with real values (no reliance on compiler inference), compares
observed against declared, and the same comparison runs on every subsequent
evaluation for free (a `NamedTuple`-type check that constant-folds away when
conformant). The rejected alternative — inference-by-evaluation as schema
authority — fails on three counts, established by walkthrough (§13.4): error
*locality* inverts (failures surface inside correct code, pointed away from the
wrong line); observed schemas are sample- and branch-dependent (the probe sees only
the initial state's branch — the §10 hazard corrupting the schema instead of a
diagnostic); and annotations have nowhere to live. Types by declaration, values by
execution, conformance by comparison.

### 13.2 The declaration inventory

```julia
struct Engine <: AbstractComponent
    ω_idle::Float64; ω_rated::Float64; J::Float64    #parameters: plain struct fields
end

#state stores: declared by initial value — types derived, nothing to drift
init_x(::Engine) = (ω = 0.0,)
init_m(::Engine) = (phase = off,)                    # off | starting | running

#input contract: bare NamedTuple of types at their Float64 faces
inputs(::Engine) = (throttle = Float64, fuel_available = Bool, M_load = Float64)

#output contract = the public interface, as a function of the sweep scalar
outputs(::Engine, ::Type{T}) where {T<:Real} = (M_shaft = T, P = T, ω = T)
#ω names a state field no stage produces → auto-published at stage 1 (§5.2)

function g_s2(eng::Engine, x, m, u, y_s1)
    M_shaft = m.phase === running ? torque_law(eng, u.throttle, x.ω) : zero(x.ω)
    return (; M_shaft, P = M_shaft * x.ω)
end

f(eng::Engine, x, m, y, u, t) = (ω = (y.M_shaft - u.M_load) / eng.J,)

#events: ordered and named — order is load-bearing (§5.2, §11.6); tier is per-event
events(::Engine) = (ignition = Event(ignition_guard, ignition_handler),)
ignition_guard(eng::Engine, x, m, y, u, t) =
    m.phase === starting && x.ω > eng.ω_idle && u.fuel_available
ignition_handler(eng::Engine, x, m, y, u, t) = (x, (; phase = running))
```

The inventory, and where each schema fact gets its authority:

- **State, modes, discrete state, workspace** (`init_x`, `init_m`, `init_z`,
  `init_workspace`): declaration *by initial value* — the type is derived from the
  value, so there is no second artifact to drift and no separate type declaration to
  check. This is the boundary of legitimate derivation: deriving from another
  declaration is sound; deriving from evaluated user code is not.
- **`inputs(::C)`**: a bare `NamedTuple` of types — zero framework vocabulary, no
  wrapper types (the last candidate, `Reduce`, died with reduce-ports, §6). Input
  declarations state the `Float64` face only: their sole job is the build-time
  wiring check, by **exact type equality** against the producer's face. (Names-only
  contracts were rejected — they lose the wiring-time type error and standalone
  checkability; subtype/pattern matching was rejected because its motivating case,
  genuinely generic consumers, is handled by `T`-parametrization below.) Inputs are
  the component's *requirements*: §8's unconnected-input error, over-wiring
  detection and did-you-mean typo messages are only definable against them.
- **`outputs(::C, ::Type{T})`** on continuous components: the declaration is a
  *function of the sweep scalar*, read literally as "this port carries whatever
  scalar this evaluation runs on" — `Float64` nominally, a `Dual` under
  linearization/trim, a tracer under §10. No substitution engine, no sentinel
  types: the author's parameter placement answers exactly the question a rewriter
  would have to guess at (which `Float64`s are the eltype and which are honest),
  and a literal `Float64`/`Bool`/enum inside a `T`-taking declaration is a visible,
  deliberate statement of non-participation. Discrete components declare plain
  `outputs(::C)` — the Tier-3 exemption (§9.2) made visible in a signature. Per
  activity, the build evaluates producer declarations at that activity's `T` to
  type the slots; during a generic sweep, gated-off discrete producers hold their
  `Float64` values, consumers gather mixed tuples, and promotion does the rest —
  semantically exact, since a frozen discrete output is a constant with zero
  partials, which is precisely what "linearize the continuous dynamics with `z`
  held" means. The tier exemption enforces itself through the type system.
- **`events(::C)`**: an ordered, named collection of guard/handler pairs with the
  Tier-2 flag as per-event annotation (§2.1). Order is semantics (§5.2 declaration
  order, §11.6 once-per-event); nothing here is inferrable.
- **No stage tags anywhere.** Which stage produces which port stays invisible in
  the contract, preserving §4.2 (moving a port between stages is non-breaking).
  Membership is *derived* with no chicken-and-egg: `g_s1` functions structurally
  take no inputs, so the build probes them first, observes their contract ports,
  assigns the remainder to stage 2, builds the graph, and probes the `g_s2` chain
  in topological order with real upstream values. The settled "decoder takes no
  inputs" property is exactly what makes the derivation well-founded.
- **Custom structs as port types** (`contact = GearContact{T}`) are first-class
  under the existing §9.2 Tier-1 rules: parametric in their real-scalar leaves,
  constructors inferring `T`, no pinned fields on the continuous path. A pinned
  field wired into the continuous path detonates at the Dual-sweep probe with an
  `InexactError` naming the offending constructor — the §9.2 CI invariant reached
  through the declaration layer with no extra machinery.

### 13.3 Visibility: the contract is the interface

**Declared = public; undeclared = private; absent `outputs()` = no outputs.** Ports
in the contract are connectable, GUI-listed and log-exported. Fields returned by
stage functions but absent from the contract are **private intermediates**: they
occupy table cells (they must survive from their computing stage to `f`/guards, and
they serve the author's own debug panels via the snapshot, §12.2), but they are
non-connectable — a build error, not a discouraged convention — and
presentation-filtered by default. Publicity is never implicit: even the minimal
component writes `outputs(::LowPassFilter, ::Type{T}) where {T} = (x = T,)`, one
line, in exchange for "public" always meaning someone wrote it down.

- **Conformance**: a declared port must be produced — by exactly one stage, or by
  **auto-publication** for declared names matching state/mode (`z`) fields that no
  stage produces (§5.2). Declared-but-unproduced and produced-by-two-stages are
  build errors; a declared name matching neither a stage product nor a state field
  errors with both lists in hand ("not produced by any stage and not a state
  field"). Undeclared returned fields are simply private. The forgotten-branch
  walkthrough survives the flip: a declared `P` missing from the taken branch's
  return fails at probe; missing from an *untaken* branch, it fails loudly at that
  branch's first execution via the always-on check.
- **Branch-shape rule**: stage returns must have the same `NamedTuple` shape on
  every branch — which Julia's type-stability discipline already demands for
  performance; the framework merely makes it a stated rule with a good error.
- **Private cell types are probe-observed** — the one place evaluation retains
  schema authority, accepted because the blast radius is structurally local: a
  private field cannot cross the component boundary, and a divergent branch fails
  at first execution. (A `Private(T)` contract entry keeping type authority total
  is the fallback if migration surfaces private-type surprises; not built.)
- **What died here**: the `unlisted` flag (§4.2 revision note) — presentational
  hiding of connectable ports — and its satellite-function representation; the
  RNG-state case that motivated it needs *nothing* now (`h` reads `z` directly,
  §5.2). The identity-publication default died with it (§9.4 step 4): publication
  driven by the contract replaces publication of everything with hiding
  annotations on top.

### 13.4 Failure walkthroughs (the error-locality grounding)

The four mistakes that decided declaration-vs-inference, with their failure sites
under this layer (each was traced under inference-by-evaluation too; in every case
the failure surfaced inside *correct* code, later, or never):

1. **Typo'd wire** (`:throtle`): build error at the connection, "no input
   `throtle`; did you mean `throttle`?" — vs. a missing-field error inside a
   correct `g_s2` at probe time, with the input set silently *defined* by the typo.
2. **Forgotten wire** (`fuel_available`, read only by a guard): §8 unconnected-input
   error at build — vs. detection contingent on the probe exercising every guard,
   framed as a missing field in event code.
3. **Forgotten branch field** (`P` returned by one branch only): probe or
   first-execution error naming the declared port — vs. a schema silently derived
   from whichever branch the initial state took, then a mid-run error (or a
   silently absent port) at the first transition.
4. **Type mismatch** (a `Float64` fraction wired into a `Bool` input): wiring-time
   error naming both endpoints and both faces — vs. a `MethodError` deep inside
   user math.

### 13.5 Assembly declaration: type-based, kind by declaration shape

An assembly is a plain struct: fields whose type is `<: AbstractComponent` are its
children (field names = path segments), all other fields are inert parameters;
substitutability and variants use ordinary parametric fields — exactly today's
`Cessna172X{K, A}` shape. Alongside it, well-known declarations:
`connections(::A)` (mandatory, even when empty), `exports(::A)`, `rates(::A)`.

**The builder is rejected** (`Assembly()` + `add!`/`connect!`, the early-sketch
spelling): the type you dispatch on and the recipe that defines its structure
become two artifacts with nothing tying them together — §13.1's drift disease at
assembly scale; it threads mutable state through declaration code; and it does not
even buy source-location capture (a called function cannot see its caller's line
any more than a returned value can). Its one real advantage, programmatic
generation, survives intact in the type-based form: a declaration is an ordinary
function body — loops and comprehensions build the returned tuple.

**No `AbstractAssembly`; one root `AbstractComponent`.** Two facts kill a kind
supertype: Julia's single inheritance is already spoken for by the domain
hierarchies (`AbstractAircraft`, engine families — a slot `E <: AbstractEngine`
must accept a primitive `PistonEngine` and a composite turbofan assembly alike),
and §13.3 says kind is an implementation detail behind the contract — encoding it
in public type identity is exactly what contract visibility exists to prevent.
Kind is instead declared by *which* well-known declarations a type defines:
`connections` (the marker, mandatory-even-if-empty — the `LowPassFilter`
precedent) makes an assembly; `outputs` + stage functions make a primitive.
Reading which declarations exist is reading declarations — the same move as the
tier-in-signature split, not §13.1's banned inference-by-evaluation. Error
taxonomy: `connections` plus `init_x`/stages on one type is a build error
(assemblies have no state of their own — §11.5's no-atomic-assemblies at
declaration time); neither marker plus component-typed fields earns a did-you-mean
("holds components but declares no `connections`").

### 13.6 Paths, wiring and exports

**Paths are slash-separated strings** — `"systems/ldg/left/trn"` — relative to the
assembly being declared, no leading slash; one canonical form, shared verbatim by
declarations, error messages, device/trace addressing (§12.3) and the HDF5 log
tree. Rejected: instance navigation (`a.ldg.left` cannot yield a *path* —
symmetric immutable siblings are `===`-identical, so path-from-instance is
unrecoverable by construction; a path-tracking proxy remains addable sugar);
tuples of symbols (structure without readability); dotted paths (a false
Julia-property affordance — the last segment is a contract port, not a field;
slashes say "named tree", which is the true model).

**`connections(::A)`** is an ordered collection of `"src/port" => "dst/port"`
pairs, strictly child-port → child-port; §8's rules apply (one wire per input,
deep paths through concretely-typed fields only, stopping at a generic child's
faces). **`exports(::A)`** is the assembly's face declaration, spelled exactly
like `connections` — an ordered collection of pairs, face name => internal
endpoint path — or a tuple of paths for an input face fanning inward
(`"trn" => ("systems/ldg/left/trn", …)`). **Face names are arbitrary strings with
two build-checked invariants: no `/` (reserved for structural paths) and
uniqueness within the assembly's face set** — every other naming choice
(separators, grouping prefixes like `"pilot.throttle_axis"`) is author
convention, not framework law; the `faces` helper's defaults (§13.8) document the
house style without legislating it. The two-notation rule this rests on: **slash
is structure** (endpoint paths walking real children and ports; snapshot and log
addressing), **face names are opaque contract tokens** — and the periphery
(devices, mappings, trace, GUI write path) speaks face names only (§12.3).
Pairs-of-strings rather than a NamedTuple also removes the `var"..."` noise that
non-identifier names would force. Direction is derived from the endpoints
(wired-to-inputs = input face, wired-from-one-output = output face; mixed or
multi-source entries are build errors). Face *types and tiers* are derived from
the internal endpoints — §13.2's blessed derivation-from-declarations — and the
derivation is forced, not merely convenient: an assembly is tier-neutral (it
exports continuous-sourced and discrete-sourced ports side by side), so neither
leaf `outputs` signature fits, and author-declared face types would need per-face
tier annotations. Rejected spellings: routing values under the leaf names
`inputs`/`outputs` (a name-level pun — a discrete leaf's exact signature with
alien value semantics, killing §13.5's name-level kind split); leaf-style *typed*
faces plus face wires inside `connections` (the tier problem above, plus
face/child namespace collisions and the weakest kind marker); routing-as-wires
with derived types and no face list (facehood implicit in wiring — publicity is
never implicit, §13.3).

**Root slots fall out with no vocabulary**: at every non-root level an exported
input face is fed by the parent's wire; at the root there is no parent, and the
face *is* the write surface (§12.3). §8's whole-tree obligation model states the
complementary error rule. An assembly never declares its external connections —
those live in the parent that instantiates it, exactly as a leaf's do.

### 13.7 Rate scopes

`rates(::A) = (fcs = 1, nav = 5)` — child name => integer multiplier `K ≥ 1`,
optional, unlisted children default to 1; §11.5 semantics unchanged (relative,
composing multiplicatively down the tree, compiled to absolute divisors). Keys are
**immediate child names only** — a deep key would edit another type's design from
outside, and the composition rule guarantees you never need to. `K` on a
continuous child is a build error (§11.5's Δt-on-continuous error at declaration
time). `Δt_base`, `h` and `n` appear in no declaration — they are deployment
decisions fixed at `Simulation` construction. Rejected: `K` carried on the child
instance, FlightCore-`Subsampled`-style — it wraps the field type (polluting
paths, dispatch and the child's contract as seen by wiring) and makes a
per-instance value of what §11.5's own rationale calls intrinsic to the design,
i.e. a fact about the assembly type.

### 13.8 Computed exports and generic boundaries

`exports` is an ordinary function evaluated at build against the concrete
instance, so it may *compute* entries from child contracts — derivation from
declarations, §13.2-blessed. The framework helper, sketched:

```julia
function faces(asm, child_path::AbstractString;
               prefix::AbstractString = child_path,   # "" → no prefixing
               sep::AbstractString = ".",
               except::Tuple = (), only::Tuple = ())  # mutually exclusive

    child = resolve(asm, child_path)      # getfield walk along "/" segments
    names = input_faces(child)            # keys(inputs(c)) for a leaf,
                                          # input entries of exports(c) for an assembly
    unknown = setdiff((except..., only...), names)
    isempty(unknown) || declaration_error(child_path, unknown, names)  # list in hand
    wanted = isempty(only) ? setdiff(names, except) : only
    label(n) = isempty(prefix) ? n : string(prefix, sep, n)
    return Tuple(label(n) => string(child_path, "/", n) for n in wanted)
end

exports(w::World) = (
    faces(w, "aircraft"; except = ("atm", "trn"))...,   # "aircraft.pilot.throttle_axis"
    faces(w, "atmosphere"; prefix = "env", sep = "_")..., # "env_wind_N"
    "view_pose" => "aircraft/pose",                      # explicit entries mix freely
)
```

The child is named by path, never passed as an instance (§13.6's `===` problem);
a face name containing dots is a legal final path segment on the right side
precisely because slash is the only structural separator. `resolve` and
`input_faces` are build-pipeline primitives needed anyway — `faces` is a thin
composition, which is what keeps it sugar rather than machinery; no `rename` hook
because `exports` is ordinary code (map over the pairs). Every error stays
first-class: an `except` face the assembly then fails to wire is an ordinary
unconnected input; a face both wired and re-exported is a two-producers error;
`except`/`only` naming a nonexistent face errors with the child's face list in
hand; a `prefix = ""` collision is caught by the build's uniqueness check like any
hand-written duplicate. The effective export list is plain printable data — the
inspectable contract of this instantiation. What computation does *not* do is
auto-bubble: the author wrote down "every input face of this child that I don't
feed, I expose under this prefix" — explicit at the type level, evaluated at
build.

**Generic holding = imposed contract.** A parent holding a child generically
constrains it exactly through the faces its wires and exports reference: build a
`World` whose concrete aircraft lacks a referenced face and the error names the
`World` entry — build-time structural typing, no new vocabulary (a formal
required-faces declaration on domain abstract types remains possible sugar).
Scalar faces make partial scripting compose: a guidance scenario component wires
`mode_req` and `EAS_ref` while the remaining faces stay exported for GUI or
defaults — impossible with a bundled face (§4.3 write-side corollary).

---

## 14. Case studies

### 14.1 `Vehicle` today → this framework

The grounding exercise that validated §5. Current `Vehicle.f_ode!`
(`aircraftbase.jl:142-170`) is a hand-woven instance of the machinery specified here:

| Today (convention) | This design (checked structure) |
|---|---|
| `kinematics.u .= dynamics.x` — velocity extracted directly from the state vector because `f_ode!(dynamics)` can't run yet | `dyn`'s stage-1 output, scheduled first by construction; the artificial loop in `VehicleDynamics` (velocity state-only, accelerations feedthrough) dissolves |
| Hand-ordered `f_ode!` body (kinematics → airdata → systems → route five `dynamics.u` assignments → dynamics last) | Build-time topological sort; wrong wiring = build error naming the cycle or dangling port |
| Velocity state duplicated (`dynamics.x` and `kinematics.u`) with manual sync, incl. `dynamics.x .= kinematics.u  #essential` in `f_init!` | One state, one owner; consumers wire to `dyn.vel` |
| `get_wr_b`/`get_mp_b`/`get_hr_b` generated tree-walk sums | Summing junctions at ownership boundaries, one explicit wire per contributor, exported totals (§6) |
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

### 14.2 Torture tests for the §5.2 interfaces: `PistonEngine` and the FCS PID cascade

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
  stage-1 port for the previous saturation (`sat_out_0` — a `z` field declared in the
  LQR's output contract, hence auto-published at stage-1 position, §13.3). The delay
  becomes an explicit property of the wiring. (The loop and its fix are
  formalism-independent; the framework's contribution is refusing to let the
  ambiguity through, and stage 1's is having the delayed value already on a port.)

Both components passed without blockers, with zero publications forced beyond current
practice — the empirical basis for §5.2's claim that derivative/output overlap is the
domain norm and the decoder matches the codebase's grain.

### 14.3 Torture test for the §12 staging shapes: filter, joystick and GUI

The exercise that selected per-device cells (§12.3) and produced the §12.5
contracts. Setup (user-level listing: `sketch_io.jl`): a first-order filter with
root inputs `u_cmd` and `τ`; a fictitious 100 Hz single-axis joystick streaming a
slow ramp onto `u_cmd` (complete writer); a 60 Hz GUI with sliders for both slots
(sparse writer); 50 Hz boundaries; pace 1. The interference on `u_cmd` is the
point. The user-level listing came out identical across the three candidate
staging shapes — ergonomics cannot discriminate between them; behavior under a
concrete interleaving did:

- **Drag against the stream** (the user grabs the `u_cmd` slider while the
  joystick streams): under per-slot cells and the batch stack, each frame's
  conflict resolves by last-store/last-push wall-clock order — with 16.7 ms
  renders against 10 ms polls, the applied input alternates between drag value
  and ramp on the cadence beat, the filter visibly wobbles, and the pattern
  differs run to run (the trace replays any given run exactly; the behavior is
  still a timing artifact). Under per-device cells the GUI stages in every drag
  frame (≥ one render per 20 ms frame, plus the active-widget contract), so it
  wins every drain by attachment order: a clean, deterministic override for
  exactly the grab duration. Same user code, qualitatively different physics.
- **Edits while paused**: per-slot cell — the user's `u_cmd` edit is overwritten
  by the still-polling joystick ~10 ms later; the knob visibly snaps back and
  the edit never applies. Batch stack — the edit is buried under newer pushes,
  and the pending chain grows at the polling rate (~10³ nodes per 10 s pause)
  with every peek walking it. Per-device cell — the edit holds in the GUI's own
  cell across the pause (the knob keeps it, §12.5 peek rule), merges with the
  `τ` edit (the sparse-accumulation case), and applies at the un-pause drain —
  for one deterministic frame before the joystick reclaims the slot, which is
  the honest semantics of one-shot editing a streamed input. The uncontested
  `τ` edit works under all three shapes.
- **Corrections the exercise forced**: the sparse-writer lost-write hazard is
  specific to one-cell-per-device layouts (per-slot cells cannot lose
  independent-slot writes — the CAS merge is per-device cells' antidote, not a
  general need); and the batch stack's conflict order is temporal, not an
  attachment-order policy — only per-device cells make precedence a rule rather
  than a race.
- **Discoveries**: the active-widget contract, and the port-resolution answer to
  panel reuse (§12.5) — prompted by asking how the filter's panel survives the
  filter becoming an embedded component with `u_cmd` driven by another component
  (the `Cessna172Xv0` → `Xv1` throttle situation).
- **v0.6/v0.7 note**: under slot exclusivity (§12.3) the contested-`u_cmd`
  scenario itself becomes an attach-time error — the drag-against-the-stream
  phase can no longer arise. The test's verdicts on cell *shapes* (atomicity,
  coalescing, pause behavior, the peek rule) stand; superseded are the
  conflict-precedence findings and the active-widget stage-every-pass contract
  (§12.5) — both correct patches for the contested-slot world this test
  examined, retired with that world.

### 14.4 The interactive C172X demo: the periphery under load (v0.6)

The full-fidelity successor to §14.3, against the real deployment:
`generic_simulation()` (`FlightApps/demos/c172_demos.jl`) —
`SimpleWorld(Cessna172Xv1, SimpleAtmosphere, HorizontalTerrain)`, GUI, joysticks,
an XPlane12 output device, ground/trim init, paced run, post-run plots. Method:
FlightCore's mechanisms are reference *behavior*, not requirements — the question
is whether the new machinery expresses the experience (move stick, plane banks),
never how to reproduce `assign_input!`. Inventory of the complete interactive
surface, with each item's home:

- **Streamed commands** (`throttle_axis`, `elevator/aileron/rudder_axis`): today
  written by joystick mappings after shaping *and* by GUI sliders on the same
  fields. Every dual-writer field in the demo is this pattern — a stream shadowed
  by a mirror, where simultaneous live writing is a bug. This adjudicated slot
  exclusivity (§12.3): claim/disable covers every case found; none needs two
  concurrent writers.
- **Edge-driven increments** (trim offsets ±5e-3 per hat release, flaps ±⅓ per
  button release): today `+=` deltas executed *inside the mappings*, accumulating
  in model `u` — the levels-never-deltas violation, live in the codebase.
  Becomes: devices stage monotonic press counters; accumulator state lives in the
  model (avionics discrete state).
- **The shaping stack**: exp curves and deadzones (defined in the aircraft
  variant module, duplicated *verbatim* across the T16000M and Gladiator
  mappings — the duplication smell), plus the `q_ref = q_sf · axis` fan-out. It
  decomposes into device conditioning (device truth), feel curves (deployment
  preference) and command semantics (FCS design); the face contract splits it —
  conditioning upstream as mapping data, semantics in-model (§12.3).
- **Mode engage** (`mode_req` + setpoint capture from current measurements — the
  GUI handler does `u.EAS_ref = EAS` read from `vehicle.y`): the one place the
  GUI composes writes from model state. Open fork: GUI peek-batch (§12.5 supports
  it as-is) vs. capture-on-engage latched inside the control laws (uniform across
  all writers, but moves behavior into the FCS).
- **Vehicle-direct and environment tunables** (engine start/stop/mixture, payload
  masses, terrain surface enum, sea-level T/p, wind NED): ordinary component
  inputs exported to root faces; GUI as sole unclaimed writer via §12.5; no
  machinery. The interactive surface is *not* one thing: pilot commands cluster
  under a prefix; environment knobs stay with their components' panels.
- **The Xv1 actuator sliders**: FlightCore's dead sliders; resolved read-only by
  §12.5. No action.
- **Outbound** (XPlane12: control-surface angles, nose-wheel steering, prop
  speed/phase, pose, `t`): a snapshot-consuming device, pure `map_output` on the
  device task (§12.2). No friction found.
- **Init/trim, pause/pace, post-run plots**: stopped-sim services (§16), control
  plane (§12.6), log/trace (§12.2–§12.3).

**Architectures examined here and rejected** (the v0.6 periphery decisions were
forced by this cast):

- **Devices as components** (a `T16000M` component wrapping SDL): replay stops
  being same-build (the trace doctrine's strongest property — same type, same
  schedule, staging fed from the recording); the GUI is irreducibly a staging
  device, so inbound uniformity is unreachable anyway (two mechanisms instead of
  one); device lifecycle would duplicate §12.9 in component vocabulary; and the
  drain stops being the single audit point for external data. The salvage: the
  *knowledge* half (a device model's semantics) is expressible as an ordinary
  in-model component wherever wanted; only the wall-clock pump stays outside. In
  an interactive, paced world the scheme is internally consistent (frame-top
  hardware sampling is well-defined) — the rejection rests on the invariants
  above, not on the §12.10 clock criterion.
- **A root-level `PilotInterface` cockpit component** assembling `pilot_commands`
  beside the physical models: its claimed jobs dissolved one by one — struct
  assembly happens in-model downstream of scalar faces (any component can gather
  and bundle), curves became mapping data, widget arbitration is §12.5 +
  exclusivity, and the stateful residue (accumulators, capture-on-engage) fits
  the avionics, where FBW stick shaping arguably belongs. What remained was a
  component with no natural place — a cockpit artifact sitting beside
  aircraft/terrain/atmosphere in `World` misstates the composition.
- **Bundled command faces** (`pilot_inputs` as one struct port): kills per-field
  claiming, liveness and trace provenance — the port is the periphery's atomic
  unit (§4.3 write-side corollary). The routing convenience the bundle bought in
  FlightCore's argument-threading world is provided here by the namespace prefix
  and `faces` (§13.8); the struct reappears legitimately downstream, assembled
  in-model by a single producer.

**Surface walkthrough (v0.7).** The demo line by line, on settled machinery:

- `SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_LOWS15))` —
  pure value construction; no `Model` wrapper (its jobs move into the build).
  `HorizontalTerrain`'s elevation is a plain field (parameter), its surface type
  an input port: the parameter/port split FlightCore kept implicit in
  `U()`-vs-field convention is now the declaration itself. The aircraft's
  `exports` block carries the `pilot.*` face group in one place, deep routes
  spanning avionics *and* systems — today's mapping writes flaps/brakes directly
  into `act`, bypassing avionics; that bypass becomes a declared route.
- `Simulation(world; algorithm = Heun(), h = 0.02, t_end = 1000)` — the entire
  build pipeline runs here: kind resolution, path validation, face derivation
  (computed exports expanded, printable), two-producers/unconnected checks,
  topological sort, probe passes, rate compilation, flat layout, slot table.
- `init!(sim, KinInit(...) | TrimParameters(...))` — stopped-sim services (§16);
  what is settled: they write `(x, m, z)`, **establish every root slot's initial
  value**, and capture the trace header. Slot initialization decisively belongs
  here, not in declarations: the trim service writes slot values it *solved for*
  (throttle, elevator) — not declaration constants.
- `attach!(sim, XPlane12Control(...), binding)` — output device: claims nothing,
  consumes snapshots via §12.8, pure `map_output` on its task. Its binding names
  snapshot paths, **validated at attach against the actual contract** — an
  aircraft substitution that breaks the binding fails at attach, not with silent
  garbage UDP (a new, cheap §12.2 obligation).
- `attach!(sim, joystick, T16000MBinding())` — the binding is a declarative
  table: axis/button → face name + conditioning params
  (`stick_y = (face = "aircraft.pilot.elevator_axis", expo = 1.0, deadzone =
  0.05)`, `button_3 = (face = "aircraft.pilot.flaps_up_count", as = :count)`).
  At attach: faces resolved against the root contract (typo → did-you-mean),
  claim set registered (second joystick on the same faces errors here). The
  Gladiator variant is the same table with different keys, zero shaping code —
  the duplication smell structurally gone.
- `run!(sim; gui = true, pace = 1)` — derived liveness with zero configuration:
  axis mirrors read-only (claimed, with provenance), mode/setpoint/mixture/
  payload/environment widgets live, actuator sliders read-only (component-fed).
  Unplug the joystick → claims release → mirrors go live at the last-held slot
  values. Post-run: `TimeSeries` over retained snapshots; the trace can re-drive
  a fresh `Simulation(world)` bit-identically — which is also the state-trajectory
  inspector (row 38 paying its way).

**Frame anatomies (v0.7).** One frame each, all on settled machinery — no
decision broke:

- *Stick motion*: device task polls, conditioning helper applies binding params,
  complete batch overwrites the cell (inter-frame polls coalesce, ZOH-correct);
  drain applies + traces; avionics tick reads the slot fresh; worst-case
  stick-to-physics latency = poll interval + frame, now by stated semantics.
- *Flaps click*: button peeks counter `k` (own-pending-else-snapshot), stages
  level `k+1` on activation; drain applies; avionics compares slot counter to its
  `z` counter, moves the detent, stores. Multi-click in one window counts via
  own-pending-first peek; repeated staging idempotent (§12.5).
- *Mode engage*: the GUI stages `mode_req` (plus optionally peek-captured
  setpoint slots); **bumpless-engage semantics live in the FCS already** — the
  current `ControlLaws` latches each controller's reference from the present
  command vector on mode transitions, so the fork dissolved: semantic capture is
  aircraft design (status quo, uniform across writers — a script engages sanely
  staging one value); the GUI peek-batch survives as display/slot-sync sugar
  only. Residual check for migration: order-sensitivity of latch vs. sync-write
  on the same boundary (believed none — both derive from the same measurements).
- *Wind slider*: sparse CAS-merge, §14.3's uncontested-`τ` case, live in the
  real cast.
- *Pause/un-pause*: control plane; GUI edits hold in its cell (peek displays),
  joystick cell coalesces bounded; un-pause drain applies both (disjoint slots —
  exclusivity makes the contested question unaskable), pacer re-anchors.
- *Window close*: §12.9 verbatim — complete boundary, final snapshot, sticky
  stopped, wake waits, unblock hooks, named-timeout joins.

Remaining open (feeding §16): the `q_sf` home (thin mapping entry vs.
avionics-internal derivation — aircraft design, not framework design).

---

## 15. Decision log (condensed)

| # | Decision | Rejected alternatives (why) |
|---|---|---|
| 1 | Hybrid causal formalism; two-tier events; projection; no DAE/SDE/per-step hook | DAEs (projection suffices), SDEs (shaping filters suffice), `f_step!` (step-size-dependent semantics) |
| 2 | Causal port-based paradigm | Acausal/MTK (fights interactivity, discrete logic, live introspection); hierarchical callables (rigor by convention — today's footguns); thin SciML library (nothing for GUI/logging/hierarchy to hang onto) |
| 3 | Taxonomy: hybrid continuous primitive + periodic discrete + assemblies; both mode factorings | Strict purity/no modes (loses reset maps; latch logic becomes wiring ceremony); uniform hybrid kind (intra-component ordering semantics murky) |
| 4 | Immutable value signals in a typed signal table | Shared mutable buffers (aliasing/staleness, concurrent-read hazards); mixed semantics (second thing to document/test) |
| 5 | Reject algebraic loops at build, explicit breaks | Implicit delays (silent math changes); per-step numerical solving (jitter, runtime failures) |
| 6 | Strict two-stage structural feedthrough (state decoder `g_s1` + all-inputs `g_s2`); component split as the refinement | N output groups with declared input subsets (declaration/validation surface for a case that never materialized); per-output declarations + tracer verification (under-declaration = silent wrongness); traced multi-pass (hot-path cost, branch unsoundness); component-atomic conservative (false loops everywhere) |
| 7 | ~~Reduce-ports with canonical fold~~ **Reversed in v0.5 → row 37** (explicit summing junctions) | Σ-junctions as default (arity/positional ceremony — objection dissolved by §13's loud declarations); contribution buses (invisible dataflow — verdict unchanged) |
| 8 | Function-valued environment signals + handle pattern | Resource injection (second composition mechanism, invisible); pre-sampling as mechanism (dependency inversion at struts) |
| 9 | Deep paths within owned assembly types only | Unrestricted deep wiring (breaks substitutability across generic boundaries); strict one-level (re-export ceremony) |
| 10 | Structured immutable state over framework-owned flat `Vector{T}` | Mutable views (aliasing, silent missing-ẋ); fully structured/no flat vector (same machinery needed anyway, loses standard integrator interface) |
| 11 | Eltype genericity on the continuous path, three-tier scoping | Float64-only + finite differences (FD noise, keeps hand-written state-space layer, no tracer) |
| 12 | Set-propagation tracers (global + sampled-local), diagnostic-only | Dual-based tracing (derivative-zero blind spot); tracing as scheduling input (soundness) |
| 13 | Immutable `z` in cells + workspace + snapshot idiom | Mutable discrete state (aliasing, snapshot cost); double-buffering (deferred; publication races) |
| 14 | Scoped allocation invariant, CI-enforced on the hot path | Blanket dogma (fights logging reality); no policy (loses the type-instability canary) |
| 15 | Fused evaluation: `f`/`h` read the fresh table (own `y` included); single computation site for derivatives and outputs. **Amended in v0.5 → row 35**: `f`/`h` additionally receive state views; single-computation-site is the rewarded idiom, not an impossibility claim | Mutable caches between `f` and outputs (S-function/FMI style — hidden state, purity violation); accepting duplicate computation (drift bug class: edited law in `f`, stale copy in `g`); derivative binding (superseded — a declaration feature subsumed by `y`-access); orthodox `f(x,u)` + atomization as standard idiom (2× components/wiring for the domain-normal overlap case) |
| 16 | Uniform component interfaces via the `g_s1` state decoder (default identity publication of state and modes); guards and handlers read the fresh boundary table, with per-event re-decode (`handler → project → g_s1 → g_s2`); `project` is the sole raw-state function (schedule-structural); unlisted-port convention for interface noise. **Amended in v0.5 → rows 34/35**: identity default and unlisted retired; `g_s1` redefined as the no-feedthrough stage; per-event re-decode unchanged | Passing state alongside `y` (double-passing; two idioms in the wild); fully private discrete state (breaks uniformity; the codebase culturally publishes state anyway); one-handler-per-component-per-boundary restriction (superseded by the cheap per-event re-decode); exposing z *without* an unlisted convention (RNG/log noise) |
| 17 | Framework-owned simulation loop; stepper seam (advance by arbitrary `h` + on-demand dense output over the last step; one-step methods only); in-house fixed-step RK4/Heun as the sole first-cut backends; `OrdinaryDiffEq` dropped from dependency to possible future extension adapter | `OrdinaryDiffEq` as substrate with `CallbackSet` choreography (semantics by convention in a foreign event loop; demonstrated churn — the `task_local_storage` regression — in exactly the interactive multi-task usage); fused loop without the seam (loses the adaptive/stiff escape hatch for ~zero savings); multistep methods (history rebuild after every handler) |
| 18 | Tier-2 localization: lazy cubic Hermite dense output + bracketed derivative-free root-finding (ITP/Brent) on guard probes that run the sweep; post-event interpolant invalidation + remainder step + bounded event budget | Newton/AD localization (guards C⁰ not C¹ — kinks and σ′ = 0 stretches; discards the bracket certificate for local guarantees; negligible savings on rare microsecond probes); re-integration probes (4× cost; σ becomes trial-h-dependent); solver-matched high-order interpolants (only matter above order 4) |
| 19 | Harmonic tick grid on step boundaries; discrete stages gated to own tick instants (ZOH by construction); assemblies virtual for execution, rate scopes for declaration (integer multipliers `K ≥ 1` composing down the tree, compiled to absolute divisors); `comp.Δt` as single source of truth, no stored `Δt`-derived parameters | Atomic assemblies, incl. opt-in (coarsened schedulable unit → §5.3 artificial loops at assembly scale; interleaving protection meaningless under the signal table; FlightCore's whole-tree atomicity was a call-tree artifact); arbitrary tick periods via time queue (variable `h`, irregular frames, no demonstrated need); absolute-period declaration as default (welds deployment rates into reusable designs; base-period variables don't compose across independently authored assemblies); re-running discrete stages every boundary (un-samples sampled-data semantics); phase offsets (no demonstrated use); `Δt` via `h`-argument only (discretized laws live in `g_s2`) |
| 20 | Boundary event phase iterates to quiescence — rounds of full re-sweep → guards → handlers (declaration order, per-event re-decode) — each event firing at most once per boundary; due `h` updates run after quiescence, outside the iteration | Single pass per boundary (cascade latency N·h — step-size-dependent semantics, the §2.2 `f_step!` footgun class, made common by §3.1 externalized FSMs); bounded-rounds cap (arbitrary K knob; livelock burns the budget then errors instead of degrading to Tier-1 granularity); event/tick fixed-point iteration (structurally unnecessary — `z⁺` is invisible until the next tick decode) |
| 21 | Pacing outside the semantics (bit-identical paced/unpaced trajectories); piecewise-affine wall-clock map, anchor re-established at pace change and un-pause (debt cleared, counted); absolute deadlines with bounded debt + re-anchor on excess; `p = ∞` as explicit pacer-off; hybrid sleep-then-spin toward `deadline − margin`, with `margin` the single knob (0 = pure sleep, ∞ = pure busy-wait = FlightCore) | Relative deadlines (permanent sim-vs-wall slip); unbounded catch-up (burst after long stalls); keeping the anchor across pace changes (retroactively reinterprets elapsed history at the new pace); `p = ∞` as arithmetic limit (perpetual-overrun diagnostics under debt accounting); dedicated busy-wait mode flag (subsumed by `margin = ∞`); separate primitive-resolution threshold (absorbed into `margin` calibration) |
| 22 | Periphery architecture: no shared mutable model — staged inputs drained at frame top + immutable snapshot published per boundary; every handoff one atomic reference op, GC as reclamation; no user code or unbounded work in framework critical sections; control on a separate atomic surface (staging cannot un-pause a drainless loop); interactive = batch + devices | Transplanted `io_lock` (loop budget hostage to arbitrary code under the lock; input timing scheduler-determined and unrecorded — replay undefinable in principle; protects a live-mutation idiom the immutable table removed); full message-passing periphery (per-device typed channels — same design with heavier ceremony) |
| 23 | Snapshot publication: build private → release-store `@atomic latest`; readers acquire-load; wait-free both ways; nothing reachable from a published snapshot ever written again; allocate per boundary; log = retained snapshot references | Preallocated snapshot rings (reintroduce the reader-liveness reclamation proof the GC already provides); `deepcopy` `SavingCallback` logging (the capture *is* the publication mechanism); mid-step publication (§11.3) |
| 24 | Inbound staging: one atomic batch cell per device; complete writers overwrite, sparse writers CAS-merge own cell (retry bounded by drain interception); drain by `atomicswap` in attachment order (conflict precedence a documented policy); levels-never-deltas doctrine; mappings pure, on the device task; device-tagged replayable input trace | Per-slot cells (conflicts by hardware store order — run-to-run behavioral variance; cross-device peek; no trace provenance; atomic-width fallback on wide slots); shared batch stack (temporal conflict order; unbounded pending under pause, taxing peeks); ordered write queue (preserves intra-frame order nothing downstream can observe) |
| 25 | One device kind: uniform handle with read (snapshot / next-boundary) + stage + control capabilities; input-only/output-only as degenerate uses; bidirectional peer = one device; GUI an ordinary device (main-thread affinity and RMW widgets its only peculiarities) | Input/output/GUI taxonomy (lock choreography artifact — blocking rules of `get_data!`/`extract_output` under `io_lock`; forces bidirectional peers into two devices sharing a socket and shutdown); special-cased GUI interface (`sync = 0` + render-under-lock ceremony, obsolete without the lock) |
| 26 | GUI write path: per-component panels name own ports; build-time resolution to root input slots; live vs first-class read-only rendering (with wiring provenance); own-pending-else-snapshot peek; active widgets stage every render pass. **Amended in v0.7 → row 47**: stage-every-pass superseded by stage-on-interaction (its motivating contest died with slot exclusivity) | Slot-naming panels (kills reuse across configurations); always-hot widgets (FlightCore's dead slider — visually live, silently overwritten); cross-device peek (re-couples devices for sub-perceptual benefit); stage-on-change only (streaming device reasserts control mid-grab) |
| 27 | Pacer coarse phase = task-yielding `sleep` (`margin` covers its overshoot); with devices attached every frame yields at least once (explicit `yield()` in unpaced/pure-spin frames); spin never yields; thread budget = sizing rule + startup warning; per-device liveness heartbeat in framework status | `Libc.systemsleep` (second knob inside `margin`; correctness re-hinges on a hard thread requirement; starves co-resident tasks silently — worse failure mode than diagnosed overruns); hard `nthreads` error (the freeze it prevented cannot reproduce: no framework thread monopolist, no stall coupling, GUI on the calling task); yielding spin (µs precision traded for scheduler noise) |
| 28 | Next-snapshot wait: monotonic frame counter + `Threads.Condition`, per-waiter predicate (`counter > last_seen && running`); newest-wins, no queues — outbound coalescing mirrors inbound ZOH; shutdown-interruptible via the predicate | `Event`-based per-frame gate (recurring signal on a latch — the reset has no correct placement under asynchronous waiters; cf. FlightCore's `io_start` reset comments); per-consumer every-boundary queues (unbounded under slow consumers; complete history is the log); polling `latest` on a timer (wasted wakeups, aliasing against the boundary rate) |
| 29 | Input trace on by default, cleared at `init!`, plain kill switch | Opt-in (the trace is primary data — the log is recomputable from it, never the reverse; the session you need replayed is the one you didn't plan to record); tying trace to the log switch (conflates primary and derived recording); rolling window/sampling (complexity without a customer) |
| 30 | Shutdown: complete the boundary → publish final snapshot → sticky stopped status → wake framework waits → `unblock!` hook (close-own-socket idiom; EOT demoted to wire courtesy) → join with named timeout; device crash = `should_close` path; loop failure runs the same protocol from the catch path | EOT as the load-bearing unblock mechanism (protocol detail doing framework work); unbounded join (one wedged device hangs `run!`); mid-frame abort (torn final snapshot; consumers observe un-swept state) |
| 31 | Mid-run mutation doctrine: root-input staging + control commands, nothing else; sim-time scripts = scenario components (clock criterion), wall-clock interaction = devices; `user_callback!` eliminated (cheap composition removed its reason to exist); manual events = slot + guard; init/trim = stopped-sim axis-7 services; mid-run intervention command = guarded addition with shape on record | Scripts as input devices (breaks unpaced — wall-clock staging against µs frames lands at scheduler-determined sim times; both demo archetypes run at `pace = Inf`); retaining `user_callback!` (the periphery's `f_step!`: unrecorded mutation, ordering by convention, invisible to replay); a raw poke API (nothing demonstrated needs it; every mid-run mutation in the codebase is a `u`-write in disguise) |
| 32 | Component declaration: declarative trait layer in plain Julia (well-known functions returning plain values; stage functions ordinary methods); schema authority — declarations define, probe evaluation checks (build probe with real values + free always-on conformance); convenience macros addable a posteriori, never load-bearing | Inference-by-evaluation as schema authority (error locality inverts — failures inside correct code; schemas sample/branch-dependent; annotations homeless); macro DSL as substrate (opaque codegen, tooling/stack-trace tax, only ever lowers to the trait layer); optional declarations with inference fallback (two idioms; the quick hacks most likely to skip are most likely to harbor branch bugs) |
| 33 | Declaration inventory: `init_*` by value (type derived — nothing to drift); `inputs(::C)` bare NamedTuple of types at `Float64` faces, exact-equality wiring check; `outputs(::C, ::Type{T})` on continuous components (functions of the sweep scalar; literal `Float64` = deliberate non-participation), plain `outputs(::C)` on discrete (Tier-3 exemption as signature); `events(::C)` ordered + per-event tier; stage membership derived (inputless `g_s1` probes first, remainder is stage 2), no stage tags | Under-the-hood `Float64→T` substitution (reflection-heavy; cannot distinguish honest `Float64`s); sentinel eltype tokens (same machinery, worse spelling); subtype/pattern matching (motivating case dissolved by `T`; abstract slots break concrete typing); names-only input contracts (lose wiring-time type errors and standalone checkability); per-stage output lists (stage membership is internal, §4.2) |
| 34 | Contract visibility: declared = public; absent `outputs()` = no outputs; undeclared stage-return fields = private intermediates (table cells, non-connectable, snapshot-visible, presentation-filtered); branch-shape-stable returns; private cell types probe-observed (blast radius structurally local) | `unlisted` presentational flag (hidden but connectable — pretends privacy without enforcing it; retired); identity-public on missing `outputs()` (implicit publicity); `Private(T)` contract entries (ceremony without a demonstrated customer — fallback on record) |
| 35 | Stores-and-views prototypes: every component function receives zero-copy views of the stores it genuinely reads — `g_s2(comp, x, m, u, y_s1)`, `f(comp, x, m, y, u, t)`, `h(comp, z, y, u, t)`, guards/handlers alike; the table holds produced signals only, never transported ones (one home per datum); identity default dies; selective auto-publication of declared state/mode fields; `g_s1` = the no-feedthrough stage | State-free evaluation prototypes with identity transport (rows 15/16 as argued: the "published anyway" camouflage fell with contract visibility; the drift-unwritability claim was overstated — `f` always had `u` plus published state); packing `u` into `y_s2` / `f(comp, y, t)` (the reductio: republishing foreign cells under local names); transition-functions-only middle position (fixes handlers but keeps hidden state transport for `f`/`h`); state cells mirroring the buffer (dead stores — no own-function reader remains) |
| 36 | Table mechanics: stage returns are NamedTuples of port values; aggregate `y` = virtual merge, gathered per call, never stored; custom structs are port values — one port, one cell, atomic in wiring; granularity guideline: bundle what shares a stage and is consumed together | Bare-struct returns with field-splatting (ambiguous, type-lossy merge, reflection-hungry); sub-field wiring (the port stops being the atomic unit; field-projection connector kept as guarded addition); per-field cells for struct internals (nested display is a lazy view, not storage) |
| 37 | Aggregation by explicit summing junctions (generic positional or named site-specific — plain components); hierarchical idiom: junctions at ownership boundaries, totals exported across generic boundaries; fold order author-visible; helper/macro sugar guarded | Reduce-ports (row 7 reversed: the declaration vocabulary's last wrapper; three-site census, all Newton–Euler, one library file; canonical-fold, multi-connection legality and identity-element machinery all retired for free; the aggregate wasn't even observable); FlightCore tree walks (silent omission — the zero-edit convenience *is* the hazard); bundled wrench/mass/momentum contribution structs (ragged contributors → identity-element noise) |
| 38 | Snapshot = boundary table (private cells included, presentation-filtered) + `t` + status — no state stores; trace header = full `(x, m, z)` at `init!` (primary data); state trajectory = derived (replay-to-inspect); checkpointing = opt-in log policy, guarded; post-run continuation reads live stores | Per-boundary full-state capture (systematically records derived data — row 29's asymmetry reversed); state wanted in logs via capture rather than declaration (publicity is the honest remedy, priced at one auto-published cell per sweep); dev auto-publish-all-state as default (a diagnostic mode, kin to workspace NaN-poisoning, not semantics) |
| 39 | Assembly declaration is type-based: plain struct, children = component-typed fields, parameters = the rest; `connections` mandatory-even-empty as the kind marker; one root `AbstractComponent`, kind by declaration shape | Builder (`add!`/`connect!`: dispatch type and structure recipe drift apart — §13.1's disease at assembly scale; mutable declaration state; doesn't even capture source locations); `AbstractAssembly` kind supertype (single inheritance is spoken for by domain hierarchies; kind is an implementation detail per §13.3); kind inferred from field types (heuristic where a declaration is wanted) |
| 40 | Slash-string paths, relative to the declaring assembly, one canonical form shared by declarations, diagnostics, devices and logs | Instance navigation (`===`-identical symmetric siblings make path-from-instance unrecoverable — proxies remain sugar); symbol tuples (structure without readability); dotted paths (false Julia-property affordance — the last segment is a port, not a field) |
| 41 | Dedicated `exports(::A)`: face => internal path(s), direction and face types/tiers derived from endpoints (assemblies are tier-neutral — derivation is forced); `connections` strictly child-to-child | Routing values under leaf `inputs`/`outputs` names (name-level pun — discrete-leaf signature with alien value semantics, kills the kind split); leaf-style typed faces + face wires in `connections` (no `outputs` signature fits a tier-neutral assembly; face/child namespace collisions; weakest kind marker); wires-only with implicit facehood (publicity never implicit) |
| 42 | `rates(::A)` optional declaration, immediate children only, `K` on a continuous child = error; `Δt_base`/`h` fixed only at `Simulation` construction | Instance wrappers (`Subsampled`-style: wraps the field type, pollutes paths/dispatch/contract; makes the type-intrinsic ratio a per-instance value); deep rate keys (edit another type's design from outside) |
| 43 | Computed exports as ordinary code + `faces(asm, path; prefix, except, only)`; root slots = the root's exported input faces; generic holding = imposed contract checked per instantiation | Auto-bubbling (forgotten wire silently promoted to a live root slot — §13.4 walkthrough 2 inverted); wildcard-export vocabulary (ordinary code suffices); `add_input!`-style root-slot declarations (second vocabulary for what exports already are) |
| 44 | Slot exclusivity: one writer per root slot at a time; device claims at attach, conflict = attach-time error, release on detach; per-device cells/CAS/drain retained for atomicity and coalescing | Cross-device attachment-order precedence as conflict *policy* (resolves races the §14.4 cast shows nobody wants — every dual writer is a stream shadowed by a mirror); FlightCore-style concurrent multi-device writing of one input (a bug surface, not a feature) |
| 45 | GUI liveness fully derived (transitive root-slot resolution ∧ slot unclaimed); faces carry writer-independent post-conditioning semantics (GUI-parity test); mappings = declarative binding data with per-axis conditioning params, on the device task; edge logic = staged counters + model-state accumulators; unexported ports unpokeable | Per-port "GUI-controlled" markings (the export chain is the marking, owned by the right author); nominally-connected + GUI override channel (second write path; breaks frame purity and trace; done right it collapses into root slots); conditioning in-model (fails GUI-parity — sliders and scripts would be deadzoned); shaping as per-device mapping code (aircraft semantics duplicated per device — today's demonstrated smell); joystick-as-component and root-level `PilotInterface` (§14.4 — replay same-build, single audit point, no natural home in `World`) |
| 46 | Face names = arbitrary strings; build invariants only no-`/` + per-assembly uniqueness; slash = structure, face names = opaque contract tokens, periphery speaks face names only; `exports` returns pairs like `connections`; `faces(asm, path; prefix, sep, except, only)` with dot-prefix *defaults* (convention, not law) | Mandated dot convention (a naming law where two invariants suffice); NamedTuple-returning `exports` (`var"..."` noise for non-identifier names; asymmetric with `connections`); slash-composed prefixes (face names would collide with structural path notation); `rename` hooks in `faces` (`exports` is ordinary code — map over the pairs) |
| 47 | Widgets stage on interaction events only (value widgets on edit, edge widgets on activation with peek-computed counter levels); levels × own-pending-first peek give idempotent repeat and correct multi-click; drain discards stale GUI entries to newly-claimed slots (warning); snapshot includes root slots (peek fallback + read-only mirrors); trace header extends to initial slot values; engage semantics stay in the FCS (the existing `ControlLaws` transition latch — uniform across writers), GUI peek-batch demoted to display-sync sugar | Stage-every-pass (motivating contest died with exclusivity; as insurance it masks invariant violations at render rate — anti-diagnostic; render-rate trace noise); held-button re-staging (auto-repeat at frame rate once the snapshot catches up); capture-on-engage as a GUI/framework obligation (already aircraft design, shipped in `ControlLaws`); slot-initial-values as export-entry defaults (trim writes slot values it *solved for* — services own initialization) |

---

## 16. Open axes

To be settled in subsequent sessions:

- **Axis 7, remainder.** The build pipeline and error-message quality (probe
  orchestration, cycle diagnostics, §13.4's walkthroughs as acceptance tests;
  wiring diagnostics name endpoint paths — the destination path uniquely
  identifies a wire, and a source-location-capturing `@wire` remains addable
  sugar; `resolve`/`input_faces` primitives per §13.8); initialization, trim and
  linearization as first-class stopped-sim services (§12.10), deleting the
  hand-written state-space mapping layer (§9.1) and replacing the per-aircraft
  NLopt trim plumbing; the services own root-slot initialization and the trace
  header's slot capture (§12.3, §14.4). Residual small forks noted in place:
  symmetric `T` on input declarations; `Private(T)` fallback; `q_sf` home
  (§14.4 — aircraft design); sketch refresh (sketches predate v0.5 and the
  v0.6–v0.7 assembly/periphery layers).
- **Migration.** Outline for FlightPhysics/FlightApps (the Tier-1 parametrization
  pass, the `KinData`-style output splits, the contributor survey feeding §6's
  aggregation chains — mechanical to extract from today's trait implementations);
  comparison criteria against FlightCore's demonstrated strengths (zero-alloc
  stepping, flexibility, interactive operation).
