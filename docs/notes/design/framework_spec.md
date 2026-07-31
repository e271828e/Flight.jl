# A Modeling & Simulation Framework for Flight.jl — Specification

This document specifies a modeling and simulation framework intended to
replace `FlightCore` as the substrate for `FlightPhysics` and `FlightApps`.
It is the normative statement of the design: what the framework *is*, in
present tense. It is derived from `framework_design.md`, which remains the
historical original and carries the full record of how the design was
reached.

All design axes are settled — the formalism, the component taxonomy, the
signal and scheduling model, time and execution, the runtime periphery, the
declaration layer, the build pipeline, error discipline and the stopped-sim
services. Only §19's items — the migration outline, the GUI panel authoring
API and the log/trace persistence deferral — remain open.

Decision rationale, including the alternatives considered and the reasons
they were rejected, lives in `framework_decisions.md`, cited throughout as
"row N".

---

## 1. Purpose and method

This document specifies the design of a modeling and simulation framework intended to
replace `FlightCore` as the substrate for `FlightPhysics` and `FlightApps`. The new
framework must match or surpass `FlightCore` in functionality, performance and
flexibility, while being more rigorous and explicit — reducing the learning curve and
the number of latent footguns for model authors.

Ground rules adopted for this design:

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

- **Continuous dynamics**: $\dot{x} = f(x, m, u, t)$ with algebraic outputs.
- **Multi-rate periodic discrete dynamics**: $z^{+} = g(z, u, t)$ at declared rates, with
  outputs held zero-order between ticks.
- **Zero-crossing events**: guard functions with handlers, in two tiers (below).
- **Post-step manifold projection**: an optional per-component hook `x ← project(x)`
  applied after each accepted step (quaternion renormalization, DCM orthonormalization,
  any manifold-valued state). This is the cheap end of the projection-methods family
  from geometric integration.
- **External inputs**: injected asynchronously by the runtime (pilot controls, network),
  under the staging rules settled in §12.

### 2.1 Events: two-tier detection

Both tiers share one declaration (guard function + handler); only the detection policy
differs:

- **Tier 1 (default, cheap):** guards are checked for sign changes at step boundaries
  only. No root-finding, no step rejection; the handler fires at the end of the step in
  which the crossing occurred. Cost: one guard evaluation per event per step. Fully
  compatible with fixed-step real-time execution.
- **Tier 2 (opt-in, per event):** localization of the crossing instant by root-finding,
  for events where timing precision genuinely matters (mechanics in §11.4). Runs
  identically in every execution mode: under real-time pacing the
  localization cost is absorbed as pacer debt like any other expensive frame
  (§11.7) — detection policy never depends on pacing, as the §11.7 invariant
  requires. (Degrading to Tier 1 under pacing was rejected: it would
  move `t*` and diverge paced from unpaced trajectories — row 80.)

This gives step-boundary logic *well-defined semantics*: the transition is defined by
the crossing; detection resolution is an execution-policy detail.

A guard defines a **condition**: a boolean predicate, or the sign of a
continuous function with **positive = condition holds** (normative; writing the
guard's condition value `σ`, holding = `σ ≥ 0`). An event fires when its
condition transitions from
not-holding to holding — edge semantics, uniform across both forms, with the
baseline bookkeeping stated in §11.6; the opposite crossing direction is
declared as a second event with the negated guard (stall entry/exit as a
pair). Tier 1 handles both forms; Tier 2 localization requires the continuous
form. This matters in practice: most transitions in FlightPhysics mix
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
  transitions, stall hysteresis latch) — for one class the mapping tightens
  semantics: level-triggered cross-component resets (the gear friction
  regulator under `!wow`) become edge-triggered events (§17.2, §19).
  Dropping the hook eliminates the footgun of
  model semantics that depend on the integrator's step size.

---

## 3. Component taxonomy

Three kinds, two of them leaves with crisp, closed semantics, one pure composition:

### 3.1 Continuous component (the hybrid primitive)

A classical hybrid automaton:

- **continuous state** `x` (isbits struct of real scalars — see §9),
- **mode variables** `m`: piecewise-constant values (enums, integers, flags) that
  parametrize the flow and change *only* through event handlers,
- **flow** $\dot{x} = f(x, m, u, t)$,
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
- **update** $z^{+} = g(z, u, t)$ at a declared rate,
- **two output stages** (feedthrough applies at update instants: a proportional path
  is direct feedthrough; a state-only output is not).

`z` influences continuous dynamics only **through signals** (outputs held zero-order
between ticks); no component ever reads another's state.

**`m` is continuous-only**: a discrete component has no
mode store — its FSM enums, flags and counters are ordinary `z` fields. `m`
exists on the continuous side because modes must change *between* flow
evaluations through handlers; on the discrete side `g` already runs at the
only instants anything may change, so a second store would duplicate `z`'s
semantics under another name.

### 3.3 Assembly

Pure composition: submodels + connections + exported ports. **No dynamics of its own.**
Hybridness emerges at the assembly level (an aircraft = continuous vehicle parts +
discrete avionics parts). The two-leaf split held under its strongest
counterexample — a strapdown IMU's periodically-reset integrators land on two
leaves with less code than the fused original (§17.5, row 56). Assemblies are flattened away for scheduling but retained as
the navigation/introspection hierarchy (GUI, logging, paths) and as declaration-level
rate scopes (§11.5).

---

## 4. Ports and signals

### 4.1 Immutable value semantics

Ports exchange **immutable values** — typically isbits structs (floats, `SVector`s,
enums, nested immutables). The framework owns a **signal table**: one concretely-typed
**cell** per output port in the flattened model. A producer's output-stage function returns
a named tuple of fresh values; the framework writes each into its cell; consumers read
cells. (*Cell* is the table entry throughout this document; *slot* is reserved for the
root input slots of §12.3.)

Consequences:

- no aliasing, ever — nothing can be mutated under a consumer's feet;
- safe concurrent reads (GUI/logging threads) by construction;
- zero allocation for isbits payloads (named tuples of isbits are isbits);
- each cell has a definite freshness tied to its producer's position in the schedule
  (unlike a monolithic `y` struct, which can be half-fresh mid-sweep with no way to
  tell).

The signal requirement, stated precisely, is **immutability plus frozen references**:
signals may reference bulk data (see §7) provided that data is read-only for the
duration of the run. `isbits` is the common case, not the rule.

### 4.2 Consumers see ports, not stages

The port is the addressable unit. A component's outputs appear to consumers, GUI and
logs as one flat namespace (`dyn.vel`, `dyn.f_c_c`, materializable lazily as a view);
which output stage computes which port is a scheduling annotation, invisible outside
the component. Moving an output between stages is non-breaking *for consumers* —
no wire, log or panel sees it; the scheduler does (the feedthrough graph and
stage membership change, §14.1).

**Visibility.** Which ports exist at all is a declaration-layer decision: the output
contract *is* the public interface, and stage-function results outside it are private
intermediates declared via `local_types` — non-connectable, presentation-filtered, see §13.3. (A presentational
*unlisted* flag — skipped in logs and GUI but still connectable — was rejected: it
pretends privacy without enforcing it, and its motivating case, RNG state feeding the
component's own update, dissolves entirely once the update function reads `z`
directly. Row 16.)

### 4.3 Table mechanics and port granularity

- **Scatter/gather is the whole protocol.** A stage function returns a named tuple;
  the framework scatters each field into that port's concretely-typed cell. Every
  reader — the next stage, `f`/`g`, guards, wired consumers, snapshot capture —
  gathers views from cells. The component's aggregate `y` is the merge of its
  stage products (`merge(y_x, y_xu)`; `merge(y_z, y_zu)` on the discrete tier) —
  declared ports and `local_types` cells alike, since both are stage products —
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
  together*. Bundling across dependency footprints is the `KinData` mistake (§17.1:
  pose is stage 1, velocity-derived quantities are stage 2 — it must split). Fan-out
  is free, so publishing both a bundle and a hot loose field (`pose` *and* `q_eb`)
  is legitimate — one extra isbits cell.
- **Write-side corollary** (from §17.4): **bundle what is written
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
# continuous component — maximal legal view set of each bundle in comments
y_x      = h_x(comp, args)          # x, m, t [, ws] — no-feedthrough stage
y_xu     = h_xu(comp, args)         # x, m, u, y_x, t [, ws]
ẋ        = f(comp, args)            # x, m, y, u, t [, ws]

# discrete component
y_z      = h_z(comp, args)          # z, t, Δt [, ws]
y_zu     = h_zu(comp, args)         # z, u, y_z, t, Δt [, ws]
z⁺       = g(comp, args)            # z, y, u, t, Δt [, ws]

# event system (continuous side only) — same fresh table, same state views:
σ        = guard(comp, args)        # x, m, y, u, t [, ws] — Bool or scalar condition value (§2.1)
(; x, m) = handler(comp, args)      # x, m, y, u, t [, ws] — keys by the return law below (§14.5)
x⁺       = project(comp, x)         # manifold projection; positional (below)
```

**The hand-off (named bundles).** Every function receives exactly two
arguments: the component and one NamedTuple bundle of zero-copy views, from
which the author **destructures by name** only what the body reads —
`f(c::LowPassFilter, (; x, u)) = ...`, `h_zu(c::PID, (; z, u, Δt)) = ...`. The
executor's call is one fixed shape, `fn(comp, args)`; unread fields are ignored
by language semantics; argument order cannot be confused because there is no
order. Naming spellings rejected, for the record: positional
signatures (dead slots written but unread, un-droppable holes mid-list, and the
`t`/`Δt` scalar pair swappable without error); keyword arguments via
`Base.kwarg_decl` reflection (a load-bearing framework seam on a binding Julia
marks internal — the §11.1 `task_local_storage` lesson); keyword arguments with
a `_...` slurp (permanent noise, and "the signature is the read-set" weakens to
"at least") — row 74. `project` alone stays positional — one store in, the same store
out, nothing to select.

**The bundle law.** A name appears in a component's bundle **iff the
corresponding store or fact exists for that component**: `x`/`m`/`z`/`ws` iff
declared (`init_x`, `init_m`, `init_z`, `workspace`), `u` iff the function kind
may see inputs **and** the component declares `input_types`, `y` iff the
component produces any table cell at all (`output_types` ∪ `local_types` ∪
auto-published), `y_x`/`y_z` iff stage-1 ports exist, `t` always, `Δt` on the
discrete tier only. Undeclared stores are *absent*, never `nothing`-filled:
destructuring a field that is not a thing for you fails at the probe inside the
§15.2 framing diagnostic, with did-you-mean against the legal field set ("`f`
of `Foo` destructures `m`, but `Foo` declares no `init_m`") — one law covering
tier facts, stage legality and declarations alike. The mechanism is structured,
not textual: destructuring an absent field throws a `FieldError` carrying the
type and the field name as data (Julia ≥ 1.12), the probe catches it *matched
against the bundle's own NamedTuple type*, and synthesizes the framing
diagnostic from the legal set — classifying the field as an undeclared store, a
wrong-tier fact, or a name illegal for this function kind. No message text is
scraped, and the bundle stays a bare NamedTuple (row 74); a getproperty-wrapper
spelling is the recorded fallback should type-matched interception prove
insufficient. The per-function name sets
are **closed**: adding one is a decision-log entry, not a convenience. The
comment in the signature block above states each function's maximal legal set;
a given component's bundle narrows it to declared reality, and the
destructuring narrows further to actual reads — a three-level funnel (stage
name ⊇ bundle ⊇ reads) worth teaching once, because a stateless component
legitimately writes `h_xu` while owning neither `x` nor `m`.

**The handler return law.** The same rule governs the return side. A handler
returns a NamedTuple carrying the stores it writes: a key is present **iff**
the corresponding store exists on the component **and** the handler updates
it. A pure FSM (modes and events, no `x`) returns `(; m = (; phase =
running))`; an `x`-only reset map returns `(; x = (; x..., ω = 0.0))`; a
handler touching both returns both. Padding forms — `((;), m⁺)`, `(x⁺, (;))`
— do not exist, for the reason row 74 gives on the argument side: a filler
value defers the error from the return to first use, with a worse message,
and a positional pair additionally lets the two stores be swapped without a
name-shaped diagnostic. Semantics per key: `x` present ⇒ the value is
complete against the state field set; `m` present ⇒ the names-subset
predicate; an unknown key ⇒ did-you-mean against `{x, m}` — the same
`FieldError`-shaped machinery §15.2 builds for bundles, now running in both
directions (§14.5).

The views themselves are unchanged in meaning: own state (`x`, `m` from the
flat buffer and mode cells; `z` from its cell), own published signals (`y`,
gathered from own table cells — `local_types` cells included alongside declared
ports, locals existing precisely for this own-later-consumption read, §13.2),
inputs (`u`, gathered from foreign cells
through the wiring's name binding), the clock (`t`, and `Δt` — see §11.5), and
scratch (`ws`, §9.3). The signal table holds only *produced* signals, never
transported ones: each datum has exactly one home — buffer for `x`, cells for
`z`/`m`, table for signals — and no store mirrors another. Every bundle field
earns its place as a view genuinely readable, and no further "simplification"
exists that does not introduce a copy: eliminating `u` would mean republishing
foreign cells under local names; eliminating `x` would mean identity transport
through the table, rejected for exactly that reason (§9.4, step 4).

**The letters**: `f` is the continuous flow, `g` the discrete
update, `h_*` the output stages — the hybrid-systems flow/jump pair
(Goebel–Sanfelice–Teel) joined to the control/estimation convention that `h` is
the output map (`y = h(x, u)`; every navigation filter's measurement function).
Bare `h` now means the integration step size only (§11), retiring a double
booking. Stage suffixes name the **dependence class**, not the argument
list: `x`/`xu` continuous, `z`/`zu` discrete — state-only versus
state-plus-input, the `y = h(x)` / `y = h(x, u)` distinction spelled in the
name, with the two tiers mirroring each other exactly. So "no `u` in the
name" *is* the no-feedthrough property, visible at every definition site.
The letters are deliberately non-exhaustive: modes fold under the state
letter (`m` is state — the objection that `h_x` omits it is answered by the
suffix's job being the feedthrough split, not an inventory; an earlier
`h_xm`/`h_xmu` spelling was traded away for the tier symmetry and the
textbook mirror), and ambient facts (`t`, `Δt`) and scratch (`ws`) ride
unnamed. A stage name on the wrong tier (`h_zu` on a continuous component)
is a build error. Rejected namings: verb names (`decode`/`compute` — encode
nothing), `Moore`/`Mealy` (exact and opaque), a tier-neutral pair (no honest
neutral state letter), `g`-for-output (forfeits both the jump-map alignment
and the step-size disambiguation).

- **`h_x`/`h_z` is the no-feedthrough stage** — defined entirely by what it
  cannot see: its bundle carries no `u`, so "no feedthrough" is unfalsifiable,
  and that structural guarantee is what its ports contribute to the schedule
  (they break would-be loops). It exists when the component has state-derived
  ports or shared state-derived intermediates for stage 2; otherwise it is
  simply absent. **A declared output that matches a state, mode or `z` field by
  name and type, and that no stage produces, is auto-published** by the framework
  from the state stores at stage-1 position (§13.3) — on the discrete tier the
  match is against `init_z` and the publication position is `h_z` — publication driven by the
  public contract, rather than the blanket identity publication of state that
  row 16 rejected.
- **`h_xu`/`h_zu` receives all wired inputs plus `y_x`/`y_z`** — its own
  stage-1 results, so shared intermediates are computed once, not re-derived —
  plus the state views; conservatively, every stage-2 output is presumed
  dependent on every wired input.
- **`f` and `g` run after the sweep**, when the full signal table — including the
  component's own stage-2 ports — is complete and fresh. The fused idiom stands:
  compute each law once, in a stage; publish it; let `f`/`g` copy from `y`. The
  interfaces *reward* single-source-of-truth (nothing ever needs computing twice)
  rather than claiming to make duplication unwritable — a claim no prototype
  could honestly make anyway, since `f` always had `u` and the published state
  (row 15).
- **Guards and handlers read the same fresh world.** At a step boundary, the order
  is *integrate → project → boundary sweep → guards*, so by guard/handler time `y`
  is a fresh decode of exactly the state being transformed, and the state views are
  that state itself. Handlers construct their `x`/`m` returns from raw state
  naturally — a reset map is `(; x = (; x..., ω = 0.0))`, no reassembly from
  published fields.
- **`project` runs between a state write and its decode** (after integration, and
  after any handler `x`-reset) — the only positions in the schedule where no fresh
  `y` of the new state can exist yet. It is not *unique* in receiving raw state
  (every function gets state views), but it is unique in that schedule position.
- All output stages must be pure (no side effects); state types make mutation
  impossible anyway (§9).

The schedule: all stage-1 functions (any order), then stage 2 in topological order, then all `f`
against the now-consistent signal table. Note the systemic consequence: *evaluating
the RHS means running the sweep* — there is no incremental `f`-only re-evaluation.
Implicit solvers, linearization and trim already work this way (seed `x`, run the
composite), so nothing is lost; §11.3/§11.4 restate it as a property of the
execution model (RHS evaluations and guard probes alike run the sweep).

**Step-boundary semantics.** At each boundary: integrate → project → boundary sweep →
evaluate **all guards once** against that sweep → for each fired event, in declaration
order: `handler → project → re-run the component's output stages`. The per-event
re-decode keeps `y` fully fresh for any subsequent handler of the *same* component
(sequential composition, no lost updates) and leaves the signal table
post-transition-consistent for whatever else the boundary does (discrete ticks,
logging). Events are rare and the re-run touches one component, so the cost is noise.
Across components, all guards and handlers read the same boundary snapshot (plus their
own component's refreshed ports); one component's transition reaches others through
the next sweep (delivered by binding time — handler-phase `u` gathers materialize at
round start, `y` binds live; mechanism and order in §11.6). Newly-enabled guards fire within the *same* boundary: the
sweep → guards → handlers phase iterates to quiescence, with each event firing at
most once per boundary (settled in §11.6).

**Departure from the orthodox formalism, stated openly.** The textbook form is
$\dot{x} = f(x, u)$, $y = g(x, u)$; this design's `f` receives the orthodox arguments
*plus* the published table: $\dot{x} = f(x, m, y, u, t)$. The composite map $x \mapsto \dot{x}$ is
mathematically identical (linearization, trim and AD are untouched); the heterodox
element is only that derivatives may read outputs. The teaching line: *"stage 1
publishes what you know from state alone; stage 2 adds what needs inputs; your
dynamics read your own published results instead of recomputing them."* The decision
was grounded in a component-by-component survey of FlightPhysics/FlightApps (§17.2):
derivative/output overlap is the *norm* in this domain (Newton–Euler, kinematics,
piston engine, gear friction, every discrete compensator), so the orthodox split
would force either systematic duplicated math (with its silent-drift bug class),
systematic component atomization, or a cache mechanism — all rejected. FlightCore's
fused `f_ode!` already embodied the same economics; this design keeps them while
adding checked scheduling.

**Shared expensive computations** are thereby solved uniformly: compute once in
stage 2, publish, and let `f`/`g` (and external consumers, e.g. an accelerometer model
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
diagnostic says so explicitly ("cycle through `systems/aero` is artificial at port
level — split the component", with the offending stage `h_xu` carried as a separate
payload field rather than dotted onto the path, §13.6/§15.2). One consequence of stage-2 conservatism worth recording: an
input consumed only by `f` (never by `h_xu`) still creates a scheduling edge if the
component has stage-2 outputs; in practice such components are integrator-shaped and
have none, and the remedy, if ever needed, is the same split.

### 5.4 Algebraic loop policy: reject at build time

A genuine cycle in the instantaneous dependency graph is a **build error** with a
diagnostic naming the full path in the canonical slash form of §13.6
(`aero/F → dyn/a → aero/α̇ → aero/F`). The user breaks
it explicitly: insert dynamics (the α-filter idiom — already standard practice in the
domain and in the current C172 model), insert an explicit unit delay (§15.7's
`UnitDelay` — note that this remedy changes the model's tier structure: the broken
signal becomes discrete, sampled at `Δt_base`, which is a modeling decision, not a
transparent wire), or restructure.

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
struct SumJunction{W, N} end        #type constructor, arity; library-provided

input_types(::SumJunction{W, N}) where {W, N} =
    NamedTuple{ntuple(i -> Symbol(:in, i), N)}(ntuple(_ -> W{Float64}, N))
output_types(::SumJunction{W, N}) where {W, N} = (; Σ = W{Float64})
h_xu(::SumJunction, (; u)) = (; Σ = +(u...))
```

(The parameter is the *unparametrized* type constructor — `SumJunction{Wrench, 3}`;
UnionAlls are legal type parameters — so both declarations derive from it at the
`Float64` nominal face; the activation leaf walk (§13.2) re-scalars the output
cell per activation. This is the same
arity-via-computed-contracts pattern §15.7 commits to for `Or{N}`.)

Wired at an ownership boundary, the junction is ordinary structure — with
`wr_sum::SumJunction{Wrench, 3}` a field of `Systems` like any other child
(§13.5):

```julia
connections(::Systems) = (
    "aero/wrench" => "wr_sum/in1",
    "pwp/wrench"  => "wr_sum/in2",
    "ldg/wrench"  => "wr_sum/in3",
    "wr_sum/Σ"    => "dynamics/wr_ext",
)
```

- **Every mistake is loud** under the declaration layer: a forgotten contributor is
  an unconnected-input error naming `in4`; a double-wired slot violates
  single-connection; a stale arity surfaces as one or the other. The bookkeeping is
  ceremony, never silence.
- **The aggregate is a first-class signal.** `wr_sum.Σ` is an ordinary port:
  loggable, GUI-visible, fanned out to a second consumer (a loads monitor) for one
  wire. (Under consumer-side folding, the total was ephemeral gather scratch.)
- **Aggregation logic is arbitrary stage-2 code** — mass-properties composition with
  its transport terms, weighted blends — not restricted to a declared
  commutative-associative binary op.
- **Fold order is author-visible**: the positional order of the junction's inputs.
  Reassigning contributors to different slots changes summation order, hence bits
  (float non-associativity) — deterministic per configuration and under author
  control, which is strictly more explicit than a framework-canonical order.
- For the handful of real sites, a **named site-specific junction**
  (`input_types(::VehicleWrenchSum) = (aero = …, ldg = …, pwp = …)`) documents the
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
guarded-addition sugar, added when migration shows the pattern repeated. The
junctions themselves — summing junctions, Bool gates — are the seed of the
standard component library committed in §15.7: ordinary components, no
framework privileges, inventory grown strictly by migration demand.

The ledger against FlightCore's tree walk, recorded: its zero-wiring convenience and
its worst failure mode were the same property — a contributor with a forgotten trait
method contributes *silently nothing* (a lighter vehicle, no diagnostic, ever).
Explicit wiring inverts every silence into a warning or error, and makes
per-contributor values and intermediate totals observable ports; the cost is that
adding a deep contributor edits one assembly level (its owner's wiring) instead of
zero. Rejected at every revision: **contribution buses** (today's mechanism
portified — dataflow invisible in the graph, scoping rules, accidental
contributions).

An earlier design specified **reduce-ports** — consumer-declared
commutative-associative folds with multi-connection legality, canonical fold order
and identity-element opt-outs — and it was abandoned (rows 7 and 37): the use-site
census never grew beyond the three Newton–Euler aggregations in one library
assembly; the mandatory typed declarations of §13 make junction mistakes loud,
dissolving the original "positional ceremony" objection; and `Reduce` was the
declaration vocabulary's last remaining wrapper — killing it leaves `input_types()` as bare
types with zero framework vocabulary and retires the canonical-fold and
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
  inside their own stage functions.
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
redesign. Substitutability behind a stable face is declared with an abstract
input entry (`terrain = AbstractTerrainField` — §13.2's structural
substitutability): the consumer wires to any concrete field type below the
bound, preserving today's `AbstractTerrain` polymorphism at the declaration
layer.

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
  hard-codes one implementation and breaks on substitution). Enforcement lives
  in the path-resolution primitive itself (`resolve`, §15.3), which walks
  declared field types alongside instances.
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
- Unconnected output ports: legal, silently (there is no build-time warning
  for them, row 84 — under mandatory `output_types` most models carry many
  observation-oriented ports no wire consumes, §12.2 blesses exactly that,
  and §12.2's path-bound readers attach after the build, so "unused" is not
  even decidable; a warning firing on every honest port poisons the sole
  warning stream, §12.5's anti-diagnostic lesson. The hazard it nominally
  guarded — a wire someone meant to draw — is caught from the consumer side,
  where the information actually lives). Unconnected input ports: build error
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

Each continuous component declares its state by value (`init_x`, §13.2): a
NamedTuple whose leaves are drawn from a **deliberately closed vocabulary —
plain real scalars and `SArray`s (static vectors and matrices) of a common
eltype `T`** — and nothing else. `Int`s/enums/`Bool`s belong in modes, and
domain wrapper types (`RQuat`, `Ranged`) are not state leaves — an attitude
state is an `SVector{4,T}`, cast where rotation semantics are wanted (below).
The framework:

- computes a **flat layout** at build time (compile-time offsets over one contiguous
  `Vector{T}` buffer it owns);
- **reconstructs** the typed immutable state value for a component at each evaluation
  and passes it to every function receiving state views (§5.2 argument rule — field
  loads at known offsets, register-level, zero cost);
- receives immutable results back: derivative functions return an `Ẋ`-typed value
  (scatter-stored into the flat `ẋ` buffer); event handlers and projection return a new
  `X` (written back).

**What `Ẋ` is.** With the leaf vocabulary closed, the answer takes one line:
`Ẋ` has exactly `X`'s shape at the activation scalar — a scalar leaf's
derivative is a `T`, an `SArray` leaf's is the same `SArray` at `T`. (This is
the vocabulary rule paying rent: an invariant-carrying leaf like a unit
quaternion has a derivative off its own type, and `Ẋ` would need a separate
derivation; here the attitude leaf is an `SVector{4,T}` and so is its rate.)
The conformance predicate is structural — *each field of `f`'s return
scatters into its field's block at `T`* (§14.5 states the check) — which
makes derivative completeness a property of the layout rather than of author
discipline. There is deliberately **no `derivative_type` hook**: a per-leaf
override would be a second register for a fact the layout already knows, and
the two could disagree.

**The buffer is authoritative; typed values are ephemeral reconstructions.** Nobody
outside the framework ever holds a mutable reference to state. "Ephemeral" is
literal: an isbits view materializes in the caller's frame (registers or spilled
stack, the compiler's business) for exactly the duration of the call and has no
existence between calls — re-materializing is the same loads, value-identical
because the value is immutable and the buffer unchanged within a sweep. Whether
repeated reads within a sweep re-materialize or reuse the loads is codegen freedom
in the literal sense: the executor is spelled rebuild-per-call, and hoisting is the
code generator's CSE, whose legality condition is exactly the
buffer-unchanged-within-a-sweep rule (§14.7). The complementary rule (§5.2): **one home per datum** —
buffer for `x`, cells for `z`/`m`, table for produced signals — and no store ever
mirrors another; in particular there are no state cells in the table beyond
contract-driven auto-published ports, which are interface, not transport.

**Why the vocabulary is closed: views must materialize without running
anyone's invariants.** Scalars and `SArray`s have invariant-free
constructors — `SVector`'s stores its tuple, `NamedTuple` construction runs
no user code, nothing normalizes or clamps — so building a view through
ordinary public construction is bit-faithful automatically:
`reconstruct(flatten(x)) == x` identically, with no constructor bypass, no
`reinterpret`, and no reliance on a custom struct's memory layout mirroring
the buffer's. Admitting invariant-carrying leaves would force one of two
readings, both rejected (row 94): run the constructor on read, and every
consumer sees a silently projected value over a buffer accumulating the raw
one — a `Ranged` leaf integrating past its bound while the clamped view hides
it and the derivative keeps pushing, `project` (§2) redundant on views yet
necessary on the buffer, §11.4's off-manifold probes (which must see the raw
interpolated state, RK-stage-like) impossible; or bypass the constructor,
which is safe under the common-eltype rule and build-time verifiable, but
plants layout-coincidence cleverness in the executor's core to save a cast
line. Domain semantics are instead an **explicit, invariant-free cast at the
point of use** — `q = RQuat(x.q, normalization = false)` — the conversion
today's `f_ode!` code performs on its raw views, now visible and chosen.
Invariants live where the design already put them: in `project` at
boundaries, and in writers — handlers build their returned values through
ordinary constructors, and the condition apply converts authored values
through ordinary `convert` methods (§16.3). Constructors run on the write
paths, never on views.

What this buys, against today's flat-`Vector` + `ComponentArrays`-views pattern:

- no aliased mutable views (who-writes-what by convention);
- derivative completeness is **structural** — the returned `Ẋ` has every field by
  construction; a forgotten `ẋ` entry is impossible rather than silently stale;
- state fields arrive as the declared scalars and `SArray`s, immutable; the domain
  wrapper, where wanted, is one explicit invariant-free cast
  (`RQuat(x.q, normalization = false)`) — the conversion the mutable-views
  pattern performed implicitly, now visible and chosen;
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

The declaration layer keeps this scoping out of the author's way: declarations
are concrete nominal types (§13.2), and cell types per activation come from the
framework's leaf walk over them — `Float64` leaves and `Real` type parameters
follow the activation scalar on the continuous side, the discrete side pins —
never from inference through user code. (Safety of the substitution rests on
the §14.5 embedding guarantee.)

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

- **Storage:** `z` and `m` live in **typed cells** (same discipline as signal cells),
  overwritten by the framework when an update/handler returns a new value. They never
  touch the integrator buffer; no arithmetic is ever done on them.
- **Type freedom:** `z` may be *any immutable value* (frozen-reference rule; isbits not
  required). Enums, integers, nested structs, RNG state (four `UInt64`s of `Xoshiro` —
  required to live in `z` for deterministic replay).
- **Snapshots are free:** copying a cell copies a reference to immutable data —
  checkpoint/replay of the entire discrete side is "copy the cell values."
- **Workspace** (for heavy algorithms, e.g. an n≈20 Kalman filter):
  component-declared mutable scratch, instantiated by the framework and arriving
  as the `ws` bundle field (§5.2) in every bundle-receiving function of the
  declaring component (`project` is positional and receives none),
  and **excluded from state semantics** — not snapshotted, not replayed, never a
  condition target (§16.1), must carry no information between calls. Checkable:
  in debug mode the framework **poisons the workspace before every call**, so
  read-before-write of stale scratch detonates immediately. The poison is scoped
  to what a store can carry: float-eltype stores get `NaN`, integer-eltype
  stores get `typemin` (which detonates on first use as an index), and element
  types with no sentinel — `Bool` flags, enums, opaque handles — are skipped,
  the skip reported once per activation. Naming the skipped stores is the point:
  a guarantee that is silently partial is exactly the anti-diagnostic pattern
  §12.5 rejects. Declared **by
  allocation**: the well-known method *is*
  the allocator — `workspace(c::KF, ::Type{T}) where {T} = (P = Matrix{T}(undef,
  c.n, c.n), x̂ = Vector{T}(undef, c.n))` on the continuous tier, plain
  `workspace(::C)` on the discrete (the `output_types` tier split) — called
  once per activation and once per scratch-store set (§16.8), sizes from the
  instance, eltypes from the activation. The `T`-signature survives here after
  row 79 retired it on the output side because the two registers ask different
  things of it: `workspace` is the by-allocation register, so the method is an
  allocator the framework *calls*, not a schema the framework *walks*. Nothing
  downstream derives from its type, and mistyped scratch detonates loudly at the
  `Dual` probe — the opposite of the silent zero-partial embedding that killed
  the output-side `T`-signature. The `undef` spelling is the
  recommended idiom: it makes "contents are meaningless" visible in the
  declaration, which is the register this store actually lives in — a
  workspace is *not* memory, so declaration-by-initial-value (§13.2) never
  legitimately covered it, and none of the §13.2 arguments for authored values
  (condition overlay base, the probe-value barrier) applies to a store that
  conditions exclude and the poison overwrites. Available on **both tiers**:
  nothing in the contract is tier-specific, and a continuous workspace simply
  joins the `T`-generic surface — under a `Dual` activation the allocator is
  called at `Dual`, and the in-place math runs through Julia's generic
  fallbacks (no BLAS; activations probe and linearize, they don't run
  marathons). The calls-per-boundary multiplicity of the continuous side (RK
  stages, Tier-2 probes, event re-sweeps) makes the no-information-between-
  calls contract *more* load-bearing there, not less — which is the poison
  check earning its keep, not a prohibition. (Rejected: `workspace_type`
  returning types with framework instantiation — array types carry no
  dimensions, sizes live in runtime fields like `kf.n` which no zero-argument
  constructor can see, and hoisting sizes into type parameters lands on the
  `MMatrix` codegen catastrophe below.)
- **Blessed idiom — zero-allocation ticks with immutable `z`:** do the in-place math
  (`mul!`, `cholesky!`, BLAS) on the workspace; at the end, snapshot into an isbits
  container (`z = KFState(SVector{20}(ws.x̂), SMatrix{20,20}(ws.P))`) and return it.
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
   the fresh signal table to `f`/`g` subsumes it: the "binding" is a one-line function
   body, no validation machinery, strictly more general (an `f` may combine published
   values with extra terms).
3. **Separate state arguments → the state decoder.** With `y` in hand, passing `x` to
   `f` duplicated information (FlightPhysics culturally publishes state in `y` anyway);
   removing it produced the uniform continuous/discrete shapes and the
   single-computation-site guarantee.
4. **Decoder-exclusive state access → stores-and-views arguments.**
   Step 3's second half was reversed (row 35). Its justification — "state is published
   anyway" — quietly depended on default identity publication; once §13.3 made
   publication a deliberate interface act, the identity decode stood revealed as
   *transport*: copying the buffer into cells so a buffer view could be replaced by
   a cell view — ceremony of exactly the kind step 2 removed elsewhere, now with
   dead stores (no own-function reads state through the table). The
   reductio that decided it: the same minimalist logic would pack `u` into the stage-2 product
   and end at `f(comp, y, t)`, republishing foreign cells under local names. The
   fixed point is §5.2's argument rule — zero-copy views of the stores a function
   genuinely reads. What survives of step 3: the uniform shapes, the fused
   economics, and the stage-1 decoder itself (today's `h_x`/`h_z`) — no longer the sole state gate, but the
   no-feedthrough stage.

Prior art, for orientation — every causal framework meets the shared-computation
problem and resolves it per its architecture: **Simulink diagrams** make integrators
explicit blocks (derivatives are ordinary wires into `1/s` — the computer/integrator
split is their native idiom); **S-functions and FMUs** use sanctioned *mutable caches*
(DWork vectors; FMI's lazy-evaluation caching) between their `mdlDerivatives`/
`mdlOutputs`-style callback pairs; **Modelica/MTK** write `der(x) = expr` natively
with symbolic CSE. The fused sweep + signal-consuming `f`/`g` is the cache-free
formulation that fits this design's purity rules — and it is also what FlightCore's
fused `f_ode!` did economically, minus the checked scheduling.

The **computer/integrator split** remains fully expressible without any framework
support (a stateless component computing derivatives as outputs, wired into a trivial
state-holding component) and is the idiom of choice when the factoring earns reuse —
e.g. one Newton–Euler solver shared across vehicle variants, or swappable kinematic
descriptors against a common integrator shape. See `sketch_decoder.jl` for the
worked example; against a split-form spelling of the same model (four components,
thirteen connections), the merged form has half the components and wiring, and everything
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
- **Logging**: amortized-zero — snapshots are records stored *inline* in a
  `Vector`, `sizehint!` for the expected duration making regrowth a non-event.
  The inline-storage claim is about the snapshot record's *fields*, not about
  everything reachable from them: a model carrying §7 field handles (heightmap
  terrain, wind grids) has a snapshot type with reference fields, and those
  ride as references to build-time-frozen data — no copy, no per-boundary
  garbage, which is what the allocation claim asserts. What the claim does not
  assert is that the snapshot is `isbits`; the per-boundary allocation cost is
  zero either way, and §7's summarize-or-skip rule governs what such a field
  contributes on export.
- **What is not recorded**: event firings. The log holds boundary snapshots and
  the trace holds staged inputs (§12.2, §12.3); neither carries a per-event
  record. Which events fired at which boundary is recovered by replay plus the
  published modes — the honest remedy of §12.2, a mode field declared public is
  in every snapshot. An event-firing stream is a guarded addition, not built.
- **Tools where garbage is unavoidable**: arena allocation (Bumper.jl-style) for scoped
  temporaries; scheduled `GC.gc(false)` at frame boundaries to move collection out of
  the critical path. (Julia has no per-object freeing; these are the honest levers.)

---

## 10. Diagnostics: feedthrough tracing

Tracing is **diagnostic only, never load-bearing**: scheduling correctness comes
exclusively from the structural two-stage split; tracing improves error messages and
verification. Triggered when the scheduler finds a cycle, to classify it (genuine →
"insert a state"; artificial → "split this component").

**Detection and naming.** A cycle surfaces as a topological-sort stall in
Stratum B. The stall's residue is not the diagnostic: it also holds the innocent
downstream cone, so it over-reports, while a single DFS back edge under-reports —
one edge of a possibly large tangle. The stalled subgraph is instead decomposed
into **strongly connected components**, and each nontrivial SCC names one cyclic
cluster exactly: one diagnostic, its members and the wires among them, presented
as one readable loop in §13.6's canonical slash form
(`aero/F → dyn/a → aero/α̇ → aero/F`).

**Classification is schedule-free.** It runs inside Stratum B's failure path,
where no schedule exists — and needs none, because each SCC member is evaluated
*once, in isolation*, at the probe point: state views from `init_*`,
out-of-cycle cells from the acyclic prefix's probe values, in-cycle cells
synthesized through §14.3's `probe_value` under tracer tags. The tracer's
product is a per-member dependence set rather than a value, so no ordering has
to be valid for the labels to come out right. The loop is **real** iff every hop
of the structural cycle survives in the traced per-member maps; **artificial**
(§5.3) iff some hop dies — the component whose stage-2 function does not in fact
route that input to that output. No Stratum C machinery is touched: no
activation, no layouts, no table. This is row 12's *local* variant — the
schedule-free per-member trace at the probe point, which is what the cycle
classifier uses; §14.4's "tracer activation" names row 12's other variant, the
global set-tracer run as an ordinary Stratum-C activation, and the two must not
be conflated.

**Caveats, carried in the diagnostic rather than assumed away.** The trace
speaks for the branch taken at the probe state (row 12's diagnostic-only
doctrine). Discrete members trace *structurally* — row 79's pinned signatures
admit no tracer scalar — which is sound as a may-depend answer but never sharp,
so the "split this component" hint is offered only for continuous members. And
if a member's evaluation itself throws, the diagnostic ships with the member
list alone: classification is a bonus on the cycle error, never its
precondition.

Two modes, degrading gracefully:

- **Global (value-blind) set-tracer** — a `Real` subtype carrying a set of input
  indices, unioned by every operation. **May-depend semantics**: sound
  over-approximation (saturated `clamp` still reports; $u^2$ at $u = 0$ still reports).
  Exact, one evaluation. Requires the traced stage-2 code to be free of
  branches/lookups on *input-tainted* values.
- **Local (primal-carrying) set-tracer at sampled states** — the fallback whenever the
  global tracer hits an undecidable branch (piecewise friction, stall blending, any
  gridded lookup at an input-tainted coordinate). Reports the dependence pattern of the
  taken paths, sampled across randomized states. Strictly dominates Dual-based tracing
  (no derivative-zero blind spot; only untaken-branch misses).

Boundaries: only *inputs* are seeded, so branching on state/modes/parameters/time never
interferes (stage-2 functions also receive state views, but neither those nor `y_x`
are ever seeded). Stage-1 functions are never traced (nothing to
seed). Derivatives, guards, handlers, projections are outside tracing's jurisdiction
entirely. Known tracer blind
spot: value-severing operations (dependence passed through a bare `Int` index, e.g.
nearest-neighbor lookup) — documented; linear/cubic interpolation is immune (dependence
flows through the fractional weights).

Both modes ride the same `T <: Real` genericity as `Dual`; Dual-cleanliness in CI
effectively guarantees traceability.

---

## 11. Time and execution

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
is what restores consistency at each accepted boundary. The rule, binding for the periphery (§12):
**external readers (GUI, logging, network output) observe the signal table only at
step boundaries.** Mid-step contents carry no meaning.

### 11.4 Tier-2 localization mechanics

Trigger: a Tier-2 guard changed sign across an accepted step $[t_n, t_{n+1}]$.

- **Interpolant (lazy).** Build the cubic Hermite continuous extension $\hat{x}(\theta)$,
  $\theta = (t - t_n)/h \in [0, 1]$, from $(x_n, \dot{x}_n, x_{n+1}, \dot{x}_{n+1})$;
  $\dot{x}_n$ is the step's first
  stage, $\dot{x}_{n+1}$ costs one sweep, paid only on trigger. Uniform accuracy $O(h^4)$ — one
  order below the discrete solution, the standard pairing, and the event time can
  only ever be as accurate as the interpolant, so nothing more expensive is worth
  probing.
- **Probes run the sweep.** Guards read `y`, so evaluating a guard at an interpolated
  state means writing $\hat{x}(\theta)$ into the state buffer and running the sweep — the same
  rule as the RHS. One sweep per probe.
- **Root-finding: bracketed and derivative-free** (ITP or Brent; bisection is an
  acceptable fallback). The observed sign change *is* an unconditional convergence
  certificate. Newton/AD localization is rejected: guards are guaranteed C⁰, not C¹
  (clamps, table knots, saturated stretches where σ′ = 0), Newton discards the
  bracket for merely local guarantees, and its superlinear convergence saves a
  handful of microsecond probes per rare event. AD earns its keep in Jacobians, not
  in root-polishing a possibly-kinked bracketed scalar.
- **Post-event.** The boundary sequence runs at `t*` (below) → **interpolant
  invalidated** (the handlers made it a lie for `t > t*`) → resume integration from
  `t*` with the remainder step targeting `tₙ₊₁` → re-check guards on the remainder,
  under a bounded per-step event budget with a chattering diagnostic. Multiple
  Tier-2 guards localizing in one step fire at the *earliest* `t*` (ties fire
  together at that boundary, declaration order within the iteration, §11.6);
  later crossings re-localize on the remainder.
- **Shared blind spot, documented:** an even number of crossings within one step
  defeats sign-change detection in both tiers; the mitigation is step size, not
  machinery.

**Budget exhaustion degrades; it does not throw.** When a step spends its event
budget, localization stops for the remainder of that frame: the remainder step
completes, and any further crossings fire in the next boundary's ordinary
iteration — Tier-1 granularity for that frame — under a warning naming the
chattering event and the localization count. The degradation is a function of
the trajectory alone, never of wall clock, so row 80's pace-independence stands
untouched and the run replays identically. A `StepError` would misclassify an
expected modeling outcome as broken machinery (§16.8's doctrine). This also
makes the asymmetry with §11.6 principled: the once-per-event bound *there* is
structural — quiescence terminates because each event fires at most once per
boundary — whereas re-localization across a frame has no structural bound, so it
takes a budget, whose exhaustion degrades exactly as that section's would.

**Projection's reach is the boundary, not the probe.** Guard probes evaluate the
raw interpolated state — the same rule RK-stage RHS evaluations already live
under, being equally off-manifold, so sweeps must tolerate near-manifold states
and already do. Authority rests with the `t*` boundary: projection runs there,
and the §11.6 iteration's edge checks read the *projected* state. If projection
moves the state back across a guard, the event does not fire and the run has
published one extra boundary — harmless, and deterministic and pace-independent
like any other localization outcome (row 80). Per-probe projection is rejected:
it puts projection cost inside the root-finder's loop to correct an interpolant
drift far below the localization tolerance.

**`t*` is a boundary, not a frame.** A *frame* is a grid step
`[tₙ, tₙ₊₁]` — the unit of scheduling: input drain at frame top (§12.3), pacer
deadlines (§11.7), tick eligibility (§11.5). A *boundary* is a published
consistency point — where the §11.6 macro-sequence completes and a snapshot
goes out. Every grid point is a boundary; `t*` and boundary zero (§16.5) are
boundaries that are not frame tops. At `t*` the full §11.6 event phase runs —
[sweep → guards → handlers] iterated to quiescence, **once-per-event
accounting scoped to this boundary** (fresh again at `tₙ₊₁`, and at a second
`t*` on the remainder) — and the settled state is **published**: snapshot,
§12.8 boundary-counter increment, `stop_on` check (§15.5) — a crash localized
at `t*` ends the run from that snapshot. What does *not* happen at `t*`:
ticks are never due (`t*` is off the harmonic grid by construction; discrete
cells ZOH-hold through the sweep), and staged inputs are not drained — input
timing is a frame fact, and replay determinism must not depend on
localization arithmetic. The publication is not separately paced: the pacer
paces frame deadlines, and a `t*` snapshot publishes when computed, mid-frame
(wall-side placement below pacer resolution; the §11.7 invariant concerns
trajectories, which are identical). Replay pointers and error messages index
boundaries by a monotonic counter with recorded `t`; the trace stays
frame-indexed — `t*` boundaries consume no inputs.

**Endpoint policy and grid integrity.** The root-finder returns the
**holding endpoint of its final bracket** — the smallest probed point where
the condition holds. Consequences: **`t* = tₙ` is structurally impossible**,
not clamped away — localization only triggers when the condition was
not-holding at `tₙ`'s quiescence (§11.6 baselines), and the interpolant
reproduces the endpoint exactly (`x̂(0) = xₙ`, probe sweeps deterministic), so
the bracket's left end stays strictly not-holding and the returned point is
strictly later than the published, immutable `tₙ` (worst rounding:
`nextfloat(tₙ)`). And the guard observably *holds* at `t*`: handlers fire in
states where their own condition is true, and the post-fire baseline records
an actual observation rather than an assumption. **`t* = tₙ₊₁` exactly is
legitimate** (a crossing at the grid point: `σ(tₙ₊₁) = 0` both triggers
detection and is the root) and **degenerates to the grid boundary**: the
localization result is discarded and the event fires inside `tₙ₊₁`'s ordinary
iteration — bitwise the Tier-1 outcome, one boundary, one snapshot, no
zero-length remainder. A near-degenerate `t*` leaves a tiny remainder step,
numerically harmless (increments are `h′`-scaled); the real hazard is
bookkeeping, killed by rule: **grid times are indexed, never accumulated** —
`tₖ = t₀ + k·h` computed from the frame index (tick gating is already
counter-modulo, §11.5), and the remainder step *targets the grid point*, with
`h′` derived at use. `t*` is a float inside a frame, never an anchor anything
else is computed from.

### 11.5 Multi-rate tick scheduling

**Harmonic grid.** Every discrete component's period is an integer multiple of a base
tick period `Δt_base`, itself an integer multiple of the continuous step
($\Delta t_{\mathrm{base}} = n \cdot h$, $n \ge 1$). Ticks therefore land on step boundaries — the only place
anything discrete ever happens. Rejected: arbitrary periods via a time-ordered tick
queue, which forces variable-length steps and irregular real-time frames for a
generality nothing demonstrated wants.

**Discrete stages are gated to tick instants.** A discrete component's `h_z`/`h_zu`
run only at its own ticks; its cells hold in between (ZOH, stated in sweep terms). The
alternative — re-running its stages at every boundary — would let outputs drift
between ticks as fresh continuous inputs flow in, silently un-sampling a sampled-data
controller. Implementation consequence, recorded: the boundary sweep is not one fixed
list — discrete stage entries are gated by tick counters, so different boundaries run
different subsets of the schedule.

**Simultaneous ticks are already well-defined** by settled machinery: all due
components run their output stages in topological order within the sweep, then all due
`g` updates run after it, in any order — each `g` reads the table and writes only its
own `z` cell. The FCS cascade's intra-tick ordering is a sweep property, not an
update-order property.

**Assemblies: virtual for execution, rate scopes for declaration.** There are no
atomic assemblies, and no opt-in variant. Execution atomicity coarsens the schedulable
unit, which is exactly how artificial algebraic loops are manufactured — §5.3 at
assembly scale, a hazard Simulink documents for its Atomic Subsystems — while the
thing it protects, non-interleaved execution, protects nothing here: the signal table
makes interleaving semantically invisible (consumers read cells whose freshness is
guaranteed by topological order, not contiguity). Its other Simulink roles — code
generation units, enabled/triggered execution — have no counterpart in the consumers;
conditional behavior ("on ground, force `direct`") is mode logic. FlightCore's
whole-tree atomicity (children's `f_periodic!` called from the parent's, hence
`Subsampled`'s parent-relative multipliers) was an artifact of call-tree execution,
not a design requirement; the signal table dissolves it.

**Rate declaration is relative, at composition.** A discrete component or sub-assembly
is instantiated with an integer multiplier $K \ge 1$ (default 1) relative to its
enclosing scope; multipliers compose multiplicatively down the tree; the root fixes
`Δt_base` in seconds. Rationale: in a layered control architecture the *ratios* are
intrinsic to the design and travel with the assembly type (inner loop at `K = 1`,
outer loops at `K = 5`, whatever the deployment), while absolute rates are deployment
decisions made once at the root. Absolute-first declaration was rejected: it welds
deployment rates into reusable definitions, and replicating relative structure with a
shared base-period variable does not compose across independently authored assemblies
(parameter-threading ceremony, cf. §7). At build, all scoping compiles away to **one
absolute tick divisor per discrete component**; the executor gates entries by counter
modulo. Recorded limitations: a child cannot tick faster than its scope ($K \ge 1$) —
soft, since assembly multipliers are declaration sugar and factors can be refactored
onto siblings; and no phase offsets in the first cut (no demonstrated use).

**`Δt` has a single source of truth: the compiled schedule.** Each discrete
component's effective period arrives read-only as the `Δt` field of every
discrete-tier bundle (§5.2) — `h_z`, `h_zu` and `g` alike, and absent from
continuous bundles, so touching it on the wrong tier is a missing-field error
rather than a rule. It must be readable in the *stages*, not just `g`: per
§17.2, the discretized laws that actually consume `Δt` — a PID's
backward-difference coefficients, a LeadLag's Tustin transform — run in
`h_zu`; `g` is a copy. (A virtual property on the component handle,
`comp.Δt`, after FlightCore's `mdl.Δt`, is the natural-looking alternative,
and it is impossible here, not merely inconvenient: `comp`
is the author's own immutable struct, and two fields holding the same
immutable value are `===`-identical — the §13.6 argument — while carrying
different `rates` keys, so the period is a fact about the *schedule
position*, with nothing on the instance to hang a property on. The value must
arrive through the call; the bundle field is where.) Author rule: **never
store `Δt`, or any `Δt`-derived coefficient, as a component parameter** —
recomputing derived coefficients per tick is a few arithmetic ops, and a
cached copy is a second thing for gain-scheduling machinery to chase.
Relative declaration structurally enforces the rule for the period
itself: under scoped multipliers a component author *cannot* know their
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

**Within-round visibility: round-start `u`, live `y`, order inert.** §5.2's
cross-component claim — all guards and handlers read the same boundary
snapshot, plus their own component's refreshed ports — is delivered by
*binding time*, not by copying. After a round's guards are evaluated, the
framework materializes the `u` gathers of every fired component **before any
handler runs** (the executor builds the same args bundle it always builds,
merely earlier), and the per-event re-decode reuses that same round-start
`u`. Because cells hold references to immutable values, loading the
references early *is* the capture — the only mid-round table writers are the
re-decodes, confined to the firing components' own cells, and the
pre-materialized bundles are beyond their reach: no shadow table, no
allocation. `y` stays bound live: it gathers own cells only, disjoint from
`u`'s foreign set, so a later event of the same component sees its own
refreshed ports (handler returns are latched per-event, before `project` and
before the next fired event runs) while foreign transitions stay invisible
until the next round's re-sweep. `u` binds round-start uniformly — even for
a port contrived to be wired from the component's own stage-1 output; the
freshness rule is about `y` only. Consequence: cross-component handler order
within a round is **semantically unobservable** (bundles carry no foreign
state, and foreign reads are round-start), so the execution order — executor
component order, declaration order within a component — is fixed only to
keep the §15.4 cursor and the diagnostics stream deterministic, never
something a trajectory depends on. Rejected: live-table reads under that
canonical order (deterministic, but model semantics would then depend on the
executor's schedule order — a rewiring that permutes it silently changes
trajectories, the same structure-the-framework-inserts class this section
already rejects); a table copy per firing round (identical semantics, paying
an allocation the pre-materialization makes unnecessary). The trade,
recorded openly: a handler cannot opt into seeing a same-round foreign
transition — same-instant sequential coupling across components is a
cascade, one round per link, deterministic; coupling tighter than that
belongs inside one component, where declaration order gives exact sequencing
(the synchronous-languages position: a micro-step sees the pre-state,
effects appear at the next micro-step).

**Why once-per-event rather than a round cap.** Termination becomes structural —
rounds are bounded by the number of declared events, with no arbitrary K knob — and a
livelock (two FSMs toggling each other) resolves deterministically into "each fired
once, quiescence, re-fire next boundary," i.e. degrades to Tier-1 granularity instead
of burning a budget and erroring. The cost: an event legitimately re-enabled within
the same boundary waits one step — accepted at the same granularity Tier 1 accepts
physical re-crossings, and flagged when it occurs (a runtime warning, §15.2).

**Guard baselines and firing semantics.** "Newly-fired" means precisely this:
an event fires at a boundary iff its condition (§2.1: predicate true, or
`σ ≥ 0`) is observed holding in some iteration round, its **baseline** — the
previous boundary's quiescent sample — was not-holding, and it has not yet
fired this boundary. The baseline lives in **loop state, per event**, updated
at each boundary's quiescence from the final post-iteration samples. It is
detection bookkeeping, not model memory — correctly *not* in `z`: not
captured, not traced, reconstructed deterministically. A condition that holds
and *keeps* holding fires once, at the boundary where it first held — sticky
flags do not re-fire every boundary. **Boundary zero establishes the baseline
as nothing-holds**, so an authored condition already in holding territory
fires at `t₀` — §16.5's behavior, derived rather than asserted —
and a warm restart (`init!` re-runs boundary zero, §16.5) re-baselines from
scratch: conditions true in the newly applied state fire again at the new
`t₀`.

**Ticks stay outside the iteration, after quiescence.** The two possible couplings
resolve asymmetrically:

- *Events → ticks: handled by machinery already in place.* Due discrete components'
  output stages (`h_z`/`h_zu`) are gated *into* the sweep (§11.5), so every
  iteration round refreshes them for free — same `z` (their `g` has not run), post-transition inputs. At
  quiescence, their published outputs reflect the settled boundary instant, which is
  exactly what "sampling at t" should mean for a logically-instantaneous cascade.
  Earlier rounds' tentative values are internal scratch, like RK stage evaluations;
  §11.3 extends naturally: external readers observe the table only after the boundary
  sequence *completes*.
- *Ticks → events: structurally impossible.* A tick's output stages contribute nothing
  guards have not already seen (they run inside the sweep, from current `z`); its `g`
  update writes `z⁺` after the sweep, and `z⁺` is first decoded at the owner's *next*
  tick, so it is invisible to every reader within the boundary — the standard
  one-sample `z⁻¹` delay of sampled-data control, here enforced by
  construction. Nothing that happens after quiescence can flip a guard, so no combined event/tick fixed point exists to
  iterate.

The boundary macro-sequence, final form (boundary zero — initialization — is
the same sequence with an empty integrate, §16.5):

> integrate → project → **[sweep → guards → handlers]** iterated to quiescence
> (once per event) → all due `g` updates → logging / I/O staging.

The mixed case — a continuous component's handler and its discrete observers' ticks
landing on one boundary (engine `starting → running` under a 50 Hz FCS) — is decided
by the sequence: the transition fires in the iteration segment, the re-sweep re-runs
the FCS's stages against `running`-mode ports, and its `g` then updates from
post-transition values.

### 11.7 Real-time pacing

**The invariant: pacing is outside the semantics.** The pacer inserts waits between
completed frames and never reorders, skips or alters the boundary sequence. A
paced and an unpaced run with identical input traces produce bit-identical
trajectories — deterministic replay (§2.2) extends over pace. Interactive runs differ
only because their *inputs* differ. Detection policy is inside the semantics:
Tier-2 localization runs identically paced or unpaced, its sweep cost
absorbed as debt like any other expensive frame (§2.1).

**Wall-clock mapping: piecewise affine, re-anchored at every knee.** The map is
$\tau(t) = \tau_{\mathrm{anchor}} + (t - t_{\mathrm{anchor}})/p$, with the anchor pair as its reference point. A
live pace change re-establishes the anchor at the current `(t, τ)` so the new slope
applies only forward; keeping the old anchor would retroactively reinterpret the
entire elapsed run at the new pace (deadlines minutes in the past or future after a
long session). Un-pause re-anchors for the same reason. Debt is cleared at re-anchor:
a deliberate user action is a natural sync point, and the counters record what was
forgiven.

**Deadline law: absolute schedule with bounded debt.** Frame deadlines come from
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

- **`margin = 0` — pure sleep:** cheapest CPU; bursty frame spacing, but the
  absolute schedule still delivers the exact *average* rate through debt repayment.
  The spin phase buys regularity, never rate correctness.
- **`margin` = a few ms — hybrid (default):** sleeps ~90% of a 20 ms budget, lands
  within µs of the deadline at a few percent of one core.
- **`margin = ∞` — pure busy-wait:** FlightCore's behavior, maximum frame
  regularity at one pinned core; the "best attempt at real time" mode is the knob's
  endpoint, not a separate mechanism.

When the frame budget is at or below the margin (e.g. `h = 0.01` at `p = 5` → 2 ms),
the hybrid degenerates to pure spin per frame by construction. Rare wake-ups past the
deadline are overruns, absorbed as debt. Which primitive the coarse phase uses —
task-yielding `sleep` vs. thread-blocking `Libc.systemsleep` — is settled in §12.7:
the coarse phase uses task-yielding `sleep`, with `margin` absorbing its overshoot.

**Diagnostics.** Overrun count, current and peak debt, forgiven-debt events, wait
statistics — published as framework status for GUI and logs (today's `SimControl`
fields are the precedent).

**Forward pointers.** The wait interval is the natural staging slot for externally
injected inputs, applied at the next boundary; the staging rules — and the
concurrency model generally, which §11.3 constrains but does not decide — belong to
§12.

---

## 12. Runtime periphery

GUI, input devices, network I/O and logging, and the concurrency model binding them
to the §11 loop: architecture and data exchange (§12.1–§12.6), loop scheduling and
thread policy (§12.7), the next-snapshot wait (§12.8), shutdown (§12.9), the
mid-run mutation doctrine (§12.10), and the run lifecycle with its partial
advance (§12.11).

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
2. **The drain:** exactly one point, at the top of each frame — never at a `t*`
   boundary (§11.4) — where the
   loop takes the staged batches and applies them to the root input slots. Between
   drains the loop owns its data exclusively — no lock is held during stepping, ever.
3. **Publication (outbound):** at the end of each boundary sequence the loop publishes
   an immutable snapshot; readers observe it without coordinating with the loop
   (§12.2).
4. **Control:** pause/pace/stop on a separate few-word atomic surface (§12.6).
5. **Task topology:** one loop, one task per rostered device except the
   GUI, all run-scoped: `run!` spawns one task per non-GUI roster entry
   after device `init!`, a mid-run `attach!` spawns immediately, and §12.9
   joins them all at every stop (§12.11). **The GUI is pinned; the loop is
   the movable piece.** CImGui ties rendering to the calling (main) task,
   so the GUI is the one rostered device *without* a spawned task: with
   `gui = true` the loop moves to a spawned task for the duration of the
   run and the calling task renders; otherwise the loop runs on the
   calling task — the batch register, what §15.4's synchronous rethrow
   presupposes, and what lets parallel batch sweeps thread `run!` inline
   with no nested task fan-out. Either way `run!` blocks its caller until
   the run ends; what varies is what the calling task spends the run
   doing. Spawn-inside-`run!` *is* the start gate — a task exists only
   once the run it serves exists — and any first-boundary synchronization a
   device needs is §12.8's counter-plus-condition predicate wait, never an
   `Event` latch: FlightCore's `io_start` gate is the once-per-run version of
   exactly the race §12.8 rejects, and inheriting it would re-import that
   race for §12.11's re-run cycle.

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
with an immutable, coherent world for as long as they like. `latest(sim)`
hands the same value to the calling task — §12.11's inspection register.
Wait-free in both
directions: a wedged reader cannot delay publication by a nanosecond; the loop cannot
tear a reader's view. Publication happens only after the boundary sequence completes
(§11.3 as extended by §11.6).

**Binding rule: nothing reachable from a published snapshot is ever written again.**
The table's immutable values (§4.1, §9) make the compiler enforce most of it; the
rule is what the soundness of lock-free reading rests on.

The captured table includes the component-local cells (`local_types`, §13.3) — the copy is
mechanical, and they serve the author's own debug panels; presentation layers (log
export, GUI listings) filter to the public contract by default. **It also includes
the root slots** (§17.4): slots are source cells of the
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

**Retention: the trace's kill switch, plus decimation.** The log takes the same
plain on/off switch the trace has (§12.3), and additionally a keep-every-kth
policy (`log_every`, Appendix B). What makes decimation admissible here and not
there is the derived/primary split (row 38): the log is recomputable from the
trace by replay, so a thinned log costs resolution in a *view*, never a record —
whereas thinning the trace would destroy the only primary account of a session
(row 29). Decimation is a retention policy only: every boundary still runs, is
still published to live readers, and still enters the trace.

**Output-device bindings are snapshot-path bindings.** An output device (telemetry, the XPlane visualizer, disk
streaming) consumes snapshots via §12.8 and addresses what it reads with
§16.4's selectors — any cell, the diagnostic register admitting `get_local`
entries and deep paths — resolved at attach against the `Build` with
did-you-mean and compiled to one gather, so `map_output` receives a labeled
NamedTuple rather than performing its own path lookups (§17.4's obligation: a
substitution that breaks a binding fails at attach, not with silent garbage
UDP). This is **diagnostic observation** (§15.5):
human-facing, no effect on run semantics — the same register as the log
retaining the full table and the GUI's deep-reading panels; local cells are
accessible (the §13.3 presentation filter is a default, not a wall — explicit
naming is intent). **A binding chooses its register**: a deep path is the
*inspection* register — zero promises, free access, right for looking at
*this* build; an exported output face (spelled `get_face(name)`, §16.4) is
the *integration* register — named,
curated, meaning-stable under substitution (§17.4's writer-independent
semantics), right for consumers that outlive the build they were configured
against. Attach validation converts *structural* drift to loud errors in both
registers; only faces protect against *semantic* drift — a substituted
aircraft publishing the same path at the same type with a different meaning, a
CG velocity under a name read as body-origin velocity — and nothing else can:
meaning is not in the schema. Semantically generic consumers (a visualizer
needs pose; every aircraft has one) should therefore bind faces, and aircraft
families should export the conventional surface such consumers need — a
library/migration deliverable (§19) — with wrapper types making face semantics
structurally checkable (`VelocityData`, its `v_eb_b` defined *at the type* as
body-origin velocity: a bare vector doesn't wire, and wrapping the wrong
quantity is a deliberate lie, not a drift).

### 12.3 Inbound: root input slots and per-device staging

**The write surface is root input slots** — and a root slot *is* the root
assembly's exported input face (§13.6): routed inward to consumers, produced by no
component, fed by the parent's wire at every non-root level — and at the root
there is no parent. (A root slot is usefully read as the output face of the
one producer the build never sees — the periphery and the services: slot
exclusivity below is a producer's one-writer right, and §16.6's totality is
its completeness obligation.) No dedicated vocabulary survives (`add_input!` in the early
sketches is dead). Slots are sources to the build-time scheduler, constants within
a frame, and the *only* thing the periphery may write (the GUI reaches them
through §12.5's resolution; control commands are not writes, §12.6); devices,
mappings, the trace and the GUI write path address them by **face name** (§13.6):
structural slash paths never cross the periphery's *write* boundary — the write
side speaks the root contract's names only. (The read side is diagnostic
observation: snapshot consumers, output devices included, address table cells
by path — §12.2/§15.5.)

**Slot exclusivity: one writer per slot at any time** (§17.4). A
device claims its slots at attach; claiming an already-claimed slot is an
attach-time error, and detaching releases the claims (a released slot's GUI
widgets re-enable, §12.5). Exclusivity replaces any cross-device conflict *policy* —
attachment-order precedence at drain, say — because such a policy resolves races the
case study shows nobody wants: every dual-writer field in the C172X demo is a joystick stream
shadowed by a GUI mirror, where simultaneous live writing is a bug. Per-device
cells, the CAS merge and the atomicswap drain all stay — they serve atomicity and
coalescing, not arbitration.

**A claim is what a device *may* write, not what it will.** Data-dependent
write-sets are ordinary: a UDP/JSON peer writes whichever subset of faces the
incoming message names, and `map_input` is arbitrary user code the framework
never inspects. Such a device therefore claims the **binding's enumerated
allowed set** — the faces the binding table lists, whether or not any given
batch touches them — and the claim is registered at attach exactly as a
joystick's is. A broad claim costs liveness: every enumerated face is claimed
for the device's whole attachment, so §12.5's derived-liveness rule renders its
GUI widget read-only even on faces the peer never writes. Narrow the binding to
narrow the claim — the enumeration *is* the interface.

**Every writer has a write surface, and the drain enforces it.** The drain
applies a batch entry **iff the named face is inside the writer's surface at
drain time**; anything else is discarded with a runtime warning (§15.2). A
surface arises in one of two ways:

- **Enumerated** — an ordinary device's surface is its claim set: static for
  the attachment, exclusively its own (claims are disjoint by construction),
  and binding-bounded even where no one else is involved — a mapping that has
  drifted onto an unenumerated face is a diagnosable anomaly
  (`OutOfClaimEntry`), never a silent write, claimed or not.
- **Derived** — the **interactive register**'s surface (the GUI, §12.5, and
  `stage!` below) is the currently-unclaimed face set: the complement of the
  union of all claims, computed rather than staked. It is shared among the
  interactive writers and guaranteed to nobody — a mid-session `attach!`
  shrinks it, and a batch staged before the claim and drained after it is
  `StaleInteractiveEntry`: no one's bug, the surface moved between staging
  and drain. §12.5's widget liveness is this surface made visible.

One rule, two surface kinds, two warning kinds. The GUI is not an exception
but the resident of the second kind — one device kind (§12.4), two *binding*
kinds (enumerated vs. derived), one drain rule — and opportunistic writing by
autonomous devices does not exist: a device that wants a face enumerates it.
Cross-writer races on one slot therefore cannot arise structurally (claimed
faces have exclusivity; unclaimed faces admit only the interactive register),
which is what keeps drain order a diagnostic fact (below) and lets a drained
GUI value simply stay (§12.5). The interactive register also has a
**task-free entry point**: `stage!(sim, "face" => value, ...)` stages a batch
from the calling task itself — the harness/REPL write path (§12.11) —
drained, traced and surface-checked exactly as the GUI's. Among interactive
writers the harness cell drains **last** — the explicit hand of code beats a
widget interaction — sequencing within one register, not the cross-device
conflict policy row 44 rejects; in practice the two rarely coexist, since
during a GUI run the calling task renders (§12.1).

**Slot initial values are owned by the init/trim services** (§17.4). Input declarations are bare types (§13.2) and carry no defaults, but a
slot unfed by any device must hold a defined value from the first frame (today's
`U()` constructors provide these: `mixture = 0.5`). Export-entry defaults were
rejected: the trim service writes slot values it *solved for* (throttle,
elevator) — not declaration constants. `init!` establishes every slot and the
trace header captures the result; totality is enforced pre-write at
`init!`/commit (§16.6).

**The roster is itself a handoff.** Devices attach and detach at any time,
running or stopped. The attached-device roster — cells, claims, attachment
order — is an immutable array republished by `attach!`/`detach!` through one
atomic reference operation (§12.1's rule applied to the periphery's own
bookkeeping); the loop acquire-loads it once at frame top and drains that
snapshot, so a device attached mid-frame is drained from the next frame — the
same granularity its staged writes already have. Attachment order is the
order within the currently published roster, carried by per-device sequence
numbers (stable across detaches). The trace tags entries with a stable device
id, never a roster index — replay provenance survives roster churn. Detaching
releases the device's claims (§12.5's widgets re-enable), and a device task
exiting through `should_close` runs the same detach path: an unplugged
joystick's claims release without any explicit call (§17.4). The §12.7
thread-budget warning re-evaluates as the roster grows.

**Staging: one atomic cell per attached device.** Each cell has a single writer — its own device task — and holds that device's latest
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

**Mappings are binding data, not shaping code** (§17.4). A mapping is
a declarative table — axis/button → slot path, plus per-axis conditioning
parameters (deadzone, expo strength) applied by a shared pure helper on the
device task. The boundary is set by the face contract: **a face's meaning is
writer-independent**, so faces carry *post-conditioning* semantics — a GUI slider
or a script writes the same command a curved stick delivers (running a mouse drag
through a deadzone would be absurd); this GUI-parity test is what places
conditioning upstream. Aircraft-semantic derivation (the C172X `q_ref = q_sf ·
axis` fan-out) must *not* ride along: it is FCS design and lives in-model — in
the avionics, or accepted as a small per-aircraft×device mapping entry (an
aircraft-design fork, §17.4). The trace records post-conditioning levels —
exactly what the model consumed, so replay is exact; raw-stick provenance (re-run
a session through *different* curves) is the known, accepted loss. Edge logic
follows the levels doctrine: devices stage monotonic press counters; accumulators
(trim offsets, flap detents) are model state, not mapping state (§17.4).

**The input trace** is the sequence of drained, device-tagged batches per frame. It
extends §11.7's determinism end-to-end: replaying a recorded interactive session —
staging fed from the recording, no devices or mappings present — reproduces the
trajectory bit-identically.

**The trace header captures the full initial state** `(x, m, z)` **plus the
initial root-slot values** at `init!` (an unfed `mixture = 0.5` never
appears in any batch, so replay is broken without them; the init/trim services own
slot initialization (§16.6), and the header capture extends naturally) — the one
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

Rejected shapes (both torture-tested in §17.3): **per-slot atomic cells** — the
simplest (no merge machinery, and a per-slot layout cannot lose independent writes)
but same-slot conflicts resolve by hardware store order, i.e. sub-frame wall-clock
phase (run-to-run behavioral variance, §17.3), peeks are cross-device, the trace
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

### 12.5 The GUI write path: port resolution, peek, staging contract

Panels remain per-component extensions in FlightCore's style — `GUI.draw!(ctx,
::LowPassFilter)`, discovered by walking the assembly — but widgets name **the
component's own ports**, never root slots. The build-time wiring answers, statically
and exactly: *is this port transitively driven by a root input slot, and which one?*
Every input port has exactly one source (§8), so the resolution is total:

- **root-driven → live widget**: peeks and stages the resolved slot through the
  GUI's own cell;
- **component-driven → read-only rendering**: displays the driven value from the
  snapshot, visually distinct, with the source as provenance ("driven by
  `avionics/throttle_cmd`" — the canonical slash form of §13.6).

This retires FlightCore's dead-slider convention (the `Cessna172Xv1` throttle: the
engine panel's slider visually live, silently overwritten by the avionics every
cycle — who commands what living in the user's head) and replaces it with checked
structure: **a widget is live exactly when the underlying input is yours to command
in this configuration.** User-commandability is a wiring decision made where
configurations are made; command-plus-manual-override is a mux component with a
root-wired select — explicit structure, not two writers racing (the same race as
§17.3's drag phase, ruled out the same way). The obligation this places on the
GUI: read-only rendering is first-class, not an error state — the author of
`input_slider!` cannot know at authoring time whether it will be live.

**Liveness is a derived property, and resolution is transitive.** A widget
is live iff its port's feed chain — walked through wires and exports across *all*
levels, not just the local assembly — terminates in a root slot, *and* that slot
is currently unclaimed by a device (§12.3 exclusivity). There is no per-port
"GUI-controlled" marking anywhere: the export chain is the marking, written by
the one author entitled to write it (a component's ports become GUI-commandable
exactly when the assemblies above surface them). The switch between "driven by
its own panel" and "driven by an external provider" is therefore automatic — at
build time by wiring archetype (a scripted `World` wires a scenario component
into the same faces the interactive `World` exports to root), at run time by
device claim state. Rejected: nominally-connected ports with a GUI
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

**Staging contract: widgets stage on interaction events only.** Value widgets (sliders, drags) stage
the new absolute level on edit; edge widgets (buttons) stage on activation, as a
level computed from the peek — a flaps button peeks the current counter `k` and
stages `k+1`. The levels doctrine makes this safe by construction: repeated staging
of the same level within a drain window is idempotent (no repeat-increment hazard),
and multi-click within one window counts correctly through the own-pending-first
peek (`k` → stage `k+1`; second click peeks pending `k+1` → stages `k+2`). Held
buttons do not re-stage — after the drain applies and the snapshot catches up,
re-staging from the peek would auto-repeat at frame rate; the activation edge is
the intent.

The alternative — active widgets staging on *every* render pass — lost its
motivation to slot exclusivity (§12.3): stage-every-pass would exist to win every
drain against a streaming device sharing the slot for the grab's duration (§17.3's
drag phase), but a slot the GUI can write is by definition unclaimed, so once staged
and drained a value simply *stays*; there is nothing to reassert against. Nor is it
worth keeping as insurance: if an
anomalous writer ever touched an unclaimed slot through a framework defect, the
slot visibly fighting is a diagnosable anomaly — continuous re-staging would mask
the invariant violation at render rate. Side benefit: staging traffic (and trace
noise) drops from render-rate-while-grabbed to actual edits. **Claim-transition
policy:** if a device claims a slot mid-interaction (attach during a drag — a
deliberate act concurrent with a held grab, vanishingly rare), the widget flips
read-only at the next render and the drain discards the stale entry —
`StaleInteractiveEntry`, §12.3's derived surface moving between staging and
drain — with a warning.

The panel-authoring calling convention — what the drawing context carries,
how widgets name their component's ports, how an assembly's panel composes
its children's — is deferred to migration (§19), where it is co-designed
against the GUI library. Its constraints are fixed here: panels name their
own ports by face-name string; resolution to root slots is baked at attach,
never performed at render; liveness and peek arrive through the
framework-supplied context, never by reaching into the loop; and assembly
panels compose children by path.

### 12.6 Control plane

Pause, un-pause, pace changes, `margin` changes, stop: a few scalar fields on a
separate atomic
surface, consulted by the loop at frame top and inside its wait and pause states.
`margin` (§11.7) rides here for the same reason `pace` does — it tunes the
wait, never the arithmetic — so retuning the coarse/spin split mid-run is safe
by construction.
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
tight for the attached population — re-checked as devices attach mid-run
(§12.3) — naming the `julia -t` remedy; sizing guidance:
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
without polling: a monotonic **boundary counter** published with the
snapshot — counting *published boundaries* (grid, `t*`, boundary zero; §11.4),
not frames, so consecutive wakes are not necessarily `h` apart — plus one
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

1. **Initiation:** `t_end` reached, a control-plane stop (GUI, device handle,
   code), or a `stop_on` face reading `true` in the just-published snapshot
   (model-detected termination, §15.5). The loop always completes the current
   boundary sequence — never stops mid-frame — publishes the final snapshot,
   then sets the sticky stopped status.
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
   (§12.7 heartbeat) and abandoned with a warning (`DeviceJoinTimeout`,
   Appendix C) rather than hanging `run!`.
   The GUI, having no spawned task (§12.1), is outside the join: its render
   loop is the calling task's own occupation of `run!`, exits by the same
   `running(handle)` predicate as any device loop, and `run!` returns after
   the joins.
6. **Device-initiated paths:** `should_close` (window ✕, peer EOT) exits the
   device's own loop; with `should_abort` set it also requests a sim stop,
   otherwise the sim continues with the device absent (its cell stops filling —
   the loop is structurally indifferent). A crashing device task is caught by the
   framework wrapper and follows the same path, logged with the device's name
   (`DeviceCrash`, Appendix C).
7. **Loop-side failure** runs (1)–(5) from the catch path — specified in
   §15.6: the failed boundary is discarded and the previous snapshot promoted
   to final (FlightCore's `SimulationTermination` catch path was the precedent;
   the exception-based termination idiom itself has no place here, §15.5) — so devices
   unwind cleanly regardless of who died.

After (5) the task set is empty — device tasks are per-run artifacts
(§12.1) — and `shutdown!` has released each device's OS resources. What
survives a stop is the roster entry: binding, claims, stable device id
(§12.3); never a task, never a live resource. The next `run!` re-runs device
`init!` — resource acquisition is per-run; FlightCore's
create-a-new-socket-each-`init!` in network.jl is the precedent — and spawns
fresh tasks against the re-armed §12.8 counter. While stopped there are no
device tasks, so `should_close` and the §12.7 liveness heartbeat are
run-scoped observables; a device unplugged while stopped surfaces as the next
run's `init!` failure, reported by name through (6)'s crash path.

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
fresh by topological order (the callback ran post-step, one boundary staler);
the pure timetable script is a one-liner reading its bundle's clock
(`h_zu(s, (; t)) = (; offset = profile(t))` — exact at its own ticks, no
latching); and
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
Initialization and trim are stopped-sim workflows (§16's first-class services),
where no concurrency perimeter exists — no loop, no devices, plain single-task
code. The guarded-addition shape is on record should demand appear: a traced,
boundary-executed intervention command applied through project → sweep → publish,
so no consumer ever observes un-decoded state.

**The doctrine, final form:** while a simulation runs, the periphery stages
root-input writes and issues control commands — nothing else, structurally.
Anything that wants to poke the model mid-run is an *input* in disguise (wire a
slot and guard), *model behavior* in disguise (add a scenario component), or a
*wall-clock interaction* (attach a device). Graceful termination follows the
same shape (§15.5): a declared condition in the model plus `stop_on` policy at
deployment — never a callback, never a thrown exception.

### 12.11 Run lifecycle and partial advance

A `Simulation` moves through four states: **built** (stores allocated,
nothing authored), **initialized** (`init!` completed boundary zero, §16.5),
**running**, and terminally **stopped** or **errored** (§15.4). `init!` is
mandatory: `run!` or `step!` on a simulation whose boundary zero has not
completed is an error in §15.2's kind set naming `init!` — distinct from
`UninitializedSlots`, which fires *inside* `init!` (§16.6). The loop runs on
the calling task unless a GUI claims it (`gui = true` moves the loop to a
spawned task, §12.1); deviceless, `run!` is fully synchronous — the batch
register (§12.1: a batch run is the same loop with empty staging), and what
§15.4's synchronous rethrow presupposes.

**Partial advance.** `step!(sim; frames = 1)` advances whole frames
synchronously through the ordinary frame sequence — drain, integrate,
boundaries, publication — and returns; a stepped simulation is bit-identical
to the same frames under `run!`. `step!(sim; t_plus = 10.0)` is the duration
spelling, mutually exclusive with `frames`: whole frames until the boundary
time first covers the duration — the migration suite's advance-by-duration
idiom. This is the test-harness register (advance,
assert, advance) and the REPL register (fly a while, inspect, continue);
neither is a script, so §12.10's scenario-component doctrine does not absorb
them.

A stepping session is **deviceless by construction**: device tasks are
per-`run!` artifacts (§12.1) and a device loop's `while running(handle)` is
false outside a run, so `attach!` while stepping only registers, exactly as
while stopped. The frame-top drain still runs — `step!` frames stay
bit-identical to `run!` frames — and what it drains is the **harness
register**: `stage!(sim, "face" => value, ...)`, §12.3's interactive
register with the calling task as writer. Staged batches are ordinary
batches — traced, so replay and bit-identity hold; applied at the next frame
top; surface-checked like any writer's (§12.3). The read half is
`latest(sim)`: the same immutable snapshot value a device handle acquires
(§12.2), navigated directly for assertions. Advance-assert-advance is
`stage!` → `step!` → `latest`. Both entry points work under `run!` too — the
harness register is not step-scoped — and an inspection accessor leaves
§15.5's rejection of closure-based termination untouched.

**Status, termination and the `run!` seam.** Between `step!` calls a simulation
reports **initialized**: no loop task exists, so `running` would lie, and
nothing is terminal, so `stopped` would too — the state reads "boundary-
consistent and ready to advance", not "sitting at boundary zero". `run!` may
therefore follow `step!`, continuing from the current boundary; so may another
`step!`. Termination policy is honored throughout, as bit-identity requires:
`t_end` reached, or a `stop_on` face holding at frame 3 of `step!(sim; frames =
10)`, ends the run there through the ordinary §12.9 tail and leaves the
simulation `stopped`. `step!` therefore returns the number of frames
**actually advanced** — the requested count in the ordinary case, fewer when
the run terminated inside the call, which is how a harness detects the
truncation without inspecting the clock.

**Re-running.** `stopped → init! → run!` is the supported cycle: `init!`
re-runs boundary zero from its condition (warm restart = `capture` → tweak →
`init!`, §16.1) and clears the trace, the log, *and* any batches still in
staging cells — the recording registers restart with the run they record, and
no stale batch survives to clobber the boundary zero it predates. Device attachments persist across
re-initialization — attachment is orthogonal to the run lifecycle (§12.3) —
and persistence means *roster* persistence: binding, claims and device id
survive; tasks and OS resources do not (§12.1's per-run topology, §12.9's
teardown). Each `run!` re-initializes every rostered device and spawns its
task; `attach!` while stopped only registers — the task appears at the next
`run!`. Run policy is re-bindable per cycle: `t_end` and `stop_on` are `Simulation`
defaults that `run!` may override for the run it starts (§15.5), so a second
run — or a `step!` register between two runs — can stop on a different clock
or a different face set without a rebuild. `errored` is terminal (row 59):
reproduction is trace replay, not resurrection.

---

## 13. The declaration layer: components and assemblies

How an author spells a component: where the structural facts live, what the build
takes as authoritative, and what is checked against what. §13.1–§13.4 cover the
component side; §13.5–§13.8 the assembly side; the build pipeline is §14 and the
stopped-sim service spellings are §16. Concrete syntax below is near-final in shape
but still illustrative in spelling. The sketches (`sketch_decoder.jl`,
`sketch_io.jl`) are written against this layer and the services spellings.

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

**Contracts are functions of the type, not of the instance.** A leaf's contract
declarations — `input_types`, `output_types`, `local_types`, `events`, and the
shapes of `init_x`/`init_m`/`init_z` — must be determined by the component's
**type**, its type parameters included, and never by its field *values*. The
value-discarding signature (`input_types(::Engine)`) is the visible form of the
rule, and the idiom for a contract that genuinely varies is the type parameter,
not the field: `SumJunction{Wrench, 3}` (§6), `Or{N}` (§15.7) — arity is spelled
in the type, at the price §6 states openly. The reason is §14.7's entry typing: a
component's bundle is a `NamedTuple` whose key set *is* its contract's, and an
executor entry carries what selects code in type parameters and what is plain
data in fields. A key set derivable only from field values would have to either
climb into the type parameters anyway — multiplying specialization and changing
§14.7's chunking cost model — or sit in fields, dissolving the static typing that
§5.1's zero runtime graph logic, §9.5's allocation invariant and §14.5's
fold-away conformance test all rest on. The build reads each declaration once,
against the concrete instance, so a value-dependent contract does not announce
itself: this is a rule authors keep, not a check the build can run. **`workspace`
is the one exception**, and explicitly so: it is the by-allocation register (row
77) — an allocator the framework *calls*, not a schema it *walks* — and it
legitimately takes sizes from the instance (`workspace(c::KF, ::Type{T})` reading
`c.n`, §9.3), because no entry type is derived from it.

### 13.2 The declaration inventory

```julia
struct Engine <: AbstractComponent
    ω_idle::Float64; ω_min::Float64; J::Float64      #parameters: plain struct fields
    ω_rated::Float64                                 #unread here; §16.2's shipped condition uses it
end

#state stores: declared by initial value — types derived, nothing to drift
init_x(::Engine) = (ω = 0.0,)
init_m(::Engine) = (phase = off,)                    # off | starting | running

#input contract: bare NamedTuple of types at their Float64 faces
input_types(::Engine) = (throttle = Float64, starter = Bool,
                         fuel_available = Bool, M_load = Float64)

#output contract = the public interface, at concrete nominal types (§13.3)
output_types(::Engine) = (M_shaft = Float64, P = Float64, ω = Float64)
#ω names a state field no stage produces → auto-published at stage 1 (§5.2)

#stage and update functions destructure their bundle by name (§5.2)
function h_xu(eng::Engine, (; x, m, u))
    M_shaft = m.phase === running ? torque_law(eng, u.throttle, x.ω) : zero(x.ω)
    return (; M_shaft, P = M_shaft * x.ω)
end

f(eng::Engine, (; x, y, u)) = (ω = (y.M_shaft - u.M_load) / eng.J,)

#events: ordered and named — order is load-bearing (§5.2, §11.6); tier by the `localize` flag
events(::Engine) = (
    start    = Event(start_guard, start_handler),                               # Tier 1 (default)
    ignition = Event(ignition_guard, ignition_handler),                         # Tier 1
    flameout = Event(flameout_guard, flameout_handler; localize = true))        # Tier 2
start_guard(::Engine, (; m, u)) =                        #manual trigger: an input (§12.10)
    m.phase === off && u.starter
start_handler(::Engine, _) = (; m = (; phase = starting))          #no `x` key: no reset
ignition_guard(eng::Engine, (; x, m, u)) =                #predicate form
    m.phase === starting && x.ω > eng.ω_idle && u.fuel_available
ignition_handler(::Engine, _) = (; m = (; phase = running))
flameout_guard(eng::Engine, (; x)) = eng.ω_min - x.ω      #continuous form: localizable
flameout_handler(::Engine, _) = (; m = (; phase = off))
```

The inventory, and where each schema fact gets its authority:

- **State, modes, discrete state** (`init_x`, `init_m`, `init_z`): declaration
  *by initial value* — the type is derived from the
  value, so there is no second artifact to drift and no separate type declaration to
  check. (The workspace is the exception to the register: it is declared *by
  allocation* — `workspace(::C, ::Type{T})` continuous, `workspace(::C)`
  discrete, the method being the allocator — because a workspace is not
  memory and none of the by-value arguments below cover it; §9.3.) This is the boundary of legitimate derivation: deriving from another
  declaration is sound; deriving from evaluated user code is not. Rejected:
  declaring types here too, `input_types`-style, with §14.3's
  `probe_value` synthesizing the initial values. The declared values are the
  condition substrate's base layer (§16.1's overlays fall back to them leaf
  by leaf, and the compiled `m`/`z` writers bake `merge(defaults, overlay)`),
  so there must be an authored value under every leaf; synthesized initial
  state would cross the probe-value barrier §16.6 makes structural (a
  fabricated zero is a fine probe input and a terrible flight condition —
  states no less than slots); and every field where synthesis picks wrong
  (modes, `Ranged` values excluding zero, trim-sensitive states) would need
  an authored default *beside* its type — the per-field two-register
  protocol §16.2 kills for `initialize` specs, aggravated by types being
  first-class values in Julia (the two registers distinguishable only by
  `isa Type`). The asymmetry against `input_types`/`output_types` is one of kind, not
  style: contracts describe table cells, recomputed from scratch every
  sweep, needing only types; `init_*` describe stores — the model's memory,
  which must have contents before the first sweep can run.
- **`input_types(::C)`**: a bare `NamedTuple` of types — zero framework vocabulary, no
  wrapper types (the last candidate, `Reduce`, died with reduce-ports, §6). Entries
  are **face constraints, not cell types**: each states a bound on the
  producer's *nominal* face, checked at wiring time as `producer_face <: entry` —
  one uniform rule, which for a concrete entry degenerates to exact equality,
  concrete types being final. Abstract entries state **structural
  substitutability** — several concrete producer types admissible behind one
  stable face, §7's field handles being the demonstrated client
  (`terrain = AbstractTerrainField`). They are never needed for eltype
  genericity: an eltype-generic producer's nominal face is concrete by
  construction, so `SVector{3, Float64}` wires to it as-is. (Names-only
  contracts were rejected — they lose the wiring-time type error and standalone
  checkability.) Inputs are the component's *requirements*: §8's
  unconnected-input error, over-wiring detection and did-you-mean typo messages
  are only definable against them. Because entries are bounds, nothing is ever
  "overwritten": cell types are single-sourced from the producer side per
  activation (§14.1), and a `Dual`-carrying cell behind a `Float64` entry is
  the design working, not a promise broken. The code-level complement is the
  **genericity obligation** — whatever scalars the wiring delivers, the
  consumer's math promotes — checked by the `Dual` probe, never declared:
  **declarations record choices; obligations are checked.** A per-leaf
  genericity marking would be a constant function across components — zero
  information (killing the envelope reading of symmetric `T`; the *predictive*
  reading is impossible outright: an input's activation type depends on
  who feeds it — a continuous producer delivers `Dual` under a `Dual`
  activation, a gated-off discrete producer `Float64` — so a symmetric
  declaration would force the consumer to state its producer's tier and break
  substitution behind the same face). **Root slots are the one place an entry
  types a cell**: produced by no component, a slot has only the consumer
  declaration to take a type from, and only a *tight* bound determines one — a
  face surfacing as a root slot must resolve to a concrete declaration (staging
  cells, the trace header and `probe_value` all need it; abstract-at-root is a
  build error). Under fan-out the slot type is the unique concrete declaration
  among its consumers — two different concrete declarations remain an error —
  and abstract co-consumers are checked against it.
- **`output_types(::C)`**: a concrete `NamedTuple` of *nominal* types, both
  tiers — same species as `input_types`. (An earlier design gave this method a
  `::Type{T}` signature, making the declaration a function of the sweep scalar; it
  was abandoned — rows 33 and 79.) Cell types per activation are
  derived by the framework's **leaf walk** over the declared type (§14.1): on
  continuous producers, `Float64` leaves — including type parameters,
  `RQuat{Float64}` → `RQuat{Dual}` — follow the activation scalar;
  `Int`/`Bool`/enum leaves pin; reference-typed fields pin (a §7 bulk-data
  handle's grid is frozen build-time data, never activation-dependent). Only
  `Real` *type* parameters follow the scalar: value parameters pin, being
  structure rather than number (`Ranged{Float64, -1, 1}` walks to
  `Ranged{Dual, -1, 1}`; the bounds are not scalars to re-type). The companion
  obligation is **constructibility at the walked type** — a Tier-1 type must be
  buildable at the activation scalar — and the `Dual` probe enforces it by
  construction: it builds real values, so a type whose constructor cannot
  accept them detonates at the probe with its own name in the message. On
  discrete producers the walk pins wholesale — the Tier-3 exemption (§9.2),
  now enforced by the typing rule rather than spelled in a signature (tier is
  known from declaration shape, §13.5). During a generic sweep, gated-off
  discrete producers hold their `Float64` values, consumers gather mixed
  tuples, and promotion does the rest — semantically exact, since a frozen
  discrete output is a constant with zero partials, which is precisely what
  "linearize the continuous dynamics with `z` held" means. What makes
  under-the-hood substitution safe — the standing objection being that it cannot
  distinguish honest `Float64`s from eltype `Float64`s (row 33) — is the **embedding
  guarantee** (§14.5): a `Float64` observed at a parametrized leaf under a
  non-nominal activation implies no `Dual` entered its computation (promotion
  is airtight; there is no lossy cast), so its true derivative along every
  seeded direction is zero and embedding it as a zero-partial constant is
  exact. The honest/eltype distinction no longer needs declaring because both
  cases get correct treatment automatically. Author-side payoffs: no
  `where {T}` ceremony, no per-leaf placement decisions, the forgotten-`T` bug
  class gone, and piecewise branches returning literal constants
  (`flow > 0 ? f(x) : 0.0`) legal as written — zero partials is the derivative
  of a locally-constant branch. Differentiation participation is selected per
  invocation by seeding (§16.10), never by typing.
  The stores walk by the same rule: the type derived from `init_x` is walked
  like a continuous producer's — real leaves and `Real` type parameters follow
  the activation scalar — while `init_m` and `init_z` pin wholesale, mirroring
  the discrete-producer rule. Declared `Float64` initial values embed as
  zero-partial constants under non-nominal activations, which is §16.3's rule
  for condition writes applied to the defaults those conditions overlay.
  `init_x`'s walk presupposes §9.1's closed leaf vocabulary (scalars and
  `SArray`s at the common eltype), so Stratum A checks it (§14.1) and reports
  a failure in the didactic register: "`init_x` field `gear_count::Int` is not
  a continuous state — integers, `Bool`s and enums belong in `init_m`";
  "`init_x` field `q_nb::RQuat` is not a state leaf — declare the
  `SVector{4}` backing and cast where rotation semantics are wanted (§9.1)".
- **`local_types(::C)`** (concrete nominal types, same species and leaf walk as
  `output_types`): the component-local
  intermediates — fields a stage returns for the component's *own* later
  consumption (`f`, guards, the author's debug panels) that are not interface.
  `local_types(::C) = (q_dyn = Float64,)` reads exactly like `output_types`,
  and its cells follow the same tier-dependent leaf walk under non-nominal
  activations, with the same embedding guarantee (§14.5) at parametrized
  leaves. **Strict**: a
  returned field declared in neither `output_types` nor `local_types` is a build error
  with did-you-mean (§13.4 walkthrough 5). Unlike `output_types` — which any
  component that publishes must declare, publicity never being implicit — an
  empty `local_types` is the framework default — absence hides nothing, because
  strictness catches every undeclared return. No auto-publication: a `local`
  naming a state field that no stage produces is simply an error (auto-pub is
  interface sugar, and local cells are not interface). The scope these are local
  *to* is the component, not a function: computed in one stage, read by the
  component's own consumers, invisible outside — cross-stage table cells, not
  the workspace (`workspace` remains the within-call scratch, §9.3).
- **`events(::C)`**: an ordered, named collection of guard/handler pairs, tier set
  per event by the `localize` flag — `Event(guard, handler; localize = true)` is
  Tier 2, and the default `false` is Tier 1 (§2.1). Order is semantics (§5.2 declaration
  order, §11.6 once-per-event); nothing here is inferrable.
- **No stage tags anywhere.** Which stage produces which port stays invisible in
  the contract, preserving §4.2 (moving a port between stages is non-breaking
  for consumers).
  Membership is *derived* with no chicken-and-egg: stage-1 functions
  (`h_x`/`h_z`) structurally receive no inputs, so the build probes them first,
  observes their contract ports, assigns the remainder to stage 2, builds the
  graph, and probes the stage-2 chain in topological order with real upstream
  values. The "decoder takes no inputs" property is exactly what makes
  the derivation well-founded. A stage name on the wrong tier (`h_zu` on a
  continuous component) is a build error in the §13.5 family.
- **Custom structs as port types** (`contact = GearContact{Float64}`) are
  first-class under the existing §9.2 Tier-1 rules: parametric in their
  real-scalar leaves, constructors inferring the scalar, no pinned fields on
  the continuous path. The leaf walk parametrizes `Real` *type parameters*
  (recursively); a struct with a hardcoded `Float64` field offers none, so its
  cell pins and any `Dual`-carrying construction detonates inside the stage
  with an `InexactError` naming the offending constructor — the §9.2 CI
  invariant reached through the declaration layer with no extra machinery.
- **Completeness of the declaration set** — four rules the build checks in
  Stratum A (§14.1), stated here because they are properties of the
  declarations, not of the wiring. **A store needs its update:** `init_x`
  without an `f` method is a build error (continuous state with no flow — the
  framework will not silently supply `ẋ = 0`, which is a model, not a default),
  and `init_z` without a `g` method likewise (a discrete store nothing updates
  is a parameter in disguise, and parameters are plain struct fields — the
  didactic register says exactly that). `init_m` carries no such obligation:
  modes are written by handlers, and a component may legitimately declare modes
  no event of its own transitions. **An event needs both halves:** an `events`
  entry whose guard or handler has no method for the component type is a build
  error, caught by method lookup at declaration-reading time rather than as a
  `MethodError` at the first firing — an event that fires only in a corner of the
  envelope would otherwise hide the omission indefinitely. **Tier is read off the
  declaration shape:** the stage names defined (`h_x`/`h_xu`/`f` versus
  `h_z`/`h_zu`/`g`), the `init_*` family declared (`init_x`/`init_m` versus
  `init_z`), and `workspace`'s arity (`(::C, ::Type{T})` versus `(::C)`) must all
  name the same tier; a stateless leaf declares no store, so its stage names
  alone decide — which is what makes a stateless `h_xu` component tier-
  transparent library material (§15.7). Members of both families, or of neither,
  are the §13.5 kind errors. **The root of a build is an assembly:** root slots
  are the root's exported input faces (§8, §12.3) and only assemblies declare
  `exports` (§13.6), so a primitive root has no faces and every input it declares
  is an unconnected-input error. Exercising a leaf alone is what §15.7's
  component test rig is for — it supplies the one-child assembly.

### 13.3 Visibility: the contract is the interface

**Declared in `output_types` = public; declared in `local_types` = private intermediate;
declared in neither = build error; absent `output_types()` = no outputs** —
visibility by declaration *site*, the same move as
kind-by-declaration-shape. Ports in the contract are connectable, GUI-listed and
log-exported. `local_types` entries occupy table cells (they must survive from their
computing stage to `f`/guards, and they serve the author's own debug panels via
the snapshot, §12.2), but they are non-connectable — a build error, not a
discouraged convention — and presentation-filtered by default. Publicity is
never implicit: even the minimal component writes
`output_types(::LowPassFilter) = (x = Float64,)`, one line, in
exchange for "public" always meaning someone wrote it down.

- **Conformance**: a declared port or local must be produced — by exactly one
  stage, or (ports only) by **auto-publication** for declared names matching
  state, mode or `z` fields that no stage produces (§5.2). Stage membership is
  derived over `output_types` ∪ `local_types` jointly (§14.1). Declared-but-unproduced
  and produced-by-two-stages are build errors, as is an `output_types`/`local_types`
  name collision; a declared port matching neither a stage product nor a state
  field errors with both lists in hand ("not produced by any stage and not a
  state field"). A *returned* field declared nowhere is a build error at probe
  with did-you-mean against `output_types` ∪ `local_types` (under an observation-based
  rule it would silently define a new private cell instead — the
  return-side analogue of §13.4 walkthrough 1). The forgotten-branch
  walkthrough holds: a declared `P` missing from the taken branch's
  return fails at probe; missing from an *untaken* branch, it fails loudly at
  that branch's first execution via the always-on check.
- **Branch-shape rule**: stage returns must have the same `NamedTuple` shape on
  every branch — which Julia's type-stability discipline already demands for
  performance; the framework merely makes it a stated rule with a good error.
- **Schema authority is total**: every field a stage returns is declared somewhere,
  and the always-on check's expected type is fully declaration-derived — return
  typos cannot silently define new cells, and every table cell traces to an
  authored fact. Protection against silently dropped partials rests on the
  embedding guarantee — promotion is airtight, so an observed `Float64` is a true
  constant, §14.5. Probe-observed expected types are rejected on two grounds:
  authority inversion, and the fact that one probe point cannot speak for
  branch-dependent types.
- **What this rules out**: the `unlisted` flag (§4.2) — presentational
  hiding of connectable ports — and its satellite-function representation; the
  RNG-state case that motivated it needs *nothing* here (`g` reads `z` directly,
  §5.2). Identity publication by default goes with it (§9.4 step 4): publication
  driven by the contract replaces publication of everything with hiding
  annotations on top. **Probe-observed private cells** and the
  `Private(T)` fallback are rejected alongside — the former by the
  authority-inversion argument above, the latter as obviated by `local_types` (a
  wrapper inside `output_types` would break "declared = public" and introduce the
  layer's first wrapper type; a separate declaration encodes visibility by site
  instead). So is the opt-in variant with a `Float64`-under-`Dual` diagnostic —
  it legislates an ambiguity that strictness dissolves. Rows 34 and 55.

### 13.4 Failure walkthroughs (the error-locality grounding)

The five mistakes that decided declaration-vs-inference, with their failure sites
under this layer (each was traced under inference-by-evaluation too; in every case
the failure surfaced inside *correct* code, later, or never):

1. **Typo'd wire** (`:throtle`): build error at the connection, "no input
   `throtle`; did you mean `throttle`?" — vs. a missing-field error inside a
   correct `h_xu` at probe time, with the input set silently *defined* by the typo.
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
5. **Typo'd return field** (`q_dny = ...` for a declared local `q_dyn`):
   probe error with did-you-mean against `output_types` ∪ `local_types`, plus the
   unproduced-`q_dyn` error with both lists in hand — vs. the typo silently
   *defining* a new private cell under observation-authority, with the intended
   name's absence surfacing later as a missing-field error inside correct
   `f`/guard code (the return-side twin of walkthrough 1).

### 13.5 Assembly declaration: type-based, kind by declaration shape

An assembly is a plain struct: fields whose type is `<: AbstractComponent` are its
children (field names = path segments), all other fields are inert parameters;
substitutability and variants use ordinary parametric fields — exactly today's
`Cessna172X{K, A}` shape. Alongside it, well-known declarations:
`connections(::A)` (mandatory, even when empty), `exports(::A)`, `rates(::A)`.

**Container children.** A field whose type is a `Tuple` or
`NamedTuple` with *every* element `<: AbstractComponent` contributes its
elements as children, path-named `"field/1"…"field/N"` (tuples) or
`"field/key"` (NamedTuples), declaration order governing layout. Containers
are **transparent grouping, not assemblies**: no contract, no `connections`,
no rate scope, no existence beyond the path segment — the elements are
children *of the parent*, whose `connections`/`exports`/`rates` address them
by element name; anything wanting its own wiring or faces declares itself an
assembly. The payoff is parametric composition: `struct Formation{NT <:
NamedTuple}; aircraft::NT; … end` holds any roster — size, names, mixed
aircraft types — per instantiation, with the declaration bodies generating
wires by comprehension over the keys (the arity-via-computed-contracts
pattern, §6's `SumJunction{W, N}`, at structure scale); §16.9's swarm worlds
and `at("aircraft/red", problem)` mounting consume it directly. Rules: a
container mixing component and non-component elements is a build error in
this section's did-you-mean family (all-component = children;
zero-component = inert parameter data); containers of containers are
rejected in the first cut (deeper grouping is what assemblies are for);
empty containers are legal (zero children — parametric code needs no special
case); abstract element types follow the same concreteness discipline as
plain fields (directly concrete or via type-parameter bounds, §13.8's
generic holding). `rates` needs no rule change: element names are immediate
child names, hence legal keys; the bare field name is sugar for a uniform
`K` across all elements.

**The builder is rejected** (`Assembly()` + `add!`/`connect!`): the type you dispatch on and the recipe that defines its structure
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
precedent) makes an **assembly**; any leaf declaration makes a **primitive** —
`init_x`/`init_m`/`init_z`, `workspace`,
`input_types`/`output_types`/`local_types`, `events`, or any stage, `f`, `g` or
`project` method. The rule is total: a `<: AbstractComponent` type declaring
neither family has no kind to read, and is a build error naming both families
rather than a silence that fails later and elsewhere.
Reading which declarations exist is reading declarations — the same move as
§13.3's visibility-by-declaration-site, not §13.1's banned
inference-by-evaluation. Kind also settles the tier fact the `T`-signature
used to spell (row 79): the leaf walk pins or parametrizes a
producer's cells by the kind its declaration shape announces. Error
taxonomy: `connections` plus `init_x`/stages on one type is a build error
(assemblies have no state of their own — §11.5's no-atomic-assemblies at
declaration time); the neither-family error sharpens into a did-you-mean when the
type has component-typed fields ("holds components but declares no
`connections`").

### 13.6 Paths, wiring and exports

**Paths are slash-separated strings** — `"systems/ldg/left/trn"` — relative to the
assembly being declared, no leading slash; one canonical form, shared verbatim by
declarations, error messages, device/trace addressing (§12.3) and the HDF5 log
tree. Container children (§13.5) add index and key segments — `"aircraft/2"`,
`"aircraft/red"` — ordinary segments, resolved against the container field. Rejected: instance navigation (`a.ldg.left` cannot yield a *path* —
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
endpoint path — or a tuple of paths for an input face routed to several
internal endpoints (fan-out through the boundary)
(`"trn" => ("systems/ldg/left/trn", …)`). **Face names are arbitrary strings with
two build-checked invariants: no `/` (reserved for structural paths) and
uniqueness within the assembly's face set** — every other naming choice
(separators, grouping prefixes like `"pilot.throttle_axis"`) is author
convention, not framework law; the `faces` helper's defaults (§13.8) document the
house style without legislating it. The two-notation rule this rests on: **slash
is structure** (endpoint paths walking real children and ports; snapshot and log
addressing), **face names are opaque contract tokens** — and the periphery's
write side (input devices, mappings, the trace, the GUI write path) speaks face
names only (§12.3); the read side addresses snapshot cells by slash path
(§12.2).
Pairs-of-strings rather than a NamedTuple also removes the `var"..."` noise that
non-identifier names would force. Direction is derived from the endpoints
(wired-to-inputs = input face, wired-from-one-output = output face; mixed or
multi-source entries are build errors). Face *types and tiers* are derived from
the internal endpoints — §13.2's blessed derivation-from-declarations — and the
derivation is forced, not merely convenient: an assembly is tier-neutral (it
exports continuous-sourced and discrete-sourced ports side by side), so
author-declared face types would need per-face tier annotations restating what
each producer's kind already fixes (§13.5; the leaf walk reads the tier from
the producer). Rejected spellings: routing values under the leaf names
`input_types`/`output_types` (a name-level pun — a discrete leaf's exact signature with
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

**A worked assembly.** §17.5's IMU, spelled in full — a mixed-tier assembly
exercising paths, exports and rates together:

```julia
struct IMU <: AbstractComponent
    integrals::IMUIntegrals    # continuous — cumulative Θ, q, Υ, V
    sampler::IMUSampler        # discrete — integrate-and-difference latches
    errors::IMUErrorModel      # discrete — scale/bias/noise on the sample
end

connections(::IMU) = (
    "integrals/Θ" => "sampler/Θ", "integrals/q" => "sampler/q",
    "integrals/Υ" => "sampler/Υ", "integrals/V" => "sampler/V",
    "sampler/sample" => "errors/sample",
)

exports(imu::IMU) = (
    "sample"      => "sampler/sample",            # ideal increments
    "sample_meas" => "errors/sample_meas",        # measured increments (the error
                                                  # model's output port)
    faces(imu, "integrals")...,                   # kinematic-truth inputs pass through
)

rates(::IMU) = (sampler = 1, errors = 1)
```

Two spellings worth reading closely: `faces` enumerates the child's **input**
faces and nothing else (§13.8 — it is inputs-only by definition, so the
pass-through of the integrals' kinematic-truth inputs (`q_eb`, `r_eb_e`,
`ω_eb_b`, `a_ib_b`, `α_ib_b`, §17.5) takes no direction argument);
and the measured-increment face re-exports `errors/sample_meas`, the error
model's *output* port, not the `errors/sample` input the sampler already
feeds — re-exporting a wired input face would be the two-producers error of
§13.8, and would derive as an input face besides.

Three facts the example carries: the assembly is tier-neutral — every face's
type and tier derive from its internal endpoint, and a `rates` key on
`integrals`, the continuous child, would be a §13.7 build error; the two
discrete children default to `K = 1` anyway, so this `rates` declaration is
declaratory — their absolute rate arrives from the enclosing scope at
deployment (§13.7); and §17.5's latch-back wire, where the integrals consume
the sampler's published latch, joins `connections` as one more ordinary pair.

### 13.7 Rate scopes

`rates(::A) = (fcs = 1, nav = 5)` — child name => integer multiplier $K \ge 1$,
optional, unlisted children default to 1; §11.5 semantics unchanged (relative,
composing multiplicatively down the tree, compiled to absolute divisors). Keys are
**immediate child names only** — a deep key would edit another type's design from
outside, and the composition rule guarantees you never need to. Container
elements (§13.5) are immediate children, so `"aircraft/red"` is a legal key;
the bare field name applies one `K` to every element. `K` on a
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
    names = input_faces(child)            # keys(input_types(c)) for a leaf,
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
because `exports` is ordinary code (map over the pairs); normative signatures
for both primitives in §15.3. Every error stays
first-class: an `except` face the assembly then fails to wire is an ordinary
unconnected input; a face both wired and re-exported is a two-producers error;
`except`/`only` naming a nonexistent face errors with the child's face list in
hand; a `prefix = ""` collision is caught by the build's uniqueness check like any
hand-written duplicate. The effective export list is plain printable data — the
inspectable contract of this instantiation. What computation does *not* do is
auto-bubble: the author wrote down "every input face of this child that I don't
feed, I expose under this prefix" — explicit at the type level, evaluated at
build.

**`faces` is inputs-only, by definition.** It reads `input_faces(child)`, and
`except`/`only` filter *face names* within that set — there is no direction
argument, and none is implied: the helper exists for the pass-through case,
where an assembly hands a child's unfed requirements up one level. Computed
*output* re-export — a `direction` keyword, or the predicate selection §15.7
records for `except`/`only` — is a **guarded addition**: plausibly wanted by the
conventional-surface work (§12.2, §19) and by test rigs, cheap to add, and not
adopted here, because every output face in the worked assemblies is an explicit
pair and no consumer has yet demonstrated the computed form.

**Generic holding = imposed contract.** A parent holding a child generically
constrains it exactly through the faces its wires and exports reference: build a
`World` whose concrete aircraft lacks a referenced face and the error names the
`World` entry — build-time structural typing, no new vocabulary (a formal
required-faces declaration on domain abstract types remains possible sugar).
Scalar faces make partial scripting compose: a guidance scenario component wires
`mode_req` and `EAS_ref` while the remaining faces stay exported for GUI or
defaults — impossible with a bundled face (§4.3 write-side corollary).

---

## 14. The build pipeline

The build consumes a root component instance and produces the runnable artifact:
resolved wires, typed signal table, evaluation schedule, absolute rate divisors,
flat state layout, root slots. §13 states what is declared and what must hold;
this section states *when* each fact is checked, against what, and with which
failure. The §13.4 walkthroughs plus §8's error rules are its acceptance tests.
Error-*reporting* policy is settled in §15.1: declarative checking passes
batch, user-code evaluation fails fast, strata are barriers — the only partial
results carried past failures are violation lists from pure checks.

### 14.1 Three strata

Three ordering constraints are forced by settled decisions: face derivation is
**bottom-up** (an assembly's `exports` evaluates against child contracts,
§13.8); the unconnected-input obligation check and cross-level two-producers
detection are **global** — decidable only at the root, after every assembly's
wires and exports are in hand (§8); and stage membership is **derived by
probing** the stage-1 functions (§13.2), so evaluation interleaves with graph construction at
exactly one blessed spot. The pipeline is therefore inherently heterogeneous,
organized as three strata:

- **Stratum A — structure.** Pure declaration reading; no user stage code
  executes (`exports`/`faces` bodies are declaration code, §13.8). Tree walk
  from the root instance: components by path, kind classification by
  declaration shape (§13.5), leaf contract collection (`input_types`,
  `output_types`, `local_types`, `init_*` values, `events`), bottom-up face
  derivation, then global wiring resolution to absolute leaf terminals:
  one-writer-per-input, typo did-you-mean against the destination's input
  list, the wiring bound check at nominal faces (`producer_face <: entry`,
  §13.2 — equality the concrete degenerate; abstract-at-root detected here),
  the whole-tree obligation
  check, root slots falling out as the root's exported input faces, and §9.1's
  closed leaf vocabulary checked on every `init_x` (the walk in §13.2 rests on
  it). Also here: §13.2's declaration-completeness rules (a store without its
  update, an event missing a guard or handler method, a leaf mixing tier
  families, a primitive at the root), and
  `rates` validation and compilation of relative multipliers into
  absolute divisors — everything except binding `Δt_base`, which is
  deployment.
- **Stratum B — schedule.** The single evaluation-feeds-structure step:
  stage-1 probes at `Float64` on `init_x`/`init_m`/`init_z` values (well-founded — the
  no-feedthrough stage takes no inputs), port classification (stage-1 /
  auto-published / stage-2 remainder), feedthrough graph from wires carrying
  stage-2 ports, topological order, §5.4 cycle rejection. The output is
  structural: names only, `T`-independent, branch-protected by the
  branch-shape rule plus the always-on check (§14.5).
- **Stratum C — activation, parametric in `T`.** Everything type-shaped:
  declared producer types (`output_types`/`local_types`) leaf-walked at the
  activation's `T` to type the cells, the probe chain run in topo order, observed compared against
  declared, flat `x` buffer and table laid out. The nominal `Float64`
  activation runs at build; other activations re-run *only this stratum*
  (§14.4).

Deployment binding (`Δt_base`, `h`, `n`, `t_end`, algorithm, harmonic-grid
validation, tick schedule instantiation) sits after all three, at `Simulation`
construction — nothing in A–C depends on it.

### 14.2 The `Build` artifact

`build(world) → Build` is a standalone entry point; `Simulation(world; ...)`
(§17.4's spelling) is the convenience that calls it and adds
deployment binding, buffers and the stopped-sim services. The `Build` is the
inspectable contract of the instantiation §13.8 gestures at — wire list, face
table, schedule, root slots as plain printable data. CI checks a model by
calling `build`; the acceptance tests target `build` errors directly;
`attach!` validates device bindings against it. Rejected: build living only
inside the `Simulation` constructor — CI would construct simulations with
dummy deployment parameters, and the phase outputs must exist as coherent
intermediate data anyway; the artifact just names them.

### 14.3 Probing and input synthesis

**Probe-everything scope.** The nominal activation probes every user function —
stages, `f`, `g`, guards, handlers, `project` — once, at the initial state,
with real values, checking shape/type conformance and discarding results (all
are pure; the cost is one evaluation each). "Fails loudly at build time where
possible" (§13.1) decides this: a malformed `f` return must not wait for the
first integrator step. Probes see only the initial state's branch — the
marginal coverage is earliness, not completeness; the always-on check (§14.5)
remains the completeness backstop.

**Probe argument sourcing.** `x`/`m`/`z` come from `init_*` declarations
(declared by value); `y_x`/`y_z` from the stage-1 probes; wired inputs from
upstream products, real values available because the stage-2 chain is probed in
topological order. Exactly one kind of terminal has no producer: **root
slots**. The build synthesizes their values via `probe_value(::Type)`:
framework methods for `Real` (`zero(T)`), `Bool` (`false`), enums (first
instance), ultimate fallback the zero-argument constructor `T()` — which is
where well-behaved constrained types already put their valid default (`RQuat()`
= identity; the `@kwdef` convention supplies it broadly). `probe_value` is
**overridable**: a type whose valid default is not reachable that way declares
its own method, which is also the seam a Tier-1 type uses to state a
constrained default. No method → build error, in the didactic register: it
names the face and the type, and asks for one of the two fixes ("no
`probe_value` for `Ranged{Float64, -1, 1}` at face `pilot.elevator_axis` —
define `probe_value(::Type{Ranged{Float64, -1, 1}})` or a zero-argument
constructor"). Synthesis never meets an abstract type:
root slots are concrete by §13.2's tight-bound rule, and abstract entries only
occur on component-fed inputs, which the probe sources from upstream products.
Physically silly values are acceptable by
construction: the probe checks types, and return types that depend on input
*values* are type instabilities (banned by the branch-shape rule); the §4.3
write-side granularity corollary keeps root slots predominantly scalar, so the
surface is small. Rejected: inputs declared by value à la `init_x` (reads as
an unwired-input default, §8's banned concept; every leaf pays for a
root-only need; and it would need its own agreement rule for fan-out through an
exported face, where the slot type already follows §13.2's rule — the unique
concrete declaration among the consumers); NaN poison
values (`Bool`/enums unpoisonable; probe values are meant to be read — poison
guards illegitimate reads, §9.3); init-service values (the build is
standalone; services post-date it).

**Probe values are strictly probe-scoped.** Everything the probe writes is
garbage once the build finishes; probe values never double as initial slot
values — that would smuggle in the default semantics rejected above. The
same doctrine covers the clock: `Δt` in seconds does not exist until
`Simulation` binds `Δt_base` (deployment post-dates the build), so
discrete-tier probes supply a placeholder period (`1.0`) in the bundle —
a fabricated, probe-scoped value like any synthesized input; the probe
checks types, not physics.
`Simulation` must not reach its first boundary with uninitialized slots;
enforcement is the pre-write `UninitializedSlots` check at `init!`/commit
(§16.6).

### 14.4 Activations: executable sets, laziness, caching

An **activation at `T`** re-runs Stratum C with a different scalar: table cells
re-typed by the leaf walk over the declared types (§13.2 — tier-dependent on
producer-fed cells, over the consumer declaration on root slots), table and
state buffers re-laid-out, workspace allocators re-invoked at `T` (a
continuous component's scratch carries the activation's scalar, §9.3), probe
chain re-run. Structure and schedule are `T`-independent by construction.

**Each activation probes exactly the function set it can execute.** A `Dual`
activation (linearization, gradient trim) evaluates the model at a frozen
instant: discrete stages are gated off holding `Float64` values (the §13.2
frozen-constant semantics), guards and handlers never run (event localization
is `Float64` sweeps by design, §11.4). Only the continuous output stages
(`h_x`/`h_xu`) and `f`
ever see a `Dual` — so only they are probed. Probing the discrete stages, `g`, or guards at `Dual`
would check code against a number type it cannot receive. One rule, no
special cases; the §10 tracer activation follows it identically. "Tracer
activation" names row 12's *global* set-tracer — a whole-model run at the
tracer scalar, an activation like any other. §10's cycle classifier is row
12's other variant, the schedule-free per-member local trace, which runs in
Stratum B's failure path and is not an activation at all.

**Lazy, with an opt-in exhaustive mode.** Non-nominal activations run at first
request, not at build: the dominant cost is compiling the continuous chain a
second time at `Dual`, pure waste for interactive fly-around use. The price,
stated openly: `build` succeeding does **not** certify the model linearizable —
a pinned `Float64` (§9.2) lurks until the first `Dual` activation detonates it
at the probe, naming the offending constructor. The repository's test suite
pins the invariant instead: `build(world; activations = (Float64, ProbeDual))`
(or a `check` entry) runs the exhaustive set in CI, catching Tier-1 genericity
violations at PR time. `ProbeDual` is the framework's exported canonical probe
scalar — `const ProbeDual = ForwardDiff.Dual{ProbeTag, Float64, 1}` — because
an activation is keyed by a *concrete* scalar type and the bare `Dual`
`UnionAll` cannot key one, be walked to, or answer `zero(T)`. Its width is
arbitrary: what CI pins is genericity, not any particular Jacobian, so one
canonical width suffices even though §16.10 chunks at whatever widths it needs.

**Caching is implementation detail, not semantics.** An activation is a pure
function of the build and the concrete scalar type; the `Build`/`Simulation`
holds (layouts, buffers, validated-flag) keyed by that type. Compiled code is
cached by Julia itself, process-wide; what the framework cache saves is probe
re-runs and re-allocation — which matters exactly in activation-reusing loops
(the envelope-grid gain-schedule case: hundreds of trim-then-linearize points
on pre-existing buffers, zero-allocation). Nothing numerical is ever cached.
Note `Dual{Tag,V,N}` carries the partial count: a different seeding width is a
different scalar type, hence a separate entry and a separate Julia compile.

### 14.5 The always-on conformance check

The probe validates each function *once*, on the initial state's branch; the
schema-authority bargain's second clause ("at first execution otherwise",
§13.1) is discharged by leaving the probe's comparison permanently in place.
At the point where the executor stores a stage return into the table, it holds
the complete expected return type at this activation — the declared types of the
names *this stage* produces, as fixed by Stratum B's stage classification: one
concrete `NamedTuple` type per (component, stage), drawn from `output_types` and
`local_types` alike. Auto-published names belong to no stage's expected type;
the framework writes those cells itself. The executor performs a single
type test against that type (conceptually `y2 isa Expected`). Type-stable conformant
code: the compiler proves the return type, decides the test at compile time,
and deletes it — zero instructions. Branch-divergent code: the test survives
as a runtime check, nanoseconds on conformant branches, a loud located error
on the divergent one at its first execution. (Type-unstable-but-conformant
code pays nanoseconds on top of the dynamic dispatch it already bought.)

**Exact match at nominal; embed-accept at parametrized leaves.** At
the nominal activation — the only one that ever runs in real time — the check
is unchanged: exact type match, no convert-on-write, one baked `isa` that
folds away (`Int` sloppiness must fail at `Float64`, not lurk until a `Dual`
activation makes "it runs" activation-dependent). The error can afford to be
didactic: "field `M_shaft`: expected `Float64`, got `Int64` — return
`zero(x.ω)`, not `0`". Under a non-nominal activation, a parametrized leaf
accepts exactly two types: the activation scalar (the fast path — the baked
`isa` unchanged) or `Float64`, which the executor **embeds** as a
zero-partial constant (`convert` through the leaf; struct-valued ports use
the standard cross-eltype constructor, a missing one failing loudly with
both types named). Nothing else is accepted. The embedding is exact, not
lenient: promotion is airtight and there is no lossy `Dual → Float64` cast,
so a `Float64` observed at a parametrized leaf means no `Dual` entered its
computation — its true derivative along every seeded direction is zero,
which is precisely what the embedded constant says. This scopes the blanket
convert-on-write rejection to the nominal check (row 53):
the bug that rejection guards against — silently zeroed partials — cannot
arise from honest code, because accidental `Float64`s from `Dual` operands
are impossible (`MethodError` at the operation site). The residual is
**deliberate stripping** (`ForwardDiff.value`): a stated intent to discard
partials, producing a silent zero in the Jacobian — the stop-gradient idiom,
occasionally legitimate (deliberately frozen couplings, opaque non-Julia
wrappers), invisible to the schema, and equally invisible to a strict
exact-match rule when applied mid-expression, so the leniency costs nothing.
Schema-visible freezing is a recorded door (§16.10).

**Uniform across all probed functions.** `f` checks against `X`'s own shape
at the activation's `T` (§9.1: a scalar leaf expects a `T`, an `SArray` leaf
the same `SArray` at `T`), the predicate being "every
field scatters into its field's block at `T`", which is what makes derivative
completeness structural rather than a matter of author discipline; guards
against their probe-derived condition form (below); `g`
against the `z` shape; handlers against the §5.2 return law, key by key.

**Handler returns, key by key.** The returned NamedTuple's key set is checked
first: an unknown key, or a key naming a store the component does not declare,
is a build error with did-you-mean against `{x, m}` narrowed to the stores that
exist — the bundle law's classification running in the return direction. Then,
per present key: `x` must be **complete** against the state field set (and
conformant at `T` like any state value), while `m` may be **partial**, checked
against a names-subset-with-matching-types predicate — still a type-level
computation that folds when inferred. An absent key is not an error and not a
no-op to diagnose: it is the handler saying it does not touch that store. The
completeness asymmetry is storage-shaped: `x` must be complete because it lives
in a flat buffer written back wholesale, while `m` may be partial because
`m` lives in per-field cells where a partial merge is the natural write.

**Guards have two admissible forms** (§2.1), so their check is form-aware
rather than a flat `isa Bool`: a predicate guard's probed return is `Bool`, a
continuous guard's is the nominal scalar — guards run only at the nominal
activation (row 52), so no parametrized-leaf case arises here. Any other probed
return type is a build error naming both admissible forms. So is a Tier-2 event
whose guard probes `Bool`: localization brackets a root, and the predicate form
offers none — return `x.ω - eng.ω_idle`, not `x.ω > eng.ω_idle`.

**Failure payload:** component path, stage, field-level diff (missing /
unexpected / per-field expected-vs-observed), simulation time. Deliberately
absent: the source branch (values carry no provenance; the diff identifies
it). The always-on input trace makes every such failure **reproducible by
replay** — the error names the boundary to replay to. At run time the failure
travels as a species of `StepError` through the single catch site (§15.4),
which adds the loop-level nonfinite-state check as its divergence sibling.

### 14.6 Stopped-sim services as Stratum-C clients

Sketched here because it grounds the strata; the services themselves are §16.
The C172 trim problem (`c172.jl`: `TrimState`, `TrimParameters`,
`θ_constraint`, the `ẋ`-reading cost) transfers near-verbatim:

- **Trim** is a write-condition → sweep → read loop on an activation — by
  default the `Dual` activation, decision variables seeded for exact residual
  Jacobians (§16.7); the derivative-free fallback runs the same loop on the
  nominal `Float64` activation with no new activation needed, and the
  always-on checks ride along either way. Decision variables stay opaque to the
  framework (only the assignment's *output* is framework vocabulary).
  `assign!` inverts from in-place mutation + self-invoked `f_ode!` to a pure
  function returning a condition value (state by path, modes, slots by face)
  that the service writes and evaluates. Domain math — the pitch constraint,
  `Kinematics.Initializer`, per-residual scalings and the equilibrium-subset
  choice — survives untouched, aircraft-side.
- **Linearization** is a `Dual` activation plus seeded sweeps: gather/scatter
  over the canonical layout replaces the hand-written
  `get_x_ss`/`assign_x_ss!` layer (§9.1's deletion discharged); root slots are
  the input surface; frozen discrete outputs are constants with zero partials,
  which is exactly "linearize with `z` held" (§13.2). Gradient-based trim —
  decision variables seeded through the `T`-generic assignment math — is the
  default (§16.7).
- The generic service loop (vectorization, optimizer setup, bounds packing,
  solved-condition write-back including root slots and the trace header's
  slot capture) replaces today's per-aircraft NLopt plumbing. A failed trim
  leaves the simulation's stores untouched — an improvement over today's
  warn-but-assign `f_init!`.

### 14.7 The compiled executor

The schedule exists in two representations at two lifecycle stages. In the
`Build` it is plain printable data (§14.2) — paths, stage names, order — the
authoring and diagnostic form. At `Simulation` construction, and per
activation, that data is compiled into the execution form: **a
concretely-typed tuple of entries over statically typed cell storage,
traversed by a compile-time-unrolled walk**. This is a forced move, not a
preference: §9.5's zero-allocation invariant, §14.5's fold-away conformance
test and §5.1's zero runtime graph logic are reachable only under full
specialization. A heterogeneous runtime vector of entries dispatches
dynamically per entry and boxes every stage return — the framework itself
would allocate, and §9.5's canary role (framework contribution exactly zero,
so any model inference regression announces itself) would die with the
invariant. An entry carries what selects code — component type, stage — in
type parameters, and what is plain data — tick divisor, layout offsets — in
fields; gating compiles to an integer test inside the specialized body.

**Phase bodies are the outer decomposition, and they are semantically
forced.** The boundary sweep's `h_x` block is order-free by definition (the
no-feedthrough stage reads no `u`); the `h_xu` block, with due discrete
stages gated in, is the only topologically ordered one; the `f` block — the
RHS body the stepper calls per stage evaluation — and the `g` block are
order-free with disjoint writes; guards and handlers are their own small
callables inside the §11.6 iteration. These bodies run at different times
and communicate only through the stores and the table, so the seams cost
nothing — no values cross them. Two doors this structure opens for free,
recorded not committed: deterministic parallel evaluation of the order-free
blocks (disjoint writes, and no floating-point reductions to reorder — §6
made every sum an ordered junction entry), and finer recompilation
granularity (editing a discrete component invalidates the boundary body, not
the RHS body).

**Chunking bounds the compile cost.** Within a large block the tuple splits
into chunks behind non-inlined but statically-typed function barriers.
Inside a chunk everything the design relies on survives — static dispatch,
inlining, view SROA, check folding, zero allocation; at the seams only
cross-entry fusion is lost, which a table-mediated dataflow barely had.
Chunk size is the implementation's *only* representation freedom (fully
fused and chunk-of-one are its endpoints), and it converts the compile cost
from superlinear in the largest body to linear in entry count. Measured
anchors (2026-07, synthetic ~15-op bodies, Apple Silicon): a fused 400-entry
sweep compiles in ~0.8 s against ~0.34 s chunked at `Float64`, the fused
curve visibly superlinear; an 8-partial `Dual` activation multiplies
instruction count ~20× and lands near ~9 s for 400 chunked entries — linear,
instruction-bound rather than structure-bound. Extrapolated to a
C172X-scale model (§17.4; roughly 200–400 entries, larger bodies), the
nominal activation sits at seconds, the `Dual` activation in the tens before
mitigation; re-measurement on the real vehicle skeleton is a §19 migration
item.

**The mitigation ladder**, in order: activations are lazy (§14.4 — a session
that never linearizes never compiles `Dual`); non-nominal activations may
compile at reduced optimizer level (their sweeps run inside service loops
where microseconds are irrelevant — a one-line per-module policy); and
activations bake into package images via ordinary precompile workloads — an
aircraft package exercising build-plus-one-sweep per activation turns TTFX
from a session tax into a CI artifact.

**Views are spelled rebuild-per-call.** Every entry constructs its bundle at
its own position; there is no framework-maintained hoisting and therefore no
cache-invalidation obligation. Hoisting belongs to the code generator: CSE
merges repeated loads exactly where no intervening store invalidates them —
which is precisely the staleness rule — and the sweep-varying bundle fields
(`u`, `y_x`) are per-call by topological necessity either way (§9.1).

**Construction is type-opaque; only the executor specializes.** Schedule
tuples are built from untyped buffers and splatted once. Generic tuple
utilities — range indexing, long `ntuple` closures, naive recursion — are
inference traps at schedule length (a 400-entry heterogeneous tuple can send
generic `getindex` inference into combinatorial collapse), so the compiled
tuple's type has exactly one consumer: the unrolled walk.

Rejected (row 86): vector-of-abstract entries (dynamic dispatch boxes
returns; union splitting rescues only toy scale and cannot close an open
method set); type-erased call tables (`FunctionWrapper`-style — the same
specialization count as chunk-of-one with extra machinery, no cross-entry
optimization, and a load-bearing seam on internal ABI, the §11.1 lesson);
framework-maintained view hoisting (lifetime obligations the compiler
already discharges, where mis-scoping means silent stale-state reads).

---

## 15. Error discipline

§13.4 fixed what must be caught and where; §14 fixed when each fact is checked.
This section fixes how failures are *reported* — the reporting policy, the
diagnostic representation, the runtime failure story, and the seam between "the
model reached a terminal state" and "the run should end". Two of FlightCore's
paid-for lessons ground it: the compact-backtrace discipline (parameterized
model types make rendered output unreadable) and the `SimulationTermination`
machinery, which §15.5 replaces.

### 15.1 Reporting policy: batch the checks, fail the evaluations fast

The fail-fast vs. compiler-style question dissolves once the build's failure
sites are split into their two populations:

- **Declarative checks over collected structure** — unconnected inputs,
  two-producers, wire typos and type mismatches, face-name uniqueness,
  `output_types`/`local_types`/state-field consistency, `rates` validation. Each is a
  pass over a list; the whole-tree obligation check literally computes *the
  set of* inputs whose obligation chain never terminates. Reporting every
  violation is the natural output of the pass — truncating to the first would
  be extra work — and these failures cluster in practice (a freshly written
  assembly has five unwired inputs; a renamed port breaks three wires).
  **These passes batch:** each returns its full violation list.
- **User-code evaluation** — `exports` bodies (Stratum A), the stage-1 probes
  (B), the probe chain (C). When user code throws there is no meaningful
  rest-of-batch: a failed `exports` leaves the parent's face derivation
  undefined; a failed stage-2 probe starves every downstream probe of its wired
  inputs (probe values flow topologically, §14.3). Continuing past these
  requires genuine compiler machinery — poisoned nodes, cascade suppression,
  dependent-check skipping — and buys little, because user-code failures are
  typically singular. **The first user-code exception aborts the phase.**

Strata are barriers: a stratum that produced any error-severity diagnostic, of
either kind, throws before the next stratum begins — probing against
unresolved wiring is meaningless. The only partial results ever carried past a
failure are violation lists from pure checking passes, so the cost that kept
this decision open — build phases having to carry partial internal results
across a failure, the machinery §14's strata would otherwise need — never
materializes.

**No cascade suppression within a stratum** — a deliberate simplification. A
typo'd wire (`:throtle`) produces both the did-you-mean error and an
unconnected-input error for the intended `throttle`; both are reported. They
render adjacently (diagnostics sort by path), the pairing is self-explanatory,
and suppression heuristics are exactly the fussy machinery this split avoids.

### 15.2 Diagnostics: structured values, one carrier exception

A diagnostic is a plain value from a **small closed set of kinds**, enumerated
normatively in **Appendix C** — kind name, payload fields, owning section,
severity — which is the artifact the §13.4 acceptance tests and the
error-message work are written against. Each kind carries its structured payload: endpoint paths, face names,
expected/observed types, the *list-in-hand* for did-you-mean rendering, and a
severity. Checking passes return diagnostics; the stratum barrier throws a
single `BuildError` wrapping the batch; `showerror` renders it compiler-style,
grouped by kind and sorted by path. A user-code exception is wrapped in a
framing diagnostic — component path, which function, the probe context
including synthesized inputs — with the original exception as `cause`, so the
didactic frame renders first and the raw throw second. One class is recognized
rather than merely framed: a `FieldError` — type and field carried as data —
matched against the bundle's own NamedTuple type becomes §5.2's bundle-law
did-you-mean, with the legal set and the undeclared-store / wrong-tier /
illegal-for-this-function classification. Nothing is recovered by reading
message text.

The §13.4 walkthroughs as acceptance tests target diagnostics: tests match on
kind plus payload fields, never on message text. Messages become pure
presentation.

Two rendering rules are doctrine, not style:

- **Strings, never instances.** Diagnostics carry paths and names as strings —
  never component instances, never model types (the `compact_backtrace`
  lesson). Expected/observed *port* types are the payload exception, and they
  are small: `Float64` vs. `Bool`, a NamedTuple field diff.
- **The didactic register is policy.** Every diagnostic states the fix or the
  lists-in-hand, not just the violation: "return `zero(x.ω)`, not `0`"; "no
  input `throtle`; did you mean `throttle`?"; the child's face list alongside
  the unknown `except` entry.

**Two warning streams, scoped separately.** The *build* diagnostic stream is
the one Appendix C's build kinds ride: warnings there carry warning severity,
render with the batch, and never trigger the throw. That stream's warning set
is **currently empty** — the unconnected-output warning having been rejected as
its sole candidate (§8, row 84) — and better an empty, trusted stream than a
noisy one; a warnings-as-errors CI switch is addable, not built. The *runtime*
status/log stream is a different channel and is not empty: it is
per-occurrence, rate-limited wherever its source can repeat, and surfaced
through the published framework status (§12.2) alongside the §11.7 pacer
diagnostics and the §12.7 liveness heartbeats — never batched, since there is
no batch to join. Nothing in row 84's argument applies to it: that decision is
about what the *build* warns on. The committed runtime warnings, in one place:

- **chattering / localization-budget exhaustion** (§11.4) — a Tier-2 event
  whose bracketing budget runs out at a boundary;
- **event re-enabled within a boundary** (§11.6) — deferred to the next
  boundary by the once-per-event rule;
- **forgiven-debt re-anchor** (§11.7) — the pacer abandoning accumulated debt
  and re-anchoring its schedule;
- **the two surface violations at the drain** (§12.3): `OutOfClaimEntry` (an
  enumerated surface's binding drift) and `StaleInteractiveEntry` (the derived
  surface moved between staging and drain, §12.5);
- **thread-budget tightness** (§12.7) — at `run!` and on each attach, as the
  roster grows;
- **device join timeout** (§12.9) — a device task exceeding the shutdown
  join timeout, abandoned by name rather than hanging `run!`;
- **device crash** (§12.9, §15.4) — a device task's failure caught by the
  framework wrapper, the sim continuing with the device absent;
- **workspace poison skip** (§9.3) — once per activation, naming the stores
  with no usable sentinel;
- **unbounded run** (Appendix B) — no finite `t_end`, no `stop_on` faces,
  `pace = Inf` at run start.

### 15.3 Build primitives: `resolve` and the face-list accessors

The §13.8 sketch's primitives, made normative:

- `resolve(asm, path::String) → AbstractComponent` — the getfield walk along
  `/`-segments. Its one non-obvious duty is enforcing §8's generic-boundary
  rule: it walks *declared field types* alongside instances, and a segment
  that traverses **past** a generically-held field (non-concrete declared
  type) is a diagnostic even though the concrete instance in hand would
  resolve it — resolving *to* a generic child is port-level access and legal.
  An unknown segment errors with the sibling field list in hand.
- `input_faces(c)` / `output_faces(c) → Vector{String}` — the keys of
  `input_types(c)` (stringified) for a leaf; the input/output entries of
  `exports(c)` for an assembly. Declaration order is preserved: deterministic
  printouts, stable diagnostics.
- `resolve_terminal(asm, path) → (component, name)` — splits a terminal
  path's final segment (unambiguous: face names may contain dots, never
  slashes, §13.6) and resolves the prefix through `resolve`. First-class
  because five clients share it — wiring resolution, condition addressing
  (§16.3), device-binding validation (§12.2), surface resolution (§16.10),
  snapshot inspection — one splitter, one did-you-mean site.

### 15.4 Runtime failures: one catch site, an execution cursor

**Where caught.** The loop wraps each execution of the boundary macro-sequence
(integrate → project → event iteration → ticks → publication) in a single
`try` — never per stage or per component, which would salt the hot path with
exception frames for no benefit. Framing information does not need to be
*caught* into existence: the executor maintains an **execution cursor**, a
plain mutable field in the loop state recording where in the compiled schedule
it is — component path (schedule index), which function
(`h_x`/`h_xu`/`h_z`/`h_zu`/`f`/`g`/guard/handler/`project`), and the boundary phase
(integration stage *k*, event round *r*, Tier-2 localization probe at trial
time, tick). One cheap store per dispatch on a single-tasked executor — no
allocation, no exception frames — and it covers every user-code surface
uniformly, including the forgettable ones: RHS evaluations at interior RK
stage points, guard evaluations inside ITP/Brent probes, environment closures.

**How handled.** The catch site wraps the original exception in `StepError` —
the runtime counterpart of `BuildError` — carrying the cursor's frame, the
boundary time, the trace boundary index (the replay pointer), and the original
exception as `cause`, rendered with compact frames per §15.2's doctrine. The
§14.5 conformance failure needs no separate path: it is thrown as its typed
diagnostic at the table-write point and arrives at the same catch site, a
species of `StepError` with the field-diff payload. Reproducibility holds by
construction: staged inputs are drained and recorded to the trace at the frame
top, *before* the boundary executes, so the failing boundary's inputs are
already in the trace when it fails — the error names the boundary to replay
to.

**Disposition.** The `Simulation` ends in a terminal status — `stopped` vs.
`errored` — with the exception retrievable. A synchronous batch run rethrows
after the shutdown tail completes, so CI fails honestly; an interactive
session logs the rendered error and surfaces the status through the control
plane and GUI.

**The nonfinite check.** Divergence is not termination: dynamics that blow up
(ground penetration, an unstable gain) produce NaNs that defeat guards — NaN
comparisons are false — so no declared condition will catch them. A loop-level
`isfinite` sweep over `x` at boundaries fails fast as a `StepError` species
naming the offending component's state block and the boundary. It catches
diverging models generally, not just post-terminal ones.

**Domain separation.** Device-side user code — mappings run on the device
task (§12.3) — fails in the device's own domain and takes the settled
per-device crash path (`should_close`, liveness heartbeat, `DeviceCrash`);
the sim keeps running. The two failure domains never mix — exactly what the
no-shared-mutable-model decision bought.

### 15.5 Termination is a state, not an exception

FlightCore's `SimulationTermination` idiom — model code throws, the loop
catches and logs it as informational — has **no counterpart here**. It sits badly
in this design: a mid-sweep throw aborts a boundary halfway, and §12.9 is built on
completing one. The discipline: **exceptions from model code are always
abnormal**; graceful termination is model *state*, reaching the loop through
declared machinery:

- **Detection** is ordinary guard/handler/mode machinery. Declare the
  condition as a Tier-2 event if the stop should be localized — touchdown
  overload is precisely a zero-crossing: the boundary is localized to the
  crossing, the handler sets `m.crashed`, and the snapshot at the crossing
  instant carries the touchdown state.
- **Publication** is an ordinary `Bool` output face, exported to the root.
  Within concretely-declared structure, deep wires gather the condition at its
  owning boundary in one visible block (§8): `Ldg` ORs its three legs through
  a junction (§6's ownership idiom, §15.7's library) and exports one `damaged`
  face; intermediate assemblies are untouched. Each *generic* seam costs one
  export entry — and that hop is the substitutability contract doing its job,
  not plumbing (§13.8's imposed contract).
- **Policy** binds at deployment: `Simulation(world; ..., stop_on = (...))`
  names root-exported `Bool` output faces, OR-combined, validated against the
  `Build`, recorded in the run metadata. After *every* published boundary —
  grid, `t*` (§11.4) and boundary zero (§16.5) alike — the
  loop reads the named faces in the snapshot it just published; the first
  `true` initiates §12.9 shutdown with *this* snapshot as the final one — the
  terminal snapshot is the terminal state, no roll-back, nothing §12.9 doesn't
  already do. `run!` therefore checks the boundary-zero snapshot before the
  first step: an authored condition already terminal ends the run at `t₀`
  with that snapshot final, integrating nothing. Default: no stop faces, run
  to `t_end` — `stop_on` is `t_end`'s model-condition sibling at the same
  declaration site.

**Both are `run!`-time overridable, with the constructor value as the
default.** `Simulation(world; t_end, stop_on)` sets the defaults for the
simulation; `run!(sim; t_end = …, stop_on = …)` binds them for **that run
only** — the `run!` argument wins where given, the constructor's value stands
where it is not, and nothing about the `Simulation` is mutated, so the next
`run!` without arguments gets the constructor's policy again. The effective
values are what the run metadata records (unchanged obligation, now stated of
the effective pair). `stop_on` face validation against the `Build` runs at
**both** binding sites, identically: an unknown or non-`Bool` face fails at
`run!` exactly as it fails at construction. This is not the root-declared stop
policy rejected below: binding moves one notch *later* along the same axis —
more deployment-flavored, not less — and §12.11's `stopped → init! → run!`
cycle and `step!` register are precisely where one `Simulation` wants different
stopping policies on different runs. The honest cost is two homes for one fact,
which the precedence rule above is what settles.

Taught contract: **stop faces are sampled at completed boundaries; declare an
event if you need the stop localized.** Both condition shapes work without
framework latching — a handler-set `m` flag is sticky by nature, a transient
stage-2 Bool is caught because the loop reacts to the first `true`. Compound
conditions compose in-model — a monitor component reading the relevant
signals and outputting one Bool — the same move §12.10 made for scripts.

Post-terminal dynamics are the model's job, and that is a feature: today
`robot2d` *throws* when it falls because it has no other way to say "my
dynamics are no longer meaningful"; here it declares the fall as an event,
switches to a frozen mode (mode-dependent `f` — machinery it already has), and
exports `fallen`. Wired, the sim ends at the fall; unwired, it integrates a
frozen robot — well-defined, unlike an uncaught throw. The discipline forces
models to have well-defined terminal states, which is better modeling.

Rejected mechanisms:

- **Predicate closures** (`stop_when = snap -> ...`): opaque (not printable
  data — against the `Build`-as-inspectable-contract doctrine), not
  serializable into run metadata, a permanent public snapshot-reading API for
  user closures, and every use of the extra expressiveness is logic that
  belongs in-model. `user_callback!` in a costume.
- **Root-type-declared stop policy**: the rates precedent read correctly cuts
  the other way — ratios travel with the design, *absolutes* bind at
  deployment, and "stop here" is absolute-flavored run policy: development
  runs past the condition to inspect, batch studies log it and continue,
  services don't step at all. (A root-declared *default* overridable at the
  ctor remains the one variant worth reopening if migration shows the ctor
  argument chronically forgotten.)
- **Blessed terminal types/names scanned from the tree, and `terminal` event
  flags**: action at a distance — a deep declaration halts the world, the
  root contract says nothing, substituting an aircraft silently changes when
  runs end, and per-deployment disabling needs masking machinery. The
  localization that terminal events appear to buy structurally is already
  available under `stop_on` via the event idiom.
- **Control-plane capability for components**: §12.6 separates the control
  plane from model semantics precisely because components live inside
  boundary semantics; a mid-sweep imperative control action is the
  `user_callback!` shape again. The device-capability precedent does not
  transfer — devices live outside the boundary loop.
- **Observation-by-path** (`stop_on` naming a deep path into any public
  output): rejected on a line worth recording as doctrine. **Diagnostic
  observation** — the log retaining the full table, GUI panels rendering a
  component's ports, replay inspection — is human-facing, has no effect on
  run semantics, and legitimately sees every public cell. **Load-bearing
  observation** — a read that changes what the run *does* — must speak the
  contract. A deep `stop_on` path makes termination semantics depend on
  internals no contract mentions: precisely the knowledge violation §8
  forbids for wires (which bans deep paths into generic children *even where
  the concrete instantiation would resolve them*), brittle under substitution
  for reasons the aircraft's contract never promised anything about, and it
  leaves the root face list lying about what can halt a run. Output devices
  are not a counterexample but the other half of the same doctrine:
  their reads are diagnostic — snapshot-path bindings, §12.2/§17.4 — while
  `stop_on` is the one read that changes what the run does, which is exactly
  why it alone must speak the contract.

The wall-clock channel — GUI stop button, device handle, code — is orthogonal
and untouched: that is the control plane's operator path. The sim-time,
model-detected channel specified here meets it at §12.9 and nowhere else.

### 15.6 Abnormal shutdown: one tail, two entries

Why a `StepError` cannot break §12.9: **the boundary is all-or-nothing outside
the sim task**. Sweeps write into table buffers, integration intermediates
live in workspace, and the only externally visible act is snapshot publication
at the very end of the sequence. A boundary that throws has published nothing;
the last *published* snapshot — a complete, consistent boundary by
construction — is still the newest thing any device, logger or waiter has
seen. The abnormal path is therefore: **discard the failed boundary, promote
the previous snapshot to final, and rejoin the ordinary tail.** The protocol
becomes one tail with two entry points — graceful entry after a *completed*
final boundary, abnormal entry after a *discarded* one — and everything
downstream of "final snapshot" runs identically: sticky stopped, waiters
woken through the boundary-counter + `Condition` path (they observe stopped
rather than a new boundary — no device task hangs), `unblock!`/close hooks,
named joins with timeout. This fills the seat §12.9's "loop failure runs the
same protocol from the catch path" reserved.

Tail hygiene: the hooks are user code too, so each is individually
caught-and-logged — shutdown runs to completion even if a device's hook
misbehaves — and the join timeout already bounds a hook that hangs rather
than throws.

What is lost is quarantined: the state stores may hold mid-boundary values (a
half-written `m`, integration intermediates). They are retained on the
errored `Simulation` for post-mortem inspection, but an errored sim is
terminally stopped, not resumable — the reproduction tool is trace replay,
not resurrection. The published record (snapshot chain, log, trace) ends at
the last consistent boundary; nothing downstream of the sim ever sees half a
boundary.

### 15.7 Tooling consequences: provenance and the component library

Computed exports gain protagonism under this section — termination chains are
their second structural customer after generic-boundary contracts — and two
commitments, a library and an idiom follow:

- **The `faces` helper family grows deliberately.** Predicate-based selection
  (an `endswith`-style filter alongside `except`/`only`) is a natural
  extension: still explicit at the declaration site, still evaluated at
  build, still printable — the blessed side of the auto-bubbling line, where
  the author writes down the *rule* and the build evaluates it into
  inspectable data.
- **The `Build` printer owes face provenance.** For every root face, the
  resolved chain down to the producing terminal (`"crashed" →
  aircraft/monitor/out ← systems/ldg/{left,right,nose}/damaged`). Once faces
  are computed rather than hand-listed, "what does this face actually reach"
  is a question the artifact must answer, not the reader — and the same
  rendering serves the wiring diagnostics, which already carry endpoint
  paths.
- **A standard component library** makes good on §6's junction promise: when
  reduce-ports were rejected, the argument leaned on explicit junctions being
  *cheap*, and a junction hand-written per arity per type is not. Starting
  inventory strictly from demonstrated need — wrench/scalar summing
  junctions, the Bool gates the termination chains use, and `UnitDelay`, the
  spelling §5.4's second loop-breaking remedy needs — growing by
  migration demand only (Simulink's library is a language; this is a
  toolbox). Doctrine: **library blocks are ordinary components** — no
  framework privileges, no special vocabulary — which keeps schema authority
  total and makes the library a permanent ergonomics torture test: if a
  three-input OR gate is painful to write under the declaration rules, the
  rules are wrong. Arity comes from a type parameter (`Or{N}` builds
  `(in1 = Bool, …, inN = Bool)` programmatically — §13.2-blessed derivation,
  and an early validation that the contract functions support parametric
  components). Tier-transparency falls out of settled semantics: a stateless
  continuous `h_xu` recomputes every sweep, so fed ZOH-held discrete signals
  its output changes only at ticks — no tier-neutral kind needed.
  `UnitDelay{V}` is a **discrete** leaf at `K = 1`: `init_z = (v = zero(V),)`,
  `h_z` publishing the stored value and `g` storing the incoming one — the
  tier's native `z⁻¹` (§11.6, the shape §17.2's `sat_out_0` hand-writes), no
  framework support. Its tier semantics are the point and must be stated
  wherever the remedy is recommended: inserting one into a *continuous* loop
  moves that signal onto the discrete tier and inserts a `Δt_base`-scale ZOH
  into the model's mathematics — a modeling decision, the delayed signal being
  genuinely sampled, not a transparent wire, which is why §5.4's diagnostic
  says so rather than offering the remedy as free. A
  migration-phase deliverable.
- **The component test rig** is the library's companion idiom: a one-child
  assembly whose `exports` surface the child's entire input face set —
  `faces(rig, "child")` verbatim (§13.8) — so any component can be built
  and simulated in isolation: every input becomes a root slot fed by
  ordinary conditions and devices, and every output is observable in the
  snapshot table. `design_world`'s little sibling (§16.9): where
  `design_world(ac)` mounts an aircraft in a minimal world for trim and
  linearization, the rig mounts one component behind a root contract for
  unit tests and open-loop probing. Deliberately ordinary machinery end to
  end — no framework support, which makes it, like the library blocks, a
  standing ergonomics test of the declaration rules.

---

## 16. Stopped-sim services

§14.6 previewed the services — initialization, trim, linearization, capture —
as Stratum-C clients. Everything they share reduces to one artifact: the
**condition value**, the datum that says "set this build to this state."
§16.1–§16.4 settle its representation, composition and application;
§16.5–§16.6 the boundary-zero sequence and slot totality; §16.7–§16.9 the
trim service in full; §16.10 linearization and `capture`.

### 16.1 Conditions are path-addressed overlays on the declared defaults

A condition may specify: continuous state fields (`x`), modes (`m`,
continuous components only, §3.2) and discrete state (`z`) — addressed by
§13.6 slash path plus field name — and root
input slots, addressed by face. Never outputs (derived data) and never
workspace (scratch). Entries are validated in the §15.1 batch register: full
list, violations collected, one `BuildError`.

**Overlay base = the declared defaults, always.** Every store has a declared
initial value (declaration-by-initial-value, §13.2), so conditions are
naturally sparse: applying one means "fresh run from the `init_*` defaults,
with these overrides." The alternative base — the stopped sim's current
stores — was rejected: it makes the result depend on run history, exactly the
hidden input the trace-header discipline exists to kill. Warm restart needs no
second semantics: a `capture` service reads the current stores **and root
slots** back *as a condition value* (capture → tweak → apply) — slot coverage
is what makes the captured condition total, hence re-applicable under §16.6 —
the same gather the trace header already needs. One mechanism, two uses.

**The mirror-tree spelling was rejected** (nested NamedTuples shaped like the
assembly): the same information, but a second spelling of structure that must
be zipped against the real one, ragged under partial specification, and
outside the path vocabulary that `connections`, diagnostics and the `Build`'s
provenance tables already speak.

**Doctrine.** This does not reopen §15.5's observation-by-path rejection:
that was *runtime* coupling — a root-authored predicate reaching through
generic seams it does not own, breaking on substitution. A condition is a
*design-time statement about a concrete build*, authored in the same register
as `connections` (which also speaks paths, about children its author owns).
§16.2's composition law makes the parallel exact.

**Pre-sweep doctrine.** Condition writes precede the first sweep by
definition; a would-be init value that depends on swept outputs is either
analytically known to the caller (trim's `α_filt = α_a`: α is a *decision
variable* — the value is known above, not computable below) or an equilibrium
constraint, i.e. a job for the trim service, not for init.

### 16.2 Fragment composition: locality without schema

The rejected alternative deserves its argument, because it names a real need.
`initialize(::C, spec)` as a schema entry — today's `f_init!` reborn
declaratively, with an assembly-level rule routing sub-specs to children —
would keep init knowledge component-local (the engine knows
`n_eng → ω = n_eng·ω_rated`). Rejected on three counts: the routing rule
makes specs a tree that must mirror the assembly tree — the two-artifact
drift that killed the assembly builder (row 39), and call-tree composition
where this design uses data composition everywhere else; spec types with
`@kwdef` defaults are a *second home* for every store's default, competing
with `init_x`; and partial overlays would need a per-field "unspecified"
protocol, while slots (root-face vocabulary) would still need the path layer
beside it — two mechanisms where one suffices. The kill shot on its
motivating example: the aero filter's "`α_filt` must equal `α`" is knowledge
about a *swept output* — a component-local `initialize` computing it would
need its inputs already swept, turning init into a third scheduled sweep.
Today's code doesn't do that either: the value is passed down from where it
is known (§16.1's pre-sweep doctrine).

What preserves the locality is an idiom, not schema: **fragment functions** —
ordinary functions, shipped beside the component, dispatched on it:

```julia
condition(eng::PistonEngine; n_eng) =
    fragment(x = (ω = n_eng * eng.ω_rated,), m = (phase = Phase.running,))
```

composed by *pull* from the structure's owner —

```julia
condition(sys::C172XSystems; n_eng, α_a, β_a) = merge(
    at("pwp/engine", condition(sys.pwp.engine; n_eng)),
    at("aero",       fragment(x = (α_filt = α_a, β_filt = β_a))))
```

— with dispatch selecting variant-specific methods (the c172s/c172x actuation
split costs no upstream edits). The three combinators are constructors of an
**inert, lazy tree** — no path arithmetic at composition:

```julia
struct Fragment{X,M,Z,S}  x::X; m::M; z::Z; slots::S  end  #self-vocabulary payloads; no paths
struct Scoped{N}  prefix::String; node::N  end             #at(prefix, node): stores, never applies
struct Merged{T<:Tuple}  nodes::T  end                     #merge(ns...): collects; order = diagnostics only
```

Every node is isbits except the interned literal prefixes, so **rebuilding
the tree per trim iteration is stack-only construction** — the zero-alloc
property of today's `assign!` loop, preserved structurally. `fragment`'s
payloads speak only about the component at the authoring point; addressing
children is exclusively `at`'s job (one way to say everything). A `slots`
payload names faces *of the authoring level's contract*; resolution walks the
export chain to the root slot and errors if the face never surfaces — an
internally-wired input has no slot and writing it is meaningless (the first
sweep overwrites), and unexported stays unpokeable for init exactly as for
the GUI (§12.5, §17.4).

**The locality law** is §8's, third instance (connections, computed exports,
now conditions): each level speaks its own fields, its declared children's
names, and its own faces; delegation by dispatch at every genericity seam;
deep `at` paths legitimate exactly where deep connections are — within an
owned concrete subtree. Absolute paths exist only in the flattened entry
list, a *compiled derivative* of the composition, as slot offsets are of
`connections`. Substituting a component invalidates precisely the fragments
its owner shipped, nothing else. Enforcement status is also §8's: convention
(ownership is a fact about who maintains the code, not build-visible),
available and idiomatic rather than machine-checked. `fragment`/`at`/`merge`
are §15.7 standard-library material — ordinary artifacts, no privileges.
Merge collisions — two entries on one leaf — are errors at resolution
reporting *both* provenance chains, and the message names the layering
combinator: "`merge` is collision-intolerant by design — use
`override(base, patch)` to layer." Last-writer-wins was rejected in the same
spirit as slot exclusivity (a silent merge is almost certainly a bug) — which
is exactly where this `merge` parts company with `Base.merge`, last-wins on
NamedTuples; the two share a name and not a semantics, and dispatch on the node
types keeps them apart mechanically. The mixed call is closed explicitly: a
`merge(::Fragment, ::NamedTuple)` — or any other blend of a condition node with
a bare NamedTuple — is an **error method**, defined so the call cannot fall
through to `Base.merge`'s last-wins semantics on a payload that looks
plausible, and its message is directive: wrap the NamedTuple in `fragment(...)`
(or `at(prefix, fragment(...))`) and merge nodes with nodes. The
explicit, *ordered* layering spelling — `override` — was deferred here and
admitted one sub-topic later, when slot totality produced its use case
(§16.6).

### 16.3 Resolution: flatten, validate, compile once

Resolution takes the root node plus a `Build`. Flattening is the only place
path strings are ever concatenated — a trivial recursion with a path
accumulator that also records each entry's **tree position** (its
`getfield`/`getindex` step tuple). The batch pass then checks each flat
entry: path resolves (did-you-mean over children), field declared in the
target's `init_x`/`init_m`/`init_z`, value type convertible to the declared
leaf type, slot faces reach root slots, no duplicate
`(path, store, field)`. The `Build` supplies two lookup families: **schema**
(the evaluated declarations: may you write this field, at what leaf type —
the authority) and **layout** (`x` backing ranges, `m`/`z` cell indices, slot
indices from the activation; face chains from Stratum A — the destination).

A valid list compiles to a plan: per leaf, a `Getter{P}` lens (the position
tuple lifted to a type parameter — type-stable navigation of the fixed tree
type), a destination offset, and a **converter baked now** (e.g. an authored
`RQuat` value → the `SVector{4}` state leaf it initializes, through the
type's ordinary `convert`/constructor methods; against a non-nominal
activation's scratch, the `Float64 → Dual` zero-partial embedding —
semantically exact for condition writes: "held at the operating point" *is*
zero partials, §16.10) — a one-time boundary decision that leaves §14.5's
nominal exact-match doctrine for table cells untouched. Converters run here
and in `capture`'s gather (§16.10) — the write paths — never on state views
(§9.1). Overlay
partiality for `m`/`z` cells is baked the same way: the writer holds
`merge(init_m_defaults, overlay)` with the base resolved at compile time
(§16.1's fork).

### 16.4 Two application registers over one plan

Execution is where the paradigm-change tax was feared, and where it
vanishes: all string work, validation and addressing are functions of the
condition's *shape*, and every hot path holds the shape fixed while varying
values. Hence resolve-once/execute-many, with two registers over one plan:

- **Specialized `apply!`** — for services that iterate (trim's per-evaluation
  write, linearization's seeding). Unrolled stores through the baked lenses
  and converters: the same machine operations as today's in-place writes,
  zero-alloc, no strings, no dispatch. The per-iteration shape check is
  §14.5's mechanism transferred: the tree type (which carries the full
  nesting, every field name and leaf type) is proven by dispatch, and a `===`
  sweep over the interned path literals closes the remainder — pointer
  compares that fold to nothing in the all-literal case, with shape drift a
  structured error, not silent corruption. Cost: Julia codegen of ~10–50 ms
  *once per condition shape* — noise against the model's own first-sweep
  warmup (seconds) and the 10³–10⁴ optimizer evaluations it amortizes over.
- **Dynamic walk** — for one-shot init. The same validated entry list,
  executed by runtime dispatch per write: microseconds total, allocation
  permitted (§9.5 — the stopped-sim path was never under the zero-alloc
  regime), and no per-shape codegen: fifty structurally different scripted
  conditions cost fifty walks, not fifty compiles.

Which register a service uses is internal, never user-facing API.

**The read-selector family is closed**: `get_state(path, field[, i])`,
`get_deriv(path, field[, i])`, `get_output(path, field[, i])`,
`get_local(path, field[, i])`, `get_slot(face)`, `get_face(name)` — one
address space for every reader of the model. The names carry a deliberate
`get_` prefix: a selector is a *deferred read* — a value describing the read
the compiled gather will perform — and the prefix both names that action and
keeps six short common nouns out of the namespace user declarations share
with domain code (a bare `local(path, field)` would not even parse — `local`
is a reserved word). `get_local` addresses the component-local cells
(`local_types`, §13.3); `get_face` addresses a root-exported output face —
§12.2's *integration* register, previously recommended but unspellable in
the family.

**A selector resolves against a source, before any client policy applies.**
The table selectors — `get_output`, `get_local`, `get_slot`, `get_face` —
resolve against a boundary snapshot. The store selectors — `get_state`,
`get_deriv` — resolve only against live stores, which only stopped-sim
service evaluations, `capture`, and post-run inspection of the live stores
(§12.2's replay-to-inspect) ever hold: the snapshot deliberately carries no
state stores (§12.2), and `ẋ` buffers are integrator scratch, not
boundary-consistent objects outside a service evaluation. A snapshot-bound
reader naming a store selector is therefore a resolution error at attach
(`ReadBindingUnresolved`), in the didactic register, with §12.2's honest
remedy: declare the field public and read the auto-published port. Client
policy rides on top — row 83's registers restated as a resolver property:

- **Load-bearing services** (trim's `reads`, linearization's surfaces) speak
  the contract: `get_state`/`get_deriv`/`get_output`/`get_slot`, within the
  scopes §8's locality law and §16.2's fragment scoping own; `get_local` is
  rejected at resolution — a service evaluation needing a private
  intermediate is the signal to export it (§16.7).
- **Diagnostic readers** (output-device bindings, GUI panels, log
  inspection) admit the whole family, within the source rule: deep paths,
  `get_local` cells and `get_face` names alike, with the store selectors
  reaching only the diagnostic clients that actually hold stores (`capture`,
  post-run inspection) — a snapshot-bound reader is barred from them by
  source, not by client.
- **`stop_on` is not a family client.** It names root-exported `Bool`
  output faces, period (§15.5, row 60): termination is run policy against
  the root contract, and no path selector reaches it.

**Compiled readers are the gather twin** over this family
and the layout tables: trim's cost read (`ẋ` and output fields), linearization's
Jacobian gather, and `capture`'s full-store readback are one primitive run in
reverse — one machinery, both directions, in the `Build`'s client kit. The
per-iteration ledger for trim: user fragment math (stack-only, the domain
computations unchanged from today) + leaf stores + folded shape check +
sweep — the sweep dominates, exactly as `f_ode!` does today. `apply!` ends at
established stores; making the model *coherent* is boundary zero, §16.5.

### 16.5 Boundary zero: an ordinary boundary with authored incoming transitions

After `apply!` establishes stores at `t₀`, the init service completes the
§11.6 macro-sequence with an empty integrate — project → [sweep → guards →
handlers]\* → due `g` updates → header capture + first snapshot — and that
parity is exact, not approximate. Piece by piece:

- **Project runs.** Authored `x` can sit off-manifold (a hand-assembled
  quaternion ulps off unit norm; a condition writing part of a constrained
  block against fresh defaults). Projection after condition writes is the
  same position it holds after any other `x` mutation, and costs nothing
  when the state is already clean.
- **The sweep runs with every tick due.** `t₀` is a grid point of every
  divisor, so all discrete output stages are gated in, publishing from the
  authored `z` — necessarily, since no earlier tick exists for a ZOH to
  hold. The `t₀` snapshot carries a fully populated table.
- **Events run.** A condition landing in guard-true territory (an authored
  stall flag, a strut authored into contact) fires visibly at `t₀` rather
  than one step later — grounded by the baseline rule (§11.6):
  boundary zero establishes the guard baseline as nothing-holds. Suppression
  would delay the identical firings while hiding the diagnostic that the
  authored condition was not quiescent — §12.5's stage-on-interaction lesson
  (insurance that masks invariant violations is anti-diagnostic). The header
  records the *authored* condition; whatever fires at boundary zero is
  deterministic under replay. (A `stop_on` face already `true` is a
  different category: nothing *fires* — the face simply reads `true` in the
  published `t₀` snapshot and the loop reacts, §15.5.)
- **Due `g` updates run.** This follows from an interval-alignment fact that
  is easy to mis-picture and is hereby a taught contract, sibling to §17.5's
  boundary-sampling line: **a boundary's `g` is the *outgoing* transition** —
  at tick `t_k` it consumes the completed boundary's samples and produces
  `z_{k+1}`, the value the next tick reads; the transition that carried `z`
  *into* `t_k` ran at `t_{k-1}`. Boundary zero is missing its incoming
  transitions on *both* tiers, and both are replaced by authorship: no
  integration over `[t_{-1}, t_0]` produced `x(0)` and no `g` at `t_{-1}`
  produced `z(0)` — the condition did. The outgoing work all runs: `g` at
  `t₀` computes `z(1)` from `t₀` samples, its only opportunity — `z(1)` must
  sit in the store before `t₁`'s gated stages read it. Skipping `g` at
  boundary zero would not preserve the authored `z(0)` — that is already
  preserved, published in the `t₀` snapshot — it would *delete the `t₀`
  sample from the discrete dynamics*: an accumulator
  $z_{k+1} = z_k + \Delta t \, e_k$ authored with $z_0 = 0$ under nonzero $e(t_0)$
  would first integrate $e(t_1)$, the whole sampled-data lattice one period
  late. (The `x`-analogue of `g`-at-`t₀` is not the empty incoming
  integrate but the first *outgoing* one, $t_0 \to t_0 + h$: both authored values
  are the published initial conditions of their outgoing transitions.)
- **`t₀` is an init-service argument** (default `0.0`), never a condition
  entry — time is not a store of any component, and the harmonic grid
  anchors at whatever `t₀` boundary zero runs at. Conditions are time-free;
  `capture` returns condition and time separately for resume-at-time.
- **Trim is untouched by all of this.** Optimizer iterations are raw
  write → sweep → read cycles on the activation — no boundaries, no events,
  no `g`; only the committed solution executes boundary zero. A guard firing
  at commit is a wanted failure signal: today's hand-written trim asserts
  (`!stall`, no weight-on-wheels, `ω > ω_idle`) become the model's own event
  logic, surfaced through the ordinary machinery instead of `@assert`.

### 16.6 Slot totality: the missing-value error and the `override` combinator

Slots are the one store family without declared defaults — §12.3's
bare-types decision, upheld here: a default inside a face declaration would
scatter condition data into the wiring contract and recreate the
competing-defaults problem that §16.2 killed for `initialize` specs. So a
slot's only source before boundary zero is the condition, and three
consequences follow.

**Totality is a precondition of starting, checked by the service.** A
condition value is legitimately partial (fragments compose; trim iterations
write subsets; capture-then-tweak patches leaves) — "every root slot
covered" is not a property of conditions but of *application at boundary
zero*. `init!` (and trim's commit, which runs the same boundary) compares
the resolved plan's slot coverage against the `Build`'s `input_faces`
before writing anything; a shortfall is one batched, declaration-ordered
diagnostic (`UninitializedSlots`, a §15.2 kind) naming every uncovered
face. Pre-write means all-or-nothing: a rejected init leaves the sim
exactly as it was, the same posture as failed trim.

**The probe-value barrier is structural.** §14.3's `probe_value` synthesis
(zero/false/first-enum/`T()`) exists so build-time probes can exercise code
with fabricated inputs; a fabricated zero is a fine probe input and a
terrible flight condition (a silently zeroed `mixture` kills the engine and
sends the user debugging aerodynamics). The services path simply contains
no call to it: a slot gets a condition value or the application errors —
no third branch. Replay likewise never synthesizes: the trace header
records every slot value, and with totality enforced its slot capture is
complete by construction (§12.3's requirement discharged).

**Baselines are aircraft-shipped condition functions, layered by
`override`.** Nobody hand-writes ~20 slot values per script; today
`SystemsInitializer`'s `@kwdef` defaults carry that load, and their
successor is ordinary user math in one authoritative home —
`ready_for_taxi(ac)`, `cold_and_dark(ac)` — returning full-coverage
conditions. But "baseline plus tweaks" collides with §16.2's
duplicate-leaf error by design: the collision *is* the intent. Hence the
fourth node kind, **`override(base, patch)`** — ordered and asymmetric
where `merge` is symmetric and collision-intolerant. At resolution a leaf
present in both takes the patch's value, with provenance recording both
sources ("patch overrode base's `throttle`"); collisions *within* one
layer remain errors; variadic layering
(`override(campaign, aircraft, todays_case)`) composes. Trim uses it on
day one: the committed condition is `override(baseline, solution)` — the
solver's handful of values over full coverage. Rejected alternative: a
service-level base keyword (`init!(sim, patch; base = ...)`) — flatter,
but it hard-codes exactly two layers and moves a composition decision out
of the condition algebra, the one place every other composition decision
lives.

### 16.7 The trim problem: NamedTuple decisions, declared reads, residual vectors

What the aircraft author ships, piece by piece against today's `c172.jl`:

- **Decision variables, initial guess and box bounds are plain, same-shaped,
  all-`Float64` NamedTuples.** The `AbstractTrimState{N}`/`FieldVector`
  supertype dies — its only job was vectorization, which is the service's
  (pack/unpack by field order, shapes checked at setup). Guess, bounds and
  the returned solution share one spelling; `Base.merge(guess, (throttle =
  0.3,))` is free warm-start tweaking; an author who wants a documented
  `@kwdef` struct keeps it privately and converts.
- **`TrimParameters` stays a plain user struct** the framework never sees;
  the assignment is the pure `trim_condition(ac, params, d)` fragment-tree
  function (§16.2), applied per iteration by the compiled plan (§16.4).
- **The read side is declared, then compiled**: `reads(name = get_deriv(path,
  field) | get_output(path, field), ...)` — `get_deriv` addresses a declared
  state field's derivative (validated against `init_x`), `get_output` a
  declared output port (validated against `output_types`); `local_types` are
  not addressable (no `get_local` in `reads`, §16.4) —
  a trim evaluation needing one is a signal the component should export it.
  The compiled reader (§16.4's gather twin) fills a stack-only NamedTuple
  per evaluation.
- **The user supplies a residual *vector*, not a scalar cost.** The
  FlightCore formulation's core — analytic elimination (`θ_constraint`
  substituting the pitch constraint, filter and actuator equilibria imposed
  by construction, the minimal 7-variable search) — is correct and survives
  verbatim as user math. What changes is the numerics: trim is a square
  root-find that FlightCore had to pose as derivative-free scalar
  minimization (BOBYQA over $\|r\|^2$ — a flat quadratic valley near the
  solution, `stopval = 1e-16` as a hand-scaled absolute threshold,
  thousands of evaluations, no per-equation diagnostics) because Jacobians
  through the mutating `f_ode!` chain and the assignment math were out of
  reach. Here the `Dual` activation seeds the decision variables through
  the `T`-generic assignment, sweep and `f` — §14.6's "open option," now
  the *default*: nonlinear least squares on $r(d)$ with exact AD Jacobians
  (trust-region/Levenberg–Marquardt register), quadratic convergence
  (~5–15 evaluations), per-residual physical tolerances, and failure
  reports naming the unbalanced equations with magnitudes. Non-squareness
  degrades gracefully (redundant actuation → weighted/minimum-norm LS;
  infeasible demands → converged nonzero residual identifying the
  impossible balance), and $\partial r / \partial d$ at the solution is free flight-physics
  data (control effectiveness) cross-checking linearization. The
  derivative-free scalar path survives as the fallback: the service squares
  the residuals — today's algorithm as the degenerate case.
- **Recorded, not built**: closed-loop sampled-data trim (append
  $g(z) - z = 0$ residuals via a nondestructive scratch evaluation of `g` —
  structurally impossible under FlightCore's mutating `f_disc!`), and
  on-ground static equilibrium (strut compressions and attitude against
  gear forces) as simply another problem value over the same service.

### 16.8 The trim service: solver seam, scratch stores, commit and report

**The backend seam.** `trim!(sim, problem; baseline, backend =
LevenbergMarquardt())`. The default is an in-house dense
Levenberg–Marquardt: for decision dimensions ~10 with exact Jacobians, the
core (damping loop, small linear solve, convergence test) is ~100 lines —
the §11.2 stepper precedent exactly (tiny needed core vs. heavy dependency),
sharpened by the fact that §16.7's per-residual physical tolerances are a
convergence test no external package spells natively. The backend contract
is documented and value-passed: given an evaluation function (packed
residual vector, optionally with Jacobian), packed guess, bounds and
tolerances, return solution, status and counts. `NLoptBackend(:LN_BOBYQA)`
lives in a package extension, squares the residuals and ignores the
Jacobian — today's algorithm one keyword away — and the framework core
carries zero optimizer dependencies. Box bounds are honored by step
projection, and a decision variable saturated *at the solution* is flagged
in the report ("converged with `elevator` at its upper bound" — the classic
CG-limit diagnostic, today inferable only from mysterious residuals).

**Scratch stores, stated without type luck.** Every `trim!` invocation
instantiates a fresh working store set — `x` backing, `m`/`z` cells, slot
and signal tables, derivative buffer — from the activation's *layout*: the
layout is the reusable compiled artifact, the buffers are per-invocation
and die with the call (stopped-sim allocation, §9.5). The `Dual` backend's
buffers being un-aliasable by type is defense in depth, not the mechanism —
a `Float64` backend (NLopt) gets equally fresh buffers. The invariant is
backend-independent: **the simulation's authoritative stores have exactly
one writer, the commit through boundary zero.** Setup applies
`override(baseline, condition(guess))` to the scratch set once (full
coverage, so sweeps see a complete world); iterations rewrite only the
problem's write-set via the compiled plan; an LM evaluation is one
Dual-seeded sweep yielding `r` (value parts) and `J` (partials) together.
No convergence → no commit → the sim is bit-for-bit untouched, including
"never initialized" — today's warn-but-assign is structurally impossible.

**The report, not an exception.** `trim!` returns a structured
`TrimReport`: converged flag, solution NamedTuple (guess-shaped —
warm-startable), final residuals with tolerances, iteration/evaluation
counts, saturated-bounds list. Non-convergence never throws — it is an
expected *outcome* (envelope-sweep data: hitting the infeasible edge is
information), per §15's exceptions-are-broken-machinery line; a malformed
problem is a `BuildError`-class failure at setup.

**The AD obligation, scoped.** The default formulation requires `Dual`
genericity of exactly: the continuous output-stage chains and `f`, plus the
user's assignment and residual math. The discrete tier's stages and `g`, and the
event system's guards and handlers, never see a `Dual` — frozen constants with
zero partials, semantically exact (§13.2). This is *not a new obligation*: it is the same activation
linearization is defined on, and AD-readiness is a build-checked property
(the Dual probe detonates pinned intermediates with a culprit-naming
`InexactError`; `build(world; activations)` puts it in CI) — robustness by
enforcement, not hope. C172 migration audit (one afternoon, one genuine
item): `Interpolations.jl` tables (propeller coefficient maps, engine
maps) must accept generic scalars — they do, but prefer cubic knots over
linear where partials matter (linear → piecewise-constant Jacobian
entries); in-model saturations (actuator limits, idle/FRC clamps) zero
Jacobian columns when active — LM damping tolerates the rank deficiency
and the report names it, and cruise trim leaves them inactive; the landing
gear is never evaluated off-zero airborne; `norm`-at-zero guards are
already in place (e5efb3a). Fallback per problem: one `backend =` keyword.

### 16.9 Mounting: problems as relocatable values

**What a `TrimProblem` is.** Not a condition, but an **implicitly specified
condition**: `condition` is a condition-*valued function* over the decision
space, `reads`/`residuals` are the equations that pin the free variables
down, `guess`/`bounds` say where to search. Solving makes the implicit
condition explicit, and the commit is then literally an init —
`override(baseline, condition(d*))` through boundary zero. The services
unify as clients of one condition algebra: `init!` applies an explicit
condition, `capture` produces one, `trim!` searches a family for the member
satisfying its equations.

**`at` lifts to problems in five lines.** Every field of a problem is
either condition-producing (path-relative) or path-free — §16.7's rule
that residual math sees only the gathered NamedTuple pays off here:

```julia
at(prefix::String, p::TrimProblem) = TrimProblem(
    p.guess, p.lower, p.upper,                 #path-free: pass through
    d -> at(prefix, p.condition(d)),           #post-compose: wrap each returned tree
    at(prefix, p.reads),                       #reads are inert selector data: same Scoped node
    p.residuals)                               #path-free: pass through
```

Resolution then needs nothing new: the flattening accumulator of §16.3
enters the `Scoped` wrapper and prefixes every entry
(`"vehicle/dynamics"` → `"wing/vehicle/dynamics"`); slot entries authored
in the aircraft's face vocabulary resolve through the export chain *from
the mount point* (`throttle` at `"wing"` → root slot `"wing.throttle"`).
An unexported face fails resolution by name — correctly: an internally
wired input (a scenario component driving the wingman's throttle) is
untrimmable from outside, and the build says so. The service compiles the
scoped condition and reads and runs the identical loop — it never knows
where its paths are mounted.

**The world wrapper dissolves.** Today's `f_init!(::Model{<:SimpleWorld})`
(initialize environment, then call the aircraft's trim) has no successor
method: the environment, the other aircraft and all slots are covered by
the `baseline` condition (§16.6), applied once at setup; the commit is
`override(baseline, at(mount, condition(d*)))`. Method nesting became value
layering.

**"Aircraft as root" is a thin world.** The aircraft is never literally the
root — its environment inputs (§7 function-valued signals) are wired from
provider components — so design tasks use a shipped rig,
`design_world(ac)` = aircraft + `SimpleAtmosphere(wind = NoWind())` +
`HorizontalTerrain`: today's ad-hoc models inside `linearize` promoted to
a named artifact. One register: the "root" case is the shallowest world,
the trim problem mounts at `"aircraft"` like anywhere else.

**Swarm doctrine.** The service solves *one problem at a time*. Sequential
independent trims (trim lead, commit, trim wing against the committed
world) cover weak/one-way coupling; a joint trim is user-side value
composition — concatenate decision NamedTuples under prefixed names, merge
the scoped condition trees, stack the residuals. If joint trims become
routine, a `product(p₁ => "lead", p₂ => "wing")` helper belongs in the
§15.7 library; recorded, not built.

### 16.10 Linearization: surface selectors, one seeded pass, a pure query

**The surface declaration.** Today's per-aircraft `XStateSpace`/
`UStateSpace`/`YStateSpace` structs plus the `get_*_ss`/`assign_*_ss!`
shuttle methods (~150 lines of bookkeeping per variant) become three
selector lists drawn from §16.4's family (`get_state`/`get_slot`/`get_output`), with the
optional component index so a vector leaf yields *named scalars* — the NamedTuple
key is the label control design slices by:
`x = (p = get_state("vehicle/dynamics", :ω_eb_b, 1), θ = get_state("vehicle/kinematics", :θ_nb), …)`,
`u = (throttle_cmd = get_slot("throttle"), …)`,
`y = (EAS = get_output("vehicle/airflow", :EAS), …)`. Validated at resolution
against `init_x`/faces/`output_types` with did-you-mean errors, compiled to
offsets once, relocatable whole via `at(prefix, surface)` — the shuttle
layer's successor is the compiled writer/reader pair, and §9.1's promised
`get_x_ss` deletion is discharged.

**The evaluation.** Instantiate a per-invocation scratch store set
(§16.8's mechanism verbatim), apply the operating-point condition, run
**one** Dual evaluation with a seed direction per `x`- and `u`-surface
entry (chunked internally). Value parts give `ẋ₀`, `y₀`; partials give
`A`, `B`, `C`, `D` simultaneously, exact to machine precision — replacing
four `FiniteDiff` jacobians, their step-size heuristics and ~4n perturbed
evaluations. Unseeded states sit constant at the operating point, and so do
unseeded slots (root-slot cells follow the activation scalar via §13.2's
leaf walk over their concrete declarations; the condition apply embeds
their `Float64` values as zero-partial constants); the discrete tier is
frozen with zero partials — precisely "linearize with `z` held" (§13.2).
Differentiation participation is a per-invocation *seeding* fact, never a
typing fact — one register for `x` and slots alike.

**A pure query, and the shape of `capture`.** Linearization is the first
service with no commit and no boundary zero: scratch buffers only, nothing
becomes authoritative, and today's restore-the-trim dance (re-`assign!`
after `FiniteDiff` dirtied the model) has no successor. Default operating
point = the sim's current committed state via `capture(sim) → (condition,
t)` — the full gather of stores *and root slots* (§16.6's totality makes slot
coverage mandatory for capture → apply): after a
`trim!` commit, `linearize(sim, surface)` is about the trim point with
nothing re-specified; an `about = <condition>` keyword linearizes anywhere
else without touching the sim.

**The returned object and `LinearizedSS`.** `linearize` returns labeled
data — `(ẋ₀, x₀, u₀, y₀, A, B, C, D)` with the surface's label sets — on
which `subsystem`/`delete_vars` survive as pure label-indexed matrix
slicing (no model involvement); the `c172x_ctl` LQR pipeline consumes it
with cosmetic changes. `LinearizedSS` the *component* survives separately
as an ordinary continuous component in the migrated library (`init_x` =
the state vector, labeled faces, the affine update in `h_xu`/`f`) — no
privileges, schema like everyone else.

**Recorded guidance.** Linearization surfaces should select
minimal-coordinate mechanizations — perturbing Euler-angle states is
meaningful where seeding quaternion components steps off the unit-norm
manifold. This is why today's code linearizes the `{NED}` variant;
`design_world(ac)` rigs it, promoting implicit practice to stated rule.
The coordinate choice belongs to the surface author, not the framework.

**Recorded, not built: the sampled-data Dual activation.** The
frozen-exact doctrine is consumer-scoped, not a capability wall: today's
services differentiate the continuous dynamics with `z` held, for which a
frozen discrete output — a ZOH constant with zero partials — is the exact
answer, enforced by the type system (§13.2). Differentiating "through" the
discrete side means differentiating a *different object*: the sampled-data
step map $\Phi : (x_k, z_k, \mathrm{slots}) \to (x_{k+1}, z_{k+1})$ (integrate one period,
then run the due ticks). The extension is additive along existing seams:
Tier-1 parametrization of `z`'s real-scalar leaves (counters/enums stay
pinned, like `m`); opt-in participation on discrete components (frozen-exact
stays the default; a participating component opts in through an explicit
trait, the leaf walk then parametrizing its cells instead of pinning them —
graceful migration, no flag day); one new §14.4 activation ("continuous chain + `f` +
discrete `h_z`/`h_zu` + `g`"); and forward sensitivities through the in-house
RK steppers for free, a payoff of owning the loop (§11.1). The honest
boundary: $\Phi$ is differentiable only where the event pattern is locally
constant — exactness across a firing needs saltation corrections — so the
scope is event-quiescent operating points (which §16.5's guards-at-commit
already makes trim points) plus a loud diagnostic if an event fires inside
a differentiated step. Consumers waiting: §16.7's closed-loop trim door
($g(z) - z = 0$ residuals currently imply the derivative-free fallback,
since frozen `g` has no Jacobian columns) and exact discrete-time
linearization of the full loop (digital design on the exact discretized
plant instead of continuous linearization + Tustin).

**Recorded, not built: declarative non-participation.** The uniform
concrete-declaration rule leaves no schema spelling for "this face is frozen
under differentiation": an opaque wrapper (an FMU, a C aerodynamic table)
strips partials deliberately (`ForwardDiff.value`, §14.5) and shows up in
Jacobians as silent zero rows; a stop-gradient coupling likewise; and no
slot can be marked never-differentiable — an unseeded slot already *is* a
held parameter, so such a marking would be pure protection, not semantics.
If these ever earn schema visibility, the vocabulary is additive, no flag
day: pinned-face declarations validated by the linearization surface
(selecting a declared-frozen output = warning), a feedthrough-graph lint
(a frozen output fed by participating inputs names the severed coupling),
forbid-seeding slot markers rejected at surface resolution. Until a real
consumer shows up, the sharp tool and the visible zero row suffice.

---

## 17. Case studies

### 17.1 `Vehicle` today → this framework

The grounding exercise that validated §5. Current `Vehicle.f_ode!`
(`aircraftbase.jl:142-170`) is a hand-woven instance of the machinery specified here:

| Today (convention) | This design (checked structure) |
|---|---|
| `kinematics.u .= dynamics.x` — velocity extracted directly from the state vector because `f_ode!(dynamics)` can't run yet | `dyn`'s stage-1 output, scheduled first by construction; the artificial loop in `VehicleDynamics` (velocity state-only, accelerations feedthrough) dissolves |
| Hand-ordered `f_ode!` body (kinematics → airdata → systems → route five `dynamics.u` assignments → dynamics last) | Build-time topological sort; wrong wiring = build error naming the cycle or dangling port |
| Velocity state duplicated (`dynamics.x` and `kinematics.u`) with manual sync, incl. `dynamics.x .= kinematics.u  #essential` in `f_init!` | One state, one owner; consumers wire to `dyn.vel` |
| `get_wr_b`/`get_mp_b`/`get_hr_b` generated tree-walk sums | Summing junctions at ownership boundaries, one explicit wire per contributor, exported totals (§6) |
| `f_step!` quaternion renorm + engine-phase/stall-latch checks | `project` hook + Tier-1 events with defined semantics |
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

### 17.2 Torture tests for the §5.2 interfaces: `PistonEngine` and the FCS PID cascade

Two components were transliterated to validate the decoder interfaces before adoption.

**`PistonEngine`** (piston.jl:310-449) — mode enum with three flow regimes, four table
lookups, two embedded continuous PI compensators, boolean transitions, argument-threaded
`fuel_available`:

- The compensator paths (`idle`, `frc`) are pure functions of the engine's own state
  (`ω`), so their complete PI laws — outputs *and* state derivatives — evaluate in
  `h_x`. (The alternative factoring, compensators as child components of an engine
  assembly, also schedules cleanly from the core's stage-1 ports.)
- `h_xu` runs the lookup chain and the mode branch once; `f` is a three-field copy
  (`ω̇`, `ẋ_idle`, `ẋ_frc`). Under the orthodox split, `f` would reproduce essentially
  the whole `f_ode!` body — four lookups and the mode branch — ×4 RK stages per step.
- `f_step!`'s transitions become Tier-1 events with mixed predicate/threshold guards
  (§2.1); `fuel_available` becomes an ordinary port (state-derived at the fuel system,
  hence stage-1 — no loop).
- Forced publications: none — everything `f` reads was already in `PistonEngineY`.

**`PID`** (control.jl:431-471) and the C172X FCS — the discrete side's representative:

- The current update entangles outputs and next state by construction (`y_i = x_i`:
  this tick's integral-path output *is* the updated integrator state). Under §5.2 the
  law runs once in `h_zu`, publishing paths, saturation and the updated states;
  `g` is a three-field copy. Under the orthodox split, `g(z, u, t)` would reproduce
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

**The supervisor slice: scheduled gains and bumpless engage.** One level
above the compensators, today's `c172x_ctl.jl` runs on two idioms the
stores-and-views rules deliberately remove: per-tick gain scheduling by
mutation (`assign!` writing `Ref`-cell parameters from EAS/altitude lookups
every 50 Hz tick, LQR matrix sets included) and mode-transition resets
(`f_init!` plus a bumpless-transfer latch, hand-ordered *before* the same
tick's `f_periodic!`). Both survive as dataflow.

*Scheduled gains are inputs.* A scheduler component owns the lookup tables as
inert parameters, reads the scheduling variables as inputs, and publishes one
gain bundle per compensator; compensators consume gains as `u`. What mutation
hid, ports expose: gain trajectories become observable in log, trace and
replay (the `Ref` writes were invisible to all three), the feedthrough graph
carries the dependency, and linearization holds unseeded gain slots constant
with zero special-casing (§16.10). One-shot design-time gains (`robot2d`'s
controller synthesis at init) are construction-time parameters or stopped-sim
service outputs — not a runtime write path.

*Resets are same-tick inputs, consumed in the output stage.* The supervisor
publishes `engage` and the latch value from its own feedthrough stage; the
compensator, topologically after it, honors them **this tick**:

```julia
h_zu(c::PI, (; z, u)) = (; u_cmd = u.engage ? u.u_latch : c.k_p*u.e + z.x_i)
g(c::PI, (; z, u, Δt)) = (; x_i = u.engage ? u.u_latch - c.k_p*u.e
                                           : z.x_i + c.k_i*Δt*u.e)
```

Honoring the reset only in `g` is legal and means something else: the state
still lands correctly at the next tick, but the *output at the engagement
tick* was already published from the stale state during the sweep — and under
ZOH the plant integrates a full step under it. That one-tick-late command is
exactly the bump that bumpless transfer exists to remove, and no diagnostic
can catch it (both spellings are meaningful designs). The update stage cannot
rescue its own boundary — republish-from-`z⁺` is rejected (row 67) — so the
output stage is the *only* same-tick path. Today's hand-ordering (`f_init!`
before `f_periodic!` in one call) is this contract enforced manually;
Appendix A carries it as the same-tick reset entry, and §12.5's
bumpless-engage answer — engage semantics live in the FCS — presupposes
exactly this spelling. One relative outside the FCS: the landing gear's
level-triggered cross-component reset (`!wow` re-initializing the friction
regulator every step) becomes an edge-triggered event owned by the regulator
— a semantic tightening recorded in §19's migration mapping.

### 17.3 Torture test for the §12 staging shapes: filter, joystick and GUI

The exercise that selected per-device cells (§12.3) and produced the §12.5
contracts. Setup (the `sketch_io.jl` listing): a first-order filter with
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
- **Note**: under slot exclusivity (§12.3) the contested-`u_cmd`
  scenario itself becomes an attach-time error — the drag-against-the-stream
  phase cannot arise. The test's verdicts on cell *shapes* (atomicity,
  coalescing, pause behavior, the peek rule) stand; its conflict-precedence
  findings and the active-widget stage-every-pass contract (§12.5) do not —
  both were patches for the contested-slot world this test examined, and that
  world is gone.

### 17.4 The interactive C172X demo: the periphery under load

The full-fidelity successor to §17.3, against the real deployment:
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
  inputs exported to root faces; GUI as the interactive register's writer via
  §12.5; no machinery. The interactive surface is *not* one thing: pilot commands cluster
  under a prefix; environment knobs stay with their components' panels.
- **The Xv1 actuator sliders**: FlightCore's dead sliders; resolved read-only by
  §12.5. No action.
- **Outbound** (XPlane12: control-surface angles, nose-wheel steering, prop
  speed/phase, pose, `t`): a snapshot-consuming device, pure `map_output` on the
  device task (§12.2). No friction found.
- **Init/trim, pause/pace, post-run plots**: stopped-sim services (§16), control
  plane (§12.6), log/trace (§12.2–§12.3).

**Architectures examined here and rejected** (the §12 periphery decisions were
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

**Surface walkthrough.** The demo line by line:

- `SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_LOWS15))` —
  pure value construction; no `Model` wrapper (its jobs move into the build).
  `HorizontalTerrain`'s elevation is a plain field (parameter), its surface type
  an input port: the parameter/port split FlightCore kept implicit in
  `U()`-vs-field convention is now the declaration itself. The aircraft's
  `exports` block carries the `pilot.*` face group in one place, deep routes
  spanning avionics *and* systems — today's mapping writes flaps/brakes directly
  into `act`, bypassing avionics; that bypass becomes a declared route.
- `Simulation(world; algorithm = RK4(), h = 0.02, n = 1, t_end = 1000)` — `n`
  binds `Δt_base = n·h` (§11.5; default 1: base tick every step). The entire
  build pipeline runs here: kind resolution, path validation, face derivation
  (computed exports expanded, printable), two-producers/unconnected checks,
  topological sort, probe passes, rate compilation, flat layout, slot table.
- `init!(sim, ready_for_taxi(ac); t0 = 0.0)` — stopped-sim services (§16;
  trim is its own service, `trim!(sim, problem; baseline, ...)`, whose commit
  runs the same boundary): they write `(x, m, z)`, **establish every root
  slot's initial value**, and capture the trace header. Slot initialization
  decisively belongs here, not in declarations: the trim service writes slot
  values it *solved for* (throttle, elevator) — not declaration constants.
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

**Frame anatomies.** One frame each:

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
- *Wind slider*: sparse CAS-merge, §17.3's uncontested-`τ` case, live in the
  real cast.
- *Pause/un-pause*: control plane; GUI edits hold in its cell (peek displays),
  joystick cell coalesces bounded; un-pause drain applies both (disjoint slots —
  exclusivity makes the contested question unaskable), pacer re-anchors.
- *Window close*: §12.9 verbatim — complete boundary, final snapshot, sticky
  stopped, wake waits, unblock hooks, named-timeout joins.

Remaining open (feeding §19): the `q_sf` home (thin mapping entry vs.
avionics-internal derivation — aircraft design, not framework design).

### 17.5 The strapdown IMU: integrate-and-dump across the tier boundary

The strongest challenge mounted against the §3 kind split, and its resolution.
The general question first: why two leaf kinds at all — why not one all-in-one
primitive carrying `x`, `m` *and* `z`, with `f`, events *and* `g`, purely
continuous or discrete components falling out by whichever facets an author
declares? (Kind is already read off declaration shape, §13.5 — the question is
whether the two declaration sets should be exclusive.)

**Why the merge buys nothing.** The split is between *time bases*, not state
kinds — the continuous primitive is already hybrid (`m`, guards, handlers, §3.1);
what separates the kinds is sweep-driven versus tick-driven execution. And the
settled rules force a merged component's two halves to communicate exactly as two
siblings do: one home per datum (§5.2), `f` has no `z` view, `g` has no `x` view,
and `z⁺` is decoded only at the owner's next tick (`g` runs last) — the very
fact that makes ticks→events structurally impossible and terminates the boundary iteration (§11.6). Cross-tier
influence inside the merged kind still routes through published table cells, so
the all-in-one component is an assembly of two primitives in a trench coat. Its
costs, meanwhile, are real: a stage cannot run both every sweep *and* only at
ticks, so the merged kind needs four tier-disambiguated stage functions; tier
stops being a component property readable off one `output_types` signature (§13.2) and
becomes per-port vocabulary; every kind-implied obligation (rate required,
`K`-on-continuous error, `Δt` bundle availability, `Dual` activation membership) becomes a
facet-conditional web; and the sampling seam — ZOH and the `z⁻¹` delay, the most
bug-prone boundary in a flight-control stack — disappears into a monolith where
the split keeps it a visible wire. (Simulink and FMI allow the fused block;
sample-time propagation confusion is the documented price.)

**The counterexample** (from a pre-design FlightCore sketch, `navsensors.jl`,
whose operative content and companion derivation note are recorded here in
full): a
strapdown IMU integrates raw increments continuously — `ẋ.ϑ_c = ω_ic_c`,
`ẋ.υ_c = f_c_c`, the coning attitude increment `q_c_cc` and the sculling integral
`ẋ.υ_c_sc = q_c_cc(f_c_c)` — and `f_disc!`, at the IMU's own `Δt`, reads the
integrals, publishes the sample, and **zeroes them**. In interval terms, the
sketch's piecewise quantities are integrals over $[t_{k-1}, t_k]$ with their
weights re-anchored at each reset: `ϑ_c` $= \int_{t_{k-1}}^{t_k} \omega^{c}_{ic} \, dt$
and `υ_c` $= \int_{t_{k-1}}^{t_k} f^{c} \, dt$ from zero,
`q_c_cc` $= q_{c_{k-1} \to c(t)}$ from identity, and
`υ_c_sc` $= \int_{t_{k-1}}^{t_k} R^{c_{k-1}}_{c} f^{c} \, dt$ with the rotation
anchored at the interval
start — exactly the forms the differencing bullets below recover from the
cumulative stores, term by term. The reset is periodic, not
condition-triggered, so events are the wrong tier; and it is a discrete-tier write
into continuous state, exactly the operation this design forbids (`g` writes only
its own `z`; handlers are the sole `x`-resetters, and they are guard-driven).
Integrate-and-dump falls squarely into the crack between the kinds:
tightly-coupled continuous and periodic dynamics in one physical device.

**The idiom: integrate-and-difference.** The reset is eliminable by algebra, not
approximation. Every interval-relative integral becomes a *cumulative* one; the
sampler differences against the previous sample, held in its `z` — the textbook
sampled-data latch, and the only new store (the memory the reset used to erase):

- *Raw increments* (linear): $\Theta(t) = \int_{t_0}^{t} \omega^{c}_{ic} \, dt$,
  $\Upsilon(t) = \int_{t_0}^{t} f^{c} \, dt$, never reset;
  $\vartheta_c = \Theta(t_k) - \Theta(t_{k-1})$,
  $\upsilon_c = \Upsilon(t_k) - \Upsilon(t_{k-1})$.
- *Coning*: cumulative $q(t) = q_{c_0 \to c(t)}$ with
  $\dot{q} = \tfrac{1}{2} \, q \otimes \omega^{c}_{ic}$ from
  identity at $t_0$; the interval rotation is $\Delta q = q(t_{k-1})' \circ q(t_k)$, exact by
  right-invariance ($\Delta q$ satisfies the same ODE with the same body rate).
- *Sculling*:
  $\int_{t_{k-1}}^{t_k} R^{c_{k-1}}_{c} f^{c} \, dt = q(t_{k-1})' \, ( V(t_k) - V(t_{k-1}) )$
  with $\dot{V} = q(t)(f^{c})$. The derivation, in two steps: re-anchor the rotation through the fixed $c_0$ frame,
  $R^{c_{k-1}}_{c} = (R^{c_0}_{c_{k-1}})^{\mathsf{T}} R^{c_0}_{c}$, so the $c_{k-1}$-dependent
  factor — constant over the interval — exits the integral; what remains is the
  cumulative integrand, and splitting its range at $t_{k-1}$ gives the
  difference of the running store:

  $$\int_{t_{k-1}}^{t_k} R^{c_{k-1}}_{c} f^{c} \, dt
  = (R^{c_0}_{c_{k-1}})^{\mathsf{T}} \left( \int_{t_0}^{t_k} R^{c_0}_{c} f^{c} \, dt
  - \int_{t_0}^{t_{k-1}} R^{c_0}_{c} f^{c} \, dt \right)
  = q(t_{k-1})' \, \big( V(t_k) - V(t_{k-1}) \big)$$

  — in code, the sampler line
  `υ_c_sc = z.q'(u.V - z.V)`. The factor leaving the integral is the **anchor
  change between two inertially-fixed frames** — constant because $t_{k-1}$ is
  in the past and latched. The physical intra-interval rotation, the thing
  sculling corrections are *about*, stays inside the integrand via $q(t)$: every
  RHS evaluation, RK stages included, applies the current cumulative attitude,
  exactly as the sketch applies the current `q_c_cc`.

**Exactness condition, stated once**: interval-relative integrals factor into
cumulative ones whenever the interval-dependence enters through a *left action by
the interval-start value of a cumulatively-integrable quantity* — the identity
action for linear integrals, right-invariance for attitude increments, constancy
of the anchor change for sculling. Two provisos: the cumulative attitude must be
integrated with the **inertial** rate, so the anchor frame is inertially fixed and
the pulled factor rigorously constant (anchoring to a rotating reference breaks
the factorization); and the equivalence survives discretization — quaternion
kinematics is linear in `q`, every RK stage composes on the right, and left
multiplication by the constant anchor commutes through, so the formulations agree
to machine precision, not merely in the continuous-time limit. Numerics of never
resetting: `q` stays unit under `project` (better conditioned than the sketch's
`normalization = false` + reset); `Θ`/`Υ`/`V` grow linearly, so differencing
loses relative precision — after an hour of flight, order $10^{-11}\ \mathrm{m/s}$ per sample
against $10^{4}\ \mathrm{m/s}$ totals, six-plus orders below any error model worth simulating.

```julia
struct IMUIntegrals <: AbstractComponent
    t_bc::FrameTransform
end
init_x(::IMUIntegrals) = (Θ = zeros(SVector{3}), q = SVector{4}(1.0, 0, 0, 0),
                          Υ = zeros(SVector{3}), V = zeros(SVector{3}))
input_types(::IMUIntegrals) = (q_eb = RQuat{Float64}, r_eb_e = SVector{3,Float64},
                          ω_eb_b = SVector{3,Float64}, a_ib_b = SVector{3,Float64},
                          α_ib_b = SVector{3,Float64})
output_types(::IMUIntegrals) =
    (Θ = SVector{3,Float64}, q = SVector{4,Float64},            # auto-published state
     Υ = SVector{3,Float64}, V = SVector{3,Float64},
     ω_ic_c = SVector{3,Float64}, f_c_c = SVector{3,Float64})   # instantaneous truth

# h_xu: the sketch's f_ode! math verbatim (lever arm, gravity, Earth rate) → (; ω_ic_c, f_c_c)
function f(imu::IMUIntegrals, (; x, y))
    q = RQuat(x.q, normalization = false)              # §9.1's explicit cast
    (Θ = y.ω_ic_c, q = SVector{4}(Attitude.dt(q, y.ω_ic_c)), Υ = y.f_c_c, V = q(y.f_c_c))
end
project(imu::IMUIntegrals, x) = (; x..., q = normalize(x.q))   # SVector normalize

struct IMUSampler <: AbstractComponent end
init_z(::IMUSampler) = (Θ = zeros(SVector{3}), q = SVector{4}(1.0, 0, 0, 0),
                        Υ = zeros(SVector{3}), V = zeros(SVector{3}))
input_types(::IMUSampler)  = (Θ = SVector{3,Float64}, q = SVector{4,Float64},
                         Υ = SVector{3,Float64}, V = SVector{3,Float64})
output_types(::IMUSampler) = (sample = IMUSample,)   # discrete kind: cells pin (frozen-exact)

function h_zu(s::IMUSampler, (; z, u, Δt))
    q_z = RQuat(z.q, normalization = false);  q_u = RQuat(u.q, normalization = false)
    ϑ_c = u.Θ - z.Θ;  υ_c = u.Υ - z.Υ
    Δq  = q_z' ∘ q_u                                   # interval rotation, exact
    υ_c_sc = q_z'(u.V - z.V)                           # constant anchor change pulled out
    (; sample = IMUSample(; ω̄_ic_c = ϑ_c / Δt, f̄_c_c = υ_c / Δt,
                            ϑ_c, ϑ_c_cc = RVec(Δq)[:], υ_c, υ_c_sc))
end
g(s::IMUSampler, (; u)) = (Θ = u.Θ, q = u.q, Υ = u.Υ, V = u.V)   # the latch
```

The `IMU` assembly wires the four integral ports across, holds the error model as
a discrete sibling consuming `sample`, and leaves the sampler at `K = 1` in its
own scope — the parent sets the IMU's rate (§13.7). `Δt` in the stage bundle is
the §11.5 single source of truth, put there for exactly this kind of discretized
law. (Initialization
consistency — sampler `z` must equal the initial integrals or the `t₀` sample is
wrong — holds by default at zeros/identity, and boundary zero discharges the
rest: its due `g` latches `z ← integrals(t₀)` for every subsequent sample, so
only the `t₀` sample itself depends on the authored `z` — a condition-authoring
obligation under trim, §16.5.)

**Why `u.V` is fresh — the line that would silently zero.** The sculling line is
correct only because a due tick samples the *completed* boundary: if `u.V` still
held the previous boundary's decode it would equal `z.V` exactly (that is the
value `g` latched), and sculling would vanish without an error anywhere. The
guarantee is the §11.6 macro-sequence, not a scheduling accident: integrate →
project → sweep, with the due sampler's stages gated *into* that sweep (§11.5) and
the integrals arriving at stage-1 position (auto-published state, §5.2) — before
any stage-2 function runs, regardless of topological placement. The rest of the
timeline closes consistently: the sampler's `h_zu` decodes `z` (the `t_{k-1}` latch)
*before* `g` runs — the `z⁻¹` semantics — and after event quiescence `g` latches
the `t_k` values for the next tick; same-boundary events re-run the gated stages
in their re-sweeps, so `g` and external readers see the settled boundary.

**Author-knowledge note** (user observation, recorded as a documentation
obligation): the clean implementation leans on the author *knowing* that "sampling
at `t_k`" means post-integration, post-projection, stage-1-fresh state. That
knowledge must be part of the framework's taught contract — §11.5/§11.6 semantics
stated in component-author documentation with this IMU as the worked example — not
internal lore. The failure mode of not knowing it is instructive: an author who
distrusts the sweep order adds a defensive one-tick delay or re-derives the
integrals in the sampler, silently degrading the model.

**When the coupling is genuinely two-way: the latch-back wire.** The IMU's
coupling is one-directional (integrals → sampler). If the flow itself needed the
interval-relative value — integrator saturation within the sampling interval,
say — the latch becomes a wire back: the sampler publishes the sample-instant
values from its *feedthrough* stage (`h_zu` reads `u`, so the latch port carries
the current tick's values, ZOH until the next; an `h_z`-published latch would be
one period stale), and the continuous `f` computes `x − u.latch`. Both cross-wires
consume the other side's ports and the schedule stays acyclic (integrals stage 1 →
sampler `h_zu`; sampler `h_zu` → the integrals' `f`-edge, §5.3). The "reset"
becomes a visible tier-crossing feedback loop — which is what it always was,
physically.

**Verdict.** The strongest counterexample landed on the two-kind taxonomy with
*less* code than the fused original — same thirteen integral scalars, same math,
minus the reset block — and three structural gains. The sampling seam became a
wire. The sketch's incidental violations became visible structure: the
`CircularBuffer` mutated inside the component struct (constants) moves to the
consumer's `z` or falls out of the log; the parent-called `f_disc!(errors)`
becomes a discrete sibling, making the truth/corrupted sample pair separately
loggable. And linearization got sane: under a `Dual` activation the discrete tier
is held (§13.2), and "integrators that never reset" *is* the cumulative
formulation — the framework's rules pushed the model into the only form its own
linearization semantics could coherently handle. Residual escape hatch, recorded
unbuilt: if interval-relative dynamics ever neither factor algebraically nor
tolerate the latch-back wire, the guarded addition is a **tick-triggered handler**
on continuous components (periodic events). Nothing surveyed needs it, and it
would be the camel's nose for the merged kind.

---

## 18. Decision log

The numbered decision rows this document cites as "row N" live in
`framework_decisions.md`: one row per settled decision, each paired with the
alternatives that were weighed against it and the reasons they were rejected.
Row numbers are stable, so a citation here always names the same row there.

---

## 19. Open axes

Still to be settled:

- **Migration.** Outline for FlightPhysics/FlightApps (the Tier-1 parametrization
  pass — whose `Ranged` rewrite targets §13.2's walk rule where `Ranged`
  survives, at ports and parameters: constructor
  discipline admitting the walked scalar with the value parameters left alone,
  plus a `probe_value` method — the `KinData`-style output splits, the
  contributor survey feeding §6's
  aggregation chains — mechanical to extract from today's trait implementations);
  comparison criteria against FlightCore's demonstrated strengths (zero-alloc
  stepping, flexibility, interactive operation); the §15.7 component library's
  starting inventory; the **conventional exported aircraft surface** for
  generic periphery consumers (§12.2's integration register): pose and
  velocity faces with wrapper types — `VelocityData`, field meaning defined at
  the type — as the `KinData` successor's periphery-facing half; the
  **supervisor seam** (§17.2): compensator gain ports plus scheduler
  components (~7 for the C172X), the same-tick reset respelling of every
  mode-transition latch, and the gear's level-triggered reset converted to
  an edge event; the **state-declaration conversion to §9.1's closed
  vocabulary** (each `RQuat` state field becomes its `SVector{4}` backing
  with the explicit `normalization = false` cast at its use sites — today's
  `Attitude.dt` already delivers the 4-wide rate — and each `Ranged` state
  field a plain scalar, its clamp respelled as dynamics or projection, never
  construction); the **exported-name surface**, decided deliberately rather than
  by accident: `condition`, `fragment`, `at`, `capture`
  and the `merge` overload (§16.2) are generic names sharing a namespace with
  FlightPhysics domain code, and `merge` in particular is a piracy surface
  whose mixed-argument methods must stay error methods — the selector
  family's `get_` prefix (§16.4) already settles this for the readers, and
  whether the condition algebra ships behind a submodule is the packaging
  question; and the **§14.7 executor compile-cost re-measurement** on
  the real vehicle skeleton — early, before the executor's shape hardens.
  Residuals: the `q_sf` home (§17.4 — aircraft design,
  belongs here); whether `stop_on` needs a root-declared overridable default
  (§15.5 — reopen only if the ctor argument proves chronically forgotten).

- **GUI panel authoring API.** The semantics are settled (§12.5: derived
  liveness, first-class read-only rendering, own-pending-else-snapshot peek,
  stage-on-interaction, claim transitions); the calling convention — context
  contents, port naming, child composition — is deferred to migration,
  co-designed against the GUI library under §12.5's four constraints.

- **Log and trace persistence.** The in-memory artifacts
  are settled — the log as retained boundary snapshots (§12.2), the
  always-on device-tagged input trace with its header of initial stores and
  slot values (§12.3, §16.5, §16.6), and the primary/derived rule (the log
  is recomputable from the trace, never the reverse) — but nothing on-disk
  is. Deferred to migration, where the consumers exist to ground the
  choices: the HDF5 export scope (whole snapshot log vs. selected
  subtrees), field-handle summarization over retained snapshots (the
  successor to `TimeSeries`'s `getproperty` navigation — today's
  post-processing entry point), and the trace file format, which doubles as
  the reproducibility carrier (§15.4's replay pointers name positions in
  it).

---

## Appendix A. Taught contracts: the author-facing index

The build pipeline enforces structure — declarations, wiring, types,
conformance. A residue of *semantic* facts is unenforceable by any check:
knowing them is what makes component and periphery code come out right, and
not knowing them produces defensive delays, duplicated math or mistimed
samples with no diagnostic firing anywhere (§17.5's author-knowledge note is
the archetype).
This appendix is an **index, not a second home** — one recall line per
contract, with the normative statement staying in the owning section (one
home per datum, applied to the document itself).

For component authors:

- **The stage funnel** (§5.2). Stage name ⊇ bundle ⊇ destructured reads:
  the stage name fixes the maximal legal view set, the component's
  declarations narrow it to the bundle, and the signature's destructuring
  narrows it to actual reads — so "no `u` in the name" *is* the
  no-feedthrough property. The teaching line: stage 1 publishes what you
  know from state alone; stage 2 adds what needs inputs; your dynamics
  read your own published results instead of recomputing them.
- **One home per datum** (§5.2, §4.3). The signal table holds *produced*
  signals only, never transported ones: buffer for `x`, cells for `z`/`m`,
  table for signals — no store mirrors another.
- **Boundary sampling** (§11.5/§11.6; worked example §17.5). "Sampling at
  `t_k`" means post-integration, post-projection, stage-1-fresh state: a
  due tick's gated stages run inside the boundary sweep and sample the
  *completed* boundary. Distrusting this — a defensive one-tick delay, a
  re-derivation inside the sampler — silently degrades the model.
- **Interval alignment** (§16.5). A boundary's `g` is the *outgoing*
  transition: at tick `t_k` it consumes the completed boundary's samples
  and produces `z_{k+1}` — the value the component's *next* tick decodes
  (the sampled-data `z⁻¹` delay, by construction). Hence `g` runs at
  boundary zero: it is the `t₀` sample's only chance.
- **Same-tick reset consumption** (§17.2). A commanded reset is an input.
  For same-tick output semantics the *output stage* consumes it — overriding
  the state-derived path — and `g` stores the matching `z⁺`; a reset honored
  only in `g` reaches the outputs one tick late (the plant integrates a full
  step under the stale command). Both spellings are legal; they mean
  different things.
- **Guard conditions, edges and baselines** (§2.1, §11.6). A guard defines
  a condition — a predicate, or a continuous condition value `σ` with
  positive = holding. Events fire on not-holding → holding *edges* against per-event
  baselines (the previous boundary's quiescent sample): a condition that
  keeps holding fires once, at the boundary where it first held. Boundary
  zero baselines as nothing-holds, so a condition authored already-holding
  fires at `t₀`. The opposite crossing direction is a second event with
  the negated guard.
- **Handler-phase visibility** (§5.2, §11.6). Your `u` is the boundary's
  world, materialized at round start: same-round foreign transitions are
  invisible and arrive through the next round's re-sweep, one round per
  causal link. Your `y` is your own fresh decode, and your `x`/`m` reflect
  earlier events of your own component in declaration order. Coupling
  tighter than one round belongs inside one component.
- **Stop-face sampling** (§15.5). Stop faces are read in completed-boundary
  snapshots; declare a Tier-2 event if the stop needs localizing.

For periphery authors and consumers:

- **Levels, never deltas** (§12.3). Staged input values are levels
  (`press_count = 17`, never `presses += 1`) — idempotent under
  coalescing; button edges ride as monotonic counters.
- **Derived liveness** (§12.5). A widget is live iff its port's feed chain
  terminates in a root slot currently unclaimed by a device; there is no
  per-port marking, and unexported ports are unpokeable.
- **The two observation registers** (§12.2, §15.5). A deep snapshot path is
  the *inspection* register: it sees everything and promises nothing
  across builds. An exported output face is the *integration* register:
  curated, writer-independent meaning — the only shield against silent
  semantic drift. Bind faces in anything meant to outlive the current
  build. The store selectors (`get_state`/`get_deriv`) belong to neither:
  they read live stores, never snapshots (§16.4's source rule).

---

## Appendix B. API synopsis: the entry points

The user-facing surface on one page — same rule as Appendix A: an index, not a
second home, with each signature normative only where its owning section settles
it. The author-side declaration surface first, then the operator surface by
lifecycle:

**Authoring** — what a component or assembly defines (§13.2, §13.5–§13.7):

- Continuous leaf: `init_x`/`init_m` (by value), `workspace(::C, ::Type{T})`
  (by allocation), `input_types`/`output_types`/`local_types` (by type),
  `events` — stages `h_x`, `h_xu`, `f`, guard/handler pairs
  (`Event(guard, handler; localize)`), `project`.
- Discrete leaf: `init_z`, `workspace(::C)`,
  `input_types`/`output_types`/`local_types` — stages `h_z`, `h_zu`, `g`.
- Assembly: `connections` (mandatory — the kind marker), `exports`, `rates`.
- Shipped conditions: `condition(::C; kw)` fragment functions (§16.2).

Bundle contents by function kind (the maximal legal sets, §5.2 — signatures
destructure less at will):

| function | bundle fields |
|---|---|
| `h_x` | `x, m, t [, ws]` |
| `h_xu` | `x, m, u, y_x, t [, ws]` |
| `f` | `x, m, y, u, t [, ws]` |
| `h_z` | `z, t, Δt [, ws]` |
| `h_zu` | `z, u, y_z, t, Δt [, ws]` |
| `g` | `z, y, u, t, Δt [, ws]` |
| guard / handler | `x, m, y, u, t [, ws]` |
| `project` | positional `(comp, x)` — no bundle |

Table footnotes, from the bundle law (§5.2) — the sets above are maximal, and
each field is present only if it exists for the component: `u` iff the function
kind may see inputs **and** the component declares `input_types`; `y` iff the
component produces any table cell (`output_types` ∪ `local_types` ∪
auto-published); `x`/`m`/`z`/`ws` iff declared; `y_x`/`y_z` iff stage-1 ports
exist; `Δt` on the discrete tier only. Returns: a stage returns a NamedTuple of
port values (§4.3); `f` returns the layout image of `X` (§9.1); a **handler
returns `(; x, m)` with each key present iff that store exists and the handler
updates it** (§5.2's return law — no padding, `x` complete, `m` partial).

**Build.**

- `build(world) → Build` — standalone; the inspectable contract artifact:
  wire list, face table with provenance, schedule, root slots (§14.2).
  `build(world; activations = (Float64, ProbeDual))` additionally pins
  activation invariants for CI (`ProbeDual` the exported canonical concrete
  probe scalar, §14.4).
- `resolve(asm, path) → AbstractComponent` — the getfield walk along `/`
  segments, enforcing §8's generic-boundary rule at the primitive (§15.3).
- `input_faces(c)` / `output_faces(c) → Vector{String}` — declaration-ordered
  face names (§15.3).
- `faces(asm, path; prefix, sep, except, only)` — the declaration-site helper
  for computed exports (§13.8).

**Deployment.**

- `Simulation(world; algorithm = RK4(), h, n = 1, t_end = Inf,
  stop_on = (), trace = true, log = true, log_every = 1, debug = false)` —
  wraps the build; `h` is required (a domain rate is not a
  framework default) and `RK4` is the default stepper (§11.2); `n` binds
  `Δt_base = n·h` (§11.5); `t_end = Inf` is the honest interactive default, and
  a run with no finite `t_end`, no `stop_on` faces and `pace = Inf` warns at
  start (an unbounded batch run is almost always an oversight); `stop_on` names
  root-exported `Bool` output faces, OR-combined, recorded in
  run metadata (§15.5; walkthrough §17.4). `t_end` and `stop_on` are
  **defaults**, overridable per run at `run!` (§15.5, §12.11). Recording:
  `trace` is the input trace's plain kill switch and `log` the snapshot log's,
  with `log_every` the log's keep-every-kth decimation — admissible on the
  derived artifact only, never on the trace (§12.2, §12.3, row 29). `debug`
  gates the workspace poison (§9.3).
- `attach!(sim, device, binding)` — write side speaks face names through a
  declarative binding table (axis/button → face + conditioning params; the
  binding's enumerated face set *is* the claim — what the device may write,
  not what it will — registered at attach with exclusivity enforced, §12.3);
  read side binds §16.4 selectors,
  validated at attach (§12.3, §12.2, §17.4). Legal at any time, running or
  stopped (§12.3).
- `detach!(sim, device)` — releases the device's claims and republishes the
  roster; the `should_close` exit path runs the same release (§12.3).
- The device handle — one kind, capabilities not taxonomy: read latest
  snapshot, wait-for-next-boundary (§12.8), stage, control access (§12.4).

**Condition algebra** (§16.1–§16.6).

- `fragment(; x, m, z, slots)` — self-vocabulary payloads at the authoring
  level; `slots` names faces of that level's contract.
- `at(prefix, node)` — scoping; stores, never applies. Also lifts whole
  `TrimProblem`s and linearization surfaces (§16.9, §16.10).
- `merge(nodes...)` — symmetric collection; duplicate leaves error with dual
  provenance; mixing a node with a bare NamedTuple is an error method, not
  `Base.merge` (§16.2).
- `override(base, patches...)` — ordered layering; patch wins, provenance
  keeps both (§16.6).
- `condition(comp; kw)` — the shipped fragment-function idiom; aircraft
  baselines (`ready_for_taxi(ac)`, `cold_and_dark(ac)`) are its
  full-coverage instances.

**Stopped-sim services** (§16).

- `init!(sim, condition; t0 = 0.0)` — slot totality checked pre-write
  (§16.6), then boundary zero: project → sweep → events → due `g` updates →
  header + first snapshot (§16.5).
- `trim!(sim, problem; baseline, backend) → TrimReport` — nonlinear least
  squares on the residual vector with exact Dual Jacobians; commit =
  `init!` with `override(baseline, solution)`; non-convergence reports, never
  throws (§16.7, §16.8).
- `capture(sim) → (condition, t)` — full-store gather including root slots;
  warm restart = capture → tweak → apply (§16.1, §16.10).
- `linearize(sim, surface) → labeled (ẋ₀, x₀, u₀, y₀, A, B, C, D)` — pure query, one
  seeded Dual pass on scratch; operating point defaults to `capture(sim)`;
  surface = `get_state`/`get_slot`/`get_output` selector lists with control-design
  labels (§16.10).

**Running.**

- `run!(sim; gui = false, pace = 1, margin = <default>, t_end = <ctor value>,
  stop_on = <ctor value>)` — paced and unpaced runs bit-identical
  (§11.7); the GUI an ordinary rostered device rendered on the calling task
  (§12.4, §12.5), and `gui = true` is sugar for attaching the standard GUI
  device — it moves the loop to a spawned task for the run (§12.1); sugar
  never activates by default. `run!` blocks until the run ends; deviceless
  it is fully synchronous on the calling task; `init!` required first
  (§12.11). `margin` is the single pacing
  knob (§11.7), `0` / a few ms / `∞` spanning the design space. `t_end` and
  `stop_on` override the constructor's defaults **for this run only**, with
  `stop_on` validated against the `Build` here exactly as at construction, and
  the effective pair recorded in the run metadata (§15.5).
- `step!(sim; frames = 1) → frames_advanced` — synchronous partial advance
  through the ordinary frame sequence, bit-identical to the same frames under
  `run!`; `t_plus = <duration>` is the mutually-exclusive duration spelling
  (whole frames until the boundary time covers it); returns the frames
  *actually* advanced, fewer than requested when
  `t_end` or a `stop_on` face ended the run inside the call. Between calls the
  simulation reports `initialized`; `run!` may follow and continues from the
  current boundary; a stepping session is deviceless — write via `stage!`,
  read via `latest` (§12.11).
- `stage!(sim, "face" => value, ...)` — task-free staging from the calling
  task into §12.3's interactive register (surface = the currently-unclaimed
  faces): traced, drained last at the next frame top, surface-checked exactly
  as the GUI's writes (§12.11's harness register; legal under `run!` and
  `step!` alike).
- `latest(sim) → snapshot` — the current published snapshot, the same
  immutable value device handles read (§12.2); the assertion/inspection
  accessor of the harness and REPL registers (§12.11).
- Control plane — pause/un-pause, pace and `margin` changes, stop on a
  separate atomic surface, never staged (§12.6; pacing sits outside the
  semantics, so both are safe to change live).
- Termination — model state via `stop_on` faces read at every published
  boundary (§15.5); shutdown completes a boundary, publishes the final
  snapshot, then joins (§12.9).
- Post-run — the log is retained snapshots; the always-on input trace
  re-drives a fresh `Simulation(world)` bit-identically and is the
  state-trajectory inspector (§12.2, §12.3; on-disk persistence deferred,
  §19).

---

## Appendix C. The diagnostic kind set

The closed set row 58 commits to, made normative: **tests match on kind plus
payload fields, never on message text** (§15.2), so this table — not any
message — is the acceptance-test contract, and adding a kind is a decision-log
entry. Every row's payload is *in addition to* what §15.2 requires of all
diagnostics: paths and names as strings, never instances; the list-in-hand
wherever a did-you-mean renders; the didactic register (state the fix). Owning
sections stay the normative home of each rule; this appendix is an index of the
values, in the manner of Appendices A and B.

Severities, in the vocabulary §15 fixes:

- **build (batch)** — a diagnostic from a declarative pass, collected with its
  siblings and thrown as one `BuildError` at the stratum barrier (§15.1);
- **build (fail-fast)** — raised while *user code* runs (an `exports` body, a
  probe); the first one aborts the phase (§15.1);
- **service** — raised by a stopped-sim service, or by
  `attach!`/`Simulation`/`run!` validating against the `Build`; batched into one
  carrier wherever the owning section says so (§16.1's register, §16.6's
  pre-write check), a single throw at the call otherwise;
- **runtime** — fail-fast during a boundary, reaching §15.4's single catch site
  as a species of `StepError`;
- **warning (runtime)** — the per-occurrence runtime stream of §15.2. The
  *build* warning severity exists and its set is currently empty (row 84).

**Declaration and wiring** (Stratum A):

| kind | payload | owner | severity |
|---|---|---|---|
| `UnknownPort` | destination path, unknown port name, the destination's port list (did-you-mean) | §8, §13.4 w1 | build (batch) |
| `UnconnectedInput` | leaf path, input name, declared entry type, the obligation chain's last level | §8, §13.4 w2 | build (batch) |
| `TwoProducers` | destination terminal, both producer terminals with provenance (sibling wire / ancestor deep route / export entry) | §8, §13.8 | build (batch) |
| `WireTypeMismatch` | both endpoint paths, both face names, declared entry type, producer face type | §8, §13.2, §13.4 w4 | build (batch) |
| `PathResolution` | path, offending segment, sibling field list; for a traversal past a generically-held field, that field's declared type | §8, §15.3 | build (batch) |
| `AbstractAtRoot` | face name, consuming leaf path, the abstract entry | §13.2 | build (batch) |
| `RootSlotTypeConflict` | face name, the consuming paths, their conflicting concrete declarations | §13.2 | build (batch) |
| `IllegalStateLeaf` | component path, `init_x` field name, leaf type, the closed vocabulary (scalar / `SArray` at the common eltype) | §9.1, §13.2 | build (batch) |
| `StoreWithoutUpdate` | component path, store (`x`/`z`), the missing function (`f`/`g`) | §13.2 | build (batch) |
| `EventHalfMissing` | component path, event name, which half, the function that has no method | §13.2 | build (batch) |
| `PrimitiveAtRoot` | root path, component type | §13.2 | build (batch) |
| `KindUnreadable` | component path, type, declarations found, both family lists; did-you-mean when the type holds component-typed fields | §13.5 | build (batch) |
| `KindMixed` | component path, the `connections` declaration and the offending leaf declarations | §13.5 | build (batch) |
| `ContainerMixed` | container field path, offending element keys/indices, their types | §13.5 | build (batch) |
| `StageOnWrongTier` | component path, stage name, the tier its declaration shape announces | §5.2, §13.5 | build (batch) |
| `FaceNameIllegal` | assembly path, face name, the violated invariant (contains `/`) | §13.6 | build (batch) |
| `FaceNameCollision` | assembly path, face name, both entries' provenance (hand-written / computed) | §13.6 | build (batch) |
| `FaceDirectionConflict` | assembly path, face name, its internal endpoints, derived directions | §13.6 | build (batch) |
| `UnknownFaceSelection` | child path, the unknown `except`/`only` names, the child's face list | §13.8 | build (batch) |
| `RatesViolation` | assembly path, offending key, reason (deep key / unknown child / `K` on a continuous child) | §11.5, §13.7 | build (batch) |
| `MissingProbeValue` | face name, type | §14.3 | build (batch) |

**Schedule and contract conformance** (Strata B and C):

| kind | payload | owner | severity |
|---|---|---|---|
| `AlgebraicCycle` | the SCC's member terminals in slash form, the wires among them, optional classification (`real`/`artificial`) with the member whose hop died | §5.4, §10 | build (batch) |
| `ProducedByTwoStages` | component path, port name, both stage names | §4.3, §13.3 | build (batch) |
| `ContractNameCollision` | component path, the name declared in both `output_types` and `local_types` | §13.3 | build (batch) |
| `DeclaredNotProduced` | component path, declared name, the stage-product list and the state-field list | §13.3 | build (batch) |
| `LocalConnected` | wire endpoints, the local cell named | §13.3 | build (batch) |
| `UndeclaredReturnField` | component path, stage, returned field name, candidates (`output_types` ∪ `local_types`) | §13.3, §13.4 w5 | build (fail-fast) |
| `ConformanceFailure` | component path, stage, field-level diff (missing / unexpected / per-field expected-vs-observed), simulation time | §14.5 | build (fail-fast) at probe; **runtime** as a `StepError` species |
| `GuardForm` | component path, event name, observed probe return type, both admissible forms | §14.5 | build (fail-fast) |
| `Tier2GuardPredicate` | component path, event name — a localizing event whose guard probes `Bool` | §14.5 | build (fail-fast) |
| `BundleFieldError` | component path, function kind, requested field, the legal field set, classification (undeclared store / wrong-tier fact / illegal for this function kind) | §5.2, §15.2 | build (fail-fast) at probe; **runtime** thereafter |
| `HandlerReturnKey` | component path, event name, offending key, the legal set `{x, m}` narrowed to the stores that exist | §5.2, §14.5 | build (fail-fast) |
| `UserCodeFraming` | component path, which function, the probe context including synthesized inputs; the original exception as `cause` | §15.2 | build (fail-fast) |

**Deployment, periphery and services:**

| kind | payload | owner | severity |
|---|---|---|---|
| `MissingInit` | the simulation's status, the entry point called (`run!`/`step!`) | §12.11 | service |
| `StopFaceInvalid` | face name, reason (unknown / not root-exported / not `Bool`), the root output-face list; the binding site (constructor or `run!`) | §15.5 | service |
| `AttachUnknownFace` | device id, binding entry, face name, the root input-face list | §12.3 | service |
| `ClaimConflict` | face name, claiming device id, incumbent device id | §12.3 | service |
| `ReadBindingUnresolved` | device id, the selector, path and field, candidates; a `reason` distinguishing an unresolved path from a store selector in a snapshot binding (§16.4's source rule) | §12.2, §16.4 | service |
| `ConditionResolution` | entry path, store, field, offending value type and declared leaf type, provenance chain; sub-kinds: unknown path, undeclared field, unconvertible value, unexported slot face | §16.2, §16.3 | service (batch) |
| `DuplicateConditionLeaf` | the leaf `(path, store, field)`, both provenance chains, the `override` advice | §16.2 | service (batch) |
| `UninitializedSlots` | every uncovered root face, in declaration order | §16.6 | service (batch), pre-write |
| `SurfaceResolution` | surface (`x`/`u`/`y`), selector kind, path, field, optional index, candidates | §16.10 | service (batch) |

**Runtime:**

| kind | payload | owner | severity |
|---|---|---|---|
| `StepError` | the carrier: cursor frame (component path, function, boundary phase — RK stage, event round, localization probe, tick), boundary time, trace boundary index (replay pointer), original exception as `cause` | §15.4 | runtime |
| `NonfiniteState` | component path, the offending state block, boundary time and index | §15.4 | runtime |
| `ChatteringBudget` | component path, event name, boundary time, the exhausted localization budget | §11.4 | warning (runtime) |
| `EventDeferred` | component path, event name, boundary time — re-enabled within the boundary, deferred by the once-per-event rule | §11.6 | warning (runtime) |
| `DebtReanchor` | forgiven debt, the new schedule anchor, boundary time | §11.7 | warning (runtime) |
| `StaleInteractiveEntry` | face name, the incumbent (claiming) device id, the discarded value, which interactive writer (GUI / harness) | §12.3, §12.5 | warning (runtime) |
| `OutOfClaimEntry` | device id, face name, the discarded value, the device's claim set; the incumbent's device id when the face is claimed elsewhere | §12.3 | warning (runtime) |
| `ThreadBudget` | thread count, device-task count, the site (`run!` or an `attach!`) | §12.7 | warning (runtime) |
| `DeviceJoinTimeout` | device id, the join timeout, boundary time and index at shutdown | §12.9 | warning (runtime) |
| `DeviceCrash` | device id, the original exception as `cause`, whether `should_abort` was set | §12.9, §15.4 | warning (runtime) |
| `PoisonSkip` | component path, the skipped workspace stores and their element types | §9.3 | warning (runtime), once per activation |
| `UnboundedRun` | the effective `t_end`, `stop_on` set and `pace` | Appendix B, §15.5 | warning (runtime), at run start |
