# A Modeling & Simulation Framework for Flight.jl — Specification

---

## Contents

- [Part I — Foundations](#part-i--foundations)
  - [1. Purpose and method](#1-purpose-and-method)
  - [2. Formalism](#2-formalism)
    - [2.1 Events: two detection policies](#21-events-two-detection-policies)
    - [2.2 Exclusions (deliberate)](#22-exclusions-deliberate)
  - [3. Component taxonomy](#3-component-taxonomy)
    - [3.1 Continuous component (the hybrid primitive)](#31-continuous-component-the-hybrid-primitive)
    - [3.2 Periodic discrete component](#32-periodic-discrete-component)
    - [3.3 Assembly](#33-assembly)
  - [4. Ports and signals](#4-ports-and-signals)
    - [4.1 Immutable value semantics](#41-immutable-value-semantics)
    - [4.2 Consumers see ports, not stages](#42-consumers-see-ports-not-stages)
    - [4.3 Table mechanics and port granularity](#43-table-mechanics-and-port-granularity)
    - [4.4 Function-valued signals: environment access](#44-function-valued-signals-environment-access)
  - [5. Evaluation order and feedthrough](#5-evaluation-order-and-feedthrough)
    - [5.1 The scheduling problem](#51-the-scheduling-problem)
    - [5.2 Two-stage outputs: signatures, bundles and the hand-off laws](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)
    - [5.3 Structural feedthrough: stage roles, schedule and step boundaries](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)
    - [5.4 Artificial loops and the escape hatch](#54-artificial-loops-and-the-escape-hatch)
    - [5.5 Algebraic loop policy: reject at build time](#55-algebraic-loop-policy-reject-at-build-time)
    - [5.6 Diagnostics: feedthrough tracing](#56-diagnostics-feedthrough-tracing)
  - [6. Composition: connections, aggregation and hierarchy](#6-composition-connections-aggregation-and-hierarchy)
    - [6.1 Connections and hierarchy](#61-connections-and-hierarchy)
    - [6.2 Aggregation: explicit summing junctions](#62-aggregation-explicit-summing-junctions)
  - [7. State and data representation](#7-state-and-data-representation)
    - [7.1 Continuous state: structured immutable, flat backing](#71-continuous-state-structured-immutable-flat-backing)
    - [7.2 Numeric genericity (eltype)](#72-numeric-genericity-eltype)
    - [7.3 Discrete state, modes, and workspace](#73-discrete-state-modes-and-workspace)
    - [7.4 The fused-evaluation lineage (prior art and how we got here)](#74-the-fused-evaluation-lineage-prior-art-and-how-we-got-here)
    - [7.5 Allocation policy: a scoped invariant](#75-allocation-policy-a-scoped-invariant)
- [Part II — Execution](#part-ii--execution)
  - [8. Time and execution](#8-time-and-execution)
    - [8.1 Loop ownership: the framework owns the simulation loop](#81-loop-ownership-the-framework-owns-the-simulation-loop)
    - [8.2 The stepper seam](#82-the-stepper-seam)
    - [8.3 Signal-table consistency is a boundary property](#83-signal-table-consistency-is-a-boundary-property)
    - [8.4 Localization mechanics](#84-localization-mechanics)
    - [8.5 Multi-rate tick scheduling](#85-multi-rate-tick-scheduling)
    - [8.6 Event iteration at boundaries: to quiescence, once per event](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)
    - [8.7 Real-time pacing](#87-real-time-pacing)
  - [9. Runtime periphery: the data plane](#9-runtime-periphery-the-data-plane)
    - [9.1 No shared mutable model: staged writes, snapshot reads](#91-no-shared-mutable-model-staged-writes-snapshot-reads)
    - [9.2 Outbound: snapshot publication](#92-outbound-snapshot-publication)
    - [9.3 Inbound: root input slots, claims and the frozen roster](#93-inbound-root-input-slots-claims-and-the-frozen-roster)
    - [9.4 Inbound: per-device staging, representation and the drain](#94-inbound-per-device-staging-representation-and-the-drain)
    - [9.5 Inbound: the input trace](#95-inbound-the-input-trace)
    - [9.6 Devices: one authoring contract, no taxonomy](#96-devices-one-authoring-contract-no-taxonomy)
    - [9.7 The GUI write path: port resolution, peek, staging contract](#97-the-gui-write-path-port-resolution-peek-staging-contract)
    - [9.8 Diagnostics and liveness: the per-writer cell](#98-diagnostics-and-liveness-the-per-writer-cell)
  - [10. Runtime periphery: lifecycle and orchestration](#10-runtime-periphery-lifecycle-and-orchestration)
    - [10.1 Control plane](#101-control-plane)
    - [10.2 Loop scheduling: wait primitive, yields, thread budget](#102-loop-scheduling-wait-primitive-yields-thread-budget)
    - [10.3 The next-snapshot wait](#103-the-next-snapshot-wait)
    - [10.4 Shutdown protocol](#104-shutdown-protocol)
    - [10.5 Scripts and the mid-run mutation doctrine](#105-scripts-and-the-mid-run-mutation-doctrine)
    - [10.6 Run lifecycle and partial advance](#106-run-lifecycle-and-partial-advance)
    - [10.7 Replay: the trace re-drives the ordinary loop](#107-replay-the-trace-re-drives-the-ordinary-loop)
- [Part III — Authoring and build](#part-iii--authoring-and-build)
  - [11. The declaration layer: components and assemblies](#11-the-declaration-layer-components-and-assemblies)
    - [11.1 Position: a declarative trait layer — plain Julia, no macros](#111-position-a-declarative-trait-layer--plain-julia-no-macros)
    - [11.2 The declaration inventory](#112-the-declaration-inventory)
    - [11.3 Visibility: the contract is the interface](#113-visibility-the-contract-is-the-interface)
    - [11.4 Failure walkthroughs (the error-locality grounding)](#114-failure-walkthroughs-the-error-locality-grounding)
    - [11.5 Assembly declaration: type-based, class by declaration shape](#115-assembly-declaration-type-based-class-by-declaration-shape)
    - [11.6 Paths, wiring and faces](#116-paths-wiring-and-faces)
    - [11.7 Rate scopes](#117-rate-scopes)
    - [11.8 Computed connections and generic boundaries](#118-computed-connections-and-generic-boundaries)
  - [12. The build pipeline](#12-the-build-pipeline)
    - [12.1 Three strata](#121-three-strata)
    - [12.2 The `Build` artifact](#122-the-build-artifact)
    - [12.3 Probing and input synthesis](#123-probing-and-input-synthesis)
    - [12.4 Activations: executable sets, laziness, caching](#124-activations-executable-sets-laziness-caching)
    - [12.5 The always-on conformance check](#125-the-always-on-conformance-check)
    - [12.6 Stopped-sim services as Stratum-C clients](#126-stopped-sim-services-as-stratum-c-clients)
    - [12.7 The compiled executor](#127-the-compiled-executor)
- [Part IV — Failure and services](#part-iv--failure-and-services)
  - [13. Error discipline](#13-error-discipline)
    - [13.1 Reporting policy: collect the checks, fail the evaluations fast](#131-reporting-policy-collect-the-checks-fail-the-evaluations-fast)
    - [13.2 Diagnostics: structured values, one carrier exception](#132-diagnostics-structured-values-one-carrier-exception)
    - [13.3 Build primitives: `resolve` and the face-list accessors](#133-build-primitives-resolve-and-the-face-list-accessors)
    - [13.4 Runtime failures: one catch site, an execution cursor](#134-runtime-failures-one-catch-site-an-execution-cursor)
    - [13.5 Termination is a state, not an exception](#135-termination-is-a-state-not-an-exception)
    - [13.6 Abnormal shutdown: one tail, two entries](#136-abnormal-shutdown-one-tail-two-entries)
    - [13.7 Tooling consequences: provenance and the component library](#137-tooling-consequences-provenance-and-the-component-library)
  - [14. Stopped-sim services](#14-stopped-sim-services)
    - [14.1 Conditions are path-addressed overlays on the declared defaults](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)
    - [14.2 Fragment composition: locality without schema](#142-fragment-composition-locality-without-schema)
    - [14.3 Resolution: flatten, validate, compile once](#143-resolution-flatten-validate-compile-once)
    - [14.4 Two application registers over one plan](#144-two-application-registers-over-one-plan)
    - [14.5 Boundary zero: an ordinary boundary with authored incoming transitions](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)
    - [14.6 Slot totality: the missing-value error and the `override` combinator](#146-slot-totality-the-missing-value-error-and-the-override-combinator)
    - [14.7 The trim problem: NamedTuple decisions, declared reads, named residuals](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)
    - [14.8 The trim service: solver seam, scratch stores, commit and report](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)
    - [14.9 Mounting: problems as relocatable values](#149-mounting-problems-as-relocatable-values)
    - [14.10 Linearization: tap selectors, one seeded pass, a pure query](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)
- [Part V — Grounding](#part-v--grounding)
  - [15. Case studies](#15-case-studies)
    - [15.1 `Vehicle` today → this framework](#151-vehicle-today--this-framework)
    - [15.2 Torture tests for the §5.2 interfaces: `PistonEngine` and the FCS PID cascade](#152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade)
    - [15.3 Torture test for the §9 staging shapes: filter, joystick and GUI](#153-torture-test-for-the-9-staging-shapes-filter-joystick-and-gui)
    - [15.4 The interactive C172X demo: the periphery under load](#154-the-interactive-c172x-demo-the-periphery-under-load)
    - [15.5 The strapdown IMU: integrate-and-dump across the tier boundary](#155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary)
  - [16. Open axes](#16-open-axes)
- [Appendix A. Taught contracts: the author-facing index](#appendix-a-taught-contracts-the-author-facing-index)
- [Appendix B. API synopsis: the entry points](#appendix-b-api-synopsis-the-entry-points)
- [Appendix C. The diagnostic kind set](#appendix-c-the-diagnostic-kind-set)
- [Appendix D. Glossary](#appendix-d-glossary)
  - [D.1 Component model and declaration layer](#d1-component-model-and-declaration-layer)
  - [D.2 Signals and data homes](#d2-signals-and-data-homes)
  - [D.3 Evaluation and scheduling](#d3-evaluation-and-scheduling)
  - [D.4 Time and events](#d4-time-and-events)
  - [D.5 Build pipeline](#d5-build-pipeline)
  - [D.6 Runtime periphery](#d6-runtime-periphery)
  - [D.7 Recording and replay](#d7-recording-and-replay)
  - [D.8 Stopped-sim services and the condition algebra](#d8-stopped-sim-services-and-the-condition-algebra)
  - [D.9 Error discipline and diagnostics](#d9-error-discipline-and-diagnostics)
  - [D.10 Meta-vocabulary](#d10-meta-vocabulary)

---

# Part I — Foundations

## 1. Purpose and method

This document specifies a modeling and simulation framework intended to
replace `FlightCore` as the substrate for `FlightPhysics` and `FlightApps`.
It is the normative statement of the design: what the framework *is*, in
present tense. It is derived from `framework_design.md`, the historical
original, which carries the full record of how the design was reached and
survives in git history. The new framework must match or surpass `FlightCore` in
functionality, performance and flexibility, while being more rigorous and
explicit — reducing the learning curve and the number of latent footguns for
model authors.

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

All design axes are settled — the formalism, the component taxonomy, the
signal and scheduling model, time and execution, the runtime periphery, the
declaration layer, the build pipeline, error discipline and the stopped-sim
services. Only [§16](#16-open-axes)'s items — the migration outline, the GUI panel authoring
API and the log/trace persistence deferral — remain open.

Decision rationale, including the alternatives considered and the reasons
they were rejected, lives in `framework_decisions.md`, cited throughout as
"row N": one row per settled decision. Row numbers are stable, so a citation
here always names the same row there.

---

## 2. Formalism

The framework simulates **hybrid causal systems**, composed of:

- **Continuous dynamics**: $\dot{x} = f(x, m, u, t)$ with algebraic outputs.
- **Multi-rate periodic discrete dynamics**: $x^{+} = g(x, u, t)$ at declared rates, with
  outputs held zero-order between ticks.
- **Zero-crossing events**: guard functions with handlers, under two detection policies (below).
- **Post-step manifold projection**: an optional per-component hook `x ← project(x)`
  applied after each accepted step (quaternion renormalization, DCM orthonormalization,
  any manifold-valued state). This is the cheap end of the projection-methods family
  from geometric integration.
- **External inputs**: injected asynchronously by the runtime (pilot controls, network),
  under the staging rules settled in [§9](#9-runtime-periphery-the-data-plane).

### 2.1 Events: two detection policies

Both policies share one declaration (guard function + handler); only detection
differs:

- **Boundary-detected (default, cheap):** guards are checked for not-holding → holding edges
  against their priors ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)) at step boundaries
  only. No root-finding, no step rejection; the handler fires at the end of the step in
  which the edge was observed. Cost: one guard evaluation per event per step. Fully
  compatible with fixed-step real-time execution.
- **Localized (opt-in, per event):** localization of the crossing instant by root-finding,
  for events where timing precision genuinely matters (mechanics in [§8.4](#84-localization-mechanics)). Runs
  identically in every execution mode: under real-time pacing the
  localization cost is absorbed as pacer debt like any other expensive frame
  ([§8.7](#87-real-time-pacing)) — detection policy never depends on pacing, as the [§8.7](#87-real-time-pacing) invariant
  requires. (Degrading to boundary detection under pacing was rejected: it would
  move `t*` and diverge paced from unpaced trajectories — row 80.)

This gives step-boundary logic *well-defined semantics*: the transition is defined by
the crossing; detection resolution is an execution-policy detail.

A guard defines a **predicate**: a `Bool`-valued form, or the sign of a
continuous function with **positive = predicate holds** (normative; writing the
guard's sign value `σ`, holding = `σ ≥ 0`). An event fires when its
predicate transitions from
not-holding to holding — edge semantics, uniform across both forms, with the
prior bookkeeping stated in [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event); the opposite crossing direction is
declared as a second event with the negated guard (stall entry/exit as a
pair). Boundary detection handles both forms; localization requires the
continuous (sign) form. This matters in practice: most transitions in FlightPhysics mix
input predicates with state thresholds (e.g. the piston engine's `starting → running`
fires on `ω > ω_idle && fuel_available`).

### 2.2 Exclusions (deliberate)

- **No DAEs / algebraic constraints.** Projection covers the actual need (state
  manifolds) at near-zero cost and zero solver complexity.
- **No SDEs / stochastic integrators.** Noise processes (Dryden/von Kármán turbulence,
  sensor noise) are modeled as ordinary RNG-driven discrete processes (shaping filters),
  which is both faithful to how they are specified and cheap. Consequence elevated to a
  framework guarantee: **deterministic replay** — RNG state lives in component discrete
  state (`x`), never in ambient globals, so same seed ⇒ bit-identical trajectory.
- **No unconditional per-step hook** (no `f_step!` equivalent). Every current use
  decomposes into projection (quaternion renorm) or boundary-detected events (engine phase
  transitions, stall hysteresis latch) — for one class the mapping tightens
  semantics: level-triggered cross-component resets (the gear friction
  regulator under `!wow`) become edge-triggered events ([§15.2](#152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade), [§16](#16-open-axes)).
  Dropping the hook eliminates the footgun of
  model semantics that depend on the integrator's step size.

---

## 3. Component taxonomy

Three classes, two of them leaves with crisp, closed semantics, one pure composition:

### 3.1 Continuous component (the hybrid primitive)

A classical hybrid automaton:

- **continuous state** `x` (isbits struct of real scalars — see [§7](#7-state-and-data-representation)),
- **mode variables** `m`: piecewise-constant values (enums, integers, flags) that
  parametrize the flow and change *only* through event handlers,
- **flow** $\dot{x} = f(x, m, u, t)$,
- **two output stages** (see [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)),
- **events**: guards + handlers (update `m`, may reset own `x`); both read the fresh
  boundary signal table ([§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)),
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

- **discrete state** `x`: any immutable value (see [§7](#7-state-and-data-representation)),
- **update** $x^{+} = g(x, u, t)$ at a declared rate,
- **two output stages** (feedthrough applies at update instants: a proportional path
  is direct feedthrough; a state-only output is not).

State carries the same letter on both tiers — a leaf is strictly one tier (row
56) and no component ever reads another's state, so `x` is never ambiguous: it
is the flow map's argument under `f` and the jump map's under `g`, standard
hybrid notation (row 173). A discrete component's `x` influences continuous
dynamics only **through signals** (outputs held zero-order between ticks).

**`m` is continuous-only**: a discrete component has no
mode store — its FSM enums, flags and counters are ordinary `x` fields. `m`
exists on the continuous side because modes must change *between* flow
evaluations through handlers; on the discrete side `g` already runs at the
only instants anything may change, so a second store would duplicate the
discrete tier's own state semantics under another name.

### 3.3 Assembly

Pure composition: submodels + child connections + boundary faces. **No dynamics of its own.**
Hybridness emerges at the assembly level (an aircraft = continuous vehicle parts +
discrete avionics parts). The two-leaf split held under its strongest
counterexample — a strapdown IMU's periodically-reset integrators land on two
leaves with less code than the fused original ([§15.5](#155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary), row 56). Assemblies are flattened away for scheduling but retained as
the navigation/introspection hierarchy (GUI, logging, paths) and as declaration-level
rate scopes ([§8.5](#85-multi-rate-tick-scheduling)).

---

## 4. Ports and signals

### 4.1 Immutable value semantics

Ports exchange **immutable values** — typically isbits structs (floats, `SVector`s,
enums, nested immutables). The framework owns a **signal table**: one concretely-typed
**cell** per output port in the flattened model. A producer's output-stage function returns
a named tuple of fresh values; the framework writes each into its cell; consumers read
cells. (Vocabulary, binding throughout this document: bare *cell* is the table
entry, and only that — the discrete-state and mode registers are *stores*, [§7.3](#73-discrete-state-modes-and-workspace), not cells;
*staging cell* is a distinct compound term, the per-device inbound register of
[§9.4](#94-inbound-per-device-staging-representation-and-the-drain), which unlike a table cell is mutated frame by frame and sits outside the
table's publish-once discipline; *slot* is reserved for the root input slots of
[§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster).)

Consequences:

- no aliasing, ever — nothing can be mutated under a consumer's feet;
- safe concurrent reads (GUI/logging threads) by construction;
- zero allocation for isbits payloads (named tuples of isbits are isbits);
- each cell has a definite freshness tied to its producer's position in the schedule
  (unlike a monolithic `y` struct, which can be half-fresh mid-sweep with no way to
  tell).

The signal requirement, stated precisely, is **immutability plus frozen references**:
signals may reference bulk data (see [§4.4](#44-function-valued-signals-environment-access)) provided that data is read-only for the
duration of the run. `isbits` is the common case, not the rule.

### 4.2 Consumers see ports, not stages

The port is the addressable unit. A component's outputs appear to consumers, GUI and
logs as one flat namespace (`dyn.vel`, `dyn.f_c_c`, materializable lazily as a view);
which output stage computes which port is a scheduling annotation, invisible outside
the component. Moving an output between stages is non-breaking *for consumers* —
no wire, log or panel sees it; the scheduler does (the feedthrough graph and
stage membership change, [§12.1](#121-three-strata)).

**Visibility.** Which ports exist at all is a declaration-layer decision: the output
contract *is* the public interface, and stage-function results outside it never enter the
table at all — they ride the stage's `w` return down to the component's own later
functions, private by construction, see [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), [§11.3](#113-visibility-the-contract-is-the-interface). (A presentational
*unlisted* flag — skipped in logs and GUI but still connectable — was rejected: it
pretends privacy without enforcing it, and its motivating case, RNG state feeding the
component's own update, dissolves entirely once the update function reads `x`
directly. Row 16.)

### 4.3 Table mechanics and port granularity

- **Scatter/gather is the whole protocol.** A stage function returns a named tuple;
  the framework scatters each field into that port's concretely-typed cell. Every
  reader — the next stage, `f`/`g`, guards, wired consumers, snapshot capture —
  gathers views from cells. The component's aggregate `y` is the merge of its
  stage products (`merge(y_x, y_xu)`, the same on either tier) —
  declared ports only, a stage's private intermediates riding its `w` return rather
  than the table ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)) —
  *semantically* but virtual *physically*: reconstructed per call from cells (field
  loads, register-level, zero cost for isbits), never stored as an object. Name
  collisions across a component's stages are a build error.
- **Stage returns are named tuples of port values, period.** A custom struct is a
  first-class port *value* — one field of the returned tuple, one declared port, one
  cell (`pose = KinPose{T}`). Nested fields get no cells of their own; GUI and logs
  drill into them lazily ([§4.2](#42-consumers-see-ports-not-stages)'s view clause). Bare-struct returns are rejected:
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
  together*. Bundling across dependency footprints is the `KinData` mistake ([§15.1](#151-vehicle-today--this-framework):
  pose is stage 1, velocity-derived quantities are stage 2 — it must split). Fan-out
  is free, so publishing both a bundle and a hot loose field (`pose` *and* `q_eb`)
  is legitimate — one extra isbits cell.
- **Write-side corollary** (from [§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)): **bundle what is written
  together.** The port is the atomic unit of the entire periphery — one cell, one
  root slot, one staged write, one device claim ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)), one trace address, one
  GUI liveness verdict ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)). A component's **ports** are its signal
  endpoints — one cell, one producer. Its **faces** are the names those ports wear
  on the component's boundary: for a leaf the two coincide; for an assembly every
  face aliases an interior port through its boundary declarations
  ([§11.6](#116-paths-wiring-and-faces)) and never creates an endpoint. The
  distinction is kind-blind — wiring and the periphery address a child's faces
  without knowing whether it is primitive or composite.
  Data written by different external writers, or at
  different cadences, must not share a port: pilot commands are scalar faces under
  a namespace prefix, and the convenient bundle is assembled *downstream*, inside
  the graph, by an ordinary component (single producer, consumed together — legal
  by the read-side rule). The two guidelines compose into one principle: a port's
  granularity is set by the finest-grained party owning either end — producers on
  the read side, external writers on the write side. Field-addressed staging (a
  lens into struct slots) stays a recorded guarded addition, unbuilt.

### 4.4 Function-valued signals: environment access

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

**The value-level constructor.** Every field-emitting component must expose the
map (component, input values) → handle as a plain, pure, exported function —
`atmospheric_field(atm; T_sl, p_sl, wind)` for the `SimpleAtmosphere`
successor — and its swept output stage must be a **one-line call to that
function**, never the other way round (the query math written into the output
stage, where only a sweep can reach it). The reason is script-side: [§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)'s
condition math must be able to construct, outside any sweep, bit-for-bit the
same handle the sweep would produce from the same slot values — one
implementation, two call sites, no drift (the silent-drift class [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries) exists to
kill). This is a *shipped component's obligation*, not something a consumer can
retrofit: the real component composes sub-models, and anyone else
reconstructing the map has re-created the drift class. For bulk-data components
the obligation is only that the query math be reachable as a plain function —
they own their resource loading, so building a handle outside a build may cost
a load, which is acceptable because condition authoring is design-time code.

Pre-sampling — a component consuming the field and a pose and emitting plain data
(`Airflow` emitting `AirData` for the whole vehicle) — is an **idiom built on top**,
used where natural; not a separate mechanism. Resource injection (declare-and-resolve
service registries) was considered and rejected for the first cut: it is a second
composition mechanism, invisible in the graph — today's argument-threading, automated.

This replaces threading `atmosphere`/`terrain` as arguments through every update
signature, and dovetails with the terrain ray-query direction of the landing-gear
redesign. Substitutability behind a stable face is declared with an abstract
input entry (`terrain = AbstractTerrainField` — [§11.2](#112-the-declaration-inventory)'s structural
substitutability): the consumer wires to any concrete field type below the
bound, preserving today's `AbstractTerrain` polymorphism at the declaration
layer.

---

## 5. Evaluation order and feedthrough

### 5.1 The scheduling problem

At every evaluation instant, all signals must be computed consistently: every consumer
reads values already produced at that instant. Build the directed graph of (a) wiring
edges and (b) intra-component feedthrough relations; if acyclic, a topological sort
yields a **static evaluation schedule**, computed once at build time. The hot loop runs
a flat list of `(component, stage)` entries — zero runtime graph logic.

### 5.2 Two-stage outputs: signatures, bundles and the hand-off laws

Every component provides exactly **two output stages**, and feedthrough is declared
**structurally, by function signature** — there are no dependency annotations anywhere
in the design:

```julia
# continuous component — maximal legal view set of each bundle in comments
y_x      = h_x(comp, args)          # x, m, t [, ws] — no-feedthrough stage
y_xu     = h_xu(comp, args)         # x, m, u, y_x, w, t [, ws]
ẋ        = f(comp, args)            # x, m, y, u, w, t [, ws]

# discrete component — same stage names, its own closed bundle sets
y_x      = h_x(comp, args)          # x, t, Δt [, ws]
y_xu     = h_xu(comp, args)         # x, u, y_x, w, t, Δt [, ws]
x⁺       = g(comp, args)            # x, y, u, w, t, Δt [, ws]

# every output stage returns y or (y, w) — the return and one-hop laws below

# event system (continuous side only) — same fresh table, same state views:
σ        = guard(comp, args)        # x, m, y, u, w, t [, ws] — Bool or scalar sign value (§2.1)
(; x, m) = handler(comp, args)      # x, m, y, u, w, t [, ws] — keys by the return law below (§12.5)
x⁺       = project(comp, x)         # manifold projection; positional (below)
```

**The hand-off (named bundles).** Every function receives exactly two
arguments: the component and one NamedTuple bundle of zero-copy views, from
which the author **destructures by name** only what the body reads —
`f(c::LowPassFilter, (; x, u)) = ...`, `h_xu(c::PID, (; x, u, Δt)) = ...`. The
executor's call is one fixed shape, `fn(comp, args)`; unread fields are ignored
by language semantics; argument order cannot be confused because there is no
order. Naming spellings rejected, for the record: positional
signatures (dead slots written but unread, un-droppable holes mid-list, and the
`t`/`Δt` scalar pair swappable without error); keyword arguments via
`Base.kwarg_decl` reflection (a load-bearing framework seam on a binding Julia
marks internal — the [§8.1](#81-loop-ownership-the-framework-owns-the-simulation-loop) `task_local_storage` lesson); keyword arguments with
a `_...` slurp (permanent noise, and "the signature is the read-set" weakens to
"at least") — row 74. `project` alone stays positional — one store in, the same store
out, nothing to select.

**The bundle law.** A name appears in a component's bundle **iff the
corresponding store or fact exists for that component**: `x`/`m`/`ws` iff
declared (`init_x`, `init_m`, `workspace`), `u` iff the function family
may see inputs **and** the component declares `input_types`, `y` iff the
component produces any table cell at all (`output_types` ∪ auto-published),
`y_x` iff stage-1 ports exist, `w` iff the stage that hands it down
returned one (the one-hop law below), `t` always, `Δt` on the
discrete tier only. **`y_x` carries the stage-1 *return*, auto-published
names excluded**: an auto-published port is the framework copying a state or
mode field into a cell at stage-1 position ([§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)), and stage 2 already
holds `x`/`m` (continuous) or `x` (discrete) directly, so carrying it in the
hand-down would be transport
for its own sake — [§7.4](#74-the-fused-evaluation-lineage-prior-art-and-how-we-got-here) step 4's rejected identity transport, in a bundle rather
than a table. The rule is what [§12.3](#123-probing-and-input-synthesis) already sources: `y_x` comes from the
stage-1 probe's return, so a component whose only stage-1 ports are
auto-published has no `y_x` in its stage-2 bundle at all (row 169). Undeclared stores are *absent*, never `nothing`-filled:
destructuring a field that is not a thing for you fails at the probe inside the
[§13.2](#132-diagnostics-structured-values-one-carrier-exception) framing diagnostic, with did-you-mean against the legal field set ("`f`
of `Foo` destructures `m`, but `Foo` declares no `init_m`") — one law covering
tier facts, stage legality and declarations alike. The mechanism is structured,
not textual: destructuring an absent field throws a `FieldError` carrying the
type and the field name as data (Julia ≥ 1.12), the probe catches it *matched
against the bundle's own NamedTuple type*, and synthesizes the framing
diagnostic from the legal set — classifying the field as an undeclared store, a
wrong-tier fact, or a name illegal for this function family. The wrong-tier
class survives the state-letter fusion on `m` and `Δt` (row 173): state itself
is legal on both tiers, `m` is continuous-only and `Δt` discrete-only, so
destructuring either on the wrong tier still lands in that bucket. No message text is
scraped, and the bundle stays a bare NamedTuple (row 74); a getproperty-wrapper
spelling is the recorded fallback should type-matched interception prove
insufficient. The per-function, **per-tier** name sets
are **closed**: adding one is a decision-log entry, not a convenience. The
comment in the signature block above states each function's maximal legal set;
a given component's bundle narrows it to declared reality, and the
destructuring narrows further to actual reads — a three-level funnel (stage
name ⊇ bundle ⊇ reads) worth teaching once, because a stateless component
legitimately writes `h_xu` while owning neither `x` nor `m`.

**The stage return law.** An output stage returns either its port NamedTuple
alone, `y`, or the pair `(y, w)`. `y` scatters into the component's declared
cells as always; `w` is a NamedTuple of **private intermediates** —
`isbits` leaves, no cell, no name in any contract, nothing to wire, list or
filter ([§11.3](#113-visibility-the-contract-is-the-interface)). A `nothing` in either slot is a probe error: the pair is a
`Tuple{NamedTuple, NamedTuple}` and padding forms do not exist, for the reason
the handler return law gives below. An empty `y = (;)` *is* legal, so the
**port-less stage** — one whose entire product is `w` — falls out of the
general law instead of needing a rule of its own: stages are discovered by
method existence and stage membership is a partition of the declared ports
that may perfectly well be empty. What is not legal is a stage that produces
neither ports nor `w`: a bare `(;)` computes nothing any consumer can read, and
is a `DeadStage` build error at the probe ([§12.3](#123-probing-and-input-synthesis)) — [§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)'s
inert-component check in the stage register.

**The one-hop law.** `w` travels exactly one hop, to the next function that
could want it, and no further. `h_x`'s `w` flows to `h_xu` if the component
defines one, and otherwise to the downstream set — `f`, guards and handlers;
`h_xu`'s `w` flows to the downstream set. The discrete tier mirrors it exactly
(`h_x` → `h_xu` → `g`), and that is the last time it is said. **Nothing flows
implicitly past its hop**: forwarding a stage-1 intermediate through stage 2 is
an explicit re-return, `(y, (; w..., extra))`, which costs a line and says in
the source that the value crosses — the alternative, a silent pass-through on a
bare `y`, would be the design's one invisible dataflow. Presence in a bundle
follows the producing stage's return under the bundle law: a bare-`y` producer
hands down no `w` key at all, so a consumer destructuring one meets the [§13.2](#132-diagnostics-structured-values-one-carrier-exception)
framing diagnostic naming its legal field set, exactly as for an undeclared
store. `w` is never persisted: the executor hands it down as an ordinary SSA
value **inside one fused pass** — a step fuses the sweep with `f`, an event
round fuses the sweep with its guards and fired handlers ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)) — so freshness
is a property of the construction rather than a rule anyone can violate, and
round fusion is thereby a design constraint on the executor, not an
optimization it may decline ([§12.7](#127-the-compiled-executor)). `w`'s types are probe-observed and its
conformance regime is [§11.3](#113-visibility-the-contract-is-the-interface)'s; the handler return law below is untouched by
any of this — handlers *receive* `w` and return stores.

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
`FieldError`-shaped machinery [§13.2](#132-diagnostics-structured-values-one-carrier-exception) builds for bundles, now running in both
directions ([§12.5](#125-the-always-on-conformance-check)).

The views themselves are unchanged in meaning: own state (on the continuous
tier `x` from the flat buffer and `m` from the mode stores; on the discrete
tier `x` from its store), own published signals (`y`,
gathered from own table cells — the declared ports, [§11.2](#112-the-declaration-inventory)), own private
intermediates (`w`, handed down by the producing stage rather than gathered from
anywhere — the one bundle field with no home at all),
inputs (`u`, gathered from foreign cells
through the wiring's name binding), the clock (`t`, and `Δt` — see [§8.5](#85-multi-rate-tick-scheduling)), and
scratch (`ws`, [§7.3](#73-discrete-state-modes-and-workspace)). The signal table holds only *produced* signals, never
transported ones: each datum has exactly one home — buffer for continuous `x`,
stores for discrete `x` and for `m`, table for signals — and no store mirrors another. Every bundle field
earns its place as a view genuinely readable, and no further "simplification"
exists that does not introduce a copy: eliminating `u` would mean republishing
foreign cells under local names; eliminating `x` would mean identity transport
through the table, rejected for exactly that reason ([§7.4](#74-the-fused-evaluation-lineage-prior-art-and-how-we-got-here), step 4).

### 5.3 Structural feedthrough: stage roles, schedule and step boundaries

**The letters**: `f` is the continuous flow, `g` the discrete
update, `h_*` the output stages — the hybrid-systems flow/jump pair
(Goebel–Sanfelice–Teel) joined to the control/estimation convention that `h` is
the output map (`y = h(x, u)`; every navigation filter's measurement function).
Bare `h` now means the integration step size only ([§8](#8-time-and-execution)), retiring a double
booking. Stage suffixes name the **dependence class**, not the argument
list: `x` versus `xu` — state-only versus
state-plus-input, the `y = h(x)` / `y = h(x, u)` distinction spelled in the
name, identically on both tiers. So "no `u` in the
name" *is* the no-feedthrough property, visible at every definition site.
The letters are deliberately non-exhaustive: modes fold under the state
letter (`m` is state — the objection that `h_x` omits it is answered by the
suffix's job being the feedthrough split, not an inventory; an earlier
`h_xm`/`h_xmu` spelling was traded away for the tier symmetry and the
textbook mirror), and ambient facts (`t`, `Δt`) and scratch (`ws`) ride
unnamed. The stage names do not distinguish the tiers at all — the state
letter is shared (row 173), so what declares a stateful leaf's tier is its
update law, `f` versus `g`, with the remaining tier-implying declarations
agreeing ([§11.2](#112-the-declaration-inventory)). Rejected namings: verb names (`decode`/`compute` — encode
nothing), `Moore`/`Mealy` (exact and opaque), `g`-for-output (forfeits both
the jump-map alignment and the step-size disambiguation). The tier-neutral
stage pair, once rejected for want of an honest neutral state letter, is what
row 173 arrives at: fusing the discrete state into `x` makes the letter
honest on both tiers, and the pair neutral with it.

- **`h_x` is the no-feedthrough stage** — defined entirely by what it
  cannot see: its bundle carries no `u`, so "no feedthrough" is unfalsifiable,
  and that structural guarantee is what its ports contribute to the schedule
  (they break would-be loops). It exists when the component has state-derived
  ports, or shared state-derived intermediates to hand down its `w` return
  ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)); otherwise it is
  simply absent — and a stage that would produce neither is the `DeadStage`
  error, an empty stage being unwritable on purpose. Guidance rather than law:
  when the component also defines `h_xu`, a `w`-only `h_x` earns nothing — fold
  the intermediates into `h_xu`, which runs exactly once per sweep just as `h_x`
  does. The port-less `h_x` earns its keep where there is no `h_xu` at all,
  its no-`u` bundle being the honest spelling of "these intermediates do not
  depend on inputs". **A declared output that matches a state or mode field by
  name and type, and that no stage produces, is auto-published** by the framework
  from the state stores at stage-1 position ([§11.3](#113-visibility-the-contract-is-the-interface)) — the match is against
  the declared stores (`init_x`, plus `init_m` on the continuous tier) and the
  publication position is `h_x` on either tier — publication driven by the
  public contract, rather than the blanket identity publication of state that
  row 16 rejected.
- **`h_xu` receives all wired inputs plus `y_x`** — its own
  stage-1 ports, and with them stage 1's `w`, so shared intermediates are
  computed once, not re-derived, whether or not they are interface —
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
  impossible anyway ([§7](#7-state-and-data-representation)).

The schedule: all stage-1 functions (any order), then stage 2 in topological order, then all `f`
against the now-consistent signal table. Note the systemic consequence: *evaluating
the RHS means running the sweep* — there is no incremental `f`-only re-evaluation.
Implicit solvers, linearization and trim already work this way (seed `x`, run the
composite), so nothing is lost; [§8.3](#83-signal-table-consistency-is-a-boundary-property)/[§8.4](#84-localization-mechanics) restate it as a property of the
execution model (RHS evaluations and guard probes alike run the *interior* sweep,
the continuous-only variant of [§8.5](#85-multi-rate-tick-scheduling) — discrete entries are absent from it by
construction, so discrete cells hold across the step).

**Step-boundary semantics.** At each boundary: integrate → project → boundary sweep →
evaluate **all guards once** against that sweep → fire the eligible events, **at most
one per component per round** (declaration order picking among a component's
simultaneously-eligible events), each firing being `handler → project`. The signal
table is written **only by sweeps**: a transition reaches the table — everyone's ports,
the firing component's own included — at the next round's re-sweep, and the round that
detects quiescence leaves the table post-transition-consistent for whatever else the
boundary does (discrete ticks, logging). Hence the epoch rule: **a handler executes
against exactly the world its guard fired on** — own `y`, foreign `u`, own `x`/`m`
alike are the firing round's sweep, so `y = h(x)` holds at every handler entry.
Same-component sequential composition happens *across* rounds, each later event
re-decided against the post-transition sweep rather than fired on a stale premise.
Newly-enabled guards fire within the *same* boundary: the
sweep → guards → handlers phase iterates to quiescence, with each event firing at
most once per boundary and each component firing at most once per round (settled in [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)).

**Departure from the orthodox formalism, stated openly.** The textbook form is
$\dot{x} = f(x, u)$, $y = g(x, u)$; this design's `f` receives the orthodox arguments
*plus* the published table: $\dot{x} = f(x, m, y, u, t)$. The composite map $x \mapsto \dot{x}$ is
mathematically identical (linearization, trim and AD are untouched); the heterodox
element is only that derivatives may read outputs. The teaching line: *"stage 1
publishes what you know from state alone; stage 2 adds what needs inputs; your
dynamics read your own published results instead of recomputing them."* The decision
was grounded in a component-by-component survey of FlightPhysics/FlightApps ([§15.2](#152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade)):
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

### 5.4 Artificial loops and the escape hatch

A component that bundles a no-feedthrough output with a feedthrough output in one
atomic evaluation unit can be **port-level acyclic yet unschedulable** (Simulink's
"artificial algebraic loop"). The canonical instance in this domain is rigid-body
dynamics: velocity out (pure state) + acceleration out (feedthrough from total force).
The two-stage split resolves it, and it is the rung that absorbs most of the class:
[§15.1](#151-vehicle-today--this-framework)'s `VehicleDynamics` instance — velocity state-only, accelerations feedthrough —
simply dissolves under it.

What survives the split is the case where a single component's stage-2 outputs
cross-couple through a neighbor (port-level acyclic, stage-level cyclic), which [§5.6](#56-diagnostics-feedthrough-tracing)'s
tracer labels **artificial**. Two remedies, in this order:

- **Re-factor the contract.** Before moving any code, re-examine the cycle's wires. An
  input the neighbor consumes *only in a fallback branch* is the archetypal false
  dependency: the neighbor is computing, on the component's behalf, a fallback whose
  semantics belong on the component's own side of the boundary. Move the branch to its
  natural owner and the wire disappears. The canonical instance is the landing gear's
  strut/steering pair: the steering model consumes the contact-point velocity azimuth
  `ψ_v` only in its disengaged (castoring) branch, but castoring is free-swiveling
  wheel physics — the strut's business, not the steering law's. Re-factoring the
  steering contract to emit `(engaged, ψ_cmd)` and computing
  `ψ_sw = engaged ? ψ_cmd : ψ_v` inside the strut deletes the backward wire outright.
  The factoring survives substitution, which is the test that it records structure
  rather than dodging the diagnostic: a stateful steering actuator produces `ψ_cmd`
  from its own state and still needs nothing from the strut ([§16](#16-open-axes) records the
  migration).
- **Split the component** — the residual remedy, when both halves genuinely belong to
  it, and it documents real structure. Its cost, stated where it bites: [§11.3](#113-visibility-the-contract-is-the-interface)'s
  visibility is binary, so every intermediate shared across the new boundary becomes
  `output_types` — public, connectable, substitution-relevant. The mitigating idiom is
  [§4.3](#43-table-mechanics-and-port-granularity)'s granularity guideline, which the split case satisfies trivially (one
  producing stage, one consumer): **one struct-valued bundle port** — a
  `StrutGeometry`-shaped value — not N loose ports. The bundle type is then contract,
  a real cost but a bounded and honest one. No visibility register is added for the
  orphaned intermediates: rows 34 and 55 (`unlisted`, `Private(T)`) stay closed, and
  the `w` channel ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)) is no exit either — it hands values between one
  component's own functions, and there is nothing for a wire to carry.

The build diagnostic offers both exits explicitly ("cycle through `systems/aero` is
artificial at port level — split the component, or narrow the neighbor's contract",
with the offending stage `h_xu` carried as a separate payload field rather than dotted
onto the path, [§11.6](#116-paths-wiring-and-faces)/[§13.2](#132-diagnostics-structured-values-one-carrier-exception)). The split is rare, and the ladder is what earns the word
rather than asserting it: the two-stage split dissolves the common shapes and the
contract re-factoring absorbs the false wires, leaving the split for cycles whose
halves really are one component's own work. One consequence of stage-2 conservatism
worth recording: an input consumed only by `f` (never by `h_xu`) still creates a
scheduling edge if the component has stage-2 outputs; in practice such components are
integrator-shaped and have none, and the remedy, if ever needed, is the same ladder.

### 5.5 Algebraic loop policy: reject at build time

A genuine cycle in the instantaneous dependency graph is a **build error** with a
diagnostic naming the full path in the canonical slash form of [§11.6](#116-paths-wiring-and-faces)
(`aero/F → dyn/a → aero/α̇ → aero/F`). The user breaks
it explicitly: insert dynamics (the α-filter idiom — already standard practice in the
domain and in the current C172 model), insert an explicit unit delay ([§13.7](#137-tooling-consequences-provenance-and-the-component-library)'s
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

### 5.6 Diagnostics: feedthrough tracing

Tracing is **diagnostic only, never load-bearing**: scheduling correctness comes
exclusively from the structural two-stage split; tracing improves error messages and
verification. Triggered when the scheduler finds a cycle, to classify it (genuine →
"insert a state"; artificial → [§5.4](#54-artificial-loops-and-the-escape-hatch)'s remedy ladder).

**Detection and naming.** A cycle surfaces as a topological-sort stall in
Stratum B. The stall's residue is not the diagnostic: it also holds the innocent
downstream cone, so it over-reports, while a single DFS back edge under-reports —
one edge of a possibly large tangle. The stalled subgraph is instead decomposed
into **strongly connected components**, and each nontrivial SCC names one cyclic
cluster exactly: one diagnostic, its members and the wires among them, presented
as one readable loop in [§11.6](#116-paths-wiring-and-faces)'s canonical slash form
(`aero/F → dyn/a → aero/α̇ → aero/F`).

**Classification is schedule-free.** It runs inside Stratum B's failure path,
where no schedule exists — and needs none, because each SCC member is evaluated
*once, in isolation*, at the probe point: state views from `init_*`,
out-of-cycle cells from the acyclic prefix's probe values, in-cycle cells
synthesized through [§12.3](#123-probing-and-input-synthesis)'s `probe_value` under tracer tags. The tracer's
product is a per-member dependence set rather than a value, so no ordering has
to be valid for the labels to come out right. The loop is **real** iff every hop
of the structural cycle survives in the traced per-member maps; **artificial**
([§5.4](#54-artificial-loops-and-the-escape-hatch)) iff some hop dies — the component whose stage-2 function does not in fact
route that input to that output. No Stratum C machinery is touched: no
activation, no layouts, no table. This is row 12's *local* variant — the
schedule-free per-member trace at the probe point, which is what the cycle
classifier uses; [§12.4](#124-activations-executable-sets-laziness-caching)'s "tracer activation" names row 12's other variant, the
global set-tracer run as an ordinary Stratum-C activation, and the two must not
be conflated.

**Caveats, carried in the diagnostic rather than assumed away.** The trace
speaks for the branch taken at the probe state (row 12's diagnostic-only
doctrine). Discrete members trace *structurally* — the discrete tier's plain,
wholesale-pinning declarations admit no tracer scalar — which is sound as a
may-depend answer but never sharp,
so the remedy hint — split this component, *or* narrow the neighbor's contract when
the dead hop's input is consumed only in a fallback branch ([§5.4](#54-artificial-loops-and-the-escape-hatch)'s ladder) — is
offered only for continuous members. And
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

## 6. Composition: connections, aggregation and hierarchy

Components become a system through wiring: connections that route signals
across the assembly hierarchy, and ordinary junction components wherever
several signals must combine into one. [§6.1](#61-connections-and-hierarchy) gives the connection and
hierarchy rules; [§6.2](#62-aggregation-explicit-summing-junctions) gives the aggregation idiom they force.

### 6.1 Connections and hierarchy

- **Deep connection paths** are allowed, with one structural rule: both endpoints must
  resolve **within the assembly type being defined**. You may deep-route into structure
  you declared (`Cessna172` routing its single `trn` input to
  `systems/ldg/{left,right,nose}/trn_field` in one visible block — no per-level
  re-exports); you may only connect port-level into submodels held **generically**
  (`World` connects `terrain/field => aircraft/trn` and knows nothing more). This kills
  the re-export ceremony where it is ceremony, and preserves substitutability where the
  boundary is load-bearing. Operationally: a path may traverse any chain of
  concretely-typed fields and stops at the first generically-held child, whose
  faces are the only things addressable beyond that point — a rule about
  the *declaration's* knowledge, not the build's (a deep path into a generic child
  is forbidden even where the concrete instantiation would resolve it, because it
  hard-codes one implementation and breaks on substitution). Enforcement lives
  in the path-resolution primitive itself (`resolve`, [§13.3](#133-build-primitives-resolve-and-the-face-list-accessors)), which walks
  declared field types alongside instances.
- Paths are validated at build time; renames break loudly.
- **Two clauses type-check a wire** ([§11.2](#112-the-declaration-inventory)). The **nominal bound check** is the
  standing rule, stated over declaration evaluations: the producer's declaration
  at `Float64` must be `<:` the consumer's entry at `Float64` — one uniform rule,
  degenerating to exact equality for a concrete entry, violated as
  `WireTypeMismatch`. Beside it, and **for a continuous consumer only**, the
  **walk-compatibility clause**: a walking producer leaf (the producer declared
  `T` there) requires a `T` entry, while a pinned producer leaf satisfies either,
  frozen values embedding upward under any activation. Both sides are declaration
  functions of `T`, so the clause is decided in Stratum A by evaluating them at a
  marker scalar — declaration reading, no user stage code ([§12.1](#121-three-strata)) — and a
  violation is `WalkingFaceAtFrozenEntry`, naming both endpoints, the leaf and
  both declared leaf types, with both remedies in the message: declare the entry
  `T` if the consumer promotes, or feed it from a non-walking source if the freeze
  is genuine. **The tier scope is load-bearing, not tidiness.** A discrete consumer
  takes the bound check alone, because its stages read exclusively at real ticks
  in the nominal world — a `Dual`-carrying cell exists only inside activations the
  discrete tier never runs in ([§12.4](#124-activations-executable-sets-laziness-caching)) — so a continuous producer feeding a
  discrete consumer is unconditionally legal; unscoped, the clause would reject
  every sensor → controller wire in the design. The clause is also what gives the
  two contract sides their **failure asymmetry**: the input-side forgotten-`T` —
  the habitual `Float64` written at an entry whose consumer really promotes —
  fails at the *first nominal build*, at the wire, with both endpoints named,
  because an input has a build-time counterparty; the output side has none, so its
  forgotten-`T` lurks (loudly, never silently) until the first `Dual` activation
  ([§11.2](#112-the-declaration-inventory)).
- Fan-out is free (one producer, many consumers). The converse is strict: every
  input port takes **exactly one** connection, no exceptions (aggregation is
  junctions, [§6.2](#62-aggregation-explicit-summing-junctions)). The rule spans levels: an input fed both inside a sub-assembly
  and by an ancestor's deep route is a two-producers build error — deep routing
  cannot silently double-feed.
- No auto-bubbling of unconnected inputs (implicit interface growth — and worse:
  the forgotten-wire error, [§11.4](#114-failure-walkthroughs-the-error-locality-grounding) walkthrough 2, would be silently *promoted to
  interface*, climbing level by level into a live root slot nobody feeds, instead
  of failing at build).
- Unconnected output ports: legal, silently (there is no build-time warning
  for them, row 84 — under mandatory `output_types` most models carry many
  observation-oriented ports no wire consumes, [§9.2](#92-outbound-snapshot-publication) blesses exactly that,
  and [§9.2](#92-outbound-snapshot-publication)'s path-bound readers attach after the build, so "unused" is not
  even decidable; a warning firing on every honest port poisons the sole
  warning stream, [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)'s anti-diagnostic lesson. The hazard it nominally
  guarded — a wire someone meant to draw — is caught from the consumer side,
  where the information actually lives). Unconnected input ports: build error
  (no silent defaults). **The check is a whole-tree property, not a
  per-declaration one**: within a single assembly declaration an unfed child input
  is simply *awaiting a claim from above* — a sibling wire, an ancestor's deep
  route, or an `input_connections` entry handing the obligation up one level ([§11.6](#116-paths-wiring-and-faces)). The
  error fires at the root build for any input whose obligation chain never
  terminates. The one legitimate terminus fed by no component is the root
  assembly's own input face — a root slot ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)).

### 6.2 Aggregation: explicit summing junctions

N-to-1 physical aggregation (total wrench, total mass properties, total internal
angular momentum — today's generated `get_wr_b`/`get_mp_b`/`get_hr_b` tree walks) is
expressed by **ordinary junction components and explicit wires**. There is no
framework aggregation mechanism: no multi-connection ports, no declared fold ops, no
identity-element opt-outs. Every input port takes exactly one connection, everywhere.

```julia
struct SumJunction{W, N} end        #type constructor, arity; library-provided

input_types(::SumJunction{W, N}, ::Type{T}) where {W, N, T <: Real} =
    NamedTuple{ntuple(i -> Symbol(:in, i), N)}(ntuple(_ -> W{T}, N))
output_types(::SumJunction{W, N}, ::Type{T}) where {W, N, T <: Real} = (; Σ = W{T})
h_xu(::SumJunction, (; u)) = (; Σ = +(u...))
```

(The parameter is the *unparametrized* type constructor — `SumJunction{Wrench, 3}`;
UnionAlls are legal type parameters — so both contracts derive their entries from
it by applying it to the activation scalar: the junction is a continuous leaf, so
its `input_types` entries are the tolerant `W{T}` a promoting consumer writes
(walking, frozen and root-slot contributors all admissible behind them) while
`output_types` re-types the output cell per activation ([§11.2](#112-the-declaration-inventory)). This is the same
arity-via-computed-contracts pattern [§13.7](#137-tooling-consequences-provenance-and-the-component-library) commits to for `Or{N}`.)

Wired at an ownership boundary, the junction is ordinary structure — with
`wr_sum::SumJunction{Wrench, 3}` a field of `Systems` like any other child
([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)):

```julia
child_connections(::Systems) = (
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
  (`input_types(::VehicleWrenchSum, ::Type{T}) where {T <: Real} = (aero = …, ldg = …, pwp = …)`) documents the
  contributor set better than generated slots, at the price of hard-coding it into a
  type; the generic positional form remains the tool for configuration-variable
  sites. Both are plain components; the framework is not involved.

**The hierarchical aggregation idiom** (what replaces the tree walk). Only physical
contributors publish these ports — a strut publishes `wr_b`, avionics publishes
nothing — and each assembly that *owns* contributors aggregates them with an internal
junction and **exports the total** ([§3.3](#33-assembly): the junction is a component inside the
assembly; the assembly exports its `Σ` port). The [§6.1](#61-connections-and-hierarchy) connection rules force this
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
standard component library committed in [§13.7](#137-tooling-consequences-provenance-and-the-component-library): ordinary components, no
framework privileges, inventory grown strictly by migration demand.

**The zero-contributor end of the same spectrum.** Ragged contributors bottom out at
none: a configuration in which a consumer's required aggregate input has *no* physical
contributors at all — the bare-propagation `Vehicle{NoVehicleSystems}`, zero
contributors to external wrench and to internal angular momentum, while
`VehicleDynamics` requires both unconditionally. There is no junction to write and no
producer to wire, [§6.1](#61-connections-and-hierarchy) bans unconnected inputs and silent defaults, and the identity
element a zero-arity junction would need was retired on purpose (row 37). The spelling
is a library `Constant` source ([§13.7](#137-tooling-consequences-provenance-and-the-component-library)) wired straight to the consumer's input —
`Constant(Wrench())` → `dynamics/wr_ext` — so the zero total becomes declared
structure, the configuration stating "external wrench ≡ 0" as a visible wire and an
observable port rather than as an identity method the framework supplies behind the
author's back. This is not [§6.1](#61-connections-and-hierarchy)'s banned default in component clothing but its
opposite: that default is silent and consumer-declared, this one is loud and
assembly-declared, the author writing the child and the wire, both inspectable.

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
assembly; the mandatory typed declarations of [§11](#11-the-declaration-layer-components-and-assemblies) make junction mistakes loud,
dissolving the original "positional ceremony" objection; and `Reduce` was the
declaration vocabulary's last remaining wrapper — killing it leaves `input_types` as bare
types with zero framework vocabulary and retires the canonical-fold and
identity-element rules wholesale.

---

## 7. State and data representation

### 7.1 Continuous state: structured immutable, flat backing

Each continuous component declares its state by value (`init_x`, [§11.2](#112-the-declaration-inventory)): a
NamedTuple whose leaves are drawn from a **deliberately closed vocabulary —
plain real scalars and `SArray`s (static vectors and matrices) of a common
eltype `T`** — and nothing else. `Int`s/enums/`Bool`s belong in modes, and
domain wrapper types (`RQuat`, `Ranged`) are not state leaves — an attitude
state is an `SVector{4,T}`, cast where rotation semantics are wanted (below).
The framework:

- computes a **flat layout** at build time (compile-time offsets over one contiguous
  `Vector{T}` buffer it owns);
- **reconstructs** the typed immutable state value for a component at each evaluation
  and passes it to every function receiving state views ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws) argument rule — field
  loads at known offsets, register-level, zero cost);
- receives immutable results back: derivative functions return an `Ẋ`-typed value
  (scatter-stored into the flat `ẋ` buffer); event handlers and projection return a new
  `X` (written back — projection at [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)'s two schedule positions).

**What `Ẋ` is.** With the leaf vocabulary closed, the answer takes one line:
`Ẋ` has exactly `X`'s shape at the activation scalar — a scalar leaf's
derivative is a `T`, an `SArray` leaf's is the same `SArray` at `T`. (This is
the vocabulary rule paying rent: an invariant-carrying leaf like a unit
quaternion has a derivative off its own type, and `Ẋ` would need a separate
derivation; here the attitude leaf is an `SVector{4,T}` and so is its rate.)
The conformance predicate is structural — *each field of `f`'s return
scatters into its field's block at `T`* ([§12.5](#125-the-always-on-conformance-check) states the check) — which
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
buffer-unchanged-within-a-sweep rule ([§12.7](#127-the-compiled-executor)). The complementary rule ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)): **one home per datum** —
buffer for continuous `x`, stores for discrete `x` and for `m`, table for produced signals — and no store ever
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
it and the derivative keeps pushing, `project` ([§2](#2-formalism)) redundant on views yet
necessary on the buffer, [§8.4](#84-localization-mechanics)'s off-manifold probes (which must see the raw
interpolated state, RK-stage-like) impossible; or bypass the constructor,
which is safe under the common-eltype rule and build-time verifiable, but
plants layout-coincidence cleverness in the executor's core to save a cast
line. Domain semantics are instead an **explicit, invariant-free cast at the
point of use** — `q = RQuat(x.q, normalization = false)` — the conversion
today's `f_ode!` code performs on its raw views, now visible and chosen.
Invariants live where the design already put them: in `project` at
boundaries, and in writers — handlers build their returned values through
ordinary constructors, and the condition apply converts authored values
through ordinary `convert` methods ([§14.3](#143-resolution-flatten-validate-compile-once)). Constructors run on the write
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

### 7.2 Numeric genericity (eltype)

The state buffer, pack/unpack machinery, and the entire **continuous evaluation path**
are generic over `T <: Real`. One design property, four consumers:

1. exact Jacobians for **linearization** (ForwardDiff duals through the whole model,
   replacing finite differences),
2. derivatives for **trim** solvers,
3. the **feedthrough tracer** ([§5.6](#56-diagnostics-feedthrough-tracing)),
4. a trivially checkable **CI invariant**: one evaluation sweep with `T = Dual` fails
   loudly (`MethodError`/`InexactError` at the offending line) on any Float64-pinning.

Consumer 1 is also where the *discrete* side's exemption stops being a
limitation and becomes the exact answer: a frozen discrete cell is a constant
with zero partials, which is what "linearize the continuous dynamics with the
discrete state held" means (`frozen_discrete_walkthrough.md` works the chain through in
detail).

The declaration layer keeps this scoping legible without putting it in the
author's way: a continuous producer's output declaration is a function of the
activation scalar ([§11.2](#112-the-declaration-inventory)), and cell types per activation are that declaration
*evaluated* at the scalar — participation authored per leaf, `T` where the leaf
follows the activation and a concrete type where it is deliberately pinned. The
state type is still derived: the framework walks the `init_x`-derived type, real
leaves and `Real` type parameters following the scalar. The discrete side stays
plain and pins wholesale. Nothing anywhere comes from inference through user
code. (Safety of the substitution rests on
the [§12.5](#125-the-always-on-conformance-check) embedding guarantee.)

Scoping (what actually needs genericity — roughly half the type inventory):

- **Walked — payload/value types constructed during evaluation (~25 structs):**
  the quaternion/attitude family (`Quaternion` becomes
  `Quaternion{N,T} <: AbstractVector{T}` — by invariance, `Float64` instances still
  match every existing `AbstractVector{Float64}` method, so existing behavior is
  untouched), `Wrench`, `FrameTransform`, `MassProperties`, `KinData`, `AirData`,
  geodesy value types, `TerrainData`, continuous output structs. Mechanical
  parametrization; constructors infer `T`, so call sites don't change; `@kwdef`
  defaults pin the no-argument case to `Float64`.
- **Pinned — parameters and definitions:** stay `Float64` (promotion handles mixing);
  no migration.
- **Exempt — the discrete side** (compensators, avionics): linearization and
  trim differentiate continuous dynamics only.

Lookups: **table data is a pinned parameter; the query coordinate is walked traffic.**
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

### 7.3 Discrete state, modes, and workspace

- **Storage:** discrete state (a discrete leaf's `x`) and modes (`m`) live in
  **typed stores** (the same immutable-value
  discipline as the table's cells, in a separate home — [§4.1](#41-immutable-value-semantics)'s vocabulary),
  overwritten by the framework when an update/handler returns a new value. They never
  touch the integrator buffer; no arithmetic is ever done on them.
- **Type freedom:** a discrete leaf's `x` may be *any immutable value* (frozen-reference rule; isbits not
  required). Enums, integers, nested structs, RNG state (four `UInt64`s of `Xoshiro` —
  required to live in discrete state for deterministic replay).
- **Snapshots are free:** copying a store copies a reference to immutable data —
  checkpoint/replay of the entire discrete side is "copy the store values."
- **Workspace** (for heavy algorithms, e.g. an n≈20 Kalman filter):
  component-declared mutable scratch, instantiated by the framework and arriving
  as the `ws` bundle field ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)) in every bundle-receiving function of the
  declaring component (`project` is positional and receives none),
  and **excluded from state semantics** — not snapshotted, not replayed, never a
  condition target ([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)), must carry no information between calls. Checkable:
  in debug mode the framework **poisons the workspace before every call**, so
  read-before-write of stale scratch detonates immediately. The poison is scoped
  to what a store can carry: float-eltype stores get `NaN`, integer-eltype
  stores get `typemin` (which detonates on first use as an index), and element
  types with no sentinel — `Bool` flags, enums, opaque handles — are skipped,
  the skip reported once per activation. Naming the skipped stores is the point:
  a guarantee that is silently partial is exactly the anti-diagnostic pattern
  [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract) rejects. Declared **by
  allocation**: the well-known method *is*
  the allocator — `workspace(c::KF, ::Type{T}) where {T} = (P = Matrix{T}(undef,
  c.n, c.n), x̂ = Vector{T}(undef, c.n))` on the continuous tier, plain
  `workspace(::C)` on the discrete (the contract declarations' tier split) — called
  once per activation and once per scratch-store set ([§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)), sizes from the
  instance, eltypes from the activation. The `T`-signature was never in doubt
  here — `workspace` is the by-allocation register, so the method is an
  allocator the framework *calls*, not a schema it reads — and it is the
  precedent [§11.2](#112-the-declaration-inventory)'s register criterion cites now that row 166 has restored the
  scalar to the by-type register as well: a `T` appears in a signature exactly
  where the author makes a choice with it. Nothing
  downstream derives from a workspace's type, and mistyped scratch detonates
  loudly at the `Dual` probe. The `undef` spelling is the
  recommended idiom: it makes "contents are meaningless" visible in the
  declaration, which is the register this store actually lives in — a
  workspace is *not* memory, so declaration-by-initial-value ([§11.2](#112-the-declaration-inventory)) never
  legitimately covered it, and none of the [§11.2](#112-the-declaration-inventory) arguments for authored values
  (condition overlay base, the probe-value barrier) applies to a store that
  conditions exclude and the poison overwrites. Available on **both tiers**:
  nothing in the contract is tier-specific, and a continuous workspace simply
  joins the `T`-generic surface — under a `Dual` activation the allocator is
  called at `Dual`, and the in-place math runs through Julia's generic
  fallbacks (no BLAS; activations probe and linearize, they don't run
  marathons). The calls-per-boundary multiplicity of the continuous side (RK
  stages, localization probes, event re-sweeps) makes the no-information-between-
  calls contract *more* load-bearing there, not less — which is the poison
  check earning its keep, not a prohibition. (Rejected: `workspace_type`
  returning types with framework instantiation — array types carry no
  dimensions, sizes live in runtime fields like `kf.n` which no zero-argument
  constructor can see, and hoisting sizes into type parameters lands on the
  `MMatrix` codegen catastrophe below.)
- **Blessed idiom — zero-allocation ticks with immutable `x`:** do the in-place math
  (`mul!`, `cholesky!`, BLAS) on the workspace; at the end, snapshot into an isbits
  container (`x = KFState(SVector{20}(ws.x̂), SMatrix{20,20}(ws.P))`) and return it.
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

### 7.4 The fused-evaluation lineage (prior art and how we got here)

The [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws) interfaces are the end point of a three-step simplification arc, recorded
here because each step replaced a mechanism with something smaller:

1. **N output groups → exactly two.** General output groups (each declaring its input
   subset) handled a cross-coupling case that never materialized in the domain; strict
   two-stage eliminates all dependency declarations, at the price of an occasional
   component split ([§5.4](#54-artificial-loops-and-the-escape-hatch)).
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
   anyway" — quietly depended on default identity publication; once [§11.3](#113-visibility-the-contract-is-the-interface) made
   publication a deliberate interface act, the identity decode stood revealed as
   *transport*: copying the buffer into cells so a buffer view could be replaced by
   a cell view — ceremony of exactly the kind step 2 removed elsewhere, now with
   dead stores (no own-function reads state through the table). The
   reductio that decided it: the same minimalist logic would pack `u` into the stage-2 product
   and end at `f(comp, y, t)`, republishing foreign cells under local names. The
   fixed point is [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)'s argument rule — zero-copy views of the stores a function
   genuinely reads. What survives of step 3: the uniform shapes, the fused
   economics, and the stage-1 decoder itself (today's `h_x`) — no longer the sole state gate, but the
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

### 7.5 Allocation policy: a scoped invariant

Not dogma — three reasons, only one about speed: (1) GC-pause jitter control for
real-time; (2) throughput for unattended runs; (3) **the canary**: an unexpected allocation
is Julia's most reliable symptom of type instability, so a zero baseline makes
`@allocated == 0` a CI-testable invariant that catches inference regressions at the
offending commit.

- **Continuous hot path** (per-stage evaluation), plus everything else that
  runs unconditionally per frame or boundary — guards (evaluated every
  boundary, firing or not) and `project` (both [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries) schedule positions):
  exactly zero, CI-enforced at the [§12.7](#127-the-compiled-executor) phase-body seam (`phase_bodies`).
- **Periodic ticks and event handlers** (episodic execution — a tick when
  due, a handler only on firing): zero by idiom (workspace + snapshot
  pattern; immutable-value returns); documented tolerance for the rare
  exception, scoped per body by the seam's granularity so it never loosens
  the continuous assertions.
- **Logging**: amortized-zero — snapshots are records stored *inline* in a
  `Vector`, `sizehint!` for the expected duration making regrowth a non-event.
  The inline-storage claim is about the snapshot record's *fields*, not about
  everything reachable from them: a model carrying [§4.4](#44-function-valued-signals-environment-access) field handles (heightmap
  terrain, wind grids) has a snapshot type with reference fields, and those
  ride as references to build-time-frozen data — no copy, no per-boundary
  garbage, which is what the allocation claim asserts. What the claim does not
  assert is that the snapshot is `isbits`; the per-boundary allocation cost is
  zero either way, and [§4.4](#44-function-valued-signals-environment-access)'s summarize-or-skip rule governs what such a field
  contributes on export.
- **What is not recorded**: event firings. The log holds boundary snapshots and
  the trace holds staged inputs ([§9.2](#92-outbound-snapshot-publication), [§9.5](#95-inbound-the-input-trace)); neither carries a per-event
  record. Which events fired at which boundary is recovered by replay plus the
  published modes — the honest remedy of [§9.2](#92-outbound-snapshot-publication), a mode field declared public is
  in every snapshot. An event-firing stream is a guarded addition, not built.
- **Tools where garbage is unavoidable**: arena allocation (Bumper.jl-style) for scoped
  temporaries; scheduled `GC.gc(false)` at frame boundaries to move collection out of
  the critical path. (Julia has no per-object freeing; these are the honest levers.)

---

# Part II — Execution

## 8. Time and execution

### 8.1 Loop ownership: the framework owns the simulation loop

The simulation loop — the [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries) boundary sequence, tick dispatch, event handling,
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

### 8.2 The stepper seam

The one delegated operation is *advance the continuous state from `t` by `h`*, behind
a narrow internal interface (the **stepper seam**). Its contract:

- **advance by arbitrary `h`** — required anyway to land on tick boundaries and to
  resume from a localized event time;
- **dense output on demand over the last completed step** — required only by event
  localization ([§8.4](#84-localization-mechanics)), constructed lazily;
- **one-step methods only** — event handlers reset state discontinuously, and a
  one-step method restarts from a new state for free; multistep methods would need
  history-rebuild machinery after every handler and are excluded.

**The seam is never entered empty.** A model with no continuous state at all is
legal — nothing in [§11.2](#112-the-declaration-inventory) requires an `x` block of anyone — and the framework
short-circuits rather than pushing the corner down the seam: with an empty `x`,
integrate degenerates to advancing `t` to the next boundary, and the stepper is
simply not called. No backend ever faces `N = 0`, and no backend contract has to
say what it would do there. This finishes structurally what [§8.1](#81-loop-ownership-the-framework-owns-the-simulation-loop) argues — the
dummy-`[0.0]` tax it charges to FlightCore is gone at the root, not just the
buffer but the step over it. Everything else about such a model is ordinary: the
boundary machinery — sweeps, events, ticks — runs unchanged.

First cut ships **in-house fixed-step RK4 and Heun** over the flat state buffer:
~a hundred lines, trivially zero-allocation (auditable for the [§7.5](#75-allocation-policy-a-scoped-invariant) CI invariant),
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
   If that day comes, [§7.2](#72-numeric-genericity-eltype) genericity supplies exact ForwardDiff Jacobians through
   the sweep for free.

### 8.3 Signal-table consistency is a boundary property

During a step, RK stages evaluate the interior sweep ([§8.5](#85-multi-rate-tick-scheduling)) at internal
stage states — the signal table is transiently **integrator scratch**. The boundary sweep in the [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries) sequence
is what restores consistency at each accepted boundary. The rule, binding for the periphery ([§9](#9-runtime-periphery-the-data-plane)):
**external readers (GUI, logging, network output) observe the signal table only at
step boundaries.** Mid-step contents carry no meaning.

### 8.4 Localization mechanics

Trigger: a localized event's predicate was not-holding at $t_n$'s quiescence (its
prior, [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)) and is holding at $t_{n+1}$ — the directional edge of [§2.1](#21-events-two-detection-policies), never a
bare sign change: a holding → not-holding transition neither fires nor localizes.

- **Interpolant (lazy).** Build the cubic Hermite continuous extension $\hat{x}(\theta)$,
  $\theta = (t - t_n)/h \in [0, 1]$, from $(x_n, \dot{x}_n, x_{n+1}, \dot{x}_{n+1})$;
  $\dot{x}_n$ is the step's first
  stage, $\dot{x}_{n+1}$ costs one sweep, paid only on trigger. Uniform accuracy $O(h^4)$ — one
  order below the discrete solution, the standard pairing, and the event time can
  only ever be as accurate as the interpolant, so nothing more expensive is worth
  probing.
- **Probes run the interior sweep.** Guards read `y`, so evaluating a guard at an interpolated
  state means writing $\hat{x}(\theta)$ into the state buffer and running the interior sweep — the same
  rule as the RHS ([§8.5](#85-multi-rate-tick-scheduling)): a probe is a mid-step evaluation, so discrete cells hold
  their tick values through localization and a guard reading a sampled output sees
  what the controller is holding. One interior sweep per probe.
- **Root-finding: bracketed and derivative-free** (ITP or Brent; bisection is an
  acceptable fallback). The observed not-holding/holding bracket *is* an unconditional
  convergence certificate. Newton/AD localization is rejected: guards are guaranteed C⁰, not C¹
  (clamps, table knots, saturated stretches where σ′ = 0), Newton discards the
  bracket for merely local guarantees, and its superlinear convergence saves a
  handful of microsecond probes per rare event. AD earns its keep in Jacobians, not
  in root-polishing a possibly-kinked bracketed scalar. **Convergence is a
  relative bracket width**: localization stops once the bracket is narrower
  than `localization_tol · h`, with `localization_tol` a `Simulation` deployment
  keyword defaulting to `1e-6` (below). Relative because an absolute-in-`t`
  tolerance is not scale-free — one number is slack at `h = 1` and unreachable at
  `h = 10⁻⁴`; `1e-6` because the event time can only ever be as accurate as the
  interpolant (`O(h⁴)`, above), so at practical `h` anything tighter buys nothing
  while every probe costs a full sweep — a handful of probes under ITP, ~20 even
  in bisection's worst case.
- **Post-event.** The boundary sequence runs at `t*` (below) → **interpolant
  invalidated** (the handlers made it a lie for `t > t*`) → resume integration from
  `t*` with the remainder step targeting `tₙ₊₁` → re-check guards on the remainder,
  under a bounded per-step event budget with a chattering diagnostic. Multiple
  events localizing in one step fire at the *earliest* `t*` (ties fire
  together at that boundary, declaration order within the iteration, [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event));
  later crossings re-localize on the remainder.
- **Shared blind spot, documented:** an even number of crossings within one step
  returns the predicate to not-holding at the boundary, so no edge is observed —
  defeating detection under both policies; the mitigation is step size, not
  machinery.

**Budget exhaustion degrades; it does not throw.** The budget is
`event_budget`, the second deployment keyword this section fixes: an integer
count of localizations permitted within one frame, defaulting to **8**. A
legitimate multi-event frame — three landing-gear struts touching down inside
one step — needs three or four; chattering needs tens; 8 bounds the pathology
without ever binding on a healthy model. When a step spends its event
budget, localization stops for the remainder of that frame: the remainder step
completes, and any further crossings fire in the next boundary's ordinary
iteration — boundary granularity for that frame — under a warning naming the
chattering event and the localization count. The degradation is a function of
the trajectory alone, never of wall clock, so row 80's pace-independence stands
untouched and the run replays identically. A `StepError` would misclassify an
expected modeling outcome as broken machinery ([§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)'s doctrine). This also
makes the asymmetry with [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) principled: the once-per-event bound *there* is
structural — quiescence terminates because each event fires at most once per
boundary — whereas re-localization across a frame has no structural bound, so it
takes a budget, whose exhaustion degrades exactly as that section's would.

**Both constants are deployment, not implementation.** `localization_tol` and
`event_budget` are `Simulation` keywords standing beside `h`, `n` and the
algorithm ([§12.1](#121-three-strata), [Appendix B](#appendix-b-api-synopsis-the-entry-points)), validated with their siblings — a positive
tolerance, an integer budget ≥ 1 — and collected into `DeploymentInvalid`
([Appendix C](#appendix-c-the-diagnostic-kind-set)). They are grid-independent, so neither enters [§8.5](#85-multi-rate-tick-scheduling)'s
harmonic-grid check. And being trajectory-determining they are **recorded**:
they ride the trace header's deployment block and join the set replay compares
up front, exactly as `h` and the algorithm do ([§9.5](#95-inbound-the-input-trace), [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)). The
replays-identically promise above is empty otherwise — a run that does not
record what its localizer was told to do cannot be re-driven through the same
localization outcomes.

**Projection's reach is the boundary, not the probe.** Guard probes evaluate the
raw interpolated state — the same rule RK-stage RHS evaluations already live
under, being equally off-manifold, so sweeps must tolerate near-manifold states
and already do. Authority rests with the `t*` boundary: projection runs there,
and the [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) iteration's edge checks read the *projected* state. If projection
moves the state back across a guard, the event does not fire and the run has
published one extra boundary — harmless, and deterministic and pace-independent
like any other localization outcome (row 80). Per-probe projection is rejected:
it puts projection cost inside the root-finder's loop to correct an interpolant
drift far below the localization tolerance.

**`t*` is a boundary, not a frame.** A *frame* is a grid step
`[tₙ, tₙ₊₁]` — the unit of scheduling: input drain at frame top ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain)), pacer
deadlines ([§8.7](#87-real-time-pacing)), tick eligibility ([§8.5](#85-multi-rate-tick-scheduling)). A *boundary* is a published
consistency point — where the [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) macro-sequence completes and a snapshot
goes out. Every grid point is a boundary; `t*` and boundary zero ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)) are
boundaries that are not frame tops. At `t*` the full [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) event phase runs —
[sweep → guards → handlers] iterated to quiescence, **once-per-event
accounting scoped to this boundary** (fresh again at `tₙ₊₁`, and at a second
`t*` on the remainder) — and the settled state is **published**: snapshot,
[§10.3](#103-the-next-snapshot-wait) boundary-counter increment, `stop_on` check ([§13.5](#135-termination-is-a-state-not-an-exception)) — a crash localized
at `t*` ends the run from that snapshot. What does *not* happen at `t*`:
ticks are never due (`t*` is off the harmonic grid by construction; discrete
cells ZOH-hold through the sweep), and staged inputs are not drained — input
timing is a frame fact, and replay determinism must not depend on
localization arithmetic. The publication is not separately paced: the pacer
paces frame deadlines, and a `t*` snapshot publishes when computed, mid-frame
(wall-side placement below pacer resolution; the [§8.7](#87-real-time-pacing) invariant concerns
trajectories, which are identical). Replay pointers and error messages index
boundaries by a monotonic counter with recorded `t`; the trace stays
frame-indexed — `t*` boundaries consume no inputs.

**Endpoint policy and grid integrity.** The root-finder returns the
**holding endpoint of its final bracket** — the smallest probed point where
the predicate holds. Consequences: **`t* = tₙ` is structurally impossible**,
not clamped away — localization only triggers when the predicate was
not-holding at `tₙ`'s quiescence ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) priors), and the interpolant
reproduces the endpoint exactly (`x̂(0) = xₙ`, probe sweeps deterministic), so
the bracket's left end stays strictly not-holding and the returned point is
strictly later than the published, immutable `tₙ` (worst rounding:
`nextfloat(tₙ)`). And the guard observably *holds* at `t*`: handlers fire in
states where their own predicate holds, and the post-fire prior records
an actual observation rather than an assumption. **`t* = tₙ₊₁` exactly is
legitimate** (a crossing at the grid point: `σ(tₙ₊₁) = 0` both triggers
detection and is the root) and **degenerates to the grid boundary**: the
localization result is discarded and the event fires inside `tₙ₊₁`'s ordinary
iteration — bitwise the boundary-detected outcome, one boundary, one snapshot, no
zero-length remainder. A near-degenerate `t*` leaves a tiny remainder step,
numerically harmless (increments are `h′`-scaled); the real hazard is
bookkeeping, killed by rule: **grid times are indexed, never accumulated** —
`tₖ = t₀ + k·h` computed from the frame index (tick gating is already
counter-modulo, [§8.5](#85-multi-rate-tick-scheduling)), and the remainder step *targets the grid point*, with
`h′` derived at use. `t*` is a float inside a frame, never an anchor anything
else is computed from.

### 8.5 Multi-rate tick scheduling

**Harmonic grid.** Every discrete component's period is an integer multiple of a base
tick period `Δt_base`, itself an integer multiple of the continuous step
($\Delta t_{\mathrm{base}} = n \cdot h$, $n \ge 1$). Ticks therefore land on step boundaries — the only place
anything discrete ever happens. Rejected: arbitrary periods via a time-ordered tick
queue, which forces variable-length steps and irregular real-time frames for a
generality nothing demonstrated wants.

**Discrete stages are gated to tick instants.** A discrete component's `h_x`/`h_xu`
run only at its own ticks; its cells hold in between (ZOH, stated in sweep terms). The
alternative — re-running its stages at every boundary — would let outputs drift
between ticks as fresh continuous inputs flow in, silently un-sampling a sampled-data
controller. Delivering that hold takes **two statically distinct sweep variants,
compiled from one entry list** — discreteness is a build-time fact, so the split is
static, not a runtime test (row 147):

- The **interior sweep** walks *continuous entries only*. RK stage evaluations
  ([§8.3](#83-signal-table-consistency-is-a-boundary-property)) and localization guard probes ([§8.4](#84-localization-mechanics)) run this variant, so the ZOH holds
  mid-step **by construction**: discrete entries are not gated out at runtime, they
  are absent from the walk at compile time, and the hot path carries no gating test
  at all. Counter-modulo gating alone cannot deliver it — at divisor 1 (the common
  $n = K = 1$ configuration) the test admits every discrete entry at every
  evaluation, so a discrete `h_xu` would re-run against interpolated `u` at each interior
  stage and the `f` reading its cell would integrate a continuously re-sampled
  controller: exactly the un-sampling this rule forbids, and invisible — such a run
  is as type-stable, allocation-free and replayable as the correct one.
- The **boundary sweep** walks the full list, with discrete entries gated by counter
  modulo against the boundary's tick index. It is the variant the [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) macro-sequence
  runs, and it is not one fixed list: different boundaries run different subsets of
  the schedule.

The split applies to **both sweep blocks**: the discrete tier's `h_x` entries are absent
from the interior stage-1 walk exactly as its `h_xu` entries are absent from the interior
stage-2 walk. The arity distinction that carries it into the phase-body surface —
zero-arg interior, tick-indexed boundary — is [§12.7](#127-the-compiled-executor)'s.

**The due set is a property of the boundary,** not of the sweep call: computed once
for the boundary and reused by every re-sweep of its quiescence iteration ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)), a
due component being at its tick instant for the whole boundary rather than for one
round of it. At a frame top it is the counter-modulo image of the frame index. At a
`t*` boundary it is **empty** — the tick counter has not advanced there, no component
is at a tick instant, and a modulo test against the unadvanced index would wrongly
re-admit the previous frame's due set. At boundary zero it is **everything**: `t₀` is
a grid point of every divisor and no earlier tick exists for a ZOH to hold ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)).

**Simultaneous ticks are already well-defined** by settled machinery: all due
components run their output stages in topological order within the sweep, then all due
`g` updates run after it, in any order — each `g` reads the table and writes only its
own `x` store. The FCS cascade's intra-tick ordering is a sweep property, not an
update-order property.

**Assemblies: virtual for execution, rate scopes for declaration.** There are no
atomic assemblies, and no opt-in variant. Execution atomicity coarsens the schedulable
unit, which is exactly how artificial algebraic loops are manufactured — [§5.4](#54-artificial-loops-and-the-escape-hatch) at
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
(parameter-threading ceremony, cf. [§4.4](#44-function-valued-signals-environment-access)). At build, all scoping compiles away to **one
absolute tick divisor per discrete component**; the boundary sweep gates entries by
counter modulo against the boundary's tick index (above), the interior sweep having
no discrete entries to gate. Recorded limitations: a child cannot tick faster than
its scope ($K \ge 1$) —
soft, since assembly multipliers are declaration sugar and factors can be refactored
onto siblings; and no phase offsets in the first cut (no demonstrated use).

**`Δt` has a single source of truth: the compiled schedule.** Each discrete
component's effective period arrives read-only as the `Δt` field of every
discrete-tier bundle ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)) — `h_x`, `h_xu` and `g` alike, and absent from
continuous bundles, so touching it on the wrong tier is a missing-field error
rather than a rule. It must be readable in the *stages*, not just `g`: per
[§15.2](#152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade), the discretized laws that actually consume `Δt` — a PID's
backward-difference coefficients, a LeadLag's Tustin transform — run in
`h_xu`; `g` is a copy. (A virtual property on the component handle,
`comp.Δt`, after FlightCore's `mdl.Δt`, is the natural-looking alternative,
and it is impossible here, not merely inconvenient: `comp`
is the author's own immutable struct, and two fields holding the same
immutable value are `===`-identical — the [§11.6](#116-paths-wiring-and-faces) argument — while carrying
different `rates` keys, so the period is a fact about the *schedule
position*, with nothing on the instance to hang a property on. The value must
arrive through the call; the bundle field is where.) Author rule: **never
store `Δt`, or any `Δt`-derived coefficient, as a component parameter** —
recomputing derived coefficients per tick is a few arithmetic ops, and a
cached copy is a second thing for gain-scheduling machinery to chase.
Relative declaration structurally enforces the rule for the period
itself: under scoped multipliers a component author *cannot* know their
absolute rate — it does not exist until composition.

### 8.6 Event iteration at boundaries: to quiescence, once per event

Resolves the question deferred in [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries). At each boundary, the event phase **iterates**:
rounds of *(re-run the boundary sweep → evaluate all guards → fire the eligible
events, at most one per component, each `handler → project`)* until a round fires
nothing, under the rule that **each declared event fires at most once per boundary**.

**Why iterate.** Under a single pass, a cascade of N logically-simultaneous
transitions (supervisor FSM → subordinate FSM → …) completes in N steps: latency N·h,
with h an execution parameter. That is model semantics depending on the integrator's
step size — the same footgun class [§2.2](#22-exclusions-deliberate) cited when killing `f_step!` — and [§3.1](#31-continuous-component-the-hybrid-primitive)'s
blessing of externalized FSM components makes cross-component cascades the expected
idiom, not a corner case. Orthodoxy concurs: hybrid automata take sequences of
instantaneous transitions at one time point; Modelica iterates events to quiescence;
Stateflow runs charts to completion within a tick. (Boundary-detection timing remains
h-dependent — that is the *resolution* at which a physical crossing is noticed; the
cascade delay would have been structure the framework inserts between transitions the
model declares simultaneous.)

**Why a full re-sweep per round.** A transition reaches the signal table *only*
through a sweep: a handler writes its component's state stores and nothing else, so
neither the transitioning component's own ports nor the downstream stage-2 chains
reading them have moved. A round therefore re-runs the whole gated schedule. Sweeps
are microseconds and rounds beyond the first require an actual cascade, so the cost
is noise.

**Within-round visibility: one writer, one epoch.** The signal table has a
single writer: **sweeps**. A handler writes nothing to it — it returns
transitions, the framework latches them into the component's state stores,
and `project` normalizes them; auto-publication is a sweep act like any
other stage-1 write ([§12.5](#125-the-always-on-conformance-check)). Nothing moves the table mid-round.
Hence the epoch rule, which is the whole of this section's content:
**a handler executes against exactly the world its guard fired on** — own
`y`, foreign `u`, own `x`/`m` are all the firing round's sweep, so
`y = h(x)` holds at every handler entry, with no bundle straddling two
epochs. Serialization is what delivers it: a component's state stores are
written only by its own handlers, and it fires **at most one event per
round**, so no same-round writer precedes any handler's entry.

A component's other eligible events are not lost, only *blocked*: each is
re-decided in the next round against the post-transition sweep. Declaration
order is therefore a **priority with re-decision** rather than a
simultaneity — an event whose premise the earlier transition falsified
simply does not fire, which under a within-round sequence it would have done
on the stale premise. Across components, handler order within a round is
**semantically unobservable** for the stronger reason that there is no
delivering mechanism at all: nothing writes the table mid-round, so there is
nothing for order to observe. The execution order — executor component
order, declaration order within a component — is fixed only to keep the
[§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor) cursor and the diagnostics stream deterministic, never
something a trajectory depends on. The natural single-pass executor, building
each handler's bundle at dispatch from the live table, is exactly correct:
no pre-materialization, no staging pass, no carrier — no shadow table, no
allocation, trivially.

The trade, recorded openly: a handler cannot opt into seeing a same-round
foreign transition — same-instant sequential coupling across components is a
cascade, one round per link, deterministic; coupling tighter than that
belongs inside one component, where declaration order gives exact sequencing
across rounds (the synchronous-languages position: a micro-step sees the
pre-state, effects appear at the next micro-step). The cost of serializing
same-component firings is one extra intra-boundary sweep per event so
serialized — microseconds, on the rare boundary that fires at all. Rejected
shapes: the per-event re-decode with a frozen round-start `u` obtained by
pre-materializing the fired components' gathers (the prior design of record,
rows 16/100/152 — two freshness classes in one bundle, whose delivery needed
a two-pass staging structure, buying only same-round multi-handler own-`y`
freshness that serialization provides free); live-table reads under the
canonical execution order (deterministic, but model semantics would then
depend on the executor's schedule order — a rewiring that permutes it
silently changes trajectories, the same structure-the-framework-inserts
class this section already rejects); a table copy per firing round
(identical semantics, paying an allocation for nothing); handlers stripped
of their own `y` (makes the incoherence unwritable, but keeps a two-epoch
bundle and the fire-on-falsified-premise hole).

**Why once-per-event rather than a round cap.** Termination becomes structural —
rounds are bounded by the number of declared events, with no arbitrary K knob — and a
livelock (two FSMs toggling each other) resolves deterministically into "each fired
once, quiescence, re-fire next boundary," i.e. degrades to boundary granularity instead
of burning a budget and erroring. The cost: an event legitimately re-enabled within
the same boundary waits one step — accepted at the same granularity boundary detection accepts
physical re-crossings, and flagged when it occurs (a runtime warning, `EventDeferred`,
[§13.2](#132-diagnostics-structured-values-one-carrier-exception); the re-arm edge that triggers it is defined below).

**Guard priors and firing semantics.** "Newly-fired" means precisely this:
an event fires at a boundary iff its predicate ([§2.1](#21-events-two-detection-policies): the `Bool` form true, or
`σ ≥ 0`) is observed holding in some iteration round, its **prior** — the
previous boundary's quiescent sample — was not-holding, and it has not yet
fired this boundary. The per-event **loop state** that decides this is three
registers, and all three are named normatively: the **prior**; a **`fired`
flag**, the once-per-event rule's register, set when the event fires and
cleared at the boundary's end; and a within-boundary **`re-arm` flag**, set
when a round observes the event already fired *and* its guard not-holding.
The re-arm flag is what makes deferral observable. `EventDeferred`
([§13.2](#132-diagnostics-structured-values-one-carrier-exception), [Appendix C](#appendix-c-the-diagnostic-kind-set)) emits when a round observes `fired ∧ re-armed ∧
holding` — a genuine intra-boundary falling-then-rising edge, an event
legitimately re-enabled within the boundary and held back by the
once-per-event rule — at most once per event per boundary. The blessed
sticky case warns nothing: an event that fires and *keeps* holding, with no
not-holding sample in between, never sets the re-arm flag, so it never emits.

The prior is updated at each boundary's quiescence from the final
post-iteration samples, **with one exception: a deferred event's prior is
recorded not-holding** — the intra-boundary re-arm edge survives the
boundary, so the deferred event genuinely fires at the next one and
"waits one step" is true by construction. (Under the bare rule the quiescent
holding sample would swallow the edge: prior holding, so no new firing, and
the deferral would never resolve.) All three registers are detection
bookkeeping, not model memory — correctly *not* in any state store: not captured, not
traced, reconstructed deterministically; the cost is two `Bool`s per event
beyond the prior. A predicate that holds and *keeps* holding fires once, at
the boundary where it first held — sticky flags do not re-fire every
boundary. **Boundary zero establishes every prior
as not-holding**, so a predicate already holding in the authored state
fires at `t₀` — [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)'s behavior, derived rather than asserted —
and a warm restart (`init!` re-runs boundary zero, [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)) resets all three
registers from scratch: predicates holding in the newly applied state fire
again at the new `t₀`.

**Ticks stay outside the iteration, after quiescence.** The two possible couplings
resolve asymmetrically:

- *Events → ticks: handled by machinery already in place.* Due discrete components'
  output stages (`h_x`/`h_xu`) are gated *into* the boundary sweep against a due set
  fixed for the whole iteration ([§8.5](#85-multi-rate-tick-scheduling)), so every
  iteration round refreshes them for free — same `x` (their `g` has not run), post-transition inputs. At
  quiescence, their published outputs reflect the settled boundary instant, which is
  exactly what "sampling at t" should mean for a logically-instantaneous cascade.
  Earlier rounds' tentative values are internal scratch, like RK stage evaluations;
  [§8.3](#83-signal-table-consistency-is-a-boundary-property) extends naturally: external readers observe the table only after the boundary
  sequence *completes*.
- *Ticks → events: structurally impossible.* A tick's output stages contribute nothing
  guards have not already seen (they run inside the sweep, from current `x`); its `g`
  update writes `x⁺` after the sweep, and `x⁺` is first decoded at the owner's *next*
  tick, so it is invisible to every reader within the boundary — the standard
  one-sample `z⁻¹` delay of sampled-data control, here enforced by
  construction. Nothing that happens after quiescence can flip a guard, so no combined event/tick fixed point exists to
  iterate.

The boundary macro-sequence, final form (boundary zero — initialization — is
the same sequence with an empty integrate, [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)):

> integrate → project → **[sweep → guards → handlers]** iterated to quiescence
> (once per event) → all due `g` updates → logging / I/O staging.

The mixed case — a continuous component's handler and its discrete observers' ticks
landing on one boundary (engine `starting → running` under a 50 Hz FCS) — is decided
by the sequence: the transition fires in the iteration segment, the re-sweep re-runs
the FCS's stages against `running`-mode ports, and its `g` then updates from
post-transition values.

### 8.7 Real-time pacing

**The invariant: pacing is outside the semantics.** The pacer inserts waits between
completed frames and never reorders, skips or alters the boundary sequence. A
paced and an unpaced run with identical input traces produce bit-identical
trajectories — deterministic replay ([§2.2](#22-exclusions-deliberate)) extends over pace. Interactive runs differ
only because their *inputs* differ. Detection policy is inside the semantics:
Event localization runs identically paced or unpaced, its sweep cost
absorbed as debt like any other expensive frame ([§2.1](#21-events-two-detection-policies)).

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
(GC, scheduler) are invisible. Debt beyond a threshold — **five frames' worth of
budget, `5·h/p`** — is forgiven by re-anchor plus warning, so long stalls
(debugger, laptop sleep) do not trigger catch-up bursts. Five: comfortably above
the ms-scale hiccups debt exists to absorb silently, and far below the
seconds-to-minutes stalls forgiveness exists for, so neither case lands near the
threshold. Rejected: relative deadlines (next = last completion +
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
purpose). **Its default is 2 ms**, the value the measurements above imply: it
covers libuv's millisecond timer granularity and `sleep`'s ≈1.4 ms median
overshoot, while anything larger merely spends more core in the spin phase. It
spans the whole design space:

- **`margin = 0` — pure sleep:** cheapest CPU; bursty frame spacing, but the
  absolute schedule still delivers the exact *average* rate through debt repayment.
  The spin phase buys regularity, never rate correctness.
- **`margin` = 2 ms — hybrid (the default):** sleeps ~90% of a 20 ms budget, lands
  within µs of the deadline at a few percent of one core.
- **`margin = ∞` — pure busy-wait:** FlightCore's behavior, maximum frame
  regularity at one pinned core; the "best attempt at real time" mode is the knob's
  endpoint, not a separate mechanism.

When the frame budget is at or below the margin (e.g. `h = 0.01` at `p = 5` → 2 ms),
the hybrid degenerates to pure spin per frame by construction. Rare wake-ups past the
deadline are overruns, absorbed as debt. Which primitive the coarse phase uses —
task-yielding `sleep` vs. thread-blocking `Libc.systemsleep` — is settled in [§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget):
the coarse phase uses task-yielding `sleep`, with `margin` absorbing its overshoot.

**Diagnostics.** Overrun count, current and peak debt, forgiven-debt events, wait
statistics — published as framework status for GUI and logs (today's `SimControl`
fields are the precedent).

**Forward pointers.** The wait interval is the natural staging slot for externally
injected inputs, applied at the next boundary; the staging rules — and the
concurrency model generally, which [§8.3](#83-signal-table-consistency-is-a-boundary-property) constrains but does not decide — belong to
[§9](#9-runtime-periphery-the-data-plane), [§10](#10-runtime-periphery-lifecycle-and-orchestration).

---

## 9. Runtime periphery: the data plane

GUI, input devices, network I/O and logging, and how data crosses between them
and the [§8](#8-time-and-execution) loop: the architecture that replaces the shared mutable model
([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)), outbound snapshot publication ([§9.2](#92-outbound-snapshot-publication)), the inbound path — root
input slots, claims and the frozen roster ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)), per-device staging and the
drain ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain)), the input trace ([§9.5](#95-inbound-the-input-trace)) — the device authoring contract
([§9.6](#96-devices-one-authoring-contract-no-taxonomy)), the GUI write path ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)), and the third cross-task channel —
runtime diagnostics and liveness ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)). The machinery that drives the loop
itself follows in [§10](#10-runtime-periphery-lifecycle-and-orchestration).

### 9.1 No shared mutable model: staged writes, snapshot reads

FlightCore's periphery is one big lock: `SimControl` and the live `Model`, guarded by
`io_lock`, with one task per attached interface reading or mutating the model under it
(sim.jl). The lock does enforce [§8.3](#83-signal-table-consistency-is-a-boundary-property)'s boundary-visibility rule — it is only ever
free between steps — but its costs are structural:

- **The loop's frame budget is hostage to its readers.** A slow GUI frame or a stalled
  `extract_output` holds the lock and the sim cannot step; blocking time is
  indistinguishable from overrun in any accounting.
- **Input timing is scheduler-determined and unrecorded.** Writes land between
  whichever boundaries the OS interleaving produced; there is no defined input trace,
  so [§8.7](#87-real-time-pacing)'s bit-identical replay is unachievable *in principle* for interactive runs.
- **It protects an idiom that no longer exists.** `assign_input!` and GUI widgets poke
  the live model; under the immutable signal table there is nothing to poke — the
  periphery needs a defined write path regardless.

The replacement has five planes. (Vocabulary, anchored here: a **frame** is one
iteration of the loop — drain, integrate, boundary sequence, publication — the
unit `step!` counts, the trace's ordinal keys, and "per frame" means throughout
this document. Distinct from the kinematic *reference frames* of the aircraft
domain, which always appear compounded: the b frame, the ECEF frame.)

1. **Staging (inbound):** devices submit pending input writes at any wall-clock
   moment, never touching live slots ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain)).
2. **The drain:** exactly one point, at the top of each frame — never at a `t*`
   boundary ([§8.4](#84-localization-mechanics)) — where the
   loop takes the staged batches and applies them to the root input slots. Between
   drains the loop owns its data exclusively — no lock is held during stepping, ever.
3. **Publication (outbound):** at the end of each boundary sequence the loop publishes
   an immutable snapshot; readers observe it without coordinating with the loop
   ([§9.2](#92-outbound-snapshot-publication)).
4. **Control:** pause/pace/stop on a separate few-word atomic surface ([§10.1](#101-control-plane)).
5. **Task topology:** one loop, one task per rostered device except the
   calling-task device, all run-scoped: `run!` spawns one task per other
   roster entry after device `init!`, and [§10.4](#104-shutdown-protocol) joins them all at every stop ([§10.6](#106-run-lifecycle-and-partial-advance)).
   `attach!` never spawns — it registers, in a stopped-sim state only
   ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)), and the task appears at the next `run!`. **The calling-task
   device is pinned; the loop is the movable piece.** Calling-task
   affinity is a device trait (`needs_calling_task`, default `false`,
   [§9.6](#96-devices-one-authoring-contract-no-taxonomy)) with at most one holder per roster ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s admission checks);
   the shipped GUI declares it — CImGui ties rendering to the calling
   (main) task. The topology is derived from the frozen roster alone —
   as device `init!` leaves it, a failed calling-task holder returning the
   loop to the calling task ([§10.4](#104-shutdown-protocol)) — and
   never from `run!`'s keywords: with a calling-task device rostered the
   loop moves to a spawned task for the duration of the run and the
   calling task runs that device's loop body — inline, inside the same
   [§9.6](#96-devices-one-authoring-contract-no-taxonomy) wrapper as any spawned device's — otherwise the loop runs on the
   calling task — the unattended register, what [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)'s synchronous rethrow
   presupposes, and what lets parallel unattended sweeps thread `run!` inline
   with no nested task fan-out (one immutable `Build` shared across the
   workers, [§12.2](#122-the-build-artifact), each `Simulation` owning its own buffers; pre-materializing
   the sweep's activations — `build(world; activations = …)`, [§12.4](#124-activations-executable-sets-laziness-caching) — leaves no
   worker synchronizing on anything). Either way `run!` blocks its caller until
   the run ends; what varies is what the calling task spends the run
   doing. Spawn-inside-`run!` *is* the start gate — a task exists only
   once the run it serves exists — and any first-boundary synchronization a
   device needs is [§10.3](#103-the-next-snapshot-wait)'s counter-plus-condition predicate wait, never an
   `Event` latch: FlightCore's `io_start` gate is the once-per-run version of
   exactly the race [§10.3](#103-the-next-snapshot-wait) rejects, and inheriting it would re-import that
   race for [§10.6](#106-run-lifecycle-and-partial-advance)'s re-run cycle.

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

Consequence, recorded because it collapses an API axis: interactive and unattended
simulation stop being different execution modes. An unattended run is the same loop with
empty staging and no snapshot readers; a replayed interactive session is the same
loop with staging fed from a recording ([§9.5](#95-inbound-the-input-trace), [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)).

### 9.2 Outbound: snapshot publication

The loop builds each snapshot — boundary-consistent signal table, `t`, framework
status — in private memory, then publishes it with a single
release-store to an `@atomic latest` reference; readers acquire-load and then work
with an immutable, coherent world for as long as they like. `latest(sim)`
hands the same value to the calling task — [§10.6](#106-run-lifecycle-and-partial-advance)'s inspection register.
Wait-free in both
directions: a wedged reader cannot delay publication by a nanosecond; the loop cannot
tear a reader's view. Publication happens only after the boundary sequence completes
([§8.3](#83-signal-table-consistency-is-a-boundary-property) as extended by [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)).

**Binding rule: nothing reachable from a published snapshot is ever written again.**
The table's immutable values ([§4.1](#41-immutable-value-semantics), [§7](#7-state-and-data-representation)) make the compiler enforce most of it; the
rule is what the soundness of lock-free reading rests on.

**The framework status is a concrete frozen value, not a window onto live
bookkeeping** — [§8.7](#87-real-time-pacing)'s pacer diagnostics, plus the per-writer diagnostic
batches, suppressed and cumulative counters and liveness timestamps the loop
takes at frame top ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)). The binding rule is what forces that shape: a
status referencing an accumulator its writers are still filling would be a
snapshot whose contents change after publication.

The captured table is the whole table — declared ports and auto-published fields,
every one of them public ([§11.3](#113-visibility-the-contract-is-the-interface)), so no presentation layer has anything to
filter. Private intermediates are not in it, never having been cells at all
([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)); the inspection path for one is **promotion to a declared output** — a line
in `output_types` and the value appears in the snapshot, the log, the GUI and the
wiring alike, its visibility an authored fact like every other. **It also includes
the root slots** ([§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)): slots are source cells of the
table, not state stores, so they ride along — and this is load-bearing, not
incidental: the [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract) peek's else-snapshot fallback is what an idle live widget
displays, and read-only mirrors of claimed slots (the axis sliders under joystick
claim) show the applied slot value from the snapshot. Slot values in the log are
derived data (recomputable from the trace), which is consistent — snapshots are
derived wholesale. The snapshot
deliberately does **not** carry the state stores (`x`, `m`): the state
trajectory is *derived* data — recomputable from the trace header plus the batches
([§9.5](#95-inbound-the-input-trace)) by bit-identical replay — and per-boundary capture would systematically
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
framework side of the [§7.5](#75-allocation-policy-a-scoped-invariant) scope (which already carved out logging); logged snapshots
are not garbage at all, unlogged ones die young. Rejected: preallocated snapshot
buffers (double/triple ring) — reuse reintroduces exactly the reader-liveness proof
the GC provides for free, to save an allocation profiling has not indicted.

**Retention: the trace's kill switch, plus decimation.** The log takes the same
plain on/off switch the trace has ([§9.5](#95-inbound-the-input-trace)), and additionally a keep-every-kth
policy (`log_every`, [Appendix B](#appendix-b-api-synopsis-the-entry-points)). What makes decimation admissible here and not
there is the derived/primary split (row 38): the log is recomputable from the
trace by replay, so a thinned log costs resolution in a *view*, never a record —
whereas thinning the trace would destroy the only primary account of a session
(row 29). Decimation is a retention policy only: every boundary still runs, is
still published to live readers, and still enters the trace.

**The bound: `log_max`.** Decimation slows the log's growth; it does not stop
it, and the default configuration — `log = true, log_every = 1` alongside
`t_end = Inf`, [Appendix B](#appendix-b-api-synopsis-the-entry-points)'s honest interactive default — grows for as long as
the session lasts, which at C172X scale and 50 Hz is gigabytes per hour and
ends in an out-of-memory nobody was warned about. So the log takes a
**retention bound** beside its switch and its stride: `log_max`, the maximum
number of retained snapshot references, default **65536** (2¹⁶), with `Inf`
the explicit opt-out. **A count, not a memory budget** — snapshots are
immutable object graphs with internal sharing ([§4.1](#41-immutable-value-semantics)), and [§4.4](#44-function-valued-signals-environment-access) field handles
ride as references to build-time-frozen data ([§7.5](#75-allocation-policy-a-scoped-invariant)), so byte accounting over
them is fuzzy and platform-dependent, while a count is exact and converts to
memory through one number the user can measure once (`Base.summarysize` of a
single snapshot). The default is **finite unconditionally**, not a modal rule
keyed on `t_end`: a finite run shorter than the bound never notices the bound,
and one number is easier to hold than two regimes. At 50 Hz and full density,
2¹⁶ boundaries is about 22 minutes before anything is dropped at all.

**When the log fills, the retention stride doubles: coverage stays global.**
The obvious policy — a rolling window, keeping the recent past and forgetting
the start of the session — optimizes for exactly the thing replay already
serves, and pays for it with the thing replay is expensive for. Instead the log
**re-decimates progressively**: after *k* generations the effective stride is
`log_every · 2^k`, so the whole run stays plottable and what coarsens is
density, never extent. That is row 38's division of labor carried through: the
log's chief consumer is the post-run plot of a session *as a whole* (nobody
plots hours at 50 Hz), while full density over any *segment* of interest is
what replay from the trace recovers ([§9.5](#95-inbound-the-input-trace), [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)). **Normative are the
guarantees, not the mechanism**: the bound is respected continuously (the
retained count never exceeds `log_max`), coverage is global at the effective
stride, and the endpoints below are kept.

Mechanism sketch, non-binding: the thinning is **amortized**. When the log
fills, the stride doubles immediately, and each subsequent retained append also
releases one predecessor of the previous generation — a cursor over the odd
indices, O(1) amortized, with physical compaction once per generation — so a
generation's thinning completes exactly when its refill does. Amortizing rather
than halving in one shot is a responsiveness choice: a one-shot halving would
make ~`log_max/2` *old-generation* snapshots unreachable at a stroke, a
major-GC burst worth a pacer-absorbed, `DebtReanchor`-visible hiccup ([§8.7](#87-real-time-pacing)) —
exponentially rare, since the wall-clock gap between generations doubles, but
needless. The amortized form drops exactly one old snapshot per retained
append: the same steady trickle a rolling window would produce, so keeping
coverage global costs nothing extra in GC pressure. The loop's own work is
pointer bookkeeping, microseconds either way and on the framework side of
[§7.5](#75-allocation-policy-a-scoped-invariant)'s scope; publication stays wait-free and readers never block — a reader
holding a released snapshot simply keeps it alive.

**The endpoints are retained unconditionally.** The boundary-zero snapshot
([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)) and the terminal snapshot ([§10.4](#104-shutdown-protocol), [§13.5](#135-termination-is-a-state-not-an-exception)) survive any `log_every`
and any `log_max`, and do not count against the bound — two extra references.
The terminal snapshot's status carries the run's final cumulative diagnostic
counters ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)), so a run's two endpoints and its complete diagnostic account
always outlive whatever retention did to the middle.

Two compositions, worth stating once. The `totals` monotonicity across logged
snapshots ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)) is untouched: re-decimation, like decimation, loses *which*
boundary within a stretch an occurrence fell on, never *how many*. And
`log_max` is a **view policy, not a trajectory-determining one** — like `log`
and `log_every` it stays out of the trace header's deployment block, and replay
neither records nor compares it ([§9.5](#95-inbound-the-input-trace), [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)). Sizing follows: [§7.5](#75-allocation-policy-a-scoped-invariant)'s
`sizehint!` for the expected duration is now naturally capped by `log_max`,
which is also what defines the hint when `t_end = Inf`.

**Output-device bindings are snapshot bindings.** An output device (telemetry, the XPlane visualizer, disk
streaming) consumes snapshots via [§10.3](#103-the-next-snapshot-wait) and addresses what it reads with
[§14.4](#144-two-application-registers-over-one-plan)'s selectors — any cell, the diagnostic register admitting deep
paths — resolved at attach against the `Build` with
did-you-mean and compiled to one gather (the output half of [§9.6](#96-devices-one-authoring-contract-no-taxonomy)'s binding
interface), so `map_output` receives a labeled
NamedTuple rather than performing its own path lookups ([§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)'s obligation: a
substitution that breaks a binding fails at attach, not with silent garbage
UDP). This is **diagnostic observation** ([§13.5](#135-termination-is-a-state-not-an-exception)):
human-facing, no effect on run semantics — the same register as the log
retaining the full table and the GUI's deep-reading panels; every cell is
reachable, the table being public throughout ([§11.3](#113-visibility-the-contract-is-the-interface)), and an intermediate a
device wants to stream is one promoted to a declared output. **A binding chooses
its register**: a deep path is the
*inspection* register — zero promises, free access, right for looking at
*this* build; an exported output face (spelled `get_face(name)`, [§14.4](#144-two-application-registers-over-one-plan)) is
the *integration* register — named,
curated, meaning-stable under substitution ([§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)'s writer-independent
semantics), right for consumers that outlive the build they were configured
against. Attach validation converts *structural* drift to loud errors in both
registers; only faces protect against *semantic* drift — a substituted
aircraft publishing the same path at the same type with a different meaning, a
CG velocity under a name read as body-origin velocity — and nothing else can:
meaning is not in the schema. Semantically generic consumers (a visualizer
needs pose; every aircraft has one) should therefore bind faces, and aircraft
families should export the conventional surface such consumers need — a
library/migration deliverable ([§16](#16-open-axes)) — with wrapper types making face semantics
structurally checkable (`VelocityData`, its `v_eb_b` defined *at the type* as
body-origin velocity: a bare vector doesn't wire, and wrapping the wrong
quantity is a deliberate lie, not a drift).

### 9.3 Inbound: root input slots, claims and the frozen roster

**The write surface is root input slots** — and a root slot *is* the root
assembly's own input face, declared through `input_connections` ([§11.6](#116-paths-wiring-and-faces)): routed inward to consumers, produced by no
component, fed by the parent's wire at every non-root level — and at the root
there is no parent. (A root slot is usefully read as the output face of the
one producer the build never sees — the periphery and the services: slot
exclusivity below is a producer's one-writer right, and [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)'s totality is
its completeness obligation.) No dedicated vocabulary survives (`add_input!` in the early
sketches is dead). Slots are sources to the build-time scheduler, constants within
a frame, and the *only* thing the periphery may write (the GUI reaches them
through [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)'s resolution; control commands are not writes, [§10.1](#101-control-plane)); devices,
mappings, the trace and the GUI write path address them by **face name** ([§11.6](#116-paths-wiring-and-faces)):
structural slash paths never cross the periphery's *write* boundary — the write
side speaks the root contract's names only. (The read side chooses per
binding: slash paths in the inspection register, face names in the
integration register and in load-bearing service reads —
[§9.2](#92-outbound-snapshot-publication)/[§13.5](#135-termination-is-a-state-not-an-exception)/[§14.4](#144-two-application-registers-over-one-plan).)

**Slot exclusivity: one writer per slot at any time** ([§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)). A
device claims its slots at attach; claiming an already-claimed slot is an
attach-time error, and detaching releases the claims (a released slot's GUI
widgets are live again from the next run, [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)). Exclusivity replaces any cross-device conflict *policy* —
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
for the device's whole attachment, so [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)'s derived-liveness rule renders its
GUI widget read-only even on faces the peer never writes. Narrow the binding to
narrow the claim — the enumeration *is* the interface.

**Every writer has a write surface, and the periphery enforces it.** A batch
entry reaches a slot **iff the named face is inside the writer's surface**;
anything else is discarded with a runtime warning ([§13.2](#132-diagnostics-structured-values-one-carrier-exception)). Because surfaces
are static per run (the roster freeze below), enforcement runs entirely at
*staging* — the earliest site, on the writer's own task — and the drain
performs no checks at all. **Every device's surface is its claim set**, and a
claim set has two *sources*:

- **Returned** — the binding enumerates the faces: it declares
  `is_input(b) = true`, `claims(b)` ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)) is
  called once at attach, and what it names is staked. Static for the
  attachment, exclusively its own (claims are disjoint by construction), and
  binding-bounded even where no one else is involved — a mapping that has
  drifted onto an unenumerated face is a diagnosable anomaly
  (`OutOfClaimEntry`), never a silent write, claimed or not.
- **Computed** — the binding declares [§9.6](#96-devices-one-authoring-contract-no-taxonomy)'s `is_greedy(b) = true` and the
  framework computes the claim at attach: all root input faces minus the
  union of the rostered claims, the unclaimed complement at that instant.
  This is the shipped GUI's claim ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)) — everything unclaimed, without
  configuration — and it is disjoint from every incumbent claim by
  construction, so exclusivity validates trivially and nothing downstream
  can tell the two sources apart.

One claim mechanism, two claim sources. The source is exhausted at the attach
point: past it, validation, roster-entry storage, [§9.4](#94-inbound-per-device-staging-representation-and-the-drain)'s shape
compilation, the drain, the trace and detach-releases-claims treat a computed
claim exactly as a returned one, which is why the GUI is not an exception but
an ordinary enumerated writer whose enumeration the framework performed.
Opportunistic writing by autonomous devices does not exist: a device that
wants a face enumerates it, and greediness is an explicit declaration, never a
default. Cross-writer races on one slot therefore cannot arise structurally —
every claim is exclusive, whatever its source — which is what keeps drain
order a diagnostic fact (below) and lets a drained GUI value simply stay
([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)).

**One framework-owned remainder: the harness register.** Beside the roster
sits a **task-free entry point** — `stage!(sim, "face" => value, ...)`, which
stages a batch from the calling task itself, the harness/REPL write path
([§10.6](#106-run-lifecycle-and-partial-advance)) — and its always-present cell, drained, traced and
surface-checked exactly as any device's. Its surface is the one thing in the
design that is *derived* rather than claimed: the unclaimed complement, the
faces no rostered device speaks for, recomputed at every stopped-sim roster
change and therefore as fixed within a run as any claim set. A `stage!` write
to a claimed face is rejected at staging (`ClaimedFaceEntry`, naming the
incumbent), and the one seam — a batch staged while stopped whose face a
subsequent `attach!` claims — is renormalized away at the attach itself
(below). The harness cell drains **last**, by convention: with every surface
disjoint the order is unobservable, so the rule exists to make the trace read
the same way every time, not to arbitrate anything.

**Slot initial values are owned by the init/trim services** ([§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)). Input declarations are bare types ([§11.2](#112-the-declaration-inventory)) and carry no defaults, but a
slot unfed by any device must hold a defined value from the first frame (today's
`U()` constructors provide these: `mixture = 0.5`). Export-entry defaults were
rejected: the trim service writes slot values it *solved for* (throttle,
elevator) — not declaration constants. `init!` establishes every slot and the
trace header captures the result; totality is enforced pre-write at every
complete-world application — `init!`, trim setup, trim commit
([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)).

**The roster is frozen per run: attach and detach are stopped-sim
operations.** `attach!`/`detach!` are legal in the `built`, `initialized`
and `stopped` states ([§10.6](#106-run-lifecycle-and-partial-advance)) and an error while `running`
(`ServiceLifecycle`, [Appendix C](#appendix-c-the-diagnostic-kind-set) — the same kind that gates the [§14](#14-stopped-sim-services)
services) — pause included: pause is a control-plane state *inside* a run
([§10.1](#101-control-plane)), and a surface that could move while paused would move mid-run. The roster —
entries, claims, attachment order — is therefore a plain immutable value
the loop reads once at `run!`, and the partition of the root face set into
per-writer surfaces plus the harness remainder is a static,
inspectable fact of the run — printable before it starts, valid until it
ends ([§13.7](#137-tooling-consequences-provenance-and-the-component-library)'s provenance register). No republication machinery exists:
no atomic roster reference, no per-frame acquire-load, no next-frame
attachment granularity, no sequence numbers — attachment order is the
roster's own order. The trace still tags entries with a stable device id,
never a roster index — ids read across runs, where the roster does change.
Attach validation, claim registration and the staging-shape compilation
(below) all run at the attach point, making `attach!`/`detach!`
stopped-sim configuration operations beside `init!` and trim ([§14](#14-stopped-sim-services)): while
a simulation runs, its configuration — build, roster, claims, surfaces — is
immutable, and [§10.5](#105-scripts-and-the-mid-run-mutation-doctrine)'s doctrine extends to its final form — the running
periphery stages writes and issues control commands, *and nothing else
changes*.

**Device identity, ids and roster admission.** Identity is the device
instance: the same object (`===`) may occupy at most one roster entry,
while two instances of the same type (two joysticks) are two devices. The
stable device id the trace, heartbeat and diagnostics speak is assigned
at `attach!` — monotonic per `Simulation`, never reused — and lives
exactly as long as the entry: across runs (roster persistence, [§10.6](#106-run-lifecycle-and-partial-advance)),
until `detach!`. Admission is a three-part check at the attach point, in
order: **identity** — an already-rostered instance is rejected
(`AlreadyAttached`, naming the entry and its binding), because rebinding
has an explicit spelling — `detach!` then `attach!`, both legal at any
stopped-sim point — and either a silent no-op or a silent rebind would
discard a binding the caller handed over; **affinity** — at most one
rostered device may declare `needs_calling_task` ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)'s topology makes
the calling task a single-slot resource; `CallerTaskConflict`, naming
both devices); **claims** — face exclusivity (`ClaimConflict`), which by
running after the identity check always names two *distinct* devices,
never a device colliding with its own earlier attachment.

**Device death does not detach.** A mid-run crash, voluntary exit or unplug
([§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§10.4](#104-shutdown-protocol)) ends the device's *task*: the cell stops filling, the [§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget)
heartbeat shows the death by name, and the roster entry — claims included —
persists to the end of the run. The orphaned claims are the accepted cost of
the freeze: the device's slots hold their last-drained values and no other
writer inherits them; [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)'s read-only widgets render the orphan visibly
("claimed by `T16000M` — task dead"), never mysteriously. Recovery is
between runs — stop, `detach!`, and either `init!` (fresh trajectory) or
`replay!`-to-end then `run!` (continuation from the interrupted boundary,
[§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)) — and the anomaly is exactly that: an anomaly, not a surface event.
One deliberate asymmetry is on record as a **guarded addition**: a pure
reader (a binding declaring `is_output` alone, [§9.6](#96-devices-one-authoring-contract-no-taxonomy) — a visualizer, a telemetry
tap) claims nothing, so attaching one mid-run would move no writer's
surface; a dynamic reader list (touching only [§10.3](#103-the-next-snapshot-wait) wakeups, the heartbeat
and the shutdown join — never the drain) is cleanly severable from the
freeze should the join-a-running-session workflow find a customer. The
[§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget) thread-budget warning runs once per `run!`, against the frozen
population.

### 9.4 Inbound: per-device staging, representation and the drain

**Staging: one atomic cell per attached device, one coalescing policy — CAS
merge, newest wins per face.** Each cell has a single writer — its own
device task — and holds that device's latest pending batch of slot writes.
Staging merges the incoming batch into the pending one: untouched faces
survive, re-staged faces take the newest level — the per-face ZOH. The CAS
can fail only because a drain intercepted the old batch, so the retry is
bounded and the failure case is precisely correct — intercepted writes are
already applied and must not be re-staged. Merge is the *only* policy
because it is always correct: for a **complete** writer (a joystick: full
write-set every poll) every incoming batch covers every face, so merge and
overwrite are provably the same operation — overwrite is a degenerate fast
path, not a second semantics — while for a **sparse** writer (the GUI, a
JSON peer: only what was touched) overwrite is a silent lost-write bug
([§15.3](#153-torture-test-for-the-9-staging-shapes-filter-joystick-and-gui)'s hazard: a pending `flaps` edit clobbered by an unrelated `gear`
message, undrained, undiagnosable). A user-facing overwrite opt-in
(`complete(binding)`, a declared promise of batch totality) was drafted and
dropped: under the fixed-shape representation below the fast path saves a
pending-read and a small tuple rebuild per staging, on the device task —
not worth a declarable promise whose false direction loses writes.

**The staged representation is fixed per attachment, compiled at attach.**
An enumerated writer's claim set and slot types are both known at attach
(`claims(binding)`, [§9.6](#96-devices-one-authoring-contract-no-taxonomy), against the root contract), so the framework
fixes the cell's content type there: a positional tuple over the claim
set, `Union{Nothing, T}` per face (isbits unions — pointer-free), with
`nothing` meaning *not touched this time*, never "reset" — the levels
doctrine is untouched, slots only ever receive the non-`nothing`
positions, and the `Union` never reaches the model. The face-name →
position schema lives in the roster entry. The consequences are each
mechanical: the merge is positional (`incoming[i] === nothing ? pending[i]
: incoming[i]`) — straight-line, union-split; the drain applies each cell
through an attach-compiled **scatter** (position → slot cell, statically
typed, `nothing` skips) — the exact mirror of [§9.2](#92-outbound-snapshot-publication)'s compiled output
gather; and authors never build the shape by hand — `map_input` returns
face ⇒ value pairs for whatever the datum touched, and `stage!` normalizes
through an attach-compiled shim (name → position, convert to the slot's
declared type, fill `nothing`), confining the residual name-shaped
dynamism to one framework-owned conversion on the device task, at the
boundary where wire-shaped data becomes system-shaped data. (Author-built
total tuples were rejected as a padding form — ten explicit `nothing`s to
say "one face touched" — the same disease row 74 and the handler return
law refuse.) A **greedy entry needs no special treatment here**: its claim
was computed at the attach point and is an ordinary claim set by the time
shapes are compiled ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)), so the GUI's cell is compiled exactly as a
joystick's. The **harness cell gets the same treatment**: under the roster
freeze its derived surface — the unclaimed complement — is as static as any
claim set, so it too is compiled to a positional shape, recompiled at each
`attach!`/`detach!` (a stopped-sim point), with the same shim, merge and
scatter; always present, it gets the compilation unasked, and it is the one
cell whose shape the framework derives rather than receives. One
representation, one mechanism; the
name-keyed dynamic path the mutable surface used to force does not exist,
and no face name is ever resolved inside the loop's frame. The
recompilation has one seam: a pending harness batch staged *before* a
stopped-sim `attach!` may hold the old shape, or name a face the new claim
covers — the attach renormalizes it (reshape, discard newly-claimed faces
with `ClaimedFaceEntry`), so the run always starts with cells matching the
run's schemas.

**Diagnostic sites follow the compilation — all of them to staging.**
Face-name validity, surface membership and value convertibility are all
static facts of the run, so every check runs in `stage!`'s normalization,
on the writer's own task: a device's out-of-claim face has no
position in the schema and is rejected (`OutOfClaimEntry` — an earlier,
better-attributed site than the drain; same kind, same payload — the GUI
included, its claim being an ordinary one); a **harness** write to a claimed
face is rejected the same way
(`ClaimedFaceEntry`, naming the incumbent device); and a value that cannot
convert to its slot's declared type is discarded with `EntryTypeMismatch`
([Appendix C](#appendix-c-the-diagnostic-kind-set)) at the same spot. Nothing remains at the drain: with surfaces
frozen for the run, there is no fact only the drain can know, and the
drain is pure application.

**Doctrine: staged values are levels, never deltas** (`press_count = 17`, never
`presses += 1`) — levels are idempotent and survive coalescing; button edges ride as
monotonic counters.

**The drain** (frame top): one `atomicswap(cell, nothing)` per device — an
indivisible take, no lost-write window — every cell applied through its
compiled scatter, **in attachment order**, retained
as a deterministic application order (under slot exclusivity, cross-device writes
to one slot no longer arise; the order matters only for diagnostics). Which
*frame* a write lands in remains wall-clock reality; what the drain guarantees is
that the frame's outcome is a pure function of the drained batches. Because
the roster is a fixed value at `run!`, the drain is fully compilable: the
cells and their scatters form a heterogeneous but *known* tuple the frame
function can specialize on — zero dynamic dispatch at frame top — the same
per-configuration compile trade [§12.7](#127-the-compiled-executor)'s executor already makes, now incurred
only at stopped-sim attach points. (The specialization is an implementation
freedom the freeze creates, not an obligation; iterating a roster array
costs a handful of dispatches per frame and remains acceptable.)

Rejected shapes (both torture-tested in [§15.3](#153-torture-test-for-the-9-staging-shapes-filter-joystick-and-gui)): **per-slot atomic cells** — the
simplest (no merge machinery, and a per-slot layout cannot lose independent writes)
but same-slot conflicts resolve by hardware store order, i.e. sub-frame wall-clock
phase (run-to-run behavioral variance, [§15.3](#153-torture-test-for-the-9-staging-shapes-filter-joystick-and-gui)), peeks are cross-device, the trace
loses provenance, and wide slot types hit Julia's atomic-width lock fallback; **a
shared lock-free batch stack** (CAS-push, swap-drain) — whole-batch atomicity and
the richest trace, but conflict order is still temporal (push order), and pending
memory is unbounded while paused, taxing every peek that must walk the chain.

**Mappings run on the device task**: today's `assign_input!(mdl, mapping, data)`
becomes pure `map_input(data, mapping) → batch`. User-extensible code thereby never
executes inside the loop's frame, and the trace consists of slot-level batches.

**Mappings are binding data, not shaping code** ([§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)). A mapping is
a declarative table — axis/button → slot path, plus per-axis conditioning
parameters (deadzone, expo strength) applied by the shipped `TableBinding`'s
generic `map_input` on the
device task ([§9.6](#96-devices-one-authoring-contract-no-taxonomy) — the shared pure helper, with an owner). The boundary is set by the face contract: **a face's meaning is
writer-independent**, so faces carry *post-conditioning* semantics — a GUI slider
or a script writes the same command a curved stick delivers (running a mouse drag
through a deadzone would be absurd); this GUI-parity test is what places
conditioning upstream. Aircraft-semantic derivation (the C172X `q_ref = q_sf ·
axis` fan-out) must *not* ride along: it is FCS design and lives in-model — in
the avionics, or accepted as a small per-aircraft×device mapping entry (an
aircraft-design fork, [§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)). The trace records post-conditioning levels —
exactly what the model consumed, so replay is exact; raw-stick provenance (re-run
a session through *different* curves) is the known, accepted loss. Edge logic
follows the levels doctrine: devices stage monotonic press counters; accumulators
(trim offsets, flap detents) are model state, not mapping state ([§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)).

### 9.5 Inbound: the input trace

**The input trace** is the sequence of drained, device-tagged batches per frame. It
extends [§8.7](#87-real-time-pacing)'s determinism end-to-end: replaying a recorded interactive session —
staging fed from the recording, no devices or mappings present — reproduces the
trajectory bit-identically.

**One record format: every batch is retained sparse.** At the drain, each
drained cell is scanned and recorded as (position ⇒ value) pairs for its
non-`nothing` entries, against the writer's face-name → position schema in
the header (below) — an O(surface-width) scan and one small allocation per
drained batch. The rule is uniform because the alternative is not
statable any more: a claim's *width* is a fact about one binding, not about
a class of writers — a greedy claim is enumerated and as wide as the root
contract ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)) — so a density dichotomy could only be re-keyed on the
claim source, which is precisely the distinction nothing downstream of
attach is allowed to see. Uniformity is what the consumers get paid in: one
record format at the trace's edge, no per-entry format flag, one decoder in
the what-if register, in disk serialization and in human inspection, and one
inverse conversion in replay — paid once, up front, off the loop
([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)). The conversion site is the drain and not the staging
shim because the drained tuple is the *coalesced* truth — a shim-side
sparse log would need its own merge.

**The costs are recorded rather than argued away.** On the wide writers the
conversion is what keeps the trace honest: a tuple as wide as the unclaimed
surface carrying one edit would otherwise make trace size track surface
width rather than information (at hundreds of faces, render-rate dragging
inflates the trace past the two-orders-below-the-log budget that justifies
trace-on-by-default, row 29). On the dense component it costs **about 2×** —
a position beside every value where the positional tuple carried the value
alone — which changes no order of magnitude and leaves row 29's budget
standing for every writer at once. The allocation is in-class with what
[§7.5](#75-allocation-policy-a-scoped-invariant)'s retention carve-out already admits and smaller per boundary
than the log's snapshot, the carve-out's standing occupant (the one
qualified exception to retains-what-was-already-allocated). And the decision
is **reversible as pure implementation**: the conversion is lossless in both
directions, so verbatim retention could return as a per-entry storage
optimization — under the same record semantics, the same header and the same
replay path — if a marathon-session measurement ever asks for it.

**The trace header captures the full initial state** `(x, m)` **plus the
initial root-slot values** at `init!` — captured **after `apply!` and the slot
writes, before the boundary-zero sequence runs** ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)). Both halves of that
placement are load-bearing: the header holds the *resolved* stores and slots
as values, never the sparse authored overlay (replay must survive edits to
declared defaults — row 38's primary-data doctrine), and never the
post-transition result — boundary zero is re-executed under replay ([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)),
so a post-sequence capture would re-fire authored-condition events on top of
already-latched state. (An unfed `mixture = 0.5` never appears in any batch,
so replay is broken without the slots; the init/trim services own slot
initialization ([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)), and the header capture extends naturally.) The header
also carries **each writer's face-name → position schema** — the run's frozen
surface partition — since positional records are meaningless without it and
replay does not reconstruct claims ([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)). And it carries the run's
**deployment block**: `t₀`, `Δt_base`, `h`, `n`, the algorithm identifier,
`localization_tol`, `event_budget` ([§8.4](#84-localization-mechanics))
and the effective `t_end`/`stop_on` pair, captured at the same instant as the
stores. The trajectory depends on these exactly as it depends on the stores —
[§12.1](#121-three-strata)'s deployment binding sits outside the `Build`, and `t₀` post-dates
even deployment ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)) — so a header without them could not back
[§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)'s bit-identity claim. This block is also the artifact **run
metadata** names ([§13.5](#135-termination-is-a-state-not-an-exception), [Appendix B](#appendix-b-api-synopsis-the-entry-points)): the recorded home of the
effective termination pair. This is the one
full-state capture in a normal run, and the other half of what "given the
initial state and the trace, the log is recomputable" requires. Header plus batches
are the *primary* record; everything else, the state trajectory included, is
derived ([§9.2](#92-outbound-snapshot-publication)).

**Trace recording is on by default** (cleared at `init!`, retrievable after the run,
plain kill switch for memory-constrained marathon sessions). The asymmetry that
decides the default: the trace is *primary* data and the log *derived* — given the
initial state and the trace, the log is recomputable (that is what bit-identical
replay means), while an untraced interactive session is unreproducible, permanently.
The cost supports it: the trace retains one small sparse record per drained
batch (above), at
drain-rate × device-count — tens of MB per hour worst case, two orders of magnitude
below the snapshot log. No sampling, no rolling window (complexity without a
customer).

### 9.6 Devices: one authoring contract, no taxonomy

FlightCore's input/output/GUI trichotomy is lock choreography, not modeling:
`get_data!` may block, so it runs outside the lock; `extract_output` must not block,
because it runs inside; the GUI breaks both rules at once and gets a third interface
(`render!` under lock, `sync = 0` so VSync cannot stall the sim, a manual framerate
sleep). With no lock, the protocol the taxonomy encoded has no referent.

#### Every attached device receives the same handle

**Every attached device receives the same handle**, carrying the two primitive
capabilities — read (latest snapshot; optionally wait-for-next-boundary, [§10.3](#103-the-next-snapshot-wait)) and
stage —
plus control access (observe running, request shutdown). **`should_abort` is an
`attach!` keyword**, defaulting to `false`: per-attachment, never a device
property — the same joystick is advisory in one deployment and load-bearing in
another — so with it clear a device's departure (loop body returning, crash, or
a failed `init!`) is reported and the run continues without it, and with it set
that departure also requests a sim stop ([§10.4](#104-shutdown-protocol)). The shipped GUI attaches with
`should_abort = true` — closing the window is the interactive session's natural
end — and `gui = true`'s run-scoped attachment states that value ([§10.6](#106-run-lifecycle-and-partial-advance),
[Appendix B](#appendix-b-api-synopsis-the-entry-points)).
Input-only and output-only devices are degenerate uses, not
framework classes; a bidirectional network peer is *one* device with one socket and one
lifecycle, not two framework devices sharing state. The GUI is an ordinary device —
the paradigm one, using every capability — with exactly two genuine peculiarities,
neither taxonomic: main-thread affinity (a launch concern) and read-modify-write
widgets ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)).

#### The authoring contract: four functions, one optional, one trait

**The authoring contract: four functions, one optional, one trait.** A device is a
user type subtyping the framework's neutral root, `MyDevice <: AbstractDevice`
— one mandatory word that costs nothing (the periphery has no competing
hierarchy to inherit from) and buys `attach!`'s dispatch gate below; the
framework asks for

```julia
init!(dev)          # per-run resource acquisition — calling task, before spawn (§10.4)
loop(dev, handle)   # the task body: owns its own wait structure
shutdown!(dev)      # per-run resource release — guaranteed on every exit path
unblock!(dev)       # optional hook, default no-op: make a blocked loop return (§10.4)
needs_calling_task(dev)   # optional trait, default false: run the loop body on the
                          # calling task (§9.1's topology; the shipped GUI's CImGui
                          # constraint). At most one holder per roster (§9.3).
```

and owns everything around them — the wrapper is the [§10.4](#104-shutdown-protocol) protocol made
structural:

```julia
init!(dev)                                   # its own bracket, pre-spawn (§10.4): a throw
                                             # here is shutdown! + DeviceCrash by name
task = Threads.@spawn try
    loop(dev, handle)
catch e
    report!(DeviceCrash, dev, e)             # §10.4(6): sim continues, device absent
finally
    shutdown!(dev)                           # any exit path: OS resources released
    mark_dead!(...)                          # heartbeat only — claims stay, §9.3
end
```

A `needs_calling_task` device runs the identical wrapper *inline* on the
calling task — the invocation site, not the contract, is its only
difference ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)'s topology, [§10.4](#104-shutdown-protocol)'s join exclusion).

**`shutdown!` must tolerate a partially initialized device.** The release
guarantee holds on the one path *outside* this wrapper too: [§10.4](#104-shutdown-protocol)'s
initialization step brackets each `init!` and hands a device that threw
half-way through acquisition straight back to `shutdown!`, so nothing it did
manage to open is leaked. The obligation that follows is "close only what is
open" — the same defensiveness `shutdown!` already owes the crash path, where a
loop body may die at any point in its own life — and `init!` is correspondingly
*not* asked to clean up after itself: the bracket does once, for every device,
what would otherwise be duplicated in each and enforced in none.

One discrimination in that wrapper: **an `InterruptException` is never a
`DeviceCrash`.** Under the spawned-loop topology the calling task is the one
running a device loop body inline — the GUI's ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)) — so an operator
Ctrl-C raises *there*, inside user code that did nothing wrong. The wrapper
forwards the control-plane stop and lets the body leave through the ordinary
`running(handle)` predicate ([§10.4](#104-shutdown-protocol)(4)): no crash report for what is not a
crash, and no `should_abort` consultation, a stop being already requested.

#### The author owns the loop body; the framework owns the bracket

**The author owns the loop body; the framework owns the bracket.** The fork
is decided by FlightCore's own history: its eight device hooks were never
sufficient alone — the framework loop calling them came in three flavors,
and the *taxonomy* carried the loop-shape information. One device contract
therefore means author-owned loop bodies. A framework-owned hook loop must
ask each device what it waits on — a poll timer, a blocking socket, the
[§10.3](#103-the-next-snapshot-wait) boundary counter — and that declaration is the taxonomy resurrected
as a trait; the bidirectional peer (one socket, one lifecycle — the
no-taxonomy headline case) needs two waits at once, which a hook loop cannot
serve without a select engine; and the GUI's render loop fits no hook set
at all. Under the author-owned body, every wait structure is ordinary user
code composed from handle primitives:

```julia
function loop(dev::T16000M, handle)              # timer-driven, full write-set
    while running(handle)
        sleep(dev.Δt_poll)
        stage!(handle, map_input(poll_axes(dev), binding(handle)))
    end
end

function loop(dev::UDPInput, handle)             # source-driven, data-dependent
    while running(handle)
        datum = recv(dev.socket)                 # blocks; unblock! closes the socket
        is_eot(datum) && return                  # voluntary exit
        try
            stage!(handle, map_input(datum, binding(handle)))
        catch e
            is_datum_error(dev, e) || rethrow()  # a bug → wrapper → DeviceCrash
            report!(handle, MalformedDatum(e))   # garbage → visible, bounded, alive
        end
    end
end

function loop(dev::Telemetry, handle)            # boundary-driven output
    while running(handle)
        snap = wait_next_snapshot(handle)        # §10.3; returns on stop
        send(dev.socket, map_output(gather(handle, snap), binding(handle)))
    end
end
```

A bidirectional peer composes both halves itself — an inner reader task
inside its own domain — rather than forcing a select engine into the
framework. Two idioms are author obligations the framework can only teach
and diagnose, never force ([Appendix A](#appendix-a-taught-contracts-the-author-facing-index)): loop on `running(handle)`, and make
blocking calls interruptible (`unblock!`, or timeouts) — a forgotten
predicate check surfaces as `DeviceJoinTimeout` with the device's name, a
stall as a stale [§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget) heartbeat (liveness timestamps ride *inside* the
handle primitives, which store them in the device's own diagnostic cell
[§9.8](#98-diagnostics-and-liveness-the-per-writer-cell), so the framework observes activity without owning the
loop). **`should_close` dissolves**: a window ✕ or peer EOT is the loop
body returning; the wrapper's exit path releases the device's OS resources,
marks it dead for the heartbeat and consults `should_abort` — claims and
roster entry persist to run end ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s freeze) — [§10.4](#104-shutdown-protocol)(6) is now
literally "the task body returned." The GUI implements the same contract; the framework calls its
`loop` inline on the calling task instead of spawning ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)'s pinning,
unchanged).

#### The binding: framework-legible by enumeration, opaque in its mappings

**The binding: framework-legible by enumeration, opaque in its mappings.**
A binding is a value subtyping `AbstractBinding` — the second mandatory
root — whose type declares which sides it has and enumerates what each side
touches. The legible
half is explicit methods returning data, called once at attach on the
calling task; the opaque half is called per datum on the device task by the
author's own loop:

```julia
struct T16000MBinding <: AbstractBinding    # the roots are mandatory: attach! dispatches
    table::NamedTuple                       # on ::AbstractDevice, ::AbstractBinding
end

is_input(::T16000MBinding)  = true     # sides are *declared*, never inferred; the root
is_output(::T16000MBinding) = false    # carries the false defaults, so silence = absent
is_greedy(::T16000MBinding) = false    # the claim source within the input side (§9.3)

claims(b::T16000MBinding)   = ...      # input side:  the enumerated face set → the claim
map_input(datum, b)                    #              datum → face ⇒ value pairs — user code
reads(b)                               # output side: §14.4 selectors → one compiled gather
map_output(nt, b)                      #              the gather's NamedTuple → wire datum
```

The framework needs no contract on the datum's shape: the datum travels
only between `loop` and `map_input`, written by the same author, and the
framework's structural knowledge comes entirely from the declared traits and
the enumeration methods — everything enumerable validates at attach,
everything opaque is
bounded at its runtime enforcement point (`map_input` by the staging
checks, [§9.4](#94-inbound-per-device-staging-representation-and-the-drain); `map_output` receives exactly the compiled gather's
NamedTuple, and what it puts on the wire is the peer's business).
`map_input`/`map_output` are, precisely, **conventions of the author-owned
loop idiom**: the framework never calls them, so they are taught
([Appendix A](#appendix-a-taught-contracts-the-author-facing-index)) and never checked — a binding whose loop calls something else
by another name is simply a binding with a different private helper.

**Sides are declared; the obligations they create are enforced both ways.**
`claims` and `reads` have **error-throwing fallbacks on the root**, so a
declared side whose method was never written fails loudly at the attach
point rather than degrading into silence, and the attach runs a
**bidirectional conformance check** over the pair (trait, method):

- `is_input && !is_greedy` ⇒ `claims(b)` is called once and its faces staked;
  the fallback firing here means "you declared an input side and wrote no
  enumeration".
- `is_input && is_greedy` ⇒ the claim is computed ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)), and a `claims`
  method defined for this binding is an error: the two sources are
  alternatives, not layers.
- `is_output` ⇒ `reads(b)` is called and the gather compiled.
- Neither side declared ⇒ attach-time error naming both traits: a binding
  that touches nothing is a configuration mistake, not a degenerate.
- `is_greedy` without `is_input` ⇒ error: greediness is a *source* within a
  side, and a source without its side is meaningless.
- A **specific** method of `claims` or `reads` defined for the binding type
  while its trait reads false ⇒ error, the converse direction of the same
  fact: a method written and never reached is exactly the drift the check
  exists to catch. Detecting it is one `which` against the fallback method —
  [§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)'s reflection class (its shadowing check is an `isdefined`/`!==`
  pair), run once at a stopped-sim service point, not inside any frame.

Every violation in that list reports `BindingContractMismatch`
([Appendix C](#appendix-c-the-diagnostic-kind-set)), naming the binding type, the trait, the method at fault and
the direction (declared-but-missing, or defined-but-undeclared). **This is
what closes the shadowing hole**: under detection-by-method-presence, a
bidirectional binding whose `claims` was written without extending the
framework's generic — [§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)'s `using`-without-`import` trap, one level
down — presented as output-only and degraded *silently*: the device
attached, staked nothing, and wrote nothing, with every diagnostic pointing
away from the missing import. With the side declared, the absent method has
something to contradict.

Greediness stays orthogonal to `reads`: a
greedy front end may also drive a
compiled output gather (legal, currently uninstantiated; the plausible
customer is a narrow-wire interactive surface — a motorized control board
whose detents must be driven back out). The binding stays an `attach!` argument, never a device field: the
same `T16000M` binds differently per aircraft, and narrowing the binding
narrows the claim ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)).

**Why the periphery gets roots where components have one.** [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape) refuses
a class supertype for two reasons, and neither reaches here: a component's
single-inheritance slot is *already spoken for* by the domain hierarchies
(`AbstractAircraft`, engine families), while a device's and a binding's are
vacant — nothing else wants them; and a component's class is
implementation detail behind its contract ([§11.3](#113-visibility-the-contract-is-the-interface)), while a binding's
**sidedness is its public contract** — the one thing every consumer of it
must know. Rejected, correspondingly: an abstract binding-type *taxonomy*
encoding the sides (`AbstractInputBinding` and friends — a bidirectional
binding cannot subtype two of them, the original objection, unchanged by the
roots); optional roots left unenforced (a signal half the ecosystem skips is
worse than none — the dispatch gate is the whole value); and a declared
`sides(b)` trait returning the side set, whose rejection stands on its
merits (redundant with methods that must exist anyway, one more thing to
drift) but is **answered rather than repeated** by the design above:
redundancy *with a cross-check* is drift detection, which is what the
bidirectional check turns the traits into — the same fact stated twice, in
two registers, with the framework paid to compare them.

**`is_greedy` is a claim source, not a device class.** The class it used to
justify — a *derived* surface elected by a marker, shared among the
writers that elected it — no longer exists: what the declaration buys is
one computation at the attach point, after which the binding holds an
ordinary claim set and every mechanism downstream ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s exclusivity
and storage, [§9.4](#94-inbound-per-device-staging-representation-and-the-drain)'s shape, shim, merge, scatter and drain, [§9.5](#95-inbound-the-input-trace)'s
trace, detach's release) is blind to where the set came from. The standing
rejection is untouched by that: opportunistic writing to unclaimed faces
"for any device" stays dead (row 44) — autonomous devices still enumerate,
and a maximal surface is now something exactly one line of a binding asks
for, in the open, checked like any other claim.

**A second greedy attach stakes the empty remainder.** The complement is
computed against the roster as it stands, so a greedy binding attached
after another has already swallowed everything gets the empty claim — which
is legal, being the honest may-write-nothing degenerate below, and useless,
which is worth saying out loud: the attach succeeds and reports
`EmptyGreedyClaim` ([Appendix C](#appendix-c-the-diagnostic-kind-set), a service warning naming the device and its
binding), the one honest reading of "you asked for what is left and nothing
was left". No singleton rule survives anywhere here. **Several interactive
front ends may be rostered at once** — a web console claiming the autopilot
faces beside a local GUI claiming the stick faces — because with explicit
claims they are simply two enumerated devices, partitioning the surface
rather than sharing it; the one thing still limited to a single holder is
`needs_calling_task` ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s affinity check), a property of the task
topology, not of interactivity.

**The shipped GUI binding is a greedy one.** It declares `is_input` and
`is_greedy` and stakes the computed claim — everything unclaimed at the moment it attaches —
defining no `claims` of its own; it declares no `reads` either, because its
read path is the handle's primitive
read — VSync-paced, it reads `latest` afresh each render ([§10.3](#103-the-next-snapshot-wait)) with an
ad-hoc, render-time read set over the whole snapshot, the inspection
register's shape ([§9.2](#92-outbound-snapshot-publication)) — so the compiled output gather has nothing to
do for it. The same GUI device type is equally attachable under a binding
that returns explicit claims: greediness is the binding's declaration, not
the device's nature. Every other interactive front end anyone might want (a
web console, a remote panel) has both spellings available — attach with the
greedy binding where the GUI would have been, or with explicit claims
beside other front ends.

**The empty enumeration is not a back door.** `is_input(b) = true` with
`claims(b) = ()` stays an
honest degenerate — a device that may write nothing, its writes
still binding-bounded, so drift onto any face is `OutOfClaimEntry` — and
there is no longer a privileged class for it to promote into: `claims`
bodies are ordinary code
([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)'s idiom, comprehensions included), and an enumeration that came
back empty by accident stays inert, exactly as written. The maximal surface
is reachable only through the explicit `is_greedy(b) = true` declaration —
the most privileged claim is the hardest to acquire by accident, and a
declared trait is deliberate authorship. For the same reason `claims` never
returns `nothing` or a sentinel to mean "compute it for me": the
enumeration contract has one meaning, the trait carries the other, and a
dual-meaning return would be exactly the ambiguity the declaration
vocabulary is built to refuse.

#### One shipped binding type; conditioning has an owner

**One shipped binding type; conditioning has an owner.** `TableBinding` is
*data-driven* — the framework writes its `map_input` once, and a table
value (axis/button entry → face, deadzone/expo parameters) is constructed
per device × aircraft pairing, where configurations are made:

```julia
TableBinding(stick_y  = (face = "elevator", deadzone = 0.05, expo = 0.6),
             throttle = (face = "throttle",),
             trigger  = (face = "brake_count",))       # levels doctrine: a counter
```

Its generic `map_input` *is* [§9.4](#94-inbound-per-device-staging-representation-and-the-drain)'s shared pure conditioning helper, now
with an owner; the entry tuple rides in the type, so the mapping
specializes per table with no dynamic dispatch. A *code-driven* binding (a
JSON telecommand peer: `claims` returns the vocabulary, `map_input` parses
bytes) looks identical to the framework. Purity note, taught in [Appendix A](#appendix-a-taught-contracts-the-author-facing-index):
cross-datum state — press counters, edge detection — lives in the device
struct, maintained by the loop, arriving *inside* the datum; `map_input`
stays pure.

#### Bad datum versus bug: two classes, two fates

**Bad datum versus bug: two classes, two fates.** A datum that cannot be
mapped for environmental reasons — a truncated datagram, malformed JSON, an
out-of-range field — is tolerated *in the loop body*: catch, stage nothing,
`report!(handle, MalformedDatum(cause))`, continue — bounded by the device's
own diagnostic cell ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)'s ring and suppressed counts, [§13.2](#132-diagnostics-structured-values-one-carrier-exception)'s
stream), visible next to a live heartbeat. Any other exception propagates,
and the wrapper turns it into `DeviceCrash` ([§10.4](#104-shutdown-protocol)). The classification is
the author's — only they know their parser — exactly as FlightCore's
`InputMappingError` docstring assigned it; what changes under the
author-owned loop is that no framework per-iteration catch site exists, so
the framework's contribution is the diagnostic channel, not the catch (a
marked exception type with no framework consumer would be vestigial and is
not provided). `report!(handle, ...)` writes device-attributed runtime
warnings into that device's diagnostic cell — the [§13.2](#132-diagnostics-structured-values-one-carrier-exception) stream's
single-writer entry point ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)) — and nothing more; it is not a general
user-diagnostics channel. Tolerating everything hides bugs as "device
attached, nothing happens"; tolerating nothing kills a live telemetry link
on its first truncated datagram — and since tasks are per-run artifacts
(row 93), kills it for the rest of the run.

### 9.7 The GUI write path: port resolution, peek, staging contract

Panels remain per-component extensions in FlightCore's style — `GUI.draw!(ctx,
::LowPassFilter)`, discovered by walking the assembly — but widgets name **the
component's own ports**, never root slots. The build-time wiring answers, statically
and exactly: *is this port transitively driven by a root input slot, and which one?*
Every input port has exactly one source ([§6.1](#61-connections-and-hierarchy)), so the resolution is total:

- **root-driven → live widget**: peeks and stages the resolved slot through the
  GUI's own staging cell;
- **component-driven → read-only rendering**: displays the driven value from the
  snapshot, visually distinct, with the source as provenance ("driven by
  `avionics/throttle_cmd`" — the canonical slash form of [§11.6](#116-paths-wiring-and-faces)).

This retires FlightCore's dead-slider convention (the `Cessna172Xv1` throttle: the
engine panel's slider visually live, silently overwritten by the avionics every
cycle — who commands what living in the user's head) and replaces it with checked
structure: **a widget is live exactly when the underlying input is yours to command
in this configuration.** User-commandability is a wiring decision made where
configurations are made; command-plus-manual-override is a mux component with a
root-wired select — explicit structure, not two writers racing (the same race as
[§15.3](#153-torture-test-for-the-9-staging-shapes-filter-joystick-and-gui)'s drag phase, ruled out the same way). The obligation this places on the
GUI: read-only rendering is first-class, not an error state — the author of
`input_slider!` cannot know at authoring time whether it will be live.

**Liveness is a derived property, and resolution is transitive.** A widget
is live iff its port's feed chain — walked through wires and boundary connections across *all*
levels, not just the local assembly — terminates in a root slot, *and* that slot
lies **inside the GUI's own claim** in the run's frozen surface partition
([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster) exclusivity) — whether that claim was computed from the unclaimed
complement under `is_greedy` or enumerated face by face by a partial-claims
binding ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)), so
"live" reads as "inside the surface I declared for". Under
the roster freeze, liveness is a static fact of the run: baked once, with the
port resolution, when the run starts — never consulted against mutable claim
state at render. There is no per-port
"GUI-controlled" marking anywhere: the export chain is the marking, written by
the one author entitled to write it (a component's ports become GUI-commandable
exactly when the assemblies above surface them). The switch between "driven by
its own panel" and "driven by an external provider" is therefore automatic — at
build time by wiring archetype (a scripted `World` wires a scenario component
into the same faces the interactive `World` exports to root), at run start by
roster claim state. Rejected: nominally-connected ports with a GUI
*override* channel — a second write path that breaks the
pure-function-of-drained-batches frame semantics, needs a parallel trace and
replay mechanism, and cannot resolve the producer conflict (either a dead widget
or a silently discarded wire); made safe — staged, frame-top, traced, exclusive —
it collapses into the root-slot mechanism it tried to bypass. The honest cost
stands: **unexported ports are unpokeable** — FlightCore's poke-any-`u` workflow
does not survive contract visibility ([§11.3](#113-visibility-the-contract-is-the-interface)), deliberately.

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
motivation to slot exclusivity ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)): stage-every-pass would exist to win every
drain against a streaming device sharing the slot for the grab's duration ([§15.3](#153-torture-test-for-the-9-staging-shapes-filter-joystick-and-gui)'s
drag phase), but a slot the GUI can write is exclusively its own claim, so once staged
and drained a value simply *stays*; there is nothing to reassert against. Nor is it
worth keeping as insurance: if an
anomalous writer ever touched a claimed slot through a framework defect, the
slot visibly fighting is a diagnosable anomaly — continuous re-staging would mask
the invariant violation at render rate. Side benefit: staging traffic (and trace
noise) drops from render-rate-while-grabbed to actual edits. No
claim-transition policy exists, because no claim transition can occur
mid-run ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s freeze); the one liveness-adjacent display rule is the
orphan case — a read-only widget whose claiming device's task has died
renders the fact in its provenance ("claimed by `T16000M` — task dead",
[§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget)'s heartbeat surfaced in place), so an orphaned slot is visible where
the user is looking, not only in the status panel.

The panel-authoring calling convention — what the drawing context carries,
how widgets name their component's ports, how an assembly's panel composes
its children's — is deferred to migration ([§16](#16-open-axes)), where it is co-designed
against the GUI library. Its constraints are fixed here: panels name their
own ports by face-name string; resolution to root slots and the liveness
verdict are baked at run start, never performed at render; liveness and peek
arrive through the framework-supplied context, never by reaching into the
loop; and assembly panels compose children by path.

### 9.8 Diagnostics and liveness: the per-writer cell

The chapter's two data channels are specified down to their memory ordering.
The runtime warning stream ([§13.2](#132-diagnostics-structured-values-one-carrier-exception)) and the liveness heartbeat ([§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget)) are
a third, and they cross the same task boundaries: they are written by the
device tasks — `OutOfClaimEntry`, `ClaimedFaceEntry` and `EntryTypeMismatch`
at staging ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain)), `MalformedDatum` from the author's loop body via
`report!(handle, …)` ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)) — and by the loop itself (`ChatteringBudget`,
`EventDeferred`, `DebtReanchor`), and read by the loop, which folds them into
the published framework status ([§9.2](#92-outbound-snapshot-publication)) and hence into every snapshot. An
unspecified structure with those writers is exactly the arbitrary shared
mutable state [§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)'s two rules exist to eliminate, so it gets the mechanism
[§9.4](#94-inbound-per-device-staging-representation-and-the-drain) already established, not one of its own.

**One diagnostic cell per writer — one per rostered device, one for the loop
itself.** Each cell has a single writer, the same ownership argument as the
staging cells: no locking, no arbitration, no new primitive. The cell holds a
**bounded accumulation** — a small ring of diagnostic values, capacity **16**,
plus a per-kind count of what the ring could not hold — and one atomic
liveness timestamp.

**That bound *is* the rate limit.** When a writer emits past the ring's
capacity within one frame, the entry is not stored and its kind's suppressed
count increments; the drop policy is earliest-in-frame retained, excess
becomes counts — the first occurrences are the ones with diagnostic content,
the hundredth is noise the count already reports. [§13.2](#132-diagnostics-structured-values-one-carrier-exception)'s "rate-limited
wherever its source can repeat" is therefore not a policy layered over the
stream but a structural property of the channel that carries it: a chattering
model or a peer flooding malformed datagrams costs at most sixteen retained
values and one integer increment per kind per frame, whatever its source
does, and no writer can starve another — the cells are disjoint.

**The drain is [§9.4](#94-inbound-per-device-staging-representation-and-the-drain)'s drain.** One `atomicswap` per cell at frame top, at
the same point and under the same indivisible-take argument as the staging
drain; what the loop swaps *in* is a shared **empty sentinel**, so a quiet
frame swaps the sentinel in and gets the sentinel back — no allocation, and
no load-only code path that goes untested on healthy runs. The take is also
what makes publication sound: the batch is exclusively the loop's before it
is ever reachable from a snapshot, so [§9.2](#92-outbound-snapshot-publication)'s binding rule — nothing
reachable from a published snapshot is ever written again — holds by
construction, and the live accumulator is never reachable from a published
value.

**The heartbeat rides in the same cell**, as an atomic timestamp field the
device task stores on every loop pass from inside the handle primitives
([§9.6](#96-devices-one-authoring-contract-no-taxonomy): the framework observes activity without owning the loop body) and
the loop acquire-loads at the drain. There is no separate liveness channel
and no second registry: a device that is alive is a device whose cell carries
a recent timestamp, and [§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget)'s 2 s staleness threshold is read against
this field. The heartbeat is not a diagnostic kind — it is a field, always
present, never enumerated in [Appendix C](#appendix-c-the-diagnostic-kind-set).

**The published framework status is a concrete frozen value.** Per writer it
carries `recent` — the ring this boundary drained, at most sixteen entries;
`suppressed` — the per-kind counts the ring refused this boundary; `totals` —
the cumulative per-writer × per-kind counts since the run began, owned
privately by the loop and *copied* into each status; and `heartbeat`. Beside
the per-writer records ride [§8.7](#87-real-time-pacing)'s pacer diagnostics, as before. Delta
plus total is what makes the status legible at any reading cadence: a GUI
panel refreshing at 60 Hz sees each occurrence once in `recent`, while a
consumer that samples occasionally still reads a complete account from
`totals` — nothing is lost by not looking.

**Presentation is where the old `maxlog` lives.** A status renderer prints a
given writer × kind up to **25** cumulative occurrences and then switches to
count-only display ("`MalformedDatum` from `UDPInput#3`: 1 482 occurrences").
That threshold is presentation policy, not channel policy: counts keep
accumulating regardless, nothing recorded depends on it, and the choice
belongs to whoever renders. The channel's own bound, above, is the one that
is normative.

**The terminal snapshot carries the run's final cumulative counters**
([§10.4](#104-shutdown-protocol), [§13.5](#135-termination-is-a-state-not-an-exception)), so an unattended run that nobody watched still
answers "what went wrong, and how often" from the value its own shutdown
published.

**Allocation.** On a quiet frame there is **zero additional heap
allocation**: the sentinel swap allocates nothing and the per-writer status
rides inline in the one per-boundary snapshot allocation [§9.2](#92-outbound-snapshot-publication) already
accepted. That requires the per-kind counters to be a **fixed-shape isbits
record, never a `Dict`** — licensed by [Appendix C](#appendix-c-the-diagnostic-kind-set)'s closed kind set, which
makes the counter layout a type rather than a lookup. On a noisy frame the
diagnostic values are allocated at emission, on the writer's own task; a
drained non-empty ring is frozen into the snapshot and can never be written
again, so the writer allocates a fresh ring lazily at its next emission —
that cost, too, landing on the writer's task — which is the same
GC-over-reuse trade [§9.2](#92-outbound-snapshot-publication) makes when it rejects preallocated snapshot
buffers. The rate limit is therefore an allocation bound as well: one ring of
sixteen entries per writer per boundary is the worst case, everything past it
an integer increment. [§7.5](#75-allocation-policy-a-scoped-invariant)'s zero-allocation invariant, scoped to the model
sweep, is untouched — the cells sit on the framework side of that scope with
publication and logging.

Composition with the log, worth stating once: because the log retains
snapshot references ([§9.2](#92-outbound-snapshot-publication)), `totals` is monotone across logged snapshots,
so `log_every` decimation loses *which* boundary within a skipped stretch an
occurrence fell on, never *how many* occurrences there were.

Rejected: **a shared queue under a lock** — cross-task contention on the
loop's own frame path, and no single-writer ownership to reason from, the
shape [§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads) refuses wholesale; **a status referencing the live
accumulator** — the cheapest thing to write and a direct violation of the
binding rule, a reader walking a structure the writer is still mutating;
**ring reuse by double-buffering** — reintroduces exactly the reader-liveness
proof the GC provides free, [§9.2](#92-outbound-snapshot-publication)'s own argument against preallocated
snapshot buffers, to save an allocation that occurs only on frames already
doing something unusual; **unbounded accumulation** — an unattended run with
no status reader grows without bound in the configuration the framework most
has to survive, and the drop policy is what a bound buys.

---

## 10. Runtime periphery: lifecycle and orchestration

Where [§9](#9-runtime-periphery-the-data-plane) fixes how data crosses the loop boundary, this chapter covers the
machinery that drives the loop itself: the control plane and the scheduling
primitives, the shutdown protocol, and the run lifecycle from `init!` through
replay.

### 10.1 Control plane

Pause, un-pause, pace changes, `margin` changes, stop: a few scalar fields on a
separate atomic
surface, consulted by the loop at frame top and inside its wait and pause states.
`margin` ([§8.7](#87-real-time-pacing)) rides here for the same reason `pace` does — it tunes the
wait, never the arithmetic — so retuning the coarse/spin split mid-run is safe
by construction. The stop's issuers are the operator's channels — GUI button,
device handle, calling code — and, in an interactive session, Ctrl-C: an
operator interrupt is caught at one of the loop's unmask points and sets exactly
this stop, no separate entry point involved ([§10.4](#104-shutdown-protocol)).
**Not staging, structurally:** staged writes apply at drains, and a paused loop
drains nothing — un-pause via staging would deadlock by construction. Riding outside
the drain/trace path is safe for determinism precisely because [§8.7](#87-real-time-pacing) put pacing
outside the semantics: control changes *when* frames execute, never what they
compute (stop merely truncates the trajectory). While paused the loop blocks on a
condition (notified on un-pause and stop), not a spin.

### 10.2 Loop scheduling: wait primitive, yields, thread budget

**Coarse phase = task-yielding `sleep`; no `systemsleep` variant.** The precision
argument for `Libc.systemsleep` (≈0.5 ms vs ≈1.4 ms median overshoot, [§8.7](#87-real-time-pacing)) is
worth ~1.5 ms of `margin` — a few percent of one core at 50 Hz. Against it: `sleep`
releases the loop's thread, making the pacer's wait slot the natural scheduling
window for co-resident device tasks (the design already spends that slot twice:
[§8.7](#87-real-time-pacing)'s staging slot, [§9.4](#94-inbound-per-device-staging-representation-and-the-drain)'s drain source); `systemsleep` turns the slot into dead
time and makes the periphery correct only when every device task has its own
thread — resurrecting FlightCore's hard `nthreads` check as a correctness
precondition. The failure asymmetry decides: `sleep`'s worst case is a late
wake-up → an overrun → absorbed as debt and *diagnosed*; `systemsleep`'s is starved
device tasks → silent functional degradation. And [§8.7](#87-real-time-pacing) committed to `margin` as
the single knob; a primitive selector is a second knob hiding inside the first (its
two settings differ by a margin recalibration). A `systemsleep` variant for
dedicated-thread hard-RT deployments is a guarded addition.

**Yield rule: with devices attached, every frame yields at least once** —
implicitly via the coarse-phase `sleep` when it runs, via an explicit `yield()`
otherwise (unpaced runs; pure-spin frames with budget ≤ margin). Zero semantic
cost: pacing, hence yielding, is outside the semantics ([§8.7](#87-real-time-pacing)). The spin phase
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
tight for the attached population — one check per run, against the frozen
roster ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)) — naming the `julia -t` remedy; sizing guidance:
one thread for the loop, the main thread for the GUI, headroom for compute-heavy or
blocking-ccall devices (libuv-backed I/O yields; raw blocking ccalls pin their
thread for the duration). No pinning, no sticky tasks.

**Liveness heartbeat.** Since starvation is survivable it must be diagnosable: the
published framework status includes per-device liveness (last-staged / last-read
wall time, task state) next to the pacer diagnostics. The mechanism is
[§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)'s per-writer cell and nothing besides — an atomic timestamp field
the device task stores on every loop pass from inside the handle primitives
([§9.6](#96-devices-one-authoring-contract-no-taxonomy)) and the loop acquire-loads at the frame-top drain, alongside that
device's diagnostics. A starved, blocked or crashed
device task shows in the GUI as a stale heartbeat with a name on it, not as
mysteriously frozen physics. **Stale means a liveness timestamp more than 2 s
behind wall clock** — deliberately loose, because the heartbeat is advisory
(a liveness display and a provenance record, never a kill trigger, never a
detach) and must tolerate a device legitimately parked in a blocking read
between rare data.

### 10.3 The next-snapshot wait

Rate-matched output devices (telemetry, disk streaming) act once per boundary
without polling: a monotonic **boundary counter** published with the
snapshot — counting *published boundaries* (grid, `t*`, boundary zero; [§8.4](#84-localization-mechanics)),
not frames, so consecutive wakes are not necessarily `h` apart — plus one
`Threads.Condition`. The loop's publication is `lock; counter += 1; notify;
unlock` — nanoseconds of framework-only code, never blocked by waiters (a waiter
parked in `wait` has released the lock as part of parking). The device-side
`wait_next_snapshot(handle)` blocks until `counter > last_seen && running` under
the canonical predicate-loop idiom, which handles waiters at different paces,
frames skipped while transmitting, and shutdown ([§10.4](#104-shutdown-protocol) wakes all waiters; each
predicate routes its owner out) with no per-frame reset. An `Event` latch is the
wrong primitive here: recurring signals require un-latching, and the reset has no
correct placement under asynchronously arriving waiters (the `io_start` reset
comments in FlightCore's sim.jl document the once-per-run version of exactly this
race). Conditions carry no facts, only "look again"; the facts (counter, running)
live in state each waiter tests privately.

**Counter home and publication order.** The boundary index is carried *in* the
snapshot (with `t`), so any holder of one — the log, an error's replay pointer
([§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)), a post-run inspector — indexes it without consulting the loop; the
loop additionally mirrors it in the state the wait predicate tests. The order
of the two publications is normative: the release-store of `latest` ([§9.2](#92-outbound-snapshot-publication))
happens **before** the counter increment under the lock, so `counter >
last_seen` implies `latest` holds at least that boundary — a waiter can never
wake onto a stale snapshot. The converse — observing a *newer* snapshot than
the increment that woke you — is expected and correct: newest-wins.

**Semantics: newest-wins, no queues.** A slow consumer skips frames and always
receives the current world. This mirrors the inbound side: coalescing to the
newest batch (in) and to the newest snapshot (out) are the same ZOH decision; no
backpressure exists in either direction, and the loop never waits on anyone.
Rejected: per-consumer every-boundary queues (unbounded under a slow consumer —
the batch stack's pause pathology again; complete history *is* the log). The GUI
does not use the wait (VSync-paced, it reads `latest` each render).

### 10.4 Shutdown protocol

1. **Initiation:** `t_end` reached, a control-plane stop (GUI, device handle,
   code, or an operator interrupt — Ctrl-C, below), or a `stop_on` face
   reading `true` in the just-published snapshot
   (model-detected termination, [§13.5](#135-termination-is-a-state-not-an-exception)). The loop always completes the current
   boundary sequence — never stops mid-frame — publishes the final snapshot,
   then sets the sticky stopped status.
   Publishing first guarantees output devices can flush the true final state,
   and that final status carries the run's cumulative diagnostic counters
   ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)) — the complete warning account of a run nobody watched. That
   terminal snapshot is retained in the log unconditionally, under any
   `log_every` and any `log_max` ([§9.2](#92-outbound-snapshot-publication)).
   **`t_end` lands on the grid:** the run ends at the first grid boundary whose
   time reaches or exceeds `t_end` — whole frames only, never a shortened final
   step, which [§8.4](#84-localization-mechanics)'s grid integrity (`tₖ = t₀ + k·h`, indexed and never
   accumulated) forbids. The final boundary may therefore overshoot `t_end` by
   up to `h`, and the termination record carries the actual final `t`
   ([§13.5](#135-termination-is-a-state-not-an-exception)). This is [§10.6](#106-run-lifecycle-and-partial-advance)'s `t_plus` spelling — whole frames until
   the boundary time first covers the duration — applied to the run's own
   clock, and it is where the two termination sources differ in kind: `t_end`
   is a grid fact, checked against boundary times on the grid, while `stop_on`
   is checked at *every* published boundary, `t*` included ([§13.5](#135-termination-is-a-state-not-an-exception)).
2. **Wake all framework waits** (next-snapshot, pause): waiters observe the
   stopped status and unwind — a stop while paused therefore works.
3. **Unblock device-specific blocking calls** via an `unblock!(device)` hook,
   default no-op; a network input's override closes its own socket, raising in
   the blocked task (caught by the framework wrapper, treated as shutdown).
   This demotes FlightCore's EOT convention from load-bearing shutdown mechanism
   to an optional wire-protocol courtesy between remote peers.
4. **Loop bodies exit:** the author's `while running(handle)` ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)'s
   contract — the predicate check and interruptible blocking are the two
   taught obligations) with all blocking points
   interruptible per (2)–(3); the wrapper's `finally shutdown!(device)` is
   guaranteed on every exit path.
5. **Join with a timeout — 5 s:** a device task exceeding it is reported *by name*
   ([§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget) heartbeat) and abandoned with a warning (`DeviceJoinTimeout`,
   [Appendix C](#appendix-c-the-diagnostic-kind-set)) rather than hanging `run!`. Five seconds is generous for
   GUI window teardown and socket closes, and short enough that an abandoned
   join reads as a diagnosed timeout rather than a hang.
   The calling-task device — the GUI — having no spawned task ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)), is
   outside the join: its loop body is the calling task's own occupation of
   `run!`, exits by the same `running(handle)` predicate as any device
   loop, and `run!` returns after the joins. One honest asymmetry: the
   abandonment protocol cannot cover it — nothing can abandon the task
   `run!` stands on — so a calling-task device that blocks past shutdown
   hangs `run!`. The trait's one authoring obligation is therefore a loop
   body that never blocks between `running` checks; the shipped GUI's
   render loop polls once per frame and never blocks.
6. **Device-initiated paths:** voluntary exit — the loop body returns
   (window ✕, peer EOT; no `should_close` hook exists, [§9.6](#96-devices-one-authoring-contract-no-taxonomy)); with
   `should_abort` set the wrapper's exit path also requests a sim stop,
   otherwise the sim continues with the device's *task* absent: its cell
   stops filling — the loop is structurally indifferent — while its roster
   entry and claims persist to run end ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s freeze: death is not
   detach; the orphaned slots hold their last-drained values, visibly,
   [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)). A crashing device task is caught by the
   framework wrapper and follows the same path, logged with the device's name
   (`DeviceCrash`, [Appendix C](#appendix-c-the-diagnostic-kind-set)).
7. **Loop-side failure** runs (1)–(5) from the catch path — specified in
   [§13.6](#136-abnormal-shutdown-one-tail-two-entries): the failed boundary is discarded and the previous snapshot promoted
   to final (FlightCore's `SimulationTermination` catch path was the precedent;
   the exception-based termination idiom itself has no place here, [§13.5](#135-termination-is-a-state-not-an-exception)) — so devices
   unwind cleanly regardless of who died.

After (5) the task set is empty — device tasks are per-run artifacts
([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)) — and `shutdown!` has released each device's OS resources. What
survives a stop is the roster entry: binding, claims, stable device id
([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)); never a task, never a live resource — a device whose task died
mid-run included, its entry indistinguishable at this point from any
other's; `stopped` is where `detach!` removes it and releases its claims.
**One roster change belongs to this tail**: a GUI attached by `run!`'s
`gui = true` is detached here, releasing its computed claim ([§10.6](#106-run-lifecycle-and-partial-advance)'s
run-scoped flag) — the only roster mutation the protocol itself performs, and
it sits in the tail precisely so that (7)'s failure path takes it too, an
everything-claim staked for one run never surviving into the next.
The next `run!` re-runs device
`init!` — resource acquisition is per-run; FlightCore's
create-a-new-socket-each-`init!` in network.jl is the precedent — and spawns
fresh tasks against the re-armed [§10.3](#103-the-next-snapshot-wait) counter. While stopped there are no
device tasks, so voluntary exit and the [§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget) liveness heartbeat are
run-scoped observables; a device unplugged while stopped surfaces as the next
run's `init!` failure, disposed of by the initialization step below.

**Initialization is a step of this protocol, taken at the top of a run.**
Before any task is spawned, the loop calls `init!` once per roster entry, in
attachment order, on the calling task — and each call sits in its own bracket:

```julia
for entry in roster                       # attachment order, calling task
    try
        init!(entry.device)
    catch e
        shutdown!(entry.device)           # release, unconditionally (§9.6)
        report!(DeviceCrash, entry.device, e)
        mark_dead!(entry)                 # from boundary zero; no task is spawned
        entry.should_abort && stop!(control)
    end
end
```

The bracket is what makes [§9.6](#96-devices-one-authoring-contract-no-taxonomy)'s "guaranteed on every exit path" true of
the path outside its wrapper: a device that throws half-way through acquisition
is handed back to `shutdown!` right there, so its partially acquired OS
resources are released rather than leaked — which is exactly why `shutdown!`
owes tolerance of a partially initialized device ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)'s taught obligation).
The report is the ordinary `DeviceCrash` ([Appendix C](#appendix-c-the-diagnostic-kind-set)), not a kind of its own:
its payload — device id, the cause exception, whether `should_abort` was set —
already carries everything an init-time failure has to say, and the name is
honest, a device that cannot acquire its resources having crashed before it
lived. No task is spawned for a failed device, so it is **dead from boundary
zero**, and that needs no machinery: its diagnostic cell never receives a
heartbeat timestamp, so it reads stale against [§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget)'s threshold from the
first frame ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)). Its **claims persist to run end** — [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s
death-is-not-detach disposition, applied one step earlier than (6)'s: the
roster is frozen for the run, and the orphaned slots hold their initial values,
well-defined by [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)'s slot totality where an orphan of (6) holds a last
drained batch.

**The run's disposition splits on `should_abort`, uniformly with (6).** Clear —
the default ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)) — the remaining entries initialize, the run starts, and
the sim runs with that device absent from frame zero: (6)'s "the sim continues
with the device's task absent", shifted to `t₀`. Set, the failure requests a
control-plane stop, and that stop is simply *already pending* when the run
reaches boundary zero — a path this protocol already has, since [§13.5](#135-termination-is-a-state-not-an-exception)'s
boundary-zero check ends a run at `t₀` with that snapshot final, integrating
nothing. No new exit protocol, therefore: the remaining entries still
initialize — every rostered device gets its `init!`/`shutdown!` pair, uniformly
— the run publishes boundary zero, and it ends `stopped` at `t₀` through this
same tail, with the termination record naming the source ([§13.5](#135-termination-is-a-state-not-an-exception)). What the
operator is left with is an ordinary stopped simulation: a terminal snapshot,
the failure named in its diagnostic account, fully serviceable by [§14](#14-stopped-sim-services) and
resumable by the next `run!` ([§10.6](#106-run-lifecycle-and-partial-advance)) once the device is plugged back in.

**Topology is derived after initialization**, not from the roster alone
([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)): a `needs_calling_task` holder whose `init!` failed returns the loop
to the calling task, which would otherwise be pinned waiting to run the loop
body of a device that does not exist. The shipped GUI attaches with
`should_abort = true`, so in practice that run ends at `t₀` anyway; the rule is
stated generally because it costs nothing.

**The operator interrupt is a stop, not a failure.** Ctrl-C in an interactive
session is a control-plane stop command issued by hand: the run completes the
current boundary, publishes the final snapshot, takes this tail like any other
stop, and ends `stopped` — boundary-consistent, fully serviceable by the [§14](#14-stopped-sim-services)
stopped-sim services, resumable by the next `run!` ([§10.6](#106-run-lifecycle-and-partial-advance)). It is the escape
from a run nothing else can end — deviceless, `t_end = Inf`, empty `stop_on`
([Appendix C](#appendix-c-the-diagnostic-kind-set)'s `UnboundedRun`) — and it needs no entry point of its own: the
control plane already carries the stop ([§10.1](#101-control-plane)), and [§13](#13-error-discipline)'s
exceptions-are-abnormal doctrine is untouched, being about *model* code, while
this is the one exception whose meaning the framework knows.
**Masking across the boundary is normative**, not an implementation hint. An
`InterruptException` is delivered asynchronously, so an interrupt landing
mid-sweep would destroy the boundary this protocol is built on completing and
leave half-written stores ([§13.6](#136-abnormal-shutdown-one-tail-two-entries)) — forcing a choice between `stopped` with
dirty stores and a terminal `errored`, and the `stopped`-with-consistent-stores
guarantee is exactly what the masking buys. The loop therefore masks delivery
across the boundary macro-sequence (Julia's `disable_sigint`; a sigatomic
counter increment, negligible per frame) and takes the deferred raise at the
unmask points — the frame top, where it already consults the control plane
([§10.1](#101-control-plane)), and inside its wait and pause blocks — all boundary-consistent.
Caught there, it sets the control-plane stop and enters this tail; [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)'s catch
site therefore never sees it. A **second interrupt during the tail** collapses
the remaining joins immediately — (5)'s abandonment path taken at once, devices
still reported by name (`DeviceJoinTimeout`) — and the run still ends `stopped`:
escalation shortens the tail, never reclassifies the run. It cannot repair (5)'s
honest asymmetry either, since nothing can abandon the task `run!` stands on.
**Interactive-session scope**, stated plainly: outside the REPL, Julia's default
(`exit_on_sigint(true)`) kills the process on SIGINT before any of this
machinery runs, and the framework flips nothing process-global — unattended runs
rely on `t_end` and `stop_on`, as they already must.

### 10.5 Scripts and the mid-run mutation doctrine

What the consumers demonstrably mutate mid-run, surveyed: FlightCore's
`user_callback!` has exactly two archetypes — the timetable script
(c172_demos.jl:290: `elevator_offset` as a function of `t`) and the synthetic
pilot (c172_demos.jl:423, 525: a phase FSM reading `y` and writing mode requests,
references, flaps, wind). Both write only `u` fields; no demo, test or GUI path
pokes `x`/`s` mid-run, and `init!`/trim appear only between construction and
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
honest `x` (visible in snapshots, logs and plots); inputs arrive same-boundary
fresh by topological order (the callback ran post-step, one boundary staler);
the pure timetable script is a one-liner reading its bundle's clock
(`h_xu(s, (; t)) = (; offset = profile(t))` — exact at its own ticks, no
latching); and
in a scenario configuration the script drives the avionics' input ports, so [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)
renders the corresponding GUI widgets read-only with provenance — today's
demo-vs-GUI dead-slider fight, resolved by the port-resolution rule.

**`user_callback!` is eliminated.** It is the periphery's `f_step!`: arbitrary
unrecorded mutation, ordered by convention, invisible to replay ([§2.2](#22-exclusions-deliberate)). Its
historical justification was FlightCore's composition cost — a supervisor required
a full `System` declaration against a ten-line closure; this framework prices a
component at roughly the closure's weight, removing the pressure. Its call sites
migrate to scenario components, not devices.

**Manual event triggering needs no mechanism:** a root input slot plus a boundary-detected
guard reading it (levels doctrine: latched commands or counters), already
expressible in settled machinery — the demos' engine start/stop buttons are
`u`-writes today.

**Mid-run re-initialization is not built, because it is not demonstrated.**
Initialization and trim are stopped-sim workflows ([§14](#14-stopped-sim-services)'s first-class services),
where no concurrency perimeter exists — no loop, no devices, plain single-task
code. The guarded-addition shape is on record should demand appear: a traced,
boundary-executed intervention command applied through project → sweep → publish,
so no consumer ever observes un-decoded state.

**The doctrine, final form:** while a simulation runs, the periphery stages
root-input writes and issues control commands — nothing else, structurally.
Anything that wants to poke the model mid-run is an *input* in disguise (wire a
slot and guard), *model behavior* in disguise (add a scenario component), or a
*wall-clock interaction* (attach a device). Graceful termination follows the
same shape ([§13.5](#135-termination-is-a-state-not-an-exception)): a declared stop face in the model plus `stop_on` policy at
deployment — never a callback, never a thrown exception.

### 10.6 Run lifecycle and partial advance

A `Simulation` moves through five states: **built** (stores allocated,
nothing authored), **initialized** (`init!` completed boundary zero, [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)),
**running**, and terminally **stopped** or **errored** ([§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)). `init!` is
mandatory: `run!` or `step!` on a simulation whose boundary zero has not
completed is an error in [§13.2](#132-diagnostics-structured-values-one-carrier-exception)'s kind set naming `init!` — distinct from
`UninitializedSlots`, which fires *inside* `init!` ([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)). (`replay!` is
the one alternative entry: it runs boundary zero from a trace header,
[§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop).) The loop runs on
the calling task unless a calling-task device — the GUI — is rostered
([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)'s roster-derived topology); deviceless, `run!` is fully
synchronous — the unattended
register ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads): an unattended run is the same loop with empty staging), and what
[§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)'s synchronous rethrow presupposes.

**Partial advance.** `step!(sim; frames = 1)` advances whole frames
synchronously through the ordinary frame sequence — drain, integrate,
boundaries, publication — and returns; a stepped simulation is bit-identical
to the same frames under `run!`. `step!(sim; t_plus = 10.0)` is the duration
spelling, mutually exclusive with `frames`: whole frames until the boundary
time first covers the duration — the migration suite's advance-by-duration
idiom. This is the test-harness register (advance,
assert, advance) and the REPL register (fly a while, inspect, continue);
neither is a script, so [§10.5](#105-scripts-and-the-mid-run-mutation-doctrine)'s scenario-component doctrine does not absorb
them.

A stepping session is **deviceless by construction**: device tasks are
per-`run!` artifacts ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)) and a device loop's `while running(handle)` is
false outside a run. Between `step!` calls the simulation is in a
stopped-sim state (`initialized`, below), so `attach!` is legal there and
does what it always does — registers ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)); the task appears at the next
`run!`. The frame-top drain still runs — `step!` frames stay
bit-identical to `run!` frames — and what it drains is the **harness
cell**: `stage!(sim, "face" => value, ...)`, [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s harness
write path with the calling task as writer. Staged batches are ordinary
batches — traced, so replay and bit-identity hold; applied at the next frame
top; surface-checked like any writer's ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)). The read half is
`latest(sim)`: the same immutable snapshot value a device handle acquires
([§9.2](#92-outbound-snapshot-publication)), navigated directly for assertions. Advance-assert-advance is
`stage!` → `step!` → `latest`. Both entry points work under `run!` too — the
harness cell is not step-scoped — and an inspection accessor leaves
[§13.5](#135-termination-is-a-state-not-an-exception)'s rejection of closure-based termination untouched.

**Status, termination and the `run!` seam.** Between `step!` calls a simulation
reports **initialized**: no loop task exists, so `running` would lie, and
nothing is terminal, so `stopped` would too — the state reads "boundary-
consistent and ready to advance", not "sitting at boundary zero". `run!` may
therefore follow `step!`, continuing from the current boundary; so may another
`step!`. Termination policy is honored throughout, as bit-identity requires:
`t_end` reached, or a `stop_on` face holding at frame 3 of `step!(sim; frames =
10)`, ends the run there through the ordinary [§10.4](#104-shutdown-protocol) tail and leaves the
simulation `stopped`. `step!` therefore returns the number of frames
**actually advanced** — the requested count in the ordinary case, fewer when
the run terminated inside the call, which is how a harness detects the
truncation without inspecting the clock.

**Re-running.** `stopped → init! → run!` is the supported cycle: `init!`
re-runs boundary zero from its condition (warm restart = `capture` → tweak →
`init!`, [§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)) and clears the trace, the log, *and* any batches still in
staging cells — the recorders restart with the run they record, and
no stale batch survives to clobber the boundary zero it predates. Device attachments persist across
re-initialization — attachment is orthogonal to the run lifecycle ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)) —
and persistence means *roster* persistence: binding, claims and device id
survive; tasks and OS resources do not ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)'s per-run topology, [§10.4](#104-shutdown-protocol)'s
teardown). Each `run!` re-initializes every rostered device and spawns its
task; `attach!` while stopped only registers — the task appears at the next
`run!`. Task topology follows the roster each time ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)): a GUI
attached *by hand* is still rostered, so the next `run!` renders it again —
loop on a spawned task — whether or not `gui = true` is repeated. **The flag
itself is run-scoped**: at run entry it attaches the standard GUI device under
the greedy binding, with `should_abort = true`, iff no GUI is rostered
([Appendix B](#appendix-b-api-synopsis-the-entry-points)), and the run's shutdown tail detaches it again ([§10.4](#104-shutdown-protocol)) —
so the roster a flagged run leaves behind is the roster it found, and a window
on every run means the flag on every run. A *persistent* GUI session is
spelled by hand — `attach!` while stopped, `detach!` when done — and against a
hand-attached GUI the flag does nothing and detaches nothing, having attached
nothing. What the scoping buys is the absence of a trap: the flag's GUI claims
everything unclaimed at attach ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s computed source), and a claim of
that shape must not outlive the run that asked for it — a joystick attached
between two runs would otherwise meet a `ClaimConflict` against an
everything-claim staked by a convenience argument nobody remembers passing.
The accepted cost is a fresh device id per run for that GUI: ids exist to be
read *across* roster changes, and each run's trace header carries its own
schemas ([§9.5](#95-inbound-the-input-trace)), so nothing that reads a completed run is affected.
Run policy is re-bindable per cycle: `t_end` and `stop_on` are `Simulation`
defaults that `run!` may override for the run it starts ([§13.5](#135-termination-is-a-state-not-an-exception)), so a second
run — or a `step!` register between two runs — can stop on a different clock
or a different face set without a rebuild. `errored` is terminal (row 59):
reproduction is trace replay ([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)), not resurrection.

### 10.7 Replay: the trace re-drives the ordinary loop

The entry point the [§9.5](#95-inbound-the-input-trace) trace exists for:

```julia
trc  = trace(sim)                     # the recorded session: header + per-frame batches
sim2 = Simulation(world)              # the same build
replay!(sim2, trc)                    # header-init, then re-drive every recorded frame
replay!(sim2, trc; to_boundary = k)   # partial: the §13.4 replay-pointer register
```

`replay!` is **the ordinary loop with exactly two substitutions** — not a
separate execution mode, which is what keeps every property proved of the
loop true of replay:

- **Boundary zero from the header.** `replay!` stands in the `init!`
  position of [§10.6](#106-run-lifecycle-and-partial-advance)'s lifecycle: it applies the header's resolved stores
  and slot values directly — no condition resolution; [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)'s totality holds
  by capture — and then executes the ordinary boundary-zero sequence
  ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)). Authored-condition events re-fire identically: the header
  predates the sequence ([§9.5](#95-inbound-the-input-trace)'s capture placement), so nothing is applied
  twice and nothing is skipped.
- **The drain reads the trace.** Each frame top applies the recording's
  batches for that **frame ordinal** instead of swapping the roster's
  staging cells. Ordinal keying is exact because the frame sequence is
  itself deterministic under replay (`t*` boundaries derive from state,
  [§8.4](#84-localization-mechanics)): frame *k* of the replay *is* frame *k* of the recording. Recorded
  batches apply **verbatim, with no surface re-check** — the write-surface
  rule ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)) ran at recording time, and claims are a live-roster fact of
  the recorded session that replay does not reconstruct.

Everything else is the loop as already specified:

- **Termination and partial replay.** Replay ends at the recording's final
  frame, or earlier at `to_boundary = k` — the consumer of [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)'s replay
  pointer, defined as running **through the frame whose execution published
  boundary `k`**, so replay always halts at a frame top. For a grid boundary
  that is exactly at `k` (the frame that publishes one ends at it), and
  [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)'s frame-entry pointer lands the same way; a localized `t*`
  boundary inside the frame is reproduced but not stoppable-at
  ([§8.4](#84-localization-mechanics)'s separation: the trace stays frame-indexed, boundaries are
  the reporting index) — or earlier still under the ordinary policies: `t_end` and
  `stop_on` overrides bind for this replay exactly as at `run!` ([§10.6](#106-run-lifecycle-and-partial-advance)); a
  termination the recorded session hit through `stop_on` reproduces itself
  anyway, deterministically.
- **Replay ends `initialized`, never `stopped`** — boundary-consistent and
  ready to advance, the same state `step!` leaves ([§10.6](#106-run-lifecycle-and-partial-advance)). This is what
  makes three promised workflows real: [§9.2](#92-outbound-snapshot-publication)'s state-trajectory inspection
  ("what was the private state at t = 37.2?" — replay there, read the live
  stores), [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)'s error reproduction (replay to `k − 1`, then `step!` the
  failing frame under instrumentation), and continuation (`run!` after
  `replay!` is a live session from the replayed boundary).
- **Replay re-records.** The trace register runs normally: the new trace
  inherits the old header and accumulates the re-drained batches — a
  bit-identical prefix — so a replayed-then-continued session leaves behind
  a complete, valid trace of *itself*, with no special stitching.
- **Pacing and the control plane are unchanged** ([§8.7](#87-real-time-pacing), [§10.1](#101-control-plane)): pacing sits
  outside the semantics, so paused, slow-motion or real-time replay is free —
  paced replay with an attached visualizer *is* session playback. Stop
  truncates, as anywhere.
- **Devices are readers.** Rostered devices init and spawn normally ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads))
  and consume snapshots ([§9.2](#92-outbound-snapshot-publication), [§10.3](#103-the-next-snapshot-wait)) — the visualizer case — but no live
  staging cell is drained while the trace feeds the loop: a batch found
  staged is discarded with a rate-limited warning
  (`ReplayDiscardedStaging`, [Appendix C](#appendix-c-the-diagnostic-kind-set)). Mixing live writes into a replay
  would destroy the property replay exists to provide; a session that wants
  live input is a continuation (`run!` after replay), not a replay.
- **Validation is loud and up front.** Before the first frame, the header is
  validated against the `Build` (store layout, slot faces) and the trace's
  batch entries against the root input-face list — attach-style, with
  did-you-mean (`ReplayHeaderMismatch`, `ReplayUnknownFace`, [Appendix C](#appendix-c-the-diagnostic-kind-set)).
  The same pass pays the trace-record conversion in reverse: every writer's
  sparse records ([§9.5](#95-inbound-the-input-trace)) are normalized to positional batches
  against the header's schemas, once, off the loop — replay has the whole
  trace in hand before frame 1 — so the replay drain applies compiled
  scatters exactly as the live drain does, and no face name is resolved
  per frame under replay either.
  *Structural* mismatch is an error; *parametric* difference is not:
  replaying against the same structure with changed parameters is the
  **what-if register** — deterministic re-driving of the recorded inputs
  through a modified model. Bit-identity is promised only against the
  identical build; the what-if register promises determinism, never
  reproduction.
  The header's deployment block ([§9.5](#95-inbound-the-input-trace)) validates in the same pass, on
  the *structural* side of that line: the six trajectory-determining
  parameters — `Δt_base`, `h`, `n`, the algorithm, `localization_tol` and
  `event_budget` ([§8.4](#84-localization-mechanics)) —
  are compared against the target `Simulation`'s own deployment binding —
  mismatch is `ReplayHeaderMismatch` with a deployment-parameter
  discriminator, never a what-if, because a deployment change moves the
  times at which the frame-ordinal batches apply: different inputs, not a
  modified model. The localization pair is compared for exactly the same
  reason the grid parameters are: it moves the trajectory, so a run that
  differs in it is not re-driving the recorded one. `t₀` is *applied*, not
  compared — replay stands in the
  `init!` position and owns the anchor, so `replay!` takes no `t0`
  argument — and the header's `t_end`/`stop_on` pair is a recorded fact of
  the recorded session, never a constraint on this one: overrides bind as
  stated above.

Rejected shapes, for the record: a `run!(sim; replay = trc)` flag (replay
replaces `init!` and swaps the drain source — it is a lifecycle *entry*, not
a run option, and folding it into `run!` muddies [§10.6](#106-run-lifecycle-and-partial-advance)'s mandatory-`init!`
rule); a synthetic playback device staging the recorded batches (wall-clock
staging cannot hit recorded frame ordinals — it would reintroduce exactly the
scheduler-determined input timing that [§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads) indicts and replay exists to
remove); replay ending `stopped` (kills all three workflows above).

---

# Part III — Authoring and build

## 11. The declaration layer: components and assemblies

How an author spells a component: where the structural facts live, what the build
takes as authoritative, and what is checked against what. [§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)–[§11.4](#114-failure-walkthroughs-the-error-locality-grounding) cover the
component side; [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)–[§11.8](#118-computed-connections-and-generic-boundaries) the assembly side; the build pipeline is [§12](#12-the-build-pipeline) and the
stopped-sim service spellings are [§14](#14-stopped-sim-services). Concrete syntax below is near-final in shape
but still illustrative in spelling. The sketches (`sketch_decoder.jl`,
`sketch_io.jl`) are written against this layer and the services spellings.

### 11.1 Position: a declarative trait layer — plain Julia, no macros

Stage functions are ordinary multiple-dispatch methods (the `GUI.draw!` precedent).
Structural facts are declared through a small set of well-known functions returning
plain values, defined alongside the methods. No macro DSL: generated code is opaque
to the debugging, tooling and comprehension workflows the charter protects ([§1](#1-purpose-and-method)), and
a macro can only ever *lower to* a layer like this one — so a convenience macro
remains addable a posteriori as pure sugar (the `@kwdef` precedent), while never
becoming load-bearing. Redundancy between declarations and function bodies is
accepted deliberately, under one non-negotiable condition: **every inconsistency
fails loudly**, at build time where possible, at first execution otherwise.
The door stays open for the declaration layer specifically: a macro generating
the well-known declarations — the `where {T <: Real}` ceremony of a continuous
`output_types` ([§11.2](#112-the-declaration-inventory)) being the obvious candidate — is admissible sugar *on
top of* the plain-Julia forms, never a replacement for them and never required
to author a component (row 166). Every rule in this part is stated over the
generated methods, so a macro that lowers to them adds convenience and no
semantics.

**Namespace and extension.** The declarations and stage functions are
**extended, not called**: authoring a component means adding methods to
framework-owned generic functions, and Julia only does that through an explicit
per-name `import` (or a qualified `Flight.f(…) = …` definition, the
`Base.show` idiom [§16](#16-open-axes) records for the extension-only periphery surface). A
component module therefore opens with

```julia
import Flight: init_x, init_m, workspace, input_types, output_types,
    events, h_x, h_xu, f, g, project, child_connections,
    input_connections, output_connections, rates
```

because `using Flight` alone is a silent trap: `f(eng::Engine, …) = …` after a
bare `using` defines a new, unrelated `MyModule.f` — no error, no warning (the
short names are deliberately unexported, precisely because `f`, `g`, `events`,
`project` are the most collision-prone identifiers in numerical Julia, so
there is no clash for the language to detect) — and the build then sees a
component with no `f` method and reports a *modeling* diagnostic
(`StoreWithoutUpdate`; `ClassUnreadable` when the whole inventory was shadowed)
for a one-line namespace mistake, pointed away from the wrong line — the [§11.4](#114-failure-walkthroughs-the-error-locality-grounding)
error-locality inversion arriving through the namespace. Two mitigations,
both normative: the import list above is authoring surface, stated wherever a
component file is first shown; and the two diagnostics run a **shadowing
check** — if the component's parent module defines a same-named function
distinct from the framework's, the message says so and names the missing
import ("`MyEngine`'s module defines its own `f`, distinct from `Flight.f` —
add `import Flight: f`"; a two-line `isdefined`/`!==` test on names the build
already looks up). A convenience macro expanding to the import list remains
addable-a-posteriori sugar per this section's macro doctrine; a re-export
submodule is *not* an alternative — `using Flight.Declarations` would carry
the identical silent-shadowing semantics, per-name `import` being the only
extension register the language provides.

**The same trap has a local-scope sibling** (row 164): written inside a `let`,
a function body or a `@testset`, `h_x(::MyComp, (; x)) = …` does not add a
method to the global `h_x` — it binds a *new local function* of that name.
Calls within the block resolve to it and look correct; the generic function the
build dispatches on never learns of the component, which therefore reads as one
declaring nothing at all. The shadowing check above cannot reach this case —
there is no parent-module binding to compare, the shadow being a local binding
that disappears with its block — so the mitigation is at the other end: **a
component that declares nothing and defines no stage is rejected at build
time**, an inert component being unwritable on purpose. That check costs a line
and catches the misspelled-declaration family with it. Test code is the
realistic victim (a fixture component defined inside its own `@testset`), and
the authoring rule is one line: declarations live at module top level. The net
under a *partially* shadowed component holds because `output_types` is still a
declaration: a component whose ports are declared but whose stage went to a
local binding reads as "declared but not produced" ([§11.3](#113-visibility-the-contract-is-the-interface)), the shadowing
note attached, rather than as a component with nothing to say.

**The schema-authority principle.** Declarations *define* the model's structure;
evaluation *checks* conformance against them — never the reverse. The build probes
user functions with real values (no reliance on compiler inference), compares
observed against declared, and the same comparison runs on every subsequent
evaluation for free (a `NamedTuple`-type check that constant-folds away when
conformant). The rejected alternative — inference-by-evaluation as schema
authority — fails on three counts, established by walkthrough ([§11.4](#114-failure-walkthroughs-the-error-locality-grounding)): error
*locality* inverts (failures surface inside correct code, pointed away from the
wrong line); observed schemas are sample- and branch-dependent (the probe sees only
the initial state's branch — the [§5.6](#56-diagnostics-feedthrough-tracing) hazard corrupting the schema instead of a
diagnostic); and annotations have nowhere to live. Types by declaration, values by
execution, conformance by comparison.

**Contracts are functions of the type, not of the instance.** A leaf's contract
declarations — `input_types`, `output_types`, `events`, and the
shapes of `init_x`/`init_m` — must be determined by the component's
**type**, its type parameters included, and never by its field *values*. The
value-discarding signature (`input_types(::Engine, ::Type{T})`) is the visible form of the
rule, and the idiom for a contract that genuinely varies is the type parameter,
not the field: `SumJunction{Wrench, 3}` ([§6.2](#62-aggregation-explicit-summing-junctions)), `Or{N}` ([§13.7](#137-tooling-consequences-provenance-and-the-component-library)) — arity is spelled
in the type, at the price [§6.2](#62-aggregation-explicit-summing-junctions) states openly. The reason is [§12.7](#127-the-compiled-executor)'s entry typing: a
component's bundle is a `NamedTuple` whose key set *is* its contract's, and an
executor entry carries what selects code in type parameters and what is plain
data in fields. A key set derivable only from field values would have to either
climb into the type parameters anyway — multiplying specialization and changing
[§12.7](#127-the-compiled-executor)'s chunking cost model — or sit in fields, dissolving the static typing that
[§5.1](#51-the-scheduling-problem)'s zero runtime graph logic, [§7.5](#75-allocation-policy-a-scoped-invariant)'s allocation invariant and [§12.5](#125-the-always-on-conformance-check)'s
fold-away conformance test all rest on. The build reads each declaration once,
against the concrete instance, so a value-dependent contract does not announce
itself: this is a rule authors keep, not a check the build can run. **`workspace`
is the one exception**, and explicitly so: it is the by-allocation register (row
77) — an allocator the framework *calls*, not a schema it *walks* — and it
legitimately takes sizes from the instance (`workspace(c::KF, ::Type{T})` reading
`c.n`, [§7.3](#73-discrete-state-modes-and-workspace)), because no entry type is derived from it.

### 11.2 The declaration inventory

```julia
struct Engine <: AbstractComponent
    ω_idle::Float64; ω_min::Float64; J::Float64      #parameters: plain struct fields
    ω_rated::Float64                                 #unread here; §14.2's shipped condition uses it
end

#state stores: declared by initial value — types derived, nothing to drift
init_x(::Engine) = (ω = 0.0,)
init_m(::Engine) = (phase = off,)                    # off | starting | running

#input contract: continuous tier ⇒ the T-form; each entry states what may arrive
input_types(::Engine, ::Type{T}) where {T <: Real} =
    (throttle = T, starter = Bool, fuel_available = Bool, M_load = T)

#output contract = the public interface (§11.3); continuous tier ⇒ the T-form, participation per leaf
output_types(::Engine, ::Type{T}) where {T <: Real} = (M_shaft = T, P = T, ω = T)
#ω names a state field no stage produces → auto-published at stage 1 (§5.3)

#stage and update functions destructure their bundle by name (§5.2)
function h_xu(eng::Engine, (; x, m, u))
    M_shaft = m.phase === running ? torque_law(eng, u.throttle, x.ω) : zero(x.ω)
    return (; M_shaft, P = M_shaft * x.ω)
end

f(eng::Engine, (; x, y, u)) = (ω = (y.M_shaft - u.M_load) / eng.J,)

#events: ordered and named — order is load-bearing (§5.3, §8.6); detection policy by the `localize` flag
events(::Engine) = (
    start    = Event(start_guard, start_handler),                               # boundary-detected (default)
    ignition = Event(ignition_guard, ignition_handler),                         # boundary-detected
    flameout = Event(flameout_guard, flameout_handler; localize = true))        # localized
start_guard(::Engine, (; m, u)) =                        #manual trigger: an input (§10.5)
    m.phase === off && u.starter
start_handler(::Engine, _) = (; m = (; phase = starting))          #no `x` key: no reset
ignition_guard(eng::Engine, (; x, m, u)) =                #predicate form
    m.phase === starting && x.ω > eng.ω_idle && u.fuel_available
ignition_handler(::Engine, _) = (; m = (; phase = running))
flameout_guard(eng::Engine, (; x)) = eng.ω_min - x.ω      #continuous form: localizable
flameout_handler(::Engine, _) = (; m = (; phase = off))
```

The inventory, and where each schema fact gets its authority:

#### State, modes, discrete state

- **State, modes, discrete state** (`init_x` on either tier, `init_m`): declaration
  *by initial value* — the type is derived from the
  value, so there is no second artifact to drift and no separate type declaration to
  check. (The workspace is the exception to the register: it is declared *by
  allocation* — `workspace(::C, ::Type{T})` continuous, `workspace(::C)`
  discrete, the method being the allocator — because a workspace is not
  memory and none of the by-value arguments below cover it; [§7.3](#73-discrete-state-modes-and-workspace).) This is the boundary of legitimate derivation: deriving from another
  declaration is sound; deriving from evaluated user code is not. Rejected:
  declaring types here too, `input_types`-style, with [§12.3](#123-probing-and-input-synthesis)'s
  `probe_value` synthesizing the initial values. The declared values are the
  condition substrate's base layer ([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)'s overlays fall back to them leaf
  by leaf, and the compiled store writers bake `merge(defaults, overlay)`),
  so there must be an authored value under every leaf; synthesized initial
  state would cross the probe-value barrier [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator) makes structural (a
  fabricated zero is a fine probe input and a terrible flight condition —
  states no less than slots); and every field where synthesis picks wrong
  (modes, `Ranged` values excluding zero, trim-sensitive states) would need
  an authored default *beside* its type — the per-field two-register
  protocol [§14.2](#142-fragment-composition-locality-without-schema) kills for `initialize` specs, aggravated by types being
  first-class values in Julia (the two registers distinguishable only by
  `isa Type`). The asymmetry against `input_types`/`output_types` is one of kind, not
  style: contracts describe table cells, recomputed from scratch every
  sweep, needing only types; `init_*` describe stores — the model's memory,
  which must have contents before the first sweep can run.
  **These declarations stay one-argument**, and the criterion is the register
  they live in (row 166), stated once here and referred to from below: a
  *by-value* declaration states nominal physics — its *types* walk by rule, [§7.1](#71-continuous-state-structured-immutable-flat-backing)
  forcing every state leaf to follow the activation scalar, so a `T` in the
  signature would record no choice its author could make, and partials enter
  through per-invocation seeding, never through initialization; a *by-type*
  declaration is a function of the activation scalar, which is why
  `input_types` and `output_types` both take it on the continuous tier; a *by-allocation* declaration
  takes the scalar too, `workspace(c, T)` being the standing precedent (row 77).
  The criterion, not uniformity, is the rule: a `T` in a signature means a choice
  was made there.

#### `input_types(::C, ::Type{T})` — continuous; `input_types(::C)` — discrete

- **`input_types`**: a bare `NamedTuple` of types — zero framework vocabulary, no
  wrapper types (the last candidate, `Reduce`, died with reduce-ports, [§6.2](#62-aggregation-explicit-summing-junctions)). On
  **continuous consumers the two-argument form is mandated**, on **discrete
  consumers the plain one**: the same tier mandate `output_types` carries, class
  read off declaration shape and the form fixed by the class ([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)), either
  violation `TierSignatureMismatch`. Entries remain **face bounds, not cell
  types**, and the reading is **permissive** (row 167): an entry states, per leaf,
  what the consumer *allows* to arrive there.
  - **`T`, alone or as a type parameter** (`SVector{3, T}`, `RQuat{T}`) — the leaf
    is **tolerant**: the activation scalar and a frozen `Float64` are both lawful
    arrivals, so a walking producer, a frozen discrete producer and a root slot
    are all admissible behind it and substitution stays intact. This is what a
    promoting consumer writes, and it is the overwhelmingly common case.
  - **`Float64`** — the leaf **demands frozen**: this input must never carry
    partials. That is the **FFI door** — a component whose internals cannot
    propagate `Dual`s (an opaque wrapper, a C table, a hand-rolled solver)
    declares it, its AD-incompatibility becomes schema-visible instead of
    folklore, and the failure moves from a `MethodError` inside user math at the
    first `Dual` probe to a named wiring error at build ([§6.1](#61-connections-and-hierarchy)).
  - **`Int`/`Bool`/enum leaves and abstract reference-typed entries** as they
    always were. Abstract entries state **structural substitutability** — several
    concrete producer types admissible behind one stable face, [§4.4](#44-function-valued-signals-environment-access)'s field
    handles being the demonstrated client (`terrain = AbstractTerrainField`) — and
    are spelled without `T`, being references rather than numbers. They are still
    never the tool for eltype genericity: that is exactly what a `T` entry is, a
    promoting consumer writing `SVector{3, T}` rather than an abstract bound.

  (Names-only contracts were rejected — they lose the wiring-time type error and
  standalone checkability.) Inputs are the component's *requirements*: [§6.1](#61-connections-and-hierarchy)'s
  unconnected-input error, over-wiring detection and did-you-mean typo messages
  are only definable against them.

  **Two clauses check a wire** ([§6.1](#61-connections-and-hierarchy)). The **nominal bound check** is unchanged in
  force, now stated over evaluations: the producer's declaration at `Float64` must
  be `<:` the entry at `Float64` — one uniform rule, degenerating to exact
  equality for a concrete entry, concrete types being final. Beside it sits the
  **tier-scoped walk-compatibility clause**: for a *continuous* consumer, a
  walking producer leaf (the producer declared `T` there) requires a `T` entry,
  while a pinned producer leaf satisfies either, frozen values embedding upward.
  Both sides are declaration functions of `T`, so the clause is decidable in
  Stratum A by evaluating them at a marker scalar — no user stage code runs
  ([§12.1](#121-three-strata)) — and a violation is `WalkingFaceAtFrozenEntry`. **Discrete consumers
  take the bound check only**, and that scope is load-bearing rather than tidy: a
  discrete stage reads exclusively at real ticks in the nominal world, a
  `Dual`-carrying cell existing only inside activations discrete stages never run
  in ([§12.4](#124-activations-executable-sets-laziness-caching)), so continuous → discrete wires are unconditionally legal — unscoped,
  the clause would reject the entire sensor → controller pattern.

  Because entries are bounds, nothing is ever "overwritten": cell types are
  single-sourced from the producer side per activation ([§12.1](#121-three-strata)), and a
  `Dual`-carrying cell behind a `T` entry is the design working, not a promise
  broken. The code-level complement is the **genericity obligation** — whatever
  scalars the wiring delivers, the consumer's math promotes — still checked by the
  `Dual` probe and never declared, now **scoped to the `T`-entries**: a
  `Float64`-entry input imposes no such obligation, that being its point. So
  **declarations record choices; obligations are checked** — the `T` entry records
  the tolerance choice, the probe checks the promotion.

  **The two readings that stay rejected, and the one that escapes them** (rows 33,
  167). The *predictive* reading — the entry saying what *will* arrive — is
  impossible outright: an input's activation type depends on who feeds it (a
  continuous producer delivers `Dual` under a `Dual` activation, a gated-off
  discrete producer `Float64`), so a predictive declaration would force the
  consumer to state its producer's tier and break substitution behind the same
  face. The *envelope* reading — the entry as a promise to promote — is a
  universal obligation every component owes anyway, hence a constant function
  across components and zero information. The permissive reading is neither: it
  predicts nothing, and it is not constant, pinned entries being rare but real.
  That third reading is what the original adjudication never had on the table, and
  it is what makes the `T` carry information here.

  **Root slots are the one place an entry types a cell**: produced by no
  component, a slot has only the consumer declaration to take a type from. The
  **slot type** is the entry evaluated at `Float64`, and only a *tight* bound
  determines one — a face surfacing as a root slot must resolve to a concrete
  declaration (staging cells, the trace header and `probe_value` all need it;
  abstract-at-root is a build error, and `AbstractAtRoot` names the remedy with
  the face: wire a concrete producer — in a test rig, a stub child, [§13.7](#137-tooling-consequences-provenance-and-the-component-library)). Under
  fan-out the slot type is the unique concrete declaration among its consumers —
  two different concrete declarations remain an error — and abstract co-consumers
  are checked against it. The **slot cells** at an activation follow it by
  evaluating that same entry at the activation's `T`, which makes **seedability
  schema-visible**: a `T`-entry slot is a lawful linearization `B`-matrix tap, and
  a `Float64`-entry slot is *declaredly* unseedable ([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)).

  **Fan-out combines tolerance by a meet, not by agreement** (row 168). Two
  consumers of one slot may agree at nominal and still differ in tolerance —
  `SVector{3, T}` and `SVector{3, Float64}` both evaluate to
  `SVector{3, Float64}`, so the slot *type* is unambiguous while the entries
  disagree about partials. That mixture is a legitimate model rather than a
  mistake: a command consumed by a promoting aerodynamics leaf and by an
  AD-opaque table is the FFI door in use, and the only remedies an error could
  offer are duplicating the slot or lying in a declaration. The slot therefore
  takes the **meet** — pinned at every activation if *any* consumer's entry
  pins, following the scalar only when every consumer tolerates. The direction
  is forced by embedding: a pinned slot cell feeds a `T` entry lawfully (frozen
  values embed upward as zero-partial constants, [§12.5](#125-the-always-on-conformance-check)), whereas a
  `Dual`-carrying cell arriving at a `Float64` entry is precisely what that
  entry forbids — so the meet is the only assignment satisfying every consumer
  at once, the mirror of [§6.1](#61-connections-and-hierarchy)'s walk-compatibility clause on the producer side.
  What the mixture costs is stated where it is paid: such a slot is unseedable,
  and a tap selecting it is rejected naming the *pinning consumer* rather than
  the face alone ([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)).

#### `output_types(::C, ::Type{T})` — continuous; `output_types(::C)` — discrete

- **`output_types`**: the public port contract, declared **by type** — same
  species as `input_types`, and carrying the activation scalar in its signature on
  the same terms, but read **literally** where the input side is read permissively:
  an entry states what the cell *carries*, not what it tolerates. On
  **continuous producers the two-argument form is
  mandated**: `output_types(::Engine, ::Type{T}) where {T <: Real} =
  (M_shaft = T, P = T, ω = T)`. On **discrete producers the plain form is
  mandated**, and it *is* the wholesale pinning of the discrete exemption
  ([§7.2](#72-numeric-genericity-eltype)) — now spelled in the signature as well as enforced by tier: class is
  read off declaration shape ([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)), the class fixes the form the declaration
  must take, and `TierSignatureMismatch` names a producer whose form and tier
  disagree in either direction. Semantics are **literal**: the cell types at an
  activation are the declaration *evaluated* at that activation's `T` — nothing
  walked, nothing inferred — so participation is authored **per leaf** and legible
  on the page:
  - **`T`, alone or as a type parameter** (`SVector{3, T}`, `RQuat{T}`,
    `MyStruct{T}`) — the leaf **participates**: its cell carries the activation
    scalar. Value parameters are structure rather than number and never take it
    (`Ranged{T, -1, 1}`; the bounds are not scalars to re-type).
  - **`Float64`** — the leaf is **deliberately pinned**, and the pin is
    schema-visible: whole-leaf freezing, declared and conformance-checked, which
    is [§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)'s recorded freeze door delivered (declare `Float64` and strip with
    `ForwardDiff.value` inside the stage — the stop-gradient stated in the
    contract instead of buried mid-expression).
  - **`Int`/`Bool`/enum leaves and reference-typed fields** pin as they always
    did (a [§4.4](#44-function-valued-signals-environment-access) bulk-data handle's grid is frozen build-time data, never
    activation-dependent).

  The companion obligation is **constructibility at `T`** — a declared type must
  be buildable at the activation scalar — and the `Dual` probe enforces it by
  construction: it builds real values, so a type whose constructor cannot accept
  them detonates at the probe with its own name in the message.
  During a generic sweep, gated-off discrete producers hold their `Float64`
  values, consumers gather mixed tuples, and promotion does the rest —
  semantically exact, since a frozen discrete output is a constant with zero
  partials, which is precisely what "linearize the continuous dynamics with the
  discrete state held" means (the frozen cell is not an AD limitation on the signal path; it is
  the true zero of an instantaneous dependence the hybrid semantics never had —
  `frozen_discrete_walkthrough.md`). What makes mixing safe — the standing
  objection being that under-the-hood substitution cannot distinguish honest
  `Float64`s from eltype `Float64`s (row 33) — is the **embedding guarantee**
  ([§12.5](#125-the-always-on-conformance-check)), now keyed on **declared-`T` leaves**: a `Float64` observed at a
  declared-`T` leaf under a non-nominal activation implies no `Dual` entered its
  computation (promotion is airtight; there is no lossy cast), so its true
  derivative along every seeded direction is zero and embedding it as a
  zero-partial constant is exact. Piecewise branches returning literal constants
  (`flow > 0 ? f(x) : 0.0`) are legal as written — zero partials is the
  derivative of a locally-constant branch — and which *invocation* carries
  partials is still chosen by seeding ([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)), never by typing: the declaration
  says which leaves *can* carry them, the seed which directions do.

  **The forgotten-`T` account, stated openly.** The whole-signature variant — a
  continuous producer declared as though it were discrete, the nastiest member
  of the class while the plain form was itself the tier marker — is extinct by
  construction: the tier mandate catches it in Stratum A, before any user code
  runs. What remains is per-leaf: an author writes `Float64` at a leaf that
  really participates. That bug **lurks, but is never silent** — no lossy
  `Dual → Float64` cast exists, so the first `Dual` activation of that component
  fails at build, and the message carries the didactic hint ("if `F`
  participates in differentiation, declare it `T`"), an observed `Dual` at a
  declared-pinned leaf having exactly one honest cause. The lurk is contained by
  policy rather than machinery: **the test suite builds a `Dual` activation of
  every component** ([§12.4](#124-activations-executable-sets-laziness-caching)'s exhaustive set — an activation is a Stratum-C
  re-run, cheap enough to make this unremarkable in CI).
  What the form buys in exchange is **reader honesty**: participation is read off
  the declaration instead of reconstructed from a framework rule carried in the
  reader's head, and a genuinely frozen leaf can say so. (History: this signature
  was the original design, abandoned in favor of a plain nominal declaration plus
  a framework leaf walk — rows 33 and 79 — and revived by row 166, whose grounds
  are recorded there: reader-honesty revalued, and two of row 79's three grounds
  dissolved independently in the meantime — class-by-declaration-shape ([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape))
  removed the tier-marker trap, and [§12.5](#125-the-always-on-conformance-check)'s embed-accept removed the
  constant-branch detonation that made the `T`-form look fragile.)

  **The stores keep their walk**; only the output side is evaluated. The type
  derived from a *continuous* leaf's `init_x` is walked exactly as a continuous
  producer's cells once
  were — real leaves and `Real` type parameters follow the activation scalar —
  while `init_m` and a *discrete* leaf's `init_x` pin wholesale, mirroring the
  discrete-producer rule. The asymmetry is the register criterion stated above under the by-value
  declarations, not an inconsistency: `init_*` declare *by value*, and [§7.1](#71-continuous-state-structured-immutable-flat-backing)
  admits no pinned state leaf for a `T` to record a choice about. Declared
  `Float64` initial values embed as zero-partial constants under non-nominal
  activations, which is [§14.3](#143-resolution-flatten-validate-compile-once)'s rule for `Float64` condition leaves applied to
  the defaults those conditions overlay.
  A continuous leaf's `init_x` walk presupposes [§7.1](#71-continuous-state-structured-immutable-flat-backing)'s closed leaf vocabulary (scalars and
  `SArray`s at the common eltype — the discrete tier's `init_x` keeps [§7.3](#73-discrete-state-modes-and-workspace)'s
  full type freedom), so Stratum A checks it ([§12.1](#121-three-strata)) and reports
  a failure in the didactic register: "`init_x` field `gear_count::Int` is not
  a continuous state — integers, `Bool`s and enums belong in `init_m`";
  "`init_x` field `q_nb::RQuat` is not a state leaf — declare the
  `SVector{4}` backing and cast where rotation semantics are wanted ([§7.1](#71-continuous-state-structured-immutable-flat-backing))".

#### `events(::C)`

- **`events(::C)`**: an ordered, named collection of guard/handler pairs, detection
  policy set per event by the `localize` flag — `Event(guard, handler; localize = true)`
  is localized, and the default `false` is boundary-detected ([§2.1](#21-events-two-detection-policies)). Order is semantics ([§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries) declaration
  order, [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) once-per-event); nothing here is inferrable.

#### No stage tags anywhere

- **No stage tags anywhere.** Which stage produces which port stays invisible in
  the contract, preserving [§4.2](#42-consumers-see-ports-not-stages) (moving a port between stages is non-breaking
  for consumers).
  Membership is *derived* with no chicken-and-egg: stage-1 functions
  (`h_x`) structurally receive no inputs, so the build probes them first,
  observes their contract ports, assigns the remainder to stage 2, builds the
  graph, and probes the stage-2 chain in topological order with real upstream
  values. The "decoder takes no inputs" property is exactly what makes
  the derivation well-founded. Stage names carry no tier information at all
  (row 173): what declares a leaf's tier is the completeness rule below.

#### Custom structs as port types

- **Custom structs as port types** (`contact = GearContact{T}`) are
  first-class under [§7.2](#72-numeric-genericity-eltype)'s scoping: parametric in their
  real-scalar leaves, constructors inferring the scalar, no pinned fields on
  the continuous path. A participating struct leaf is declared with the scalar
  in its parameter position (`GearContact{T}`, recursively for nested
  parameters); a struct with a hardcoded `Float64` field offers no such
  position, so it can only be declared bare — a pinned leaf, honestly spelled —
  and any `Dual`-carrying construction detonates inside the stage
  with an `InexactError` naming the offending constructor — the [§7.2](#72-numeric-genericity-eltype) CI
  invariant reached through the declaration layer with no extra machinery.

#### Completeness of the declaration set

- **Completeness of the declaration set** — four rules the build checks in
  Stratum A ([§12.1](#121-three-strata)), stated here because they are properties of the
  declarations, not of the wiring. **A store needs its update:** `init_x` with
  neither an `f` nor a `g` method is a build error — continuous state with no
  flow, or a discrete store nothing updates. The
  framework will not silently supply `ẋ = 0`, which is a model, not a default;
  and an unupdated discrete store
  is a parameter in disguise, parameters being plain struct fields — the
  didactic register says exactly that. `init_m` carries no such obligation:
  modes are written by handlers, and a component may legitimately declare modes
  no event of its own transitions. **An event needs both halves:** an `events`
  entry whose guard or handler has no method for the component type is a build
  error, caught by method lookup at declaration-reading time rather than as a
  `MethodError` at the first firing — an event that fires only in a corner of the
  envelope would otherwise hide the omission indefinitely. **Tier is declared by
  the update law:** for a **stateful** leaf, `f` marks continuous and `g` marks
  discrete — the flow/jump pair carries the tier now that the stage names are
  shared (row 173) — and the remaining tier-implying declarations must agree:
  `init_m` and `events` are continuous-only (the event system is
  continuous-side only: [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), [§3.2](#32-periodic-discrete-component), [§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)),
  `workspace`'s arity splits the tiers (`(::C, ::Type{T})` versus `(::C)`), and
  so do `output_types`' and `input_types`' arities (rows 166–167). Disagreement is
  `DeclarationOnWrongTier` ([Appendix C](#appendix-c-the-diagnostic-kind-set) — the offending
  declaration, with the tier the leaf's other declarations announce), surviving
  on a reduced member set with its job unchanged: declaring both `f` and `g`, or
  a `g` beside a two-argument `output_types`. The wrong-letter class — an `h_z`
  written on a continuous leaf — dies with the distinction it policed.
  A **stateless** leaf declares no store and no update law, so its tier is
  decided by its contract arities — `output_types`, mandatory hence always the
  decider, with `input_types`, where declared, agreeing; and the arity is no
  mere marker — it *is* the tier's semantics (rows 166–167): the two-argument
  forms declare cells and tolerances at the activation scalar, walking with it,
  where the plain forms declare the pinned discrete world. The tier-transparency point stands unchanged: a
  stateless `h_xu` component is tier-transparent library material
  ([§13.7](#137-tooling-consequences-provenance-and-the-component-library)). Members of both families, or of neither,
  are the [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape) class errors. **The root of a build is an assembly:** root slots
  are the root's input faces declared through `input_connections` ([§6.1](#61-connections-and-hierarchy), [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)), and only
  assemblies declare boundary connections ([§11.6](#116-paths-wiring-and-faces)), so a primitive root has no root
  slots — its faces are just its own port names — and every input it declares
  is an unconnected-input error. Exercising a leaf alone is what [§13.7](#137-tooling-consequences-provenance-and-the-component-library)'s
  component test rig is for — it supplies the one-child assembly.

### 11.3 Visibility: the contract is the interface

**Declared in `output_types` = public; returned in `w` = private *by
construction*; returned in `y` and declared nowhere = build error; absent
`output_types()` = no outputs** — visibility by *where the value goes*, the
same move as class-by-declaration-shape. Ports in the contract are connectable,
GUI-listed, snapshot-carried and log-exported; the table is public throughout,
every cell a declared port or an auto-published one, so nothing anywhere needs
a presentation filter. Private intermediates are not filtered, they are simply
not there: `w` is handed from stage to consumer as a value ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)'s one-hop
law), with no cell, no name in a contract and nothing for a wire, a listing or
a log to reach. The inspection path for an intermediate is **promotion**: one
line in `output_types` makes it public, checked and visible everywhere at once
(row 165 — FlightCore's precedent, where an intermediate was inspected by
putting it in the `Model` output and no other way). Publicity is
never implicit: even the minimal component writes
`output_types(::LowPassFilter, ::Type{T}) where {T <: Real} = (x = T,)`, one
line, in exchange for "public" always meaning someone wrote it down.

- **Conformance**: a declared port must be produced — by exactly one
  stage, or by **auto-publication** for declared names matching
  state or mode fields that no stage produces ([§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)). Stage membership is
  derived over `output_types` alone ([§12.1](#121-three-strata)). Declared-but-unproduced
  and produced-by-two-stages are build errors; a declared port matching neither
  a stage product nor a state
  field errors with both lists in hand ("not produced by any stage and not a
  state field"). A *returned port field* declared nowhere is a build error at probe
  with did-you-mean against `output_types` (under an observation-based
  rule it would silently define a new cell instead — the
  return-side analogue of [§11.4](#114-failure-walkthroughs-the-error-locality-grounding) walkthrough 1). The forgotten-branch
  walkthrough holds: a declared `P` missing from the taken branch's
  return fails at probe; missing from an *untaken* branch, it fails loudly at
  that branch's first execution via the always-on check.
- **`w`'s regime is probe-observed, and that is sound precisely because `w` is
  not a cell.** A cell needs a fixed type per activation, which only a
  declaration can supply; a value flowing between two functions in one fused
  pass has no type contract to violate, and mixed branches are handled exactly
  by promotion. So the probe takes `w` as it finds it: it checks that the
  second return slot is a `NamedTuple` at all, and it checks the *consumer's*
  reads against the observed field set, a destructured name that is not there
  failing inside [§13.2](#132-diagnostics-structured-values-one-carrier-exception)'s framing diagnostic with did-you-mean from the
  actual fields. That is weaker than a declaration-backed message — it can say
  "`f` of `Foo` reads `w.q_dny`; the producing stage returned `q_dyn`" but
  cannot say which spelling was intended — and it is located, name-shaped and
  costs no declaration.
- **Branch-shape rule**: stage returns must have the same `NamedTuple` shape on
  every branch — which Julia's type-stability discipline already demands for
  performance; the framework merely makes it a stated rule with a good error.
  `w` is inside the rule: a stage whose `w` changes shape between branches is
  as broken as one whose `y` does.
- **The always-on check covers `w` at the nominal activation only.** Beside the
  `y` test ([§12.5](#125-the-always-on-conformance-check)) sits one baked `isa` against the type the *nominal probe
  observed* — folding to nothing while the stage is stable, and converting the
  unintended-branch-divergence class, which otherwise announces itself only as
  an allocation in somebody's benchmark, into a loud located field-naming error
  at the divergent branch's first execution. The blame text says what it
  honestly knows: "expected what the nominal probe observed". This complements
  [§7.5](#75-allocation-policy-a-scoped-invariant)'s canary — the canary detects, this localizes. **Non-nominal
  activations run no `w` check at all**, and deliberately: there is no
  declaration to anchor a branch-independent expectation, no store whose typing
  the check would be protecting, and a strict probe-observed expectation would
  fire on the legal constant-branch idiom (the one-probe-point argument, which
  is exactly what kills observation authority for cells). Correctness needs no
  guard there: a `Float64` in `w` under a `Dual` activation is an honest
  zero-partial constant by the embedding guarantee ([§12.5](#125-the-always-on-conformance-check)), and its
  downstream promotion is exact. Walking the nominal observation to synthesize
  non-nominal expectations was rejected as machinery kept alive for a check
  that catches nothing the nominal one misses.
- **Schema authority is total over the table**: every *cell* traces to an
  authored declaration, the always-on check's expected type for `y` is fully
  declaration-derived, and return typos cannot silently define new cells.
  Protection against silently dropped partials rests on the
  embedding guarantee — promotion is airtight, so an observed `Float64` is a true
  constant, [§12.5](#125-the-always-on-conformance-check). Probe-observed expected types remain rejected *for
  cells* on two grounds — authority inversion, and the fact that one probe
  point cannot speak for branch-dependent types — and neither ground reaches
  `w`, which declares nothing and types nothing.
- **What this rules out**: the `unlisted` flag ([§4.2](#42-consumers-see-ports-not-stages)) — presentational
  hiding of connectable ports — and its satellite-function representation; the
  RNG-state case that motivated it needs *nothing* here (`g` reads `x` directly,
  [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)). Identity publication by default goes with it ([§7.4](#74-the-fused-evaluation-lineage-prior-art-and-how-we-got-here) step 4): publication
  driven by the contract replaces publication of everything with hiding
  annotations on top. **Probe-observed private cells** and the
  `Private(T)` fallback are rejected alongside — the former by the
  authority-inversion argument above, the latter as obviated by `w` (a
  wrapper inside `output_types` would break "declared = public" and introduce the
  layer's first wrapper type, where the return channel encodes privacy by
  where the value travels, with no declaration to write at all). So is the
  opt-in variant with a `Float64`-under-`Dual` diagnostic —
  it legislates an ambiguity that strictness dissolves. Rows 34, 55 and 165.

### 11.4 Failure walkthroughs (the error-locality grounding)

The five mistakes that decided declaration-vs-inference, with their failure sites
under this layer (each was traced under inference-by-evaluation too; in every case
the failure surfaced inside *correct* code, later, or never):

1. **Typo'd wire** (`:throtle`): build error at the connection, "no input
   `throtle`; did you mean `throttle`?" — vs. a missing-field error inside a
   correct `h_xu` at probe time, with the input set silently *defined* by the typo.
2. **Forgotten wire** (`fuel_available`, read only by a guard): [§6.1](#61-connections-and-hierarchy) unconnected-input
   error at build — vs. detection contingent on the probe exercising every guard,
   framed as a missing field in event code.
3. **Forgotten branch field** (`P` returned by one branch only): probe or
   first-execution error naming the declared port — vs. a schema silently derived
   from whichever branch the initial state took, then a mid-run error (or a
   silently absent port) at the first transition.
4. **Type mismatch** (a `Float64` fraction wired into a `Bool` input): wiring-time
   error naming both endpoints and both faces — vs. a `MethodError` deep inside
   user math.
5. **Typo'd return field**, in its two registers. A typo'd *port*
   (`P_shft = ...` for a declared `P_shaft`) keeps the full strength of the
   declaration: a probe error with did-you-mean against `output_types`, plus
   the unproduced-`P_shaft` error with both the stage-product and state-field
   lists in hand — vs. the typo silently *defining* a new cell under
   observation-authority, with the intended name's absence surfacing later as a
   missing-field error inside correct `f`/guard code (the return-side twin of
   walkthrough 1). A typo *inside* `w` (`q_dny` where the consumer reads
   `q_dyn`) has no declaration to be checked against and surfaces one hop
   later, at the consumer, as [§13.2](#132-diagnostics-structured-values-one-carrier-exception)'s framing diagnostic carrying the
   producing stage's observed field set ("`f` of `Foo` reads `w.q_dyn`; the
   stage returned `q_dny`") — weaker than declaration-backed, since the
   framework cannot know which of the two spellings was meant, but located at
   the pair of lines that disagree and still name-shaped. That is the price of
   the private channel, paid where no interface is at stake, and the remedy for
   an intermediate worth stronger checking is the one [§11.3](#113-visibility-the-contract-is-the-interface) names: declare it
   an output.

### 11.5 Assembly declaration: type-based, class by declaration shape

An assembly is a plain struct: fields whose type is `<: AbstractComponent` are its
children (field names = path segments), all other fields are inert parameters;
substitutability and variants use ordinary parametric fields — exactly today's
`Cessna172X{K, A}` shape. Alongside it, well-known declarations:
`child_connections(::A)` (mandatory, even when empty), `input_connections(::A)`,
`output_connections(::A)`, `rates(::A)`.

**Container children.** A field whose type is a `Tuple` or
`NamedTuple` with *every* element `<: AbstractComponent` contributes its
elements as children, path-named `"field/1"…"field/N"` (tuples) or
`"field/key"` (NamedTuples), declaration order governing layout. Containers
are **transparent grouping, not assemblies**: no contract, no `child_connections`,
no rate scope, no existence beyond the path segment — the elements are
children *of the parent*, whose `child_connections`/`input_connections`/
`output_connections`/`rates` address them
by element name; anything wanting its own wiring or faces declares itself an
assembly. The payoff is parametric composition: `struct Formation{NT <:
NamedTuple}; aircraft::NT; … end` holds any roster — size, names, mixed
aircraft types — per instantiation, with the declaration bodies generating
wires by comprehension over the keys (the arity-via-computed-contracts
pattern, [§6.2](#62-aggregation-explicit-summing-junctions)'s `SumJunction{W, N}`, at structure scale); [§14.9](#149-mounting-problems-as-relocatable-values)'s swarm worlds
and `at("aircraft/red", problem)` mounting consume it directly. Rules: a
container mixing component and non-component elements is a build error in
this section's did-you-mean family (all-component = children;
zero-component = inert parameter data); containers of containers are
rejected in the first cut (deeper grouping is what assemblies are for);
empty containers are legal (zero children — parametric code needs no special
case); abstract element types follow the same concreteness discipline as
plain fields (directly concrete or via type-parameter bounds, [§11.8](#118-computed-connections-and-generic-boundaries)'s
generic holding). `rates` needs no rule change: element names are immediate
child names, hence legal keys; the bare field name is sugar for a uniform
`K` across all elements.

**The builder is rejected** (`Assembly()` + `add!`/`connect!`): the type you dispatch on and the recipe that defines its structure
become two artifacts with nothing tying them together — [§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)'s drift disease at
assembly scale; it threads mutable state through declaration code; and it does not
even buy source-location capture (a called function cannot see its caller's line
any more than a returned value can). Its one real advantage, programmatic
generation, survives intact in the type-based form: a declaration is an ordinary
function body — loops and comprehensions build the returned tuple.

**No `AbstractAssembly`; one root `AbstractComponent`.** Two facts kill a class
supertype: Julia's single inheritance is already spoken for by the domain
hierarchies (`AbstractAircraft`, engine families — a slot `E <: AbstractEngine`
must accept a primitive `PistonEngine` and a composite turbofan assembly alike),
and [§11.3](#113-visibility-the-contract-is-the-interface) says class is an implementation detail behind the contract — encoding it
in public type identity is exactly what contract visibility exists to prevent.
Class is instead declared by *which* well-known declarations a type defines:
`child_connections` (the marker, mandatory-even-if-empty — the `LowPassFilter`
precedent) makes an **assembly**; any leaf declaration makes a **primitive** —
`init_x`/`init_m`, `workspace`,
`input_types`/`output_types`, `events`, or any stage, `f`, `g` or
`project` method. The rule is total: a `<: AbstractComponent` type declaring
neither family has no class to read, and is a build error naming both families
rather than a silence that fails later and elsewhere.
Reading which declarations exist is reading declarations — the same move as
[§11.3](#113-visibility-the-contract-is-the-interface)'s visibility-by-declaration-site, not [§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)'s banned
inference-by-evaluation. Class also **mandates the shape of the contract
signatures** rather than merely being read from them (rows 166, 167): **both**
contract declarations follow the tier — a continuous leaf's `input_types` and
`output_types` must take the two-argument form
`input_types(::C, ::Type{T}) where {T <: Real}`,
`output_types(::C, ::Type{T}) where {T <: Real}`, a discrete leaf's both the plain
one-argument form, and either violation — a continuous declaration missing the
`T`-form, a discrete declaration carrying one — is `TierSignatureMismatch`
([Appendix C](#appendix-c-the-diagnostic-kind-set)), reported with the component path, the declaration at fault, the tier its
other declarations announce and the form found versus the form mandated. The
check is Stratum A and collected: declaration shape is read, nothing is
evaluated. So the tier fact is spelled in the signature *and* fixed by the
class, the two kept in agreement by a check rather than by convention — which
is what makes the whole-signature forgotten-`T` bug (row 79's worst case, when
the plain form silently *was* the tier marker) unwritable. Error
taxonomy: `child_connections` plus `init_x`/stages on one type is a build error
(assemblies have no state of their own — [§8.5](#85-multi-rate-tick-scheduling)'s no-atomic-assemblies at
declaration time); the neither-family error sharpens into a did-you-mean when the
type has component-typed fields ("holds components but declares no
`child_connections`").

### 11.6 Paths, wiring and faces

**Paths are slash-separated strings** — `"systems/ldg/left/trn"` — relative to the
assembly being declared, no leading slash; one canonical form, shared verbatim by
declarations, error messages, device/trace addressing ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)) and the HDF5 log
tree. Container children ([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)) add index and key segments — `"aircraft/2"`,
`"aircraft/red"` — ordinary segments, resolved against the container field. Rejected: instance navigation (`a.ldg.left` cannot yield a *path* —
symmetric immutable siblings are `===`-identical, so path-from-instance is
unrecoverable by construction; a path-tracking proxy remains addable sugar);
tuples of symbols (structure without readability); dotted paths (a false
Julia-property affordance — the last segment is a contract port, not a field;
slashes say "named tree", which is the true model).

**`child_connections(::A)`** is an ordered collection of `"src/port" => "dst/port"`
pairs, strictly child-port → child-port; [§6.1](#61-connections-and-hierarchy)'s rules apply (one wire per input,
deep paths through concretely-typed fields only, stopping at a generic child's
faces). The assembly's **boundary** is declared by two further methods, one per
direction. **`input_connections(::A)`** is an ordered collection of pairs, face
name => internal endpoint path — or a tuple of paths for an input face routed to
several internal endpoints (fan-out through the boundary)
(`"trn" => ("systems/ldg/left/trn", …)`). **`output_connections(::A)`** runs the
other way, internal source path => face name
(`"aircraft/pose" => "view_pose"`), so that its pairs, like every other pair in
the three declarations, read along the flow. **Face names are arbitrary strings with
two build-checked invariants: no `/` (reserved for structural paths) and
uniqueness across the union of the two boundary declarations' face names** — every other naming choice
(separators, grouping prefixes like `"pilot.throttle_axis"`) is author
convention, not framework law; the `input_passthrough` helper's defaults ([§11.8](#118-computed-connections-and-generic-boundaries)) document the
house style without legislating it. The two-notation rule this rests on is directional —
structure vs. contract, not read vs. write: **slash is structure** (endpoint
paths walking real children and ports; the inspection register's snapshot and
log addressing), **face names are opaque contract tokens** — the periphery's
write side (input devices, mappings, the trace, the GUI write path) speaks
face names exclusively ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)), and the read side speaks them wherever it
wants meaning that outlives the build: integration bindings ([§9.2](#92-outbound-snapshot-publication)'s
`get_face`) and load-bearing service reads ([§14.4](#144-two-application-registers-over-one-plan)).
Pairs-of-strings rather than a NamedTuple also removes the `var"..."` noise that
non-identifier names would force.

One invariant spans all three declarations: every pair's arrow points the way the
signal flows — the left side is a producer or entry point, the right side a
consumer — and every right side is fed exactly once. **Direction is therefore
declared by the method**, not inferred: the resolved endpoints only *cross-check*
it, and an entry whose endpoint resolves to a port of the wrong direction is a
build error naming the method, the entry and the resolved port's actual direction.
A mixed entry is no longer expressible, so that error class disappears with the
single list that made it possible; two entries producing the same output face
remain the ordinary two-producers error. Face *types and tiers* are derived from
the internal endpoints — [§11.2](#112-the-declaration-inventory)'s blessed derivation-from-declarations — and the
derivation is forced, not merely convenient: an assembly is tier-neutral (it
exports continuous-sourced and discrete-sourced ports side by side), so
author-declared face types would need per-face tier annotations restating what
each producer's class already fixes ([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape); a face's cells follow the
producer's own declaration, evaluated at the activation scalar on the
continuous tier and pinned on the discrete). Rejected spellings: routing values under the leaf names
`input_types`/`output_types` (a name-level pun — a discrete leaf's exact signature with
alien value semantics, killing [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)'s name-level class split); leaf-style *typed*
faces plus face wires inside `child_connections` (the tier problem above, plus
face/child namespace collisions and the weakest class marker); routing-as-wires
with derived types and no face list (facehood implicit in wiring — publicity is
never implicit, [§11.3](#113-visibility-the-contract-is-the-interface)).

**Root slots fall out with no vocabulary**: at every non-root level an input face
declared through `input_connections` is fed by the parent's wire; at the root there
is no parent, and the
face *is* the write surface ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)). [§6.1](#61-connections-and-hierarchy)'s whole-tree obligation model states the
complementary error rule. An assembly never declares its external connections —
those live in the parent that instantiates it, exactly as a leaf's do.

**A worked assembly.** [§15.5](#155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary)'s IMU, spelled in full — a mixed-tier assembly
exercising paths, faces and rates together:

```julia
struct IMU <: AbstractComponent
    integrals::IMUIntegrals    # continuous — cumulative Θ, q, Υ, V
    sampler::IMUSampler        # discrete — integrate-and-difference latches
    errors::IMUErrorModel      # discrete — scale/bias/noise on the sample
end

child_connections(::IMU) = (
    "integrals/Θ" => "sampler/Θ", "integrals/q" => "sampler/q",
    "integrals/Υ" => "sampler/Υ", "integrals/V" => "sampler/V",
    "sampler/sample" => "errors/sample",
)

input_connections(imu::IMU) = (
    input_passthrough(imu, "integrals")...,       # kinematic-truth inputs pass through
)

output_connections(::IMU) = (
    "sampler/sample"     => "sample",             # ideal increments
    "errors/sample_meas" => "sample_meas",        # measured increments (the error
                                                  # model's output port)
)

rates(::IMU) = (sampler = 1, errors = 1)
```

Two spellings worth reading closely: `input_passthrough` enumerates the child's
**input** faces and nothing else ([§11.8](#118-computed-connections-and-generic-boundaries)), which is why the pass-through of the
integrals' kinematic-truth inputs (`q_eb`, `r_eb_e`, `ω_eb_b`, `a_ib_b`,
`α_ib_b`, [§15.5](#155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary)) is a bare splat with nothing to say about direction;
and the measured-increment face sources `errors/sample_meas`, the error
model's *output* port, not the `errors/sample` input the sampler already
feeds — listing `errors/sample` in `output_connections` would fail the direction
cross-check, and listing it in `input_connections` while it is wired is the
two-producers error of [§11.8](#118-computed-connections-and-generic-boundaries).

Three facts the example carries: the assembly is tier-neutral — every face's
type and tier derive from its internal endpoint, and a `rates` key on
`integrals`, the continuous child, would be a [§11.7](#117-rate-scopes) build error; the two
discrete children default to `K = 1` anyway, so this `rates` declaration is
declaratory — their absolute rate arrives from the enclosing scope at
deployment ([§11.7](#117-rate-scopes)); and [§15.5](#155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary)'s latch-back wire, where the integrals consume
the sampler's published latch, joins `child_connections` as one more ordinary pair.

### 11.7 Rate scopes

`rates(::A) = (fcs = 1, nav = 5)` — child name => integer multiplier $K \ge 1$,
optional, unlisted children default to 1; [§8.5](#85-multi-rate-tick-scheduling) semantics unchanged (relative,
composing multiplicatively down the tree, compiled to absolute divisors). Keys are
**immediate child names only** — a deep key would edit another type's design from
outside, and the composition rule guarantees you never need to. Container
elements ([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)) are immediate children, so `"aircraft/red"` is a legal key;
the bare field name applies one `K` to every element. `K` on a
continuous child is a build error ([§8.5](#85-multi-rate-tick-scheduling)'s Δt-on-continuous error at declaration
time). `Δt_base`, `h` and `n` appear in no declaration — they are deployment
decisions fixed at `Simulation` construction. Rejected: `K` carried on the child
instance, FlightCore-`Subsampled`-style — it wraps the field type (polluting
paths, dispatch and the child's contract as seen by wiring) and makes a
per-instance value of what [§8.5](#85-multi-rate-tick-scheduling)'s own rationale calls intrinsic to the design,
i.e. a fact about the assembly type.

### 11.8 Computed connections and generic boundaries

`input_connections` and `output_connections` are ordinary functions evaluated at
build against the concrete
instance, so they may *compute* entries from child contracts — derivation from
declarations, [§11.2](#112-the-declaration-inventory)-blessed. The framework helper, sketched:

```julia
function input_passthrough(asm, child_path::AbstractString;
                     prefix::AbstractString = child_path,   # "" → no prefixing
                     sep::AbstractString = ".",
                     except::Tuple = (), only::Tuple = ())  # mutually exclusive

    child = resolve(asm, child_path)      # getfield walk along "/" segments
    names = input_faces(child)            # the leaf's input_types keys,
                                          # entries of input_connections(c) for an assembly
    isempty(except) || isempty(only) ||
        declaration_error(child_path, :both_given)  # exclusivity enforced, not documented
    unknown = setdiff((except..., only...), names)
    isempty(unknown) || declaration_error(child_path, unknown, names)  # list in hand
    wanted = isempty(only) ? setdiff(names, except) : only
    label(n) = isempty(prefix) ? n : string(prefix, sep, n)
    return Tuple(label(n) => string(child_path, "/", n) for n in wanted)
end

input_connections(w::World) = (
    input_passthrough(w, "aircraft"; except = ("atm", "trn"))...,   # "aircraft.pilot.throttle_axis"
    input_passthrough(w, "atmosphere"; prefix = "env", sep = "_")..., # "env_wind_N"
)

output_connections(w::World) = (
    "aircraft/pose" => "view_pose",
)
```

The child is named by path, never passed as an instance ([§11.6](#116-paths-wiring-and-faces)'s `===` problem);
a face name containing dots is a legal final path segment on the internal-endpoint
side precisely because slash is the only structural separator, and computed entries mix
freely with hand-written ones in either declaration. `resolve` and
`input_faces` are build-pipeline primitives needed anyway — `input_passthrough` is a thin
composition, which is what keeps it sugar rather than machinery; no `rename` hook
because the boundary declarations are ordinary code (map over the pairs); normative signatures
for both primitives in [§13.3](#133-build-primitives-resolve-and-the-face-list-accessors). Every error stays
first-class: an `except` face the assembly then fails to wire is an ordinary
unconnected input; a face both wired and passed through is a two-producers error;
`except`/`only` naming a nonexistent face errors with the child's face list in
hand; a `prefix = ""` collision is caught by the build's uniqueness check like any
hand-written duplicate. The effective face list is plain printable data — the
inspectable contract of this instantiation. What computation does *not* do is
auto-bubble: the author wrote down "every input face of this child that I don't
feed, I expose under this prefix" — explicit at the type level, evaluated at
build.

**The name carries the direction.** `input_passthrough` reads
`input_faces(child)` and `except`/`only` filter *face names* within that set;
the helper exists for the pass-through case, where an assembly hands a child's
unfed requirements up one level. Computed *output* re-export — a sibling
`output_passthrough` splatted into `output_connections`, with the predicate
selection [§13.7](#137-tooling-consequences-provenance-and-the-component-library)
records for `except`/`only` — is a **guarded addition**: plausibly wanted by the
conventional-surface work ([§9.2](#92-outbound-snapshot-publication), [§16](#16-open-axes)) and by test rigs, cheap to add, and not
adopted here, because every output face in the worked assemblies is an explicit
pair and no consumer has yet demonstrated the computed form. It could not be a
keyword on one helper in any case: after the boundary split a single call cannot
emit entries into two different declarations.

**One authored list, two declarations.** The `World` example's two-entry
`except` understates the real shape. Every level of a realistic tree is a
generic seam, and an assembly that feeds some of a child's input faces while
passing the rest up must name the fed ones in `except` — at C172X scale, four
seams and roughly ten names at the innermost one, restating in each `except`
tuple the wire list sitting in the same assembly's `child_connections`. That is
[§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)'s
"structure kept in two artifacts" (row 39), the shape this design refuses
elsewhere. It needs no vocabulary: declaration bodies are ordinary code
([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)), so
the author writes the feed list *once* and both declarations compute their
share of it.

```julia
# one authored artifact: actuator output face => destination child input face
const ACT_FEEDS = (
    "e"          => "aero/e",
    "a"          => "aero/a",
    "r"          => "aero/r",
    "brake_left" => "ldg/left.brake",
    …                                            # ~10 entries for the C172X
)

# the face names of `child` the feed list targets
fed_faces(feeds, child) = Tuple(chopprefix(dst, child * "/")
                                for (_, dst) in feeds
                                if startswith(dst, child * "/"))

child_connections(::Systems) = (
    (("act/" * src) => dst for (src, dst) in ACT_FEEDS)...,
    "aero/wrench" => "wr_sum/in1",               # non-feed wires unchanged
    …
)

input_connections(s::Systems) = (
    input_passthrough(s, "aero"; except = fed_faces(ACT_FEEDS, "aero"))...,
    input_passthrough(s, "ldg";  except = fed_faces(ACT_FEEDS, "ldg"))...,
    …
)
```

Adding an actuator channel is then one edit: the new pair simultaneously
creates the wire and removes the face from the input face surface. The two
declarations cannot drift, because neither holds the shared names — both are
projections of the authored list, so the drift class is removed rather than
detected. Every misspelling stays loud: a mistyped destination is an
unknown-face error with the child's face list in hand, whether the wire or
the `except` entry meets it first. One honest asymmetry: a pair *omitted*
from the list is not an error but a structural change — the face leaves the
`except` set and joins the input face surface, ultimately a root slot for
conditions to cover
([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)).
What the idiom preserves, and the helper below surrenders, is that the feed
statement exists to be reviewed: an omission is legible in one authored
artifact, not defined away as the complement of the wire list.

**The line not to cross** is deriving `except` from `child_connections` itself — a
helper spelled `except = fed(s, "aero")`, reading the assembly's own wire
list. That is auto-bubbling under another name (row 43): the author's explicit
statement of which faces are fed would vanish, so a *forgotten* wire would no
longer be a build-time unconnected-input error but a silent promotion of the
face to a live root slot — caught at best later as an `UninitializedSlots`
deployment error of misleading shape
([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)),
at worst not at all, once a GUI or a condition writes it.
[§11.4](#114-failure-walkthroughs-the-error-locality-grounding)'s walkthrough 2,
inverted. The single source must be **authored data, never inferred
structure**.

**Generic holding = imposed contract.** A parent holding a child generically
constrains it exactly through the faces its wires and boundary connections reference: build a
`World` whose concrete aircraft lacks a referenced face and the error names the
`World` entry — build-time structural typing, no new vocabulary (a formal
required-faces declaration on domain abstract types remains possible sugar).
Scalar faces make partial scripting compose: a guidance scenario component wires
`mode_req` and `EAS_ref` while the remaining faces stay exported for GUI or
defaults — impossible with a bundled face ([§4.3](#43-table-mechanics-and-port-granularity) write-side corollary).

---

## 12. The build pipeline

The build consumes a root component instance and produces the runnable artifact:
resolved wires, typed signal table, evaluation schedule, absolute rate divisors,
flat state layout, root slots. [§11](#11-the-declaration-layer-components-and-assemblies) states what is declared and what must hold;
this section states *when* each fact is checked, against what, and with which
failure. The [§11.4](#114-failure-walkthroughs-the-error-locality-grounding) walkthroughs plus [§6.1](#61-connections-and-hierarchy)'s error rules are its acceptance tests.
Error-*reporting* policy is settled in [§13.1](#131-reporting-policy-collect-the-checks-fail-the-evaluations-fast): declarative checking passes
collect, user-code evaluation fails fast, strata are barriers — the only partial
results carried past failures are violation lists from pure checks.

### 12.1 Three strata

Three ordering constraints are forced by settled decisions: face derivation is
**bottom-up** (an assembly's boundary connections evaluate against child contracts,
[§11.8](#118-computed-connections-and-generic-boundaries)); the unconnected-input obligation check and cross-level two-producers
detection are **global** — decidable only at the root, after every assembly's
wires and faces are in hand ([§6.1](#61-connections-and-hierarchy)); and stage membership is **derived by
probing** the stage-1 functions ([§11.2](#112-the-declaration-inventory)), so evaluation interleaves with graph construction at
exactly one blessed spot. The pipeline is therefore inherently heterogeneous,
organized as three strata:

- **Stratum A — structure.** Pure declaration reading; no user stage code
  executes (`input_connections`/`output_connections`/`input_passthrough` bodies are
  declaration code, [§11.8](#118-computed-connections-and-generic-boundaries)). Tree walk
  from the root instance: components by path, class read off
  declaration shape ([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)), leaf contract collection (`input_types`,
  `output_types`, `init_*` values, `events`), bottom-up face
  derivation, then global wiring resolution to absolute leaf terminals:
  one-writer-per-input, typo did-you-mean against the destination's input
  list, the two wiring type clauses ([§6.1](#61-connections-and-hierarchy), [§11.2](#112-the-declaration-inventory)) — the **bound check** at
  nominal faces (the producer's declaration at `Float64` `<:` the entry at
  `Float64`, equality the concrete degenerate; abstract-at-root detected here)
  and, for continuous consumers, the **walk-compatibility clause**, decided by
  evaluating both declarations at a marker scalar and comparing per leaf
  (`WalkingFaceAtFrozenEntry`) — which stays inside this stratum's charter
  because both sides are declaration functions of `T`: declarations are
  evaluated, no user stage code runs —
  the whole-tree obligation
  check, root slots falling out as the root's input faces, and [§7.1](#71-continuous-state-structured-immutable-flat-backing)'s
  closed leaf vocabulary checked on every `init_x` (the walk in [§11.2](#112-the-declaration-inventory) rests on
  it). Also here: [§11.2](#112-the-declaration-inventory)'s declaration-completeness rules (a store without its
  update, an event missing a guard or handler method, a leaf mixing tier
  families, a contract signature whose form contradicts the leaf's tier
  ([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)), a primitive at the root), and
  `rates` validation and compilation of relative multipliers into
  absolute divisors — everything except binding `Δt_base`, which is
  deployment.
- **Stratum B — schedule.** The single evaluation-feeds-structure step:
  workspace allocation at the probing scalar (sound this early: the allocator
  reads only the instance and the scalar, row 77 — no layout dependence),
  stage-1 probes at `Float64` on `init_x`/`init_m` values (well-founded — the
  no-feedthrough stage takes no inputs), port classification over
  `output_types` alone (stage-1 /
  auto-published / stage-2 remainder — a stage's `w` classifies nothing, being
  no part of the contract, [§11.3](#113-visibility-the-contract-is-the-interface)), feedthrough graph from wires carrying
  stage-2 ports, topological order, [§5.5](#55-algebraic-loop-policy-reject-at-build-time) cycle rejection. The output is
  structural: names only, `T`-independent, branch-protected by the
  branch-shape rule plus the always-on check ([§12.5](#125-the-always-on-conformance-check)).
- **Stratum C — activation, parametric in `T`.** Everything type-shaped:
  the producers' output declarations **evaluated** at the activation's `T` to
  type the cells ([§11.2](#112-the-declaration-inventory)'s literal semantics — a continuous producer's
  two-argument declaration called at `T`, a discrete producer's plain one read
  once and pinned), the `init_x`-derived state type walked as before, the probe chain run in topo order —
  threading each stage's `w` to its one-hop consumers ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), [§12.3](#123-probing-and-input-synthesis)) — observed compared against
  declared, flat `x` buffer and table laid out. The nominal `Float64`
  activation runs at build; other activations re-run *only this stratum*
  ([§12.4](#124-activations-executable-sets-laziness-caching)).

Deployment binding (`Δt_base`, `h`, `n`, `t_end`, algorithm,
`localization_tol`, `event_budget`, harmonic-grid
validation, tick schedule instantiation) sits after all three, at `Simulation`
construction — nothing in A–C depends on it. Its validation is collected like
its declarative siblings ([§13.1](#131-reporting-policy-collect-the-checks-fail-the-evaluations-fast)): a nonpositive `h`, an `n < 1`, a
harmonic-grid violation, an algorithm the stepper seam does not know, a
nonpositive `localization_tol`, an `event_budget` that is not an integer ≥ 1
are collected and reported as `DeploymentInvalid` ([Appendix C](#appendix-c-the-diagnostic-kind-set) — parameter,
value, the violated constraint). The localization pair validates on its own
terms only — grid-independent, so it takes no part in the harmonic-grid
check ([§8.4](#84-localization-mechanics)).

### 12.2 The `Build` artifact

`build(world) → Build` is a standalone entry point; `Simulation(world; ...)`
([§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)'s spelling) is the convenience that calls it and adds
deployment binding, buffers and the stopped-sim services. The constructor is
two entry points, not one: `Simulation(build::Build; ...)` accepts the
artifact directly, and `Simulation(world; ...)` is *defined as*
`Simulation(build(world); ...)` — the build CI checked, an acceptance test
targeted or a face-provenance table was printed from is the one deployed,
never an assumed-equal reconstruction (computed boundary-connection bodies are ordinary
user code re-evaluated on every build, so equality between two builds of the
same world is an assumption the factorization removes). Deployment binding
still happens only at `Simulation` construction, whichever entry point runs.
**The `Build` is immutable and may back any number of `Simulation`s,
concurrently** — true by construction once buffers are single-owner ([§12.4](#124-activations-executable-sets-laziness-caching)): each
`Simulation` materializes its own from the shared layouts, so nothing writable
is shared. The one mutable thing on the artifact is the lazily populated
activation cache, whose insertion [§12.4](#124-activations-executable-sets-laziness-caching) makes torn-state-free.
The `Build` is the
inspectable contract of the instantiation [§11.8](#118-computed-connections-and-generic-boundaries) gestures at — wire list, face
table, schedule, root slots as plain printable data. CI checks a model by
calling `build`; the acceptance tests target `build` errors directly;
`attach!` validates device bindings against it. Rejected: build living only
inside the `Simulation` constructor — CI would construct simulations with
dummy deployment parameters, and the phase outputs must exist as coherent
intermediate data anyway; the artifact just names them.

### 12.3 Probing and input synthesis

**Probe-everything scope.** The nominal activation probes every user function —
stages, `f`, `g`, guards, handlers, `project` — once, at the initial state,
with real values, checking shape/type conformance and discarding results (all
are pure; the cost is one evaluation each). "Fails loudly at build time where
possible" ([§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)) decides this: a malformed `f` return must not wait for the
first integrator step. Probes see only the initial state's branch — the
marginal coverage is earliness, not completeness; the always-on check ([§12.5](#125-the-always-on-conformance-check))
remains the completeness backstop.

**Probe argument sourcing.** `x`/`m` come from `init_*` declarations
(declared by value); `y_x` from the stage-1 probes' *returns* (an
auto-published name is a framework write, never a probe product, so it is
absent from the hand-down — [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)); wired inputs from
upstream products, real values available because the stage-2 chain is probed in
topological order. **`w` threads through the probe in stage order**, by the
[§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws) one-hop law: the `h_x` probe's second return slot enters `h_xu`'s probe
bundle if the component defines one, otherwise the downstream bundles, and the
`w` that survives the last output stage enters the `f`, guard and handler
probes — so every consumer is probed against the same value it will receive at
run time. Three checks ride the same pass: the return's shape (a stage
returning something other than a `NamedTuple` or a `NamedTuple` pair fails
here, `nothing` in either slot included), the consumer's `w` reads against the
observed field set ([§11.3](#113-visibility-the-contract-is-the-interface)'s did-you-mean, through the [§13.2](#132-diagnostics-structured-values-one-carrier-exception) framing
diagnostic), and the dead-stage rule — a stage returning bare `(;)`, neither
ports nor `w`, is `DeadStage`, fail-fast. The bundle law's two remaining fields ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)): `t` is
probe-scoped `0.0` — deployment binds no clock and `t₀` post-dates even
deployment ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)), so like `Δt` below it is a fabricated, probe-scoped
value; `ws` comes from invoking the component's `workspace` allocator at the
probing scalar, which reads only the instance and the scalar (row 77),
deriving nothing from layouts, so it runs before the Stratum B probes that
need it. Exactly one kind of terminal has no producer: **root
slots**. The build synthesizes their values via `probe_value(::Type)`:
framework methods for `Real` (`zero(T)`), `Bool` (`false`), enums (first
instance), ultimate fallback the zero-argument constructor `T()` — which is
where well-behaved constrained types already put their valid default (`RQuat()`
= identity; the `@kwdef` convention supplies it broadly). `probe_value` is
**overridable**: a type whose valid default is not reachable that way declares
its own method, which is also the seam a walked type uses to state a
constrained default. No method → build error, in the didactic register: it
names the face and the type, and asks for one of the two fixes ("no
`probe_value` for `Ranged{Float64, -1, 1}` at face `pilot.elevator_axis` —
define `probe_value(::Type{Ranged{Float64, -1, 1}})` or a zero-argument
constructor"). Synthesis never meets an abstract type:
root slots are concrete by [§11.2](#112-the-declaration-inventory)'s tight-bound rule (the slot type is the consuming
entry evaluated at `Float64`), and abstract entries only
occur on component-fed inputs, which the probe sources from upstream products.
Physically silly values are acceptable by
construction: the probe checks types, and return types that depend on input
*values* are type instabilities (banned by the branch-shape rule); the [§4.3](#43-table-mechanics-and-port-granularity)
write-side granularity corollary keeps root slots predominantly scalar, so the
surface is small. Rejected: inputs declared by value à la `init_x` (reads as
an unwired-input default, [§6.1](#61-connections-and-hierarchy)'s banned concept; every leaf pays for a
root-only need; and it would need its own agreement rule for fan-out through an
exported face, where the slot type already follows [§11.2](#112-the-declaration-inventory)'s rule — the unique
concrete declaration among the consumers); NaN poison
values (`Bool`/enums unpoisonable; probe values are meant to be read — poison
guards illegitimate reads, [§7.3](#73-discrete-state-modes-and-workspace)); init-service values (the build is
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
enforcement is the pre-write `UninitializedSlots` check carried by every
complete-world application — `init!`, trim setup, trim commit
([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)).

**The author's side of that bargain.** Silly values are acceptable *because* the
author is obliged to accept them: **stage code must be total over type-valid
inputs** — every probed user function (stages, `f`, `g`, guards, handlers,
`project`) evaluates without throwing on any input satisfying its declared
types. The domain is type-validity, not the probe's particular synthesized
values: the branch-shape rule already bans value-dependent return types, so
types are the only domain the framework can speak of, and the probe is the
enforcement moment, not the reason — asking whether one can detect being probed
reads the contract backwards. Two consequence sites, the same throw at both: at
build it is a `UserCodeFraming`-wrapped build failure ([§13.1](#131-reporting-policy-collect-the-checks-fail-the-evaluations-fast)) whose diagnostic
points at code that is "correct" on every trajectory it has ever seen; at
runtime it is a `StepError` and the run ends `errored` ([§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)) — exceptions from
model code are always abnormal ([§13.5](#135-termination-is-a-state-not-an-exception)). Three habits of shipped code have
sanctioned spellings. A *plausibility* check meaning "stop the run" — a strut
throwing on a touchdown overload — is a published `Bool` output face plus
`stop_on` ([§13.5](#135-termination-is-a-state-not-an-exception)), machinery already there. A *self-consistency* assert — an
author checking that their own contact algebra cancels a velocity component to a
hard tolerance — is a regression test about that algebra, and its home is the
test suite; it is also the most probe-fragile of the three, since a
near-degenerate synthesized geometry can keep the cancellation algebraically
exact while missing an absolute tolerance in floating point. And a *defensive
exhaustiveness* branch — an `else error("unrecognized surface type")` over a
closed enum, or a coefficient constructor asserting an ordering of its
arguments when that constructor runs per step inside a stage — is not banned
validation but **mislocated** validation: totality over a closed enum means
handling every instance (an `else error` is an admission that the function is
partial), and parameter validation belongs where user-controlled data enters —
the constructors of parameter and instance values, which run before the build,
where asserts are perfectly legitimate — never inside a stage, on probe-fed
data.

### 12.4 Activations: executable sets, laziness, caching

An **activation at `T`** re-runs Stratum C with a different scalar:
producer-fed cells re-typed by *evaluating* the producing component's output
declaration at `T` ([§11.2](#112-the-declaration-inventory) — a continuous producer's two-argument
declaration called at the scalar, a discrete producer's plain declaration
pinning), root-slot cells by *evaluating* the consuming `input_types` entry at
`T` ([§11.2](#112-the-declaration-inventory)'s permissive entries — a `T` entry follows the activation, a
`Float64` entry stays frozen) and the
state type by the walk over `init_x`'s, table and
state buffers re-laid-out, workspace allocators re-invoked at `T` (re-invoked,
not introduced — the first invocation precedes the Stratum B probes,
[§12.1](#121-three-strata)/[§12.3](#123-probing-and-input-synthesis); a continuous component's scratch carries the activation's
scalar, [§7.3](#73-discrete-state-modes-and-workspace)), probe chain re-run. Structure and schedule are `T`-independent by construction.

**Each activation probes exactly the function set it can execute.** A `Dual`
activation (linearization, gradient trim) evaluates the model at a frozen
instant: discrete stages are gated off holding `Float64` values (the [§11.2](#112-the-declaration-inventory)
frozen-constant semantics), guards and handlers never run (event localization
is `Float64` sweeps by design, [§8.4](#84-localization-mechanics)). Only the continuous output stages
(`h_x`/`h_xu`) and `f`
ever see a `Dual` — so only they are probed. Probing the discrete stages, `g`, or guards at `Dual`
would check code against a number type it cannot receive. One rule, no
special cases; the [§5.6](#56-diagnostics-feedthrough-tracing) tracer activation follows it identically. "Tracer
activation" names row 12's *global* set-tracer — a whole-model run at the
tracer scalar, an activation like any other. [§5.6](#56-diagnostics-feedthrough-tracing)'s cycle classifier is row
12's other variant, the schedule-free per-member local trace, which runs in
Stratum B's failure path and is not an activation at all.

**Lazy, with an opt-in exhaustive mode.** Non-nominal activations run at first
request, not at build: the dominant cost is compiling the continuous chain a
second time at `Dual`, pure waste for interactive fly-around use. The price,
stated openly: `build` succeeding does **not** certify the model linearizable —
a pinned `Float64` ([§7.2](#72-numeric-genericity-eltype)), whether hidden in a constructor or written into an
output declaration at a leaf that really participates ([§11.2](#112-the-declaration-inventory)'s per-leaf
forgotten-`T`), lurks until the first `Dual` activation detonates it
at the probe, naming the offending constructor or leaf. The repository's test suite
pins the invariant instead, and row 166 makes it policy rather than advice:
**every component gets a `Dual` activation built in CI** —
`build(world; activations = (Float64, ProbeDual))`
(or a `check` entry) runs the exhaustive set, catching both genericity
violations and forgotten-`T` leaves at PR time, at the cost of a Stratum-C
re-run per component. The same keyword is the recommended idiom for [§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)'s
parallel-sweep register: pre-materialize the activations the sweep will need and
the shared `Build` is a fully immutable artifact, with no synchronization on any
path. `ProbeDual` is the framework's exported canonical probe
scalar — `const ProbeDual = ForwardDiff.Dual{ProbeTag, Float64, 1}` — because
an activation is keyed by a *concrete* scalar type and the bare `Dual`
`UnionAll` cannot key one, be walked to, or answer `zero(T)`. Its width is
arbitrary: what CI pins is genericity, not any particular Jacobian, so one
canonical width suffices even though [§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query) chunks at whatever widths it needs.

**Caching is implementation detail, not semantics.** An activation is a pure
function of the build and the concrete scalar type — so the cache is the
`Build`'s, and it holds (layouts, compiled plans, validated-flag) keyed by that
type: immutable once constructed, hence freely shareable. **Buffers are never
cached**, because every buffer set has exactly one owner. The `Simulation` owns
its nominal activation's buffers — materialized from the cached layouts at
construction, what the loop's zero-allocation stepping runs on — and every
service invocation owns the scratch set it instantiates from those same layouts;
[§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report) states this for `trim!`, and it is the general rule, not a trim-local one.
Compiled code is cached by Julia itself, process-wide; what the framework cache
saves is the expensive part — probe re-runs, layout construction, and Julia's
compilation of the `Dual` chain — which is what actually amortizes in
activation-reusing loops (the envelope-grid gain-schedule case: hundreds of
trim-then-linearize points paying those costs once). What does not amortize is
the per-point allocation of a working store set, O(model size) and trivial
against the solve it feeds: [§7.5](#75-allocation-policy-a-scoped-invariant)'s zero-allocation invariant is scoped to the
stepping loop, and the services were always allocation-tolerant. Nothing
numerical is ever cached.
Note `Dual{Tag,V,N}` carries the partial count: a different seeding width is a
different scalar type, hence a separate entry and a separate Julia compile.
**Lazy materialization is torn-state-free**, normatively: concurrent first
requests for the same activation must never expose partially populated cache
state. The mechanism is unspecified — a guard around insertion suffices, paid at
service time and never on the hot path — and since an activation is a pure
function of build and scalar, the worst benign race is duplicated work; torn
state is excluded by contract, not by luck.

### 12.5 The always-on conformance check

The probe validates each function *once*, on the initial state's branch; the
schema-authority bargain's second clause ("at first execution otherwise",
[§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)) is discharged by leaving the probe's comparison permanently in place.
At the point where the executor stores a stage return into the table, it holds
the complete expected return type at this activation — the declared types of the
names *this stage* produces, as fixed by Stratum B's stage classification: one
concrete `NamedTuple` type per (component, stage), drawn from
`output_types(c, T)` on a continuous producer and from `output_types(c)` on a
discrete one ([§11.2](#112-the-declaration-inventory)).
Auto-published names belong to no stage's expected type;
the framework writes those cells itself. The executor canonicalizes the
observed return to that type's field order by a type-level reorder
(`NamedTuple{names(Expected)}(y2)`) and performs a single
type test against that type (conceptually `y2 isa Expected`). Type-stable conformant
code: the compiler proves the return type, decides the test at compile time,
and deletes it — zero instructions. Branch-divergent code: the test survives
as a runtime check, nanoseconds on conformant branches, a loud located error
on the divergent one at its first execution. (Type-unstable-but-conformant
code pays nanoseconds on top of the dynamic dispatch it already bought.)

**The names are the pairing; field order carries no semantics.** `Expected`'s
order is an internal fact — derived from `output_types`,
stage-filtered, auto-published names removed, an order no single declaration
shows the author — and the author never reproduces it: a return spelling the
right names at the right types conforms in any order, `(; P = M*ω, M_shaft = M)`
and `(; M_shaft = M, P = M*ω)` being the same return. This is the general rule
at every author↔framework `NamedTuple` seam ([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals) states it for the trim
problem's decisions and residuals), and it is what downstream consumption
already assumes: the scatter writes each returned field into its own *named*
cell ([§4.3](#43-table-mechanics-and-port-granularity)), so order-sensitivity in the check would be incidental
strictness rather than protection. The canonicalizing reorder costs nothing:
a compile-time permutation of an already-typed value — register shuffling SROA
deletes — folding exactly where the test folds, so row 53's economics stand
unamended (one whole-type test, never per-field checks) and [§7.5](#75-allocation-policy-a-scoped-invariant)'s canary
verifies the fold empirically rather than by assertion. What remains an error is
what always was: a key-set mismatch or a per-field type mismatch, reported by
the unchanged payload below. A permutation is not an error at all — which is
equally why that diff never has to express one.

**Exact match at nominal; embed-accept at declared-`T` leaves.** At
the nominal activation — the only one that ever runs in real time — the check
is unchanged: exact type match, no convert-on-write, one baked `isa` that
folds away (`Int` sloppiness must fail at `Float64`, not lurk until a `Dual`
activation makes "it runs" activation-dependent). The error can afford to be
didactic: "field `M_shaft`: expected `Float64`, got `Int64` — return
`zero(x.ω)`, not `0`". Under a non-nominal activation the two leaf kinds
[§11.2](#112-the-declaration-inventory)'s declaration distinguishes are checked differently. A **declared-`T`
leaf** — the author wrote `T` there — accepts exactly two types: the activation
scalar (the fast path — the baked
`isa` unchanged) or `Float64`, which the executor **embeds** as a
zero-partial constant (`convert` through the leaf; struct-valued ports use
the standard cross-eltype constructor, a missing one failing loudly with
both types named). Nothing else is accepted. A **declared-pinned leaf** — the
author wrote a concrete type, `Float64` at the head of the list — takes the
nominal-style exact check at *every* activation, its declaration having said the
leaf never carries partials; an observed `Dual` there is the per-leaf
forgotten-`T` error, with the didactic hint attached ("if `F` participates in
differentiation, declare it `T`"), that being the one honest cause. The
embedding is exact, not
lenient: promotion is airtight and there is no lossy `Dual → Float64` cast,
so a `Float64` observed at a declared-`T` leaf means no `Dual` entered its
computation — its true derivative along every seeded direction is zero,
which is precisely what the embedded constant says. This scopes the blanket
convert-on-write rejection to the nominal check (row 53):
the bug that rejection guards against — silently zeroed partials — cannot
arise from honest code, because accidental `Float64`s from `Dual` operands
are impossible (`MethodError` at the operation site). The residual is
**deliberate stripping** (`ForwardDiff.value`): a stated intent to discard
partials, producing a silent zero in the Jacobian — the stop-gradient idiom,
occasionally legitimate (deliberately frozen couplings, opaque non-Julia
wrappers), and equally invisible to a strict
exact-match rule when applied mid-expression, so the leniency costs nothing.
What it is no longer is invisible to the schema *by necessity*: **the
declared-pinned leaf is the schema-visible freeze** — an author who means to
strip declares the leaf `Float64` and strips inside the stage, and the check
above holds the freeze to its word at every activation. Stripping mid-expression
at a leaf still declared `T` remains legal and remains unseen, as the sharp tool
it is.

**`w` is checked at the nominal activation, and nowhere else.** A stage that
returns the `(y, w)` pair ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)) gets a second baked `isa` beside the first,
against the type the **nominal probe observed** — no declaration exists to draw
an expectation from, and none is wanted, `w` being a value in flight rather
than a cell to type. It folds exactly as its sibling does when the stage is
stable, and when it does not fold it earns its nanoseconds: the branch that
quietly returns a `w` of a different shape is otherwise invisible until it
shows up as an allocation in someone's benchmark, and here it is a located
error naming the offending field at that branch's first execution. The message
cites its authority honestly — expected *what the nominal probe observed*, not
what anybody declared. Under **non-nominal activations the `w` test is
absent**: with no declaration there is no branch-independent anchor, there is
no cell whose typing the check would protect, and a probe-observed expectation
would reject the constant-branch idiom that [§11.2](#112-the-declaration-inventory) keeps legal on the
output side. Nothing is lost in correctness — a `Float64` arriving in `w` under
a `Dual` activation is a true zero-partial constant by the embedding guarantee
above, and promotion at its first use with a `Dual` operand is exact. The
reasoning is [§11.3](#113-visibility-the-contract-is-the-interface)'s, recorded in row 165.

**Uniform across all probed functions.** `f` checks against `X`'s own shape
at the activation's `T` ([§7.1](#71-continuous-state-structured-immutable-flat-backing): a scalar leaf expects a `T`, an `SArray` leaf
the same `SArray` at `T`), the predicate being "every
field scatters into its field's block at `T`", which is what makes derivative
completeness structural rather than a matter of author discipline; guards
against their probe-derived predicate form (below); `g`
against its leaf's `x` shape; handlers against the [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws) return law, key by key;
`project` against `X`'s own shape at `T`, **complete** — the same
predicate as a handler's `x` key, since its result is written back to the
buffer wholesale at both of [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)'s schedule positions, and a projection
with a mode-dependent branch first executes its second branch at run
time.

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
`m` lives in per-field stores where a partial merge is the natural write.

**Guards have two admissible forms** ([§2.1](#21-events-two-detection-policies)), so their check is form-aware
rather than a flat `isa Bool`: a `Bool`-form guard's probed return is `Bool`, a
sign-form guard's is the nominal scalar — guards run only at the nominal
activation (row 52), so no parametrized-leaf case arises here. Any other probed
return type is a build error naming both admissible forms. So is a localized event
whose guard probes `Bool`: localization brackets a root, and the `Bool` form
offers none — return `x.ω - eng.ω_idle`, not `x.ω > eng.ω_idle`.

**Failure payload:** component path, function, field-level diff (missing /
unexpected / per-field expected-vs-observed), simulation time. Deliberately
absent: the source branch (values carry no provenance; the diff identifies
it). The always-on input trace makes every such failure **reproducible by
replay** — the error names the boundary to replay to ([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)'s
`to_boundary`). At run time the failure
travels as a species of `StepError` through the single catch site ([§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)),
which adds the loop-level nonfinite-state check as its divergence sibling.

### 12.6 Stopped-sim services as Stratum-C clients

Sketched here because it grounds the strata; the services themselves are [§14](#14-stopped-sim-services).
The C172 trim problem (`c172.jl`: `TrimState`, `TrimParameters`,
`θ_constraint`, the `ẋ`-reading cost) transfers near-verbatim:

- **Trim** is a write-condition → sweep → read loop on an activation — by
  default the `Dual` activation, decision variables seeded for exact residual
  Jacobians ([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)); the derivative-free fallback runs the same loop on the
  nominal `Float64` activation with no new activation needed, and the
  always-on checks ride along either way. Decision variables stay opaque to the
  framework (only the assignment's *output* is framework vocabulary).
  `assign!` inverts from in-place mutation + self-invoked `f_ode!` to a pure
  function returning a condition value (state by path, modes, slots by face)
  that the service writes and evaluates. Domain math — the pitch constraint,
  `Kinematics.Initializer`, per-residual scalings and the equilibrium-subset
  choice — survives aircraft-side, with one respelling: the initializer's
  `atmosphere::Model` argument becomes a field handle ([§4.4](#44-function-valued-signals-environment-access)), built at value
  level by the atmosphere's value-level constructor or held directly as a rig
  slot value ([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults), [§14.9](#149-mounting-problems-as-relocatable-values)).
- **Linearization** is a `Dual` activation plus seeded sweeps: gather/scatter
  over the canonical layout replaces the hand-written
  `get_x_ss`/`assign_x_ss!` layer ([§7.1](#71-continuous-state-structured-immutable-flat-backing)'s deletion discharged); root slots are
  the input surface; frozen discrete outputs are constants with zero partials,
  which is exactly "linearize with the discrete state held" ([§11.2](#112-the-declaration-inventory)). Gradient-based trim —
  decision variables seeded through the `T`-generic assignment math — is the
  default ([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)).
- The generic service loop (vectorization, optimizer setup, bounds packing,
  solved-condition write-back including root slots and the trace header's
  slot capture) replaces today's per-aircraft NLopt plumbing. A failed trim
  leaves the simulation's stores untouched — an improvement over today's
  warn-but-assign `f_init!`.

### 12.7 The compiled executor

The schedule exists in two representations at two lifecycle stages. In the
`Build` it is plain printable data ([§12.2](#122-the-build-artifact)) — paths, stage names, order — the
authoring and diagnostic form. At `Simulation` construction, and per
activation, that data is compiled into the execution form: **a
concretely-typed tuple of entries over statically typed cell storage,
traversed by a compile-time-unrolled walk**. This is a forced move, not a
preference: [§7.5](#75-allocation-policy-a-scoped-invariant)'s zero-allocation invariant, [§12.5](#125-the-always-on-conformance-check)'s fold-away conformance
test and [§5.1](#51-the-scheduling-problem)'s zero runtime graph logic are reachable only under full
specialization. A heterogeneous runtime vector of entries dispatches
dynamically per entry and boxes every stage return — the framework itself
would allocate, and [§7.5](#75-allocation-policy-a-scoped-invariant)'s canary role (framework contribution exactly zero,
so any model inference regression announces itself) would die with the
invariant. An entry carries what selects code — component type, stage — in
type parameters, and what is plain data — tick divisor, layout offsets — in
fields; gating compiles to an integer test inside the specialized *boundary*
body, the interior bodies holding no discrete entries to test ([§8.5](#85-multi-rate-tick-scheduling)).

**Cells are stored per element type, not per cell.** The signal table is one
contiguous buffer per element type — [§7.1](#71-continuous-state-structured-immutable-flat-backing)'s construction pointed at signals
rather than state — and a cell address is a build-time offset into it, carried
in an entry *field* with the port type as the address's own parameter; gathers
reconstruct and scatters flatten through the same leaf walk, so the closed
vocabulary earns its keep twice. This is the entry rule above paying rent: two
instances of one component type then differ only in field values, share an
entry type, and compile to **one** body — where a store enumerating every cell
in its own type, addressed by index in the type domain, compiles one body per
instance and grows the store type with the model. Measured rather than argued
(row 162, `prototypes/cellstore_bench`): at 400 identical instances the two
representations cost 1.1 s and 56 s of cold compile at equal zero allocation,
and on the 8-partial `Dual` activation the per-element-type store's compile
cost *saturates* — bounded by chunk-type count, not model size — where the
type-domain one keeps climbing.

**Phase bodies are the outer decomposition, and they are semantically
forced.** The boundary sweep's `h_x` block is order-free by definition (the
no-feedthrough stage reads no `u`); the `h_xu` block, with due discrete
stages gated in, is the only topologically ordered one; the `f` block — the
RHS body the stepper calls per stage evaluation — and the `g` block are
order-free with disjoint writes; guards and handlers are their own small
callables inside the [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) iteration. **Each sweep block compiles in two
arities off one entry list** ([§8.5](#85-multi-rate-tick-scheduling)'s interior/boundary split): the zero-arg
`sweep_hx()`/`sweep_hxu()` are the interior variants, over continuous entries
only — which is what makes `@ballocated(sweep_hxu()) == 0` a well-defined
measurement *of the interior path*, rather than of whichever tick phase the
simulation happens to be sitting in — while `sweep_hx(tick)`/`sweep_hxu(tick)`
are the boundary variants, gating their discrete entries by modulo against the
passed index, symmetric with `ticks(tick)`; `rhs` takes no index (row 147,
amending row 116). Across passes these bodies communicate only
through the stores and the table, so the seams between passes cost
nothing — no values cross them. **Within** a pass one kind of value does: a
stage's `w` ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)) is handed to its one-hop consumers as an ordinary SSA
argument, across block and chunk seams alike, never through storage. That is
what makes the fusion a design constraint rather than an optimization: a step's
sweep and its `f` block compile into one pass, and an event round's sweep, its
guards and its fired handlers into another ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)), because that is the scope
over which a private intermediate is fresh by construction. Two doors this structure opens for free,
recorded not committed: deterministic parallel evaluation of the order-free
blocks (disjoint writes, and no floating-point reductions to reorder — [§6.2](#62-aggregation-explicit-summing-junctions)
made every sum an ordered junction entry), and finer recompilation
granularity (editing a discrete component invalidates the boundary body, not
the RHS body — literal under the two-arity split, discrete entries existing
only in the boundary variants).

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
C172X-scale model ([§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load); roughly 200–400 entries, larger bodies), the
nominal activation sits at seconds, the `Dual` activation in the tens before
mitigation; re-measurement on the real vehicle skeleton is a [§16](#16-open-axes) migration
item.

**The mitigation ladder**, in order: activations are lazy ([§12.4](#124-activations-executable-sets-laziness-caching) — a session
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
(`u`, `y_x`) are per-call by topological necessity either way ([§7.1](#71-continuous-state-structured-immutable-flat-backing)).

**Construction is type-opaque; only the executor specializes.** Schedule
tuples are built from untyped buffers and splatted once. Generic tuple
utilities — range indexing, long `ntuple` closures, naive recursion — are
inference traps at schedule length (a 400-entry heterogeneous tuple can send
generic `getindex` inference into combinatorial collapse), so the compiled
tuple's type has exactly one consumer: the unrolled walk.

**The phase bodies are the [§7.5](#75-allocation-policy-a-scoped-invariant) measurement seam.** `phase_bodies(sim)`
returns the compiled bodies of the nominal activation as named callables
bound over the simulation's own buffers: the four blocks (`rhs` — the `f`
block — `sweep_hx` and `sweep_hxu`, each in both arities, and `ticks`, which
takes the tick index its entries gate on), plus the per-event guards and
handlers and the per-component `project` callables, keyed by the model's own roster.
The four-body roster is fixed and total: the accessor returns all of it
always, whatever the model happens to declare. A model with no discrete
components, no events or no continuous state at all still gets every body;
the empty ones are legal, compile to no-ops, and their `@ballocated`
assertion passes vacuously — which is the point, because consumers then
iterate the roster uniformly, with no existence checks and no per-model
branching in the measurement code. One
promise, in the diagnostic register ([§13.5](#135-termination-is-a-state-not-an-exception)): **these are the bodies the loop
runs** — not re-derivations, which is what makes the measurement honest, and
why each callable carries the real in-loop argument types by construction
(the thing a hand-built standalone test cannot reproduce, and the reason
per-component tests cannot discharge the invariant: they never call the
executor whose contribution the claim quantifies over). CI is
warm-then-assert over the roster — one call compiles, then
`@ballocated(body()) == 0` — at per-body granularity, each sweep arity asserted
in its own right (the interior call bare, the boundary call at a due index), so
a documented [§7.5](#75-allocation-policy-a-scoped-invariant)
tolerance loosens exactly one assertion; this is the successor of the
migration suite's `@ballocated f_ode!`/`f_step!`/`f_periodic!` idiom and the
seam the [§16](#16-open-axes) FlightCore comparison measures through. Publication is not a
phase body — [§7.5](#75-allocation-policy-a-scoped-invariant)'s carve-out made structural: what the accessor exposes is
exactly what the invariant claims is zero. Invoking bodies in isolation
mutates the simulation's buffers outside any frame sequence (a tick entry
advances discrete state with no clock advance), leaving them valid but
off-trajectory; a session that wants to continue meaningfully re-runs
`init!`.

Rejected (row 86): vector-of-abstract entries (dynamic dispatch boxes
returns; union splitting rescues only toy scale and cannot close an open
method set); type-erased call tables (`FunctionWrapper`-style — the same
specialization count as chunk-of-one with extra machinery, no cross-entry
optimization, and a load-bearing seam on internal ABI, the [§8.1](#81-loop-ownership-the-framework-owns-the-simulation-loop) lesson);
framework-maintained view hoisting (lifetime obligations the compiler
already discharges, where mis-scoping means silent stale-state reads).

---

# Part IV — Failure and services

## 13. Error discipline

[§11.4](#114-failure-walkthroughs-the-error-locality-grounding) fixed what must be caught and where; [§12](#12-the-build-pipeline) fixed when each fact is checked.
This section fixes how failures are *reported* — the reporting policy, the
diagnostic representation, the runtime failure story, and the seam between "the
model reached a terminal state" and "the run should end". Two of FlightCore's
paid-for lessons ground it: the compact-backtrace discipline (parameterized
model types make rendered output unreadable) and the `SimulationTermination`
machinery, which [§13.5](#135-termination-is-a-state-not-an-exception) replaces.

### 13.1 Reporting policy: collect the checks, fail the evaluations fast

The fail-fast vs. compiler-style question dissolves once the build's failure
sites are split into their two populations:

- **Declarative checks over collected structure** — unconnected inputs,
  two-producers, wire typos and type mismatches, face-name uniqueness,
  `output_types`/state-field consistency, `rates` validation. Each is a
  pass over a list; the whole-tree obligation check literally computes *the
  set of* inputs whose obligation chain never terminates. Reporting every
  violation is the natural output of the pass — truncating to the first would
  be extra work — and these failures cluster in practice (a freshly written
  assembly has five unwired inputs; a renamed port breaks three wires).
  **These passes collect:** each returns its full violation list.
- **User-code evaluation** — boundary-connection bodies (Stratum A), the stage-1 probes
  (B), the probe chain (C). When user code throws there is no meaningful
  rest-of-collection: a failed `input_connections` leaves the parent's face derivation
  undefined; a failed stage-2 probe starves every downstream probe of its wired
  inputs (probe values flow topologically, [§12.3](#123-probing-and-input-synthesis)). Continuing past these
  requires genuine compiler machinery — poisoned nodes, cascade suppression,
  dependent-check skipping — and buys little, because user-code failures are
  typically singular. **The first user-code exception aborts the phase.**

Strata are barriers: a stratum that produced any error-severity diagnostic, of
either kind, throws before the next stratum begins — probing against
unresolved wiring is meaningless. The only partial results ever carried past a
failure are violation lists from pure checking passes, so the cost that kept
this decision open — build phases having to carry partial internal results
across a failure, the machinery [§12](#12-the-build-pipeline)'s strata would otherwise need — never
materializes.

**No cascade suppression within a stratum** — a deliberate simplification. A
typo'd wire (`:throtle`) produces both the did-you-mean error and an
unconnected-input error for the intended `throttle`; both are reported. They
render adjacently (diagnostics sort by path), the pairing is self-explanatory,
and suppression heuristics are exactly the fussy machinery this split avoids.

### 13.2 Diagnostics: structured values, one carrier exception

A diagnostic is a plain value from a **small closed set of kinds**, enumerated
normatively in **[Appendix C](#appendix-c-the-diagnostic-kind-set)** — kind name, payload fields, owning section,
severity — which is the artifact the [§11.4](#114-failure-walkthroughs-the-error-locality-grounding) acceptance tests and the
error-message work are written against. Each kind carries its structured payload: endpoint paths, face names,
expected/observed types, the *list-in-hand* for did-you-mean rendering, and a
severity. Checking passes return diagnostics; the stratum barrier throws a
single `BuildError` wrapping the collection; `showerror` renders it compiler-style,
grouped by kind and sorted by path. A user-code exception is wrapped in a
framing diagnostic — component path, which function, the probe context
including synthesized inputs — with the original exception as `cause`, so the
didactic frame renders first and the raw throw second. One class is recognized
rather than merely framed: a `FieldError` — type and field carried as data —
matched against the bundle's own NamedTuple type becomes [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)'s bundle-law
did-you-mean, with the legal set and the undeclared-store / wrong-tier /
illegal-for-this-function classification. Nothing is recovered by reading
message text.

The [§11.4](#114-failure-walkthroughs-the-error-locality-grounding) walkthroughs as acceptance tests target diagnostics: tests match on
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
the one [Appendix C](#appendix-c-the-diagnostic-kind-set)'s build kinds ride: warnings there carry warning severity,
render with the collection, and never trigger the throw. That stream's warning set
is **currently empty** — the unconnected-output warning having been rejected as
its sole candidate ([§6.1](#61-connections-and-hierarchy), row 84) — and better an empty, trusted stream than a
noisy one; a warnings-as-errors CI switch is addable, not built. The *runtime*
status/log stream is a different channel and is not empty: it is
per-occurrence, carried by [§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)'s per-writer diagnostic cells — which is
where the rate limit lives, as a structural bound (a bounded ring plus
per-kind suppressed counts, drained at frame top) rather than a policy
layered over the stream, so "rate-limited wherever its source can repeat"
holds of every kind below without any kind arranging it — and surfaced
through the published framework status ([§9.2](#92-outbound-snapshot-publication)) alongside the [§8.7](#87-real-time-pacing) pacer
diagnostics and the [§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget) liveness heartbeats, which ride in the same cells
— never collected, since there is
no collection to join. Nothing in row 84's argument applies to it: that decision is
about what the *build* warns on. A *service* warning (`TrimCommitEvents`,
`TrimCommitResiduals` — [Appendix C](#appendix-c-the-diagnostic-kind-set)) is
neither stream: a synchronous per-call annotation, emitted once at a
stopped-sim service call's return beside the value it returns, its payload
duplicated as plain report fields — no carrier cell, no collection, no rate
limit to arrange. The committed runtime warnings, in one place:

- **chattering / localization-budget exhaustion** ([§8.4](#84-localization-mechanics)) — a localized event
  whose bracketing budget runs out at a boundary;
- **event re-enabled within a boundary** ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)) — deferred to the next
  boundary by the once-per-event rule;
- **forgiven-debt re-anchor** ([§8.7](#87-real-time-pacing)) — the pacer abandoning accumulated debt
  and re-anchoring its schedule;
- **the write-surface and entry violations** ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)), all at staging — the
  drain checks nothing: `OutOfClaimEntry` (an
  enumerated surface's binding drift — no position in
  the attach-compiled schema), `ClaimedFaceEntry` (a harness write to
  a face claimed in the run's frozen partition, naming the incumbent; also
  fired by the stopped-sim attach renormalizing a pending batch) and
  `EntryTypeMismatch` (a value unconvertible to its
  slot's declared type, rejected at staging for every writer);
- **a tolerated device-side datum failure** ([§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)):
  `MalformedDatum` — emitted by the author's loop via `report!(handle, …)`
  into the device's own cell ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)), the `InputMappingError` successor;
- **staging discarded during replay** ([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)): `ReplayDiscardedStaging` —
  a live batch found staged while the trace feeds the drain;
- **thread-budget tightness** ([§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget)) — once per `run!`, against the frozen
  roster;
- **device join timeout** ([§10.4](#104-shutdown-protocol)) — a device task exceeding the shutdown
  join timeout, abandoned by name rather than hanging `run!`;
- **device crash** ([§10.4](#104-shutdown-protocol), [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)) — a device task's failure caught by the
  framework wrapper, the sim continuing with the device absent;
- **workspace poison skip** ([§7.3](#73-discrete-state-modes-and-workspace)) — once per activation, naming the stores
  with no usable sentinel;
- **unbounded run** ([Appendix B](#appendix-b-api-synopsis-the-entry-points)) — no finite `t_end`, no `stop_on` faces,
  `pace = Inf` at run start.

### 13.3 Build primitives: `resolve` and the face-list accessors

The [§11.8](#118-computed-connections-and-generic-boundaries) sketch's primitives, made normative:

- `resolve(asm, path::String) → AbstractComponent` — the getfield walk along
  `/`-segments. Its one non-obvious duty is enforcing [§6.1](#61-connections-and-hierarchy)'s generic-boundary
  rule: it walks *declared field types* alongside instances, and a segment
  that traverses **past** a generically-held field (non-concrete declared
  type) is a diagnostic even though the concrete instance in hand would
  resolve it — resolving *to* a generic child is port-level access and legal.
  An unknown segment errors with the sibling field list in hand.
  **The duty is register-scoped** — row 83's load-bearing/diagnostic line
  carried into resolution, client policy riding on one primitive exactly as
  in [§14.4](#144-two-application-registers-over-one-plan):
    - *Structural* (wiring resolution, Stratum A): the strict rule above,
      verbatim — the register [§6.1](#61-connections-and-hierarchy)'s law lives in.
    - *Load-bearing* (condition entries, trim `reads`, taps — [§14.3](#143-resolution-flatten-validate-compile-once), [§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals),
      [§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)): strict, evaluated **at the authoring or mount level**. The
      locality law is an authoring-level law ([§14.2](#142-fragment-composition-locality-without-schema): absolute paths are a
      compiled derivative), and a mount prefix is checked by the mount
      itself ([§14.9](#149-mounting-problems-as-relocatable-values) validates the problem's declared type against the mount
      point's contract) — this register checks the authored path below it.
    - *Diagnostic* (device read bindings, GUI panels, snapshot and log
      inspection — [§9.2](#92-outbound-snapshot-publication), [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)): the instance walk. A generic seam is not an
      error for a client that never claimed substitutability — "what is in
      *this* build" is the inspection register's defining question — and
      drift stays loud: an unknown path is an attach-time
      `ReadBindingUnresolved` with did-you-mean.
  Which register a client resolves under is internal framework fact, never
  user-facing API — the same status as [§14.4](#144-two-application-registers-over-one-plan)'s two `apply!` registers.
- `input_faces(c)` / `output_faces(c) → Vector{String}` — the stringified keys of
  a leaf's `input_types` (the key set is `T`-independent); the entries of
  `input_connections(c)` / `output_connections(c)` for an assembly. Declaration order is preserved: deterministic
  printouts, stable diagnostics.
- `resolve_terminal(asm, path) → (component, name)` — splits a terminal
  path's final segment (unambiguous: face names may contain dots, never
  slashes, [§11.6](#116-paths-wiring-and-faces)) and resolves the prefix through `resolve`. First-class
  because five clients share it across the three registers — wiring
  resolution (structural); condition addressing ([§14.3](#143-resolution-flatten-validate-compile-once)) and tap
  resolution ([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)) (load-bearing); device-binding validation ([§9.2](#92-outbound-snapshot-publication))
  and snapshot inspection (diagnostic) — one splitter, one did-you-mean
  site.

### 13.4 Runtime failures: one catch site, an execution cursor

**Where caught.** The loop wraps each execution of the boundary macro-sequence
(integrate → project → event iteration → ticks → publication) in a single
`try` — never per stage or per component, which would salt the hot path with
exception frames for no benefit. Framing information does not need to be
*caught* into existence: the executor maintains an **execution cursor**, a
plain mutable field in the loop state recording where in the compiled schedule
it is — component path (schedule index), which function
(`h_x`/`h_xu`/`f`/`g`/guard/handler/`project`), and the boundary phase
(integration stage *k*, event round *r*, localization probe at trial
time, tick). One cheap store per dispatch on a single-tasked executor — no
allocation, no exception frames — and it covers every user-code surface
uniformly, including the forgettable ones: RHS evaluations at interior RK
stage points, guard evaluations inside ITP/Brent probes, environment closures.

**How handled.** The catch site wraps the original exception in `StepError` —
the runtime counterpart of `BuildError` — carrying the cursor's frame, the
boundary time, the **frame-entry boundary index** (the replay pointer: the
frame-top boundary — grid or boundary zero — at which the failing frame
began, always a legal replay halt, [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)), and the original
exception as `cause`, rendered with compact frames per [§13.2](#132-diagnostics-structured-values-one-carrier-exception)'s doctrine. The
[§12.5](#125-the-always-on-conformance-check) conformance failure needs no separate path: it is thrown as its typed
diagnostic at the table-write point and arrives at the same catch site, a
species of `StepError` with the field-diff payload. Reproducibility holds by
construction: staged inputs are drained and recorded to the trace at the frame
top, *before* the boundary executes, so the failing boundary's inputs are
already in the trace when it fails — the error names the frame-entry
boundary `k` to replay to, and `replay!(sim2, trc; to_boundary = k)` halts
exactly at that frame top; `step!` then re-executes the failing frame under
instrumentation, localized boundaries included ([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)).

**The one exception never wrapped.** An `InterruptException` is not model code
failing but the operator's stop command ([§10.4](#104-shutdown-protocol)), so the catch site discriminates
it and routes it to the stop path: the run takes the ordinary graceful tail and
ends `stopped`, never `errored` under a `StepError`. With [§10.4](#104-shutdown-protocol)'s boundary
masking in force this branch is unreachable in practice — the interrupt is
deferred to a frame-top or wait unmask point and never raises inside the guarded
sequence — and it is kept defensively because the cost of being wrong about that
is a terminally errored session in place of a clean stop.

**Disposition.** The `Simulation` ends in a terminal status — `stopped` vs.
`errored` — with the exception retrievable. A synchronous unattended run rethrows
after the shutdown tail completes, so CI fails honestly; an interactive
session logs the rendered error and surfaces the status through the control
plane and GUI.

**The nonfinite check.** Divergence is not termination: dynamics that blow up
(ground penetration, an unstable gain) produce NaNs that defeat guards — NaN
comparisons are false — so no declared condition will catch them. A loop-level
`isfinite` sweep over `x` at boundaries fails fast as a `StepError` species
naming the offending component's state block and the boundary. It catches
diverging models generally, not just post-terminal ones.

*Placement is the whole value.* The sweep is the boundary's **first act** —
immediately after integrate returns, before `project` and before the boundary
sweep. Run there, `NonfiniteState` names the component whose own block
diverged. Run later, the NaN has already propagated: it reaches an innocent
downstream component through the ordinary signal path and surfaces as that
component's lookup-table `DomainError` or an `InexactError` in its
conversion — [§11.4](#114-failure-walkthroughs-the-error-locality-grounding)'s error-locality inversion, designed out of the build
tier and quietly reintroduced at runtime. One `isfinite` pass over a flat
buffer is cheap enough that placement, not cost, decides.

*Scope: `ẋ` does not participate.* A nonfinite derivative contaminates its own
state block's step result within that very step, so the `x` check at the next
boundary is the same detection with identical component attribution — the
`ẋ` sweep would report nothing new, one boundary earlier. And `ẋ` buffers are
integrator scratch: written per stage, meaningful only inside a step, and not
boundary-consistent in the sense the check is stated over.

**Domain separation.** Device-side user code — loop bodies and mappings run
on the device task ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain), [§9.6](#96-devices-one-authoring-contract-no-taxonomy)) — fails in the device's own domain, and in
two classes: a genuine bug takes the per-device crash path (liveness
heartbeat, `DeviceCrash`) while the sim keeps running; an unmappable datum
is not a failure at all — the loop body tolerates and reports it
(`MalformedDatum`, [§9.6](#96-devices-one-authoring-contract-no-taxonomy)). The two failure domains never mix — exactly what
the no-shared-mutable-model decision bought.

### 13.5 Termination is a state, not an exception

FlightCore's `SimulationTermination` idiom — model code throws, the loop
catches and logs it as informational — has **no counterpart here**. It sits badly
in this design: a mid-sweep throw aborts a boundary halfway, and [§10.4](#104-shutdown-protocol) is built on
completing one. The discipline: **exceptions from model code are always
abnormal**; graceful termination is model *state*, reaching the loop through
declared machinery:

- **Detection** is ordinary guard/handler/mode machinery. Declare the
  predicate as a localized event if the stop should be localized — touchdown
  overload is precisely a zero-crossing: the boundary is localized to the
  crossing, the handler sets `m.crashed`, and the snapshot at the crossing
  instant carries the touchdown state.
- **Publication** is an ordinary `Bool` output face, exported to the root.
  Within concretely-declared structure, deep wires gather the condition at its
  owning boundary in one visible block ([§6.1](#61-connections-and-hierarchy)): `Ldg` ORs its three legs through
  a junction ([§6.2](#62-aggregation-explicit-summing-junctions)'s ownership idiom, [§13.7](#137-tooling-consequences-provenance-and-the-component-library)'s library) and exports one `damaged`
  face; intermediate assemblies are untouched. Each *generic* seam costs one
  output connection entry — and that hop is the substitutability contract doing its job,
  not plumbing ([§11.8](#118-computed-connections-and-generic-boundaries)'s imposed contract).
- **Policy** binds at deployment: `Simulation(world; ..., stop_on = (...))`
  names root-exported `Bool` output faces, OR-combined, validated against the
  `Build`, recorded in the run metadata (the trace header's deployment block,
  [§9.5](#95-inbound-the-input-trace)). After *every* published boundary —
  grid, `t*` ([§8.4](#84-localization-mechanics)) and boundary zero ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)) alike — the
  loop reads the named faces in the snapshot it just published; the first
  `true` initiates [§10.4](#104-shutdown-protocol) shutdown with *this* snapshot as the final one — the
  terminal snapshot is the terminal state, no roll-back, nothing [§10.4](#104-shutdown-protocol) doesn't
  already do (and its status carries the run's final cumulative diagnostic
  counters, [§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)). `run!` therefore checks the boundary-zero snapshot before the
  first step: an authored condition already terminal ends the run at `t₀`
  with that snapshot final, integrating nothing. Default: no stop faces, run
  to `t_end` — `stop_on` is `t_end`'s model-declared sibling at the same
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
more deployment-flavored, not less — and [§10.6](#106-run-lifecycle-and-partial-advance)'s `stopped → init! → run!`
cycle and `step!` register are precisely where one `Simulation` wants different
stopping policies on different runs. The honest cost is two homes for one fact,
which the precedence rule above is what settles.

**The termination record names the source.** Where run metadata carries the
effective *policy*, the run's termination record carries its *outcome*: which of
the four sources ended the run — `t_end` reached, a named `stop_on` face reading
`true`, a control-plane stop (GUI button, device handle, calling code), or an
operator interrupt ([§10.4](#104-shutdown-protocol)) — so a `stopped` simulation answers "why did it
stop?" without its consumer reconstructing the answer from the clock. The
interrupt is a tag on an ordinary stop, not a diagnostic kind of its own
([Appendix C](#appendix-c-the-diagnostic-kind-set) gains nothing here): nothing failed.

Taught contract: **stop faces are sampled at completed boundaries; declare an
event if you need the stop localized.** Both stop-flag shapes work without
framework latching — a handler-set `m` flag is sticky by nature, a transient
stage-2 Bool is caught because the loop reacts to the first `true`. Compound
stop logic composes in-model — a monitor component reading the relevant
signals and outputting one Bool — the same move [§10.5](#105-scripts-and-the-mid-run-mutation-doctrine) made for scripts.

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
  runs past the condition to inspect, unattended studies log it and continue,
  services don't step at all. (A root-declared *default* overridable at the
  ctor remains the one variant worth reopening if migration shows the ctor
  argument chronically forgotten.)
- **Blessed terminal types/names scanned from the tree, and `terminal` event
  flags**: action at a distance — a deep declaration halts the world, the
  root contract says nothing, substituting an aircraft silently changes when
  runs end, and per-deployment disabling needs masking machinery. The
  localization that terminal events appear to buy structurally is already
  available under `stop_on` via the event idiom.
- **Control-plane capability for components**: [§10.1](#101-control-plane) separates the control
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
  internals no contract mentions: precisely the knowledge violation [§6.1](#61-connections-and-hierarchy)
  forbids for wires (which bans deep paths into generic children *even where
  the concrete instantiation would resolve them*), brittle under substitution
  for reasons the aircraft's contract never promised anything about, and it
  leaves the root face list lying about what can halt a run. Output devices
  are not a counterexample but the other half of the same doctrine:
  their reads are diagnostic — snapshot-path bindings, [§9.2](#92-outbound-snapshot-publication)/[§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load) — while
  `stop_on` is the one read that changes what the run does, which is exactly
  why it alone must speak the contract.

The wall-clock channel — GUI stop button, device handle, code — is orthogonal
and untouched: that is the control plane's operator path. The sim-time,
model-detected channel specified here meets it at [§10.4](#104-shutdown-protocol) and nowhere else.

### 13.6 Abnormal shutdown: one tail, two entries

Why a `StepError` cannot break [§10.4](#104-shutdown-protocol): **the boundary is all-or-nothing outside
the sim task**. Sweeps write into table buffers, integration intermediates
live in framework-owned integrator buffers (never any component's `workspace`
in [§7.3](#73-discrete-state-modes-and-workspace)'s sense), and the only externally visible act is snapshot publication
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
named joins with timeout. This fills the seat [§10.4](#104-shutdown-protocol)'s "loop failure runs the
same protocol from the catch path" reserved.

Tail hygiene: the hooks are user code too, so each is individually
caught-and-logged — shutdown runs to completion even if a device's hook
misbehaves — and the join timeout already bounds a hook that hangs rather
than throws.

What is lost is quarantined: the state stores may hold mid-boundary values (a
half-written `m`, integration intermediates). They are retained on the
errored `Simulation` for post-mortem inspection, but an errored sim is
terminally stopped, not resumable — the reproduction tool is trace replay,
not resurrection, and the stopped-sim services enforce it by refusing an
errored simulation outright (`ServiceLifecycle`, [§14](#14-stopped-sim-services)). The published record (snapshot chain, log, trace) ends at
the last consistent boundary; nothing downstream of the sim ever sees half a
boundary.

### 13.7 Tooling consequences: provenance and the component library

Computed boundary connections gain protagonism under this section — termination chains are
their second structural customer after generic-boundary contracts — and two
commitments, a library and an idiom follow:

- **The `input_passthrough` helper family grows deliberately.** Predicate-based selection
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
- **A standard component library** makes good on [§6.2](#62-aggregation-explicit-summing-junctions)'s junction promise: when
  reduce-ports were rejected, the argument leaned on explicit junctions being
  *cheap*, and a junction hand-written per arity per type is not. Starting
  inventory strictly from demonstrated need — wrench/scalar summing
  junctions, the Bool gates the termination chains use, `UnitDelay`, the
  spelling [§5.5](#55-algebraic-loop-policy-reject-at-build-time)'s second loop-breaking remedy needs, and
  `Constant{V}`, the source block — growing by migration demand only
  (Simulink's library is a language; this is a toolbox). Doctrine: **library blocks are ordinary components** — no
  framework privileges, no special vocabulary — which keeps schema authority
  total and makes the library a permanent ergonomics torture test: if a
  three-input OR gate is painful to write under the declaration rules, the
  rules are wrong. Arity comes from a type parameter (`Or{N}` builds
  `(in1 = Bool, …, inN = Bool)` programmatically — [§11.2](#112-the-declaration-inventory)-blessed derivation,
  and an early validation that the contract functions support parametric
  components). Tier-transparency falls out of settled semantics: a stateless
  continuous `h_xu` recomputes every sweep, so fed ZOH-held discrete signals
  its output changes only at ticks — no tier-neutral class needed.
  `UnitDelay{V}` is a **discrete** leaf at `K = 1`: `init_x = (v = zero(V),)`,
  `h_x` publishing the stored value and `g` storing the incoming one — the
  tier's native `z⁻¹` ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event), the shape [§15.2](#152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade)'s `sat_out_0` hand-writes), no
  framework support. Its tier semantics are the point and must be stated
  wherever the remedy is recommended: inserting one into a *continuous* loop
  moves that signal onto the discrete tier and inserts a `Δt_base`-scale ZOH
  into the model's mathematics — a modeling decision, the delayed signal being
  genuinely sampled, not a transparent wire, which is why [§5.5](#55-algebraic-loop-policy-reject-at-build-time)'s diagnostic
  says so rather than offering the remedy as free. `Constant{V}` is the
  **source block**: no inputs, no state,
  `output_types(::Constant{V}, ::Type{T}) where {V, T <: Real} = (out = V,)`,
  and a stage-1 body returning
  the value the instance holds — a stateless continuous leaf, so the
  tier-transparency argument above already covers discrete consumers and no
  discrete variant is needed. The declaration takes the activation scalar and
  ignores it, which is the point: the block's output *is* its stored value, so
  the leaf is **deliberately pinned** at that value's own type — a
  `Constant{Float64}` declares a `Float64` port and means it, and
  [§11.2](#112-the-declaration-inventory)'s embedding turns it into the zero-partial constant it already was
  under any `Dual` activation. The honest pin is now spelled rather than
  inferred. Two demonstrated needs admit it under this
  bullet's charter: [§6.2](#62-aggregation-explicit-summing-junctions)'s zero-contributor configurations, where a
  required aggregate input has no physical contributor and the zero total
  must be spelled as a wire, and the rig stub below. Its value is instance
  data, like junction arity — not an overridable default: a configuration
  wanting an externally settable source uses a root slot ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)), which
  keeps the block from drifting into a back-door input default. A
  migration-phase deliverable.
- **The component test rig** is the library's companion idiom: a one-child
  assembly whose `input_connections` surface the child's entire input face set —
  `input_passthrough(rig, "child")` verbatim ([§11.8](#118-computed-connections-and-generic-boundaries)) — so any component can be built
  and simulated in isolation: every input becomes a root slot fed by
  ordinary conditions and devices, and every output is observable in the
  snapshot table. One qualification, from [§11.2](#112-the-declaration-inventory)'s root-slot rule: an
  *abstract* input entry (`terrain = AbstractTerrainField`, [§4.4](#44-function-valued-signals-environment-access)) cannot
  surface as a root slot — abstract-at-root is a build error — so the rig
  satisfies it *inside* the rig: a concrete stub child (a
  `SampleTerrainField` provider) wired to the face, and the concrete
  remainder exposed via `input_passthrough(rig, "strut"; except = ("terrain",))`. That
  stub child is typically just a `Constant` holding the test handle — the
  source block's first shipped instance — with bespoke stubs remaining
  ordinary components wherever the double must compute something.
  Zero new machinery — wiring and `except` already exist — and it is the
  substitutability contract doing its job: an abstract entry declares that
  a substitute must be chosen, and the rig choosing its test double
  explicitly, as ordinary inspectable code, is precisely the isolation the
  rig exists to provide. `design_world`'s little sibling ([§14.9](#149-mounting-problems-as-relocatable-values)): where
  `design_world(ac)` mounts an aircraft in a minimal world for trim and
  linearization, the rig mounts one component behind a root contract for
  unit tests and open-loop probing. Deliberately ordinary machinery end to
  end — no framework support, which makes it, like the library blocks, a
  standing ergonomics test of the declaration rules.

---

## 14. Stopped-sim services

[§12.6](#126-stopped-sim-services-as-stratum-c-clients) previewed the services — initialization, trim, linearization, capture —
as Stratum-C clients. Everything they share reduces to one artifact: the
**condition value**, the datum that says "set this build to this state."
[§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)–[§14.4](#144-two-application-registers-over-one-plan) settle its representation, composition and application;
[§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)–[§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator) the boundary-zero sequence and slot totality; [§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)–[§14.9](#149-mounting-problems-as-relocatable-values) the
trim service in full; [§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query) linearization and `capture`.

**Lifecycle preconditions.** Every service requires a non-running
simulation — pause included, by the roster freeze's own doctrine ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster),
[§10.1](#101-control-plane)): while a run exists the loop owns the stores between drains, and a
service reading or writing them would race it. Within the stopped-sim
states, legality follows each service's inputs: `capture` requires
`initialized` or `stopped` (it reads committed, boundary-consistent
stores); `init!` and `trim!` are additionally legal from `built` — their
inputs are authored conditions (`trim!`'s scratch world is built from
`override(baseline, condition(guess))`, [§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report), never from the sim's
stores); `linearize` inherits `capture`'s precondition when its operating
point defaults to `capture(sim)`, and `init!`'s legality with an explicit
`about` ([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)). **`errored` is terminal for all four** (row 59): the
errored sim's stores may hold mid-boundary values ([§13.6](#136-abnormal-shutdown-one-tail-two-entries)), and a captured
condition is indistinguishable from a healthy one once produced —
`capture(errored) → init!` would be resurrection with extra steps, and
`linearize` about a half-written point the warn-but-assign failure mode
[§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report) exists to kill. Post-mortem inspection of an errored sim's stores,
log and trace stays available as a diagnostic read; it may not become a
condition value. A violation is `ServiceLifecycle` ([Appendix C](#appendix-c-the-diagnostic-kind-set) — the
operation, the current status, the legal statuses), the same kind
`attach!`/`detach!` raise while `running` ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)): one register for "this
operation is illegal in the current lifecycle state," distinct from
`MissingInit`, which names a missing prior step.

### 14.1 Conditions are path-addressed overlays on the declared defaults

A condition may specify: state fields (`x`, on either tier) and modes (`m`,
continuous components only, [§3.2](#32-periodic-discrete-component)) — addressed by
[§11.6](#116-paths-wiring-and-faces) slash path plus field name — and root
input slots, addressed by face. Never outputs (derived data) and never
workspace (scratch). Entries are validated in the [§13.1](#131-reporting-policy-collect-the-checks-fail-the-evaluations-fast) collecting register: full
list, violations collected, one `BuildError`.

**Overlay base = the declared defaults, always.** Every store has a declared
initial value (declaration-by-initial-value, [§11.2](#112-the-declaration-inventory)), so conditions are
naturally sparse: applying one means "fresh run from the `init_*` defaults,
with these overrides." The alternative base — the stopped sim's current
stores — was rejected: it makes the result depend on run history, exactly the
hidden input the trace-header discipline exists to kill. Warm restart needs no
second semantics: a `capture` service reads the current stores **and root
slots** back *as a condition value* (capture → tweak → apply) — slot coverage
is what makes the captured condition total, hence re-applicable under [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator) —
the same gather the trace header already needs. One mechanism, two uses.

**The mirror-tree spelling was rejected** (nested NamedTuples shaped like the
assembly): the same information, but a second spelling of structure that must
be zipped against the real one, ragged under partial specification, and
outside the path vocabulary that `child_connections`, diagnostics and the `Build`'s
provenance tables already speak.

**Doctrine.** This does not reopen [§13.5](#135-termination-is-a-state-not-an-exception)'s observation-by-path rejection:
that was *runtime* coupling — a root-authored predicate reaching through
generic seams it does not own, breaking on substitution. A condition is a
*design-time statement about a concrete build*, authored in the same register
as `child_connections` (which also speaks paths, about children its author owns).
[§14.2](#142-fragment-composition-locality-without-schema)'s composition law makes the parallel exact.

**Pre-sweep doctrine.** Condition writes precede the first sweep by
definition; a would-be init value that depends on swept outputs is either
analytically known to the caller (trim's `α_filt = α_a`: α is a *decision
variable* — the value is known above, not computable below) or an equilibrium
constraint, i.e. a job for the trim service, not for init.
"Caller-computable" reaches past closed-form knowledge to **environment
queries**: a condition needing one constructs the same handle the sweep will
produce — through [§4.4](#44-function-valued-signals-environment-access)'s value-level constructor, applied to the same values
its `baseline` writes into the environment component's slots, or, in a rig
where the handle itself is a root slot value, by simply holding the value it
wrote there — and then calls the same query function the consuming component
calls. One implementation of the field math, evaluated one level up: no
pre-sweep, no new mechanism. And where closed-form enforcement of a target is
not wanted at all, the second escape already covers the case — promote the
eliminated state coordinates to decision variables and enforce the targets as
residuals on swept outputs, which needs no environment access at condition time
whatsoever.

### 14.2 Fragment composition: locality without schema

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
is known ([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)'s pre-sweep doctrine).

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
struct Fragment{X,M,S}  x::X; m::M; slots::S  end          #self-vocabulary payloads; no paths
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
the GUI ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract), [§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)).

**The locality law** is [§6.1](#61-connections-and-hierarchy)'s, third instance (child connections, computed boundary connections,
now conditions): each level speaks its own fields, its declared children's
names, and its own faces; delegation by dispatch at every genericity seam;
deep `at` paths legitimate exactly where deep connections are — within an
owned concrete subtree. Absolute paths exist only in the flattened entry
list, a *compiled derivative* of the composition, as slot offsets are of
`child_connections`. Substituting a component invalidates precisely the fragments
its owner shipped, nothing else. Enforcement status is also [§6.1](#61-connections-and-hierarchy)'s: convention
(ownership is a fact about who maintains the code, not build-visible),
available and idiomatic rather than machine-checked. `fragment`/`at`/`merge`
are [§13.7](#137-tooling-consequences-provenance-and-the-component-library) standard-library material — ordinary artifacts, no privileges.
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
(or `at(prefix, fragment(...))`) and merge nodes with nodes. The rejection
carries a kind like every other (`ConditionNodeMisuse`, [Appendix C](#appendix-c-the-diagnostic-kind-set) — the
offending argument's type, the node kinds in hand): it is raised at
composition time, before any resolution pass or provenance chain exists,
which is why it is its own kind and not a `ConditionResolution` sub-kind
([§14.3](#143-resolution-flatten-validate-compile-once)). The
explicit, *ordered* layering spelling — `override` — was deferred here and
admitted one sub-topic later, when slot totality produced its use case
([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)).

### 14.3 Resolution: flatten, validate, compile once

Resolution takes the root node plus a `Build`. Flattening is the only place
path strings are ever concatenated — a trivial recursion with a path
accumulator that also records each entry's **tree position** (its
`getfield`/`getindex` step tuple). The collecting pass then checks each flat
entry: path resolves (did-you-mean over children), field declared in the
target's `init_x`/`init_m`, value type convertible to the declared
leaf type, slot faces reach root slots, no duplicate
`(path, store, field)`. The `Build` supplies two lookup families: **schema**
(the evaluated declarations: may you write this field, at what leaf type —
the authority) and **layout** (`x` backing ranges, store indices for `m` and for discrete `x`, slot
indices from the activation; face chains from Stratum A — the destination).

A valid list compiles to a plan: per leaf, a `Getter{P}` lens (the position
tuple lifted to a type parameter — type-stable navigation of the fixed tree
type), a destination offset, and a **converter baked now**, *selected per
leaf from that leaf's type in the resolved shape* — leaf types are shape
facts, carried by the tree type along with the full nesting and every field
name ([§14.4](#144-two-application-registers-over-one-plan)), so selection consults no runtime fact and stays a
resolution-time bake. Two cases. A leaf **already at the activation's scalar
type** is decision-descended — under a `Dual`-seeded evaluation of a
type-stable `trim_condition(d)` every decision-dependent leaf is `Dual`-typed
in the shape ([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)) — and takes the type's ordinary `convert`/constructor
methods *at that eltype*: an authored `RQuat` of `Dual`s → the `SVector{4}`
state leaf at `Dual`, partials flowing through untouched, which is what makes
the seeded decisions reach the sweep at all; at the nominal activation the
same rule is the ordinary `Float64` conversion (an authored `RQuat` value →
the `SVector{4}` state leaf it initializes). A plain **`Float64` leaf against
a non-nominal activation's scratch** is a held constant and takes the
`Float64 → Dual` zero-partial embedding — semantically exact there, and
exactly there: "held at the operating point" *is* zero partials, the whole of
a linearization operating-point condition, authored decision-free
([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)). The
selection is a one-time boundary decision that leaves [§12.5](#125-the-always-on-conformance-check)'s
nominal exact-match doctrine for table cells untouched. Converters run here
and in `capture`'s gather ([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)) — the write paths — never on state views
([§7.1](#71-continuous-state-structured-immutable-flat-backing)). Overlay
partiality for the `m` and discrete-`x` stores is baked the same way: the writer holds
`merge(init_m_defaults, overlay)` with the base resolved at compile time
([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)'s fork).

### 14.4 Two application registers over one plan

Execution is where the paradigm-change tax was feared, and where it
vanishes: all string work, validation and addressing are functions of the
condition's *shape*, and every hot path holds the shape fixed while varying
values. Hence resolve-once/execute-many, with two registers over one plan:

- **Specialized `apply!`** — for services that iterate (trim's per-evaluation
  write, linearization's seeding). Unrolled stores through the baked lenses
  and converters: the same machine operations as today's in-place writes,
  zero-alloc, no strings, no dispatch. The per-iteration shape check is
  [§12.5](#125-the-always-on-conformance-check)'s mechanism transferred: the tree type (which carries the full
  nesting, every field name and leaf type) is proven by dispatch, and a `===`
  sweep over the interned path literals closes the remainder — pointer
  compares that fold to nothing in the all-literal case, with shape drift a
  structured error, not silent corruption. Cost: Julia codegen of ~10–50 ms
  *once per condition shape* — noise against the model's own first-sweep
  warmup (seconds) and the 10³–10⁴ optimizer evaluations it amortizes over.
- **Dynamic walk** — for one-shot init. The same validated entry list,
  executed by runtime dispatch per write: microseconds total, allocation
  permitted ([§7.5](#75-allocation-policy-a-scoped-invariant) — the stopped-sim path was never under the zero-alloc
  regime), and no per-shape codegen: fifty structurally different scripted
  conditions cost fifty walks, not fifty compiles.

Which register a service uses is internal, never user-facing API.

**The read-selector family is closed**: `get_state(path, field[, i])`,
`get_deriv(path, field[, i])`, `get_output(path, field[, i])`,
`get_slot(face)`, `get_face(name)` — one
address space for every reader of the model. The names carry a deliberate
`get_` prefix: a selector is a *deferred read* — a value describing the read
the compiled gather will perform — and the prefix both names that action and
keeps five short common nouns out of the namespace user declarations share
with domain code. There is no selector for a component's private
intermediates, and cannot be: they are values in flight, not cells anything
could address ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)), so a reader that wants one is a component that should
declare it an output ([§11.3](#113-visibility-the-contract-is-the-interface)). `get_face` addresses a root-exported output face —
[§9.2](#92-outbound-snapshot-publication)'s *integration* register, previously recommended but unspellable in
the family.

**A selector resolves against a source, before any client policy applies.**
The table selectors — `get_output`, `get_slot`, `get_face` —
resolve against a *table source*: a boundary snapshot, or the scratch tables
a service evaluation instantiates ([§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)) — the axis separates
table-borne values from store-borne ones, not snapshots from services. The
store selectors — `get_state`,
`get_deriv` — resolve only against live stores, which only stopped-sim
service evaluations, `capture`, and post-run inspection of the live stores
([§9.2](#92-outbound-snapshot-publication)'s replay-to-inspect) ever hold: the snapshot deliberately carries no
state stores ([§9.2](#92-outbound-snapshot-publication)), and `ẋ` buffers are integrator scratch, not
boundary-consistent objects outside a service evaluation. A snapshot-bound
reader naming a store selector is therefore a resolution error at attach
(`ReadBindingUnresolved`), in the didactic register, with [§9.2](#92-outbound-snapshot-publication)'s honest
remedy: declare the field public and read the auto-published port. Client
policy rides on top — row 83's registers restated as a resolver property:

- **Load-bearing services** (trim's `reads`, linearization's taps) speak
  the contract: `get_state`/`get_deriv`/`get_output`/`get_slot`/`get_face`,
  within the scopes [§6.1](#61-connections-and-hierarchy)'s locality law and [§14.2](#142-fragment-composition-locality-without-schema)'s fragment scoping own.
  `get_face` is the set's seam-crossing member: it resolves through export
  chains exactly as [§14.9](#149-mounting-problems-as-relocatable-values)'s mounting resolves slot faces — the read
  side mirroring the write side — so an equilibrium equation reaching
  behind a generically-held child binds the curated face register instead
  of a path the locality law forbids. A service evaluation needing a private
  intermediate has one remedy, and it is the same at every register: the
  component exports it ([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)).
- **Diagnostic readers** (output-device bindings, GUI panels, log
  inspection) admit the whole family, within the source rule: deep paths
  and `get_face` names alike, with the store selectors
  reaching only the diagnostic clients that actually hold stores (`capture`,
  post-run inspection) — a snapshot-bound reader is barred from them by
  source, not by client.
- **`stop_on` is not a family client.** It names root-exported `Bool`
  output faces, period ([§13.5](#135-termination-is-a-state-not-an-exception), row 60): termination is run policy against
  the root contract, and no path selector reaches it.

**Compiled readers are the gather twin** over this family
and the layout tables: trim's cost read (`ẋ` and output fields), linearization's
Jacobian gather, and `capture`'s full-store readback are one primitive run in
reverse — one machinery, both directions, in the `Build`'s client kit. The
per-iteration ledger for trim: user fragment math (stack-only, the domain
computations unchanged from today) + leaf stores + folded shape check +
sweep — the sweep dominates, exactly as `f_ode!` does today. `apply!` ends at
established stores; making the model *coherent* is boundary zero, [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions).

### 14.5 Boundary zero: an ordinary boundary with authored incoming transitions

After `apply!` establishes stores at `t₀` — and the trace header captures
them, together with the slot values, *before anything below runs* ([§9.5](#95-inbound-the-input-trace)'s
capture placement; a post-sequence capture would hand replay
already-transitioned state) — the init service completes the
[§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) macro-sequence with an empty integrate — project → [sweep → guards →
handlers]\* → due `g` updates → first snapshot — and that
parity is exact, not approximate. Piece by piece:

- **Project runs.** Authored `x` can sit off-manifold (a hand-assembled
  quaternion ulps off unit norm; a condition writing part of a constrained
  block against fresh defaults). Projection after condition writes is the
  same position it holds after any other `x` mutation, and costs nothing
  when the state is already clean.
- **The sweep runs with every tick due.** `t₀` is a grid point of every
  divisor, so all discrete output stages are gated in, publishing from the
  authored `x` — necessarily, since no earlier tick exists for a ZOH to
  hold. The `t₀` snapshot carries a fully populated table.
- **Events run.** A condition landing a guard predicate in holding territory
  (an authored stall flag, a strut authored into contact) fires visibly at `t₀` rather
  than one step later — grounded by the prior rule ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)):
  boundary zero establishes every guard prior as not-holding. Suppression
  would delay the identical firings while hiding the diagnostic that the
  authored condition was not quiescent — [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)'s stage-on-interaction lesson
  (insurance that masks invariant violations is anti-diagnostic). The header
  records the *resolved pre-sequence* stores and slots ([§9.5](#95-inbound-the-input-trace)), so replay
  re-executes boundary zero from the same starting point and whatever fires
  at `t₀` fires again identically ([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)) — the firings are recomputed,
  never recorded. (A `stop_on` face already `true` is a
  different category: nothing *fires* — the face simply reads `true` in the
  published `t₀` snapshot and the loop reacts, [§13.5](#135-termination-is-a-state-not-an-exception).)
- **Due `g` updates run.** This follows from an interval-alignment fact that
  is easy to mis-picture and is hereby a taught contract, sibling to [§15.5](#155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary)'s
  boundary-sampling line: **a boundary's `g` is the *outgoing* transition** —
  at tick `t_k` it consumes the completed boundary's samples and produces
  `x_{k+1}`, the value the next tick reads; the transition that carried `x`
  *into* `t_k` ran at `t_{k-1}`. Boundary zero is missing its incoming
  transitions on *both* tiers, and both are replaced by authorship: no
  integration over `[t_{-1}, t_0]` produced a continuous leaf's `x(0)`, and no
  `g` at `t_{-1}`
  produced a discrete leaf's `x(0)` — the condition did. The outgoing work all runs: `g` at
  `t₀` computes `x(1)` from `t₀` samples, its only opportunity — `x(1)` must
  sit in the store before `t₁`'s gated stages read it. Skipping `g` at
  boundary zero would not preserve the authored `x(0)` — that is already
  preserved, published in the `t₀` snapshot — it would *delete the `t₀`
  sample from the discrete dynamics*: an accumulator
  $x_{k+1} = x_k + \Delta t \, e_k$ authored with $x_0 = 0$ under nonzero $e(t_0)$
  would first integrate $e(t_1)$, the whole sampled-data lattice one period
  late. (The continuous-tier analogue of `g`-at-`t₀` is not the empty incoming
  integrate but the first *outgoing* one, $t_0 \to t_0 + h$: both authored values
  are the published initial conditions of their outgoing transitions.)
- **`t₀` is an init-service argument** (default `0.0`), never a condition
  entry — time is not a store of any component, and the harmonic grid
  anchors at whatever `t₀` boundary zero runs at. Both init-service entry
  points carry it: `init!(sim, condition; t0)` and `trim!`'s commit
  (`trim!(sim, problem; baseline, t0, backend)`, [§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)), same default.
  Conditions are time-free;
  `capture` returns condition and time separately for resume-at-time —
  the returned `t` passed back as `t0`.
- **Trim is untouched by all of this.** Optimizer iterations are raw
  write → sweep → read cycles on the activation — no boundaries, no events,
  no `g`; only the committed solution executes boundary zero. A guard firing
  at commit is a wanted failure signal: today's hand-written trim asserts
  (`!stall`, no weight-on-wheels, `ω > ω_idle`) become the model's own event
  logic, surfaced through the ordinary machinery instead of `@assert`. And
  the signal has a channel: boundary zero reports the fired set on the
  `TrimReport` and raises `TrimCommitEvents` when it is non-empty
  ([§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report), [Appendix C](#appendix-c-the-diagnostic-kind-set)) — a handler that fires at commit moves the
  committed stores off the solved point, and saying nothing would be
  warn-but-assign relocated.
  A commit-fired handler is not the only mover, and the second one is
  unconditional: boundary zero's *first* act is `project`, so the committed
  `x` is `project(x*)`, not the solver's `x*` — an attitude quaternion
  renormalized by a few ulps is the canonical case. Legitimate, wanted, and
  usually invisible in the residuals; but the point the stores sit at is no
  longer the point the verdict was read at, and the doctrine is the same for
  both movers. So is the remedy: the `TrimReport` carries the
  committed-state residuals beside the solved-point ones, and a converged
  solve whose committed-state residuals fail the box test raises
  `TrimCommitResiduals` ([§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)) — the move made visible rather than silent.

### 14.6 Slot totality: the missing-value error and the `override` combinator

Slots are the one initialized datum without declared defaults — [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s
bare-types decision, upheld here: a default inside a face declaration would
scatter condition data into the wiring contract and recreate the
competing-defaults problem that [§14.2](#142-fragment-composition-locality-without-schema) killed for `initialize` specs. So a
slot's only source before boundary zero is the condition, and three
consequences follow.

**Totality is a precondition of starting, checked by the service.** A
condition value is legitimately partial (fragments compose; trim iterations
write subsets; capture-then-tweak patches leaves) — "every root slot
covered" is not a property of conditions but of *every application that
establishes a complete world over virgin stores*. That principle, not an
enumeration, names the sites: `init!`, trim's setup application to freshly
allocated scratch stores, and trim's commit through boundary zero
([§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)) — one class, one mechanism, one kind. Each compares the resolved
plan's slot coverage against the `Build`'s `input_faces` before writing
anything; a shortfall is one collected, declaration-ordered diagnostic
(`UninitializedSlots`, a [§13.2](#132-diagnostics-structured-values-one-carrier-exception) kind) naming every uncovered face.
Coverage is a *plan-level fact* — both operands are resolution-time data —
so the check is one comparison and runs before any evaluation, not merely
before any write. Pre-write means all-or-nothing: a rejected init leaves
the sim exactly as it was, the same posture as failed trim.

**The probe-value barrier is structural.** [§12.3](#123-probing-and-input-synthesis)'s `probe_value` synthesis
(zero/false/first-enum/`T()`) exists so build-time probes can exercise code
with fabricated inputs; a fabricated zero is a fine probe input and a
terrible flight condition (a silently zeroed `mixture` kills the engine and
sends the user debugging aerodynamics). The services path simply contains
no call to it: a slot gets a condition value or the application errors —
no third branch. Replay likewise never synthesizes: the trace header
records every slot value, and with totality enforced its slot capture is
complete by construction ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s requirement discharged).

**Baselines are aircraft-shipped condition functions, layered by
`override`.** Nobody hand-writes ~20 slot values per script; today
`SystemsInitializer`'s `@kwdef` defaults carry that load, and their
successor is ordinary user math in one authoritative home —
`ready_for_taxi(ac)`, `cold_and_dark(ac)` — returning full-coverage
conditions. But "baseline plus tweaks" collides with [§14.2](#142-fragment-composition-locality-without-schema)'s
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

### 14.7 The trim problem: NamedTuple decisions, declared reads, named residuals

What the aircraft author ships, piece by piece against today's `c172.jl`:

- **Decision variables, initial guess and box bounds are plain, same-*named*,
  all-`Float64` NamedTuples.** The `AbstractTrimState{N}`/`FieldVector`
  supertype dies — its only job was vectorization, which is the service's
  (pack/unpack by field order, that order being the `guess` NamedTuple's own —
  [§12.5](#125-the-always-on-conformance-check)'s rule at this seam: **the names are the pairing, order carries
  no semantics**). `lower` and `upper` are checked at setup for key-set equality
  with `guess` and `Float64` fields, then canonicalized to `guess`'s field order
  by the same type-level reorder (`NamedTuple{keys(guess)}(lower)`), so a
  permuted bound spelling is a non-event rather than `α`'s bound silently
  applied to `throttle`; a key-set or field-type mismatch is
  `TrimProblemInvalid`. Guess, bounds and
  the returned solution share one spelling; `Base.merge(guess, (throttle =
  0.3,))` is free warm-start tweaking; an author who wants a documented
  `@kwdef` struct keeps it privately and converts.
- **`TrimParameters` stays a plain user struct** the framework never sees;
  the assignment is the pure `trim_condition(ac, params, d)` fragment-tree
  function ([§14.2](#142-fragment-composition-locality-without-schema)), applied per iteration by the compiled plan ([§14.4](#144-two-application-registers-over-one-plan)).
- **The read side is declared, then compiled**: `reads(name = get_state(path,
  field) | get_deriv(path, field) | get_output(path, field) | get_slot(face) |
  get_face(name), ...)` — [§14.4](#144-two-application-registers-over-one-plan)'s load-bearing set. `get_state` and
  `get_deriv` address a declared state field and its derivative (validated
  against `init_x`), `get_output` a declared output port (validated against
  `output_types`), `get_slot` and `get_face` a root input and output face
  (validated against the root face lists). The path selectors reach only
  through [§6.1](#61-connections-and-hierarchy)'s locality scopes; an equilibrium equation crossing a
  generic seam reads a face. A component's private intermediates are
  not addressable at all ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), [§14.4](#144-two-application-registers-over-one-plan)) —
  a trim evaluation needing one is a signal the component should export it —
  and a derivative wanted across a contract boundary takes the same remedy:
  publish it as an ordinary output port computed in `h_xu` ([§7.4](#74-the-fused-evaluation-lineage-prior-art-and-how-we-got-here) step 2's
  one-line binding, made contract), leaving `get_deriv` scoped to owned
  concrete subtrees.
  The compiled reader ([§14.4](#144-two-application-registers-over-one-plan)'s gather twin) fills a stack-only NamedTuple
  per evaluation.
- **The user supplies a residual *system*, not a scalar cost** — authored as a
  NamedTuple of named equations, packed to the solver's vector by field order —
  `tolerances`' field order here, the declared side again, the residual return
  canonicalized to it (the decisions rule, symmetric on both ends of the
  seam: names pair, order never does). The
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
  the `T`-generic assignment, sweep and `f`; the seeds survive the condition
  write boundary because [§14.3](#143-resolution-flatten-validate-compile-once) selects the baked converter per leaf from the
  shape — a decision-descended leaf is `Dual`-typed there and takes the
  structural conversion, the zero-partial embedding staying on the held
  `Float64` leaves — [§12.6](#126-stopped-sim-services-as-stratum-c-clients)'s "open option," now
  the *default*: nonlinear least squares on $r(d)$ with exact AD Jacobians
  (trust-region/Levenberg–Marquardt register), quadratic convergence
  (~5–15 evaluations), per-residual physical tolerances — the convergence
  verdict itself service-owned and backend-independent
  ([§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)) —
  and failure reports naming the unbalanced equations with magnitudes. Non-squareness
  degrades gracefully (redundant actuation → weighted/minimum-norm LS;
  infeasible demands → converged nonzero residual identifying the
  impossible balance), and $\partial r / \partial d$ at the solution is free flight-physics
  data (control effectiveness) cross-checking linearization. The
  derivative-free scalar path survives as the fallback: the service squares
  *and normalizes* the residuals against the tolerances
  ($\sum (r_i/\mathit{tol}_i)^2$ at `stopval = 1`,
  [§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)),
  which is where the hand-scaled absolute threshold is repaired — today's
  algorithm as the degenerate case.
- **Recorded, not built**: closed-loop sampled-data trim (append
  $g(x) - x = 0$ residuals via a nondestructive scratch evaluation of `g` —
  structurally impossible under FlightCore's mutating `f_disc!`), and
  on-ground static equilibrium (strut compressions and attitude against
  gear forces) as simply another problem value over the same service.

**The field set, normative** — [§14.9](#149-mounting-problems-as-relocatable-values)'s lift is field-by-field, so this list
is closed: `guess`, `lower`, `upper` (same-named all-`Float64` NamedTuples),
`condition` (the condition-valued function over decisions), `reads` (the
declared read set), `residuals` (the residual function), and `tolerances` —
an all-`Float64` NamedTuple, same-named as the residual function's return,
carried *in the problem* because a relocated problem must carry its own
convergence test (`at` passes it through untouched). The residual signature:

```julia
residuals(reads::NamedTuple, d::NamedTuple) → NamedTuple
```

The rule: **what the solver varies is passed; what is fixed per problem is
closed over.** The gathered reads and the decision NamedTuple arrive as
arguments — `d` is the one value that *cannot* be closed over — while
`TrimParameters` stays behind the closure, exactly as `condition` already
holds it (the framework never sees it). Being user-shaped, that record is also
where any environment handles the condition math needs conventionally ride
([§4.4](#44-function-valued-signals-environment-access)'s value-level constructor, [§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)'s pre-sweep doctrine): the problem
*receives* the environment, and never writes it ([§14.9](#149-mounting-problems-as-relocatable-values)).
The returned NamedTuple's names are
the equation names the report and the failure messages use; the service packs
residuals and tolerances by field order — `tolerances`' order, canonical for
the r-side, each residual return reordered to it (`NamedTuple{keys(tolerances)}(r)`)
before packing, so an equation list spelled in a different order inside the
lambda pairs correctly and costs nothing. Names and types are checked at setup:
the guess evaluation the service performs anyway observes the residual key set,
and a `tolerances` key-set mismatch — or any field-type disagreement, on either
side of the seam — is `TrimProblemInvalid` ([Appendix C](#appendix-c-the-diagnostic-kind-set)), with the offending
field and the names or types in hand. Order is never a mismatch.
Worked, the C172 cruise case reduced to
its three-equation core (the real problem is the same shape with the full
7-variable search; the θ-constraint elimination survives inside
`trim_condition`):

```julia
cruise = TrimProblem(
    guess      = (throttle = 0.6, elevator = -0.05, α = 0.06),
    lower      = (throttle = 0.0, elevator = -1.0,  α = -0.1),
    upper      = (throttle = 1.0, elevator =  1.0,  α =  0.3),
    condition  = d -> trim_condition(ac, params, d),    #params closed over (§14.2)
    reads      = reads(v̇_b = get_deriv("vehicle/dynamics", :v_eb_b),
                       ω̇_b = get_deriv("vehicle/dynamics", :ω_eb_b)),
    residuals  = (r, d) -> (axial_force  = r.v̇_b[1],
                            normal_force = r.v̇_b[3],
                            pitch_moment = r.ω̇_b[2]),
    tolerances = (axial_force = 1e-3, normal_force = 1e-3, pitch_moment = 1e-4))
```

### 14.8 The trim service: solver seam, scratch stores, commit and report

**The backend seam.** `trim!(sim, problem; baseline, t0 = 0.0, backend =
LevenbergMarquardt())`. The default is an in-house dense
Levenberg–Marquardt: for decision dimensions ~10 with exact Jacobians, the
core (damping loop, small linear solve, convergence test) is ~100 lines —
the [§8.2](#82-the-stepper-seam) stepper precedent exactly (tiny needed core vs. heavy dependency),
sharpened by the fact that [§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)'s per-residual physical tolerances are a
convergence test no external package spells natively — which is precisely
why the *service*, not the backend, applies it. The backend contract is a
**pinned signature**, value-passed — one required method per backend:

```julia
solve(backend, eval!, d0, lower, upper, tol) -> (; d, status, nevals, niters)
```

`eval!(r, J, d)` is in-place and always fills `r`, the residual vector packed
in `tolerances`' field order. It fills `J` **iff `J !== nothing`**: the
Jacobian is requested by argument, so a Jacobian-free backend simply always
passes `nothing` and there is still exactly one evaluation method to
implement. `d0`, `lower` and `upper` are packed `Vector{Float64}` in `guess`'s
field order ([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals): the declared side is canonical on both ends, and the
service has canonicalized before packing), with `±Inf` meaning unbounded — a
backend that ignores bounds ignores two vectors, not a missing argument.
`tol` is a `Vector{Float64}` in `tolerances`' field order: data the backend
*may* stop on (that is the service's per-register translation, below) and
decisive of nothing. Returned: `d`, the solution; `status::Symbol`, from a
deliberately **open** set, recorded verbatim in the report — the verdict is
the service's, so a closed enum would buy nothing and would launder foreign
solvers' vocabularies back into per-backend meaning; and `nevals`/`niters`,
diagnostic counts. The name `solve` is subject to the [§16](#16-open-axes) naming audit like
every other API spelling. The backend sees vectors and never names, so the
solution it returns unpacks by the same order it was given.

**The convergence verdict is the service's, uniformly.** `converged` means
`all(abs.(rᵢ) .≤ tolᵢ)` — [§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)'s per-residual box test in its own
physical units — evaluated **by the service at the backend's returned
point**, one residual evaluation, noise against the solve that produced it.
That verdict, and nothing else, gates the commit and fills
`TrimReport.converged`; the backend's returned `status` is recorded in the
report as diagnostic data and is authoritative over nothing.

**The tolerance translation is the service's too**, per register. In the
least-squares register the tolerances *are* the stopping criterion: they
feed the per-residual test directly, LM's damping loop testing exactly what
the service will re-test. For the derivative-free scalar fallback —
`NLoptBackend(:LN_BOBYQA)` in a package extension, passing `nothing` for the
Jacobian, today's algorithm one keyword away, with the framework core carrying zero
optimizer dependencies — the service squares *and normalizes*: the objective
minimized is $\sum_i (r_i/\mathit{tol}_i)^2$ with `stopval = 1`,
dimensionless where FlightCore's threshold was hand-scaled and absolute
([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)), and a well-scaled valley where a raw $\|r\|^2$ sums forces
against moments. The two criteria cannot disagree in the dangerous
direction: $\sum_i (r_i/\mathit{tol}_i)^2 \le 1$ implies every
$(r_i/\mathit{tol}_i)^2 \le 1$, so the `stopval` sphere is *inscribed* in
the tolerance box and a fallback stopping at `stopval` necessarily passes
the service's box test. The converse disagreement — a backend stopping
early and reporting an optimistic status — is caught by the re-check, which
remains the single authority. What is *not* claimed is point identity:
different backends may land on different solutions, an algorithmic
difference and a legitimate one. What is eliminated is per-backend meanings
of `converged`.

Box bounds are honored by step
projection, and a decision variable saturated *at the solution* is flagged
in the report ("converged with `elevator` at its upper bound" — the classic
CG-limit diagnostic, today inferable only from mysterious residuals).

**Scratch stores, stated without type luck.** Every `trim!` invocation
instantiates a fresh working store set — `x` backing, `m` and discrete-`x` stores, slot
and signal tables, derivative buffer — from the activation's *layout*: the
layout is the reusable compiled artifact, the buffers are per-invocation
and die with the call (stopped-sim allocation, [§7.5](#75-allocation-policy-a-scoped-invariant)). The `Dual` backend's
buffers being un-aliasable by type is defense in depth, not the mechanism —
a `Float64` backend (NLopt) gets equally fresh buffers. The invariant is
backend-independent: **the simulation's authoritative stores have exactly
one writer, the commit through boundary zero.** Setup applies
`override(baseline, condition(guess))` to the scratch set once, its full
coverage *checked here* — [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)'s comparison of the resolved plan against
the `Build`'s `input_faces`, one plan-level comparison before the first
evaluation — so sweeps see a complete world. Raw instantiation is sound
exactly because of that check: every slot is written before any read, and
an incomplete `baseline` is one declaration-ordered `UninitializedSlots`
at setup rather than a whole solve against undefined cells. Since the
commit applies the same composite over the same `baseline`
([§14.9](#149-mounting-problems-as-relocatable-values)'s `override(baseline, condition(d*))`), its coverage is setup's:
commit's totality check is structurally unfailable through the trim path
and stands as the shared `init!`-boundary defense — a converged solve is
always committable, so `TrimReport` carries no committed flag and the
no-throw doctrine needs no exception. Iterations rewrite only the
problem's write-set via the compiled plan; an LM evaluation is one
Dual-seeded sweep yielding `r` (value parts) and `J` (partials) together.
No convergence — the service's box test failing at the returned point,
whatever status the backend attached to it — → no commit → the sim is
bit-for-bit untouched, including
"never initialized" — today's warn-but-assign is structurally impossible. The
same structure covers an interrupt: Ctrl-C during a long solve unwinds an
ordinary Julia call operating on per-invocation scratch stores, no commit has
happened, and the simulation is bit-for-bit untouched exactly as a
non-converged solve leaves it — the services need no counterpart to the loop's
boundary masking ([§10.4](#104-shutdown-protocol)).

**The commit, in full.** The committed solution is applied as an `init!`
in every respect: `override(baseline, solution)` through boundary zero —
[§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator) pre-write slot totality, [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)'s sequence, guards at commit — with
the harmonic grid anchored at `trim!`'s own `t0` argument (default `0.0`,
the same default as `init!`'s: one rule for both init-service entry
points), and the recorders cleared exactly as [§10.6](#106-run-lifecycle-and-partial-advance) states for
`init!` — trace, log, and any batches still in staging cells. A fresh
recording starting at its own anchor is the unattended register's natural
shape; fly-then-retrim keeps continuity explicitly — `capture` returns
`(condition, t)` separately for exactly this ([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)), so the resumed
spelling is `trim!(sim, problem; baseline = c, t0 = t)`.

**The report, not an exception.** `trim!` returns a structured
`TrimReport`: the `converged` flag — the service's own box test at the
returned point, never a backend's opinion — solution NamedTuple
(guess-shaped — warm-startable), the **solved-point residuals** with their
tolerances (the very numbers the verdict is read off, gathered at the
backend's returned point), the **committed-state residuals** (the same
residuals re-gathered from the boundary-zero world after the commit — nearly
free, since that boundary's sweep has already run and the residuals'
declared reads need only gather from it: the numbers describing the state
the simulation is actually *in*, which is the point a `capture`-defaulted
`linearize` reads), the backend's returned status together
with its iteration/evaluation counts (diagnostic throughout: informative
about *how* the solve went, decisive about nothing),
saturated-bounds list, and the commit's fired events (component
paths and event names, empty when boundary zero ran quiet — [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions); a
non-empty set also raises `TrimCommitEvents`, [Appendix C](#appendix-c-the-diagnostic-kind-set): the committed
stores then sit at the post-handler point, not the reported solution, and a
`capture`-defaulted `linearize` ([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)) reads that point).
The two residual sets are what make the moved point auditable: a converged
solve whose *committed-state* residuals violate the box test raises
`TrimCommitResiduals` ([Appendix C](#appendix-c-the-diagnostic-kind-set)), naming the offending residuals with
their committed values and tolerances — the move (`project`, or a
commit-fired handler, [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)) surfaced rather than left silent. The verdict
itself is not re-litigated: it gated the commit, at the solved point, and
row 150's numbers stand as reported. Non-convergence never throws — it is an
expected *outcome* (envelope-sweep data: hitting the infeasible edge is
information), per [§13](#13-error-discipline)'s exceptions-are-broken-machinery line; a malformed
problem is a `BuildError`-class failure at setup (`TrimProblemInvalid`,
[Appendix C](#appendix-c-the-diagnostic-kind-set) — a guess/bounds key-set or field-type disagreement, an
unknown `reads` selector, a
`tolerances`/residual key-set mismatch observed at the setup guess evaluation:
the offending field with the names or types in hand, collected, mirroring
linearization's `TapResolution`); a permuted spelling is none of these
([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)).
An *empty* one is none of them either: `TrimProblem(guess = (;), …)` is
legal, not `TrimProblemInvalid`. With zero decision variables the solver is
bypassed outright — nothing to pack, no seeded activation, no backend call —
and the service simply evaluates the residuals once at the baseline, the
ordinary box test deciding `converged` and the commit as usual. The
degenerate problem is the "is this operating point an equilibrium?" probe:
evaluate this condition's equations and report, useful in its own right and
free.

**The AD obligation, scoped.** The default formulation requires `Dual`
genericity of exactly: the continuous output-stage chains and `f`, plus the
user's assignment and residual math. The discrete tier's stages and `g`, and the
event system's guards and handlers, never see a `Dual` — frozen constants with
zero partials, semantically exact ([§11.2](#112-the-declaration-inventory)). This is *not a new obligation*: it is the same activation
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

### 14.9 Mounting: problems as relocatable values

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
either condition-producing (path-relative) or path-free — [§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)'s rule
that residual math sees only the gathered NamedTuple pays off here:

```julia
at(prefix::String, p::TrimProblem) = TrimProblem(
    p.guess, p.lower, p.upper,                 #path-free: pass through
    d -> at(prefix, p.condition(d)),           #post-compose: wrap each returned tree
    at(prefix, p.reads),                       #reads are inert selector data: same Scoped node
    p.residuals, p.tolerances)                 #path-free: pass through
```

Resolution then needs nothing new: the flattening accumulator of [§14.3](#143-resolution-flatten-validate-compile-once)
enters the `Scoped` wrapper and prefixes every entry
(`"vehicle/dynamics"` → `"wing/vehicle/dynamics"`); slot entries authored
in the aircraft's face vocabulary resolve through the export chain *from
the mount point* (`throttle` at `"wing"` → root slot `"wing.throttle"`).
An unexported face fails resolution by name — correctly: an internally
wired input (a scenario component driving the wingman's throttle) is
untrimmable from outside, and the build says so. The service compiles the
scoped condition and reads and runs the identical loop — it never knows
where its paths are mounted.

**A problem never authors the environment.** The environment — a sibling
component's slots in a full world, a handle-valued root slot in a thin rig —
is the world's and the `baseline`'s business; the problem *receives* its
handles through the user parameter record ([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)) and only queries them. The
reason is the resolution rule just stated: a condition entry naming a wired
input fails by name, correctly — so a problem writing an environment face
would be applicable only to those rigs where that face happens to be
unconnected, and the relocatability this section exists to guarantee would be
lost.

**The world wrapper dissolves.** Today's `f_init!(::Model{<:SimpleWorld})`
(initialize environment, then call the aircraft's trim) has no successor
method: the environment, the other aircraft and all slots are covered by
the `baseline` condition ([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)), applied once at setup; the commit is
`override(baseline, at(mount, condition(d*)))`. Method nesting became value
layering.

**"Aircraft as root" is a thin world.** By default the aircraft is not
literally the root — its environment inputs ([§4.4](#44-function-valued-signals-environment-access) function-valued
signals) are wired from provider components — so design tasks use a shipped rig,
`design_world(ac)` = aircraft + `SimpleAtmosphere(wind = NoWind())` +
`HorizontalTerrain`: today's ad-hoc models inside `linearize` promoted to
a named artifact. One register: the "root" case is the shallowest world,
the trim problem mounts at `"aircraft"` like anywhere else. Leaving an
environment face *unconnected* is legal by construction, though: the face
becomes an ordinary root slot holding the handle **value**, written by the
`baseline` like any other slot — the test-rig register, the function-valued
sibling of a constant source, zero ceremony for a frozen environment. For
design tasks the shipped rig stays `design_world(ac)`, which keeps the
environment's tunables in the slot vocabulary that conditions, `capture`,
linearization's input surface and the trace header already speak.

**Swarm doctrine.** The service solves *one problem at a time*. Sequential
independent trims (trim lead, commit, trim wing against the committed
world) cover weak/one-way coupling; a joint trim is user-side value
composition — concatenate decision NamedTuples under prefixed names, merge
the scoped condition trees, stack the residuals. If joint trims become
routine, a `product(p₁ => "lead", p₂ => "wing")` helper belongs in the
[§13.7](#137-tooling-consequences-provenance-and-the-component-library) library; recorded, not built.

### 14.10 Linearization: tap selectors, one seeded pass, a pure query

**The tap declaration.** Today's per-aircraft `XStateSpace`/
`UStateSpace`/`YStateSpace` structs plus the `get_*_ss`/`assign_*_ss!`
shuttle methods (~150 lines of bookkeeping per variant) become three
selector lists drawn from [§14.4](#144-two-application-registers-over-one-plan)'s family (`get_state`/`get_slot`/`get_output`), with the
optional component index so a vector leaf yields *named scalars* — the NamedTuple
key is the label control design slices by:
`x = (p = get_state("vehicle/dynamics", :ω_eb_b, 1), θ = get_state("vehicle/kinematics", :θ_nb), …)`,
`u = (throttle_cmd = get_slot("throttle"), …)`,
`y = (EAS = get_output("vehicle/airflow", :EAS), …)`. Validated at resolution
against `init_x`/faces/`output_types` with did-you-mean errors, compiled to
offsets once, relocatable whole via `at(prefix, taps)` — the shuttle
layer's successor is the compiled writer/reader pair, and [§7.1](#71-continuous-state-structured-immutable-flat-backing)'s promised
`get_x_ss` deletion is discharged.

**The evaluation.** Instantiate a per-invocation scratch store set
([§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)'s mechanism verbatim), apply the operating-point condition, run
**one** Dual evaluation with a seed direction per `x`- and `u`-tap
entry (chunked internally). Value parts give `ẋ₀`, `y₀`; partials give
`A`, `B`, `C`, `D` simultaneously, exact to machine precision — replacing
four `FiniteDiff` jacobians, their step-size heuristics and ~4n perturbed
evaluations. Unseeded states sit constant at the operating point, and so do
unseeded slots (a root-slot cell follows the activation scalar by *evaluating*
its consuming `input_types` entry at that scalar, [§11.2](#112-the-declaration-inventory); the condition apply embeds
their `Float64` values as zero-partial constants); the discrete tier is
frozen with zero partials — precisely "linearize with the discrete state held" ([§11.2](#112-the-declaration-inventory)).
Differentiation participation is a per-invocation *seeding* fact for every
slot the schema leaves seedable — one register for `x` and slots alike — with
one declared exception now visible in the schema: a slot whose entry is
declared `Float64` is **declaredly unseedable**, its cell frozen at every
activation, so selecting it as a `B`-matrix tap is rejected at tap resolution
with the offending entry in hand rather than silently yielding a zero column
(row 167). Under fan-out the rejection names the **pinning consumer**, not
the face alone: a slot is unseedable whenever any one of its consumers demands
frozen ([§11.2](#112-the-declaration-inventory)'s meet, row 168), and the author's next move — promote that leaf
to a tolerant entry, or route the tap around it — depends on knowing which leaf
froze the slot.

**A pure query, and the shape of `capture`.** Linearization is the first
service with no commit and no boundary zero: scratch buffers only, nothing
becomes authoritative, and today's restore-the-trim dance (re-`assign!`
after `FiniteDiff` dirtied the model) has no successor. Default operating
point = the sim's current committed state via `capture(sim) → (condition,
t)` — the full gather of stores *and root slots* ([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)'s totality makes slot
coverage mandatory for capture → apply): after a
`trim!` commit, `linearize(sim, taps)` is about the trim point with
nothing re-specified; an `about = <condition>` keyword linearizes anywhere
else without touching the sim.

**The returned object and `LinearizedSS`.** `linearize` returns labeled
data — `(ẋ₀, x₀, u₀, y₀, A, B, C, D)` with the taps' label sets — on
which `subsystem`/`delete_vars` survive as pure label-indexed matrix
slicing (no model involvement); the `c172x_ctl` LQR pipeline consumes it
with cosmetic changes. `LinearizedSS` the *component* survives separately
as an ordinary continuous component in the migrated library (`init_x` =
the state vector, labeled faces, the affine update in `h_xu`/`f`) — no
privileges, schema like everyone else.

**Recorded guidance.** Linearization taps should select
minimal-coordinate mechanizations — perturbing Euler-angle states is
meaningful where seeding quaternion components steps off the unit-norm
manifold. This is why today's code linearizes the `{NED}` variant;
`design_world(ac)` rigs it, promoting implicit practice to stated rule.
The coordinate choice belongs to the tap author, not the framework.

**Recorded, not built: the sampled-data Dual activation.** The
frozen-exact doctrine is consumer-scoped, not a capability wall: today's
services differentiate the continuous dynamics with the discrete state held, for which a
frozen discrete output — a ZOH constant with zero partials — is the exact
answer, enforced by the type system ([§11.2](#112-the-declaration-inventory)). Stated once, because the
question recurs: **the frozen discrete cell is not an AD limitation on the
signal path; it is the true zero of an instantaneous dependence that the hybrid
semantics never had** — the dataflow through a discrete component is temporal,
not instantaneous, and AD follows actual dataflow (`frozen_discrete_walkthrough.md`
works the three-component chain through). Differentiating "through" the
discrete side means differentiating a *different object*: the sampled-data
step map $\Phi : (x_k, \mathrm{slots}) \to x_{k+1}$ over the model's *whole*
state, continuous and discrete leaves alike (integrate one period,
then run the due ticks). The extension is additive along existing seams:
walked-leaf parametrization of the discrete tier's real-scalar state leaves (counters/enums stay
pinned, like `m`); opt-in participation on discrete components (frozen-exact
stays the default; a participating component opts in through an explicit
trait, which **brings the two-argument `T`-form of `output_types` with it** —
the trait flips the leaf's mandated declaration shape from the plain form to
the continuous one ([§11.2](#112-the-declaration-inventory), [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)), so participation stays authored per leaf on
that tier too, and the hinge is recorded here so the two forms stay compatible
— graceful migration, no flag day); one new [§12.4](#124-activations-executable-sets-laziness-caching) activation ("continuous chain + `f` +
discrete `h_x`/`h_xu` + `g`"); and forward sensitivities through the in-house
RK steppers for free, a payoff of owning the loop ([§8.1](#81-loop-ownership-the-framework-owns-the-simulation-loop)). The honest
boundary: $\Phi$ is differentiable only where the event pattern is locally
constant — exactness across a firing needs saltation corrections — so the
scope is event-quiescent operating points (which [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)'s guards-at-commit
already makes trim points) plus a loud diagnostic if an event fires inside
a differentiated step. Consumers waiting: [§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)'s closed-loop trim door
($g(x) - x = 0$ residuals currently imply the derivative-free fallback,
since frozen `g` has no Jacobian columns) and exact discrete-time
linearization of the full loop (digital design on the exact discretized
plant instead of continuous linearization + Tustin).

**Delivered, and what is still only recorded: declarative non-participation.**
**Both halves of this door are now built.** The output half (row 166): a
continuous producer's declaration is per-leaf, so "this port is frozen under
differentiation" has a spelling — declare the leaf `Float64` and strip with
`ForwardDiff.value` inside the stage ([§11.2](#112-the-declaration-inventory), [§12.5](#125-the-always-on-conformance-check)). An opaque wrapper (an FMU,
a C aerodynamic table) and a deliberately severed coupling can therefore both
say so in the schema instead of showing up in Jacobians as unexplained zero
rows, and the conformance check holds them to it at every activation. The input
half (row 167): a consumer's entries are per-leaf too, so a `Float64` entry
declares "never hand me partials" — the AD-incompatible component's own
statement, enforced at the wire ([§6.1](#61-connections-and-hierarchy)) — and, at a root slot, *is* the
forbid-seeding marker this section once recorded as speculative. The marker
therefore exists, and it carries semantics rather than mere protection: it types
the slot cell at every activation — so an unseeded slot is a *choice* where a
`Float64`-entry slot is a *declaration* — and tap resolution rejects the latter
with the offending entry in hand instead of returning a silent zero column. What
stays recorded is only the remaining **tooling** over that visibility:
pinned-face validation by the tap declaration (selecting a declared-frozen
output = warning) and a feedthrough-graph lint (a frozen output fed by
participating inputs names the severed coupling). Additive when a consumer shows
up, no flag day; until then the declared pins and the visible zero rows suffice.

---

# Part V — Grounding

## 15. Case studies

### 15.1 `Vehicle` today → this framework

The grounding exercise that validated [§5](#5-evaluation-order-and-feedthrough). Current `Vehicle.f_ode!`
(`aircraftbase.jl:142-170`) is a hand-woven instance of the machinery specified here:

| Today (convention) | This design (checked structure) |
|---|---|
| `kinematics.u .= dynamics.x` — velocity extracted directly from the state vector because `f_ode!(dynamics)` can't run yet | `dyn`'s stage-1 output, scheduled first by construction; the artificial loop in `VehicleDynamics` (velocity state-only, accelerations feedthrough) dissolves |
| Hand-ordered `f_ode!` body (kinematics → airdata → systems → route five `dynamics.u` assignments → dynamics last) | Build-time topological sort; wrong wiring = build error naming the cycle or dangling port |
| Velocity state duplicated (`dynamics.x` and `kinematics.u`) with manual sync, incl. `dynamics.x .= kinematics.u  #essential` in `f_init!` | One state, one owner; consumers wire to `dyn.vel` |
| `get_wr_b`/`get_mp_b`/`get_hr_b` generated tree-walk sums | Summing junctions at ownership boundaries, one explicit wire per contributor, exported totals ([§6.2](#62-aggregation-explicit-summing-junctions)) |
| `f_step!` quaternion renorm + engine-phase/stall-latch checks | `project` hook + boundary-detected events with defined semantics |
| `Aircraft.f_ode!` runs avionics before the vehicle → continuous avionics reads one-stage-stale `vehicle.y` (implicit delay) | Avionics scheduled inside the sweep after the stage-1 outputs it consumes — no delay; or declared periodic and samples post-step by stated semantics |
| `atmosphere`/`terrain` threaded as arguments through every signature | Field-handle signals through ordinary ports ([§4.4](#44-function-valued-signals-environment-access)) |

The migration cost surfaced by the same exercise: today's monolithic `KinData` splits
into `pose` (stage 1: `q_eb`, `r_eb_e`, `ϕ_λ_h`, ...) and `kin_vel` (stage 2:
`v_eb_n`, `v_gnd`, `χ`, `γ`, ...) because the parts genuinely have different
dependencies. The recurring trade, stated once: the framework asks authors to write
down structure previously kept in their heads, and pays them back by never letting it
silently rot.

The genuine algebraic loop in the domain — α̇-dependent aerodynamics — is already
broken in the current C172 model by a filter state, exactly the explicit break [§5.5](#55-algebraic-loop-policy-reject-at-build-time)
prescribes. Evidence that reject-loops matches domain practice rather than fighting it.

### 15.2 Torture tests for the §5.2 interfaces: `PistonEngine` and the FCS PID cascade

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
- `f_step!`'s transitions become boundary-detected events with mixed predicate/threshold guards
  ([§2.1](#21-events-two-detection-policies)); `fuel_available` becomes an ordinary port (state-derived at the fuel system,
  hence stage-1 — no loop).
- Forced publications: none — everything `f` reads was already in `PistonEngineY`.

**`PID`** (control.jl:431-471) and the C172X FCS — the discrete side's representative:

- The current update entangles outputs and next state by construction (`y_i = x_i`:
  this tick's integral-path output *is* the updated integrator state). Under [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries) the
  law runs once in `h_xu`, publishing paths, saturation and the updated states;
  `g` is a three-field copy. Under the orthodox split, `g(x, u, t)` would reproduce
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
  stage-1 port for the previous saturation (`sat_out_0` — an `x` field declared in the
  LQR's output contract, hence auto-published at stage-1 position, [§11.3](#113-visibility-the-contract-is-the-interface)). The delay
  becomes an explicit property of the wiring. (The loop and its fix are
  formalism-independent; the framework's contribution is refusing to let the
  ambiguity through, and stage 1's is having the delayed value already on a port.)

Both components passed without blockers, with zero publications forced beyond current
practice — the empirical basis for [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)'s claim that derivative/output overlap is the
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
with zero special-casing ([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)). One-shot design-time gains (`robot2d`'s
controller synthesis at init) are construction-time parameters or stopped-sim
service outputs — not a runtime write path.

*Resets are same-tick inputs, consumed in the output stage.* The supervisor
publishes `engage` and the latch value from its own feedthrough stage; the
compensator, topologically after it, honors them **this tick**:

```julia
h_xu(c::PI, (; x, u)) = (; u_cmd = u.engage ? u.u_latch : c.k_p*u.e + x.x_i)
g(c::PI, (; x, u, Δt)) = (; x_i = u.engage ? u.u_latch - c.k_p*u.e
                                           : x.x_i + c.k_i*Δt*u.e)
```

Honoring the reset only in `g` is legal and means something else: the state
still lands correctly at the next tick, but the *output at the engagement
tick* was already published from the stale state during the sweep — and under
ZOH the plant integrates a full step under it. That one-tick-late command is
exactly the bump that bumpless transfer exists to remove, and no diagnostic
can catch it (both spellings are meaningful designs). The update stage cannot
rescue its own boundary — republish-from-`x⁺` is rejected (row 67) — so the
output stage is the *only* same-tick path. Today's hand-ordering (`f_init!`
before `f_periodic!` in one call) is this contract enforced manually;
[Appendix A](#appendix-a-taught-contracts-the-author-facing-index) carries it as the same-tick reset entry, and [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)'s
bumpless-engage answer — engage semantics live in the FCS — presupposes
exactly this spelling. One relative outside the FCS: the landing gear's
level-triggered cross-component reset (`!wow` re-initializing the friction
regulator every step) becomes an edge-triggered event owned by the regulator
— a semantic tightening recorded in [§16](#16-open-axes)'s migration mapping. There the
respelling is not a stylistic one: the continuous tier admits no input
spelling at all, because only handlers write `x` ([§3.1](#31-continuous-component-the-hybrid-primitive)), so the event is
necessity rather than taste ([Appendix A](#appendix-a-taught-contracts-the-author-facing-index) carries the contract), and the
reimplemented `PIVector`'s optional reset face ([§16](#16-open-axes)) is sugar over
exactly this event.

### 15.3 Torture test for the §9 staging shapes: filter, joystick and GUI

The exercise that selected per-device cells ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain)) and produced the [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)
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
  cell across the pause (the knob keeps it, [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract) peek rule), merges with the
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
  panel reuse ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)) — prompted by asking how the filter's panel survives the
  filter becoming an embedded component with `u_cmd` driven by another component
  (the `Cessna172Xv0` → `Xv1` throttle situation).
- **Note**: under slot exclusivity ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)) the contested-`u_cmd`
  scenario itself becomes an attach-time error — the drag-against-the-stream
  phase cannot arise. The test's verdicts on cell *shapes* (atomicity,
  coalescing, pause behavior, the peek rule) stand; its conflict-precedence
  findings and the active-widget stage-every-pass contract ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)) do not —
  both were patches for the contested-slot world this test examined, and that
  world is gone.

### 15.4 The interactive C172X demo: the periphery under load

The full-fidelity successor to [§15.3](#153-torture-test-for-the-9-staging-shapes-filter-joystick-and-gui), against the real deployment:
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
  exclusivity ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)): claim/disable covers every case found; none needs two
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
  conditioning upstream as mapping data, semantics in-model ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain)).
- **Mode engage** (`mode_req` + setpoint capture from current measurements — the
  GUI handler does `u.EAS_ref = EAS` read from `vehicle.y`): the one place the
  GUI composes writes from model state. Open fork: GUI peek-batch ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract) supports
  it as-is) vs. capture-on-engage latched inside the control laws (uniform across
  all writers, but moves behavior into the FCS).
- **Vehicle-direct and environment tunables** (engine start/stop/mixture, payload
  masses, terrain surface enum, sea-level T/p, wind NED): ordinary component
  inputs exported to root faces; the GUI writes them under its greedy claim via
  [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract); no machinery. The interactive surface is *not* one thing: pilot commands cluster
  under a prefix; environment knobs stay with their components' panels.
- **The Xv1 actuator sliders**: FlightCore's dead sliders; resolved read-only by
  [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract). No action.
- **Outbound** (XPlane12: control-surface angles, nose-wheel steering, prop
  speed/phase, pose, `t`): a snapshot-consuming device, pure `map_output` on the
  device task ([§9.2](#92-outbound-snapshot-publication)). No friction found.
- **Init/trim, pause/pace, post-run plots**: stopped-sim services ([§14](#14-stopped-sim-services)), control
  plane ([§10.1](#101-control-plane)), log/trace ([§9.2](#92-outbound-snapshot-publication), [§9.5](#95-inbound-the-input-trace)).

#### Architectures examined here and rejected

**Architectures examined here and rejected** (the [§9](#9-runtime-periphery-the-data-plane), [§10](#10-runtime-periphery-lifecycle-and-orchestration) periphery decisions were
forced by this cast):

- **Devices as components** (a `T16000M` component wrapping SDL): replay stops
  being same-build (the trace doctrine's strongest property — same type, same
  schedule, staging fed from the recording); the GUI is irreducibly a staging
  device, so inbound uniformity is unreachable anyway (two mechanisms instead of
  one); device lifecycle would duplicate [§10.4](#104-shutdown-protocol) in component vocabulary; and the
  drain stops being the single audit point for external data. The salvage: the
  *knowledge* half (a device model's semantics) is expressible as an ordinary
  in-model component wherever wanted; only the wall-clock pump stays outside. In
  an interactive, paced world the scheme is internally consistent (frame-top
  hardware sampling is well-defined) — the rejection rests on the invariants
  above, not on the [§10.5](#105-scripts-and-the-mid-run-mutation-doctrine) clock criterion.
- **A root-level `PilotInterface` cockpit component** assembling `pilot_commands`
  beside the physical models: its claimed jobs dissolved one by one — struct
  assembly happens in-model downstream of scalar faces (any component can gather
  and bundle), curves became mapping data, widget arbitration is [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract) +
  exclusivity, and the stateful residue (accumulators, capture-on-engage) fits
  the avionics, where FBW stick shaping arguably belongs. What remained was a
  component with no natural place — a cockpit artifact sitting beside
  aircraft/terrain/atmosphere in `World` misstates the composition.
- **Bundled command faces** (`pilot_inputs` as one struct port): kills per-field
  claiming, liveness and trace provenance — the port is the periphery's atomic
  unit ([§4.3](#43-table-mechanics-and-port-granularity) write-side corollary). The routing convenience the bundle bought in
  FlightCore's argument-threading world is provided here by the namespace prefix
  and `input_passthrough` ([§11.8](#118-computed-connections-and-generic-boundaries)); the struct reappears legitimately downstream, assembled
  in-model by a single producer.

#### Surface walkthrough

**Surface walkthrough.** The demo line by line:

- `SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_LOWS15))` —
  pure value construction; no `Model` wrapper (its jobs move into the build).
  `HorizontalTerrain`'s elevation is a plain field (parameter), its surface type
  an input port: the parameter/port split FlightCore kept implicit in
  `U()`-vs-field convention is now the declaration itself. The aircraft's
  `input_connections` block carries the `pilot.*` face group in one place, deep routes
  spanning avionics *and* systems — today's mapping writes flaps/brakes directly
  into `act`, bypassing avionics; that bypass becomes a declared route.
- `Simulation(world; algorithm = RK4(), h = 0.02, n = 1, t_end = 1000)` — `n`
  binds `Δt_base = n·h` ([§8.5](#85-multi-rate-tick-scheduling); default 1: base tick every step). The entire
  build pipeline runs here: class resolution, path validation, face derivation
  (computed boundary connections expanded, printable), two-producers/unconnected checks,
  topological sort, probe passes, rate compilation, flat layout, slot table.
- `init!(sim, ready_for_taxi(ac); t0 = 0.0)` — stopped-sim services ([§14](#14-stopped-sim-services);
  trim is its own service, `trim!(sim, problem; baseline, ...)`, whose commit
  runs the same boundary): they write `(x, m)`, **establish every root
  slot's initial value**, and capture the trace header. Slot initialization
  decisively belongs here, not in declarations: the trim service writes slot
  values it *solved for* (throttle, elevator) — not declaration constants.
- `attach!(sim, XPlane12Control(...), binding)` — output device: claims nothing,
  consumes snapshots via [§10.3](#103-the-next-snapshot-wait), pure `map_output` on its task. Its binding names
  snapshot paths, **validated at attach against the actual contract** — an
  aircraft substitution that breaks the binding fails at attach, not with silent
  garbage UDP (a new, cheap [§9.2](#92-outbound-snapshot-publication) obligation).
- `attach!(sim, joystick, T16000MBinding())` — the binding is a declarative
  table: axis/button → face name + conditioning params
  (`stick_y = (face = "aircraft.pilot.elevator_axis", expo = 1.0, deadzone =
  0.05)`, `button_3 = (face = "aircraft.pilot.flaps_up_count", as = :count)`).
  At attach: faces resolved against the root contract (typo → did-you-mean),
  claim set registered (second joystick on the same faces errors here). The
  Gladiator variant is the same table with different keys, zero shaping code —
  the duplication smell structurally gone.
- `run!(sim; gui = true, pace = 1)` — a greedy claim over every unclaimed face
  and liveness with zero configuration, both settled at run start against the
  frozen roster ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster); the flag's attachment lasts exactly this run,
  [§10.6](#106-run-lifecycle-and-partial-advance)):
  axis mirrors read-only (claimed, with provenance), mode/setpoint/mixture/
  payload/environment widgets live, actuator sliders read-only (component-fed).
  Unplug the joystick → its task exits, the mirrors stay read-only with the
  death in their provenance ("claimed by `T16000M` — task dead") and the axes
  hold their last-drained values — the accepted orphan anomaly ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster));
  recovery is between runs: stop, `detach!`, then `init!` for a fresh
  trajectory or `replay!`-to-end + `run!` to continue the interrupted one
  ([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)). Post-run: `TimeSeries` over retained snapshots; the trace can re-drive
  a fresh `Simulation(world)` bit-identically — which is also the state-trajectory
  inspector (row 38 paying its way).

#### Frame anatomies

**Frame anatomies.** One frame each:

- *Stick motion*: device task polls, conditioning helper applies binding params,
  complete batch overwrites the cell (inter-frame polls coalesce, ZOH-correct);
  drain applies + traces; avionics tick reads the slot fresh; worst-case
  stick-to-physics latency = poll interval + frame, now by stated semantics.
- *Flaps click*: button peeks counter `k` (own-pending-else-snapshot), stages
  level `k+1` on activation; drain applies; avionics compares slot counter to its
  `x` counter, moves the detent, stores. Multi-click in one window counts via
  own-pending-first peek; repeated staging idempotent ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)).
- *Mode engage*: the GUI stages `mode_req` (plus optionally peek-captured
  setpoint slots); **bumpless-engage semantics live in the FCS already** — the
  current `ControlLaws` latches each controller's reference from the present
  command vector on mode transitions, so the fork dissolved: semantic capture is
  aircraft design (status quo, uniform across writers — a script engages sanely
  staging one value); the GUI peek-batch survives as display/slot-sync sugar
  only. Residual check for migration: order-sensitivity of latch vs. sync-write
  on the same boundary (believed none — both derive from the same measurements).
- *Wind slider*: sparse CAS-merge, [§15.3](#153-torture-test-for-the-9-staging-shapes-filter-joystick-and-gui)'s uncontested-`τ` case, live in the
  real cast.
- *Pause/un-pause*: control plane; GUI edits hold in its cell (peek displays),
  joystick cell coalesces bounded; un-pause drain applies both (disjoint slots —
  exclusivity makes the contested question unaskable), pacer re-anchors.
- *Window close*: [§10.4](#104-shutdown-protocol) verbatim — complete boundary, final snapshot, sticky
  stopped, wake waits, unblock hooks, named-timeout joins.

Remaining open (feeding [§16](#16-open-axes)): the `q_sf` home (thin mapping entry vs.
avionics-internal derivation — aircraft design, not framework design).

### 15.5 The strapdown IMU: integrate-and-dump across the tier boundary

The strongest challenge mounted against the [§3](#3-component-taxonomy) class split, and its resolution.
The general question first: why two leaf classes at all — why not one all-in-one
primitive carrying continuous state, modes *and* discrete state, with `f`, events *and* `g`, purely
continuous or discrete components falling out by whichever facets an author
declares? (Class is already read off declaration shape, [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape) — the question is
whether the two declaration sets should be exclusive.)

#### Why the merge buys nothing

**Why the merge buys nothing.** The split is between *time bases*, not state
classes — the continuous primitive is already hybrid (`m`, guards, handlers, [§3.1](#31-continuous-component-the-hybrid-primitive));
what separates the classes is sweep-driven versus tick-driven execution. And the
settled rules force a merged component's two halves to communicate exactly as two
siblings do: one home per datum ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)), `f` sees only the continuous
state and `g` only the discrete one,
and `x⁺` is decoded only at the owner's next tick (`g` runs last) — the very
fact that makes ticks→events structurally impossible and terminates the boundary iteration ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)). Cross-tier
influence inside the merged class still routes through published table cells, so
the all-in-one component is an assembly of two primitives in a trench coat. Its
costs, meanwhile, are real: a stage cannot run both every sweep *and* only at
ticks, so the merged class needs four tier-disambiguated stage functions; tier
stops being a component property readable off one `output_types` signature ([§11.2](#112-the-declaration-inventory)) and
becomes per-port vocabulary; every class-implied obligation (rate required,
`K`-on-continuous error, `Δt` bundle availability, `Dual` activation membership) becomes a
facet-conditional web; and the sampling seam — ZOH and the `z⁻¹` delay, the most
bug-prone boundary in a flight-control stack — disappears into a monolith where
the split keeps it a visible wire. (Simulink and FMI allow the fused block;
sample-time propagation confusion is the documented price.)

#### The counterexample

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
its own `x`; handlers are the sole resetters of continuous state, and they are
guard-driven).
Integrate-and-dump falls squarely into the crack between the classes:
tightly-coupled continuous and periodic dynamics in one physical device.

#### The idiom: integrate-and-difference

**The idiom: integrate-and-difference.** The reset is eliminable by algebra, not
approximation. Every interval-relative integral becomes a *cumulative* one; the
sampler differences against the previous sample, held in its `x` — the textbook
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
  `υ_c_sc = x.q'(u.V - x.V)`. The factor leaving the integral is the **anchor
  change between two inertially-fixed frames** — constant because $t_{k-1}$ is
  in the past and latched. The physical intra-interval rotation, the thing
  sculling corrections are *about*, stays inside the integrand via $q(t)$: every
  RHS evaluation, RK stages included, applies the current cumulative attitude,
  exactly as the sketch applies the current `q_c_cc`.

#### Exactness condition, stated once

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
input_types(::IMUIntegrals, ::Type{T}) where {T <: Real} =
    (q_eb = RQuat{T}, r_eb_e = SVector{3,T},
     ω_eb_b = SVector{3,T}, a_ib_b = SVector{3,T}, α_ib_b = SVector{3,T})
output_types(::IMUIntegrals, ::Type{T}) where {T <: Real} =
    (Θ = SVector{3,T}, q = SVector{4,T},                        # auto-published state
     Υ = SVector{3,T}, V = SVector{3,T},
     ω_ic_c = SVector{3,T}, f_c_c = SVector{3,T})               # instantaneous truth

# h_xu: the sketch's f_ode! math verbatim (lever arm, gravity, Earth rate) → (; ω_ic_c, f_c_c)
function f(imu::IMUIntegrals, (; x, y))
    q = RQuat(x.q, normalization = false)              # §7.1's explicit cast
    (Θ = y.ω_ic_c, q = SVector{4}(Attitude.dt(q, y.ω_ic_c)), Υ = y.f_c_c, V = q(y.f_c_c))
end
project(imu::IMUIntegrals, x) = (; x..., q = normalize(x.q))   # SVector normalize

struct IMUSampler <: AbstractComponent end
init_x(::IMUSampler) = (Θ = zeros(SVector{3}), q = SVector{4}(1.0, 0, 0, 0),
                        Υ = zeros(SVector{3}), V = zeros(SVector{3}))
input_types(::IMUSampler)  = (Θ = SVector{3,Float64}, q = SVector{4,Float64},  # discrete class: plain
                         Υ = SVector{3,Float64}, V = SVector{3,Float64})       # form, bound check only
output_types(::IMUSampler) = (sample = IMUSample,)   # discrete class: cells pin (frozen-exact)

function h_xu(s::IMUSampler, (; x, u, Δt))
    q_x = RQuat(x.q, normalization = false);  q_u = RQuat(u.q, normalization = false)
    ϑ_c = u.Θ - x.Θ;  υ_c = u.Υ - x.Υ
    Δq  = q_x' ∘ q_u                                   # interval rotation, exact
    υ_c_sc = q_x'(u.V - x.V)                           # constant anchor change pulled out
    (; sample = IMUSample(; ω̄_ic_c = ϑ_c / Δt, f̄_c_c = υ_c / Δt,
                            ϑ_c, ϑ_c_cc = RVec(Δq)[:], υ_c, υ_c_sc))
end
g(s::IMUSampler, (; u)) = (Θ = u.Θ, q = u.q, Υ = u.Υ, V = u.V)   # the latch
```

The `IMU` assembly wires the four integral ports across, holds the error model as
a discrete sibling consuming `sample`, and leaves the sampler at `K = 1` in its
own scope — the parent sets the IMU's rate ([§11.7](#117-rate-scopes)). `Δt` in the stage bundle is
the [§8.5](#85-multi-rate-tick-scheduling) single source of truth, put there for exactly this kind of discretized
law. (Initialization
consistency — the sampler's `x` must equal the initial integrals or the `t₀` sample is
wrong — holds by default at zeros/identity, and boundary zero discharges the
rest: its due `g` latches `x ← integrals(t₀)` for every subsequent sample, so
only the `t₀` sample itself depends on the authored `x` — a condition-authoring
obligation under trim, [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions).)

#### Why `u.V` is fresh — the line that would silently zero

**Why `u.V` is fresh — the line that would silently zero.** The sculling line is
correct only because a due tick samples the *completed* boundary: if `u.V` still
held the previous boundary's decode it would equal `x.V` exactly (that is the
value `g` latched), and sculling would vanish without an error anywhere. The
guarantee is the [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) macro-sequence, not a scheduling accident: integrate →
project → sweep, with the due sampler's stages gated *into* that sweep ([§8.5](#85-multi-rate-tick-scheduling)) and
the integrals arriving at stage-1 position (auto-published state, [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)) — before
any stage-2 function runs, regardless of topological placement. The rest of the
timeline closes consistently: the sampler's `h_xu` decodes `x` (the `t_{k-1}` latch)
*before* `g` runs — the `z⁻¹` semantics — and after event quiescence `g` latches
the `t_k` values for the next tick; same-boundary events re-run the gated stages
in their re-sweeps, so `g` and external readers see the settled boundary.

#### Author-knowledge note

**Author-knowledge note** (user observation, recorded as a documentation
obligation): the clean implementation leans on the author *knowing* that "sampling
at `t_k`" means post-integration, post-projection, stage-1-fresh state. That
knowledge must be part of the framework's taught contract — [§8.5](#85-multi-rate-tick-scheduling)/[§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) semantics
stated in component-author documentation with this IMU as the worked example — not
internal lore. The failure mode of not knowing it is instructive: an author who
distrusts the sweep order adds a defensive one-tick delay or re-derives the
integrals in the sampler, silently degrading the model.

**When the coupling is genuinely two-way: the latch-back wire.** The IMU's
coupling is one-directional (integrals → sampler). If the flow itself needed the
interval-relative value — integrator saturation within the sampling interval,
say — the latch becomes a wire back: the sampler publishes the sample-instant
values from its *feedthrough* stage (`h_xu` reads `u`, so the latch port carries
the current tick's values, ZOH until the next; an `h_x`-published latch would be
one period stale), and the continuous `f` computes `x − u.latch`. Both cross-wires
consume the other side's ports and the schedule stays acyclic (integrals stage 1 →
sampler `h_xu`; sampler `h_xu` → the integrals' `f`-edge, [§5.4](#54-artificial-loops-and-the-escape-hatch)). The "reset"
becomes a visible tier-crossing feedback loop — which is what it always was,
physically.

**Verdict.** The strongest counterexample landed on the two-class taxonomy with
*less* code than the fused original — same thirteen integral scalars, same math,
minus the reset block — and three structural gains. The sampling seam became a
wire. The sketch's incidental violations became visible structure: the
`CircularBuffer` mutated inside the component struct (constants) moves to the
consumer's `x` or falls out of the log; the parent-called `f_disc!(errors)`
becomes a discrete sibling, making the truth/corrupted sample pair separately
loggable. And linearization got sane: under a `Dual` activation the discrete tier
is held ([§11.2](#112-the-declaration-inventory)), and "integrators that never reset" *is* the cumulative
formulation — the framework's rules pushed the model into the only form its own
linearization semantics could coherently handle. Residual escape hatch, recorded
unbuilt: if interval-relative dynamics ever neither factor algebraically nor
tolerate the latch-back wire, the guarded addition is a **tick-triggered handler**
on continuous components (periodic events). Nothing surveyed needs it, and it
would be the camel's nose for the merged class.

---

## 16. Open axes

Still to be settled:

- **Migration.** Outline for FlightPhysics/FlightApps (the walked-leaf parametrization
  pass — whose `Ranged` rewrite targets [§11.2](#112-the-declaration-inventory)'s walk rule where `Ranged`
  survives, at ports and parameters: constructor
  discipline admitting the walked scalar with the value parameters left alone,
  plus a `probe_value` method — the `KinData`-style output splits, the
  contributor survey feeding [§6.2](#62-aggregation-explicit-summing-junctions)'s
  aggregation chains — mechanical to extract from today's trait implementations);
  comparison criteria against FlightCore's demonstrated strengths (zero-alloc
  stepping — measured through [§12.7](#127-the-compiled-executor)'s `phase_bodies` seam, apples-to-apples
  with today's `@ballocated f_ode!` suites — flexibility, interactive
  operation); the [§13.7](#137-tooling-consequences-provenance-and-the-component-library) component library's
  starting inventory; the **conventional exported aircraft surface** for
  generic periphery consumers ([§9.2](#92-outbound-snapshot-publication)'s integration register): pose and
  velocity faces with wrapper types — `VelocityData`, field meaning defined at
  the type — as the `KinData` successor's periphery-facing half; the
  **supervisor seam** ([§15.2](#152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade)): compensator gain ports plus scheduler
  components (~7 for the C172X), the same-tick reset respelling of every
  mode-transition latch, and the gear's level-triggered reset converted to
  an edge event — which lands on the *library* side: the reimplemented
  `PIVector` gains a **flag-gated reset face**, `PIVector(; reset = true)`
  adding a `Bool` input face plus the event, the default omitting both
  (declarations are ordinary functions of the instance, [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape) — the honest
  version of Simulink's checkbox), under one fixed policy: rising edge →
  reset to the declared `init_x` values, implemented internally as an
  ordinary guard/handler event ([Appendix A](#appendix-a-taught-contracts-the-author-facing-index)'s continuous-reset contract in its
  worked instance). Falling-edge consumers wire a NOT gate ([§13.7](#137-tooling-consequences-provenance-and-the-component-library)'s Bool
  gates); level-pinning and reset-to-an-external-value are different blocks
  (tracking), not options. The gear then wires `strut.wow → frc.reset`: the
  **touchdown** edge ([§2.1](#21-events-two-detection-policies)'s not-holding → holding semantics), fresh
  regulator state per contact episode. The liftoff edge (`!wow`) was
  rejected — its equivalence to today's level reset rests on the strut's
  airborne zero-default (`v_ec_xy = [0,0]` in the no-contact branch) plus the
  integrator leak, a cross-component dependency the touchdown edge dissolves
  rather than documents. Boundary-detected policy suffices (the regulator's
  input ramps from zero at touchdown; localization buys nothing), and a sim
  initialized on ground fires the reset at boundary zero harmlessly (declared
  inits are zero, and boundary-zero priors are not-holding, [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)). The
  engine's two `PIVector` instances (`PistonEngine`'s `idle` and `frc`)
  migrate **unchanged, flag off** — verified reset-free in today's code,
  where windup across unused phases is already handled by the saturation
  bounds and `int_halted`; their `f_init!` gain writes become
  construction-time parameters (row 89, as `Contact`'s do). The PI *law* is
  shared as plain pure functions called by the block's stages (row 139's
  laws-as-plain-functions pattern), and `sat_ext` poses the same
  always-on-vs-flag-gated face question, to be decided at reimplementation
  time on the same axis; the **steering contract re-factoring** ([§5.4](#54-artificial-loops-and-the-escape-hatch)'s middle rung,
  worked on the shipped instance): `AbstractSteering` moves from "give me the
  angle" to `(engaged, ψ_cmd)`, with the castoring fallback
  (`ψ_sw = engaged ? ψ_cmd : ψ_v`) computed inside `Strut`, which deletes the
  strut → steering → strut artificial loop that stage-2 conservatism would
  otherwise manufacture — beside [§15.1](#151-vehicle-today--this-framework)'s `VehicleDynamics` instance, which
  dissolves under the two-stage split alone; splitting `Strut`, its shared
  geometry crossing the new boundary as one `StrutGeometry` bundle port, is the
  residual remedy, recorded and not taken (an aircraft-library call — a
  component's own contract — recorded here, not framework vocabulary); the
  **state-declaration conversion to [§7.1](#71-continuous-state-structured-immutable-flat-backing)'s closed
  vocabulary** (each `RQuat` state field becomes its `SVector{4}` backing
  with the explicit `normalization = false` cast at its use sites — today's
  `Attitude.dt` already delivers the 4-wide rate — and each `Ranged` state
  field a plain scalar, its clamp respelled as dynamics or projection, never
  construction); the **exported-name surface**, decided deliberately rather than
  by accident: `condition`, `fragment`, `at`, `capture`
  and the `merge` overload ([§14.2](#142-fragment-composition-locality-without-schema)) are generic names sharing a namespace with
  FlightPhysics domain code, and `merge` in particular is a piracy surface
  whose mixed-argument methods must stay error methods — the selector
  family's `get_` prefix ([§14.4](#144-two-application-registers-over-one-plan)) already settles this for the readers, and
  whether the condition algebra ships behind a submodule is the packaging
  question. The audit is a full-surface sweep (per user, 2026-08-01): every
  API method name is either specific enough to export or gets renamed or
  left unexported — with *unexported* the preferred disposition for
  extension-only surface: the declaration and stage family of [§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)'s import
  list — the larger half of the question, on every component file's first
  line, settled there — plus the [§9.6](#96-devices-one-authoring-contract-no-taxonomy) binding interface `claims`/`reads` and
  the side traits `is_input`/`is_output`/`is_greedy` (with `map_input`/`map_output`
  outside the question, being loop-idiom conventions the framework never calls)
  and the device contract `init!`/`loop`/`shutdown!`/
  `unblock!`/`needs_calling_task`, which authors extend by `import` or qualified name, `Base.show`-style,
  rather than call every day. Its criterion is the **four-register naming
  convention** (row 144): declarations the author defines and the framework
  calls are bare nouns or `init_*`/`_types` (`child_connections`,
  `input_connections`/`output_connections`,
  `events`, `input_types`, `workspace`, the stage letters, [§9.6](#96-devices-one-authoring-contract-no-taxonomy)'s `claims(b)`);
  value selectors called against `reads` and snapshots carry `get_`
  ([§14.4](#144-two-application-registers-over-one-plan)); lifecycle and mutating actions are verbs, `!` when they mutate;
  build primitives ([§13.3](#133-build-primitives-resolve-and-the-face-list-accessors)) are plain verbs. A name in the wrong register is a
  rename candidate on that ground alone — which is what the `passthrough`
  rename ([§11.8](#118-computed-connections-and-generic-boundaries), now `input_passthrough`, row 171) already applied, retiring a caller-side helper that wore
  declaration dress and collided with `get_face`; the [§9.6](#96-devices-one-authoring-contract-no-taxonomy) binding-method
  renames (`faces` → `claims`, `selectors` → `reads`, row 146) then applied the
  convention's **semantic axis** — right register, wrong noun: both were correctly
  bare-noun declarations but named their *content*, where `claims` and `reads`
  name the *consequence* the declaration has. (That axis's earlier exemplar was
  `exports`, retired by row 170's split; the `*_connections` family that replaced
  it names content deliberately, for authoring transparency — a recorded choice,
  not register drift.) Five items are flagged
  for the sweep and deliberately not settled now: `input_faces`/`output_faces`
  (noun accessors punning on the `_types` declarations, mitigated by being
  framework-facing); `workspace` (a declaration whose bare noun reads as an
  accessor — every candidate replacement is clunkier, so lean keep); `loop`
  ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)'s device contract — a
  mutating task body spelled as a bare noun among its verb-`!` siblings
  `init!`/`shutdown!`/`unblock!`, with `run!` taken and the "loop body" prose
  entrenched, so it needs the audit's whole-surface view); the bare-noun
  accessor family `trace(sim)`, `latest(sim)`, `binding(handle)`,
  `phase_bodies(sim)` (value selectors outside register (2)'s `get_` rule,
  `trace` the sharpest — the constructor kill-switch `trace = false` and the
  post-run accessor `trace(sim)` are one name in two senses, the overload
  pattern rows 122/144 retire); and whether register (1) needs an explicit
  exemption for predicate traits (`is_greedy`, `needs_calling_task`) —
  boundary cases drawn after row 144's list was fixed,
  not defects. The
  **[§12.7](#127-the-compiled-executor) executor compile-cost re-measurement** runs on
  the real vehicle skeleton — early, before the executor's shape hardens.
  Residuals: the `q_sf` home ([§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load) — aircraft design,
  belongs here); whether `stop_on` needs a root-declared overridable default
  ([§13.5](#135-termination-is-a-state-not-an-exception) — reopen only if the ctor argument proves chronically forgotten).

- **GUI panel authoring API.** The semantics are settled ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract): derived
  liveness, first-class read-only rendering, own-pending-else-snapshot peek,
  stage-on-interaction, orphan display); the calling convention — context
  contents, port naming, child composition — is deferred to migration,
  co-designed against the GUI library under [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)'s four constraints.

- **Log and trace persistence.** The in-memory artifacts
  are settled — the log as retained boundary snapshots ([§9.2](#92-outbound-snapshot-publication)), the
  always-on device-tagged input trace with its header of initial stores and
  slot values ([§9.5](#95-inbound-the-input-trace), [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions), [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)), and the primary/derived rule (the log
  is recomputable from the trace, never the reverse) — but nothing on-disk
  is. Deferred to migration, where the consumers exist to ground the
  choices: the HDF5 export scope (whole snapshot log vs. selected
  subtrees), field-handle summarization over retained snapshots (the
  successor to `TimeSeries`'s `getproperty` navigation — today's
  post-processing entry point), and the trace file format, which doubles as
  the reproducibility carrier ([§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)'s replay pointers name positions in
  it).

---

## Appendix A. Taught contracts: the author-facing index

The build pipeline enforces structure — declarations, wiring, types,
conformance. A residue of *semantic* facts is unenforceable by any check:
knowing them is what makes component and periphery code come out right, and
not knowing them produces defensive delays, duplicated math or mistimed
samples with no diagnostic firing anywhere ([§15.5](#155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary)'s author-knowledge note is
the archetype).
This appendix is an **index, not a second home** — one recall line per
contract, with the normative statement staying in the owning section (one
home per datum, applied to the document itself).

For component authors:

- **The stage funnel** ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)). Stage name ⊇ bundle ⊇ destructured reads:
  the stage name fixes the maximal legal view set, the component's
  declarations narrow it to the bundle, and the signature's destructuring
  narrows it to actual reads — so "no `u` in the name" *is* the
  no-feedthrough property. The teaching line: stage 1 publishes what you
  know from state alone; stage 2 adds what needs inputs; your dynamics
  read your own published results instead of recomputing them.
- **One home per datum** ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), [§4.3](#43-table-mechanics-and-port-granularity)). The signal table holds *produced*
  signals only, never transported ones: buffer for continuous `x`, stores for
  discrete `x` and for `m`, table for signals — no store mirrors another.
- **The value-level constructor** ([§4.4](#44-function-valued-signals-environment-access)). A field-emitting component ships
  the map (component, input values) → handle as a plain exported function,
  and its output stage merely calls it: [§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)'s condition math must be able
  to produce the sweep's exact handle outside any sweep, and only the
  component's author can write that function without re-creating the
  drift class.
- **Boundary sampling** ([§8.5](#85-multi-rate-tick-scheduling)/[§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event); worked example [§15.5](#155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary)). "Sampling at
  `t_k`" means post-integration, post-projection, stage-1-fresh state: a
  due tick's gated stages run inside the boundary sweep and sample the
  *completed* boundary. Distrusting this — a defensive one-tick delay, a
  re-derivation inside the sampler — silently degrades the model.
- **Interval alignment** ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)). A boundary's `g` is the *outgoing*
  transition: at tick `t_k` it consumes the completed boundary's samples
  and produces `x_{k+1}` — the value the component's *next* tick decodes
  (the sampled-data `z⁻¹` delay, by construction). Hence `g` runs at
  boundary zero: it is the `t₀` sample's only chance.
- **Same-tick reset consumption** ([§15.2](#152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade)) — *discrete tier*. A commanded reset
  of a discrete component's `x` is an input. For same-tick output semantics
  the *output stage* consumes it — overriding
  the state-derived path — and `g` stores the matching `x⁺`; a reset honored
  only in `g` reaches the outputs one tick late (the plant integrates a full
  step under the stale command). Both spellings are legal; they mean
  different things. The continuous tier has no such choice — next entry.
- **A continuous component's state reset is an event** ([§3.1](#31-continuous-component-the-hybrid-primitive), [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event), [§15.2](#152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade)).
  Only handlers write `x`, so even a *commanded* reset — the condition
  arriving as an ordinary `Bool` input, weight-on-wheels being the shipped
  instance — is spelled as an event whose guard reads that input; the
  discrete tier's input spelling does not transfer. The reason is semantic,
  not stylistic: only the discrete tier's update stage is already a jump map,
  so a reset there is just another value for `x⁺`, whereas a continuous state
  jump must be solver-visible, applied *between* integration segments — the
  flow/jump split every hybrid tool converges on (Simulink applies its reset
  ports through zero-crossing events plus a solver restart; Modelica's
  `reinit` is syntactically legal only inside a `when`). And there is no
  stale-output hazard to manage: [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) re-sweeps outputs to quiescence after
  handlers, so a continuous edge-reset is same-boundary by construction.
- **Guard predicates, edges and priors** ([§2.1](#21-events-two-detection-policies), [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)). A guard defines
  a predicate — a `Bool` form, or a sign value `σ` with
  positive = holding. Events fire on not-holding → holding *edges* against per-event
  priors (the previous boundary's quiescent sample): a predicate that
  keeps holding fires once, at the boundary where it first held. Boundary
  zero sets every prior to not-holding, so a predicate already holding in
  the authored state fires at `t₀`. The opposite crossing direction is a
  second event with the negated guard.
- **Handler-phase visibility** ([§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries), [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)). A handler executes
  against exactly the world its guard fired on: own `y`, foreign `u` and own
  `x`/`m` are all the firing round's sweep, and a component fires at most one
  event per round — later own events are re-decided against the next round's
  sweep, one round per causal link, within and across components alike. The
  table is written only by sweeps.
- **Stage totality** ([§12.3](#123-probing-and-input-synthesis); [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor), [§13.5](#135-termination-is-a-state-not-an-exception)). Stage code is total over
  type-valid inputs: the probe evaluates every user function against values
  chosen for their types alone, and a value-level throw is a build failure
  there and a `StepError` at runtime. Physical plausibility is a published
  `Bool` and `stop_on`; self-consistency asserts belong in tests; parameter
  validation belongs at instance construction, not inside a stage.
- **Stop-face sampling** ([§13.5](#135-termination-is-a-state-not-an-exception)). Stop faces are read in completed-boundary
  snapshots; declare a localized event if the stop needs localizing.

For periphery authors and consumers:

- **Levels, never deltas** ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain)). Staged input values are levels
  (`press_count = 17`, never `presses += 1`) — idempotent under
  coalescing; button edges ride as monotonic counters. Cross-datum state
  (press counters, edge detection) lives in the device struct, maintained
  by the loop, arriving *inside* the datum — `map_input` is pure ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)).
- **The device loop idioms** ([§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§10.4](#104-shutdown-protocol)). Loop on `running(handle)`;
  make every blocking call interruptible (an `unblock!` override, or
  timeouts); voluntary exit is returning. Three canonical shapes:
  timer-poll (sleep, poll, stage), source-driven (block on your socket;
  `unblock!` closes it), boundary-driven (`wait_next_snapshot`, gather,
  send). A forgotten predicate check surfaces as `DeviceJoinTimeout` with
  your device's name; a stall as a stale heartbeat.
- **`shutdown!` closes only what is open** ([§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§10.4](#104-shutdown-protocol)). The framework
  runs it on every exit path, your own `init!`'s failure included: a throw
  half-way through acquisition hands the half-built device straight back to
  you, so guard each release (`isopen`, a `nothing` handle) rather than
  assuming initialization completed. The converse is a burden you do *not*
  carry: `init!` owes no cleanup of its own.
- **Binding traits are declarations, mappings are your own idiom** ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)).
  Keep `is_input`/`is_output`/`is_greedy` trivial — a literal, or a flag read
  off a field fixed at the constructor call — because the framework calls them
  once, at attach, and cross-checks each against the enumeration method it
  implies. `map_input`/`map_output` are the other kind of thing: conventions of
  the loop idiom, called only by your own `loop`, never by the framework — the
  names are worth keeping for readers, and nothing enforces them.
- **Bad datum vs. bug** ([§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)). Catch what your parser can throw,
  `report!(handle, MalformedDatum(cause))`, stage nothing, continue; let
  everything else propagate — the wrapper makes it `DeviceCrash`.
  Tolerating everything hides bugs as "device attached, nothing happens";
  tolerating nothing kills a live link on its first truncated datagram.
- **Derived liveness** ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)). A widget is live iff its port's feed chain
  terminates in a root slot inside the GUI's own claim in the run's frozen
  partition; there is no per-port marking, and unexported ports are
  unpokeable.
- **The two observation registers** ([§9.2](#92-outbound-snapshot-publication), [§13.5](#135-termination-is-a-state-not-an-exception)). A deep snapshot path is
  the *inspection* register: it sees everything and promises nothing
  across builds. An exported output face is the *integration* register:
  curated, writer-independent meaning — the only shield against silent
  semantic drift. Bind faces in anything meant to outlive the current
  build. The store selectors (`get_state`/`get_deriv`) belong to neither:
  they read live stores, never snapshots ([§14.4](#144-two-application-registers-over-one-plan)'s source rule).

---

## Appendix B. API synopsis: the entry points

The user-facing surface on one page — same rule as [Appendix A](#appendix-a-taught-contracts-the-author-facing-index): an index, not a
second home, with each signature normative only where its owning section settles
it. The author-side declaration surface first, then the operator surface by
lifecycle:

**Authoring** — what a component or assembly defines ([§11.2](#112-the-declaration-inventory), [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)–[§11.7](#117-rate-scopes)):

- Continuous leaf: `init_x`/`init_m` (by value), `workspace(::C, ::Type{T})`
  (by allocation), `input_types(::C, ::Type{T})` and
  `output_types(::C, ::Type{T})` (by type),
  `events` — stages `h_x`, `h_xu`, `f`, guard/handler pairs
  (`Event(guard, handler; localize)`), `project`.
- Discrete leaf: `init_x`, `workspace(::C)`,
  `input_types`/`output_types` — stages `h_x`, `h_xu`, `g`.
- Assembly: `child_connections` (mandatory — the class marker),
  `input_connections`, `output_connections`, `rates`.
- Shipped conditions: `condition(::C; kw)` fragment functions ([§14.2](#142-fragment-composition-locality-without-schema)).

Bundle contents by function family (the maximal legal sets, [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws) — signatures
destructure less at will):

| function | bundle fields |
|---|---|
| `h_x` (continuous) | `x, m, t [, ws]` |
| `h_xu` (continuous) | `x, m, u, y_x, w, t [, ws]` |
| `f` | `x, m, y, u, w, t [, ws]` |
| `h_x` (discrete) | `x, t, Δt [, ws]` |
| `h_xu` (discrete) | `x, u, y_x, w, t, Δt [, ws]` |
| `g` | `x, y, u, w, t, Δt [, ws]` |
| guard / handler | `x, m, y, u, w, t [, ws]` |
| `project` | positional `(comp, x)` — no bundle |

Table footnotes, from the bundle law ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)) — the sets above are maximal, and
each field is present only if it exists for the component: `u` iff the function
family may see inputs **and** the component declares `input_types`; `y` iff the
component produces any table cell (`output_types` ∪
auto-published); `x`/`m`/`ws` iff declared; `y_x` iff the stage-1
*return* is non-empty (auto-published names excluded — [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), row 169);
`w` iff the stage handing it down returned one ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)'s one-hop law);
`Δt` on the discrete tier only. Returns: a stage returns a NamedTuple of
port values, or the pair `(y, w)` adding its private intermediates
([§4.3](#43-table-mechanics-and-port-granularity), [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)); `f` returns the layout image of `X` ([§7.1](#71-continuous-state-structured-immutable-flat-backing)); a **handler
returns `(; x, m)` with each key present iff that store exists and the handler
updates it** ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)'s return law — no padding, `x` complete, `m` partial).

**Build.**

- `build(world) → Build` — standalone; the inspectable contract artifact:
  wire list, face table with provenance, schedule, root slots ([§12.2](#122-the-build-artifact)).
  `build(world; activations = (Float64, ProbeDual))` additionally pins
  activation invariants for CI (`ProbeDual` the exported canonical concrete
  probe scalar, [§12.4](#124-activations-executable-sets-laziness-caching)), and pre-materializes activations so a parallel
  sweep shares a fully immutable `Build` ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads), [§12.4](#124-activations-executable-sets-laziness-caching)).
- `resolve(asm, path) → AbstractComponent` — the getfield walk along `/`
  segments, enforcing [§6.1](#61-connections-and-hierarchy)'s generic-boundary rule at the primitive ([§13.3](#133-build-primitives-resolve-and-the-face-list-accessors)).
- `input_faces(c)` / `output_faces(c) → Vector{String}` — declaration-ordered
  face names ([§13.3](#133-build-primitives-resolve-and-the-face-list-accessors)).
- `input_passthrough(asm, path; prefix, sep, except, only)` — the declaration-site
  helper for computed input connections ([§11.8](#118-computed-connections-and-generic-boundaries)).

**Deployment.**

- `Simulation(world; algorithm = RK4(), h, n = 1, t_end = Inf,
  stop_on = (), localization_tol = 1e-6, event_budget = 8,
  trace = true, log = true, log_every = 1, log_max = 65536, debug = false)` —
  wraps the build (`Simulation(world; ...) = Simulation(build(world); ...)`;
  the `Build` overload takes the same deployment keywords and deploys an
  inspected artifact directly, [§12.2](#122-the-build-artifact)); `h` is required (a domain rate is not a
  framework default) and `RK4` is the default stepper ([§8.2](#82-the-stepper-seam)); `n` binds
  `Δt_base = n·h` ([§8.5](#85-multi-rate-tick-scheduling)); `t_end = Inf` is the honest interactive default —
  open-ended in time but bounded in memory, `log_max` being what keeps such a
  session from growing without limit ([§9.2](#92-outbound-snapshot-publication)) — and
  a run with no finite `t_end`, no `stop_on` faces and `pace = Inf` warns at
  start (an unbounded unattended run is almost always an oversight); `stop_on`
  names root-exported `Bool` output faces, OR-combined, recorded in
  run metadata — the trace header's deployment block ([§9.5](#95-inbound-the-input-trace), [§13.5](#135-termination-is-a-state-not-an-exception);
  walkthrough [§15.4](#154-the-interactive-c172x-demo-the-periphery-under-load)). `t_end` and `stop_on` are
  **defaults**, overridable per run at `run!` ([§13.5](#135-termination-is-a-state-not-an-exception), [§10.6](#106-run-lifecycle-and-partial-advance)); a run ends
  at the first grid boundary reaching or exceeding `t_end`, whole frames only
  ([§10.4](#104-shutdown-protocol)). `localization_tol` is the root-finder's relative bracket-width
  convergence test (`localization_tol · h`) and `event_budget` the per-frame
  localization allowance ([§8.4](#84-localization-mechanics)): trajectory-determining like their siblings,
  hence validated with them (`DeploymentInvalid`) and recorded in the
  deployment block, where replay compares them ([§9.5](#95-inbound-the-input-trace), [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)). Recording:
  `trace` is the input trace's plain kill switch and `log` the snapshot log's,
  with `log_every` the log's keep-every-kth decimation — admissible on the
  derived artifact only, never on the trace ([§9.2](#92-outbound-snapshot-publication), [§9.5](#95-inbound-the-input-trace), row 29) — and
  `log_max` the maximum number of retained snapshots, finite by default with
  `Inf` the opt-out: when the log fills, the retention stride doubles, so the
  whole run stays covered at coarsening density, the boundary-zero and terminal
  snapshots being retained unconditionally and outside the bound ([§9.2](#92-outbound-snapshot-publication)). Both
  are view policies, not trajectory-determining: neither enters the deployment
  block, and replay neither records nor compares them. `debug`
  gates the workspace poison ([§7.3](#73-discrete-state-modes-and-workspace)).
- `attach!(sim, dev::AbstractDevice, binding::AbstractBinding; should_abort = false)`
  — the roots are mandatory and the signature is the gate; **sides are declared**
  by the Bool traits `is_input`/`is_output` (framework defaults false on
  `AbstractBinding`), with a conformance check at attach pairing each trait
  against its method — error fallbacks for a declared side whose `claims`/`reads`
  was never written, `which`-against-the-fallback for a method defined under a
  false trait, both `BindingContractMismatch` ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)):
  `claims(b)`/`map_input(datum, b)` is the input side (the enumerated face
  set *is* the claim — what the device may write,
  not what it will — registered with exclusivity enforced, the staged
  shape and normalization shim compiled, [§9.4](#94-inbound-per-device-staging-representation-and-the-drain)); the trait
  `is_greedy(b) = true` switches the claim's *source* — the framework
  computes the unclaimed complement at attach instead of calling `claims`,
  everything downstream being identical, an empty remainder legal and
  reported (`EmptyGreedyClaim`), and `is_greedy` without `is_input` an error
  ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster), [§9.6](#96-devices-one-authoring-contract-no-taxonomy)); `reads(b)`/
  `map_output(nt, b)` is the output side ([§14.4](#144-two-application-registers-over-one-plan) selectors validated and
  compiled to one gather, [§9.2](#92-outbound-snapshot-publication)); `TableBinding` is the shipped
  data-driven binding, the standard GUI binding the shipped greedy
  one ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)). A stopped-sim operation — legal in `built`,
  `initialized` and `stopped`, an error while `running`
  (`ServiceLifecycle`, [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s roster freeze); admission checks identity
  (`AlreadyAttached` — one roster entry per instance, rebinding =
  `detach!` + `attach!`), calling-task affinity (`CallerTaskConflict` —
  at most one holder) and claims (`ClaimConflict`), [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster); registers
  only, the task appears at the next `run!`. `should_abort` is the
  per-attachment failure policy: set, the device's departure — loop body
  returning, crash, or a failed `init!` — also requests a sim stop; clear (the
  default), the run continues with the device absent and its claims held to run
  end ([§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§10.4](#104-shutdown-protocol)).
- `detach!(sim, device)` — removes the roster entry and releases the
  device's claims; stopped-sim only, like `attach!`. A loop body's
  voluntary exit or crash mid-run does *not* detach: the task dies, the
  claims persist to run end ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster), [§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§10.4](#104-shutdown-protocol)).
- The device contract — `MyDevice <: AbstractDevice` plus `init!(dev)` /
  `loop(dev, handle)` /
  `shutdown!(dev)` / optional `unblock!(dev)` / optional trait
  `needs_calling_task(dev) = false` ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)'s topology: at most one per
  roster, its loop body runs inline on the calling task): per-run `init!`
  on the calling task — bracketed, so a throw there is `shutdown!` plus
  `DeviceCrash` by name and the device is dead from boundary zero
  ([§10.4](#104-shutdown-protocol)) — the author-owned task body inside the framework's
  try/catch/finally wrapper, voluntary exit = return ([§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§10.4](#104-shutdown-protocol)).
- The device handle — one type, capabilities not taxonomy: `running`,
  `latest`, `wait_next_snapshot` ([§10.3](#103-the-next-snapshot-wait)), `stage!`, `binding`, `gather`,
  `report!` ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)).

**Condition algebra** ([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)–[§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)).

- `fragment(; x, m, slots)` — self-vocabulary payloads at the authoring
  level; `slots` names faces of that level's contract.
- `at(prefix, node)` — scoping; stores, never applies. Also lifts whole
  `TrimProblem`s and linearization tap sets ([§14.9](#149-mounting-problems-as-relocatable-values), [§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)).
- `merge(nodes...)` — symmetric collection; duplicate leaves error with dual
  provenance; mixing a node with a bare NamedTuple is an error method, not
  `Base.merge` ([§14.2](#142-fragment-composition-locality-without-schema)).
- `override(base, patches...)` — ordered layering; patch wins, provenance
  keeps both ([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)).
- `condition(comp; kw)` — the shipped fragment-function idiom; aircraft
  baselines (`ready_for_taxi(ac)`, `cold_and_dark(ac)`) are its
  full-coverage instances.

**Stopped-sim services** ([§14](#14-stopped-sim-services)).

- `init!(sim, condition; t0 = 0.0)` — slot totality checked pre-write
  ([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)), then boundary zero: project → sweep → events → due `g` updates →
  header + first snapshot ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)).
- `trim!(sim, problem; baseline, t0 = 0.0, backend) → TrimReport` —
  nonlinear least squares on the packed residuals with exact Dual
  Jacobians, against the problem's own `tolerances`
  (`residuals(reads, d) → NamedTuple`, packed in `tolerances`' field order as
  decisions pack in `guess`'s — names pair, order is the declared side's,
  [§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)'s
  closed seven-field problem); setup and commit both carry [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)'s slot-totality
  check; commit = `init!` with `override(baseline, solution)` —
  boundary zero anchored at `t0`, recordings cleared ([§10.6](#106-run-lifecycle-and-partial-advance)); resume-at-
  time = `capture`'s returned `t` as `t0`; `converged` = the service's
  per-residual box test at the backend's returned point, backend-independent
  and the commit's gate, with the backend's status and counts recorded
  diagnostically; the backend seam a pinned one-method signature,
  `solve(backend, eval!, d0, lower, upper, tol) → (; d, status, nevals, niters)`
  — in-place `eval!(r, J, d)` filling `J` only when it is not `nothing`,
  packed vectors in the declared orders, `status` an open `Symbol` recorded
  verbatim; non-convergence reports, never
  throws ([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals), [§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)).
- `capture(sim) → (condition, t)` — full-store gather including root slots;
  warm restart = capture → tweak → apply ([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults), [§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)).
- `linearize(sim, taps) → labeled (ẋ₀, x₀, u₀, y₀, A, B, C, D)` — pure query, one
  seeded Dual pass on scratch; operating point defaults to `capture(sim)`;
  taps = `get_state`/`get_slot`/`get_output` selector lists with control-design
  labels ([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)).

**Running.**

- `run!(sim; gui = false, pace = 1, margin = 0.002, t_end = <ctor value>,
  stop_on = <ctor value>)` — paced and unpaced runs bit-identical
  ([§8.7](#87-real-time-pacing)); the GUI an ordinary rostered device rendered on the calling task
  ([§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)); `gui = true` is **run-scoped attachment** — at run entry it
  attaches the standard GUI device under the standard greedy binding, with
  `should_abort = true`, **iff no GUI is already rostered**, and the shutdown
  tail detaches it again ([§10.4](#104-shutdown-protocol)), the error path included, so a
  hand-attached GUI makes the flag a
  no-op rather than an admission error and nothing the flag did survives the
  run; a persistent GUI session is spelled `attach!`/`detach!` by hand.
  Placement follows the roster, not the flag: a rostered GUI
  moves the loop to a spawned task for as long as it is rostered
  ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads), [§10.6](#106-run-lifecycle-and-partial-advance)); sugar
  never activates by default. `run!` blocks until the run ends; deviceless
  it is fully synchronous on the calling task; `init!` required first
  ([§10.6](#106-run-lifecycle-and-partial-advance)). `margin` is the single pacing
  knob ([§8.7](#87-real-time-pacing)), in seconds and defaulting to 2 ms — the sleep primitive's
  granularity plus its measured overshoot — with `0` / 2 ms / `∞` spanning the
  design space. `t_end` and
  `stop_on` override the constructor's defaults **for this run only**, with
  `stop_on` validated against the `Build` here exactly as at construction, and
  the effective pair recorded in the run metadata ([§13.5](#135-termination-is-a-state-not-an-exception)).
- `step!(sim; frames = 1) → frames_advanced` — synchronous partial advance
  through the ordinary frame sequence, bit-identical to the same frames under
  `run!`; `t_plus = <duration>` is the mutually-exclusive duration spelling
  (whole frames until the boundary time covers it); returns the frames
  *actually* advanced, fewer than requested when
  `t_end` or a `stop_on` face ended the run inside the call. Between calls the
  simulation reports `initialized`; `run!` may follow and continues from the
  current boundary; a stepping session is deviceless — write via `stage!`,
  read via `latest` ([§10.6](#106-run-lifecycle-and-partial-advance)).
- `stage!(sim, "face" => value, ...)` — task-free staging from the calling
  task into [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s harness register (surface = the currently-unclaimed
  faces): traced, drained last at the next frame top, surface-checked exactly
  as the GUI's writes ([§10.6](#106-run-lifecycle-and-partial-advance)'s harness cell; legal under `run!` and
  `step!` alike).
- `latest(sim) → snapshot` — the current published snapshot, the same
  immutable value device handles read ([§9.2](#92-outbound-snapshot-publication)); the assertion/inspection
  accessor of the harness and REPL registers ([§10.6](#106-run-lifecycle-and-partial-advance)).
- `phase_bodies(sim) → named callables` — the compiled phase bodies of the
  nominal activation, bound over the simulation's own buffers: the four
  blocks (`rhs`, `sweep_hx`, `sweep_hxu`, `ticks` — the sweeps in both
  arities, zero-arg interior and tick-indexed boundary; `ticks` takes the tick
  index) plus per-event guards/handlers and per-component `project`, keyed
  by the model's roster. The [§7.5](#75-allocation-policy-a-scoped-invariant) allocation seam: warm, then
  `@ballocated(body()) == 0` per body; diagnostic register, the one promise
  being identity with what the loop runs; isolated invocation leaves buffers
  valid but off-trajectory — re-run `init!` to continue ([§12.7](#127-the-compiled-executor)).
- Control plane — pause/un-pause, pace and `margin` changes, stop on a
  separate atomic surface, never staged ([§10.1](#101-control-plane); pacing sits outside the
  semantics, so both are safe to change live).
- Termination — model state via `stop_on` faces read at every published
  boundary ([§13.5](#135-termination-is-a-state-not-an-exception)); shutdown completes a boundary, publishes the final
  snapshot, then joins ([§10.4](#104-shutdown-protocol)).
- Post-run — the log is retained snapshots; `trace(sim) → trc` retrieves the
  always-on input trace, and `replay!(sim2, trc; to_boundary = k)` re-drives
  a fresh `Simulation(world)` bit-identically through the ordinary loop
  (boundary zero from the trace header, drain fed by frame ordinal), ending
  `initialized` — inspect via `latest`/live stores, advance via `step!`,
  continue via `run!`; the state-trajectory inspector and the `StepError`
  reproduction tool ([§9.2](#92-outbound-snapshot-publication), [§9.5](#95-inbound-the-input-trace), [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop); on-disk persistence deferred,
  [§16](#16-open-axes)).

---

## Appendix C. The diagnostic kind set

The closed set row 58 commits to, made normative: **tests match on kind plus
payload fields, never on message text** ([§13.2](#132-diagnostics-structured-values-one-carrier-exception)), so this table — not any
message — is the acceptance-test contract, and adding a kind is a decision-log
entry. Every row's payload is *in addition to* what [§13.2](#132-diagnostics-structured-values-one-carrier-exception) requires of all
diagnostics: paths and names as strings, never instances; the list-in-hand
wherever a did-you-mean renders; the didactic register (state the fix). Owning
sections stay the normative home of each rule; this appendix is an index of the
values, in the manner of Appendices A and B.

Severities, in the vocabulary [§13](#13-error-discipline) fixes:

- **build (collected)** — a diagnostic from a declarative pass, collected with its
  siblings and thrown as one `BuildError` at the stratum barrier ([§13.1](#131-reporting-policy-collect-the-checks-fail-the-evaluations-fast));
- **build (fail-fast)** — raised while *user code* runs (a boundary-connection body, a
  probe); the first one aborts the phase ([§13.1](#131-reporting-policy-collect-the-checks-fail-the-evaluations-fast));
- **service** — raised by a stopped-sim service, or by
  `attach!`/`Simulation`/`run!` validating against the `Build`; collected into one
  carrier wherever the owning section says so ([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)'s register, [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)'s
  pre-write check), a single throw at the call otherwise;
- **runtime** — fail-fast during a boundary, reaching [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)'s single catch site
  as a species of `StepError`;
- **warning (runtime)** — the per-occurrence runtime stream of [§13.2](#132-diagnostics-structured-values-one-carrier-exception),
  carried by [§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)'s per-writer diagnostic cells and rate-limited by them:
  every kind in this severity is bounded per writer per boundary (a ring of
  sixteen retained values, the excess becoming per-kind suppressed counts).
  The per-row qualifiers below record where that bound is load-bearing — a
  source that can repeat within a frame — and where the source itself fires
  once. The
  *build* warning severity exists and its set is currently empty (row 84).
- **warning (service)** — raised by a stopped-sim service call that
  *completed*: emitted at the call site through the standard logging backend,
  beside the returned value, never thrown, part of no collection; no rate
  limit — each kind fires at most once per call, and its payload is drawn from
  the report the call returns ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions), [§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)).

**Declaration and wiring** (Stratum A):

| kind | payload | owner | severity |
|---|---|---|---|
| `UnknownPort` | the wire end (`source`/`destination`), that end's path, the unknown port name, that end's port list (did-you-mean) | [§6.1](#61-connections-and-hierarchy), [§11.4](#114-failure-walkthroughs-the-error-locality-grounding) w1 | build (collected) |
| `UnconnectedInput` | leaf path, input name, declared entry type, the obligation chain's last level | [§6.1](#61-connections-and-hierarchy), [§11.4](#114-failure-walkthroughs-the-error-locality-grounding) w2 | build (collected) |
| `TwoProducers` | destination terminal, both producer terminals with provenance (sibling wire / ancestor deep route / boundary connection entry) | [§6.1](#61-connections-and-hierarchy), [§11.8](#118-computed-connections-and-generic-boundaries) | build (collected) |
| `WireTypeMismatch` | both endpoint paths, both face names, declared entry type, producer face type | [§6.1](#61-connections-and-hierarchy), [§11.2](#112-the-declaration-inventory), [§11.4](#114-failure-walkthroughs-the-error-locality-grounding) w4 | build (collected) |
| `WalkingFaceAtFrozenEntry` | consumer path and entry name, producer path and face name, the offending leaf, both declared leaf types; both remedies in the message ("declare the entry `T` if the consumer promotes; feed it from a non-walking source if the freeze is genuine") | [§6.1](#61-connections-and-hierarchy), [§11.2](#112-the-declaration-inventory) | build (collected) |
| `PathResolution` | path, offending segment, sibling field list; for a traversal past a generically-held field, that field's declared type | [§6.1](#61-connections-and-hierarchy), [§13.3](#133-build-primitives-resolve-and-the-face-list-accessors) | build (collected) |
| `AbstractAtRoot` | face name, consuming leaf path, the abstract entry; remedy hint (wire a concrete producer — in a rig, a stub child, [§13.7](#137-tooling-consequences-provenance-and-the-component-library)) | [§11.2](#112-the-declaration-inventory) | build (collected) |
| `RootSlotTypeConflict` | face name, the consuming paths, their conflicting concrete declarations at nominal (a tolerance difference is not a conflict — [§11.2](#112-the-declaration-inventory)'s meet) | [§11.2](#112-the-declaration-inventory) | build (collected) |
| `IllegalStateLeaf` | component path, `init_x` field name, leaf type, the closed vocabulary (scalar / `SArray` at the common eltype) | [§7.1](#71-continuous-state-structured-immutable-flat-backing), [§11.2](#112-the-declaration-inventory) | build (collected) |
| `StoreWithoutUpdate` | component path, the `init_x` store, the missing update (neither `f` nor `g` has a method); shadowing note when the parent module defines its own `f`/`g` ([§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)) | [§11.2](#112-the-declaration-inventory) | build (collected) |
| `EventHalfMissing` | component path, event name, which half, the function that has no method | [§11.2](#112-the-declaration-inventory) | build (collected) |
| `PrimitiveAtRoot` | root path, component type | [§11.2](#112-the-declaration-inventory) | build (collected) |
| `ClassUnreadable` | component path, type, declarations found, both family lists; did-you-mean when the type holds component-typed fields; shadowing note when the parent module defines same-named declaration functions ([§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)) | [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape) | build (collected) |
| `ClassMixed` | component path, the `child_connections` declaration and the offending leaf declarations | [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape) | build (collected) |
| `ContainerMixed` | container field path, offending element keys/indices, their types | [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape) | build (collected) |
| `DeclarationOnWrongTier` | component path, the offending declaration (`f`/`g`, `events`, `init_m`, or a `workspace`/`output_types` arity — stage names carry no tier since row 173), the tier the leaf's other declarations announce | [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), [§11.2](#112-the-declaration-inventory), [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape) | build (collected) |
| `TierSignatureMismatch` | component path, the declaration at fault (`input_types` or `output_types`), the leaf's tier, the signature form found versus the form mandated (two-argument `(::C, ::Type{T})` on the continuous tier, plain `(::C)` on the discrete); stateful leaves only — on a stateless leaf `output_types`' arity *is* the tier ([§11.2](#112-the-declaration-inventory)), so there is nothing to mismatch | [§11.2](#112-the-declaration-inventory), [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape) | build (collected) |
| `FaceNameIllegal` | assembly path, face name, the violated invariant (contains `/`) | [§11.6](#116-paths-wiring-and-faces) | build (collected) |
| `FaceNameCollision` | assembly path, face name, both entries' provenance (hand-written / computed) | [§11.6](#116-paths-wiring-and-faces) | build (collected) |
| `FaceDirectionConflict` | assembly path, the declaring method, the offending entry, the resolved port's actual direction | [§11.6](#116-paths-wiring-and-faces) | build (collected) |
| `UnknownFaceSelection` | child path, reason (unknown names / both `except` and `only` given), the offending names, the child's face list | [§11.8](#118-computed-connections-and-generic-boundaries) | build (collected) |
| `RatesViolation` | assembly path, offending key, reason (deep key / unknown child / `K` on a continuous child) | [§8.5](#85-multi-rate-tick-scheduling), [§11.7](#117-rate-scopes) | build (collected) |
| `MissingProbeValue` | face name, type | [§12.3](#123-probing-and-input-synthesis) | build (collected) |

**Schedule and contract conformance** (Strata B and C):

| kind | payload | owner | severity |
|---|---|---|---|
| `AlgebraicCycle` | the SCC's member terminals in slash form, the wires among them, optional classification (`real`/`artificial`) with the member whose hop died | [§5.5](#55-algebraic-loop-policy-reject-at-build-time), [§5.6](#56-diagnostics-feedthrough-tracing) | build (collected) |
| `ProducedByTwoStages` | component path, port name, both stage names | [§4.3](#43-table-mechanics-and-port-granularity), [§11.3](#113-visibility-the-contract-is-the-interface) | build (collected) |
| `DeclaredNotProduced` | component path, declared name, the stage-product list and the state-field list | [§11.3](#113-visibility-the-contract-is-the-interface) | build (collected) |
| `UndeclaredReturnField` | component path, stage, returned field name, candidates (`output_types`) | [§11.3](#113-visibility-the-contract-is-the-interface), [§11.4](#114-failure-walkthroughs-the-error-locality-grounding) w5 | build (fail-fast) |
| `DeadStage` | component path, stage — a stage method returning bare `(;)`, producing neither ports nor `w` | [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), [§12.3](#123-probing-and-input-synthesis) | build (fail-fast) at probe |
| `ConformanceFailure` | component path, function, field-level diff (missing / unexpected / per-field expected-vs-observed — order-insensitive, the return having been canonicalized first), simulation time | [§12.5](#125-the-always-on-conformance-check) | build (fail-fast) at probe; **runtime** as a `StepError` species |
| `GuardForm` | component path, event name, observed probe return type, both admissible forms | [§12.5](#125-the-always-on-conformance-check) | build (fail-fast) |
| `LocalizedGuardForm` | component path, event name — a localized event whose guard probes `Bool` (localization requires the sign form) | [§12.5](#125-the-always-on-conformance-check) | build (fail-fast) |
| `BundleFieldError` | component path, function family, requested field, the legal field set, classification (undeclared store / wrong-tier fact / illegal for this function family) | [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), [§13.2](#132-diagnostics-structured-values-one-carrier-exception) | build (fail-fast) at probe; **runtime** thereafter |
| `HandlerReturnKey` | component path, event name, offending key, the legal set `{x, m}` narrowed to the stores that exist | [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), [§12.5](#125-the-always-on-conformance-check) | build (fail-fast) |
| `UserCodeFraming` | component path, which function, the probe context including synthesized inputs; the original exception as `cause` | [§13.2](#132-diagnostics-structured-values-one-carrier-exception) | build (fail-fast) |

**Deployment, periphery and services:**

| kind | payload | owner | severity |
|---|---|---|---|
| `MissingInit` | the simulation's status, the entry point called (`run!`/`step!`) | [§10.6](#106-run-lifecycle-and-partial-advance) | service |
| `ServiceLifecycle` | the operation (`attach!`/`detach!`/`init!`/`trim!`/`capture`/`linearize`), the current status, the legal statuses | [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster), [§14](#14-stopped-sim-services) | service |
| `StopFaceInvalid` | face name, reason (unknown / not root-exported / not `Bool`), the root output-face list; the binding site (constructor or `run!`) | [§13.5](#135-termination-is-a-state-not-an-exception) | service |
| `DeploymentInvalid` | the deployment parameter (`h`, `n`, algorithm, `localization_tol`, `event_budget`, the harmonic-grid relation), the value in hand, the violated constraint | [§12.1](#121-three-strata) | service (collected) |
| `AttachUnknownFace` | device id, binding entry, face name, the root input-face list | [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster) | service |
| `AlreadyAttached` | the device id of the existing roster entry, its binding | [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster) | service |
| `CallerTaskConflict` | both device ids — the rostered `needs_calling_task` holder and the candidate | [§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads), [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster) | service |
| `ClaimConflict` | face name, claiming device id, incumbent device id | [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster) | service |
| `EmptyGreedyClaim` | the greedy device's id and its binding — the computed complement was empty, every root input face being claimed already | [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster), [§9.6](#96-devices-one-authoring-contract-no-taxonomy) | warning (service) |
| `BindingContractMismatch` | the binding type, the trait and the method at fault, and the direction: a declared side whose enumeration method is missing (`is_input`/`is_output` true, the root's error fallback reached), or a `claims`/`reads` method defined under a false trait (detected by `which` against the fallback); `is_greedy` without `is_input`, and `claims` defined on a greedy binding, report here too | [§9.6](#96-devices-one-authoring-contract-no-taxonomy) | service |
| `ReadBindingUnresolved` | device id, the selector, path and field, candidates; a `reason` distinguishing an unresolved path from a store selector in a snapshot binding ([§14.4](#144-two-application-registers-over-one-plan)'s source rule) | [§9.2](#92-outbound-snapshot-publication), [§14.4](#144-two-application-registers-over-one-plan) | service |
| `ConditionResolution` | entry path, store, field, offending value type and declared leaf type, provenance chain; sub-kinds: unknown path, undeclared field, unconvertible value, unexported slot face | [§14.2](#142-fragment-composition-locality-without-schema), [§14.3](#143-resolution-flatten-validate-compile-once) | service (collected) |
| `DuplicateConditionLeaf` | the leaf `(path, store, field)`, both provenance chains, the `override` advice | [§14.2](#142-fragment-composition-locality-without-schema) | service (collected) |
| `ConditionNodeMisuse` | the offending argument's type, the node kinds in hand | [§14.2](#142-fragment-composition-locality-without-schema) | service |
| `UninitializedSlots` | every uncovered root face, in declaration order | [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator) | service (collected), pre-write |
| `TapResolution` | tap set (`x`/`u`/`y`), selector kind, path, field, optional index, candidates; for a declaredly-unseedable slot, the pinning consumer's path and its `input_types` entry | [§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query) | service (collected) |
| `TrimProblemInvalid` | the offending `TrimProblem` field, the names or types in hand (a key-set or field-type mismatch; never a field-order difference) | [§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals), [§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report) | service (collected) |
| `TrimCommitEvents` | the events fired at boundary zero: component paths and event names; the same list rides the `TrimReport` | [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions), [§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report) | warning (service) |
| `TrimCommitResiduals` | the offending residual names with committed-state values and tolerances — a converged solve whose committed-state residuals violate the box test | [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions), [§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report) | warning (service) |
| `ReplayHeaderMismatch` | the mismatch, discriminated: a store or slot (component path, store, expected vs. found layout/type) or a deployment parameter (`Δt_base`/`h`/`n`/algorithm/`localization_tol`/`event_budget`, recorded vs. bound value); the build's and the trace's provenance | [§9.5](#95-inbound-the-input-trace), [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop) | service |
| `ReplayUnknownFace` | face name, frame ordinal, the trace's device tag, the root input-face list | [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop) | service (collected) |

**Runtime:**

| kind | payload | owner | severity |
|---|---|---|---|
| `StepError` | the carrier: cursor frame (component path, function, boundary phase — RK stage, event round, localization probe, tick), boundary time, frame-entry boundary index (replay pointer), original exception as `cause` | [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor) | runtime |
| `NonfiniteState` | component path, the offending state block, boundary time and index | [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor) | runtime |
| `ChatteringBudget` | component path, event name, boundary time, the exhausted `event_budget` and the frame's localization count | [§8.4](#84-localization-mechanics) | warning (runtime) |
| `EventDeferred` | component path, event name, boundary time — re-enabled within the boundary, deferred by the once-per-event rule | [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) | warning (runtime) |
| `DebtReanchor` | forgiven debt, the new schedule anchor, boundary time | [§8.7](#87-real-time-pacing) | warning (runtime) |
| `ClaimedFaceEntry` | face name, the incumbent (claiming) device id, the discarded value; the site (staging, or a stopped-sim attach's renormalization). Harness-register only — a device's out-of-surface entry is `OutOfClaimEntry` | [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster), [§9.4](#94-inbound-per-device-staging-representation-and-the-drain) | warning (runtime) |
| `OutOfClaimEntry` | device id, face name, the discarded value, the device's claim set; the incumbent's device id when the face is claimed elsewhere | [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster) | warning (runtime) |
| `ThreadBudget` | thread count, device-task count | [§10.2](#102-loop-scheduling-wait-primitive-yields-thread-budget) | warning (runtime), at `run!` |
| `DeviceJoinTimeout` | device id, the join timeout, boundary time and index at shutdown | [§10.4](#104-shutdown-protocol) | warning (runtime) |
| `DeviceCrash` | device id, the original exception as `cause`, whether `should_abort` was set; also the init-time failure, reported pre-spawn from the initialization bracket after its `shutdown!` | [§10.4](#104-shutdown-protocol), [§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor) | warning (runtime) |
| `ReplayDiscardedStaging` | device id, the discarded batch's face names, frame ordinal | [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop) | warning (runtime), repeating source — rate-limited per writer ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)) |
| `MalformedDatum` | device id, the cause exception; emitted by the author's loop body via `report!(handle, …)` | [§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor) | warning (runtime), repeating source — rate-limited per writer ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)) |
| `EntryTypeMismatch` | writer id, face name, the offending value's type, the slot's declared type, the discarded value | [§9.4](#94-inbound-per-device-staging-representation-and-the-drain) | warning (runtime) |
| `PoisonSkip` | component path, the skipped workspace stores and their element types | [§7.3](#73-discrete-state-modes-and-workspace) | warning (runtime), once per activation |
| `UnboundedRun` | the effective `t_end`, `stop_on` set and `pace`; the remedy names both, and — interactively — the operator interrupt as the sanctioned escape from the configuration warned about ([§10.4](#104-shutdown-protocol)) | [Appendix B](#appendix-b-api-synopsis-the-entry-points), [§13.5](#135-termination-is-a-state-not-an-exception) | warning (runtime), at run start |

---

## Appendix D. Glossary

*Non-normative. Each entry compresses the meaning its owning section fixes
and cites that section; where an entry and its owning section diverge, **the
owning section wins** — the same precedence rule the companion walkthrough
explainers carry. The glossary's job is to route a reader to the normative
text and to make drift visible, never to be a second source of truth. Entries
are grouped by subject and alphabetical within each group; a term appears
once, in the group that owns it, with a "not to be confused with" clause
wherever a neighboring term is genuinely close.*

### D.1 Component model and declaration layer

**abstract entry** — an `input_types` entry whose declared type is abstract,
stating **structural substitutability**: any concrete producer face below the
bound wires to it ([§4.4](#44-function-valued-signals-environment-access)'s field handles are the demonstrated client). Never
needed for eltype genericity, and illegal where the face surfaces as a root
slot (`AbstractAtRoot`) ([§11.2](#112-the-declaration-inventory)).

**assembly** — pure composition: component-typed fields as children, plus
`child_connections` (mandatory, the class marker), `input_connections`,
`output_connections` and `rates`, with no
dynamics of its own; flattened away for scheduling, retained as the
navigation hierarchy and as declaration-level rate scopes ([§3.3](#33-assembly), [§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)).

**auto-published port** — a declared output that matches a state or mode field
by name and type and that no stage produces: the framework publishes it from
the store at stage-1 position (`h_x`, either tier; the match is against
`init_x`, plus `init_m` on the continuous tier). Contract-driven — row 16
rejected blanket identity publication of state — a framework write, never a
probe product, and excluded from the stage-1 hand-down `y_x` ([§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries), [§11.3](#113-visibility-the-contract-is-the-interface),
[§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), row 169).

**class** — a component's primitive-vs-assembly status, read off *which*
well-known declarations its type defines: `child_connections` ⇒ assembly, any leaf
declaration ⇒ primitive, neither ⇒ `ClassUnreadable` ([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)). Not to be
confused with *tier* (continuous vs. discrete, [§D.4](#d4-time-and-events)) or with a diagnostic
*kind* ([§D.9](#d9-error-discipline-and-diagnostics)).

**component** — the unit of modeling: a leaf (continuous or periodic discrete
primitive) or an assembly of components; "primitive" and "leaf" are used
interchangeably for the non-assembly classes ([§3](#3-component-taxonomy)).

**container children** — a `Tuple`/`NamedTuple` field whose elements are all
components, contributing them as children path-named `"field/1"` or
`"field/key"`. Transparent grouping, not an assembly: no contract, no
`child_connections`, no rate scope ([§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)).

**continuous component** — the hybrid primitive: continuous state `x`, modes
`m`, flow `f`, two output stages, events (guards + handlers) and optional
`project`; any facet may be empty, so a state-free instance is an FSM ([§3.1](#31-continuous-component-the-hybrid-primitive)).

**contract** — a component's declared interface: `input_types` (its
requirements, read permissively — what each entry *allows* to arrive) and
`output_types` (its public ports, read literally — what each cell *carries*).
Both take the two-argument `T`-form on the continuous tier and the plain form on
the discrete one. Declared in
`output_types` = public, returned in `w` = private by construction, returned in
`y` and declared nowhere = build error ([§11.2](#112-the-declaration-inventory), [§11.3](#113-visibility-the-contract-is-the-interface)).

**declaration inventory** — the closed set of well-known functions a component
or assembly defines — `init_x`/`init_m`, `workspace`,
`input_types`/`output_types`, `events`, the stages, `f`/`g`/
`project`, and `child_connections`/`input_connections`/`output_connections`/`rates` — each declared in a stated
register of authority: by value, by type, by allocation ([§11.2](#112-the-declaration-inventory)).

**function family** — which bundle fields a given function may legally
receive: `h_x`/`h_xu`/`f`/`g`/guard/handler/`project` (the `h_*` sets being
per-tier, [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)), with
[§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)'s comment block stating each family's maximal legal set and
`BundleFieldError` classifying a read as illegal for the family ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)).
Not a diagnostic *kind* ([§D.9](#d9-error-discipline-and-diagnostics)).

**generic holding** — a parent holding a child through a non-concrete field
type; the child is opaque below its faces, and the wires and boundary connections
referencing those faces *are* the imposed contract, checked per instantiation
([§11.8](#118-computed-connections-and-generic-boundaries), [§6.1](#61-connections-and-hierarchy)).

**hybrid causal system** — what the framework simulates: continuous flow with
algebraic outputs, multi-rate periodic discrete dynamics, zero-crossing
events, post-step manifold projection, and externally injected inputs ([§2](#2-formalism)).

**the letters** — `f` the continuous flow, `g` the discrete update, `h_*` the
output stages (suffix = dependence class, not argument list), `x` the state on
either tier and `m` the continuous-only mode store, `u` wired inputs, `y` own
published signals, `ws` the
workspace. Bare `h` means the integration step size only ([§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries), [§8](#8-time-and-execution));
bare `z` means only the shift operator `z⁻¹` (row 173).

**periodic discrete component** — a leaf with state `x`, update `g`
at a declared rate, and two output stages whose cells hold zero-order between
ticks; it has no `m` store, and its `x` reaches others only through signals
([§3.2](#32-periodic-discrete-component)).

**rate scope** — an assembly's `rates` declaration: immediate child name ⇒
integer multiplier `K ≥ 1` relative to the enclosing scope, composing
multiplicatively down the tree and compiled to one absolute tick divisor per
discrete component ([§11.7](#117-rate-scopes), [§8.5](#85-multi-rate-tick-scheduling)).

**schema authority** — the principle that declarations *define* structure and
evaluation only *checks* conformance against them, never the reverse; types by
declaration, values by execution, conformance by comparison ([§11.1](#111-position-a-declarative-trait-layer--plain-julia-no-macros)).

**stage function / two-stage outputs** — every component provides exactly two
output stages, `h_x` (no `u` in the bundle, hence structurally no
feedthrough) and `h_xu`; feedthrough is thereby declared by signature,
with no dependency annotations anywhere ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)).

**workspace** — component-declared mutable scratch, declared *by allocation*
(`workspace(::C, ::Type{T})` continuous, `workspace(::C)` discrete), arriving
as the `ws` bundle field; excluded from state semantics, never a condition
target, and poisoned before every call in debug mode ([§7.3](#73-discrete-state-modes-and-workspace)).

### D.2 Signals and data homes

**buffer** — the framework-owned contiguous `Vector{T}` backing all continuous
state, laid out at build time; authoritative, with typed state values as
ephemeral reconstructions of it ([§7.1](#71-continuous-state-structured-immutable-flat-backing)). [§13.6](#136-abnormal-shutdown-one-tail-two-entries)'s integration intermediates live
in framework-owned integrator buffers, never in a component's workspace.

**bundle** — the single `NamedTuple` of zero-copy views a component function
receives beside the component itself. Under the bundle law a name is present
**iff** the corresponding store or fact exists for that component; undeclared
stores are absent, never `nothing`-filled ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)).

**cell** — one concretely-typed entry of the signal table, one per output port
of the flattened model, written by its producing
stage and read by every gatherer ([§4.1](#41-immutable-value-semantics)). Every cell is public, private
intermediates never being cells ([§11.3](#113-visibility-the-contract-is-the-interface)). Bare "cell" is only this — see
*staging cell* ([§D.6](#d6-runtime-periphery)) and *store*.

**constant source** — an ordinary library component with no inputs and no
state, publishing a value its instance holds (`Constant{V}`); the spelling for
an aggregate input with zero contributors and for the rig stub feeding an
abstract face. Its value is instance data, never a default ([§13.7](#137-tooling-consequences-provenance-and-the-component-library), [§6.2](#62-aggregation-explicit-summing-junctions)).

**entry** — never used bare: the spec's compounds are table entry (a cell),
input entry ([§11.2](#112-the-declaration-inventory)), executor entry ([§12.7](#127-the-compiled-executor)), roster entry ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)), batch entry
([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)) and condition entry ([§14.3](#143-resolution-flatten-validate-compile-once)), each a different thing.

**face** — the name a port wears on its component's boundary: for a leaf the
port's own name, for an assembly a name declared in `input_connections` or
`output_connections`, aliasing an interior port. An opaque token with
two build-checked invariants (no `/`, unique within the assembly), its type
derived from its internal endpoint and its direction declared by the method
that names it. The periphery's write side
speaks face names only; the read side speaks them wherever it wants contract
rather than structure ([§11.6](#116-paths-wiring-and-faces), [§9.2](#92-outbound-snapshot-publication), [§14.4](#144-two-application-registers-over-one-plan)).

**feedthrough** — an instantaneous input→output dependence. **Structural
feedthrough** is this design's version: fixed by which stage produces a port
rather than annotated, with every stage-2 output conservatively presumed
dependent on every wired input ([§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)).

**field handle / function-valued signal** — an immutable query object carried
on an ordinary port (`ISAField`, `TerrainField`) that consumers evaluate at
arguments of their own choosing; bulk data rides as build-time-frozen
references, never as mutable caches ([§4.4](#44-function-valued-signals-environment-access)).

**immutable value semantics** — the signal rule, stated precisely as
immutability *plus frozen references* (`isbits` is the common case, not the
rule): no aliasing, safe concurrent reads, and a definite per-cell freshness
tied to the producer's schedule position ([§4.1](#41-immutable-value-semantics)).

**one home per datum** — buffer for continuous `x`, stores for discrete `x` and for `m`, table for
produced signals; no store mirrors another, and the table never holds
transported data ([§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws), [§7.1](#71-continuous-state-structured-immutable-flat-backing)).

**port** — the addressable unit of the model: one declared name, one cell, one
root slot, one staged write, one device claim, one trace address, one GUI
liveness verdict. Wiring is port-granular, and which stage computes a port is
invisible outside the component ([§4.2](#42-consumers-see-ports-not-stages), [§4.3](#43-table-mechanics-and-port-granularity)).

**signal table** — the framework-owned collection of cells holding every
produced signal of the flattened model; consumers gather views from it, and
its consistency is a boundary property, transiently integrator scratch within
a step ([§4.1](#41-immutable-value-semantics), [§8.3](#83-signal-table-consistency-is-a-boundary-property)).

**slot** — a root input slot: the root assembly's own input face,
produced by no component, constant within a frame, and the only thing the
periphery may write ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster), [§11.6](#116-paths-wiring-and-faces)).

**staging cell** — the per-device atomic holding place where a device's pending
write batch waits between drains; mutated frame by frame, hence outside the
table's publish-once discipline ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain)). Not a table cell ([§4.1](#41-immutable-value-semantics)).

**store** — the typed home of `m` and of a discrete leaf's `x`: overwritten by the framework when a
handler or update returns a new value, never arithmetic-touched, snapshot-free
to copy. Never called a cell — root slots, by contrast, *are* source cells of
the table ([§7.3](#73-discrete-state-modes-and-workspace), [§4.1](#41-immutable-value-semantics), [§9.2](#92-outbound-snapshot-publication)).

**summing junction** — an ordinary library component performing N-to-1
aggregation through explicit wires (`SumJunction{W, N}` or a named
site-specific variant); there is no framework aggregation mechanism, and fold
order is the junction's positional input order ([§6.2](#62-aggregation-explicit-summing-junctions)).

**value-level constructor** — the plain exported function (component, input
values) → field handle that every field-emitting component is obliged to
provide, its own swept output stage being a one-line call to it; the device by
which [§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults) condition math queries the environment before any sweep exists
([§4.4](#44-function-valued-signals-environment-access)).

**view** — a zero-copy reconstruction of a store handed to a function through
its bundle; it materializes in the caller's frame for the duration of the call
and is value-identical on re-materialization within a sweep ([§7.1](#71-continuous-state-structured-immutable-flat-backing), [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)).

**`w` (private intermediates)** — the optional second slot of an output stage's
return, `(y, w)`: an `isbits`-leaf NamedTuple of values the component computes
for its own later use. Private by construction — no cell, no contract entry,
nothing to wire, list, filter or address — travelling exactly one hop by
[§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)'s law, SSA-passed within a fused pass, probe-observed in type and
checked at the nominal activation only ([§11.3](#113-visibility-the-contract-is-the-interface), [§12.5](#125-the-always-on-conformance-check)). The inspection path
for one is promotion to a declared output.

### D.3 Evaluation and scheduling

**algebraic loop** — a genuine cycle in the instantaneous dependency graph: a
build error naming the path in canonical slash form, broken by the author with
inserted dynamics, an explicit `UnitDelay` or restructuring ([§5.5](#55-algebraic-loop-policy-reject-at-build-time)). Not to be
confused with an **artificial loop**, port-level acyclic but stage-level
cyclic, whose remedy is a ladder — the two-stage split, contract re-factoring,
and as residual a component split ([§5.4](#54-artificial-loops-and-the-escape-hatch)).

**flow / RHS** — `f`, the continuous derivative function. Evaluating the RHS
means running the whole sweep, since `f` reads the fresh table: there is no
incremental `f`-only re-evaluation ([§3.1](#31-continuous-component-the-hybrid-primitive), [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)).

**frame** — one iteration of the loop — drain, integrate, boundary sequence,
publication — the unit `step!` counts and the trace's ordinal key ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)).
Distinct from a *boundary* ([§D.4](#d4-time-and-events)), and from the kinematic reference frames of
the aircraft domain, which always appear compounded ("the b frame").

**projection** — the optional per-component hook `x ← project(x)`, run in the
only two schedule positions between a state write and its decode (after
integration, after a handler's `x`-reset); the cheap end of geometric
integration's projection methods ([§2](#2-formalism), [§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)).

**schedule** — the static evaluation order computed once at build time from
wiring edges plus intra-component feedthrough: all stage-1 functions in any
order, stage 2 in topological order, then `f`. The hot loop runs a flat list of
`(component, stage)` entries, with zero runtime graph logic ([§5.1](#51-the-scheduling-problem)).

**sweep** — one execution of that schedule against the current state, in one of
two statically distinct variants compiled from the same entry list: the
**interior sweep** walks continuous entries only — what RK stage evaluations and
localization guard probes run, so discrete cells hold ZOH mid-step by
construction — and the **boundary sweep** walks the full list, with the
boundary's due discrete entries gated in by counter modulo. Mid-step sweeps are
integrator scratch; the boundary sweep restores table consistency, and the event
phase re-runs whole boundary sweeps, against that boundary's fixed due set,
until quiescence ([§5.3](#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries), [§8.3](#83-signal-table-consistency-is-a-boundary-property), [§8.5](#85-multi-rate-tick-scheduling), [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)).

### D.4 Time and events

**boundary** — a published consistency point: where the [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event) macro-sequence
completes and a snapshot goes out. Every grid point is a boundary, but `t*`
and boundary zero are boundaries that are not frame tops ([§8.4](#84-localization-mechanics)). *Boundary
zero* ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions), [§D.8](#d8-stopped-sim-services-and-the-condition-algebra)) is a hyponym — it is an ordinary boundary whose incoming
transitions are authored rather than computed.

**boundary-detected** — the default detection policy: guards are checked for
not-holding → holding edges against their priors at step boundaries only, with
no root-finding and no step rejection, the handler firing at the end of the
step in which the edge was observed ([§2.1](#21-events-two-detection-policies)).

**chattering / event budget** — the bounded per-step localization allowance,
`event_budget`, a `Simulation` deployment keyword defaulting to 8;
exhaustion *degrades* rather than throws — localization stops for the rest of
the frame and further crossings fire at the next boundary, under a
`ChatteringBudget` warning naming the event ([§8.4](#84-localization-mechanics)).

**`Δt_base`** — the base tick period, an integer multiple `n·h` of the
continuous step and fixed at `Simulation` construction; every discrete
component's period is an integer multiple of it ([§8.5](#85-multi-rate-tick-scheduling)).

**due** — a discrete component is due at a boundary when its compiled absolute
divisor admits that boundary's tick index; due components' output stages are
gated into the *boundary* sweep (never the interior one) and their `g` updates
run after quiescence. The due set is a property of the boundary, fixed for its
whole event iteration: the modulo image of the frame index at a frame top, empty
at `t*`, everything at boundary zero ([§8.5](#85-multi-rate-tick-scheduling), [§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)).

**edge semantics / holding** — an event fires on a not-holding → holding
transition of its predicate, never on a bare sign change; the opposite
crossing direction is declared as a second event with the negated guard ([§2.1](#21-events-two-detection-policies),
[§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)).

**guard** — the declared function defining an event's predicate, evaluated
against the fresh boundary table and paired with a handler in an ordered,
named `events` collection; its detection policy is set per event by the
`localize` flag ([§2.1](#21-events-two-detection-policies), [§11.2](#112-the-declaration-inventory)).

**harmonic grid** — the rule that every discrete period is an integer multiple
of `Δt_base`, itself an integer multiple of `h`, so ticks land only on step
boundaries; grid times are indexed from the frame count, never accumulated
([§8.5](#85-multi-rate-tick-scheduling), [§8.4](#84-localization-mechanics)).

**interpolant** — the lazily built cubic Hermite continuous extension over the
last completed step, from which localization probes read the states they sweep;
invalidated at `t*`, where the handlers have made it a lie ([§8.4](#84-localization-mechanics)).

**localized** — the opt-in per-event detection policy: the crossing instant is
bracketed by derivative-free root-finding over probe sweeps of interpolated
states, to a bracket narrower than `localization_tol · h` (a deployment
keyword, default `1e-6`). Requires the continuous (sign) guard form, and runs
identically paced or unpaced ([§2.1](#21-events-two-detection-policies), [§8.4](#84-localization-mechanics)).

**once per event** — the rule bounding the boundary event iteration: each
declared event fires at most once per boundary, which makes termination
structural; an event legitimately re-enabled within the boundary waits one
step, with an `EventDeferred` warning ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)).

**pacing / pacer debt** — the pacer inserts waits between completed frames and
never alters the boundary sequence; a frame exceeding its wall budget leaves
**debt** that later frames repay, with excess forgiven by re-anchor plus
warning ([§8.7](#87-real-time-pacing)).

**predicate** — what a guard defines: a `Bool`-valued form, or the sign of a
continuous function with positive = holding (writing the sign value `σ`,
holding = `σ ≥ 0`) ([§2.1](#21-events-two-detection-policies)). Not to be confused with [§14](#14-stopped-sim-services)'s *condition*, the
value that sets a build's state ([§D.8](#d8-stopped-sim-services-and-the-condition-algebra)).

**prior** — the per-event stored sample of its predicate at the previous
boundary's quiescence (a deferred event's prior is recorded not-holding, so
the deferral resolves at the next boundary), held in loop state and never in a state store; "newly fired"
is defined against it, and boundary zero establishes every prior as
not-holding ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)).

**quiescence** — the fixed point of the boundary event phase: rounds of
[sweep → guards → handlers] iterate until a round fires nothing, after which
the priors are updated and due `g` updates run ([§8.6](#86-event-iteration-at-boundaries-to-quiescence-once-per-event)).

**remainder step** — the integration from `t*` to the original grid target
after a localized event, with `h′` derived at use; guards are re-checked on it
under the event budget ([§8.4](#84-localization-mechanics)).

**`t*`** — the localized event time: the holding endpoint of the root-finder's
final bracket, structurally strictly later than `tₙ`. A full boundary runs
there, but no ticks are due and no staged inputs are drained ([§8.4](#84-localization-mechanics)).

**tick** — an instant at which a discrete component's stages and update run,
gated by counter modulo against the harmonic grid inside the boundary sweep;
different boundaries therefore run different subsets of the schedule ([§8.5](#85-multi-rate-tick-scheduling)).

**tier** — the continuous or discrete side of the hybrid formalism, read off a
leaf's declaration shape (`DeclarationOnWrongTier` names a violation) ([§11.2](#112-the-declaration-inventory),
[§11.5](#115-assembly-declaration-type-based-class-by-declaration-shape)). Bare "tier" means only this: the genericity classes are *walked /
pinned / exempt* ([§D.5](#d5-build-pipeline)) and the detection policies *boundary-detected /
localized*.

### D.5 Build pipeline

**activation** — a re-run of Stratum C at a given scalar type `T`: cells
re-typed (producer-fed ones by evaluating the producer's output declaration at
`T`, root slots by evaluating the consuming `input_types` entry at `T`, the
state type by the leaf walk), buffers re-laid-out,
workspace allocators re-invoked, probe chain re-run. Structure and schedule are `T`-independent;
non-nominal activations are lazy, with an opt-in exhaustive set for CI ([§12.4](#124-activations-executable-sets-laziness-caching)).

**always-on conformance check** — the probe's comparison left permanently in
place: one type test of a stage return against the declaration-derived
expected `NamedTuple` at the table-write point, preceded by a type-level
reorder to that type's field order (the names pair; order carries no
semantics), both folding to zero instructions
when the return type is proven ([§12.5](#125-the-always-on-conformance-check)).

**`Build`** — the artifact `build(world)` produces: wire list, face table with
provenance, schedule and root slots as plain printable data — the inspectable
contract of the instantiation, and what `attach!`, `stop_on`, replay and
condition resolution all validate against ([§12.2](#122-the-build-artifact)).

**chunking** — splitting a large phase body's entry tuple into statically
typed chunks behind non-inlined function barriers; the implementation's only
representation freedom, converting compile cost from superlinear in body size
to linear in entry count ([§12.7](#127-the-compiled-executor)).

**executable set** — the function set an activation can actually run, hence
exactly what it probes: a `Dual` activation sees only the continuous output
stages and `f` — never the discrete stages, `g`, guards or handlers ([§12.4](#124-activations-executable-sets-laziness-caching)).

**executor** — the compiled execution form of the schedule: a concretely-typed
tuple of entries over statically typed cell storage, traversed by a
compile-time-unrolled walk, with code-selecting facts in type parameters and
plain data in fields ([§12.7](#127-the-compiled-executor)).

**leaf walk** — the framework's derivation of per-activation types from a
declared nominal type: real leaves and `Real` type parameters follow the
activation scalar, everything else pins. It survives on the **state** side alone
(the type derived from a continuous leaf's `init_x`; `init_m` and a discrete
leaf's `init_x` pin wholesale). **Cells are not
walked**: an output cell comes from evaluating the producer's `output_types` at
the activation scalar (row 166) and a root-slot cell from evaluating the
consuming `input_types` entry at it (row 167), participation and tolerance
authored per leaf in both ([§11.2](#112-the-declaration-inventory); applied in [§12.1](#121-three-strata)'s
Stratum C).

**lens (`Getter`)** — the compiled navigation step of a condition entry: its
tree position tuple lifted to a type parameter, giving type-stable access to
the authored value at apply time ([§14.3](#143-resolution-flatten-validate-compile-once)).

**measurement seam / phase bodies** — `phase_bodies(sim)` returns the compiled
bodies of the nominal activation bound over the simulation's own buffers
(`rhs`, `sweep_hx`, `sweep_hxu` — the sweeps in both arities, zero-arg interior
and tick-indexed boundary — `ticks`, plus per-event guards and handlers
and per-component `project`). Its one promise is identity with what the loop
runs, which is what makes [§7.5](#75-allocation-policy-a-scoped-invariant)'s allocation assertions honest ([§12.7](#127-the-compiled-executor)).

**nominal** — the `Float64` activation, and of a declaration its `Float64`
face (for a continuous producer's output declaration, its evaluation at
`Float64`); the only activation that runs in real time, and the one where the
conformance check demands exact type match ([§11.2](#112-the-declaration-inventory), [§12.4](#124-activations-executable-sets-laziness-caching), [§12.5](#125-the-always-on-conformance-check)).

**probe** — the build's single evaluation of a user function with real values,
checking shape and type conformance and discarding the result. Every user
function is probed once, at the initial state; probes see only that state's
branch ([§12.3](#123-probing-and-input-synthesis)).

**probe value / input synthesis** — `probe_value(::Type)` fabricates values
for the one kind of terminal with no producer, root slots
(`zero(T)`/`false`/first enum/`T()`, overridable). Strictly probe-scoped:
never an initial slot value, which [§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator) makes a structural barrier ([§12.3](#123-probing-and-input-synthesis)).

**`ProbeDual`** — the framework's exported canonical concrete probe scalar
(`ForwardDiff.Dual{ProbeTag, Float64, 1}`), which keys the CI activation
pinning walked-leaf genericity; its width is arbitrary, since what CI pins is
genericity, not a particular Jacobian ([§12.4](#124-activations-executable-sets-laziness-caching)).

**schema vs. layout** — the two lookup families the `Build` supplies to
condition resolution: *schema* is the evaluated declarations (may you write
this field, at what leaf type — the authority), *layout* is where it
physically lives (buffer ranges, store and slot indices) ([§14.3](#143-resolution-flatten-validate-compile-once)).

**stratum** — one of the build's three phases: A structure (pure declaration
reading), B schedule (the single evaluation-feeds-structure step), C
activation (everything type-shaped). Strata are barriers — a stratum that
produced any error-severity diagnostic throws before the next begins ([§12.1](#121-three-strata),
[§13.1](#131-reporting-policy-collect-the-checks-fail-the-evaluations-fast)).

**walked / pinned / exempt** — the eltype-genericity classes: walked
payload/value types follow the activation scalar, pinned parameters and
definitions stay `Float64`, and the discrete side is exempt. Enforced by the
leaf walk on the state side and stated per leaf in a continuous leaf's contract
declarations on the cell side — `output_types` for what a producer's cells
carry, `input_types` for what a consumer's entries tolerate ([§7.2](#72-numeric-genericity-eltype), [§11.2](#112-the-declaration-inventory)).

### D.6 Runtime periphery

**bad datum** — a datum unmappable for environmental reasons (truncated
datagram, malformed JSON, out-of-range field): tolerated *in the loop body* —
catch, stage nothing, `report!(handle, MalformedDatum(cause))`, continue —
while any other exception propagates and becomes `DeviceCrash`. The
classification is the device author's ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)).

**batch** — a device's staged set of face ⇒ value writes, coalesced in its
staging cell and applied whole at the next drain ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain)). The word means only
this; error reporting *collects* ([§D.9](#d9-error-discipline-and-diagnostics)).

**binding** — the value passed at `attach!` that makes a device
framework-legible: a subtype of `AbstractBinding` declaring its sides by the
Bool traits `is_input`/`is_output` (false by default on the root), with
`is_greedy` switching the input side's claim source from returned to computed
(the unclaimed complement, in place of `claims`). `claims` and `reads` carry
error fallbacks on the root, and attach cross-checks each trait against its
method in both directions (`BindingContractMismatch`); `map_input`/`map_output`
are loop-idiom conventions the framework never calls. Every input-side binding
stakes a claim; `TableBinding` is the shipped
data-driven one ([§9.6](#96-devices-one-authoring-contract-no-taxonomy), [§9.4](#94-inbound-per-device-staging-representation-and-the-drain)).

**boundary counter** — the monotonic count of *published boundaries* carried
in the snapshot and mirrored in the loop state the wait predicate tests;
incremented after the `latest` release-store, so a waking waiter can never see
a stale snapshot ([§10.3](#103-the-next-snapshot-wait)).

**calling task** — the task that invoked `run!`. It runs the loop itself (the
unattended register) unless a `needs_calling_task` device is rostered, in
which case it runs that device's loop body inline and the loop moves to a
spawned task ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)).

**claim** — the set of faces a device *may* write, registered at attach —
either returned by its binding's `claims` or computed as the unclaimed
complement under `is_greedy` — and released at detach; claiming an
already-claimed face is an attach-time error (`ClaimConflict`), and a broad
claim costs GUI liveness ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)).

**coalescing** — the CAS merge keeping one pending batch per device:
untouched faces survive, re-staged faces take the newest level (the per-face
ZOH). Its outbound mirror is newest-wins snapshot delivery ([§9.4](#94-inbound-per-device-staging-representation-and-the-drain), [§10.3](#103-the-next-snapshot-wait)).

**control plane** — the separate few-word atomic surface carrying pause,
un-pause, pace, `margin` and stop, consulted at frame top and inside the
loop's wait and pause states; structurally not staging, since a paused loop
drains nothing ([§10.1](#101-control-plane)).

**derived liveness** — the rule that a GUI widget is live iff its port's feed
chain terminates in a root slot inside the GUI's own claim in the run's frozen
surface partition; baked once at run start, with no per-port "GUI-controlled"
marking anywhere ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)).

**device** — any attached participant in the periphery: a subtype of
`AbstractDevice` under one authoring
contract (`init!`/`loop`/`shutdown!`, optional `unblock!` and
`needs_calling_task`) and one handle; input-only and output-only are
degenerate uses, and the GUI is an ordinary device ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)).

**diagnostic cell** — the single-writer cell each rostered device and the loop
itself owns for runtime diagnostics and liveness: a bounded ring (capacity 16)
of diagnostic values plus per-kind suppressed counts — the bound being the
rate limit itself — and an atomic heartbeat timestamp, taken by the loop with
`atomicswap` at the frame-top drain and frozen into the published status
([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)).

**drain** — the single point at the top of each frame where the loop takes
each staging cell by `atomicswap` and applies it through the attach-compiled
scatter, in attachment order; never at a `t*` boundary, and under the roster
freeze it performs no checks at all ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads), [§9.4](#94-inbound-per-device-staging-representation-and-the-drain)). The diagnostic
cells are taken at the same point ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)).

**framework status** — the concrete frozen value each snapshot carries beside
the signal table: [§8.7](#87-real-time-pacing)'s pacer diagnostics plus, per writer, this
boundary's drained diagnostics (`recent`), the counts the ring refused
(`suppressed`), the loop's cumulative per-writer × per-kind counters copied in
(`totals`) and the liveness timestamp ([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell), [§9.2](#92-outbound-snapshot-publication)).

**greedy claim** — the claim a binding declaring `is_greedy` receives: the
unclaimed complement computed by the framework at attach instead of returned
by `claims`, ordinary in every respect afterwards; an empty remainder is legal
and reported (`EmptyGreedyClaim`), and the shipped GUI binding is the shipped
instance ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster), [§9.6](#96-devices-one-authoring-contract-no-taxonomy)).

**harness cell** — the always-present staging cell of the harness register,
written by `stage!(sim, "face" => value, …)` from the calling task itself:
ordinary batches, traced and surface-checked, drained last by convention
([§10.6](#106-run-lifecycle-and-partial-advance), [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)).

**harness register** — the framework-owned write path of the calling task —
`stage!(sim, …)` and its cell — and the design's sole *derived* surface: the
unclaimed complement, the faces no rostered device claims, recomputed at every
stopped-sim roster change; a write to a claimed face is `ClaimedFaceEntry`
([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster), [§10.6](#106-run-lifecycle-and-partial-advance)).

**`latest`** — the `@atomic` reference a published snapshot is release-stored
into and readers acquire-load; `latest(sim)` hands the calling task the same
immutable value a device handle gets ([§9.2](#92-outbound-snapshot-publication)).

**next-snapshot wait** — `wait_next_snapshot(handle)`: the boundary counter
plus one `Threads.Condition` under the canonical predicate loop
(`counter > last_seen && running`), newest-wins, no queues, no per-frame reset
([§10.3](#103-the-next-snapshot-wait)).

**operator interrupt** — Ctrl-C in an interactive session, read as a
control-plane stop rather than a failure: delivery is masked across the boundary
macro-sequence (`disable_sigint`) and raised at a frame-top or wait unmask
point, so the run takes the ordinary graceful tail and ends `stopped`. A second
one collapses the device joins; outside the REPL, SIGINT still kills the process
([§10.4](#104-shutdown-protocol), [§10.1](#101-control-plane)).

**orphaned claims** — the claims of a device whose task died mid-run. Death is
not detach: the roster entry and claims persist to run end, the slots hold
their last-drained values, and the GUI renders the fact in the widget's
provenance; recovery is between runs ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)).

**peek** — the GUI display rule: a widget shows its own pending write if any,
else the snapshot value. Own-cell only, which is what makes multi-click
counting and paused editing correct ([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)).

**periphery** — everything outside the loop that exchanges data with it — GUI,
input devices, network I/O, logging — together with the concurrency model
binding them: staged writes inbound, snapshot reads outbound, control on its
own surface ([§9](#9-runtime-periphery-the-data-plane), [§10](#10-runtime-periphery-lifecycle-and-orchestration)).

**roster** — the list of attached device entries (binding, claims, stable
device id, attachment order): a plain immutable value the loop reads once at
`run!`, since `attach!`/`detach!` are stopped-sim operations ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)).

**scenario component** — the home of a sim-time script under the mid-run
mutation doctrine: an ordinary periodic discrete component executed
synchronously in the loop, deterministic paced or unpaced and replayed by
recomputation. The clock is the criterion — wall-clock interactions are
devices ([§10.5](#105-scripts-and-the-mid-run-mutation-doctrine)).

**selector (read-selector family)** — the closed set of deferred reads
`get_state`/`get_deriv`/`get_output`/`get_slot`/`get_face`, each
resolving against a source (table sources — a boundary snapshot or a service
evaluation's scratch tables — vs. live stores) before any client policy
applies ([§14.4](#144-two-application-registers-over-one-plan)).

**`should_abort`** — the per-attachment failure policy, an `attach!` keyword
defaulting to `false`: set, a device's departure — loop body returning, crash,
or a failed `init!` — also requests a control-plane stop; clear, the run
continues with the device absent and its claims held to run end. An attachment
fact, never a device property, the same device being advisory in one deployment
and load-bearing in another; the shipped GUI attaches with `true` ([§9.6](#96-devices-one-authoring-contract-no-taxonomy),
[§10.4](#104-shutdown-protocol)).

**snapshot** — the immutable per-boundary publication: boundary-consistent
signal table (root slots included), `t`, boundary index and
framework status. It deliberately carries no state stores — the state
trajectory is derived data ([§9.2](#92-outbound-snapshot-publication)).

**stage-on-interaction** — the GUI staging contract: value widgets stage the
new level on edit, edge widgets on activation as a level computed from the
peek; held buttons do not re-stage, and no widget stages per render pass
([§9.7](#97-the-gui-write-path-port-resolution-peek-staging-contract)).

**unattended run** — a run with empty staging and no snapshot readers: the
same loop, fully synchronous on the calling task, rethrowing after the
shutdown tail so CI fails honestly ([§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads), [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)).

**write surface** — the set of faces a writer's batch entries may reach: a
device's claim set, whether returned by `claims` or computed under
`is_greedy` ([§9.6](#96-devices-one-authoring-contract-no-taxonomy)), and for the harness register the derived
unclaimed complement. Static per run and enforced entirely at
staging — `OutOfClaimEntry` for a device, `ClaimedFaceEntry` for the harness
([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)).

### D.7 Recording and replay

**decimation** — the log's keep-every-kth retention policy (`log_every`),
admissible on the log alone because it is derived data; every boundary still
runs, publishes to live readers and enters the trace. Bounded by `log_max`, the
maximum number of retained snapshot references (default finite, `Inf` the
opt-out): when the log fills, the effective stride doubles — *progressive
re-decimation*, so coverage stays global at `log_every · 2^k` instead of
collapsing to a rolling window — with the boundary-zero and terminal snapshots
retained unconditionally and outside the bound. A view policy throughout, never
trajectory-determining ([§9.2](#92-outbound-snapshot-publication)).

**frame ordinal** — the trace's key: replay applies the recording's batches
for frame *k* at frame *k*, exact because the frame sequence is itself
deterministic under replay ([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop), [§9.1](#91-no-shared-mutable-model-staged-writes-snapshot-reads)).

**log** — the retained sequence of published snapshots (the same objects, no
copies), with a plain kill switch and `log_every` decimation; derived data,
recomputable from the trace by replay ([§9.2](#92-outbound-snapshot-publication)).

**recorders** — the trace and the log jointly, cleared together at `init!` and
at a trim commit so they restart with the run they record ([§10.6](#106-run-lifecycle-and-partial-advance), [§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)).

**replay** — the ordinary loop with exactly two substitutions: boundary zero
from the trace header, and a drain reading the trace by frame ordinal. It
re-records, ends `initialized`, and validates the header — stores, slot
faces and deployment block — up front, applying the header's `t₀`
([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)).

**run metadata** — the trace header's deployment block: `t₀`, `Δt_base`,
`h`, `n`, the algorithm identifier, `localization_tol`, `event_budget` and the
effective `t_end`/`stop_on`
pair ([§9.5](#95-inbound-the-input-trace), [§13.5](#135-termination-is-a-state-not-an-exception)).

**trace** — the primary record of a session: the sequence of drained,
device-tagged batches per frame, plus its header. On by default, because the
log is recomputable from the trace and never the reverse ([§9.5](#95-inbound-the-input-trace)).

**trace header** — the trace's preamble: the resolved initial stores
`(x, m)`, the initial root-slot values, each writer's face-name →
position schema, and the deployment block — captured after `apply!` and the
slot writes, before the boundary-zero sequence runs ([§9.5](#95-inbound-the-input-trace), [§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)).

**trace record** — the retained form of a drained batch, uniform for every
writer: (position ⇒ value) pairs for the non-`nothing` entries, converted at
the drain against the header's schema, so trace size tracks information
rather than surface width and consumers meet one format and one replay path
([§9.5](#95-inbound-the-input-trace), row 176).

**what-if register** — replaying a trace against the same structure with
changed parameters: deterministic re-driving of the recorded inputs through a
modified model, promising determinism but never bit-identical reproduction
([§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)).

### D.8 Stopped-sim services and the condition algebra

**`at` / `Scoped`** — the scoping combinator: `at(prefix, node)` stores a
prefix beside a condition node and applies nothing, path concatenation
happening once at resolution. It also lifts whole `TrimProblem`s and
linearization tap sets ([§14.2](#142-fragment-composition-locality-without-schema), [§14.9](#149-mounting-problems-as-relocatable-values)).

**baseline** — an aircraft-shipped, full-coverage condition function
(`ready_for_taxi(ac)`, `cold_and_dark(ac)`) layered under tweaks by
`override`, and the `baseline` keyword `init!`/`trim!` take ([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)). Not to be
confused with an event *prior* ([§D.4](#d4-time-and-events)).

**boundary zero** — the initialization boundary: the ordinary macro-sequence
with an empty integrate — project → [sweep → guards → handlers]\* → due `g`
updates → header and first snapshot — run at `t₀` once `apply!` has
established the stores ([§14.5](#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions)).

**capture** — the service reading the current committed stores *and* root
slots back as a condition value, returning `(condition, t)`; the gather twin
of `apply!`, and what makes warm restart need no second semantics ([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults),
[§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)).

**component test rig** — a one-child assembly exporting the child's entire
input face set, so any component can be built and simulated in isolation; an
abstract entry is satisfied *inside* the rig by a concrete **stub child**
wired to that face ([§13.7](#137-tooling-consequences-provenance-and-the-component-library)).

**condition** — the datum that says "set this build to this state": a
path-addressed sparse overlay on the declared defaults, covering `x` and `m`
fields plus root slots by face — never outputs, never workspace ([§14.1](#141-conditions-are-path-addressed-overlays-on-the-declared-defaults)).
[§14](#14-stopped-sim-services) owns the word; a guard defines a *predicate* ([§D.4](#d4-time-and-events)).

**`design_world`** — the shipped thin world (aircraft +
`SimpleAtmosphere(wind = NoWind())` + `HorizontalTerrain`) that mounts an
aircraft for trim and linearization; "aircraft as root" is the shallowest
world, not a special case ([§14.9](#149-mounting-problems-as-relocatable-values)).

**fragment** — the leaf node of the condition algebra: `fragment(; x, m,
slots)` payloads speaking only about the component at the authoring point
(**self-vocabulary**), with addressing left entirely to `at` ([§14.2](#142-fragment-composition-locality-without-schema)).

**fragment tree** — the inert, lazy composition of `Fragment`/`Scoped`/
`Merged`/override nodes; isbits but for the interned prefixes, so rebuilding
it per trim iteration is stack-only construction ([§14.2](#142-fragment-composition-locality-without-schema)).

**merge** — the symmetric, collision-intolerant combinator over condition
nodes: a duplicate leaf is an error naming both provenance chains, and mixing
a node with a bare `NamedTuple` is an error method, not `Base.merge` ([§14.2](#142-fragment-composition-locality-without-schema)).

**mounting** — relocating a whole problem or tap set with `at(prefix, …)`:
every field is either condition-producing (path-relative, post-composed) or
path-free, so the service never knows where its paths sit ([§14.9](#149-mounting-problems-as-relocatable-values)).

**override** — the ordered, asymmetric layering combinator: on a shared leaf
the patch wins and provenance keeps both sources, while collisions *within* a
layer remain errors; variadic ([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)).

**service lifecycle** — the `Simulation` states `built` / `initialized` /
`running` / `stopped` / `errored` ([§10.6](#106-run-lifecycle-and-partial-advance)) and each service's legality against
them; a violation is `ServiceLifecycle`, and `errored` is terminal for all
four services ([§14](#14-stopped-sim-services)).

**slot totality** — the pre-write requirement that an application establishing
a complete world over virgin stores — `init!`, trim setup, trim commit —
cover every root slot. Conditions themselves are legitimately partial; a
shortfall is `UninitializedSlots`, collected and declaration-ordered, leaving
the simulation untouched ([§14.6](#146-slot-totality-the-missing-value-error-and-the-override-combinator)).

**taps** — the three selector lists (`x`, `u`, `y`) declaring what
linearization seeds and reports, with an optional component index so a vector
leaf yields named scalars; validated at resolution (`TapResolution`) and
relocatable via `at` ([§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)).

**`TrimProblem`** — the closed seven-field value
`guess`/`lower`/`upper`/`condition`/`reads`/`residuals`/`tolerances`: an
*implicitly specified* condition, solved as a square root-find over named
residuals and committed as an `init!` of `override(baseline, solution)`
([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals), [§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)).

### D.9 Error discipline and diagnostics

**carrier exception** — the single exception a set of diagnostics travels in:
`BuildError` thrown at a stratum barrier, `StepError` at the runtime catch
site. Diagnostics themselves are plain values ([§13.2](#132-diagnostics-structured-values-one-carrier-exception), [§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)).

**collect the checks, fail the evaluations fast** — the reporting policy:
declarative passes over collected structure return their full violation list,
while the first user-code exception aborts the phase; strata are barriers, and
the site column spells the collecting case "build (collected)" ([§13.1](#131-reporting-policy-collect-the-checks-fail-the-evaluations-fast),
[Appendix C](#appendix-c-the-diagnostic-kind-set)).

**did-you-mean** — the required shape of any name-shaped failure: the
offending name plus the list-in-hand it should have come from, carried as
payload rather than baked into message text ([§13.2](#132-diagnostics-structured-values-one-carrier-exception)).

**error locality** — the property the declaration layer buys: a mistake fails
at the site of the mistake, not later and inside correct code. [§11.4](#114-failure-walkthroughs-the-error-locality-grounding)'s five
walkthroughs are its grounding cases and the acceptance tests ([§11.4](#114-failure-walkthroughs-the-error-locality-grounding)).

**execution cursor** — the plain mutable field recording where in the compiled
schedule execution is (component path, function, boundary phase); one cheap
store per dispatch, so a runtime failure gets its frame without exception
frames in the hot path ([§13.4](#134-runtime-failures-one-catch-site-an-execution-cursor)).

**feedthrough tracer** — the set-propagation instrument (global value-blind,
or local primal-carrying at sampled states) used to classify a rejected cycle
as real or artificial; diagnostic only, never an input to scheduling ([§5.6](#56-diagnostics-feedthrough-tracing)).

**kind** — a diagnostic's identity in the closed set enumerated normatively in
[Appendix C](#appendix-c-the-diagnostic-kind-set), with payload fields, owning section and severity; tests match on
kind plus payload, never on message text ([§13.2](#132-diagnostics-structured-values-one-carrier-exception)). Not a component *class*
([§D.1](#d1-component-model-and-declaration-layer)) or a *function family* ([§D.1](#d1-component-model-and-declaration-layer)).

**payload** — the structured data a diagnostic carries beside its kind: paths
and names as strings (never instances or model types), expected/observed port
types, the list-in-hand, the severity ([§13.2](#132-diagnostics-structured-values-one-carrier-exception), [Appendix C](#appendix-c-the-diagnostic-kind-set)).

**poison** — the debug-mode overwrite of a workspace before every call (`NaN`
for float-eltype stores, `typemin` for integer ones) so read-before-write of
stale scratch detonates immediately; stores with no sentinel are skipped and
named once per activation (`PoisonSkip`) ([§7.3](#73-discrete-state-modes-and-workspace)).

**`stop_on` / termination is a state** — graceful termination is model state,
never an exception: detection is ordinary event machinery, publication an
ordinary root-exported `Bool` output face, and `stop_on` the deployment policy
naming the faces the loop reads after every published boundary ([§13.5](#135-termination-is-a-state-not-an-exception)).

**warning streams** — two, scoped separately: the *build* stream, whose
warning set is deliberately empty, and the *runtime* stream — per-occurrence,
carried by the per-writer diagnostic cells that structurally rate-limit it
([§9.8](#98-diagnostics-and-liveness-the-per-writer-cell)), surfaced through published
framework status, with its committed inventory listed in [§13.2](#132-diagnostics-structured-values-one-carrier-exception).

### D.10 Meta-vocabulary

**blessed** — the spec's marker for a practice it explicitly sanctions where a
neighboring one is forbidden: derivation from other declarations ([§11.2](#112-the-declaration-inventory)), the
one spot where evaluation feeds structure ([§12.1](#121-three-strata)), the workspace-plus-snapshot
idiom for zero-allocation ticks ([§7.3](#73-discrete-state-modes-and-workspace)).

**the freeze** — the roster freeze: `attach!`/`detach!` are stopped-sim
operations, so the roster, its claims and the run's partition of the root face
set into write surfaces are static, inspectable facts of each run ([§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster), rows
106–107).

**guarded addition** — a capability the design admits but does not build,
weighed against Flight.jl's fundamental strengths and recorded with its shape
so adoption stays additive ([§1](#1-purpose-and-method); e.g. [§4.3](#43-table-mechanics-and-port-granularity)'s field-addressed staging, [§9.3](#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s
mid-run reader attach).

**normative / index, not a second home** — the spec is the normative statement
of the design, and its appendices are indices: each recall line's normative
statement stays in the owning section ([Appendix A](#appendix-a-taught-contracts-the-author-facing-index), and the same rule for
Appendices B, C and D). The design directory's walkthrough explainers are
non-normative companions by their own preambles.

**recorded, not built** — the disposition of a worked-out extension
deliberately left unimplemented, with its seams named so adoption is additive
([§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals)'s closed-loop trim, [§14.10](#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query)'s sampled-data `Dual` activation and
declarative non-participation).

**register** — the spec's word for a mode or idiom in which something is done,
always compounded: the didactic register ([§13.2](#132-diagnostics-structured-values-one-carrier-exception)), the inspection and
integration registers ([§9.2](#92-outbound-snapshot-publication)), the by-allocation register ([§11.2](#112-the-declaration-inventory)), the
harness, unattended and what-if registers ([§10.6](#106-run-lifecycle-and-partial-advance), [§10.7](#107-replay-the-trace-re-drives-the-ordinary-loop)). Reserved for this
sense — the recording artifacts are the *recorders* ([§D.7](#d7-recording-and-replay)).

**row** — a numbered entry of `framework_decisions.md`, cited throughout as
"row N": one settled decision with the alternatives weighed against it. Row
numbers are stable and never reused, and each row states its *current*
position, a superseded one being demoted into its rejected-alternatives column
rather than rewritten ([§1](#1-purpose-and-method)).

**seam** — a narrow, named interface kept deliberately thin so what sits
behind it can be replaced or measured: the stepper seam ([§8.2](#82-the-stepper-seam)), the backend
seam ([§14.8](#148-the-trim-service-solver-seam-scratch-stores-commit-and-report)), the measurement seam ([§12.7](#127-the-compiled-executor)), the phase-body seams of the
compiled executor ([§12.7](#127-the-compiled-executor)).

**torture test** — an existing, maximally awkward artifact transliterated
against a proposed mechanism to validate it before adoption: `PistonEngine`
and the FCS cascade against [§5.2](#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws) ([§15.2](#152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade)), filter/joystick/GUI against [§9](#9-runtime-periphery-the-data-plane)'s
staging shapes ([§15.3](#153-torture-test-for-the-9-staging-shapes-filter-joystick-and-gui)), the strapdown IMU against the leaf split ([§15.5](#155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary)); the
standard component library is the standing ergonomics one ([§13.7](#137-tooling-consequences-provenance-and-the-component-library)).

**worked (example)** — a full spelling of a mechanism against a real artifact,
carried in the spec rather than left to the reader: the worked assembly of
[§11.6](#116-paths-wiring-and-faces), the worked C172 cruise problem of [§14.7](#147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals), and [§15.5](#155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary)'s IMU as the
boundary-sampling example [Appendix A](#appendix-a-taught-contracts-the-author-facing-index) points at.
