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
    - [8.6 Event iteration at boundaries: to quiescence, budgeted](#86-event-iteration-at-boundaries-to-quiescence-budgeted)
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
It is the [normative statement](#g-normative) of the design: what the framework *is*, in
present tense. The new framework must match or surpass `FlightCore` in
functionality, performance and flexibility, while being more rigorous and
explicit — reducing the learning curve and the number of latent footguns for
model authors.

Ground rules adopted for this design:

- **Capability grounding, not interface grounding.** Requirements are derived from what
  `FlightPhysics` and `FlightApps` demonstrably *do* (code, unit tests, demos). Every
  `FlightCore` call site in the consumers is read as evidence of a capability the
  substrate must provide, never as a prescription for how it should be spelled.
- **No interface compatibility.** The new framework need not be source-compatible with
  the current consumers. A non-trivial migration of `FlightPhysics` and `FlightApps` is
  expected and accepted.
- **[Guarded additions](#g-guarded-addition).** Whenever the design admits functionality beyond what the
  consumers demonstrate, it must be weighed against the fundamental strengths of
  Flight.jl: zero-allocation stepping, type stability, real-time interactive operation,
  live introspection (GUI), and compositional flexibility.

All design axes are settled — the formalism, the [component](#g-component) taxonomy, the
signal and scheduling model, time and execution, the runtime [periphery](#g-periphery), the
declaration layer, the build pipeline, error discipline and the stopped-sim
services. Only the items ([§16][s16]) — the migration outline, the GUI panel authoring
API and the log/[trace](#g-trace) persistence deferral — remain open.

Decision rationale, including the alternatives considered and the reasons
they were rejected, lives in `framework_decisions.md`, cited throughout as
"row N": one row per settled decision. Row numbers are stable, so a citation
here always names the same row there.

---

## 2. Formalism

The framework simulates **[hybrid causal systems](#g-hybrid-causal-system)**, composed of:

- **Continuous dynamics**: $\dot{x} = f(x, m, u, t)$ with algebraic outputs.
- **Multi-rate periodic discrete dynamics**: $x^{+} = g(x, u, t)$ at declared rates, with
  outputs held zero-order between [ticks](#g-tick).
- **Zero-crossing events**: [guard](#g-guard) functions with handlers, under two detection policies (below).
- **Post-step manifold [projection](#g-projection)**: an optional per-[component](#g-component) hook `x ← project(x)`
  applied after each accepted step (quaternion renormalization, DCM orthonormalization,
  any manifold-valued state). This is the cheap end of the projection-methods family
  from geometric integration.
- **External inputs**: injected asynchronously by the runtime (pilot controls, network),
  under the staging rules settled in [§9][s9].

### 2.1 Events: two detection policies

Both policies share one declaration ([guard](#g-guard) function + handler); only detection
differs:

- **[Boundary-detected](#g-boundary-detected) (cheap):** guards are checked for not-holding → [holding](#g-edge-semantics) edges
  against their [priors](#g-prior) ([§8.6][s8-6]) at step [boundaries](#g-boundary)
  only. No root-finding, no step rejection; the handler fires at the end of the step in
  which the edge was observed. Cost: one guard evaluation per event per step. Fully
  compatible with fixed-step real-time execution.
- **[Localized](#g-localized) (precise):** localization of the crossing instant by root-finding,
  for events where timing precision genuinely matters (mechanics in [§8.4][s8-4]).
  Detection policy never depends on real-time [pacing](#g-pacing) ([§8.7][s8-7]).

This arrangement gives step-boundary logic *well-defined semantics*: the transition is
defined by the crossing; detection resolution is an execution-policy detail.

A guard defines a **[predicate](#g-predicate)**: either a `Bool`-valued form, or the
sign of a continuous function under the normative convention **positive = predicate
holds**. Writing the guard's sign value `σ`, holding = `σ ≥ 0`.

An event fires when its predicate transitions from not-holding to holding —
[edge semantics](#g-edge-semantics), uniform across both forms. The prior bookkeeping
is stated in [§8.6][s8-6]. The opposite crossing direction is declared as a second
event with the negated guard (stall entry/exit as a pair).

Which form an author declares is not a free choice: **the guard's return type is the
declared policy**. A `Bool` return declares boundary-detected; the sign form declares
localized (row 179). [§8.4][s8-4] states the rule, the exactness result that
motivates the `Bool` form, and the gate idiom for localizing mixed
predicates.

### 2.2 Exclusions (deliberate)

- **No DAEs / algebraic constraints.** [Projection](#g-projection) is what covers the
  actual need: state manifolds (row 1).
- **No SDEs / stochastic integrators.** Noise processes such as Dryden/von Kármán
  turbulence and sensor noise are modeled as ordinary RNG-driven discrete processes
  (shaping filters). That modeling is both faithful to how they are specified and cheap.
  One consequence is elevated to a framework guarantee: **deterministic
  [replay](#g-replay)**. RNG state lives in [component](#g-component) discrete state
  (`x`), never in ambient globals, so same seed ⇒ bit-identical trajectory.
- **No unconditional per-step hook** (no `f_step!` equivalent). Every current use
  decomposes into one of two mechanisms. Projection covers quaternion renorm;
  [boundary-detected](#g-boundary-detected) events (checked for edges at step boundaries
  only, no root-finding) cover engine phase transitions and the stall hysteresis latch.
  For one class the mapping tightens semantics: level-triggered cross-component resets
  become edge-triggered events. The gear friction regulator under `!wow` is one such
  reset ([§15.2][s15-2], [§16][s16], row 1).

---

## 3. Component taxonomy

Three classes, two of them leaves with crisp, closed semantics, one pure composition:

### 3.1 Continuous component (the hybrid primitive)

A classical hybrid automaton:

- **continuous state** `x` (isbits struct of real scalars — see [§7][s7]),
- **mode variables** `m`: piecewise-constant values (enums, integers, flags) that
  parametrize the [flow](#g-flow) and change *only* through event handlers,
- **flow** $\dot{x} = f(x, m, u, t)$,
- **two output stages** (see [§5.2][s5-2]),
- **events**: [guards](#g-guard) + handlers (update `m`, may reset own `x`); both read the fresh
  [boundary](#g-boundary) [signal table](#g-signal-table) ([§5.3][s5-3]),
- optional **[projection](#g-projection)**.

Any facet may be empty. In particular, a [component](#g-component) with *no* continuous state — only
modes, events and mode-valued outputs — is an FSM. Both factorings of mode logic are
therefore supported with a single primitive:

- **internal modes** for tightly-coupled cases (stall hysteresis inside the aero
  component) — preserves cohesion and enables reset maps;
- **external FSM component** feeding modes through a [port](#g-port) — maximal purity, independent
  testability, swappable supervision logic.

Rule of meaning: a supervisor *commanding* a mode change is an ordinary **input**; a
component *detecting* its own transition is an **event**. Two mechanisms, two meanings,
no overlap.

### 3.2 Periodic discrete component

- **discrete state** `x`: any immutable value (see [§7][s7]),
- **update** $x^{+} = g(x, u, t)$ at a declared rate,
- **two output stages**, with [feedthrough](#g-feedthrough) applying at update instants: a
  proportional path is direct feedthrough; a state-only output is not.

State carries the same letter on both [tiers](#g-tier): `x` is the argument of the
[flow](#g-flow) map under `f` and of the jump map under `g`. One letter for both maps
is standard hybrid notation (row 173). The letter is never ambiguous, because
a leaf is strictly one tier (row 56) and no [component](#g-component) ever reads
another's state. A discrete component's `x` influences continuous dynamics only
**through signals**. Those outputs are held zero-order between [ticks](#g-tick).

**`m` is continuous-only.** A discrete component has no mode store: its FSM
enums, flags and counters are ordinary `x` fields.

**Why.** `m` exists on the continuous side because modes must change *between*
flow evaluations, which is what handlers do. On the discrete side `g` already
runs at the only instants at which anything may change, so a second store would
duplicate the discrete tier's own state semantics under another name.

### 3.3 Assembly

Pure composition: submodels + child connections + [boundary](#g-boundary) [faces](#g-face). **No dynamics of its own.**
Hybridness emerges at the [assembly](#g-assembly) level (an aircraft = continuous vehicle parts +
discrete avionics parts). The two-leaf split was upheld against the
integrate-and-dump challenge ([§15.5][s15-5], row 56). Assemblies are flattened away for scheduling but retained as
the navigation/introspection hierarchy (GUI, logging, paths) and as declaration-level
[rate scopes](#g-rate-scope) ([§8.5][s8-5]).

---

## 4. Ports and signals

### 4.1 Immutable value semantics

[Ports](#g-port) exchange **immutable values** — typically isbits structs (floats, `SVector`s,
enums, nested immutables). The framework owns a **[signal table](#g-signal-table)**: one
concretely-typed **[cell](#g-cell)** per output port in the flattened model. A producer's
output-[stage function](#g-stage-function) (`h_x` or `h_xu`, the two output stages every component
provides) returns a named tuple of fresh values. The framework writes each of those
values into its cell, and consumers read cells.

**Vocabulary.** These names are binding throughout this document:

- bare *cell* — the table entry, and only that;
- *store* — the discrete-state and mode registers ([§7.3][s7-3]), which are not cells;
- *[staging cell](#g-staging-cell)* — a distinct compound term, the per-[device](#g-device)
  inbound register of [§9.4][s9-4]; unlike a table cell it is mutated frame by frame and
  sits outside the table's publish-once discipline;
- *[slot](#g-slot)* — reserved for the root input slots of [§9.3][s9-3].

The signal requirement, stated precisely, is **immutability plus frozen references**:
signals may reference bulk data (see [§4.4][s4-4]) provided that data is read-only for the
duration of the run. `isbits` is the common case, not the rule.

Consequences:

- no aliasing, ever — nothing can be mutated under a consumer's feet;
- safe concurrent reads (GUI/logging threads) by construction;
- zero allocation for isbits payloads (named tuples of isbits are isbits);
- each cell has a definite freshness tied to its producer's position in the [schedule](#g-schedule)
  (row 4).

### 4.2 Consumers see ports, not stages

The [port](#g-port) is the addressable unit. A [component](#g-component)'s outputs appear to consumers, GUI and
logs as one flat namespace (`dyn.vel`, `dyn.f_c_c`), materializable lazily as a view.
Which output stage computes which port is a scheduling annotation, invisible outside
the component. Moving an output between stages is non-breaking *for consumers*: no
wire, log or panel sees it. The scheduler does see it — the [feedthrough](#g-feedthrough) graph and
stage membership change ([§12.1][s12-1]).

**Visibility.** Which ports exist at all is a declaration-layer decision: the output
[contract](#g-contract) *is* the public interface. Stage-function results outside that contract
never enter the table at all. They ride the stage's `w` return down to the component's
own later functions, private by construction ([§5.2][s5-2], [§11.3][s11-3]). A presentational
*unlisted* flag — skipped in logs and GUI but still connectable — is closed (row 16).

### 4.3 Table mechanics and port granularity

[§4.2][s4-2] fixed the [port](#g-port) as the addressable unit. Four questions remain:
what a port is on a [component](#g-component)'s [boundary](#g-boundary), how values travel into cells and
out of them, what a port may hold, and how much a single port should carry. The
last question has two answers, one owed to the parties that read a port and one
to the parties that write it.

#### Ports and faces

A component's **ports** are its signal endpoints — one [cell](#g-cell), one producer. Its
**[faces](#g-face)** are the names those ports wear on the component's boundary. For a leaf
the two coincide. For an [assembly](#g-assembly) every face aliases an interior port through
its boundary declarations ([§11.6][s11-6]) and never creates an endpoint. The
distinction is kind-blind — wiring and the [periphery](#g-periphery) address a child's faces
without knowing whether it is primitive or composite.

**Rule.** The port is the atomic unit of the entire periphery — one cell, one
root [slot](#g-slot), one staged write, one [device](#g-device) [claim](#g-claim) ([§9.3][s9-3]), one [trace](#g-trace)
address, one GUI liveness verdict ([§9.7][s9-7]).

#### Scatter and gather

**Scatter/gather is the whole protocol.** A [stage function](#g-stage-function) (`h_x` or `h_xu`,
the two output stages every component provides) returns a named tuple. The
framework scatters each field into that port's concretely-typed cell. Every
reader — the next stage, `f`/`g`, [guards](#g-guard), wired consumers, [snapshot](#g-snapshot)
capture — gathers views from cells.

**The aggregate `y` is a merge semantically and virtual physically.**
Semantically, a component's `y` is the merge of its stage products
(`merge(y_x, y_xu)`, the same on either [tier](#g-tier)). It carries declared ports only;
a stage's private intermediates ride its `w` return rather than the table
([§5.2][s5-2]). Physically no such object exists: `y` is reconstructed per call from
cells — field loads, register-level, zero cost for isbits — and never stored as
an object.

Name collisions across a component's stages are a build error.

#### What a port may hold

**Stage returns are named tuples of port values, period.** A custom struct is a
first-class port *value* — one field of the returned tuple, one declared port, one
cell (`pose = KinPose{T}`). Nested fields get no cells of their own; GUI and logs
drill into them lazily (the view clause, [§4.2][s4-2]). Bare-struct returns are rejected
(row 36).

#### Granularity, read side

**Rule.** Wiring is port-granular: no sub-field connections. A consumer that
wants less than a bundle asks the producer for a loose port, or takes the bundle
and destructures. A field-projection connector is a [guarded addition](#g-guarded-addition) (a
capability the design admits but does not build). Its shape is obvious, and it
is not built.

**Granularity guideline** for authors: bundle what *shares a stage* *and is
consumed together*. The first criterion is trivially enforced, because each port
has exactly one producing function. Bundling across dependency footprints is the
`KinData` mistake ([§15.1][s15-1]). Pose is stage 1, velocity-derived quantities are
stage 2 — it must split. Fan-out is free, so publishing both a bundle and a hot
loose field (`pose` *and* `q_eb`) is legitimate — one extra isbits cell.

**Example.** The bundle, the loose hot field and a face aliasing an interior
port, in declaration form:

```julia
#a continuous leaf: one bundle port and the hot field published loose — two cells
output_types(::Kinematics, ::Type{T}) where {T <: Real} =
    (pose = KinPose{T}, q_eb = RQuat{T})

#the enclosing assembly: each face aliases an interior port, creating no endpoint
output_connections(::Vehicle) = ("kin/pose" => "pose", "kin/q_eb" => "q_eb")
```

#### Granularity, write side

**Write-side corollary** (from [§15.4][s15-4]): **bundle what is written together.**

**Rule.** Data written by different external writers, or at different cadences,
must not share a port.

Pilot commands are the case in point: scalar faces under a namespace prefix,
with the convenient bundle assembled *downstream*, inside the graph, by an
ordinary component (single producer, consumed together — legal by the read-side
rule).

The two guidelines compose into one principle: a port's granularity is set by
the finest-grained party owning either end — producers on the read side,
external writers on the write side. Field-addressed staging (a [lens](#g-lens) into
struct slots) stays a recorded guarded addition, unbuilt.

### 4.4 Function-valued signals: environment access

Atmosphere and terrain are **query-shaped**: consumers evaluate them at arguments of
their own choosing (each gear strut at its own contact point; airflow at the vehicle
pose). They are therefore carried by ordinary [ports](#g-port) as **immutable query objects**
("[field handles](#g-field-handle)"):

- An environment [component](#g-component) emits a field value (`ISAField(T_sl, p_sl,
  wind)`, `TerrainField(...)`); consumers receive it through ordinary input ports.
  Inside their own [stage functions](#g-stage-function) (`h_x` or `h_xu`, the two
  output stages every component provides) they call query functions on it:
  `airdata(field, pos, vel)`, `ray_intersect(field, p, u)`.
- **Parametric models are isbits** (ISA, uniform wind, horizontal terrain). **Bulk-data
  models use the handle pattern**: an immutable struct combining isbits parameters with
  references to bulk data (heightmaps, wind grids, the geoid undulation grid) loaded at
  build time and frozen. Handles are rebuilt per evaluation allocation-free — immutable
  structs with existing references. Never `Ref`s, whose mutable cell allocates.
- **No mutable caches inside field objects** (memoizing interpolators, lazy loaders):
  concurrent consumers and the GUI thread would race. Caches belong in the consumer's
  state, or the interpolant is restructured to be pure.
- Loggers treat field-handle signals specially (skip or summarize).

**The [value-level constructor](#g-value-level-constructor).** Every field-emitting
component must expose the map (component, input values) → handle as a plain,
pure, exported function — `atmospheric_field(atm; T_sl, p_sl, wind)` for the
`SimpleAtmosphere` successor. The field-emitting component's swept output stage
must be a **one-line call to that function**, never the other way round. The
other way round puts the query math in the output stage, where only a
[sweep](#g-sweep) can reach it.

The reason is script-side: the condition math ([§14.1][s14-1]) must be able to
construct, outside any sweep, bit-for-bit the same handle the sweep would
produce from the same [slot](#g-slot) values. One implementation, two call
sites, no drift — and the drift avoided here is the silent-drift class that
[§5.3][s5-3] exists to kill. This is a *shipped component's obligation*, not
something a consumer can retrofit. The real component composes sub-models, and
anyone else reconstructing the map has re-created the drift class.

For bulk-data components the obligation is only that the query math be
reachable as a plain function. They own their resource loading, so building a
handle outside a build may cost a load. That cost is acceptable, because
condition authoring is design-time code.

**Example.** The map as a plain function, and the stage that does nothing but
call it:

```julia
#the map: plain, pure, exported — callable outside any sweep
atmospheric_field(atm; T_sl, p_sl, wind) = ISAField(…)

#the swept output stage: one line, nothing but the call
h_xu(atm, args) = (; … = atmospheric_field(atm; T_sl = …, p_sl = …, wind = …))
```

Pre-sampling — a component consuming the field and a pose and emitting plain data
(`Airflow` emitting `AirData` for the whole vehicle) — is an **idiom built on top**,
used where natural; not a separate mechanism. Resource injection (declare-and-resolve
service registries) is closed for the first cut (row 8).

The field-handle mechanism replaces threading `atmosphere`/`terrain` as arguments
through every update signature, and dovetails with the terrain ray-query direction
of the landing-gear redesign. Substitutability behind a stable [face](#g-face) is
declared with an abstract input entry — `terrain = AbstractTerrainField`, structural
substitutability ([§11.2][s11-2]). The consumer wires to any concrete field type
below the bound, preserving today's `AbstractTerrain` polymorphism at the
declaration layer.

---

## 5. Evaluation order and feedthrough

### 5.1 The scheduling problem

At every evaluation instant, all signals must be computed consistently: every consumer
reads values already produced at that instant. Build the directed graph of (a) wiring
edges and (b) intra-[component](#g-component) [feedthrough](#g-feedthrough) relations; if acyclic, a topological sort
yields a **static evaluation [schedule](#g-schedule)**, computed once at build time. The hot loop runs
a flat list of `(component, stage)` entries — zero runtime graph logic.

### 5.2 Two-stage outputs: signatures, bundles and the hand-off laws

Every [component](#g-component) provides exactly **two output stages**, and [feedthrough](#g-feedthrough) is declared
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

The laws that govern that surface follow: how a function receives its
arguments, which names it receives, what an output stage returns, how far a
private intermediate travels, and what a handler returns.

#### The hand-off: one component, one bundle

**Rule.** Every function receives exactly two arguments: the component and one
NamedTuple [bundle](#g-bundle) of zero-copy views. From that bundle the author
**destructures by name** only what the body reads:
`f(c::LowPassFilter, (; x, u)) = ...`, `h_xu(c::PID, (; x, u, Δt)) = ...`.

**Why.** The [executor](#g-executor) (the compiled execution form of the
schedule) issues one fixed call shape, `fn(comp, args)`. Unread fields are
ignored by language semantics. Argument order cannot be confused, because there
is no order.

Positional, `kwarg_decl`-reflected and slurping-keyword spellings are all
closed (row 74). `project` alone stays positional: one store in, the same store
out, nothing to select.

#### The bundle law: which names a component receives

**Rule.** Under the [bundle law](#g-bundle), a name appears in a component's
bundle **iff the corresponding store or fact exists for that component**.

| bundle field | present iff |
|---|---|
| `x` | the component declares `init_x` |
| `m` | the component declares `init_m` |
| `ws` | the component declares `workspace` |
| `u` | the [function family](#g-function-family) (which bundle fields a given function may legally receive) may see inputs **and** the component declares `input_types` |
| `y` | the component produces any table [cell](#g-cell) at all (`output_types` ∪ auto-published) |
| `y_x` | stage-1 [ports](#g-port) exist |
| `w` | the stage that hands it down returned one (the one-hop law below) |
| `t` | always |
| `Δt` | the component is on the discrete [tier](#g-tier) |

**`y_x` carries the stage-1 *return*, auto-published names excluded.** An
[auto-published port](#g-auto-published-port) is the framework copying a state
or mode field into a cell at stage-1 position ([§5.3][s5-3]), and stage 2
already holds `x`/`m` (continuous) or `x` (discrete) directly. The rule is what
[§12.3][s12-3] already sources: `y_x` comes from the stage-1
[probe](#g-probe)'s return. So a component whose only stage-1 ports are
auto-published has no `y_x` in its stage-2 bundle at all (row 169).

Undeclared stores are *absent*, never `nothing`-filled. Destructuring a field
that is not a thing for you fails at the probe inside the [§13.2][s13-2]
framing diagnostic, with [did-you-mean](#g-did-you-mean) (the offending name
plus the list-in-hand it should have matched) against the legal field set:
"`f` of `Foo` destructures `m`, but `Foo` declares no `init_m`". One law covers
tier facts, stage legality and declarations alike.

The mechanism is structured, not textual. Destructuring an absent field throws
a `FieldError` carrying the type and the field name as data (Julia ≥ 1.12). The
probe catches it *matched against the bundle's own NamedTuple type*, and
synthesizes the framing diagnostic from the legal set — classifying the field
as an undeclared store, a wrong-tier fact, or a name illegal for this function
family. No message text is scraped, and the bundle stays a bare NamedTuple
(row 74); a getproperty-wrapper spelling is the recorded fallback should
type-matched interception prove insufficient.

The wrong-tier class holds under the shared state letter (row 173): state
itself is legal on both tiers, `m` is continuous-only and `Δt` discrete-only,
so destructuring either on the wrong tier still lands in that bucket.

The per-function, **per-tier** name sets are **closed**: adding one is a
decision-log entry, not a convenience. The comment in the signature block above
states each function's maximal legal set; a given component's bundle narrows it
to declared reality; the destructuring narrows further to actual reads. That
three-level funnel (stage name ⊇ bundle ⊇ reads) is worth teaching once,
because a stateless component legitimately writes `h_xu` while owning neither
`x` nor `m`.

#### The stage return law

**Rule.** An output stage returns either its port NamedTuple alone, `y`, or the
pair `(y, w)`.

`y` scatters into the component's declared cells as always. `w` is a NamedTuple
of **private intermediates** — `isbits` leaves, no cell, no name in any
[contract](#g-contract), nothing to wire, list or filter ([§11.3][s11-3]). A
`nothing` in either slot is a probe error: the pair is a
`Tuple{NamedTuple, NamedTuple}` and padding forms do not exist, for the reason
the handler return law gives below.

An empty `y = (;)` *is* legal, so the **port-less stage** — one whose entire
product is `w` — falls out of the general law instead of needing a rule of its
own. Stages are discovered by method existence, and stage membership is a
partition of the declared ports that may perfectly well be empty. What is not
legal is a stage that produces neither ports nor `w`: a bare `(;)` computes
nothing any consumer can read, and is a `DeadStage` build error at the probe
([§12.3][s12-3]) — the inert-component check in the stage register
([§11.1][s11-1]).

#### The one-hop law

**Rule.** `w` travels exactly one hop, to the next function that could want it,
and no further.

`h_x`'s `w` flows to `h_xu` if the component defines one, and otherwise to the
downstream set — `f`, [guards](#g-guard) and handlers. `h_xu`'s `w` flows to
the downstream set. The discrete tier mirrors it exactly (`h_x` → `h_xu` →
`g`), and that is the last time it is said.

**Nothing flows implicitly past its hop.** Forwarding a stage-1 intermediate
through stage 2 is an explicit re-return, `(y, (; w..., extra))`.

**Why.** The re-return costs a line and says in the source that the value
crosses; a silent pass-through on a bare `y` was rejected (row 165).

```julia
# one hop, on a component defining both stages
h_x(c::Foo,  (; x))       = (y, w)   # w is produced here
h_xu(c::Foo, (; x, u, w)) = y        # it arrives here, and stops
f(c::Foo,    (; x, u))    = ẋ        # no w key: h_xu returned a bare y

# the same stage, forwarding across the hop — f, guards and handlers see it now
h_xu(c::Foo, (; x, u, w)) = (y, (; w..., extra))
```

Presence in a bundle follows the producing stage's return under the bundle law:
a bare-`y` producer hands down no `w` key at all, so a consumer destructuring
one meets the [§13.2][s13-2] framing diagnostic naming its legal field set,
exactly as for an undeclared store.

`w` is never persisted. The executor hands it down as an ordinary SSA value
**inside one fused pass** — a step fuses the [sweep](#g-sweep) with `f`, an
event round fuses the sweep with its guards and fired handlers ([§8.6][s8-6]).
So freshness is a property of the construction rather than a rule anyone can
violate, and round fusion is thereby a design constraint on the executor, not
an optimization it may decline ([§12.7][s12-7]).

`w`'s types are probe-observed, and [§11.3][s11-3] states the conformance
regime that governs them. The handler return law below is untouched by any of
this — handlers *receive* `w` and return stores.

#### The handler return law

**Rule.** A handler returns a NamedTuple carrying the stores it writes: a key
is present **iff** the corresponding store exists on the component **and** the
handler updates it.

That is the bundle law's *iff* shape, now governing the return side. A pure FSM
(modes and events, no `x`) returns `(; m = (; phase = running))`; an `x`-only
reset map returns `(; x = (; x..., ω = 0.0))`; a handler touching both returns
both. Padding forms — `((;), m⁺)`, `(x⁺, (;))` — do not exist
(row 90, on row 74's argument-side ground).

Semantics per key: `x` present ⇒ the value is complete against the state field
set; `m` present ⇒ the names-subset [predicate](#g-predicate); an unknown key ⇒
did-you-mean against `{x, m}` — the same `FieldError`-shaped machinery
[§13.2][s13-2] builds for bundles, now running in both directions
([§12.5][s12-5]).

#### What the views are

The views themselves are unchanged in meaning:

- own state — on the continuous tier `x` from the flat [buffer](#g-buffer) and
  `m` from the mode [stores](#g-store); on the discrete tier `x` from its store;
- own published signals — `y`, gathered from own table cells (the declared
  ports, [§11.2][s11-2]);
- own private intermediates — `w`, handed down by the producing stage rather
  than gathered from anywhere, the one bundle field with no home at all;
- inputs — `u`, gathered from foreign cells through the wiring's name binding;
- the clock — `t`, and `Δt` (see [§8.5][s8-5]);
- scratch — `ws` ([§7.3][s7-3]).

The [signal table](#g-signal-table) holds only *produced* signals, never
transported ones: each datum has exactly one home — buffer for continuous `x`,
stores for discrete `x` and for `m`, table for signals — and no store mirrors
another. Every bundle field earns its place as a view genuinely readable, and
no minimization of the set survives without introducing a copy (row 35).

### 5.3 Structural feedthrough: stage roles, schedule and step boundaries

[§5.2][s5-2] fixes the two-stage surface and the laws that govern it. What remains is
the reading. A stage's role has two halves: what its name asserts, and what it may see.
The rest of the section orders the stages into a schedule, then puts that schedule
inside a step [boundary](#g-boundary).

#### Stage roles: the letters

**The letters**: `f` is the continuous [flow](#g-flow), `g` the discrete update, `h_*`
the output stages. The scheme joins two traditions: `f` and `g` are the hybrid-systems
flow/jump pair (Goebel–Sanfelice–Teel), and `h` is the control/estimation convention
that the output map is `y = h(x, u)` — every navigation filter's measurement function.
Bare `h` denotes the integration step size only ([§8][s8]).

**Rule.** Stage suffixes name the **dependence class**, not the argument list. `x`
versus `xu` is state-only versus state-plus-input — the `y = h(x)` / `y = h(x, u)`
distinction spelled in the name, identically on both [tiers](#g-tier). So "no `u` in the
name" *is* the no-[feedthrough](#g-feedthrough) property, visible at every definition
site.

The letters are deliberately non-exhaustive. Modes fold under the state letter: `m` is
state, and the suffix names the [feedthrough](#g-feedthrough) split rather than an
argument inventory (row 75). Ambient facts (`t`, `Δt`) and scratch (`ws`) ride unnamed.

The stage names do not distinguish the tiers at all — the state letter is shared
(row 173). So what declares a stateful leaf's tier is its update law, `f` versus `g`,
with the remaining tier-implying declarations agreeing ([§11.2][s11-2]).

#### Stage roles: what each stage may see

**`h_x` is the no-feedthrough stage**, defined entirely by what it cannot see. Its
[bundle](#g-bundle) (the NamedTuple of zero-copy views a component function receives)
carries no `u`, so "no feedthrough" cannot be violated by construction. That structural
guarantee is what `h_x` [ports](#g-port) contribute to the [schedule](#g-schedule): they
break would-be loops.

`h_x` exists when the [component](#g-component) has state-derived ports, or shared
state-derived intermediates to hand down its `w` return ([§5.2][s5-2]); otherwise it is
simply absent. A stage that would produce neither is the `DeadStage` error, an empty
stage being unwritable on purpose.

Guidance rather than law: when the component also defines `h_xu`, a `w`-only `h_x`
earns nothing. Fold the intermediates into `h_xu`, which runs exactly once per
[sweep](#g-sweep) just as `h_x` does. The port-less `h_x` earns its keep where there is
no `h_xu` at all, its no-`u` bundle being the honest spelling of "these intermediates do
not depend on inputs".

**Rule.** A declared output that matches a state or mode field by name and type, and
that no stage produces, is auto-published by the framework from the state stores at
stage-1 position ([§11.3][s11-3]). The match is against the declared stores — `init_x`,
plus `init_m` on the continuous tier — and the publication position is `h_x` on either
tier. Publication is driven by the public [contract](#g-contract) (row 16).

**`h_xu` receives all wired inputs plus `y_x`** — its own stage-1 ports. With them it
receives stage 1's `w`, so shared intermediates are computed once, not re-derived,
whether or not they are interface. It receives the [state views](#g-view) too.
Conservatively, every stage-2 output is presumed dependent on every wired input.

**`f` and `g` run after the sweep**, when the full [signal table](#g-signal-table) is
complete and fresh — the component's own stage-2 ports included. The fused idiom stands:
compute each law once, in a stage; publish it; let `f`/`g` copy from `y`.

**Why.** The interfaces *reward* single-source-of-truth rather than making duplication
unwritable: nothing ever needs computing twice (rows 15, 35).

All output stages must be pure (no side effects); state types make mutation impossible
anyway ([§7][s7]).

#### The schedule

**Rule.** The schedule runs all stage-1 functions in any order, then stage 2 in
topological order, then all `f` against the now-consistent signal table.

Note the systemic consequence: *evaluating the [RHS](#g-flow) means running the sweep*.
There is no incremental `f`-only re-evaluation. Nothing is lost by that, because
implicit solvers, linearization and trim already work this way — seed `x`, run the
composite. [§8.3][s8-3]/[§8.4][s8-4] restate that consequence as a property of the
execution model. RHS evaluations and guard [probes](#g-probe) alike run the *interior*
sweep, the continuous-only variant of [§8.5][s8-5]; discrete entries are absent from it
by construction, so discrete [cells](#g-cell) hold across the step.

#### Step boundaries

**[Guards](#g-guard) and handlers read the same fresh world.** At a step boundary the
order is *integrate → project → [boundary sweep](#g-sweep) → guards*, so by
guard/handler time `y` is a fresh decode of exactly the state being transformed, and the
state views are that state itself. Handlers construct their `x`/`m` returns from raw
state naturally: a reset map is `(; x = (; x..., ω = 0.0))`, no reassembly from published
fields.

**`project` runs between a state write and its decode** — after integration, and after
any handler `x`-reset. Those are the only positions in the schedule where no fresh `y`
of the new state can exist yet. `project` is not *unique* in receiving raw state (every
function gets state views), but it is unique in that schedule position.

**The boundary sequence.** At each boundary: integrate → project → boundary sweep →
evaluate **all guards once** against that sweep → fire the eligible events, each firing
being `handler → project`. The sweep → guards → handlers phase then iterates to
[quiescence](#g-quiescence) (the fixed point where a round of handlers fires nothing),
so newly-enabled guards fire within the *same* boundary:

> integrate → project → **[ sweep → guards → handlers ]** iterated to quiescence

The signal table is written **only by sweeps**, so a transition reaches the table at the
next round's re-sweep ([§8.6][s8-6]). The round that detects quiescence leaves the table
post-transition-consistent for whatever else the boundary does — discrete
[ticks](#g-tick), logging.

**Rule.** Hence the [epoch rule](#g-input-epoch): a handler executes against exactly the
world its guard fired on. Own `y`, foreign `u`, own `x`/`m` alike are the firing round's
sweep, so `y = h(x)` holds at every handler entry.

[§8.6][s8-6] settles the iteration itself: how far it runs, and how often each event may
fire under the [firing budget](#g-firing-budget) (the per-boundary cap on how often each
event fires). Two of its rules matter here. Within a round each component
fires at most one event, declaration order picking among that component's
simultaneously-eligible events. And same-component sequential composition happens
*across* rounds, each later event re-decided against the post-transition sweep rather
than fired on a stale premise.

#### Why derivatives may read outputs

**Departure from the orthodox formalism, stated openly.** The textbook form is
$\dot{x} = f(x, u)$, $y = g(x, u)$; this design's `f` receives the orthodox arguments
*plus* the published table: $\dot{x} = f(x, m, y, u, t)$. The composite map
$x \mapsto \dot{x}$ is mathematically identical (linearization, trim and AD are
untouched). The heterodox element is only that derivatives may read outputs.

The teaching line: *"stage 1 publishes what you know from state alone; stage 2 adds what
needs inputs; your dynamics read your own published results instead of recomputing
them."*

**Why.** The decision was grounded in a component-by-component survey of
FlightPhysics/FlightApps ([§15.2][s15-2]). Derivative/output overlap is the *norm* in
this domain — Newton–Euler, kinematics, piston engine, gear friction, every discrete
compensator — which is what makes the orthodox split expensive here (row 15).
FlightCore's fused `f_ode!` already embodied the same economics; this design keeps them
while adding checked scheduling.

**Shared expensive computations** are thereby solved uniformly: compute once in stage 2,
publish, and let `f`/`g` consume the ports. External consumers read the same ports — an
accelerometer model reading `f_c_c`, for instance. The **computer/integrator split**
remains fully expressible without framework support ([§7.4][s7-4] carries the full
statement, including when the factoring earns its keep). Purity rules forbid the third
classic resolution, mutable caching, by design.

### 5.4 Artificial loops and the escape hatch

A [component](#g-component) that bundles a no-[feedthrough](#g-feedthrough) output with a feedthrough output in one
atomic evaluation unit can be **[port](#g-port)-level acyclic yet unschedulable** — Simulink's
"artificial [algebraic loop](#g-algebraic-loop)". The canonical instance in this domain is rigid-body
dynamics: velocity out is pure state, acceleration out is feedthrough from total
force. The two-stage split resolves it, and it is the rung that absorbs most of the
class. The `VehicleDynamics` instance ([§15.1][s15-1]) is velocity state-only with
accelerations feedthrough, and it simply dissolves under the split.

What survives the split is the case where a single component's stage-2 outputs
cross-couple through a neighbor: port-level acyclic, stage-level cyclic. The tracer
([§5.6][s5-6]) labels that case **artificial**. Two remedies apply, in this order:

- **Re-factor the [contract](#g-contract).** Before moving any code, re-examine the cycle's wires. An
  input the neighbor consumes *only in a fallback branch* is the archetypal false
  dependency: the neighbor is computing, on the component's behalf, a fallback whose
  semantics belong on the component's own side of the [boundary](#g-boundary). Move the branch to its
  natural owner and the wire disappears.

  The canonical instance is the landing gear's strut/steering pair. The steering
  model consumes the contact-point velocity azimuth `ψ_v` only in its disengaged,
  castoring branch. Castoring, however, is free-swiveling wheel physics: the strut's
  business, not the steering law's.

  ```
  # before
  strut ──ψ_v──▶ steering        # ψ_v consumed only by the castoring branch
  steering ─────▶ strut          # the pair of wires is the cycle

  # after
  steering ──(engaged, ψ_cmd)──▶ strut
  strut:  ψ_sw = engaged ? ψ_cmd : ψ_v    # nothing flows back
  ```

  Re-factoring the steering contract to emit `(engaged, ψ_cmd)`, and computing
  `ψ_sw = engaged ? ψ_cmd : ψ_v` inside the strut, deletes the backward wire
  outright. The factoring survives substitution, which is the test that it records
  structure rather than dodging the diagnostic: a stateful steering actuator produces
  `ψ_cmd` from its own state and still needs nothing from the strut ([§16][s16]
  records the migration).
- **Split the component.** This is the residual remedy, taken when both halves
  genuinely belong to the component and the split documents real structure. Its cost
  is stated where it bites: visibility ([§11.3][s11-3]) is binary. Every intermediate
  shared across the new boundary therefore becomes
  `output_types` — public, connectable, substitution-relevant. The mitigating
  idiom is the granularity guideline ([§4.3][s4-3]), which the split case
  satisfies trivially: one producing stage, one consumer. The idiom spells out as
  **one struct-valued bundle port**, a `StrutGeometry`-shaped value, not N loose
  ports. The bundle type is then contract — a real cost, but a bounded and honest
  one. No visibility register is added for the orphaned intermediates: rows 34 and 55
  (`unlisted`, `Private(T)`) stay closed. The `w` channel ([§5.2][s5-2]) is no exit
  either: it hands values between one component's own functions, so there is
  nothing for a wire to carry (row 165).

The build diagnostic offers both exits explicitly: "cycle through `systems/aero` is
artificial at port level — split the component, or narrow the neighbor's contract".
The offending stage `h_xu` is carried as a separate [payload](#g-payload) field rather than dotted
onto the path ([§11.6][s11-6]/[§13.2][s13-2]).

The split is rare, and the ladder is what earns that word rather than asserting it:
the two-stage split dissolves the common shapes, and the contract re-factoring
absorbs the false wires. What is left for the split is cycles whose halves really are
one component's own work.

One consequence of stage-2 conservatism is worth recording. An input consumed only
by `f`, never by `h_xu`, still creates a scheduling edge if the component has stage-2
outputs. In practice such components are integrator-shaped and have no stage-2
outputs, and the remedy, if ever needed, is the same ladder.

### 5.5 Algebraic loop policy: reject at build time

A genuine cycle in the instantaneous dependency graph is a **build error**. The
diagnostic names the full path in the canonical slash form of [§11.6][s11-6]:
`aero/F → dyn/a → aero/α̇ → aero/F`.

The user breaks it explicitly, by one of three routes: insert dynamics (the α-filter
idiom), insert an explicit unit delay (`UnitDelay`, [§13.7][s13-7]), or restructure.
The α-filter idiom is already standard practice in the domain and in the current C172
model. The unit delay carries a caveat: it changes the model's [tier](#g-tier)
structure. The broken signal becomes discrete, sampled at [`Δt_base`](#g-dt_base) (the
base tick period, an integer multiple `n·h`). That is a modeling decision, not a
transparent wire. Implicit delays and per-step numerical loop solving are both closed
(row 5).

Implicit *algebraic balances* inside a [component](#g-component) (e.g. a turbomachinery operating
point) remain the component author's business: local, owned, bounded. Rejecting
framework-level loops does not forbid such models.

### 5.6 Diagnostics: feedthrough tracing

Tracing is **diagnostic only, never load-bearing**: scheduling correctness comes
exclusively from the structural two-stage split; tracing improves error messages and
verification. Triggered when the scheduler finds a cycle, to classify it (genuine →
"insert a state"; artificial → the remedy ladder, [§5.4][s5-4]).

**Detection and naming.** A cycle surfaces as a topological-sort stall in
[Stratum](#g-stratum) B (one of the build's three phases: structure, schedule, activation).
The stalled subgraph is decomposed into **strongly connected [components](#g-component)**, and
each nontrivial SCC names one cyclic cluster exactly: one diagnostic, its
members and the wires among them, presented as one readable loop in the
canonical slash form ([§11.6][s11-6], `aero/F → dyn/a → aero/α̇ → aero/F`). Neither the
raw stall residue nor a single back edge names the cluster correctly (row 12).

**Classification is [schedule](#g-schedule)-free.** It runs inside Stratum B's failure path,
where no schedule exists — and needs none, because each SCC member is evaluated
*once, in isolation*, at the [probe](#g-probe) point: state views from `init_*`,
out-of-cycle [cells](#g-cell) from the acyclic prefix's probe values, in-cycle cells
synthesized through `probe_value` ([§12.3][s12-3]) under tracer tags. The tracer's
product is a per-member dependence set rather than a value, so no ordering has
to be valid for the labels to come out right. The loop is **real** iff every hop
of the structural cycle survives in the traced per-member maps; **artificial**
([§5.4][s5-4]) iff some hop dies — the component whose stage-2 function does not in fact
route that input to that output. No Stratum C machinery is touched: no
[activation](#g-activation), no layouts, no table. This is the *local* variant (row 12) — the
schedule-free per-member trace at the probe point, which is what the cycle
classifier uses; the "tracer activation" ([§12.4][s12-4]) names the other variant (row 12), the
global set-tracer run as an ordinary Stratum-C activation, and the two must not
be conflated.

**Caveats, carried in the diagnostic rather than assumed away.** The trace
speaks for the branch taken at the probe state (the diagnostic-only doctrine,
row 12). Discrete members trace *structurally* — the discrete [tier](#g-tier)'s
plain,
wholesale-pinning declarations admit no tracer scalar — which is sound as a
may-depend answer but never sharp,
so the remedy hint — split this component, *or* narrow the neighbor's [contract](#g-contract) when
the dead hop's input is consumed only in a fallback branch (the ladder, [§5.4][s5-4]) — is
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
  taken paths, sampled across randomized states; its only misses are untaken branches
  (row 12).

[Boundaries](#g-boundary): only *inputs* are seeded, so branching on state/modes/parameters/time never
interferes (stage-2 functions also receive state views, but neither those nor `y_x`
are ever seeded). Stage-1 functions are never traced (nothing to
seed). Derivatives, [guards](#g-guard), handlers, [projections](#g-projection) are outside tracing's jurisdiction
entirely. Known tracer blind
spot: value-severing operations (dependence passed through a bare `Int` index, e.g.
nearest-neighbor lookup) — documented; linear/cubic interpolation is immune (dependence
flows through the fractional weights).

Both modes ride the same `T <: Real` genericity as `Dual`; Dual-cleanliness in CI
effectively guarantees traceability.

---

## 6. Composition: connections, aggregation and hierarchy

[Components](#g-component) become a system through wiring: connections that route signals
across the [assembly](#g-assembly) hierarchy, and ordinary junction components wherever
several signals must combine into one. [§6.1][s6-1] gives the connection and
hierarchy rules; [§6.2][s6-2] gives the aggregation idiom they force.

### 6.1 Connections and hierarchy

- **Deep connection paths** are allowed, with one structural rule: both endpoints must
  resolve **within the [assembly](#g-assembly) type being defined**. You may deep-route into structure
  you declared (`Cessna172` routing its single `trn` input to
  `systems/ldg/{left,right,nose}/trn_field` in one visible block — no per-level
  re-exports); you may only connect [port](#g-port)-level into submodels held **generically**
  (`World` connects `terrain/field => aircraft/trn` and knows nothing more). This kills
  the re-export ceremony where it is ceremony, and preserves substitutability where the
  [boundary](#g-boundary) is load-bearing. Operationally: a path may traverse any chain of
  concretely-typed fields and stops at the first generically-held child, whose
  [faces](#g-face) are the only things addressable beyond that point — a rule about
  the *declaration's* knowledge, not the build's (a deep path into a generic child
  is forbidden even where the concrete instantiation would resolve it, because it
  hard-codes one implementation and breaks on substitution). Enforcement lives
  in the path-resolution primitive itself (`resolve`, [§13.3][s13-3]), which walks
  declared field types alongside instances.
- Paths are validated at build time; renames break loudly.
- **Two clauses type-check a wire** ([§11.2][s11-2]). The **nominal bound check** is
  stated over declaration evaluations: the producer's declaration
  at `Float64` must be `<:` the consumer's entry at `Float64` — one uniform rule,
  degenerating to exact equality for a concrete entry, violated as
  `WireTypeMismatch`. Beside it, and **for a continuous consumer only**, the
  **walk-compatibility clause**: a walking producer leaf (the producer declared
  `T` there) requires a `T` entry, while a [pinned](#g-walked) producer leaf satisfies either,
  frozen values embedding upward under any [activation](#g-activation). Both sides are declaration
  functions of `T`, so the clause is decided in [Stratum](#g-stratum) A by evaluating them at a
  marker scalar — declaration reading, no user stage code ([§12.1][s12-1]) — and a
  violation is `WalkingFaceAtFrozenEntry`, naming both endpoints, the leaf and
  both declared leaf types, with both remedies in the message: declare the entry
  `T` if the consumer promotes, or feed it from a non-walking source if the freeze
  is genuine. **The [tier](#g-tier) scope is load-bearing, not tidiness.** A discrete consumer
  takes the bound check alone, because its stages read exclusively at real [ticks](#g-tick)
  in the [nominal](#g-nominal) world — a `Dual`-carrying [cell](#g-cell) exists only inside activations the
  discrete tier never runs in ([§12.4][s12-4]) — so a continuous producer feeding a
  discrete consumer is unconditionally legal (row 167). The clause is also what gives the
  two [contract](#g-contract) sides their **failure asymmetry**: the input-side forgotten-`T` —
  the habitual `Float64` written at an entry whose consumer really promotes —
  fails at the *first nominal build*, at the wire, with both endpoints named,
  because an input has a build-time counterparty; the output side has none, so its
  forgotten-`T` lurks (loudly, never silently) until the first `Dual` activation
  ([§11.2][s11-2]).
- Fan-out is free (one producer, many consumers). The converse is strict: every
  input port takes **exactly one** connection, no exceptions (aggregation is
  junctions, [§6.2][s6-2]). The rule spans levels: an input fed both inside a sub-assembly
  and by an ancestor's deep route is a two-producers build error — deep routing
  cannot silently double-feed.
- No auto-bubbling of unconnected inputs (row 43).
- Unconnected output ports: legal, silently, with no build-time warning (row 84).
  Unconnected input ports: build error
  (no silent defaults). **The check is a whole-tree property, not a
  per-declaration one**: within a single assembly declaration an unfed child input
  is simply *awaiting a claim from above* — a sibling wire, an ancestor's deep
  route, or an `input_connections` entry handing the obligation up one level ([§11.6][s11-6]). The
  error fires at the root build for any input whose obligation chain never
  terminates. The one legitimate terminus fed by no [component](#g-component) is the root
  assembly's own input face — a root slot ([§9.3][s9-3]).

### 6.2 Aggregation: explicit summing junctions

N-to-1 physical aggregation (total wrench, total mass properties, total internal
angular momentum — today's generated `get_wr_b`/`get_mp_b`/`get_hr_b` tree walks) is
expressed by **ordinary junction [components](#g-component) and explicit wires**. There is no
framework aggregation mechanism: no multi-connection [ports](#g-port), no declared fold ops, no
identity-element opt-outs. Every input port takes exactly one connection, everywhere.

```julia
struct SumJunction{W, N} end        #type constructor, arity; library-provided

input_types(::SumJunction{W, N}, ::Type{T}) where {W, N, T <: Real} =
    NamedTuple{ntuple(i -> Symbol(:in, i), N)}(ntuple(_ -> W{T}, N))
output_types(::SumJunction{W, N}, ::Type{T}) where {W, N, T <: Real} = (; Σ = W{T})
h_xu(::SumJunction, (; u)) = (; Σ = +(u...))
```

(The parameter is the *unparametrized* type constructor — `SumJunction{Wrench, 3}`;
UnionAlls are legal type parameters — so both [contracts](#g-contract) derive their entries from
it by applying it to the [activation](#g-activation) scalar: the junction is a continuous leaf, so
its `input_types` entries are the tolerant `W{T}` a promoting consumer writes
(walking, frozen and root-[slot](#g-slot) contributors all admissible behind them) while
`output_types` re-types the output [cell](#g-cell) per activation ([§11.2][s11-2]). This is the same
arity-via-computed-contracts pattern [§13.7][s13-7] commits to for `Or{N}`.)

Wired at an ownership [boundary](#g-boundary), the junction is ordinary structure — with
`wr_sum::SumJunction{Wrench, 3}` a field of `Systems` like any other child
([§11.5][s11-5]):

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
  wire.
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
nothing — and each [assembly](#g-assembly) that *owns* contributors aggregates them with an internal
junction and **exports the total** ([§3.3][s3-3]: the junction is a component inside the
assembly; the assembly exports its `Σ` port). The [§6.1][s6-1] connection rules force this
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
junctions themselves — [summing junctions](#g-summing-junction), Bool gates — are the seed of the
standard component library committed in [§13.7][s13-7]: ordinary components, no
framework privileges, inventory grown strictly by migration demand.

**The zero-contributor end of the same spectrum.** Ragged contributors bottom out at
none: a configuration in which a consumer's required aggregate input has *no* physical
contributors at all — the bare-propagation `Vehicle{NoVehicleSystems}`, zero
contributors to external wrench and to internal angular momentum, while
`VehicleDynamics` requires both unconditionally. There is no junction to write and no
producer to wire, [§6.1][s6-1] bans unconnected inputs and silent defaults, and the identity
element a zero-arity junction would need is deliberately absent (row 37). The spelling
is a library `Constant` source ([§13.7][s13-7]) wired straight to the consumer's input —
`Constant(Wrench())` → `dynamics/wr_ext` — so the zero total becomes declared
structure, the configuration stating "external wrench ≡ 0" as a visible wire and an
observable port rather than as an identity method the framework supplies behind the
author's back. This is not the banned default ([§6.1][s6-1]) in component clothing but its
opposite: that default is silent and consumer-declared, this one is loud and
assembly-declared, the author writing the child and the wire, both inspectable.

The cost of explicit wiring, recorded: adding a deep contributor edits one assembly
level (its owner's wiring) instead of zero, and buys in exchange per-contributor
values and intermediate totals as observable ports, with every silence inverted into
a warning or error (rows 7, 37).

Consumer-declared folds with multi-connection legality are closed (rows 7 and 37).

---

## 7. State and data representation

### 7.1 Continuous state: structured immutable, flat backing

Each [continuous component](#g-continuous-component) declares its state by value (`init_x`, [§11.2][s11-2]): a
NamedTuple whose leaves are drawn from a **deliberately closed vocabulary —
plain real scalars and `SArray`s (static vectors and matrices) of a common
eltype `T`** — and nothing else. `Int`s/enums/`Bool`s belong in modes, and
domain wrapper types (`RQuat`, `Ranged`) are not state leaves — an attitude
state is an `SVector{4,T}`, cast where rotation semantics are wanted (below).
The framework:

- computes a **flat layout** at build time (compile-time offsets over one contiguous
  `Vector{T}` [buffer](#g-buffer) it owns);
- **reconstructs** the typed immutable state value for a [component](#g-component) at each evaluation
  and passes it to every function receiving state views ([§5.2][s5-2] argument rule — field
  loads at known offsets, register-level, zero cost);
- receives immutable results back: derivative functions return an `Ẋ`-typed value
  (scatter-stored into the flat `ẋ` buffer); event handlers and [projection](#g-projection) return a new
  `X` (written back — projection at the two [schedule](#g-schedule) positions, [§5.3][s5-3]).

**What `Ẋ` is.** With the leaf vocabulary closed, the answer takes one line:
`Ẋ` has exactly `X`'s shape at the [activation](#g-activation) scalar — a scalar leaf's
derivative is a `T`, an `SArray` leaf's is the same `SArray` at `T`. (This is
the vocabulary rule paying rent: an invariant-carrying leaf like a unit
quaternion has a derivative off its own type, and `Ẋ` would need a separate
derivation; here the attitude leaf is an `SVector{4,T}` and so is its rate.)
The conformance [predicate](#g-predicate) is structural — *each field of `f`'s return
scatters into its field's block at `T`* ([§12.5][s12-5] states the check) — which
makes derivative completeness a property of the layout rather than of author
discipline. There is deliberately **no `derivative_type` hook** (row 190).

**The buffer is authoritative; typed values are ephemeral reconstructions.** Nobody
outside the framework ever holds a mutable reference to state. "Ephemeral" is
literal: an isbits view materializes in the caller's frame (registers or spilled
stack, the compiler's business) for exactly the duration of the call and has no
existence between calls — re-materializing is the same loads, value-identical
because the value is immutable and the buffer unchanged within a [sweep](#g-sweep). Whether
repeated reads within a sweep re-materialize or reuse the loads is codegen freedom
in the literal sense: the [executor](#g-executor) (the compiled execution form of the
schedule) is spelled rebuild-per-call, and hoisting is the code generator's
CSE, whose legality condition is exactly the buffer-unchanged-within-a-sweep
rule ([§12.7][s12-7]). The complementary rule ([§5.2][s5-2]): **[one home per datum](#g-one-home-per-datum)** — buffer
for continuous `x`, stores for discrete `x` and for `m`, table for produced
signals — and no store ever mirrors another; in particular there are no state
[cells](#g-cell) in the table beyond [contract](#g-contract)-driven [auto-published ports](#g-auto-published-port), which are
interface, not transport.

**Why the vocabulary is closed: views must materialize without running
anyone's invariants.** Scalars and `SArray`s have invariant-free
constructors — `SVector`'s stores its tuple, `NamedTuple` construction runs
no user code, nothing normalizes or clamps — so building a view through
ordinary public construction is bit-faithful automatically:
`reconstruct(flatten(x)) == x` identically, with no constructor bypass, no
`reinterpret`, and no reliance on a custom struct's memory layout mirroring
the buffer's. Invariant-carrying leaves are closed (row 94). Domain semantics
are instead an **explicit, invariant-free cast at the
point of use** — `q = RQuat(x.q, normalization = false)` — the conversion
today's `f_ode!` code performs on its raw views, now visible and chosen.
Invariants live where the design already put them: in `project` at
[boundaries](#g-boundary), and in writers — handlers build their returned values through
ordinary constructors, and the condition apply converts authored values
through ordinary `convert` methods ([§14.3][s14-3]). Constructors run on the write
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

The state [buffer](#g-buffer), pack/unpack machinery, and the entire **continuous evaluation path**
are generic over `T <: Real`. One design property, four consumers:

1. exact Jacobians for **linearization** (ForwardDiff duals through the whole model,
   replacing finite differences),
2. derivatives for **trim** solvers,
3. the **[feedthrough tracer](#g-feedthrough-tracer)** ([§5.6][s5-6]),
4. a trivially checkable **CI invariant**: one evaluation [sweep](#g-sweep) with `T = Dual` fails
   loudly (`MethodError`/`InexactError` at the offending line) on any Float64-pinning.

Consumer 1 is also where the *discrete* side's exemption stops being a
limitation and becomes the exact answer: a frozen discrete [cell](#g-cell) is a constant
with zero partials, which is what "linearize the continuous dynamics with the
discrete state held" means (`frozen_discrete_walkthrough.md` works the chain through in
detail).

The declaration layer keeps this scoping legible without putting it in the
author's way: a continuous producer's output declaration is a function of the
[activation](#g-activation) scalar ([§11.2][s11-2]), and cell types per activation are that declaration
*evaluated* at the scalar — participation authored per leaf, `T` where the leaf
follows the activation and a concrete type where it is deliberately [pinned](#g-walked). The
state type is still derived: the framework walks the `init_x`-derived type, real
leaves and `Real` type parameters following the scalar. The discrete side stays
plain and pins wholesale. Nothing anywhere comes from inference through user
code. (Safety of the substitution rests on
the [§12.5][s12-5] embedding guarantee.)

Scoping (what actually needs genericity — roughly half the type inventory):

- **[Walked](#g-walked) — [payload](#g-payload)/value types constructed during evaluation (~25 structs):**
  the quaternion/attitude family (`Quaternion` becomes
  `Quaternion{N,T} <: AbstractVector{T}` — by invariance, `Float64` instances still
  match every existing `AbstractVector{Float64}` method, so existing behavior is
  untouched), `Wrench`, `FrameTransform`, `MassProperties`, `KinData`, `AirData`,
  geodesy value types, `TerrainData`, continuous output structs. Mechanical
  parametrization; constructors infer `T`, so call sites don't change; `@kwdef`
  defaults pin the no-argument case to `Float64`.
- **Pinned — parameters and definitions:** stay `Float64` (promotion handles mixing);
  no migration.
- **[Exempt](#g-walked) — the discrete side** (compensators, avionics): linearization and
  trim differentiate continuous dynamics only.

Lookups: **table data is a pinned parameter; the query coordinate is walked traffic.**
Interpolations.jl evaluates generically over the coordinate (`itp(x::Dual)` works
through the `BSpline`/`scale`/`extrapolate` compositions in use). Caveats: `Linear()`
[interpolants](#g-interpolant) have kinked derivatives at knots (no regression vs. finite differences;
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
  discipline as the table's [cells](#g-cell), in a separate home — [§4.1][s4-1]'s vocabulary),
  overwritten by the framework when an update/handler returns a new value. They never
  touch the integrator [buffer](#g-buffer); no arithmetic is ever done on them.
- **Type freedom:** a discrete leaf's `x` may be *any immutable value* (frozen-reference rule; isbits not
  required). Enums, integers, nested structs, RNG state (four `UInt64`s of `Xoshiro` —
  required to live in discrete state for deterministic [replay](#g-replay)).
- **[Snapshots](#g-snapshot) are free:** copying a store copies a reference to immutable data —
  checkpoint/replay of the entire discrete side is "copy the store values."
- **[Workspace](#g-workspace)** (for heavy algorithms, e.g. an n≈20 Kalman filter):
  [component](#g-component)-declared mutable scratch, instantiated by the framework and arriving
  as the `ws` [bundle](#g-bundle) field ([§5.2][s5-2]) in every bundle-receiving function of the
  declaring component (`project` is positional and receives none),
  and **excluded from state semantics** — not snapshotted, not replayed, never a
  condition target ([§14.1][s14-1]), must carry no information between calls. The framework
  **never inspects or mutates a workspace** — the workspace is an opaque,
  opt-in escape hatch from value semantics, used at the author's own risk, and
  its rules are [contract](#g-contract), not checks: at call entry, contents are unspecified
  beyond the structure the allocator itself established (a plan or
  factorization configured at allocation is valid from then on); scratch is
  garbage until written this call; nothing a previous call left behind may be
  relied upon. No poisoning of scratch is attempted (row 183).
  Declared **by allocation**: the well-known method *is* the
  allocator — `workspace(c::KF, ::Type{T}) where {T} = (P = Matrix{T}(undef,
  c.n, c.n), x̂ = Vector{T}(undef, c.n))` on the continuous [tier](#g-tier), plain
  `workspace(::C)` on the discrete (the contract declarations' tier split) — called
  once per [activation](#g-activation) and once per scratch-store set ([§14.8][s14-8]), sizes from the
  instance, eltypes from the activation. The `T`-signature was never in doubt
  here — `workspace` is the by-allocation register, so the method is an
  allocator the framework *calls*, not a schema it reads — and it is the
  precedent the register criterion ([§11.2][s11-2]) cites now that row 166 has restored the
  scalar to the by-type register as well: a `T` appears in a signature exactly
  where the author makes a choice with it. Nothing
  downstream derives from a workspace's type, and mistyped scratch detonates
  loudly at the `Dual` [probe](#g-probe). The `undef` spelling is the
  recommended idiom and the sole visible marker that contents are meaningless:
  it puts that fact in the declaration, which is the register this store
  actually lives in — declaration by allocation, never by initial value (row 77).
  Available
  on **both tiers**: nothing in the contract is tier-specific, and a continuous
  workspace simply joins the `T`-generic surface — under a `Dual` activation
  the allocator is called at `Dual`, and the in-place math runs through Julia's
  generic fallbacks (no BLAS; activations probe and linearize, they don't run
  marathons). The calls-per-[boundary](#g-boundary) multiplicity of the continuous side (RK
  stages, localization probes, event re-[sweeps](#g-sweep)) makes the no-information-between-
  calls contract *more* load-bearing there, not less.
- **Blessed idiom — zero-allocation [ticks](#g-tick) with immutable `x`:** do the in-place math
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
- **Double-buffered mutable state**: a possible future extension only,
  deferred (row 13).

### 7.4 The fused-evaluation lineage (prior art and how we got here)

The [§5.2][s5-2] interfaces are the end point of a three-step simplification arc, recorded
here because each step replaced a mechanism with something smaller:

1. **N output groups → exactly two** (row 6), at the price of an occasional
   [component](#g-component) split ([§5.4][s5-4]).
2. **Derivative binding → own-output access** (row 15): passing the fresh
   [signal table](#g-signal-table) to `f`/`g` subsumes the declaration feature, and the
   "binding" becomes a one-line function body.
3. **Separate state arguments → the state decoder** (rows 16, 35) — the step
   later reversed by step 4.
4. **Decoder-exclusive state access → stores-and-views arguments.** Step 3's
   second half was reversed (row 35): once [§11.3][s11-3] made publication a
   deliberate interface act, the identity decode stood revealed as *transport* —
   copying the [buffer](#g-buffer) into [cells](#g-cell) so a buffer view could be replaced by a
   cell view. The fixed point is the argument rule ([§5.2][s5-2]) — zero-copy views of
   the stores a function genuinely reads. What survives of step 3: the uniform
   shapes, the fused economics, and the stage-1 decoder itself (today's `h_x`) —
   no longer the sole state gate, but the no-[feedthrough](#g-feedthrough) stage.

Prior art, for orientation — every causal framework meets the shared-computation
problem and resolves it per its architecture: **Simulink diagrams** make integrators
explicit blocks (derivatives are ordinary wires into `1/s` — the computer/integrator
split is their native idiom); **S-functions and FMUs** use sanctioned *mutable caches*
(DWork vectors; FMI's lazy-evaluation caching) between their `mdlDerivatives`/
`mdlOutputs`-style callback pairs; **Modelica/MTK** write `der(x) = expr` natively
with symbolic CSE. The fused [sweep](#g-sweep) + signal-consuming `f`/`g` is the cache-free
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
real-time; (2) throughput for [unattended runs](#g-unattended-run); (3) **the canary**: an unexpected allocation
is Julia's most reliable symptom of type instability, so a zero baseline makes
`@allocated == 0` a CI-testable invariant that catches inference regressions at the
offending commit.

- **Continuous hot path** (per-stage evaluation), plus everything else that
  runs unconditionally per frame or [boundary](#g-boundary) — [guards](#g-guard) (evaluated every
  boundary, firing or not) and `project` (both [§5.3][s5-3] [schedule](#g-schedule) positions):
  exactly zero, CI-enforced at the [§12.7][s12-7] phase-body [seam](#g-seam) (`phase_bodies`).
- **Periodic [ticks](#g-tick) and event handlers** (episodic execution — a tick when
  due, a handler only on firing): zero by idiom ([workspace](#g-workspace) + [snapshot](#g-snapshot)
  pattern; immutable-value returns); documented tolerance for the rare
  exception, scoped per body by the seam's granularity so it never loosens
  the continuous assertions.
- **Logging**: amortized-zero — snapshots are records stored *inline* in a
  `Vector`, `sizehint!` for the expected duration making regrowth a non-event.
  The inline-storage claim is about the snapshot record's *fields*, not about
  everything reachable from them: a model carrying [§4.4][s4-4] [field handles](#g-field-handle) (heightmap
  terrain, wind grids) has a snapshot type with reference fields, and those
  ride as references to build-time-frozen data — no copy, no per-boundary
  garbage, which is what the allocation claim asserts. What the claim does not
  assert is that the snapshot is `isbits`; the per-boundary allocation cost is
  zero either way, and the summarize-or-skip rule ([§4.4][s4-4]) governs what such a field
  contributes on export.
- **What is not recorded**: event firings. The [log](#g-log) holds boundary snapshots and
  the [trace](#g-trace) holds staged inputs ([§9.2][s9-2], [§9.5][s9-5]); neither carries a per-event
  record. Which events fired at which boundary is recovered by [replay](#g-replay) plus the
  published modes — the honest remedy of [§9.2][s9-2], a mode field declared public is
  in every snapshot. An event-firing stream is a [guarded addition](#g-guarded-addition), not built.
- **Tools where garbage is unavoidable**: arena allocation (Bumper.jl-style) for scoped
  temporaries; scheduled `GC.gc(false)` at frame boundaries to move collection out of
  the critical path. (Julia has no per-object freeing; these are the honest levers.)

---

# Part II — Execution

## 8. Time and execution

### 8.1 Loop ownership: the framework owns the simulation loop

Six activities make up the simulation loop: the [§5.3][s5-3] [boundary](#g-boundary)
sequence, [tick](#g-tick) dispatch, event handling, logging, input staging, and
[pacing](#g-pacing) (waits inserted between completed frames, never altering the
boundary sequence).

**Rule.** All six are **framework code, unconditionally**. The loop is written
here, not assembled out of callbacks registered with a third-party solver.

**Why.** The step-boundary [contract](#g-contract) is this design's central
invariant, and only a loop the framework owns can enforce it by construction
rather than by convention. Choreographing the same sequence as an ordered
`CallbackSet` inside a foreign event loop is rejected on exactly that ground
(row 17).

`OrdinaryDiffEq` is therefore **dropped as a dependency** of the new core (row 17).

### 8.2 The stepper seam

Loop ownership stops at one operation: *advance the continuous state from `t` by
`h`*. That operation is delegated across a narrow internal interface, the
**[stepper seam](#g-seam)**, so that the integration method can be replaced
without the loop changing.

#### What the seam requires of a backend

The seam [contract](#g-contract) has three clauses.

- **Advance by arbitrary `h`.** This is required anyway: the loop lands on
  [tick](#g-tick) [boundaries](#g-boundary), and it resumes from a
  [localized](#g-localized) event time (the crossing instant bracketed by
  root-finding over probe sweeps).
- **Dense output on demand over the last completed step.** Only event
  localization needs it ([§8.4][s8-4]), so it is constructed lazily.
- **One-step methods only.** Event handlers reset state discontinuously, and a
  one-step method restarts from a new state for free. Multistep methods are
  excluded (row 17).

#### Models with no continuous state

A model with no continuous state at all is legal: nothing in [§11.2][s11-2]
requires an `x` block of anyone. Such a model still has to be run.

**The seam is never entered empty.** The framework short-circuits rather than
pushing the corner down the seam. With an empty `x`, integrate degenerates to
advancing `t` to the next boundary, and the stepper is simply not called. No
backend ever faces `N = 0`, and no backend contract has to say what it would do
there.

The ownership rule of [§8.1][s8-1] pays off structurally here, rather than as an
argument: the dummy-`[0.0]` tax an empty state pays under a foreign solver loop
(row 17) is gone at the root, not just the [buffer](#g-buffer) but the step over
it. Everything else about such a model is ordinary. The boundary machinery —
[sweeps](#g-sweep), events, ticks — runs unchanged.

#### The first-cut backends

First cut ships **in-house fixed-step RK4 and Heun** over the flat state buffer.
Together they are ~a hundred lines: trivially zero-allocation, hence auditable
against the CI invariant of [§7.5][s7-5], and trivially `T`-generic. Genericity is
not even required of the stepper, since linearization and the tracer drive the
*sweep*, never the integrator.

An `OrdinaryDiffEq`-backed stepper can exist later as a package extension, if an
offline study genuinely demands adaptive or stiff methods. Per the
guarded-additions rule it is not built until then.

#### Why fixed-step low-order suffices

The domain argument is recorded here because it is load-bearing for the whole
axis.

1. **The closed-loop tick cap.** Every application beyond bare propagation runs
   periodic avionics (50 Hz today), whose commands are zero-order-held signals.
   Integrating past a tick with stale commands is wrong, so the integrator must
   land on every tick boundary regardless of method. Adaptive and high-order
   methods pay off exactly when steps can stretch, and the execution model
   forbids the stretch by construction.
2. **A piecewise-smooth [RHS](#g-flow) starves high order.** Linearly-interpolated
   lookup tables (C¹-kinked at every knot), clamps, friction blends and mode
   branches deny high-order error estimators and implicit-solver Newton
   iterations the smoothness they assume. RK4 at 50 Hz already puts integration
   error orders of magnitude below the model uncertainty of a coefficient-table
   aircraft model.
3. **Stiffness has a remedy ladder.** The fastest continuous dynamics in the
   current codebase — actuator poles ~31 rad/s, gear damper decay, friction
   compensators — sit inside RK4's stability region at `h = 0.02`, and the
   crosswind-landing demo is the empirical proof. If a future model exceeds that
   region, the ladder runs in order: first shrink `h` (the RHS costs
   microseconds, and 500 Hz real-time is unremarkable), then subcycle the stepper
   against the tick grid, and only then reach for an implicit method through the
   adapter. If that day comes, the eltype genericity of [§7.2][s7-2] supplies
   exact ForwardDiff Jacobians through the sweep for free.

### 8.3 Signal-table consistency is a boundary property

During a step, RK stages evaluate the [interior sweep](#g-sweep) ([§8.5][s8-5]) at internal
stage states — the [signal table](#g-signal-table) is transiently **integrator scratch**. The [boundary sweep](#g-sweep) in the [§5.3][s5-3] sequence
is what restores consistency at each accepted [boundary](#g-boundary). The rule, binding for the [periphery](#g-periphery) ([§9][s9]):
**external readers (GUI, logging, network output) observe the signal table only at
step boundaries.** Mid-step contents carry no meaning.

### 8.4 Localization mechanics

A [guard](#g-guard)'s [predicate](#g-predicate) can cross inside an integration step, strictly
between two grid points. Two answers are available. The crossing can be noticed
at the step's end, at grid resolution; or the crossing instant itself can be
found and published. This section fixes which guards get which answer, and
spells out the machinery of the second.

Two words carry the section, and they are not synonyms.

- A **[frame](#g-frame)** is a grid step `[tₙ, tₙ₊₁]`, the unit of scheduling.
  Three things are keyed to it: the input [drain](#g-drain) (the frame-top swap
  that publishes staged device inputs into the root slots, [§9.4][s9-4]), pacer
  deadlines ([§8.7][s8-7]) and [tick](#g-tick) eligibility ([§8.5][s8-5]).
- A **[boundary](#g-boundary)** is a published consistency point, where the
  [§8.6][s8-6] macro-sequence completes and a [snapshot](#g-snapshot) goes out.

Every grid point is a boundary. The converse fails: `t*`, the localized event
time, and [boundary zero](#g-boundary-zero) (the initialization boundary at `t₀`,
[§14.5][s14-5]) are boundaries that are not frame tops.

A frame in which one event localizes runs like this, boundaries in bold:

> **tₙ** → integrate → arrival sweep at tₙ₊₁ → trigger → θ = 0 probe → bracket
> → root-find → **t\*** → remainder step → **tₙ₊₁**

That chain is the order of operations, not a walk along the time axis. The
arrival sweep at tₙ₊₁ is what raises the trigger, and integration then resumes
from `t*`, which lies behind it.

#### Which guards localize: the form is the policy

**Rule.** Detection policy is declared by the guard's return type, not by a
flag:

- a guard returning `Bool` is **[boundary-detected](#g-boundary-detected)**
  (checked for edges at step boundaries only, no root-finding);
- a guard returning the nominal scalar, the continuous sign form, is
  **[localized](#g-localized)** (the crossing instant bracketed by root-finding
  over probe sweeps).

The build reads the policy off the [probe](#g-probe) it already runs
([§12.3][s12-3], nominal [activation](#g-activation)), so `Event(guard, handler)`
carries no detection keyword (row 179).

**Why.** Localization brackets a root, and only the sign form offers one. The
form *is* the policy, and the illegal pairing is thereby unrepresentable rather
than merely diagnosed.

**De-localizing is a one-line rewrite with no semantic cost.** Cast the guard to
its predicate: return `σ ≥ 0` instead of `σ`. That cast is the definition
([§2.1][s2-1]), so what remains is the same predicate and the same edges,
observed at boundary resolution.

#### Boundary detection is exact for guards over `u` and `m` alone

**Rule.** For a guard reading only `u` and `m`, boundary detection is *exact*.

**Why.** Such a predicate is piecewise frame-constant: `u` steps only at the
frame-top drain ([§9.4][s9-4]), `m` only through handlers, at boundaries. It
cannot cross mid-step, so there is no interior instant for a root-finder to
find. The boundary is not a resolution limit here but the crossing itself, and
localization is machinery that class has no use for.

#### Mixed predicates: the gate idiom

This matters in practice, because most transitions in FlightPhysics mix input
predicates with state thresholds. The piston engine's `starting → running`
fires on `ω > ω_idle && fuel_available`.

**Rule.** The blessed spelling when such a transition wants localizing is the
gate form `(gate) ? σ : -one(σ)`. The `Bool` factors ride the branch, the
continuous factor rides the value.

**Why.** The idiom is honest, not a trick to defeat the check. Probes vary only
θ, with `u` and `m` fixed through a localization, so the gates are constant over
the bracket. σ restricted to the bracket is then the continuous atom,
bracketable as such.

#### The trigger

**Rule.** A localized event triggers when its predicate was
not-[holding](#g-edge-semantics) at tₙ's [quiescence](#g-quiescence) (the fixed
point where a round of handlers fires nothing) and is holding at tₙ₊₁. The tₙ
sample is the event's **[prior](#g-prior)**, its stored predicate sample from the
previous boundary ([§8.6][s8-6]).

This is the directional edge of [§2.1][s2-1], never a bare sign change: a
holding → not-holding transition neither fires nor localizes.

**The trigger check runs against the arrival [sweep](#g-sweep) at tₙ₊₁**, the
sweep that closes the integration step. It therefore runs *before* the due-gated
[boundary sweep](#g-sweep) refreshes any discrete [cell](#g-cell). The ZOH clause
below already forces that ordering, since probes must see the values the frame
actually held. Stating it here fixes the sequencing at the top: every `t*`
firing precedes tₙ₊₁'s boundary sequence entirely.

#### The localization loop

One localized event within one frame, sketched; every step is normed below.

```julia
# θ = (t − tₙ)/h, and every σ(θ) is a probe: write the state, run the
# interior sweep, evaluate the guard.
σ₁ = σ(1)                              # already computed by the arrival sweep
(prior not-holding && σ₁ holding) || return    # no trigger, nothing to localize

σ₀ = σ(0)                              # the θ = 0 validation: x̂(0) = xₙ, no interpolant
σ₀ holding && return                   # epoch-caused edge: fall through to tₙ₊₁

build x̂ from (xₙ, ẋₙ, xₙ₊₁, ẋₙ₊₁)      # one sweep for ẋₙ₊₁, paid only here
lo, hi = 0, 1                          # the not-holding / holding bracket, in θ
while hi - lo > localization_tol       # relative: bracket width (hi − lo)·h vs. tol·h
    θ = next bracketed guess in (lo, hi)       # ITP, Brent or bisection
    σ(θ) holding ? (hi = θ) : (lo = θ)
end
t* = tₙ + hi·h                         # the holding endpoint of the final bracket
```

**The θ = 0 validation.** On trigger the *first* act is a probe at the left end:
write xₙ into the state [buffer](#g-buffer), run one [interior sweep](#g-sweep),
and evaluate the guard to get σ₀. Nothing new is kept, xₙ being already retained
by the stepper for the [interpolant](#g-interpolant) (the cubic Hermite
continuous extension over the last completed step). The probe runs before any
interpolant cost, because x̂(0) = xₙ identically and no interpolant is needed to
place it.

It pays for itself twice over. σ₀ is the **left bracket value** the value-based
root-finders need: σ₁ = σ(tₙ₊₁) is retained from the arrival evaluation, while
σ₀ was otherwise unsourced. And σ₀ **discriminates the edge's cause**.

The discrimination needs one term. An **[input epoch](#g-input-epoch)** is a
maximal span of constant `u`, delimited by frame-top drains ([§9.4][s9-4]).
Within an epoch a guard can change only through the trajectory; at a seam it can
jump without crossing anything.

**Why the discriminator is conclusive.** `u` is the *only* thing that can differ
between the prior's evaluation context and this probe. `m` changes only via
handlers at boundaries, and priors are sampled at quiescence, after them.
Discrete cells ZOH-hold, and the interior sweep excludes discrete entries
([§8.5][s8-5]). `t = tₙ` exactly, by the indexed-grid rule below. And sweeps are
deterministic. So under the honest [priors](#g-prior) of [§8.6][s8-6] the
frame-top drain is the sole possible source of disagreement.

- σ₀ **not-holding** ⇒ a **trajectory-caused** edge, a genuine in-frame
  crossing. Pay ẋₙ₊₁'s sweep, build the interpolant and root-find on the bracket
  $(t_n, \sigma_0)$/$(t_{n+1}, \sigma_1)$.
- σ₀ **holding** ⇒ an **epoch-caused** edge. The drain flipped the guard at the
  frame top, and σ holds at both ends, so no in-frame crossing exists to find.

**An epoch-caused edge is discarded, not degraded.** The localization is
**discarded** and the event fires inside tₙ₊₁'s ordinary iteration.
Mechanically, *not localizing is the action*: fall through, and the boundary
iteration detects and fires it like any boundary-detected event. This path costs
one interior sweep, never ẋₙ₊₁ and never an interpolant, and it consumes no
`localization_budget`. It also **warns nothing**: input timing is a frame fact,
the same doctrine that forbids draining at `t*` below, and boundary detection is
*exact* for a `u`-caused edge (above; row 179), so boundary firing is the
correct semantics rather than a degradation. It is the left-end mirror of the
`t* = tₙ₊₁` degeneracy below.

**The interpolant is lazy.** It is the cubic Hermite continuous extension
$\hat{x}(\theta)$, $\theta = (t - t_n)/h \in [0, 1]$, built from
$(x_n, \dot{x}_n, x_{n+1}, \dot{x}_{n+1})$. $\dot{x}_n$ is the step's first
stage; $\dot{x}_{n+1}$ costs one sweep, paid only on a validated trigger. The
θ = 0 probe precedes it, so an epoch-caused edge never pays it at all. Uniform
accuracy is $O(h^4)$, one order below the discrete solution and the standard
pairing. The event time can only ever be as accurate as the interpolant, so
nothing more expensive is worth probing.

**Probes run the interior sweep.** Guards read `y`, so evaluating a guard at an
interpolated state means writing $\hat{x}(\theta)$ into the state buffer and
running the interior sweep. This is the rule the [RHS](#g-flow) already lives
under ([§8.5][s8-5]): a probe is a mid-step evaluation. Discrete cells therefore
hold their [tick](#g-tick) values through localization, and a guard reading a
sampled output sees what the controller is holding. One interior sweep per
probe.

**Root-finding is bracketed and derivative-free** (ITP or Brent; bisection is an
acceptable fallback). The observed not-holding/holding bracket *is* an
unconditional convergence certificate. Newton/AD localization is rejected
(row 18).

**Convergence is a relative bracket width.** Localization stops once the bracket
is narrower than `localization_tol · h`, with `localization_tol` a `Simulation`
deployment keyword defaulting to `1e-6`. Relative, because an absolute-in-`t`
tolerance is not scale-free (row 133). `1e-6`, because the event time can only
ever be as accurate as the interpolant — `O(h⁴)`, above — so at practical `h`
anything tighter buys nothing while every probe costs a full sweep. The bill is
a handful of probes under ITP, ~20 even in bisection's worst case.

**Post-event.** The boundary sequence runs at `t*` (below) → **interpolant
invalidated** → resume integration from `t*` with the
[remainder step](#g-remainder-step) targeting tₙ₊₁ → re-check guards on the
remainder. The interpolant is invalidated because the handlers made it a lie for
`t > t*`. The re-check runs under the per-frame
[localization budget](#g-chattering), with a chattering diagnostic.

**Multiple events localizing in one step fire at the *earliest* `t*`.** Ties
fire together at that boundary, in declaration order within the iteration
([§8.6][s8-6]). Later crossings re-localize on the remainder.

**A shared blind spot, documented.** An even number of crossings within one step
returns the predicate to not-holding at the boundary, so no edge is observed.
That defeats detection under both policies. The mitigation is step size, not
machinery.

#### Endpoint policy and grid integrity

**Rule.** The root-finder returns the **holding endpoint of its final bracket**,
the smallest probed point where the predicate holds.

**Consequence: `t* = tₙ` is structurally impossible**, not clamped away.

**Why.** The argument rests on **observation, not testimony**. Root-finding is
entered only after the θ = 0 validation has *measured* σ₀ not-holding under the
frame's own `u`. The bracket's left end is therefore strictly not-holding by the
same kind of evidence as its right end, and the returned point is strictly later
than the published, immutable tₙ (worst rounding: `nextfloat(tₙ)`). This holds
unconditionally, with no appeal to the prior and no residual epoch hole: the
case where the prior's testimony and the frame's `u` disagree is exactly the
epoch-caused edge, which never reaches the root-finder.

The guard also observably *holds* at `t*`. Handlers therefore fire in states
where their own predicate holds, and the post-fire prior records an actual
observation rather than an assumption.

**`t* = tₙ₊₁` exactly is legitimate**: a crossing at the grid point, where
σ(tₙ₊₁) = 0 both triggers detection and is the root. It **degenerates to the
grid boundary** — the localization result is discarded and the event fires
inside tₙ₊₁'s ordinary iteration. That outcome is bitwise the boundary-detected
one: one boundary, one snapshot, no zero-length remainder.

**Grid times are indexed, never accumulated.** A near-degenerate `t*` leaves a
tiny remainder step. That is numerically harmless, increments being `h′`-scaled;
the real hazard is bookkeeping, and this rule kills it. `tₖ = t₀ + k·h` is
computed from the frame index (tick gating is already counter-modulo,
[§8.5][s8-5]), and the remainder step *targets the grid point*, with `h′`
derived at use. `t*` is a float inside a frame, never an anchor anything else is
computed from.

#### What a `t*` boundary does, and does not, do

At `t*` the full [§8.6][s8-6] event phase runs: [sweep → guards → handlers]
iterated to quiescence, with **firing-budget accounting scoped to this
boundary** — fresh again at tₙ₊₁, and at a second `t*` on the remainder. The
settled state is then **published**: snapshot, [§10.3][s10-3] boundary-counter
increment, [`stop_on`](#g-stop_on) check ([§13.5][s13-5]). A crash localized at
`t*` ends the run from that snapshot.

**What does not happen at `t*`.** Ticks are never due there: `t*` is off the
[harmonic grid](#g-harmonic-grid) (every discrete period an integer multiple of
`Δt_base`) by construction, and discrete cells ZOH-hold through the sweep.
Staged inputs are not drained either, for two reasons. Input timing is a frame
fact, and replay determinism must not depend on localization arithmetic.

**The `t*` publication is not separately paced.** The pacer paces frame
deadlines, and a `t*` snapshot publishes when computed, mid-frame. Where that
lands in wall-clock time is below what pacing resolves, and the [§8.7][s8-7]
invariant is about trajectories, which are identical either way.

Replay pointers and error messages index boundaries by a monotonic counter with
recorded `t`. The trace stays frame-indexed — `t*` boundaries consume no
inputs.

#### Projection's reach is the boundary, not the probe

**Rule.** Guard probes evaluate the raw interpolated state. Authority rests with
the `t*` boundary: [projection](#g-projection) runs there, and the [§8.6][s8-6]
iteration's edge checks read the *projected* state.

**Why.** RK-stage RHS evaluations already live under that same rule, being
equally off-manifold, so sweeps must tolerate near-manifold states and already
do. Per-probe projection is rejected (row 18).

If projection moves the state back across a guard, the event does not fire and
the run has published one extra boundary. That is harmless, and deterministic
and pace-independent like any other localization outcome (row 80).

#### Budget exhaustion degrades; it does not throw

**Rule.** `localization_budget` is an integer count of localizations permitted
within one frame, defaulting to **8**. It is the second deployment keyword this
section fixes.

**Why 8.** A legitimate multi-event frame — three landing-gear struts touching
down inside one step — needs three or four; [chattering](#g-chattering) needs
tens. 8 bounds the pathology without ever binding on a healthy model.

**When a step spends its event budget**, localization stops for the remainder of
that frame. The remainder step completes, and any further crossings fire in the
next boundary's ordinary iteration — boundary granularity for that frame — under
a `ChatteringBudget` warning ([Appendix C][sC]) naming the chattering event and
the localization count.

The degradation is a function of the trajectory alone, never of wall clock, so
the pace-independence guarantee (row 80) stands untouched and the run replays
identically. A `StepError` would misclassify an expected modeling outcome as
broken machinery (the doctrine, [§14.8][s14-8]).

**The same doctrine governs [§8.6][s8-6].** Neither the boundary iteration nor
cross-frame re-localization has a structural bound, so each takes a budget:
`firing_budget` there, `localization_budget` here. Both degrade loudly rather
than erroring, under a warning that names the offending event. They differ only
in *what* exhaustion sheds. Localization sheds root-finding precision,
preserving every firing at boundary granularity; the firing budget sheds
firings, which is precisely what bounds the iteration.

#### Both constants are deployment, not implementation

`localization_tol` and `localization_budget` are `Simulation` keywords standing
beside `h`, `n` and the algorithm ([§12.1][s12-1], [Appendix B][sB]). They are
validated with their siblings — a positive tolerance, an integer budget ≥ 1 —
and collected into `DeploymentInvalid` ([Appendix C][sC]). The `firing_budget`
([§8.6][s8-6]) stands beside them in every one of these lists: same validation,
same `DeploymentInvalid`, same trace header, same replay comparison. Both are
grid-independent, so neither enters the harmonic-grid check ([§8.5][s8-5]).

**Being trajectory-determining, both are recorded.** They ride the
[trace header](#g-trace-header)'s deployment block and join the set
[replay](#g-replay) compares up front, exactly as `h` and the algorithm do
([§9.5][s9-5], [§10.7][s10-7]).

**Why.** The replays-identically promise above is empty otherwise. A run that
does not record what its localizer was told to do cannot be re-driven through
the same localization outcomes.

### 8.5 Multi-rate tick scheduling

A model runs several clocks at once. The integrator advances on the continuous
step `h`, an inner control loop samples at one rate, an outer loop at another,
a receiver at a third — and each of those must hold its outputs steady between
its own firings. Three things have to be fixed for that to be well-defined: the
time lattice every rate shares, the test that decides which components run at a
given [boundary](#g-boundary), and the surface an author declares a rate on.
This section fixes all three, in that order.

#### The base grid, and the pair every rate compiles to

**Rule.** Every discrete [component](#g-component)'s period is an integer
multiple of a base [tick](#g-tick) period `Δt_base`, and `Δt_base` is itself an
integer multiple of the continuous step ($\Delta t_{\mathrm{base}} = n \cdot h$,
$n \ge 1$). That is the **[harmonic grid](#g-harmonic-grid)**. Ticks therefore
land on step boundaries — the only place anything discrete ever happens
(row 19).

**One pair per component.** However an author declares a rate, and however
deeply the declaration is nested, the build compiles it to two integers per
discrete component: a divisor `D`, the component's period in base ticks, and a
[phase](#g-phase) `Φ`, its offset in base ticks. The pair is kept in the
canonical residue `0 ≤ Φ < D`, so the component's ticks fall at base-tick
indices `Φ`, `Φ + D`, `Φ + 2D`, and so on.

**The gate.** A component is **[due](#g-due)** at a boundary when
`(idx − Φ) % D == 0`, where `idx` is the boundary's tick index. That
subtraction and remainder are the whole admission test: one subtraction more
than a phase-free test would cost, over a lattice fixed at build time. Where a
component's `(D, Φ)` comes from is the declaration surface below.

#### Discrete stages run only at their own ticks

**Rule.** A discrete component's `h_x`/`h_xu` run only at its own ticks, and
its [cells](#g-cell) hold in between — ZOH, stated in [sweep](#g-sweep) terms.

**Why.** Re-running a discrete component's stages at every boundary would
un-sample a sampled-data controller (row 19).

Delivering that hold takes **two statically distinct sweep variants, compiled
from one entry list**. Discreteness is a build-time fact, so the split is
static rather than a runtime test (row 147).

- The **[interior sweep](#g-sweep)** walks *continuous entries only*. RK stage
  evaluations ([§8.3][s8-3]) and localization [guard](#g-guard)
  [probes](#g-probe) ([§8.4][s8-4]) run this variant. The ZOH therefore holds
  mid-step **by construction**: discrete entries are not gated out at runtime,
  they are absent from the walk at compile time. The hot path carries no gating
  test at all.
- The **[boundary sweep](#g-sweep)** walks the full list, with discrete entries
  gated by `(idx − Φ) % D` against the boundary's tick index. It is the variant
  the [§8.6][s8-6] macro-sequence runs. It is not one fixed list either:
  different boundaries run different subsets of the [schedule](#g-schedule).

The split applies to **both sweep blocks**. The discrete [tier](#g-tier)'s `h_x`
entries are absent from the interior stage-1 walk exactly as its `h_xu` entries
are absent from the interior stage-2 walk. The two sweep variants surface in
the phase-body signatures: interior bodies take no arguments, boundary bodies
take the tick index ([§12.7][s12-7]).

#### The due set is a property of the boundary

**Rule.** The due set is computed once for the boundary and reused by every
re-sweep of its [quiescence](#g-quiescence) iteration (the fixed point where a
round of handlers fires nothing, [§8.6][s8-6]). It is a property of the
boundary, not of the sweep call.

**Why.** A due component is at its tick instant for the whole boundary, not for
one round of it.

Three kinds of boundary, three due sets:

- At a **frame top**, the due set is the gate's image of the frame index.
- At a **`t*` boundary**, it is **empty**. The tick counter has not advanced
  there, so no component is at a tick instant. A modulo test against the
  unadvanced index would wrongly re-admit the previous frame's due set.
- At **[boundary zero](#g-boundary-zero)** (the initialization boundary: the
  ordinary macro-sequence with an empty integrate), it is **everything with
  `Φ = 0`**. At `idx = 0` the gate reads `(0 − Φ) % D == 0`, which under the
  canonical residue `0 ≤ Φ < D` holds if and only if `Φ = 0`. The rule is
  implemented by nothing: it falls out of the ordinary gate.

An offset component's first tick is at `Φ·Δt_base`. Until then its cells hold
the values the build probe populated ([§12.3][s12-3], [§14.5][s14-5]), which is
a coherent ZOH story: those are exactly what a tick at `t₀⁻` would have
produced. In a phase-free model every `Φ` is 0, so "at boundary zero everything
is due" is the degenerate case of the same identity.

#### Simultaneous ticks are already well-defined

Several components can be due at one boundary, and settled machinery already
orders them. All due components run their output stages in topological order
within the sweep. All due `g` updates run after it, in any order — each `g`
reads the table and writes only its own `x` store. The FCS cascade's intra-tick
ordering is thereby a sweep property, not an update-order property.

#### Coincidence and stagger are modeling choices with observable consequences

Coincident ticks give a consumer fresh same-instant reads via topological
order, the idealized synchronous-sampling picture. A phase stagger makes the
same reads pipelined and deterministically aged instead. That is the structural
expression of an acquisition pipeline's latency, obtained with no delay blocks.

A stagger is also a load-shaping tool under real-time [pacing](#g-pacing)
(waits inserted between completed frames, never altering the boundary
sequence). Staggered stacks never share a [frame](#g-frame), so worst-case
frame cost is a `max` rather than a sum ([§8.7][s8-7]).

Both patterns are worked in `sample_time_proposal.md`, together with how
silently an offset edit rewires a coincidence structure. The
[bound schedule](#g-bound-schedule) (the printable per-component `(D, Φ, Δt)`
artifact deployment binding produces) and its hyperperiod chart
([§12.2][s12-2]) are how a user audits which pattern a model actually has.

#### Assemblies: virtual for execution, rate scopes for declaration

A [rate scope](#g-rate-scope) is an assembly's `sample_times` declaration
against the enclosing scope. There are no atomic [assemblies](#g-assembly), and
no opt-in variant (row 19).

**Why no coarsening is needed:** the [signal table](#g-signal-table).
Interleaving is semantically invisible under it, consumers reading cells whose
freshness is guaranteed by topological order rather than by contiguity.

#### Declaring a sample time: two registers, one concept

**Rule.** A discrete component or sub-assembly is scheduled by a `sample_times`
entry in its enclosing assembly ([§11.7][s11-7]). The entry declares one
(period, phase) pair in one of two unit systems, and the wrapper type names the
unit system.

| entry | unit system | tick instants | constraints |
|---|---|---|---|
| `Relative(K, Φ = 0)` | scope ticks | every `K`-th tick of the enclosing scope, starting from its `Φ`-th | `K ≥ 1`, `0 ≤ Φ < K` |
| `Absolute(q, τ = 0)` | seconds | `t = τ + k·T`, with `T = period(q)` | `T > 0`, `0 ≤ τ < T` |

`K = 1` therefore admits no stagger. Two same-rate siblings are staggered one
level down instead: declare the scope at twice their rate, then `Relative(2, 0)`
and `Relative(2, 1)`.

`q` is a quantity value, `Period(1//50)` or `Hz(50)` — a spelling choice,
normalized to the rational period at construction. Every period and offset is
an exact `Rational{Int}`, because grid derivation is GCD arithmetic and
ill-defined over floats. A float argument throws the teaching error naming the
exact spelling (`Period(1//50)`, or `Hz(1//2)` for 0.5 Hz).

The wrappers are the whole vocabulary: a bare integer or bare quantity is a
declaration error. An unlisted discrete child defaults to `Relative(1)`, so the
common case costs nothing and a multiplied or anchored child always appears
explicitly.

**Validation belongs to [Stratum](#g-stratum) A**, the build's
declaration-validation stratum, and is collected with path attribution
([§12.1][s12-1], [§13.1][s13-1]). It covers `K ≥ 1`, `0 ≤ Φ < K`, `T > 0`,
`0 ≤ τ < T`, and keys naming discrete or scope children. The constructors
themselves are plain data carriers, with no checks of their own (row 185).

#### The relative register composes affinely and stays on the scope grid

**Rule.** Multipliers compose multiplicatively and phases affinely down the
tree. Under a scope compiled to divisor and phase `(D_s, Φ_s)` in base ticks, a
child declared `Relative(K, φ)` compiles to `D = K·D_s` and `Φ = Φ_s + φ·D_s`.

Composition preserves the canonical residue `0 ≤ Φ < D`, for which
`sample_time_proposal.md` carries the one-line induction. All scoping therefore
compiles away at build to **one `(D, Φ)` pair per discrete component**, and the
boundary sweep gates on that pair with the `(idx − Φ) % D == 0` test above. The
lattice stays static, and the interior sweep still holds no discrete entries to
gate.

**Why relative is the default register.** In a layered control architecture the
*ratios* are intrinsic to the design and travel with the assembly type: inner
loop at `Relative(1)`, outer loops at `Relative(5)`, whatever the deployment.
The convention that keeps `K ≥ 1` livable is that **a scope's base rate is its
fastest relative member**, the member that then gets `K = 1`.

**Two structural properties confine grid cost to the other register.** A
relative phase selects among scope ticks that already exist, so it never
refines the base grid. And it cannot place a tick *between* scope ticks:
staggering off-grid means declaring the offset in seconds, or declaring the
scope base finer than its fastest member so unused [slots](#g-slot) exist.

#### An `Absolute` entry detaches its child from the scope's grid

**Rule.** An `Absolute` entry may appear in any scope's `sample_times`, not
only the root's. The `(T, τ)` pair it establishes is an **[anchor](#g-anchor)**,
and the child hangs from that anchor: it is severed from the enclosing scope's
grid, with no relation to the scope's ticks remaining.

Three corollaries follow:

- `K ≥ 1` reads "a child cannot tick faster than the scope it is *relative*
  to", so an anchored child may tick faster than its scope.
- The fastest-member convention counts relative members only.
- Phase relationships between an anchored child and its relative siblings are
  **deployment-emergent**: whether their ticks ever coincide depends on how the
  grid derivation works out. That is the question the printable bound schedule
  ([§12.2][s12-2]) exists to answer.

Relative children *of* an anchored subtree compose against the anchor exactly
as against the root grid, and a nested anchor simply severs again (the fold,
[§12.1][s12-1]).

**Absolute periods and nonzero offsets jointly constrain the base grid.** They
join the deployment-time constraint pool ([§12.1][s12-1]). This is the subtlety
with teeth: an offset of `T/2` can cost a 2× finer grid, and `T/1000` a 1000×
one. The cost is relational, incurred against everything else declared, which
is why attribution is the engine's job ([§12.2][s12-2]).

#### A declaration, its compiled pairs, and one hyperperiod

Three discrete components under two scopes, at a deployment that binds
`Δt_base = 2 ms` ([§12.1][s12-1]):

```julia
# Root scope: (D_s, Φ_s) = (1, 0).
sample_times(::Vehicle) = (fcs  = Relative(1),        # → (D, Φ) = (1, 0)
                           gnss = Absolute(Hz(50)))   # → (10, 0), anchored at T = 20 ms

# fcs is the enclosing scope of its own children, at (D_s, Φ_s) = (1, 0).
sample_times(::FCS) = (inner = Relative(1),           # → (1, 0)
                       outer = Relative(5, 2))        # → (5, 2), D = 5·1, Φ = 0 + 2·1
```

Those pairs are what the boundary gate reads at run time. One hyperperiod is
`lcm(Dᵢ) = 10` base ticks, and because the gate is pure modulo arithmetic that
one hyperperiod is the complete truth rather than a sample ([§12.2][s12-2]):

```
base tick k:  0  1  2  3  4  5  6  7  8  9 | 0  1  …
inner         •  •  •  •  •  •  •  •  •  • | •  •       (D, Φ) = (1, 0)
outer               •              •       |            (5, 2)
gnss          •                            | •          (10, 0)
```

`outer` is due where `(k − 2) % 5 == 0`, `gnss` where `k % 10 == 0`. Should
`outer` read a `gnss` cell, the ZOH makes that read two base ticks old at
`k = 2` and seven at `k = 7` — the deterministic aging of a stagger, in this
model's numbers.

#### Where the doctrinal line falls on mid-tree anchors

Absolute-first declaration as the default register is rejected (rows 19, 186).
What mid-tree anchors legitimize is narrower.

**Rule.** An absolute declaration inside a library type is legitimate when the
rate is **a fact about the modeled system, not a preference about the
simulation**.

A GPS receiver emitting at 1 Hz, a bus schedule, an ADC pipeline's fixed
conversion offset are as intrinsic to the assembly as its wiring. Forcing them
to the root breaks encapsulation, since the root would have to know
[device](#g-device) internals to re-declare them.

"Run the controller at 400 Hz in this study" remains a deployment choice, and
the existing idiom remains its answer: the assembly exposes its multiplier as a
constructor parameter. Absolute pinning *from outside* a subtree's
[contract](#g-contract) stays rejected as action at a distance.

The framework cannot police the distinction. It is authoring doctrine, recorded
here.

**Anchoring leaves the never-cache-`Δt` argument below fully intact.** The
pinning happens in the enclosing assembly's `sample_times`, the same site where
the multiplier lives. The component type itself therefore stays rate-agnostic,
and still consumes the `Δt` of its [bundle](#g-bundle) (the NamedTuple of
zero-copy views a component function receives).

#### `Δt` has a single source of truth: the compiled schedule

**Rule.** Each discrete component's effective period arrives read-only as the
`Δt` field of every discrete-tier bundle ([§5.2][s5-2]) — `h_x`, `h_xu` and `g`
alike. The field is absent from continuous bundles, so touching it on the wrong
tier is a missing-field error rather than a rule.

**It must be readable in the *stages*, not just in `g`.** Per [§15.2][s15-2],
the discretized laws that actually consume `Δt` — a PID's backward-difference
coefficients, a LeadLag's Tustin transform — run in `h_xu`; `g` is a copy.

The value must arrive through the call, and the bundle field is where. A
`comp.Δt` virtual property is impossible here, not merely inconvenient
(row 19).

**Author rule: never store `Δt`, or any `Δt`-derived coefficient, as a
component parameter.** Recomputing derived coefficients per tick is a few
arithmetic ops. A cached copy is a second thing for gain-scheduling machinery
to chase.

**Relative declaration structurally enforces that rule for the period itself.**
Under scoped multipliers a component author *cannot* know their absolute rate:
it does not exist until composition.

**Phases change none of this.** The bundle's `Δt` is still `D·Δt_base`, an
offset shifting firing instants and never the period, so the discretized laws
([§15.2][s15-2]) are unaffected by staggering.

### 8.6 Event iteration at boundaries: to quiescence, budgeted

This section settles what [§5.3][s5-3] leaves open: how far the event phase runs
at a [boundary](#g-boundary), and how often each event may fire while it does.
The phase **iterates**. One round re-runs the [boundary sweep](#g-sweep),
evaluates all [guards](#g-guard) against it, and fires the eligible events — at
most one per [component](#g-component), each firing being `handler → project`.
Rounds continue to **[quiescence](#g-quiescence)** (the fixed point where a round
of handlers fires nothing).

**Rule.** An event fires in an iteration round iff three things are true: its
[predicate](#g-predicate) is observed holding in that round, the sample observed
before it was not-holding, and the event's firing count for this boundary is
below `firing_budget`. That is the whole of "newly-fired". The predicate is the
one [§2.1][s2-1] defines: the `Bool` form true, or `σ ≥ 0`. `firing_budget` is a
`Simulation` deployment keyword, an integer ≥ 1 defaulting to **4**, and it caps
how many times each declared event may fire at one boundary.

The per-event **loop state** that decides the rule is three registers, all named
normatively:

- the **[prior](#g-prior)** — the previous boundary's quiescent sample, which is
  what the boundary's *first* round tests against;
- the **last-observed sample** — initialized from the prior when the boundary
  opens and overwritten by every round's evaluation, which is what every later
  round tests against; the one exception is an event that was eligible but
  blocked (below), whose sample stands — blocking defers its edge rather than
  consuming it;
- the boundary's **firing count** — incremented at each firing and reset when the
  boundary ends.

Eligibility inside a boundary is therefore an [edge](#g-edge-semantics) like any
other: a not-holding → holding transition, never a bare sign change. What differs
is the reference sample — the edge is read against the last-observed sample, not
against the boundary's entry prior. Two consequences follow directly. Sticky
predicates need no special case: an event that fires and *keeps* holding presents
no further not-holding → holding edge, so it fires once, at the boundary where it
first held. And a predicate genuinely falsified and re-enabled inside the
boundary — its effect reverted by another handler's cascade — fires again, at
this boundary, against a fresh sweep (row 181).

One boundary's iteration, sketched:

```julia
# entering the boundary, per event:  last ← prior,  count ← 0
while the previous round fired something   # the first round always runs
    boundary sweep                         # whole gated schedule, due set fixed for the boundary
    per event:      eligible ← last not-holding && now holding && count < firing_budget
    per component:  firing ← its first eligible event, in declaration order
    per event:      last ← now, unless eligible and not firing   # a blocked edge stays unconsumed
    fire the firing events                 # handler → project, count += 1
end                                        # the exit condition is quiescence
per event:  prior ← last                   # the settled boundary's honest sample
```

**The prior is updated at each boundary's quiescence, from the final
post-iteration samples.** That update is unconditional, with no exception. Every
prior is therefore an honest observation of a settled boundary, which is what
makes the θ = 0 discriminator ([§8.4][s8-4]) conclusive: the
[frame-top drain](#g-drain) is the sole possible source of disagreement between
the prior and the probed left end.

All three registers are detection bookkeeping, not model memory. They are
correctly *not* in any state store — not captured, not traced, reconstructed
deterministically. Beyond the prior, the cost is one `Bool` and one small counter
per event.

**[Boundary zero](#g-boundary-zero) establishes every prior as not-holding.** A
predicate already holding in the authored state therefore fires at `t₀`, and that
behavior ([§14.5][s14-5]) is derived rather than asserted. A warm restart resets
all three registers from scratch, `init!` re-running boundary zero
([§14.5][s14-5]): predicates holding in the newly applied state fire again at the
new `t₀`.

**Why iterate.** Under a single pass, a cascade of N logically-simultaneous
transitions (supervisor FSM → subordinate FSM → …) completes in N steps, at
latency N·h. That makes model semantics depend on the integrator's step size,
`h` being an execution parameter — the same footgun class [§2.2][s2-2] cited when
killing `f_step!`. Cascades are not a corner case either: the blessing of
externalized FSM components ([§3.1][s3-1]) makes cross-component ones the
expected idiom. Orthodoxy concurs. Hybrid automata take sequences of
instantaneous transitions at one time point, Modelica iterates events to
quiescence, and Stateflow runs charts to completion within a [tick](#g-tick).

Boundary-detection timing remains h-dependent, but that is a different quantity:
the *resolution* at which a physical crossing is noticed. The cascade delay would
have been structure the framework inserts between transitions the model declares
simultaneous.

**Why a full re-sweep per round.** A transition reaches the
[signal table](#g-signal-table) *only* through a sweep. A handler writes its
component's state stores and nothing else, so neither the transitioning
component's own [ports](#g-port) nor the downstream stage-2 chains reading them
have moved. A round therefore re-runs the whole gated [schedule](#g-schedule).
The cost is noise: sweeps are microseconds, and rounds beyond the first require
an actual cascade.

**Within-round visibility: the signal table has a single writer, and it is the
sweep.** A handler writes nothing to the table. It returns transitions, the
framework latches them into the component's state stores, and `project`
normalizes them; auto-publication is a sweep act like any other stage-1 write
([§12.5][s12-5]). Nothing moves the table mid-round.

Hence the [epoch rule](#g-input-epoch), which is the whole of this section's
content: **a handler executes against exactly the world its guard fired on.** Own
`y`, foreign `u`, own `x`/`m` are all the firing round's sweep, so `y = h(x)`
holds at every handler entry. No [bundle](#g-bundle) — the NamedTuple of
zero-copy views a component function receives — ever straddles two epochs.
Serialization is what delivers it: a component's state stores are written only by
its own handlers, and it fires **at most one event per round**, so no same-round
writer precedes any handler's entry.

**A component's other eligible events are blocked, not lost.** Each is re-decided
in the next round, against the post-transition sweep. Declaration order is
therefore a **priority with re-decision** rather than a simultaneity. An event
whose premise the earlier transition falsified simply does not fire; under a
within-round sequence it would have fired on the stale premise. Blocking is
register-visible (row 191): an eligible-but-blocked event is the one case whose
last-observed sample is *not* overwritten, so the edge it presented stands
unconsumed. A guard that keeps holding therefore fires in the next round on that
same edge; one the transition falsified records not-holding as usual, and any
later re-rise is a fresh edge. The prior stays honest for free: the quiescent
round fires nothing, hence blocks nothing, so every sample takes its final
update before the prior is written.

**Across components, handler order within a round is semantically unobservable**
for the stronger reason that there is no delivering mechanism at all: nothing
writes the table mid-round, so there is nothing for order to observe. Execution
order is fixed all the same — [executor](#g-executor) component order, then
declaration order within a component — to keep the [§13.4][s13-4] cursor and the
diagnostics stream deterministic. No trajectory depends on it. The natural
single-pass executor is thereby exactly correct: it builds each handler's bundle
at dispatch, from the live table. No pre-materialization, no staging pass, no
carrier — no shadow table, no allocation, trivially.

**The trade, recorded openly: a handler cannot opt into seeing a same-round
foreign transition.** Same-instant sequential coupling across components is a
cascade, one round per link, deterministic. Coupling tighter than that belongs
inside one component, where declaration order gives exact sequencing across
rounds. That is the synchronous-languages position: a micro-step sees the
pre-state, and effects appear at the next micro-step. The cost of serializing
same-component firings is one extra intra-boundary sweep per event so serialized,
which is microseconds on the rare boundary that fires at all. The rejected shapes
are recorded in row 154: the per-event re-decode with a frozen round-start `u`
(rows 16, 100, 152), live-table reads under the canonical execution order, a
table copy per firing round (row 100), and handlers stripped of their own `y`.

**Why a per-event budget.** A re-enabled event fires at its *true* boundary,
against a fresh sweep; the deferral design and the per-round cap are both
rejected (rows 20, 181). Priors stay honest as a consequence: every prior is a
sample actually taken, never a value recorded to make a rule work out.

Termination is then budget-bounded rather than structural. For `E` declared
events a boundary admits at most `firing_budget · E` firings, hence a bounded
number of rounds, deterministically and independently of pace. A livelock — two
FSMs toggling each other — does not resolve silently. Each toggler spends its
budget and warns (below), and the run proceeds and replays identically:
degradation, not an error, per the doctrine ([§8.4][s8-4]). The warning names the
actual chatterer, while every other event's iteration continues untouched.

The trade is recorded openly. Termination being budget-bounded rather than
structural, the arbitrary-K objection lives on in `firing_budget`. What that buys
is the absence of deferral machinery — no manufactured prior, no re-arm flag, no
`EventDeferred` warning, no one-step artifact, and no collapse of several
intra-boundary re-arms into a single later firing.

**Budget exhaustion degrades; it does not throw.** When an event has fired
`firing_budget` times at a boundary, its further edges there are **lost** for the
remainder of that boundary: it is skipped by the eligibility test while every
other event iterates normally. Exhaustion emits a `FiringBudget` warning
([§13.2][s13-2], [Appendix C][sC]), at most once per event per boundary. The
warning carries the component path, the event name, the boundary time and the
exhausted budget beside the boundary's firing count.

The default of **4** is chosen the way [§8.4][s8-4]'s is: a legitimate re-enable
is one or two firings deep, a toggling FSM pair chatters without bound, and 4
separates them without ever binding on a healthy model. Like every other
degradation here it is a function of the trajectory alone, so the run replays
identically.

**Ticks stay outside the iteration, after quiescence.** The two possible
couplings resolve asymmetrically:

- *Events → ticks: handled by machinery already in place.* Due discrete
  components' output stages (`h_x`/`h_xu`) are gated *into* the boundary sweep
  against a due set fixed for the whole iteration ([§8.5][s8-5]). Every iteration
  round therefore refreshes them for free, against the same `x` — their `g` has
  not run — and post-transition inputs. At quiescence, their published outputs
  reflect the settled boundary instant, which is exactly what "sampling at t"
  should mean for a logically-instantaneous cascade. Earlier rounds' tentative
  values are internal scratch, like RK stage evaluations. [§8.3][s8-3] extends
  naturally: external readers observe the table only after the boundary sequence
  *completes*.
- *Ticks → events: structurally impossible.* A tick's output stages contribute
  nothing guards have not already seen, since they run inside the sweep, from
  current `x`. Its `g` update writes `x⁺` after the sweep, and `x⁺` is first
  decoded at the owner's *next* tick, so it is invisible to every reader within
  the boundary — the standard one-sample `z⁻¹` delay of sampled-data control,
  here enforced by construction. Nothing that happens after quiescence can flip a
  guard, so no combined event/tick fixed point exists to iterate.

The boundary macro-sequence, final form (boundary zero — initialization — is
the same sequence with an empty integrate, [§14.5][s14-5]):

> integrate → project → **[sweep → guards → handlers]** iterated to quiescence
> (under the [firing budget](#g-firing-budget)) → all due `g` updates → logging / I/O staging.

The sequence decides the mixed case: a
[continuous component](#g-continuous-component)'s handler and its discrete
observers' ticks landing on one boundary. Take an engine's `starting → running`
transition under a 50 Hz FCS. The transition fires in the iteration segment, the
re-sweep re-runs the FCS's stages against `running`-mode ports, and its `g` then
updates from post-transition values.

### 8.7 Real-time pacing

**The invariant: [pacing](#g-pacing) is outside the semantics.** The pacer inserts waits between
completed [frames](#g-frame) and never reorders, skips or alters the [boundary](#g-boundary) sequence. A
paced and an unpaced run with identical input [traces](#g-trace) produce bit-identical
trajectories — deterministic [replay](#g-replay) ([§2.2][s2-2]) extends over pace. Interactive runs differ
only because their *inputs* differ. Detection policy is inside the semantics:
event localization runs identically paced or unpaced, its [sweep](#g-sweep) cost
absorbed as debt like any other expensive frame ([§8.4][s8-4]); degrading to
boundary detection under pacing was rejected (row 80).

**Wall-clock mapping: piecewise affine, re-anchored at every knee.** The map is
$\tau(t) = \tau_{\mathrm{anchor}} + (t - t_{\mathrm{anchor}})/p$, with the anchor pair as its reference point. A
live pace change re-establishes the anchor at the current `(t, τ)` so the new slope
applies only forward (row 21); un-pause re-anchors for the same reason. Debt is
cleared at re-anchor: a deliberate user action is a natural sync point, and the
counters record what was forgiven.

**Deadline law: absolute [schedule](#g-schedule) with bounded [debt](#g-pacing)** (row 21). Frame deadlines come from
the map; a frame exceeding its wall budget `h/p` leaves debt that subsequent frames
repay by running short or waitless — the long-run rate is exact and ms-scale hiccups
(GC, scheduler) are invisible. Debt beyond a threshold — **five frames' worth of
budget, `5·h/p`** — is forgiven by re-anchor plus warning, so long stalls
(debugger, laptop sleep) do not trigger catch-up bursts. Five frames' worth sits
comfortably above the ms-scale hiccups debt exists to absorb silently, and far
below the seconds-to-minutes stalls forgiveness exists for, so neither case lands
near the threshold.

**`p = ∞` is pacer-off, not a limit value** (row 21). Unpaced mode is the explicit
*absence* of deadlines — no waits, no debt, no warnings; by the invariant, the same
execution with the waits deleted.

**Wait mechanism: hybrid sleep-then-spin, one knob.** Non-realtime OSes guarantee
only a lower bound on sleep: the thread becomes runnable no earlier than requested;
the wake-up is best-effort (timer granularity, scheduler load, macOS timer
coalescing), with no hard upper bound. Measured on the dev machine, idle,
against 2 ms requests (2026-07): Julia `sleep` overshoots by ≈ 1.4 ms median, and
`Libc.systemsleep` by ≈ 0.5 ms. Behind the `sleep` figure are libuv's
millisecond-granularity timers; sub-ms requests are accepted and rounded up.
Spikes under load are unbounded. The pacer therefore sleeps toward
`deadline − margin` and spins the remainder:

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
task-yielding `sleep` vs. thread-blocking `Libc.systemsleep` — is settled in [§10.2][s10-2]:
the coarse phase uses task-yielding `sleep`, with `margin` absorbing its overshoot.

**Diagnostics.** Overrun count, current and peak debt, forgiven-debt events and wait
statistics are published as [framework status](#g-framework-status) (the frozen
diagnostics value each snapshot carries beside the table) for GUI and logs.
Today's `SimControl` fields are the precedent.

**Forward pointers.** The wait interval is the natural staging [slot](#g-slot) for externally
injected inputs, applied at the next boundary; the staging rules — and the
concurrency model generally, which [§8.3][s8-3] constrains but does not decide — belong to
[§9][s9], [§10][s10].

---

## 9. Runtime periphery: the data plane

This chapter covers GUI, input [devices](#g-device), network I/O and logging, and how
data crosses between them and the [§8][s8] loop:

- the architecture that replaces the shared mutable model ([§9.1][s9-1]);
- outbound [snapshot](#g-snapshot) publication ([§9.2][s9-2]);
- the inbound path:
  - root input [slots](#g-slot), [claims](#g-claim) and the [frozen roster](#g-roster) ([§9.3][s9-3]);
  - per-device staging and the [drain](#g-drain) ([§9.4][s9-4]);
  - the input [trace](#g-trace) ([§9.5][s9-5]);
- the device authoring [contract](#g-contract) ([§9.6][s9-6]);
- the GUI write path ([§9.7][s9-7]);
- the third cross-task channel, runtime diagnostics and liveness ([§9.8][s9-8]).

The machinery that drives the loop itself follows in [§10][s10].

### 9.1 No shared mutable model: staged writes, snapshot reads

Everything in the [periphery](#g-periphery) — GUI, input [devices](#g-device), network I/O, logging —
runs on tasks the loop does not control, and every one of them wants the data
the loop is stepping. FlightCore's answer is one big lock: `SimControl` and the
live `Model`, guarded by `io_lock`, with one task per attached interface
reading or mutating the model under it (sim.jl). That lock does enforce the
[boundary](#g-boundary)-visibility rule ([§8.3][s8-3]), being only ever free between steps.
Transplanting it here was nevertheless rejected, on three structural costs
(row 22). One of the three is load-bearing for everything below: under a lock,
input timing is scheduler-determined and unrecorded. There is then no defined
input [trace](#g-trace), and bit-identical [replay](#g-replay) ([§8.7][s8-7]) is unachievable *in principle*
for interactive runs.

The replacement has five planes. Its unit of account is the [frame](#g-frame) (one
iteration of the loop: [drain](#g-drain), integrate, boundary sequence, publication), the
grid step [§8.4][s8-4] names. A frame is what `step!` counts, and the ordinal that
keys the trace; "per frame" means this throughout the document. The kinematic
*reference frames* of the aircraft domain are a different word, and always
appear compounded: the b frame, the ECEF frame.

1. **Staging (inbound):** devices submit pending input writes at any wall-clock
   moment, never touching live [slots](#g-slot) ([§9.4][s9-4]).
2. **The drain:** exactly one point in each frame, at its top, where the loop
   takes the staged batches and applies them to the root input slots. Never at
   a `t*` boundary ([§8.4][s8-4]). Between drains the loop owns its data
   exclusively, and no lock is held during stepping, ever.
3. **Publication (outbound):** at the end of each boundary sequence the loop
   publishes an immutable [snapshot](#g-snapshot). Readers observe it without coordinating
   with the loop ([§9.2][s9-2]).
4. **Control:** pause, pace and stop on a separate few-word atomic surface
   ([§10.1][s10-1]).
5. **Task topology:** one loop, one task per rostered device except the
   calling-task device, all run-scoped.

The fifth plane carries the most machinery, and the rest of this section spells
it out.

**Device tasks are run-scoped.** `run!` spawns one task per other [roster](#g-roster) entry
after device `init!`, and [§10.4][s10-4] joins them all at every stop ([§10.6][s10-6]).
`attach!` never spawns: it registers, in a stopped-sim state only ([§9.3][s9-3]), and
the task appears at the next `run!`.

**The calling-task device is pinned; the loop is the movable piece.**
Calling-task affinity is a device trait — `needs_calling_task`, default `false`
([§9.6][s9-6]) — held by at most one device per roster (the admission checks,
[§9.3][s9-3]). The shipped GUI declares it, because CImGui ties rendering to the
calling (main) task. With such a device rostered, the loop moves to a spawned
task for the duration of the run, and the [calling task](#g-calling-task) (the task that invoked
`run!`) runs that device's loop body. It runs it inline, inside the same
[§9.6][s9-6] wrapper that brackets any spawned device's loop body. With no such
device rostered, the loop runs on the calling task.

| | no calling-task device | calling-task device rostered |
|---|---|---|
| the calling task runs | the loop | that device's loop body, inline |
| spawned tasks | one per rostered device | the loop, plus one per other rostered device |

The loop-on-the-calling-task case is the unattended register. It is what the synchronous
rethrow ([§13.4][s13-4]) presupposes, and what lets parallel unattended [sweeps](#g-sweep) thread
`run!` inline with no nested task fan-out: one immutable [`Build`](#g-build) is shared
across the workers ([§12.2][s12-2]), and each `Simulation` owns its own [buffers](#g-buffer).
Pre-materializing the sweep's [activations](#g-activation) — its per-eltype
[executable sets](#g-executable-set), `build(world; activations = …)` ([§12.4][s12-4]) — then leaves no
worker synchronizing on anything. Either way `run!` blocks its caller until the
run ends; what varies is what the calling task spends the run doing.

**Topology is derived after initialization**, from the [frozen roster](#g-roster) plus the
outcomes of device `init!`, and never from `run!`'s keywords. [§10.4][s10-4] carries
the rule and the case that motivates it: a calling-task holder whose `init!`
failed returns the loop to the calling task.

**Spawn-inside-`run!` is the start gate.** A task exists only once the run it
serves exists. Any first-boundary synchronization a device needs is the
counter-plus-condition [predicate](#g-predicate) wait ([§10.3][s10-3]), never an `Event` latch
(row 93).

Two rules bind the implementation:

- **Every handoff is one atomic reference operation.** Both shared structures
  reduce their mutable surface to single words: release/acquire `@atomic`
  fields. The GC is the reclamation mechanism, since an object a reader still
  holds is reachable and therefore never recycled. That dissolves the
  reclamation problem (hazard pointers, epochs, RCU grace periods) which makes
  these patterns hard in non-GC languages. Deep immutability of the exchanged
  objects is what makes this sound.
- **No user code, no unbounded work, ever, inside a framework critical section.**
  FlightCore's pathologies are all "arbitrary code under a shared lock". Here
  mappings run on device tasks, rendering runs against snapshots, and the loop's
  frame contains framework and model code only. A stalled device produces stale
  snapshots and late staging; it cannot stall the loop. Two exposures are
  residual and pre-existing: GC pauses and OS scheduling, for which the pacer's
  debt absorption is the mitigation.

One consequence is recorded here because it collapses an API axis: interactive
and unattended simulation stop being different execution modes. An
[unattended run](#g-unattended-run) is the same loop with empty staging and no snapshot readers. A
replayed interactive session is the same loop with staging fed from a recording
([§9.5][s9-5], [§10.7][s10-7]).

### 9.2 Outbound: snapshot publication

Every consumer of model state — GUI panels, output [devices](#g-device), the log — needs a
picture of the model that holds still while it is read, and none of them should
have to make the loop wait to get one. Publication is the answer: one immutable
value per [boundary](#g-boundary), handed out through a single atomic reference. Three
things follow from it, in that order: how a [snapshot](#g-snapshot) is built and published,
how the log retains published snapshots under a bound, and how an output device
addresses what it reads.

#### The snapshot and its publication

**Rule.** The loop builds each snapshot in private memory, then publishes it with
a single release-store to an [`@atomic latest`](#g-latest) reference. A snapshot carries
the boundary-consistent [signal table](#g-signal-table), `t` and the framework status. Readers
acquire-load that reference and then work with an immutable, coherent world for
as long as they like. The [calling task](#g-calling-task) (the task that invoked `run!`) reads the
same value through `latest(sim)` — the inspection register ([§10.6][s10-6]).

The exchange is wait-free in both directions: a wedged reader cannot delay
publication by a nanosecond, and the loop cannot tear a reader's view.
Publication happens only after the boundary sequence completes ([§8.3][s8-3] as
extended by [§8.6][s8-6]).

**Binding rule: nothing reachable from a published snapshot is ever written
again.** The table's immutable values ([§4.1][s4-1], [§7][s7]) make the compiler enforce most
of it. **Why.** The soundness of lock-free reading rests on this rule.

**The [framework status](#g-framework-status) is a concrete frozen value, not a window onto live
bookkeeping.** It carries the pacer diagnostics ([§8.7][s8-7]), plus the per-writer
diagnostic batches, suppressed and cumulative counters and liveness timestamps
the loop takes at frame top ([§9.8][s9-8]). **Why.** The binding rule forces that
shape: a status referencing an accumulator its writers are still filling would
be a snapshot whose contents change after publication.

**The captured table is the whole table** — declared [ports](#g-port) and auto-published
fields, every one of them public ([§11.3][s11-3]). No presentation layer has anything to
filter. Private intermediates are not in it, never having been [cells](#g-cell) at all
([§5.2][s5-2]). The inspection path for one is **promotion to a declared output**: a
line in `output_types`, and the value appears in the snapshot, the log, the GUI
and the wiring alike. Its visibility is then an authored fact like every other.

**The captured table also includes the root [slots](#g-slot)** ([§15.4][s15-4]). Slots are source
cells of the table, not state stores, so they ride along. That is load-bearing,
not incidental. The [§9.7][s9-7] [peek](#g-peek) (showing a widget's own pending write, else the
snapshot value) falls back to the snapshot, and that fallback is what an idle
live widget displays. Read-only mirrors of claimed slots — the axis sliders
under joystick [claim](#g-claim) — show the applied slot value from the snapshot. Slot
values in the log are derived data, recomputable from the [trace](#g-trace), which is
consistent: snapshots are derived wholesale.

**The snapshot deliberately does *not* carry the state stores (`x`, `m`).** Two
reasons. The state trajectory is *derived* data, recomputable from the
[trace header](#g-trace-header) plus the batches ([§9.5][s9-5]) by bit-identical [replay](#g-replay). And
per-boundary capture would systematically record derived data — the same
asymmetry the trace-default decision (row 29) refuses in the other direction.

"What was the private state at t = 37.2?" is answered by replaying to 37.2 and
inspecting the live stores. A state field wanted in logs or GUI has the honest
remedy of being declared public, at a cost of one auto-published cell per
[sweep](#g-sweep). Post-run continuation reads the live stores directly. Periodic full-state checkpoints — warm restart without
replay-from-zero — are a [guarded addition](#g-guarded-addition) shaped as an opt-in log policy, and
a dev-mode flag auto-publishing all state fields is a possible future
diagnostic.

#### The log: retained snapshots under a bound

**Logging dissolves into publication.** The log is a vector of retained snapshot
references: the same objects, zero extra copies. The per-step `deepcopy` detour
of the `SavingCallback` disappears. The cost is one snapshot allocation per
boundary, on the framework side of the [§7.5][s7-5] scope, which already carved out
logging. Logged snapshots are not garbage at all; unlogged ones die young.
Preallocated snapshot [buffers](#g-buffer) are rejected (row 23).

**Retention: the trace's kill switch, plus [decimation](#g-decimation).** The log takes the
same plain on/off switch the trace has ([§9.5][s9-5]), and additionally a
keep-every-kth retention policy: `log_every` ([Appendix B][sB]). **Why decimation is
admissible here and not there.** The derived/primary split (row 38): the log is
recomputable from the trace by replay, so a thinned log costs resolution in a
*view*, never a record. Thinning the trace would destroy the only primary
account of a session (row 29). Decimation is a retention policy only — every
boundary still runs, is still published to live readers, and still enters the
trace.

Decimation slows the log's growth; it does not stop it. The default
configuration is `log = true, log_every = 1` alongside `t_end = Inf`, the honest
interactive default ([Appendix B][sB]), and it grows for as long as the session lasts.
At C172X scale and 50 Hz that is gigabytes per hour, and it ends in an
out-of-memory nobody was warned about.

**Rule.** The log takes a retention bound beside its switch and its stride:
`log_max`, the maximum number of retained snapshot references, default **65536**
(2¹⁶), with `Inf` the explicit opt-out.

**A count, not a memory budget.** Snapshots are immutable object graphs with
internal sharing ([§4.1][s4-1]), and a [field handle](#g-field-handle) (an immutable query object
consumers evaluate at their own arguments, [§4.4][s4-4]) rides as a reference to
build-time-frozen data ([§7.5][s7-5]). Byte accounting over them is therefore fuzzy and
platform-dependent. A count is exact, and converts to memory through one number
the user can measure once: `Base.summarysize` of a single snapshot.

**The default is finite unconditionally**, not a modal rule keyed on `t_end`. A
finite run shorter than the bound never notices the bound, and one number is
easier to hold than two regimes. At 50 Hz and full density, 2¹⁶ boundaries is
about 22 minutes before anything is dropped at all.

**When the log fills, the retention stride doubles: coverage stays global.** A
rolling window — recent past at full density, the start of the session forgotten
— was rejected (row 137). Instead the log **re-decimates progressively**: after
*k* generations the effective stride is `log_every · 2^k`. The whole run stays
plottable, and what coarsens is density, never extent.

**Why.** That is the division of labor (row 38) carried through. The log's chief
consumer is the post-run plot of a session *as a whole*, and nobody plots hours
at 50 Hz. Full density over any *segment* of interest is what replay from the
trace recovers ([§9.5][s9-5], [§10.7][s10-7]).

**Normative are the guarantees, not the mechanism.** Three of them: the bound is
respected continuously, the retained count never exceeding `log_max`; coverage
is global at the effective stride; the endpoints below are kept.

Mechanism sketch, non-binding: the thinning is **amortized**. When the log
fills, the stride doubles immediately, and each subsequent retained append also
releases one predecessor of the previous generation — a cursor over the odd
indices, O(1) amortized, with physical compaction once per generation. A
generation's thinning therefore completes exactly when its refill does.

> fill at stride `log_every` → **full** → stride doubles → thin one predecessor
> per retained append, while refilling → **full** at stride `log_every · 2` → ⋯
> → generation *k* at `log_every · 2^k`

Amortizing rather than halving in one shot is a responsiveness choice (row 137).
The amortized form drops exactly one old snapshot per retained append, the same
steady trickle a rolling window would produce, so keeping coverage global costs
nothing extra in GC pressure. The loop's own work is pointer bookkeeping,
microseconds either way and on the framework side of the scope ([§7.5][s7-5]).
Publication stays wait-free and readers never block: a reader holding a released
snapshot simply keeps it alive.

**The endpoints are retained unconditionally.** The boundary-zero snapshot
([§14.5][s14-5]) and the terminal snapshot ([§10.4][s10-4], [§13.5][s13-5]) survive any `log_every` and
any `log_max`, and do not count against the bound — two extra references. The
terminal snapshot's status carries the run's final cumulative diagnostic
counters ([§9.8][s9-8]). A run's two endpoints and its complete diagnostic account
therefore always outlive whatever retention did to the middle.

Two compositions are worth stating once. The `totals` monotonicity across logged
snapshots ([§9.8][s9-8]) is untouched: re-decimation, like decimation, loses *which*
boundary within a stretch an occurrence fell on, never *how many*. And `log_max`
is a **view policy, not a trajectory-determining one** — like `log` and
`log_every` it stays out of the trace header's deployment block, and replay
neither records nor compares it ([§9.5][s9-5], [§10.7][s10-7]). Sizing follows: the `sizehint!`
for the expected duration ([§7.5][s7-5]) is now naturally capped by `log_max`, which is
also what defines the hint when `t_end = Inf`.

#### Output-device bindings

**Output-device bindings are snapshot bindings.** An output device — telemetry,
the XPlane visualizer, disk streaming — consumes snapshots via [§10.3][s10-3]. It
addresses what it reads with the [selectors](#g-selector) (the closed family of deferred reads
resolving against a source, [§14.4][s14-4]), which reach any cell, the diagnostic
register admitting deep paths. A binding is resolved at attach against the
`Build` with [did-you-mean](#g-did-you-mean) (the offending name plus the list-in-hand it should
have matched), and compiled to one gather — the output half of the binding
interface ([§9.6][s9-6]). `map_output` therefore receives a labeled NamedTuple instead
of performing its own path lookups. That discharges the obligation stated in
[§15.4][s15-4]: a substitution that breaks a binding fails at attach, not with silent
garbage UDP.

This is **diagnostic observation** ([§13.5][s13-5]): human-facing, with no effect on run
semantics. It is the same register as the log retaining the full table and the
GUI's deep-reading panels. Every cell is reachable, the table being public
throughout ([§11.3][s11-3]), and an intermediate a device wants to stream is one
promoted to a declared output.

**A binding chooses its register.** A deep path is the *inspection* register:
zero promises, free access, right for looking at *this* build. An exported
output [face](#g-face) — spelled `get_face(name)` ([§14.4][s14-4]) — is the *integration* register:
named, curated, meaning-stable under substitution, right for consumers that
outlive the build they were configured against. What makes a face meaning-stable
is writer-independent semantics ([§15.4][s15-4]).

**Why the choice matters.** Attach validation converts *structural* drift to
loud errors in both registers. Only faces protect against *semantic* drift — a
substituted aircraft publishing the same path at the same type with a different
meaning, a CG velocity under a name read as body-origin velocity. Nothing else
can: meaning is not in the schema.

Semantically generic consumers should therefore bind faces: a visualizer needs
pose, and every aircraft has one. Aircraft families should export the
conventional surface such consumers need, a library/migration deliverable
([§16][s16]). Wrapper types make face semantics structurally checkable, as in
`VelocityData` with its `v_eb_b` defined *at the type* as body-origin velocity:
a bare vector doesn't wire, and wrapping the wrong quantity is a deliberate lie,
not a drift.

### 9.3 Inbound: root input slots, claims and the frozen roster

**The [write surface](#g-write-surface) (the set of faces a writer's batch entries may reach)
is root input [slots](#g-slot).** A root slot *is* the root [assembly](#g-assembly)'s own input
[face](#g-face), declared through `input_connections` ([§11.6][s11-6]): routed inward to
consumers, produced by no [component](#g-component). At every non-root level an input face
is fed by the parent's wire, and at the root there is no parent. No dedicated
vocabulary is needed.

A root slot is usefully read as the output face of the one producer the build
never sees: the [periphery](#g-periphery) and the services. On that reading, slot
exclusivity (below) is that producer's one-writer right, and the totality
([§14.6][s14-6]) is its completeness obligation.

Slots are sources to the build-time scheduler, constants within a frame, and
the *only* thing the periphery may write. The GUI reaches them through the
resolution ([§9.7][s9-7]), and control commands are not writes ([§10.1][s10-1]).
[Devices](#g-device), mappings, the [trace](#g-trace) and the GUI write path address slots by
**face name** ([§11.6][s11-6]). Structural slash paths never cross the periphery's
*write* [boundary](#g-boundary): the write side speaks the root [contract](#g-contract)'s names only.
The read side chooses per binding — slash paths in the inspection register,
face names in the integration register and in load-bearing service reads
([§9.2][s9-2]/[§13.5][s13-5]/[§14.4][s14-4]).

**Slot exclusivity: one writer per slot at any time** ([§15.4][s15-4]). A
device [claims](#g-claim) its slots at attach, and claiming an already-claimed slot is
an attach-time error. Detaching releases the claims: a released slot's GUI
widgets are live again from the next run ([§9.7][s9-7]). Exclusivity replaces any
cross-device conflict *policy* — attachment-order precedence at
[drain](#g-drain), say (row 44). Per-device [cells](#g-cell), the CAS merge and the
atomicswap drain all stay; they serve atomicity and [coalescing](#g-coalescing), not
arbitration.

**A claim is what a device *may* write, not what it will.** Data-dependent
write-sets are ordinary: a UDP/JSON peer writes whichever subset of faces the
incoming message names, and `map_input` is arbitrary user code the framework
never inspects. Such a device therefore claims the **binding's enumerated
allowed set** — the faces the binding table lists, whether or not any given
batch touches them. The claim is registered at attach exactly as a joystick's
is. A broad claim costs liveness: every enumerated face is claimed for the
device's whole attachment. The derived-liveness rule ([§9.7][s9-7]) therefore renders
the device's GUI widget read-only even on faces the peer never writes. Narrow
the binding to narrow the claim — the enumeration *is* the interface.

**Every writer has a write surface, and the periphery enforces it.** A batch
entry reaches a slot **iff the named face is inside the writer's surface**;
anything else is discarded with a runtime warning ([§13.2][s13-2]). Because surfaces
are static per run (the [roster](#g-roster) freeze, below), enforcement runs entirely at
*staging* — the earliest site, on the writer's own task. The drain performs no
checks at all. **Every device's surface is its claim set**, and a claim set has
two *sources*:

- **Returned** — the binding enumerates the faces: it declares
  `is_input(b) = true`, `claims(b)` ([§9.6][s9-6]) is
  called once at attach, and what it names is staked. Such a claim is static
  for the attachment and exclusively its own, since claims are disjoint by
  construction. It is binding-bounded even where no one else is involved: a
  mapping that has drifted onto an unenumerated face is a diagnosable anomaly
  (`OutOfClaimEntry`), never a silent write, claimed or not.
- **Computed** — the binding declares `is_greedy(b) = true` ([§9.6][s9-6]) and the
  framework computes the claim at attach: all root input faces minus the
  union of the rostered claims, the unclaimed complement at that instant.
  This is the shipped GUI's claim ([§9.7][s9-7]) — everything unclaimed, without
  configuration. It is disjoint from every incumbent claim by construction, so
  exclusivity validates trivially and nothing downstream can tell the two
  sources apart.

One claim mechanism, two claim sources. The source is exhausted at the attach
point: past it, validation, roster-entry storage, shape compilation
([§9.4][s9-4]), the drain, the trace and detach-releases-claims treat a computed
claim exactly as a returned one. The GUI is therefore not an exception but an
ordinary enumerated writer whose enumeration the framework performed.
Opportunistic writing by autonomous devices does not exist: a device that
wants a face enumerates it, and greediness is an explicit declaration, never a
default. Cross-writer races on one slot therefore cannot arise structurally:
every claim is exclusive, whatever its source. That is what keeps drain order
a diagnostic fact (below) and lets a drained GUI value simply stay
([§9.7][s9-7]).

**One framework-owned remainder: the [harness register](#g-harness-register).** Beside the roster
sits a **task-free entry point**, `stage!(sim, "face" => value, ...)`, the
harness/REPL write path ([§10.6][s10-6]). It stages a batch from the
[calling task](#g-calling-task) itself (the task that invoked `run!`). Its always-present
cell is drained, traced and surface-checked exactly as any device's. The
register's surface is the one thing in the design that is *derived* rather
than claimed: the unclaimed complement, the faces no rostered device speaks
for. That surface is recomputed at every stopped-sim roster change, and is
therefore as fixed within a run as any claim set. A `stage!` write to a
claimed face is rejected at staging (`ClaimedFaceEntry`, naming the
incumbent). The one seam — a batch staged while stopped whose face a
subsequent `attach!` claims — is renormalized away at the attach itself
(below). The [harness cell](#g-harness-cell) (the always-present staging cell of the harness
register) drains **last**, by convention: with every surface disjoint the order
is unobservable, so the rule exists to make the trace read the same way every
time, not to arbitrate anything.

**Slot initial values are owned by the init/trim services** ([§15.4][s15-4]).
Input declarations are bare types ([§11.2][s11-2]) and carry no defaults, yet a
slot unfed by any device must hold a defined value from the first frame.
Today's `U()` constructors provide these (`mixture = 0.5`). Export-entry
defaults were rejected (row 47). `init!` establishes every slot, and the
[trace header](#g-trace-header) captures the result. Totality is enforced pre-write at every
complete-world application — `init!`, trim setup, trim commit
([§14.6][s14-6]).

**The roster is frozen per run: attach and detach are stopped-sim
operations.** `attach!`/`detach!` are legal in the `built`, `initialized`
and `stopped` states ([§10.6][s10-6]) and an error while `running`. That error is
`ServiceLifecycle` ([Appendix C][sC]), the same kind that gates the [§14][s14]
services. The prohibition includes pause: pause is a control-plane state
*inside* a run ([§10.1][s10-1]), and a surface that could move while paused would
move mid-run. The roster — entries, claims, attachment order — is therefore a
plain immutable value the loop reads once at `run!`. The partition of the root
face set into per-writer surfaces plus the harness remainder is a static,
inspectable fact of the run: printable before it starts, valid until it ends
(the provenance register, [§13.7][s13-7]). No republication machinery exists —
no atomic roster reference, no per-frame acquire-load, no next-frame
attachment granularity, no sequence numbers. Attachment order is the roster's
own order. The trace still tags entries with a stable device id, never a
roster index, because ids read across runs, where the roster does change.
Attach validation, claim registration and the staging-shape compilation
(below) all run at the attach point, which makes `attach!`/`detach!`
stopped-sim configuration operations beside `init!` and trim ([§14][s14]). While
a simulation runs, its configuration — build, roster, claims, surfaces — is
immutable. The doctrine ([§10.5][s10-5]) extends to its final form: the running
periphery stages writes and issues control commands, *and nothing else
changes*.

**Device identity, ids and roster admission.** Identity is the device
instance: the same object (`===`) may occupy at most one roster entry. Two
instances of the same type — two joysticks — are two devices. The stable
device id the trace, heartbeat and diagnostics speak is assigned at `attach!`,
monotonic per `Simulation` and never reused. It lives exactly as long as the
entry: across runs (roster persistence, [§10.6][s10-6]), until `detach!`.

Admission is a three-part check at the attach point, in order:

```
# each line: the condition that rejects the attach → the diagnostic it raises
identity   this instance is already rostered
               → AlreadyAttached      (names the entry and its binding)
affinity   this device declares needs_calling_task, and a rostered
           device already declares it
               → CallerTaskConflict   (names both devices)
claims     face exclusivity: this device's claim set meets a rostered claim
               → ClaimConflict        (names two distinct devices)
```

An already-rostered instance is rejected rather than silently absorbed because
rebinding has an explicit spelling: `detach!` then `attach!`, both legal at
any stopped-sim point. Either a silent no-op or a silent rebind would discard
a binding the caller handed over. The affinity check admits at most one
rostered device declaring `needs_calling_task`, because the topology
([§9.1][s9-1]) makes the calling task a single-slot resource. Running the claims
check after the identity check is what makes `ClaimConflict` always name two
*distinct* devices, never a device colliding with its own earlier attachment.

**Device death does not detach.** A mid-run crash, voluntary exit or unplug
([§9.6][s9-6], [§10.4][s10-4]) ends the device's *task*: the cell stops filling, the [§10.2][s10-2]
heartbeat shows the death by name, and the roster entry — claims included —
persists to the end of the run. The [orphaned claims](#g-orphaned-claims) are the accepted cost of
[the freeze](#g-the-freeze): the device's slots hold their last-drained values and no other
writer inherits them. The read-only widgets ([§9.7][s9-7]) render the orphan visibly
("claimed by `T16000M` — task dead"), never mysteriously. Recovery is between
runs: stop, `detach!`, and either `init!` (fresh trajectory) or
`replay!`-to-end then `run!` (continuation from the interrupted boundary,
[§10.7][s10-7]). The death is an anomaly, not a surface event.

One deliberate asymmetry is on record as a **[guarded addition](#g-guarded-addition)** (a
capability the design admits but does not build). A pure reader — a binding
declaring `is_output` alone ([§9.6][s9-6]), a visualizer or a telemetry tap — claims
nothing, so attaching one mid-run would move no writer's surface. A dynamic
reader list would touch only [§10.3][s10-3] wakeups, the heartbeat and the shutdown
join, never the drain, and is cleanly severable from the freeze should the
join-a-running-session workflow find a customer. The [§10.2][s10-2] thread-budget
warning runs once per `run!`, against the frozen population.

### 9.4 Inbound: per-device staging, representation and the drain

A device produces writes on its own task, whenever its hardware or its peer
hands it a datum; the loop consumes them at frame top. Between those two rates
something has to hold each device's pending writes and hand the loop a value
it can apply. That something is the staging cell. This section fixes the
policy under which a device writes into its cell, the shape the cell holds,
and the [drain](#g-drain) (the frame-top swap that publishes staged device
inputs into the root slots) that empties it.

**Rule.** Staging keeps one atomic [cell](#g-staging-cell) per attached
[device](#g-device) under one [coalescing](#g-coalescing) policy: CAS merge,
newest wins per [face](#g-face). Each cell has a single writer, its own device
task, and holds that device's latest pending [batch](#g-batch) of
[slot](#g-slot) writes. Staging merges the incoming batch into the pending
one. Untouched faces survive; re-staged faces take the newest level, the
per-face ZOH. The CAS can fail only because a drain intercepted the old batch,
so the retry is bounded. The failure case is precisely correct as well:
intercepted writes are already applied, and must not be re-staged.

**Why.** Merge is the *only* policy because it is always correct. A
**complete** writer covers every face in every batch it stages — a joystick
delivering its full write-set every poll is the type case. For such a writer
merge and overwrite are provably the same operation, which makes overwrite a
degenerate fast path rather than a second semantics. A **sparse** writer — the
GUI, a JSON peer, each staging only what was touched — loses writes silently
under overwrite instead. [§15.3][s15-3] works that hazard through: a pending
`flaps` edit clobbered by an unrelated `gear` message, undrained and
undiagnosable. A user-facing overwrite opt-in (`complete(binding)`) is closed
(row 104).

**Rule.** The staged representation is fixed per attachment, compiled at
attach. An enumerated writer's [claim](#g-claim) set and slot types are both
known at attach, `claims(binding)` ([§9.6][s9-6]) read against the root
[contract](#g-contract). So the framework fixes the cell's content type there:
a positional tuple over the claim set, `Union{Nothing, T}` per face. Those
unions are isbits, hence pointer-free. `nothing` means *not touched this
time*, never "reset". The levels doctrine is therefore untouched, slots only
ever receive the non-`nothing` positions, and the `Union` never reaches the
model. The face-name → position schema lives in the [roster](#g-roster) entry.

In sketch form:

```julia
# claim set enumerated f₁ … fₙ; Tᵢ = declared type of the slot behind fᵢ

# attach fixes the cell's content type:
Tuple{Union{Nothing,T₁}, …, Union{Nothing,Tₙ}}

# stage!, on the device task: the shim turns the author's face ⇒ value pairs
# into `incoming` — name → position, convert to Tᵢ, fill nothing

# staging merges into the pending batch, newest wins per face:
merged[i] = incoming[i] === nothing ? pending[i] : incoming[i]

# the drain scatters, nothing skipped:
batch[i] === nothing || (the slot at position i receives batch[i])
```

The consequences are each mechanical. The merge is positional
(`incoming[i] === nothing ? pending[i] : incoming[i]`), so it compiles
straight-line and union-split. The drain applies each cell through an
attach-compiled **scatter**: position → slot cell, statically typed, `nothing`
skipped. That scatter is the exact mirror of the compiled output gather
([§9.2][s9-2]).

Authors never build the shape by hand. `map_input` returns face ⇒ value pairs
for whatever the datum touched, and `stage!` normalizes those pairs through an
attach-compiled shim. The shim does three things: name → position, convert to
the slot's declared type, fill `nothing`. It thereby confines the residual
name-shaped dynamism to one framework-owned conversion on the device task, at
the [boundary](#g-boundary) where wire-shaped data becomes system-shaped data.
Author-built total tuples are rejected as a padding form, the same disease
row 74 and the handler return law refuse (row 104).

**A greedy entry needs no special treatment here.** Its claim was computed at
the attach point, and by the time shapes are compiled it is an ordinary claim
set ([§9.3][s9-3]). The GUI's cell is compiled exactly as a joystick's.

**The [harness cell](#g-harness-cell) (the always-present staging cell of the
harness register) gets the same treatment.** Under the roster freeze its
derived surface — the unclaimed complement — is as static as any claim set, so
it too is compiled to a positional shape. That shape is recompiled at each
`attach!`/`detach!`, both stopped-sim points, and it carries the same shim,
merge and scatter. Being always present, the harness cell gets the
compilation unasked. It is also the one cell whose shape the framework derives
rather than receives.

One representation, one mechanism. The name-keyed dynamic path the mutable
surface used to force does not exist, and no face name is ever resolved inside
the loop's frame.

The recompilation has one seam. A pending harness batch staged *before* a
stopped-sim `attach!` may hold the old shape, or may name a face the new claim
covers. The attach renormalizes that batch: it is reshaped, and newly-claimed
faces are discarded with `ClaimedFaceEntry`. So the run always starts with
cells matching the run's schemas.

**Rule.** Diagnostic sites follow the compilation, all of them to staging.
Face-name validity, surface membership and value convertibility are all static
facts of the run. Every check therefore runs in `stage!`'s normalization, on
the writer's own task. A device's out-of-claim face has no position in the
schema and is rejected with `OutOfClaimEntry`. Staging is an earlier,
better-attributed site than the drain for that rejection, and the kind and
[payload](#g-payload) are the same either way. The GUI is included, its claim
being an ordinary one. A **harness** write to a claimed face is rejected the
same way, with `ClaimedFaceEntry` naming the incumbent device. And a value
that cannot convert to its slot's declared type is discarded with
`EntryTypeMismatch` ([Appendix C][sC]), at the same spot.

Nothing remains at the drain. With surfaces frozen for the run, there is no
fact only the drain can know, and the drain is pure application.

**Doctrine: staged values are levels, never deltas** (`press_count = 17`,
never `presses += 1`). Levels are idempotent and survive coalescing. Button
edges ride as monotonic counters.

**Rule.** At frame top the drain takes each device's cell with one
`atomicswap(cell, nothing)`. The swap is an indivisible take, so there is no
lost-write window. Each taken cell is then applied through its compiled
scatter, **in attachment order**. Attachment order is retained as a
deterministic application order. Under slot exclusivity, cross-device writes
to one slot cannot arise, so the order matters only for diagnostics.

Which *frame* a write lands in remains wall-clock reality. What the drain
guarantees is that the frame's outcome is a pure function of the drained
batches.

Because the roster is a fixed value at `run!`, the drain is fully compilable.
The cells and their scatters form a heterogeneous but *known* tuple the frame
function can specialize on, which means zero dynamic dispatch at frame top.
That is the same per-configuration compile trade the [executor](#g-executor)
([§12.7][s12-7]) already makes, now incurred only at stopped-sim attach
points. The specialization is an implementation freedom
[the freeze](#g-the-freeze) creates, not an obligation. Iterating a roster
array costs a handful of dispatches per frame and remains acceptable.

Two shapes were rejected, both torture-tested in [§15.3][s15-3]: per-slot
atomic cells, and a shared lock-free [batch](#g-batch) stack (row 24).

**Mappings run on the device task.** Today's
`assign_input!(mdl, mapping, data)` becomes the pure
`map_input(data, mapping) → batch`. User-extensible code thereby never
executes inside the loop's frame, and the trace consists of slot-level
batches.

**Mappings are binding data, not shaping code** ([§15.4][s15-4]). A mapping is
a declarative table: axis/button → slot path, plus per-axis conditioning
parameters (deadzone, expo strength). The shipped `TableBinding` applies those
parameters in its generic `map_input`, on the device task ([§9.6][s9-6] — the
shared pure helper, with an owner).

The boundary is set by the face contract: **a face's meaning is
writer-independent**. Faces therefore carry *post-conditioning* semantics. A
GUI slider or a script writes the same command a curved stick delivers, and
running a mouse drag through a deadzone would be absurd. This GUI-parity test
is what places conditioning upstream.

Aircraft-semantic derivation must *not* ride along, the C172X
`q_ref = q_sf · axis` fan-out being the case in point. It is FCS design and
lives in-model, in the avionics. Alternatively it is accepted as a small
per-aircraft×device mapping entry, an aircraft-design fork ([§15.4][s15-4]).

The [trace records](#g-trace-record) post-conditioning levels. Those are
exactly what the model consumed, so [replay](#g-replay) is exact. Raw-stick
provenance — re-running a session through *different* curves — is the known,
accepted loss. Edge logic follows the levels doctrine: devices stage monotonic
press counters. Accumulators (trim offsets, flap detents) are model state, not
mapping state ([§15.4][s15-4]).

### 9.5 Inbound: the input trace

**The input [trace](#g-trace)** is the sequence of drained,
[device](#g-device)-tagged batches per frame. It extends the determinism
([§8.7][s8-7]) end-to-end. Replaying a recorded interactive session reproduces
the trajectory bit-identically, with staging fed from the recording and no
devices or mappings present.

**One record format: every batch is retained sparse.** At the
[drain](#g-drain) — the frame-top swap that publishes staged device inputs
into the root slots — each drained [cell](#g-staging-cell) is scanned. Its
non-`nothing` entries are recorded as (position ⇒ value) pairs, against the
writer's [face](#g-face)-name → position schema in the header (below). That is
an O(surface-width) scan and one small allocation per drained batch.

**Why.** The rule is uniform because a [claim](#g-claim)'s *width* is a fact
about one binding, not about a class of writers. A
[greedy claim](#g-greedy-claim) is enumerated and as wide as the root
[contract](#g-contract) ([§9.3][s9-3]), so keying retention by claim source is
rejected (row 176). Every consumer then handles one format instead of two:
one record format at the trace's edge, no per-entry format flag, one decoder
in the [what-if register](#g-what-if-register) (replay with edited inputs), in
disk serialization and in human inspection, and one inverse conversion in
[replay](#g-replay). That work is paid once, up front, off the loop
([§10.7][s10-7]). The conversion site is the drain and not the staging shim
because the drained tuple is the *coalesced* truth: a shim-side sparse log
would need its own merge.

**The costs are recorded rather than argued away.** On the wide writers the
conversion is what keeps the trace honest. A tuple as wide as the unclaimed
surface carrying one edit would otherwise make trace size track surface
width rather than information; at hundreds of faces, render-rate dragging
inflates the trace past the two-orders-below-the-log budget that justifies
trace-on-by-default (row 29). On the dense [component](#g-component) it costs
**about 2×** — a position beside every value where the positional tuple
carried the value alone. That changes no order of magnitude, and it leaves
the budget (row 29) standing for every writer at once. The allocation is
in-class with what the retention carve-out ([§7.5][s7-5]) already admits, and
per [boundary](#g-boundary) it is smaller than the log's [snapshot](#g-snapshot),
the carve-out's standing occupant — the one qualified exception to
retains-what-was-already-allocated. And the decision is **reversible as pure
implementation**: the conversion is lossless in both directions, so verbatim
retention could return as a per-entry storage optimization if a
marathon-session measurement ever asks for it. Such a return would leave the
record semantics, the header and the replay path exactly as they are.

**The [trace header](#g-trace-header) captures the full initial state** `(x, m)` **plus the
initial root-[slot](#g-slot) values** at `init!`. The capture happens **after `apply!`
and the slot writes, before the boundary-zero sequence runs** ([§14.5][s14-5]). Both
halves of that placement are load-bearing:

- the header holds the *resolved* stores and slots as values, never the sparse
  authored overlay — replay must survive edits to declared defaults, the
  primary-data doctrine (row 38);
- and it never holds the post-transition result, since
  [boundary zero](#g-boundary-zero) is re-executed under replay ([§10.7][s10-7]):
  a post-sequence capture would re-fire authored-condition events on top of
  already-latched state.

An unfed `mixture = 0.5` never appears in any batch, so replay is broken
without the slots; the init/trim services
own slot initialization ([§14.6][s14-6]), and the header capture extends
naturally. The header carries two further things:

- **each writer's face-name → position schema** — the run's frozen surface
  partition — since positional records are meaningless without it and replay
  does not reconstruct claims ([§10.7][s10-7]);
- **the run's deployment block**: `t₀`, `Δt_base`, `h`, `n`, the algorithm
  identifier, `localization_tol`, `localization_budget` ([§8.4][s8-4]),
  `firing_budget` ([§8.6][s8-6]) and the effective `t_end`/`stop_on` pair,
  captured at the same instant as the stores.

The trajectory depends on the deployment block exactly as it depends on the
stores. The deployment binding ([§12.1][s12-1]) sits outside the `Build`, and
`t₀` post-dates even deployment ([§14.5][s14-5]). A header without them could
therefore not back the bit-identity claim ([§10.7][s10-7]). This block is also
what the artifact's **run metadata** names ([§13.5][s13-5], [Appendix B][sB]) —
the recorded home of the effective termination pair. The header capture is the
one full-state capture in a normal run, and the other half of what "given the
initial state and the trace, the log is recomputable" requires. Header plus
batches are the *primary* record; everything else, the state trajectory
included, is derived ([§9.2][s9-2]).

**Trace recording is on by default.** The trace is cleared at `init!` and
retrievable after the run, and a plain kill switch covers memory-constrained
marathon sessions. The asymmetry that decides the default is that the trace is
*primary* data and the log *derived*: given the initial state and the trace,
the log is recomputable, which is what bit-identical replay means. An untraced
interactive session, by contrast, is unreproducible, permanently. The cost
supports the default. The trace retains one small sparse record per drained
batch (above), at drain-rate × device-count — tens of MB per hour worst case,
two orders of magnitude below the snapshot log. No sampling, no rolling window
(row 29).

### 9.6 Devices: one authoring contract, no taxonomy

FlightCore's input/output/GUI trichotomy is lock choreography, not modeling:
with no lock, the protocol the taxonomy encoded has no referent (row 25).

#### Every attached device receives the same handle

The handle carries the two primitive capabilities, read and stage, plus
control access (observe running, request shutdown). Read returns the latest
[snapshot](#g-snapshot), optionally waiting for the next
[boundary](#g-boundary) ([§10.3][s10-3]).

**[`should_abort`](#g-should_abort) is an `attach!` keyword**, defaulting to
`false`. It is per-attachment, never a device property: the same joystick is
advisory in one deployment and load-bearing in another. With it clear, a
device's departure is reported and the run continues without it; with it set,
that departure also requests a sim stop ([§10.4][s10-4]). A departure is the
loop body returning, a crash, or a failed `init!`. The shipped GUI attaches
with `should_abort = true`, since closing the window is the interactive
session's natural end, and `gui = true`'s run-scoped attachment states that
value ([§10.6][s10-6], [Appendix B][sB]).

Input-only and output-only devices are degenerate uses, not framework classes.
A bidirectional network peer is *one* device with one socket and one
lifecycle, not two framework devices sharing state. The GUI is an ordinary
device — the paradigm one, using every capability. It has exactly two genuine
peculiarities, neither taxonomic: main-thread affinity (a launch concern) and
read-modify-write widgets ([§9.7][s9-7]).

#### The authoring contract: four functions, one optional, one trait

A [device](#g-device) is a user type subtyping the framework's neutral root:
`MyDevice <: AbstractDevice`. That is one mandatory word, and it costs nothing
— the [periphery](#g-periphery) has no competing hierarchy to inherit from.
What it buys is `attach!`'s dispatch gate below. The framework asks for

```julia
init!(dev)          # per-run resource acquisition — calling task, before spawn (§10.4)
loop(dev, handle)   # the task body: owns its own wait structure
shutdown!(dev)      # per-run resource release — guaranteed on every exit path
unblock!(dev)       # optional hook, default no-op: make a blocked loop return (§10.4)
needs_calling_task(dev)   # optional trait, default false: run the loop body on the
                          # calling task (§9.1's topology; the shipped GUI's CImGui
                          # constraint). At most one holder per roster (§9.3).
```

The framework owns everything around them. The wrapper is the shutdown
protocol ([§10.4][s10-4]) made structural:

```julia
init!(dev)                                   # its own bracket, pre-spawn (§10.4): a throw
                                             # here is shutdown! + DeviceCrash by name
task = Threads.@spawn try
    loop(dev, handle)
catch e
    report!(handle, DeviceCrash(e))          # §10.4(6): sim continues, device absent
finally
    shutdown!(dev)                           # any exit path: OS resources released
    mark_dead!(...)                          # heartbeat only — claims stay, §9.3
end
```

A `needs_calling_task` device runs the identical wrapper *inline* on the
[calling task](#g-calling-task). The invocation site, not the contract, is its
only difference (the topology, [§9.1][s9-1]; the join exclusion,
[§10.4][s10-4]).

**`shutdown!` must tolerate a partially initialized device.** The release
guarantee holds on the one path *outside* this wrapper too. The initialization
step ([§10.4][s10-4]) brackets each `init!`, and a device that threw half-way
through acquisition goes straight back to `shutdown!`, so nothing it did
manage to open is leaked. The obligation that follows is "close only what is
open" — the same defensiveness `shutdown!` already owes the crash path, where a
loop body may die at any point in its own life. `init!` is correspondingly
*not* asked to clean up after itself: the bracket does once, for every device,
what would otherwise be duplicated in each and enforced in none.

One discrimination in that wrapper: **an `InterruptException` is never a
`DeviceCrash`.** Under the spawned-loop topology the calling task is the one
running a device loop body inline — the GUI's ([§9.1][s9-1]). An operator
Ctrl-C therefore raises *there*, inside user code that did nothing wrong. The
wrapper forwards the control-plane stop and lets the body leave through the
ordinary `running(handle)` [predicate](#g-predicate) ([§10.4][s10-4](4)). There
is no crash report for what is not a crash, and no `should_abort`
consultation, a stop being already requested.

#### The author owns the loop body; the framework owns the bracket

One device
[contract](#g-contract) means author-owned loop bodies: a framework-owned hook
loop would have to ask each device what it waits on, which is the rejected
taxonomy resurrected as a trait (row 102). Under the author-owned body, every
wait structure is ordinary user code composed from handle primitives:

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

A bidirectional peer composes both halves itself, with an inner reader task
inside its own domain, rather than forcing a select engine into the
framework. Two idioms are author obligations the framework can only teach
and diagnose, never force ([Appendix A][sA]): loop on `running(handle)`, and make
blocking calls interruptible (`unblock!`, or timeouts). A forgotten
[predicate](#g-predicate) check surfaces as `DeviceJoinTimeout` with the
device's name; a stall surfaces as a stale heartbeat ([§10.2][s10-2]).
Liveness timestamps ride *inside* the handle primitives, which store them in
the device's own [diagnostic cell](#g-diagnostic-cell) ([§9.8][s9-8]), so the
framework observes activity without owning the loop.

**`should_close` dissolves**: a window ✕ or peer EOT is the loop body
returning. The wrapper's exit path releases the device's OS resources, marks
it dead for the heartbeat and consults `should_abort`; [claims](#g-claim) and
[roster](#g-roster) entry persist to run end (the freeze, [§9.3][s9-3]).
[§10.4][s10-4](6) is literally "the task body returned." The GUI implements
the same contract; the framework calls its `loop` inline on the
[calling task](#g-calling-task) instead of spawning (the pinning,
[§9.1][s9-1]).

#### The binding: framework-legible by enumeration, opaque in its mappings

A binding is a value subtyping `AbstractBinding`, the second mandatory root.
Its type declares which sides it has and enumerates what each side touches.
The legible half is explicit methods returning data, called once at attach on
the [calling task](#g-calling-task) (the task that invoked `run!`). The opaque
half is called per datum on the [device](#g-device) task by the author's own
loop:

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

The framework needs no [contract](#g-contract) on the datum's shape. The datum
travels only between `loop` and `map_input`, written by the same author, and
the framework's structural knowledge comes entirely from the declared traits
and the enumeration methods. Everything enumerable validates at attach;
everything opaque is bounded at its runtime enforcement point. `map_input` is
bounded by the staging checks ([§9.4][s9-4]). `map_output` receives exactly the
compiled gather's NamedTuple, and what it puts on the wire is the peer's
business. `map_input`/`map_output` are, precisely, **conventions of the
author-owned loop idiom**: the framework never calls them, so they are taught
([Appendix A][sA]) and never checked. A binding whose loop calls something else
by another name is simply a binding with a different private helper.

**Sides are declared; the obligations they create are enforced both ways.**
`claims` and `reads` have **error-throwing fallbacks on the root**, so a
declared side whose method was never written fails loudly at the attach
point rather than degrading into silence, and the attach runs a
**bidirectional conformance check** over the pair (trait, method):

- `is_input && !is_greedy` ⇒ `claims(b)` is called once and its [faces](#g-face) staked;
  the fallback firing here means "you declared an input side and wrote no
  enumeration".
- `is_input && is_greedy` ⇒ the [claim](#g-claim) is computed ([§9.3][s9-3]), and a `claims`
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
  the reflection class ([§11.1][s11-1]), where the shadowing check is an
  `isdefined`/`!==` pair. The `which` detection runs once at a stopped-sim
  service point, not inside any frame.

Every violation in that list reports `BindingContractMismatch`
([Appendix C][sC]). The report names the binding type, the trait, the method
at fault and the direction (declared-but-missing, or defined-but-undeclared).

**This is what closes the shadowing hole.** Under the rejected alternative,
detection by method presence (row 177), a bidirectional binding whose `claims`
was written without extending the framework's generic presented as output-only
and degraded *silently*. That omission is the `using`-without-`import` trap
([§11.1][s11-1]), one level down. With the side declared, the absent method has
something to contradict.

Greediness stays orthogonal to `reads`: a greedy front end may also drive a
compiled output gather. That combination is legal and currently
uninstantiated; its plausible customer is a narrow-wire interactive surface, a
motorized control board whose detents must be driven back out. The binding
stays an `attach!` argument, never a device field: the same `T16000M` binds
differently per aircraft, and narrowing the binding narrows the claim
([§9.3][s9-3]).

**Why the [periphery](#g-periphery) gets roots where [components](#g-component) have one.** [§11.5][s11-5] refuses
a class supertype for two reasons, and neither reaches here. First, a
component's single-inheritance [slot](#g-slot) is *already spoken for* by the
domain hierarchies (`AbstractAircraft`, engine families), while a device's and
a binding's are vacant — nothing else wants them. Second, a component's class
is implementation detail behind its contract ([§11.3][s11-3]), while a
binding's **sidedness is its public contract**, the one thing every consumer of
it must know.

Rejected, correspondingly (row 177): an abstract binding-type *taxonomy*
encoding the sides, optional roots left unenforced, and a declared `sides(b)`
trait returning the side set. The last of the three is **answered rather than
repeated** by the design above, since redundancy *with a cross-check* is drift
detection. That is what the bidirectional check turns the traits into: the
same fact stated twice, in two registers, with the framework paid to compare
them.

**`is_greedy` is a claim source, not a device class.** What the declaration
buys is one computation at the attach point. After it the binding holds an
ordinary claim set, and every mechanism downstream is blind to where the set
came from: exclusivity and storage ([§9.3][s9-3]); shape, shim, merge, scatter
and [drain](#g-drain) ([§9.4][s9-4]); [trace](#g-trace) ([§9.5][s9-5]);
detach's release. There is no derived surface shared among the writers that
elected it, and no device class hanging off the marker. The standing
rejection is untouched: opportunistic writing to unclaimed faces
"for any device" stays dead (row 44). Autonomous devices still enumerate,
and a maximal surface is what exactly one line of a binding asks
for, in the open, checked like any other claim.

**A second greedy attach stakes the empty remainder.** The complement is
computed against the [roster](#g-roster) as it stands. A greedy binding attached
after another has already swallowed everything therefore gets the empty claim.
That claim is legal, being the honest may-write-nothing degenerate below, and
it is useless, which is worth saying out loud. The attach succeeds and reports
`EmptyGreedyClaim` ([Appendix C][sC], a service warning naming the device and
its binding) — the one honest reading of "you asked for what is left and
nothing was left".

**Several interactive front ends may be rostered at once**, a web console
claiming the autopilot faces beside a local GUI claiming the stick faces. With
explicit claims they are simply two enumerated devices, partitioning the
surface rather than sharing it. The one thing still limited to a single holder
is `needs_calling_task` (the affinity check, [§9.3][s9-3]), which is a
property of the task topology, not of interactivity.

**The shipped GUI binding is a greedy one.** It declares `is_input` and
`is_greedy` and stakes the computed claim — everything unclaimed at the moment
it attaches — and defines no `claims` of its own. It declares no `reads`
either, because its read path is the handle's primitive read. VSync-paced, it
reads `latest` afresh each render ([§10.3][s10-3]), with an ad-hoc,
render-time read set over the whole [snapshot](#g-snapshot) — the inspection
register's shape ([§9.2][s9-2]). The compiled output gather therefore has
nothing to do for it. The same GUI device type is equally attachable under a binding
that returns explicit claims: greediness is the binding's declaration, not
the device's nature. Every other interactive front end anyone might want (a
web console, a remote panel) has both spellings available — attach with the
greedy binding where the GUI would have been, or with explicit claims
beside other front ends.

**The empty enumeration is not a back door.** `is_input(b) = true` with
`claims(b) = ()` stays an honest degenerate: a device that may write nothing.
Its writes are still binding-bounded, so drift onto any face is
`OutOfClaimEntry`. There is no privileged class for it to promote into either.
`claims` bodies are ordinary code (the idiom, [§11.5][s11-5]; comprehensions
included), and an enumeration that came back empty by accident stays inert,
exactly as written. The maximal surface is reachable only through the explicit
`is_greedy(b) = true` declaration: the most privileged claim is the hardest to
acquire by accident, and a declared trait is deliberate authorship. For the
same reason, `claims` never returns `nothing` or a sentinel to mean "compute
it for me". The enumeration contract has one meaning and the trait carries the
other; a dual-meaning return would be exactly the ambiguity the declaration
vocabulary is built to refuse.

#### One shipped binding type; conditioning has an owner

`TableBinding` is *data-driven*. The framework writes its `map_input` once,
and a table value (axis/button entry → [face](#g-face), deadzone/expo
parameters) is constructed per [device](#g-device) × aircraft pairing, where
configurations are made:

```julia
TableBinding(stick_y  = (face = "elevator", deadzone = 0.05, expo = 0.6),
             throttle = (face = "throttle",),
             trigger  = (face = "brake_count",))       # levels doctrine: a counter
```

Its generic `map_input` *is* the shared pure conditioning helper
([§9.4][s9-4]), and its owner. The entry tuple rides in the type, so the
mapping specializes per table with no dynamic dispatch. A *code-driven*
binding looks identical to the framework: a JSON telecommand peer whose
`claims` returns the vocabulary and whose `map_input` parses bytes. Purity
note, taught in [Appendix A][sA]: cross-datum state — press counters, edge
detection — lives in the device struct, maintained by the loop, and arrives
*inside* the datum. `map_input` stays pure.

#### Bad datum versus bug: two classes, two fates

A datum that cannot be mapped for environmental reasons — a truncated
datagram, malformed JSON, an out-of-range field — is tolerated *in the loop
body*: catch, stage nothing, `report!(handle, MalformedDatum(cause))`,
continue. That tolerance is bounded by the [device](#g-device)'s own
[diagnostic cell](#g-diagnostic-cell) (the ring and suppressed counts,
[§9.8][s9-8]; the stream, [§13.2][s13-2]), and what it records is visible next
to a live heartbeat. Any other exception propagates, and the wrapper turns it
into `DeviceCrash` ([§10.4][s10-4]).

The classification is the author's — only they know their parser — exactly as
FlightCore's `InputMappingError` docstring assigned it. What changes under the
author-owned loop is that no framework per-iteration catch site exists, so the
framework's contribution is the diagnostic channel, not the catch; a marked
exception type is not provided (row 105). `report!(handle, ...)` writes
device-attributed runtime warnings into that device's diagnostic
[cell](#g-cell), the single-writer entry point into the runtime warning stream
([§13.2][s13-2], [§9.8][s9-8]), and nothing more. It is not a general
user-diagnostics channel. Tolerating everything hides bugs as "device
attached, nothing happens"; tolerating nothing kills a live telemetry link on
its first truncated datagram, and since tasks are per-run artifacts (row 93),
kills it for the rest of the run.

### 9.7 The GUI write path: port resolution, peek, staging contract

Panels remain per-[component](#g-component) extensions in FlightCore's style — `GUI.draw!(ctx,
::LowPassFilter)`, discovered by walking the [assembly](#g-assembly) — but widgets name **the
component's own [ports](#g-port)**, never root [slots](#g-slot). The build-time wiring answers, statically
and exactly: *is this port transitively driven by a root input slot, and which one?*
Every input port has exactly one source ([§6.1][s6-1]), so the resolution is total:

- **root-driven → live widget**: [peeks](#g-peek) and stages the resolved slot through the
  GUI's own [staging cell](#g-staging-cell);
- **component-driven → read-only rendering**: displays the driven value from the
  [snapshot](#g-snapshot), visually distinct, with the source as provenance ("driven by
  `avionics/throttle_cmd`" — the canonical slash form of [§11.6][s11-6]).

This retires FlightCore's dead-slider convention and replaces it with checked
structure: **a widget is live exactly when the underlying input is yours to
command in this configuration.** The dead slider is the `Cessna172Xv1`
throttle: the engine panel's slider is visually live and the avionics silently
overwrite it every cycle, so who commands what lives in the user's head.
User-commandability is a wiring decision made where configurations are made.
Command-plus-manual-override is a mux component with a root-wired select —
explicit structure, not two writers racing (the same race as the drag phase,
[§15.3][s15-3], ruled out the same way). The obligation this places on the
GUI: read-only rendering is first-class, not an error state — the author of
`input_slider!` cannot know at authoring time whether it will be live.

**Liveness is a [derived property](#g-derived-liveness), and resolution is transitive.** A widget
is live iff two things hold: its port's feed chain terminates in a root slot,
*and* that slot lies **inside the GUI's own [claim](#g-claim)** in the run's
frozen surface partition (slot exclusivity, [§9.3][s9-3]). The feed chain is
walked through wires and [boundary](#g-boundary) connections across *all*
levels, not just the local assembly. The claim may have been computed from the
unclaimed complement under `is_greedy`, or enumerated [face](#g-face) by face
by a partial-claims binding ([§9.6][s9-6]); either way, "live" reads as
"inside the surface I declared for".

Under the [roster](#g-roster) freeze, liveness is a static fact of the run: baked once, with the
port resolution, when the run starts — never consulted against mutable claim
state at render. There is no per-port "GUI-controlled" marking anywhere. The
export chain is the marking, written by the one author entitled to write it: a
component's ports become GUI-commandable exactly when the assemblies above
surface them. The switch between "driven by its own panel" and "driven by an
external provider" is therefore automatic. At build time it follows the wiring
archetype — a scripted `World` wires a [scenario component](#g-scenario-component)
into the same faces the interactive `World` exports to root — and at run start
it follows roster claim state. Rejected: nominally-connected ports with a GUI
*override* channel (row 45). The honest cost stands: **unexported ports are
unpokeable** — FlightCore's poke-any-`u` workflow does not survive
[contract](#g-contract) visibility ([§11.3][s11-3]), deliberately.

**Peek rule:** a widget displays its **own pending write if any, else the snapshot
value**. Own-[cell](#g-cell) only: another [device](#g-device)'s pending write is invisible
by design. Its applied value arrives via the snapshot one frame later, and
cross-device peek is rejected (row 26). While paused, staged edits display
indefinitely and apply at the un-pause [drain](#g-drain) (the frame-top swap that publishes
staged device inputs into the root slots). Fan-out is consistent for free:
widgets on ports resolving to the same slot peek the same pending value.

**Staging contract: widgets stage on interaction events only.** Value widgets (sliders, drags) stage
the new absolute level on edit. Edge widgets (buttons) stage on [activation](#g-activation), as a
level computed from the peek: a flaps button peeks the current counter `k` and
stages `k+1`. The levels doctrine makes this safe by construction. Repeated
staging of the same level within a drain window is idempotent, so there is no
repeat-increment hazard. Multi-click within one window counts correctly through
the own-pending-first peek (`k` → stage `k+1`; second click peeks pending
`k+1` → stages `k+2`). Held buttons do not re-stage: after the drain applies
and the snapshot catches up, re-staging from the peek would auto-repeat at
frame rate. The activation edge is the intent.

The alternative — active widgets staging on *every* render pass — is rejected
(row 26): under slot exclusivity ([§9.3][s9-3]) it has no motivation. Side
benefit: staging traffic (and trace noise) drops from
render-rate-while-grabbed to actual edits.

No claim-transition policy exists, because no claim transition can occur
mid-run (the freeze, [§9.3][s9-3]). The one liveness-adjacent display rule is
the orphan case: a read-only widget whose claiming device's task has died
renders the fact in its provenance — "claimed by `T16000M` — task dead", the
heartbeat surfaced in place ([§10.2][s10-2]). An orphaned slot is therefore
visible where the user is looking, not only in the status panel.

The panel-authoring calling convention — what the drawing context carries,
how widgets name their component's ports, how an assembly's panel composes
its children's — is deferred to migration ([§16][s16]), where it is co-designed
against the GUI library. Its constraints are fixed here: panels name their
own ports by face-name string; resolution to root slots and the liveness
verdict are baked at run start, never performed at render; liveness and peek
arrive through the framework-supplied context, never by reaching into the
loop; and assembly panels compose children by path.

### 9.8 Diagnostics and liveness: the per-writer cell

The chapter's two data channels are specified down to their memory ordering.
The runtime warning stream ([§13.2][s13-2]) and the liveness heartbeat ([§10.2][s10-2]) are
a third, and they cross the same task [boundaries](#g-boundary): they are written by the
[device](#g-device) tasks — `OutOfClaimEntry`, `ClaimedFaceEntry` and `EntryTypeMismatch`
at staging ([§9.4][s9-4]), `MalformedDatum` from the author's loop body via
`report!(handle, …)` ([§9.6][s9-6]) — and by the loop itself (`ChatteringBudget`,
`FiringBudget`, `DebtReanchor`), and read by the loop, which folds them into
the published [framework status](#g-framework-status) ([§9.2][s9-2]) and hence into every [snapshot](#g-snapshot). An
unspecified structure with those writers is exactly the arbitrary shared
mutable state the two rules ([§9.1][s9-1]) exist to eliminate, so it gets the mechanism
[§9.4][s9-4] already established, not one of its own.

**One [diagnostic cell](#g-diagnostic-cell) per writer — one per rostered device, one for the loop
itself.** Each [cell](#g-cell) has a single writer, the same ownership argument as the
[staging cells](#g-staging-cell): no locking, no arbitration, no new primitive. The cell holds a
**bounded accumulation** — a small ring of diagnostic values, capacity **16**,
plus a per-kind count of what the ring could not hold — and one atomic
liveness timestamp.

**That bound *is* the rate limit.** When a writer emits past the ring's
capacity within one frame, the entry is not stored and its kind's suppressed
count increments; the drop policy is earliest-in-frame retained, excess
becomes counts — the first occurrences are the ones with diagnostic content,
the hundredth is noise the count already reports. [§13.2][s13-2]'s "rate-limited
wherever its source can repeat" is therefore not a policy layered over the
stream but a structural property of the channel that carries it: a [chattering](#g-chattering)
model or a peer flooding malformed datagrams costs at most sixteen retained
values and one integer increment per kind per frame, whatever its source
does, and no writer can starve another — the cells are disjoint.

**The [drain](#g-drain) is [§9.4][s9-4]'s drain.** One `atomicswap` per cell at frame top, at
the same point and under the same indivisible-take argument as the staging
drain; what the loop swaps *in* is a shared **empty sentinel**, so a quiet
frame swaps the sentinel in and gets the sentinel back — no allocation, and
no load-only code path that goes untested on healthy runs. The take is also
what makes publication sound: the batch is exclusively the loop's before it
is ever reachable from a snapshot, so the binding rule ([§9.2][s9-2]) — nothing
reachable from a published snapshot is ever written again — holds by
construction, and the live accumulator is never reachable from a published
value.

**The heartbeat rides in the same cell**, as an atomic timestamp field the
device task stores on every loop pass from inside the handle primitives
([§9.6][s9-6]: the framework observes activity without owning the loop body) and
the loop acquire-loads at the drain. There is no separate liveness channel
and no second registry: a device that is alive is a device whose cell carries
a recent timestamp, and the 2 s staleness threshold ([§10.2][s10-2]) is read against
this field. The heartbeat is not a diagnostic kind — it is a field, always
present, never enumerated in [Appendix C][sC].

**The published framework status is a concrete frozen value.** Per writer it
carries `recent` — the ring this boundary drained, at most sixteen entries;
`suppressed` — the per-kind counts the ring refused this boundary; `totals` —
the cumulative per-writer × per-kind counts since the run began, owned
privately by the loop and *copied* into each status; and `heartbeat`. Beside
the per-writer records ride the pacer diagnostics ([§8.7][s8-7]). Delta
plus total is what makes the status legible at any reading cadence: a GUI
panel refreshing at 60 Hz sees each occurrence once in `recent`, while a
consumer that samples occasionally still reads a complete account from
`totals` — nothing is lost by not looking.

**Presentation is where `maxlog` lives.** A status renderer prints a
given writer × kind up to **25** cumulative occurrences and then switches to
count-only display ("`MalformedDatum` from `UDPInput#3`: 1 482 occurrences").
That threshold is presentation policy, not channel policy: counts keep
accumulating regardless, nothing recorded depends on it, and the choice
belongs to whoever renders. The channel's own bound, above, is the one that
is normative.

**The terminal snapshot carries the run's final cumulative counters**
([§10.4][s10-4], [§13.5][s13-5]), so an [unattended run](#g-unattended-run) (a run with empty staging and no snapshot
readers) that nobody watched still answers "what went wrong, and how often"
from the value its own shutdown published.

**Allocation.** On a quiet frame there is **zero additional heap
allocation**: the sentinel swap allocates nothing and the per-writer status
rides inline in the one per-boundary snapshot allocation [§9.2][s9-2] already
accepted. That requires the per-kind counters to be a **fixed-shape isbits
record, never a `Dict`** — licensed by the closed kind set ([Appendix C][sC]), which
makes the counter layout a type rather than a lookup. On a noisy frame the
diagnostic values are allocated at emission, on the writer's own task; a
drained non-empty ring is frozen into the snapshot and can never be written
again, so the writer allocates a fresh ring lazily at its next emission —
that cost, too, landing on the writer's task — which is the same
GC-over-reuse trade [§9.2][s9-2] makes when it rejects preallocated snapshot
[buffers](#g-buffer). The rate limit is therefore an allocation bound as well: one ring of
sixteen entries per writer per boundary is the worst case, everything past it
an integer increment. The zero-allocation invariant ([§7.5][s7-5]), scoped to the model
[sweep](#g-sweep), is untouched — the cells sit on the framework side of that scope with
publication and logging.

Composition with the log, worth stating once: because the log retains
snapshot references ([§9.2][s9-2]), `totals` is monotone across logged snapshots,
so `log_every` [decimation](#g-decimation) (the log's keep-every-kth retention policy) loses
*which* boundary within a skipped stretch an occurrence fell on, never *how
many* occurrences there were.

Rejected: a shared queue under a lock, a status referencing the live
accumulator, ring reuse by double-buffering, and unbounded accumulation
(row 136).

---

## 10. Runtime periphery: lifecycle and orchestration

Where [§9][s9] fixes how data crosses the loop [boundary](#g-boundary), this chapter covers the
machinery that drives the loop itself: the [control plane](#g-control-plane) and the scheduling
primitives, the shutdown protocol, and the run lifecycle from `init!` through
[replay](#g-replay).

### 10.1 Control plane

Pause, un-pause, pace changes, `margin` changes, stop: a few scalar fields on a
separate atomic
surface, consulted by the loop at frame top and inside its wait and pause states.
`margin` ([§8.7][s8-7]) rides here for the same reason `pace` does — it tunes the
wait, never the arithmetic — so retuning the coarse/spin split mid-run is safe
by construction. The stop's issuers are the operator's channels — GUI button,
[device](#g-device) handle, calling code — and, in an interactive session, Ctrl-C: an
[operator interrupt](#g-operator-interrupt) is caught at one of the loop's unmask points and sets exactly
this stop, no separate entry point involved ([§10.4][s10-4]).
**Not staging, structurally:** staged writes apply at [drains](#g-drain), and a paused loop
drains nothing — un-pause via staging would deadlock by construction. Riding outside
the drain/[trace](#g-trace) path is safe for determinism precisely because [§8.7][s8-7] put [pacing](#g-pacing)
outside the semantics: control changes *when* frames execute, never what they
compute (stop merely truncates the trajectory). While paused the loop blocks on a
condition (notified on un-pause and stop), not a spin.

### 10.2 Loop scheduling: wait primitive, yields, thread budget

[§8.7][s8-7] fixed the shape of the pacer's wait, hybrid sleep-then-spin, but
left the coarse phase's primitive open. That choice is a scheduling decision
rather than an arithmetic one: it settles what else can run while a frame
waits. It is made here, together with the two questions that trail it —
whether a frame is guaranteed to yield at all, and how many threads a session
needs.

**Rule.** The coarse phase uses task-yielding `sleep`; there is no
`systemsleep` variant (row 27).

**Why.** What the choice buys is the wait [slot](#g-slot). `sleep` releases the
loop's thread, which makes the pacer's wait the natural scheduling window for
co-resident [device](#g-device) tasks. The design already spends that slot
twice, as the staging slot ([§8.7][s8-7]) and as the [drain](#g-drain) source
([§9.4][s9-4]).

A `systemsleep` variant for dedicated-thread hard-RT deployments is a
[guarded addition](#g-guarded-addition) (a capability the design admits but does not build).

**Rule.** With devices attached, every frame yields at least once.

**Why.** The rule is semantically free: [pacing](#g-pacing), and hence
yielding, is outside the semantics ([§8.7][s8-7]).

The yield is implicit in the coarse-phase `sleep` whenever that phase runs. An
explicit `yield()` covers the frames where it does not run: unpaced runs, and
pure-spin frames with budget ≤ margin. The spin phase itself never yields.
Yielding there would trade its µs precision for scheduler noise.

The consequence is a bound on thread occupancy: the loop holds a thread for at
most one frame before the scheduler can run anyone else. Julia's
cooperative-scheduler freeze requires a thread monopolist — a never-yielding
task that holds its thread forever — and that precondition is structurally
absent from framework tasks.

**Rule.** The thread budget is a documented sizing rule and a startup warning,
not a hard error (row 27).

The freeze FlightCore's `nthreads` error prevented cannot reproduce here, for
three reasons. The loop yields every frame. Nothing couples a stall to anyone
else, the GUI least of all: it waits on nothing, ever — a
[snapshot](#g-snapshot) acquire-load, its own [staging cell](#g-staging-cell),
atomic control. And the GUI runs on the *calling* task, so it cannot fail to
be scheduled. Under any starvation, then, the window keeps rendering and the
stop button keeps working. Undersized sessions degrade to laggy inputs and
stale snapshots, which are visible, recoverable states.

`run!` warns when `Threads.nthreads()` is tight for the attached population,
naming the `julia -t` remedy. That is one check per run, against the frozen
[roster](#g-roster) ([§9.3][s9-3]). The sizing guidance behind it: one thread
for the loop, the main thread for the GUI, and headroom for compute-heavy or
blocking-ccall devices. libuv-backed I/O yields; raw blocking ccalls pin their
thread for the duration. No pinning, no sticky tasks.

**Liveness heartbeat.** Since starvation is survivable, it must be
diagnosable. The record is the published
[framework status](#g-framework-status), the frozen diagnostics value each
snapshot carries beside the table. It includes per-device liveness — last-staged
and last-read wall time, task state — next to the pacer diagnostics. The
mechanism is the per-writer [cell](#g-cell) and nothing besides, specified in
full by [§9.8][s9-8]. A starved, blocked or crashed device task shows in the
GUI as a stale heartbeat with a name on it, not as mysteriously frozen physics.

**Stale means a liveness timestamp more than 2 s behind wall clock.** The
threshold is deliberately loose, because the heartbeat is advisory: a liveness
display and a provenance record, never a kill trigger, never a detach. It must
also tolerate a device legitimately parked in a blocking read between rare
data.

### 10.3 The next-snapshot wait

Rate-matched output [devices](#g-device) — telemetry, disk streaming — act once
per [boundary](#g-boundary). What they need from the framework is a way to learn
that a boundary has happened without polling for it.

**Rule.** Two artifacts provide it: a monotonic
**[boundary counter](#g-boundary-counter)** published with the
[snapshot](#g-snapshot), plus one `Threads.Condition`.

The counter counts *published boundaries* — grid, `t*`,
[boundary zero](#g-boundary-zero) ([§8.4][s8-4]) — not frames. Consecutive wakes
are therefore not necessarily `h` apart.

The loop's publication is `lock; counter += 1; notify; unlock`, nanoseconds of
framework-only code. Waiters never block it: one parked in `wait` has released
the lock as part of parking.

The device side is `wait_next_snapshot(handle)`, which blocks until
`counter > last_seen && running` under the canonical
[predicate](#g-predicate)-loop idiom. That idiom handles waiters at different
paces, frames skipped while transmitting, and shutdown, all with no per-frame
reset. Shutdown works because [§10.4][s10-4] wakes all waiters and each
predicate then routes its owner out.

An `Event` latch is the wrong primitive here (row 28). Conditions carry no
facts, only "look again"; the facts that matter, the counter and `running`,
live in state each waiter tests privately.

**Counter home and publication order.** The boundary index is carried *in* the
snapshot, with `t`. Any holder of one therefore indexes it without consulting
the loop — the log, an error's [replay](#g-replay) pointer ([§13.4][s13-4]), a
post-run inspector. The loop additionally mirrors the index in the state the
wait predicate tests.

**Rule.** The order of the two publications is normative: the release-store of
`latest` ([§9.2][s9-2]) happens **before** the counter increment under the lock.

```julia
# the loop, at every published boundary
@atomic latest = snap    # 1. release-store: the snapshot becomes reachable (§9.2)
lock(cond) do            # 2. only then the counter, under the lock
    counter += 1
    notify(cond)         #    parked waiters wake and re-test their predicate
end
```

`counter > last_seen` therefore implies that `latest` holds at least that
boundary, and a waiter can never wake onto a stale snapshot. The converse —
observing a *newer* snapshot than the increment that woke you — is expected and
correct: newest-wins.

**Semantics: newest-wins, no queues.** A slow consumer skips frames and always
receives the current world. This mirrors the inbound side: [coalescing](#g-coalescing) to the
newest batch (in) and to the newest snapshot (out) are the same ZOH decision.
No backpressure exists in either direction, and the loop never waits on anyone.
Rejected: per-consumer every-boundary queues (row 28).

The GUI does not use the wait. Being VSync-paced, it reads `latest` at each
render.

### 10.4 Shutdown protocol

A run ends when its time is up, when someone stops it, or when the model itself
says so. Whatever the cause, the same work has to happen: the loop has to stop
on a [boundary](#g-boundary) rather than mid-frame, every device task has to be
woken out of whatever it was blocked on, and every resource acquired for the run
has to be released before `run!` returns. This section specifies that sequence,
which it calls the tail throughout. It also specifies the two things that
bracket the tail. One is the pre-spawn initialization at the top of a run, which
is a step of this same protocol. The other is the operator interrupt, which
enters the tail rather than bypassing it.

#### The tail: the ordered sequence every stop takes

**Rule.** Steps (1) through (5) below run in that order on every stop, whatever
initiated it. Steps (6) and (7) specify what happens when a device task or the
loop itself ends first.

1. **Initiation.** Three events start a shutdown: `t_end` is reached, a
   control-plane stop is issued, or a `stop_on` [face](#g-face) reads `true` in
   the just-published [snapshot](#g-snapshot). The stop's issuers are the GUI, a
   [device](#g-device) handle, code, or an
   [operator interrupt](#g-operator-interrupt) — Ctrl-C, treated below. The
   third event is model-detected termination ([§13.5][s13-5]). The loop always
   completes the current boundary sequence and never stops mid-frame. It then
   publishes the final snapshot. Only then does it set the sticky stopped
   status.
2. **Wake all framework waits.** The waits are the next-snapshot wait and the
   pause. Each waiter observes the stopped status and unwinds. A stop issued
   while paused therefore works.
3. **Unblock device-specific blocking calls.** The hook is `unblock!(device)`,
   default no-op. A network input's override closes its own socket, which raises
   in the blocked task. The framework wrapper catches that raise and treats it
   as shutdown. This demotes FlightCore's EOT convention from load-bearing
   shutdown mechanism to an optional wire-protocol courtesy between remote
   peers.
4. **Loop bodies exit.** The exit is the author's own `while running(handle)`
   loop, the [contract](#g-contract) of [§9.6][s9-6]. That contract teaches two
   obligations: the [predicate](#g-predicate) check and interruptible blocking.
   Steps (2) and (3) are what make every blocking point interruptible. The
   wrapper's `finally shutdown!(device)` is guaranteed on every exit path.
5. **Join with a timeout of 5 s.** A device task exceeding the timeout is
   reported *by name*, through the [§10.2][s10-2] heartbeat. It is then
   abandoned with a `DeviceJoinTimeout` warning ([Appendix C][sC]) rather than
   left to hang `run!`.
6. **Device-initiated paths.** A device exits voluntarily when its loop body
   returns, at a window ✕ or a peer EOT. No `should_close` hook exists
   ([§9.6][s9-6]). With `should_abort` set, the wrapper's exit path also
   requests a sim stop. Otherwise the sim continues with the device's *task*
   absent: its [cell](#g-cell) stops filling, and the loop is structurally
   indifferent. Its [roster](#g-roster) entry and its [claims](#g-claim) persist
   to run end, because [§9.3][s9-3] freezes the roster for the run and death is
   not detach. The orphaned [slots](#g-slot) hold their last-drained values,
   visibly ([§9.7][s9-7]). A crashing device task is caught by the framework
   wrapper and follows the same path, logged with the device's name
   (`DeviceCrash`, [Appendix C][sC]).
7. **Loop-side failure.** A failure on the loop's own side runs steps (1)
   through (5) from the catch path, specified in [§13.6][s13-6]. The failed
   boundary is discarded and the previous snapshot is promoted to final.
   FlightCore's `SimulationTermination` catch path was the precedent, though the
   exception-based termination idiom itself has no place here ([§13.5][s13-5]).
   Devices therefore unwind cleanly regardless of who died.

The ordered part, in one line:

> **(1)** stop observed → current boundary completed → final snapshot published
> → sticky `stopped` set → **(2)** framework waits woken → **(3)** `unblock!`
> per device → **(4)** loop bodies return, each through its wrapper's
> `finally shutdown!` → **(5)** join under the 5 s cap → `run!` returns

**Why the final snapshot goes out before the status is set.** Publishing first
guarantees that output devices can flush the true final state. The status that
follows carries the run's cumulative diagnostic counters ([§9.8][s9-8]) — the
complete warning account of a run nobody watched.

**Rule.** That terminal snapshot is retained in the log unconditionally, under
any `log_every` and any `log_max` ([§9.2][s9-2]).

**`t_end` lands on the grid.** The run ends at the first grid boundary whose
time reaches or exceeds `t_end`. Whole frames only, never a shortened final
step, which grid integrity forbids ([§8.4][s8-4]; `tₖ = t₀ + k·h`, indexed and
never accumulated). The final boundary may therefore overshoot `t_end` by up to
`h`. The termination record carries the actual final `t` ([§13.5][s13-5]). This
is the `t_plus` spelling ([§10.6][s10-6]) applied to the run's own clock: whole
frames until the boundary time first covers the duration.

**The two termination sources differ in kind.** `t_end` is a grid fact, checked
against boundary times on the grid. `stop_on` is checked at *every* published
boundary, `t*` included ([§13.5][s13-5]).

**Why five seconds.** It is generous for GUI window teardown and socket closes.
It is short enough that an abandoned join reads as a diagnosed timeout rather
than a hang.

**The calling-task device sits outside the join.** The
[calling task](#g-calling-task) is the task that invoked `run!`. The device it
hosts is the GUI, which has no spawned task ([§9.1][s9-1]). That device's loop
body is the calling task's own occupation of `run!`. It exits by the same
`running(handle)` predicate as any device loop, and `run!` returns after the
joins. One honest asymmetry follows: the abandonment path of (5) cannot cover
it, because nothing can abandon the task `run!` stands on. A calling-task device
that blocks past shutdown therefore hangs `run!`. The trait's one authoring
obligation is a loop body that never blocks between `running` checks. The
shipped GUI's render loop polls once per frame and never blocks.

**What survives the tail.** After (5) the task set is empty, device tasks being
per-run artifacts ([§9.1][s9-1]), and `shutdown!` has released each device's OS
resources. What survives a stop is the roster entry: binding, claims, stable
device id ([§9.3][s9-3]). Never a task, never a live resource. That holds for a
device whose task died mid-run too, its entry being indistinguishable at this
point from any other's. `stopped` is where `detach!` removes an entry and
releases its claims.

**One roster change belongs to this tail.** A GUI attached by `run!`'s
`gui = true` is detached here, releasing its computed claim (the run-scoped
flag, [§10.6][s10-6]). It is the only roster mutation the protocol itself
performs. It sits in the tail precisely so that (7)'s failure path takes it too,
an everything-claim staked for one run never surviving into the next.

**The next run re-acquires everything.** The next `run!` re-runs device `init!`,
resource acquisition being per-run; FlightCore's
create-a-new-socket-each-`init!` in network.jl is the precedent. It also spawns
fresh tasks against the re-armed [§10.3][s10-3] counter. While stopped there are
no device tasks at all, so voluntary exit and the [§10.2][s10-2] liveness
heartbeat are run-scoped observables. A device unplugged while stopped surfaces
as the next run's `init!` failure, disposed of by the initialization bracket
below.

#### Initialization: the pre-spawn bracket

**Initialization is a step of this protocol, taken at the top of a run.**

**Rule.** Before any task is spawned, the loop calls `init!` once per roster
entry, in attachment order, on the calling task. Each call sits in its own
bracket.

```julia
for entry in roster                       # attachment order, calling task
    try
        init!(entry.device)
    catch e
        shutdown!(entry.device)           # release, unconditionally (§9.6)
        report!(entry, DeviceCrash(e))    # pre-spawn: the entry is the address, no handle yet
        mark_dead!(entry)                 # from boundary zero; no task is spawned
        entry.should_abort && stop!(control)
    end
end
```

**Why the bracket.** It is what makes "guaranteed on every exit path"
([§9.6][s9-6]) true of the path outside that wrapper. A device that throws
half-way through acquisition is handed back to `shutdown!` right there, so its
partially acquired OS resources are released rather than leaked. That is exactly
why `shutdown!` owes tolerance of a partially initialized device, a taught
obligation of [§9.6][s9-6].

**The report is the ordinary `DeviceCrash`, not a kind of its own**
([Appendix C][sC]). Its [payload](#g-payload) already carries everything an
init-time failure has to say: the device id, the cause exception, and whether
`should_abort` was set. The name is honest, a device that cannot acquire its
resources having crashed before it lived.

**Rule.** The report is written through the ordinary
`report!(address, diagnostic)` entry point, addressed by the roster entry rather
than by a handle.

There is no device task to hold a handle before the spawn. The address supplies
the device identity either way, which is why no call passes a device id
([§9.8][s9-8]).

**Rule.** No task is spawned for a failed device, so it is dead from
[boundary zero](#g-boundary-zero) (the initialization boundary: the ordinary
macro-sequence with an empty integrate).

That needs no machinery. Its [diagnostic cell](#g-diagnostic-cell) (the
single-writer cell each writer owns for diagnostics and heartbeat) never
receives a heartbeat timestamp. The cell therefore reads stale against the
[§10.2][s10-2] threshold from the first frame ([§9.8][s9-8]).

**Rule.** The claims of a failed device persist to run end.

This is the death-is-not-detach disposition ([§9.3][s9-3]), applied one step
earlier than (6)'s. The roster is frozen for the run, and the orphaned slots
hold their initial values, well-defined by [slot totality](#g-slot-totality)
([§14.6][s14-6]; every root slot must hold a value). An orphan of (6) holds a
last drained batch instead.

**The run's disposition splits on `should_abort`, uniformly with (6).** With the
flag clear, which is the default ([§9.6][s9-6]), the remaining entries
initialize, the run starts, and the sim runs with that device absent from frame
zero. That is (6)'s "the sim continues with the device's *task* absent", shifted
to `t₀`.

With the flag set, the failure requests a control-plane stop, and that stop is
simply *already pending* when the run reaches boundary zero. This protocol
already has that path: the boundary-zero check ([§13.5][s13-5]) ends a run at
`t₀` with that snapshot final, integrating nothing. No new exit protocol is
needed, therefore.
The remaining entries still initialize, every rostered device getting its
`init!`/`shutdown!` pair uniformly. The run publishes boundary zero. It ends
`stopped` at `t₀` through this same tail, with the termination record naming the
source ([§13.5][s13-5]). What the operator is left with is an ordinary stopped
simulation: a terminal snapshot, with the failure named in its diagnostic
account. It is fully serviceable by [§14][s14] and resumable by the next `run!`
([§10.6][s10-6]) once the device is plugged back in.

**Topology is derived after initialization**, not from the roster alone
([§9.1][s9-1]): a `needs_calling_task` holder whose `init!` failed returns the
loop to the calling task, which would otherwise be pinned waiting to run the
loop body of a device that does not exist. The shipped GUI attaches with
`should_abort = true`, so in practice that run ends at `t₀` anyway. The rule is
stated generally because it costs nothing.

#### The operator interrupt

**The operator interrupt is a stop, not a failure.** Ctrl-C in an interactive
session is a control-plane stop command issued by hand. The run completes the
current boundary, publishes the final snapshot, takes this tail like any other
stop, and ends `stopped`. The result is boundary-consistent. It is fully
serviceable by the [§14][s14] stopped-sim services and resumable by the next
`run!` ([§10.6][s10-6]).

The interrupt is the escape from a run nothing else can end: deviceless,
`t_end = Inf`, empty `stop_on` (`UnboundedRun`, [Appendix C][sC]). It needs no
entry point of its own. The stop already rides on the
[control plane](#g-control-plane), the separate atomic surface carrying pause,
pace and stop ([§10.1][s10-1]). The exceptions-are-abnormal doctrine
([§13][s13]) is untouched. That doctrine is about *model* code, while this is
the one exception whose meaning the framework knows.

**Rule.** Masking across the boundary is normative, not an implementation hint.

**Why.** An `InterruptException` is delivered asynchronously. An interrupt
landing mid-[sweep](#g-sweep) would therefore destroy the boundary this protocol
is built on completing, and would leave half-written stores ([§13.6][s13-6]).
That forces a choice between `stopped` with dirty stores and a terminal
`errored`. The `stopped`-with-consistent-stores guarantee is exactly what the
masking buys.

The loop masks delivery across the boundary macro-sequence, using Julia's
`disable_sigint` — a sigatomic counter increment, negligible per frame. It takes
the deferred raise at the unmask points: the frame top, where it already
consults the control plane ([§10.1][s10-1]), and inside its wait and pause
blocks. All of those points are boundary-consistent. Caught at one of them, the
interrupt sets the control-plane stop and enters this tail. The catch site
([§13.4][s13-4]) therefore never sees it.

**A second interrupt during the tail** collapses the remaining joins
immediately. That is (5)'s abandonment path taken at once, with devices still
reported by name (`DeviceJoinTimeout`). The run still ends `stopped`: escalation
shortens the tail, never reclassifies the run. Nor can a second interrupt repair
(5)'s honest asymmetry, since nothing can abandon the task `run!` stands on.

**Interactive-session scope, stated plainly.** Outside the REPL, Julia's default
(`exit_on_sigint(true)`) kills the process on SIGINT before any of this
machinery runs. The framework flips nothing process-global.
[Unattended runs](#g-unattended-run), those with empty staging and no snapshot
readers, rely on `t_end` and `stop_on`, as they already must.

### 10.5 Scripts and the mid-run mutation doctrine

What the consumers demonstrably mutate mid-run, surveyed: FlightCore's
`user_callback!` has exactly two archetypes. The first is the timetable script
(c172_demos.jl:290: `elevator_offset` as a function of `t`). The second is the
synthetic pilot (c172_demos.jl:423, 525: a phase FSM reading `y` and writing
mode requests, references, flaps, wind). Both write only `u` fields. No demo,
test or GUI path pokes `x`/`s` mid-run, and `init!`/trim appear only between
construction and `run!` (c172_demos.jl:303).

**Sim-time scripts are model behavior, so they become
[scenario components](#g-scenario-component)** (ordinary periodic discrete
components holding a sim-time script). Both archetypes are clocked by *sim
time*: `t`, the trajectory. Mapping them to [devices](#g-device) is rejected
(row 31). The clock is the criterion.

**Rule.** A sim-time script becomes a source or supervisor
[component](#g-component); a wall-clock interaction becomes a device.

A script mapped to a component is periodic discrete, with `K = 1` for today's
`dt = 0.02` callbacks. It executes synchronously in the loop, deterministic
paced or unpaced. It is replayed by recomputation, with no [trace](#g-trace). A
device, by contrast, is traced and replayed from the trace.

**The component mapping is strictly richer than the callback it replaces.**

- The `Ref(:init)` phase closure becomes honest `x`, visible in
  [snapshots](#g-snapshot), logs and plots.
- Inputs arrive same-[boundary](#g-boundary) fresh by topological order; the
  callback ran post-step, one boundary staler.
- The pure timetable script is a one-liner reading the clock out of its
  [bundle](#g-bundle) (the NamedTuple of zero-copy views a component function
  receives). That one-liner is `h_xu(s, (; t)) = (; offset = profile(t))`,
  exact at its own [ticks](#g-tick), with no latching.
- In a scenario configuration the script drives the avionics' input
  [ports](#g-port). [§9.7][s9-7] therefore renders the corresponding GUI widgets
  read-only with provenance — today's demo-vs-GUI dead-slider fight, resolved
  by the port-resolution rule.

**`user_callback!` is eliminated** (row 31). It is the
[periphery](#g-periphery)'s `f_step!`, and cheap composition leaves it without
justification. Its call sites migrate to scenario components, not devices.

**Manual event triggering needs no mechanism.** It takes a root input
[slot](#g-slot) plus a [boundary-detected](#g-boundary-detected)
[guard](#g-guard) reading that slot — an edge check at step boundaries only,
with no root-finding. That is already expressible in settled machinery. The
levels doctrine applies: latched commands or counters. The demos' engine
start/stop buttons are `u`-writes today.

**Mid-run re-initialization is not built, because it is not demonstrated.**
Initialization and trim are stopped-sim workflows (first-class services,
[§14][s14]). No concurrency perimeter exists there: no loop, no devices, plain
single-task code. The guarded-addition shape is on record should demand appear.
It would be a traced, boundary-executed intervention command applied through
project → [sweep](#g-sweep) → publish, so that no consumer ever observes
un-decoded state.

**The doctrine, final form.** While a simulation runs, the periphery stages
root-input writes and issues control commands. Structurally, it does nothing
else. Anything that wants to poke the model mid-run is one of three things. It
is an *input* in disguise, so wire a slot and a guard. It is *model behavior*
in disguise, so add a scenario component. Or it is a *wall-clock interaction*,
so attach a device. Graceful termination follows the same shape
([§13.5][s13-5]): a declared stop [face](#g-face) in the model, plus `stop_on`
policy at deployment. Never a callback, and never a thrown exception.

### 10.6 Run lifecycle and partial advance

A `Simulation` moves through five states: **built**, **initialized**,
**running**, and terminally **stopped** or **errored** ([§13.4][s13-4]).
**Built** is stores allocated, nothing authored. **Initialized** is `init!`
completed [boundary zero](#g-boundary-zero), the initialization boundary run as
the ordinary macro-sequence with an empty integrate ([§14.5][s14-5]).

**Rule.** `init!` is mandatory.

`run!` or `step!` on a simulation whose [boundary](#g-boundary) zero has not
completed is an error in the kind set ([§13.2][s13-2]) naming `init!`. That is
distinct from `UninitializedSlots`, which fires *inside* `init!`
([§14.6][s14-6]). `replay!` is the one alternative entry: it runs boundary zero
from a [trace header](#g-trace-header) ([§10.7][s10-7]).

**Where the loop runs.** The loop runs on the [calling task](#g-calling-task),
the task that invoked `run!`, unless a calling-task [device](#g-device) is
rostered. That device is the GUI, and the topology is derived from the
[roster](#g-roster) ([§9.1][s9-1]). Deviceless, `run!` is fully synchronous.
That is the unattended register: an [unattended run](#g-unattended-run) is the
same loop with empty staging ([§9.1][s9-1]). It is also what the synchronous
rethrow presupposes ([§13.4][s13-4]).

**Partial advance.** `step!(sim; frames = 1)` advances whole frames
synchronously through the ordinary frame sequence — [drain](#g-drain),
integrate, boundaries, publication — and returns. A stepped simulation is
bit-identical to the same frames under `run!`. `step!(sim; t_plus = 10.0)` is
the duration spelling, mutually exclusive with `frames`. It advances whole
frames until the boundary time first covers the duration, which is the
migration suite's advance-by-duration idiom.

Partial advance is the test-[harness register](#g-harness-register): advance,
assert, advance. It is equally the REPL register: fly a while, inspect,
continue. Neither is a script, so the scenario-[component](#g-component)
doctrine does not absorb them ([§10.5][s10-5]).

**A stepping session is deviceless by construction.** Device tasks are
per-`run!` artifacts ([§9.1][s9-1]), and a device loop's `while running(handle)`
is false outside a run. Between `step!` calls the simulation is in a stopped-sim
state — `initialized`, below — so `attach!` is legal there and does what it
always does: it registers ([§9.3][s9-3]). The task appears at the next `run!`.

**The [frame-top drain](#g-drain) still runs**, `step!` frames staying
bit-identical to `run!` frames. What it drains is the **harness
[cell](#g-cell)**. The harness write path is
`stage!(sim, "face" => value, ...)` ([§9.3][s9-3]), with the calling task as
writer. Staged batches are ordinary batches. They are traced, so
[replay](#g-replay) and bit-identity hold. They are applied at the next frame
top. They are surface-checked like any writer's ([§9.3][s9-3]).

The read half is `latest(sim)`. It hands back the same immutable
[snapshot](#g-snapshot) value a device handle acquires ([§9.2][s9-2]), navigated
directly for assertions. Advance-assert-advance is `stage!` → `step!` →
`latest`. Both entry points work under `run!` too: the
[harness cell](#g-harness-cell), the always-present staging cell of the harness
register, is not step-scoped. An inspection accessor leaves the rejection of
closure-based termination ([§13.5][s13-5]) untouched.

**Status, termination and the `run!` [seam](#g-seam).** Between `step!` calls a
simulation reports **initialized**. No loop task exists, so `running` would lie;
nothing is terminal, so `stopped` would lie too. The state reads
"boundary-consistent and ready to advance", not "sitting at boundary zero".
`run!` may therefore follow `step!`, continuing from the current boundary; so
may another `step!`.

Termination policy is honored throughout, as bit-identity requires. `t_end`
reached, or a `stop_on` [face](#g-face) [holding](#g-edge-semantics) at frame 3
of `step!(sim; frames = 10)`, ends the run there through the ordinary
[§10.4][s10-4] tail and leaves the simulation `stopped`. `step!` therefore
returns the number of frames **actually advanced** — the requested count in the
ordinary case, fewer when the run terminated inside the call. That return is how
a harness detects the truncation without inspecting the clock.

**Re-running: `stopped → init! → run!` is the supported cycle.** `init!` re-runs
boundary zero from its condition, the warm restart being `capture` → tweak →
`init!` ([§14.1][s14-1]). It clears the [trace](#g-trace), the log, *and* any
batches still in [staging cells](#g-staging-cell). The [recorders](#g-recorders)
restart with the run they record, and no stale batch survives to clobber the
boundary zero it predates.

**Device attachments persist across re-initialization**, attachment being
orthogonal to the run lifecycle ([§9.3][s9-3]). Persistence means *roster*
persistence: binding, [claims](#g-claim) and device id survive. Tasks and OS
resources do not (the per-run topology, [§9.1][s9-1]; the teardown,
[§10.4][s10-4]). Each `run!` re-initializes every rostered device and spawns its
task. `attach!` while stopped only registers, the task appearing at the next
`run!`.

**Task topology follows the roster each time** ([§9.1][s9-1]). A GUI attached
*by hand* is still rostered, so the next `run!` renders it again — loop on a
spawned task — whether or not `gui = true` is repeated.

**The `gui = true` flag itself is run-scoped.** At run entry it attaches the
standard GUI device under the greedy binding, with `should_abort = true`, iff
no GUI is rostered ([Appendix B][sB]). The run's shutdown tail detaches it again
([§10.4][s10-4]). So the roster a flagged run leaves behind is the roster it
found, and a window on every run means the flag on every run. A *persistent* GUI
session is spelled by hand: `attach!` while stopped, `detach!` when done.
Against a hand-attached GUI the flag does nothing and detaches nothing, having
attached nothing.

**What the scoping buys is the absence of a trap.** The flag's GUI claims
everything unclaimed at attach (the computed source, [§9.3][s9-3]), and a claim
of that shape must not outlive the run that asked for it. A joystick attached
between two runs would otherwise meet a `ClaimConflict` against an
everything-claim staked by a convenience argument nobody remembers passing.

**The accepted cost is a fresh device id per run for that GUI.** Ids exist to be
read *across* roster changes, and each run's trace header carries its own
schemas ([§9.5][s9-5]). Nothing that reads a completed run is therefore
affected.

**Run policy is re-bindable per cycle.** `t_end` and `stop_on` are `Simulation`
defaults that `run!` may override for the run it starts ([§13.5][s13-5]). A
second run — or a `step!` register between two runs — can therefore stop on a
different clock or a different face set without a rebuild.

**`errored` is terminal** (row 59). Reproduction is trace replay
([§10.7][s10-7]), not resurrection.

### 10.7 Replay: the trace re-drives the ordinary loop

The entry point the [§9.5][s9-5] [trace](#g-trace) exists for:

```julia
trc  = trace(sim)                     # the recorded session: header + per-frame batches
sim2 = Simulation(world)              # the same build
replay!(sim2, trc)                    # header-init, then re-drive every recorded frame
replay!(sim2, trc; to_boundary = k)   # partial: the §13.4 replay-pointer register
```

`replay!` is **the ordinary loop with exactly two substitutions**, not a
separate execution mode. That is what keeps every property proved of the loop
true of [replay](#g-replay):

- **[Boundary zero](#g-boundary-zero) from the header** — the initialization
  boundary, the ordinary macro-sequence with an empty integrate. `replay!`
  stands in the `init!` position of the lifecycle ([§10.6][s10-6]): it applies
  the header's resolved stores and [slot](#g-slot) values directly. There is no
  condition resolution, and the totality ([§14.6][s14-6]) holds by capture. It
  then executes the ordinary [boundary](#g-boundary)-zero sequence
  ([§14.5][s14-5]). Authored-condition events re-fire identically: the header
  predates the sequence (the capture placement, [§9.5][s9-5]), so nothing is
  applied twice and nothing is skipped.
- **The [drain](#g-drain) reads the trace.** Each frame top applies the
  recording's batches for that **[frame ordinal](#g-frame-ordinal)**, the frame
  index a batch replays at. It does not swap the [roster](#g-roster)'s
  [staging cells](#g-staging-cell), where a device's pending write batch waits
  between drains. Ordinal keying is exact because the frame sequence is itself
  deterministic under replay (`t*` boundaries derive from state, [§8.4][s8-4]):
  frame *k* of the replay *is* frame *k* of the recording. Recorded batches
  apply **verbatim, with no surface re-check**. The write-surface rule
  ([§9.3][s9-3]) ran at recording time, and [claims](#g-claim) are a live-roster
  fact of the recorded session that replay does not reconstruct.

Everything else is the loop as already specified:

- **Termination and partial replay.** Replay ends at the recording's final
  frame, or earlier at `to_boundary = k` — the consumer of the replay pointer
  ([§13.4][s13-4]). That spelling is defined as running **through the frame
  whose execution published boundary `k`**, so replay always halts at a frame
  top. For a grid boundary the halt is exactly at `k`, the frame that publishes
  one ending at it. The frame-entry pointer ([§13.4][s13-4]) lands the same way.
  A [localized](#g-localized) `t*` boundary inside the frame — the crossing
  instant bracketed by root-finding over probe sweeps — is reproduced but not
  stoppable-at. [§8.4][s8-4] separates the two indices: the trace stays
  frame-indexed, and boundaries are the reporting index. Replay may also end
  earlier still under the ordinary policies, since `t_end` and `stop_on`
  overrides bind for this replay exactly as at `run!` ([§10.6][s10-6]). A
  termination the recorded session hit through `stop_on` reproduces itself
  anyway, deterministically.
- **Replay ends `initialized`, never `stopped`** — boundary-consistent and
  ready to advance, the same state `step!` leaves ([§10.6][s10-6]). This is what
  makes three promised workflows real. State-trajectory inspection asks "what
  was the private state at t = 37.2?", and answers it by replaying there and
  reading the live stores ([§9.2][s9-2]). Error reproduction replays to
  `k − 1`, then `step!`s the failing frame under instrumentation
  ([§13.4][s13-4]). Continuation is `run!` after `replay!`, a live session from
  the replayed boundary.
- **Replay re-records.** The trace register runs normally: the new trace
  inherits the old header and accumulates the re-drained batches, a
  bit-identical prefix. A replayed-then-continued session therefore leaves
  behind a complete, valid trace of *itself*, with no special stitching.
- **[Pacing](#g-pacing) and the [control plane](#g-control-plane) are
  unchanged** ([§8.7][s8-7], [§10.1][s10-1]). Pacing (waits inserted between
  completed frames, never altering the boundary sequence) sits outside the
  semantics, so paused, slow-motion or real-time replay is free. Paced replay
  with an attached visualizer *is* session playback. Stop truncates, as
  anywhere.
- **[Devices](#g-device) are readers.** Rostered devices init and spawn normally
  ([§9.1][s9-1]) and consume [snapshots](#g-snapshot) ([§9.2][s9-2],
  [§10.3][s10-3]) — the visualizer case. But no live staging [cell](#g-cell) is
  drained while the trace feeds the loop: a batch found staged is discarded with
  a rate-limited warning (`ReplayDiscardedStaging`, [Appendix C][sC]). Mixing
  live writes into a replay would destroy the property replay exists to provide.
  A session that wants live input is a continuation (`run!` after replay), not a
  replay.
- **Validation is loud and up front.** Before the first frame, the header is
  validated against the `Build` (store layout, slot [faces](#g-face)), and the
  trace's batch entries against the root input-face list. The checks are
  attach-style, and a failure reports [did-you-mean](#g-did-you-mean): the
  offending name plus the list-in-hand it should have matched. The kinds are
  `ReplayHeaderMismatch` and `ReplayUnknownFace` ([Appendix C][sC]).

  The same pass pays the trace-record conversion in reverse. Every writer's
  sparse records ([§9.5][s9-5]) are normalized to positional batches against the
  header's schemas, once, off the loop — replay has the whole trace in hand
  before frame 1. The replay drain therefore applies compiled scatters exactly
  as the live drain does, and no face name is resolved per frame under replay
  either.

  *Structural* mismatch is an error; *parametric* difference is not. Replaying
  against the same structure with changed parameters is the
  **[what-if register](#g-what-if-register)** — deterministic re-driving of the
  recorded inputs through a modified model. Bit-identity is promised only
  against the identical build; the what-if register promises determinism, never
  reproduction.

  The header's deployment block ([§9.5][s9-5]) validates in the same pass, on
  the *structural* side of that line. The seven trajectory-determining
  parameters are `Δt_base`, `h`, `n`, the algorithm, `localization_tol`,
  `localization_budget` ([§8.4][s8-4]) and `firing_budget` ([§8.6][s8-6]). All
  seven are compared against the target `Simulation`'s own deployment binding.
  Mismatch is `ReplayHeaderMismatch` with a deployment-parameter discriminator,
  never a what-if. A deployment change moves the times at which the
  frame-ordinal batches apply: different inputs, not a modified model. The
  event trio — the localization pair and `firing_budget` — is compared for
  exactly the same reason the grid parameters are. It moves the trajectory,
  so a run that differs in it is not re-driving the recorded one.

  `t₀` is *applied*, not compared. Replay stands in the `init!` position and
  owns the anchor, so `replay!` takes no `t0` argument. The header's
  `t_end`/`stop_on` pair is a recorded fact of the recorded session, never a
  constraint on this one; overrides bind as stated above.

The dispositions, by header content:

| header content | disposition |
|---|---|
| store layout, slot faces | compared against the `Build` |
| the deployment block's seven trajectory-determining parameters: `Δt_base`, `h`, `n`, the algorithm, `localization_tol`, `localization_budget`, `firing_budget` | compared against the target `Simulation`'s own deployment binding |
| resolved stores, slot values | applied directly at boundary zero |
| `t₀` | applied; `replay!` takes no `t0` argument |
| `t_end`, `stop_on` | neither compared nor applied: a recorded fact of the recorded session |

Rejected shapes, for the record (row 101): a `run!(sim; replay = trc)` flag, a
synthetic playback device staging the recorded batches, and replay ending
`stopped`.

---

# Part III — Authoring and build

## 11. The declaration layer: components and assemblies

How an author spells a [component](#g-component): where the structural facts live, what the build
takes as authoritative, and what is checked against what. [§11.1][s11-1]–[§11.4][s11-4] cover the
component side; [§11.5][s11-5]–[§11.8][s11-8] the [assembly](#g-assembly) side; the build pipeline is [§12][s12] and the
stopped-sim service spellings are [§14][s14]. Concrete syntax below is near-final in shape
but still illustrative in spelling. The sketches (`sketch_decoder.jl`,
`sketch_io.jl`) are written against this layer and the services spellings.

### 11.1 Position: a declarative trait layer — plain Julia, no macros

A component is authored in ordinary Julia. Its
[stage functions](#g-stage-function) — `h_x` or `h_xu`, the two output stages
every component provides — are ordinary multiple-dispatch methods, on the
`GUI.draw!` precedent, and its structural facts are declared through a small set
of well-known functions returning plain values, defined alongside those methods.
Four questions about that layer are settled here: what it rules out and what it
still admits (macros); how an author's methods reach the framework's generic
functions, and how they can silently fail to (the namespace); which of
declaration and evaluation is authoritative (the schema); and what a component's
contract may depend on (the type).

#### Plain Julia, not a macro DSL

**Rule.** There is no macro DSL.

**Why.** The charter's debugging, tooling and comprehension workflows ([§1][s1])
decide it (row 32).

Redundancy between declarations and function bodies is accepted deliberately,
under one non-negotiable condition: **every inconsistency fails loudly**, at
build time where possible, at first execution otherwise.

A macro can only ever *lower to* a layer like this one. A convenience macro
therefore remains addable a posteriori as pure sugar — the `@kwdef` precedent —
while never becoming load-bearing.

The door stays open for the declaration layer specifically. A macro generating
the well-known declarations is admissible sugar *on top of* the plain-Julia
forms, never a replacement for them and never required to author a
[component](#g-component) (row 166). The obvious candidate is the
`where {T <: Real}` ceremony of a continuous `output_types` ([§11.2][s11-2]).
Every rule in this part is stated over the generated methods, so a macro that
lowers to them adds convenience and no semantics.

#### The namespace: declarations are extended, not called

**Rule.** The declarations and stage functions are extended, not called:
authoring a component means adding methods to framework-owned generic functions.

Julia admits that only through an explicit per-name `import`, or through a
qualified `Flight.f(…) = …` definition — the `Base.show` idiom [§16][s16]
records for the extension-only [periphery](#g-periphery) surface. A component
module therefore opens with

```julia
import Flight: init_x, init_m, workspace, input_types, output_types,
    events, h_x, h_xu, f, g, project, child_connections,
    input_connections, output_connections, sample_times
```

**Why the explicit list: `using Flight` alone is a silent trap.** After a bare
`using`, `f(eng::Engine, …) = …` defines a new, unrelated `MyModule.f` — no
error, no warning. The short names are deliberately unexported, precisely
because `f`, `g`, `events` and `project` are the most collision-prone
identifiers in numerical Julia, so there is no clash for the language to
detect. The build then
sees a component with no `f` method and reports a *modeling* diagnostic:
`StoreWithoutUpdate`, or `ClassUnreadable` when the whole inventory was
shadowed. A one-line namespace mistake is thereby reported away from the wrong
line. That is the [§11.4][s11-4] inversion of [error locality](#g-error-locality)
— the property that a mistake fails at the site of the mistake — arriving
through the namespace.

Two mitigations, both normative. The first: the import list above is
authoring surface, stated wherever a component file is first shown. The second:
the two diagnostics run a **shadowing check**. If the component's parent module
defines a same-named function distinct from the framework's, the message says so
and names the missing import: "`MyEngine`'s module defines its own `f`, distinct
from `Flight.f` — add `import Flight: f`". The check is a two-line
`isdefined`/`!==` test on names the build already looks up.

A convenience macro expanding to the import list remains addable-a-posteriori
sugar per this section's macro doctrine. A re-export submodule is not an
alternative: per-name `import` is the only extension register the language
provides (row 117).

**The same trap has a local-scope sibling** (row 164). Written inside a `let`, a
function body or a `@testset`, `h_x(::MyComp, (; x)) = …` does not add a method
to the global `h_x`; it binds a *new local function* of that name. Calls within
the block resolve to it and look correct. The generic function the build
dispatches on never learns of the component, which therefore reads as one
declaring nothing at all.

```julia
@testset "mycomp" begin
    h_x(::MyComp, (; x)) = …   #a NEW local h_x, not a method of Flight.h_x
    …                          #calls here resolve to it, and look correct
end
#outside: Flight.h_x still has no MyComp method
```

The shadowing check above cannot reach this case: there is no parent-module
binding to compare, the shadow being a local binding that disappears with its
block. So the mitigation is at the other end — **a component that declares
nothing and defines no stage is rejected at build time**, an inert component
being unwritable on purpose. That check costs a line, and it catches the
misspelled-declaration family with it.

Test code is the realistic victim: a fixture component defined inside its own
`@testset`. The authoring rule is one line — declarations live at module top
level.

The net holds under a *partially* shadowed component too, because `output_types`
is still a declaration. A component whose [ports](#g-port) are declared but whose
stage went to a local binding reads as "declared but not produced"
([§11.3][s11-3]), the shadowing note attached, rather than as a component with
nothing to say.

#### Declarations are the schema authority

**Rule.** Declarations *define* the model's structure; evaluation *checks*
conformance against them — never the reverse.

The build [probes](#g-probe) user functions with real values, with no reliance on
compiler inference, and compares observed against declared. The same comparison
then runs on every subsequent evaluation for free: a `NamedTuple`-type check
that constant-folds away when conformant.

Inference-by-evaluation as schema authority is rejected on three counts,
established by walkthrough ([§11.4][s11-4]) and litigated in row 32. Types by
declaration, values by execution, conformance by comparison.

#### Contracts are functions of the type, not of the instance

**Rule.** A leaf's [contract](#g-contract) declarations — `input_types`,
`output_types`, `events`, and the shapes of `init_x`/`init_m` — must be
determined by the component's **type**, its type parameters included, and never
by its field *values*.

The value-discarding signature `input_types(::Engine, ::Type{T})` is the visible
form of the rule. The idiom for a contract that genuinely varies is the type
parameter, not the field: `SumJunction{Wrench, 3}` ([§6.2][s6-2]), `Or{N}`
([§13.7][s13-7]). Arity is spelled in the type, at the price [§6.2][s6-2] states
openly.

**Why.** The entry typing decides it ([§12.7][s12-7]). A component's
[bundle](#g-bundle) is the `NamedTuple` of zero-copy views a component function
receives, and its key set *is* its contract's. An entry of the
[executor](#g-executor), the compiled execution form of the schedule, carries
what selects code in type parameters and what is plain data in fields. A key set
derivable only from field values would therefore have to go one of two ways. It
could climb into the type parameters anyway, multiplying specialization and
changing the cost model ([§12.7][s12-7]) of [chunking](#g-chunking), the
splitting of a large phase body into statically typed chunks. Or it could sit in
fields, dissolving the static typing that the zero runtime graph logic
([§5.1][s5-1]), the allocation invariant ([§7.5][s7-5]) and the fold-away
conformance test ([§12.5][s12-5]) all rest on.

The build reads each declaration once, against the concrete instance, so a
value-dependent contract does not announce itself. This is a rule authors keep,
not a check the build can run.

**`workspace` is the one exception**, and explicitly so. It is the by-allocation
[register](#g-register) (row 77): an allocator the framework *calls*, not a
schema it *walks*. It legitimately takes sizes from the instance
(`workspace(c::KF, ::Type{T})` reads `c.n`, [§7.3][s7-3]), because no entry type
is derived from it.

### 11.2 The declaration inventory

One continuous primitive, declared end to end:

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

#events: ordered and named — order is load-bearing (§5.3, §8.6); detection policy by the guard's return type (§2.1)
events(::Engine) = (
    start    = Event(start_guard, start_handler),        # boundary-detected: Bool guard
    ignition = Event(ignition_guard, ignition_handler),  # boundary-detected: Bool guard
    flameout = Event(flameout_guard, flameout_handler))  # localized: sign-form guard
start_guard(::Engine, (; m, u)) =                        #manual trigger: an input (§10.5)
    m.phase === off && u.starter
start_handler(::Engine, _) = (; m = (; phase = starting))          #no `x` key: no reset
ignition_guard(eng::Engine, (; x, m, u)) =                #predicate form
    m.phase === starting && x.ω > eng.ω_idle && u.fuel_available
ignition_handler(::Engine, _) = (; m = (; phase = running))
flameout_guard(eng::Engine, (; x)) = eng.ω_min - x.ω      #continuous form: localizable
flameout_handler(::Engine, _) = (; m = (; phase = off))
```

The blocks below take that inventory declaration by declaration, and record
where each schema fact gets its authority.

#### State, modes, discrete state

**Rule.** `init_x` on either [tier](#g-tier), and `init_m`, declare *by initial
value*: the type is derived from the value.

There is consequently no second artifact to drift and no separate type
declaration to check. The [workspace](#g-workspace) (component-declared mutable
scratch arriving as the `ws` bundle field) is the exception to that
[register](#g-register): it is declared *by allocation*,
`workspace(::C, ::Type{T})` on the continuous tier and `workspace(::C)` on the
discrete one, the method itself being the allocator. A workspace earns the
exception because it is not memory and none of the by-value arguments below
cover it ([§7.3][s7-3]).

This is the [boundary](#g-boundary) of legitimate derivation: deriving from
another declaration is sound, deriving from evaluated user code is not.
Declaring types here too, `input_types`-style, with `probe_value`
([§12.3][s12-3]) synthesizing the initial values, was rejected (row 73).

**Why.** The declared values are the base layer of the
[condition](#g-condition) substrate — a condition being the path-addressed
sparse overlay that sets a build's state. The overlays ([§14.1][s14-1]) fall
back to them leaf by leaf, and the compiled store writers bake
`merge(defaults, overlay)`, so there must be an authored value under every
leaf.

The asymmetry against `input_types`/`output_types` is one of kind, not style.
[Contracts](#g-contract) describe table [cells](#g-cell), recomputed from
scratch every [sweep](#g-sweep), and so need only types. `init_*` describe
[stores](#g-store) — the model's memory, which must have contents before the
first sweep can run.

**These declarations stay one-argument**, and the criterion is the register
they live in (row 166). It is stated once here, and the blocks below refer back
to it. A *by-value* declaration states nominal physics, and its *types*
[walk by rule](#g-leaf-walk) — the derivation of per-activation types from a
declared nominal type. [§7.1][s7-1] forces every state leaf to follow the
[activation](#g-activation) scalar (a re-run of Stratum C at a given scalar
type), so a `T` in the signature would record no choice its author could make.
Partials enter through per-invocation seeding, never through initialization. A
*by-type* declaration is a function of the activation scalar, which is why
`input_types` and `output_types` both take it on the continuous tier. A
*by-allocation* declaration takes the scalar too, `workspace(c, T)` being the
standing precedent (row 77). The criterion, not uniformity, is the rule: a `T`
in a signature means a choice was made there.

#### `input_types(::C, ::Type{T})` — continuous; `input_types(::C)` — discrete

An `input_types` declaration is a bare `NamedTuple` of types: zero framework
vocabulary, no wrapper types. On **continuous consumers the two-argument form
is mandated**, on **discrete consumers the plain one**. That is the same
[tier](#g-tier) mandate `output_types` carries. The [class](#g-class) — a
component's primitive-vs-assembly status — is read off declaration shape, and
the class fixes the form the declaration must take ([§11.5][s11-5]); either
violation is `TierSignatureMismatch`.

Entries are **[face](#g-face) bounds, not [cell](#g-cell) types**, and the
reading is **permissive** (row 167): an entry states, per leaf, what the
consumer *allows* to arrive there. Entries come in three forms:

| entry | the leaf is | what may lawfully arrive |
|---|---|---|
| `T`, alone or as a type parameter (`SVector{3, T}`, `RQuat{T}`) | **tolerant** | the [activation](#g-activation) scalar or a frozen `Float64` |
| `Float64` | **demanding frozen** | never partials |
| `Int`/`Bool`/enum leaves, abstract reference-typed entries | as it always was | what the declared bound admits |

A `T` entry is what a promoting consumer writes, and it is the overwhelmingly
common case: a walking producer, a frozen discrete producer and a root
[slot](#g-slot) are all admissible behind it, so substitution stays intact.

A `Float64` entry is the **FFI door**: this input must never carry partials. A
[component](#g-component) whose internals cannot propagate `Dual`s — an opaque
wrapper, a C table, a hand-rolled solver — declares it, and its
AD-incompatibility becomes schema-visible instead of folklore. The failure then
moves from a `MethodError` inside user math at the first `Dual`
[probe](#g-probe) to a named wiring error at build ([§6.1][s6-1]).

`Int`/`Bool`/enum leaves and abstract reference-typed entries stand as they
always were. [Abstract entries](#g-abstract-entry) state **structural
substitutability**: several concrete producer types admissible behind one
stable face. The field handles ([§4.4][s4-4]) are the demonstrated client —
`terrain = AbstractTerrainField`. They are spelled without `T`, being
references rather than numbers, and they are still never the tool for eltype
genericity: that is exactly what a `T` entry is, a promoting consumer writing
`SVector{3, T}` rather than an abstract bound.

Names-only [contracts](#g-contract) were rejected (row 33). Inputs are the
component's *requirements*, and only against them are the unconnected-input
error ([§6.1][s6-1]), over-wiring detection and
[did-you-mean](#g-did-you-mean) typo messages definable at all. A did-you-mean
message is the offending name plus the list-in-hand it should have matched.

**Two clauses check a wire** ([§6.1][s6-1]). The **nominal bound check** is
stated over evaluations: the producer's declaration at `Float64` must be `<:`
the entry at `Float64`. It is one uniform rule, degenerating to exact equality
for a concrete entry, concrete types being final. Beside it sits the
**tier-scoped walk-compatibility clause**: for a *continuous* consumer, a
walking producer leaf — the producer declared `T` there — requires a `T` entry,
while a [pinned](#g-walked) producer leaf satisfies either, frozen values
embedding upward. Both sides are declaration functions of `T`, so the clause is
decidable in [Stratum](#g-stratum) A (one of the build's three phases:
structure, schedule, activation) by evaluating them at a marker scalar. No user
stage code runs ([§12.1][s12-1]), and a violation is
`WalkingFaceAtFrozenEntry`.

**Discrete consumers take the bound check only**, and that scope is
load-bearing rather than tidy.

**Why.** A discrete stage reads exclusively at real [ticks](#g-tick) in the
nominal world, and a `Dual`-carrying cell exists only inside activations
discrete stages never run in ([§12.4][s12-4]). Continuous → discrete wires are
therefore unconditionally legal. The unscoped variant is rejected in row 167.

Because entries are bounds, nothing is ever "overwritten": cell types are
single-sourced from the producer side per activation ([§12.1][s12-1]), and a
`Dual`-carrying cell behind a `T` entry is the design working, not a promise
broken. The code-level complement is the **genericity obligation** — whatever
scalars the wiring delivers, the consumer's math promotes — still checked by
the `Dual` probe, never declared, and **scoped to the `T`-entries**. A
`Float64`-entry input imposes no such obligation, that being its point. So
**declarations record choices; obligations are checked**: the `T` entry records
the tolerance choice, the probe checks the promotion.

**The permissive reading is the operative one, and the two readings it escapes
are rejected** (rows 33, 54, 167). The *predictive* reading has the entry
saying what *will* arrive; the *envelope* reading has it as a promise to
promote. The permissive reading predicts nothing, and it is not constant —
pinned entries are rare but real — which is what makes the `T` carry
information here.

**Root slots are the one place an entry types a cell.** Produced by no
component, a slot has only the consumer declaration to take a type from. The
**slot type** is the entry evaluated at `Float64`, and only a *tight* bound
determines one. A face surfacing as a root slot must therefore resolve to a
concrete declaration, which [staging cells](#g-staging-cell), the
[trace header](#g-trace-header) and `probe_value` all need. Abstract-at-root is
a build error, and `AbstractAtRoot` names the remedy with the face: wire a
concrete producer, in a test rig a stub child ([§13.7][s13-7]). Under fan-out
the slot type is the unique concrete declaration among its consumers — two
different concrete declarations remain an error — and abstract co-consumers are
checked against it. The **slot cells** at an activation follow the slot type by
evaluating that same entry at the activation's `T`, which makes **seedability
schema-visible**: a `T`-entry slot is a lawful linearization `B`-matrix tap,
and a `Float64`-entry slot is *declaredly* unseedable ([§14.10][s14-10]).

**Fan-out combines tolerance by a meet, not by agreement** (row 168): the slot
pins at every activation if *any* consumer's entry pins, and follows the scalar
only when every consumer tolerates. Two consumers of one slot may agree at
nominal and still differ in tolerance. `SVector{3, T}` and
`SVector{3, Float64}` both evaluate to `SVector{3, Float64}`, so the slot
*type* is unambiguous while the entries disagree about partials. That mixture
is a legitimate model rather than a mistake: a command consumed by a promoting
aerodynamics leaf and by an AD-opaque table is the FFI door in use.

**Why.** The direction of the meet is forced by embedding. A pinned slot cell
feeds a `T` entry lawfully, frozen values embedding upward as zero-partial
constants ([§12.5][s12-5]), whereas a `Dual`-carrying cell arriving at a
`Float64` entry is precisely what that entry forbids. The meet is therefore the
only assignment satisfying every consumer at once, the mirror of the
walk-compatibility clause ([§6.1][s6-1]) on the producer side.

What the mixture costs is stated where it is paid: such a slot is unseedable,
and a tap selecting it is rejected naming the *pinning consumer* rather than
the face alone ([§14.10][s14-10]).

#### `output_types(::C, ::Type{T})` — continuous; `output_types(::C)` — discrete

`output_types` declares the public [port](#g-port) [contract](#g-contract), and
declares it **by type**. It is the same species as `input_types`, carrying the
[activation](#g-activation) scalar in its signature on the same terms. Where
the input side is read permissively, though, this one is read **literally**: an
entry states what the [cell](#g-cell) *carries*, not what it tolerates.

On **continuous producers the two-argument form is mandated**, the spelling
being `output_types(::Engine, ::Type{T}) where {T <: Real} =
(M_shaft = T, P = T, ω = T)`. On **discrete producers the plain form is
mandated**, and it *is* the wholesale pinning of the discrete exemption
([§7.2][s7-2]) — spelled in the signature as well as enforced by
[tier](#g-tier). Class is read off declaration shape ([§11.5][s11-5]), the
class fixes the form the declaration must take, and `TierSignatureMismatch`
names a producer whose form and tier disagree in either direction.

Semantics are **literal**: the cell types at an activation are the declaration
*evaluated* at that activation's `T` — nothing [walked](#g-walked), nothing
inferred. Participation is therefore authored **per leaf** and legible on the
page:
- **`T`, alone or as a type parameter** (`SVector{3, T}`, `RQuat{T}`,
  `MyStruct{T}`) — the leaf **participates**: its cell carries the activation
  scalar. Value parameters are structure rather than number and never take it
  (`Ranged{T, -1, 1}`; the bounds are not scalars to re-type).
- **`Float64`** — the leaf is **deliberately [pinned](#g-walked)**, and the pin
  is schema-visible: whole-leaf freezing, declared and conformance-checked.
  That is the recorded freeze door ([§14.10][s14-10]) delivered — declare
  `Float64` and strip with `ForwardDiff.value` inside the stage, the
  stop-gradient stated in the contract instead of buried mid-expression.
- **`Int`/`Bool`/enum leaves and reference-typed fields** pin as they always
  did (a [§4.4][s4-4] bulk-data handle's grid is frozen build-time data, never
  activation-dependent).

The companion obligation is **constructibility at `T`**: a declared type must
be buildable at the activation scalar. The `Dual` [probe](#g-probe) enforces it
by construction — it builds real values, so a type whose constructor cannot
accept them detonates at the probe with its own name in the message.

During a generic [sweep](#g-sweep), gated-off discrete producers hold their
`Float64` values, consumers gather mixed tuples, and promotion does the rest.
That is semantically exact: a frozen discrete output is a constant with zero
partials, which is precisely what "linearize the continuous dynamics with the
discrete state held" means. The frozen cell is not an AD limitation on the
signal path; it is the true zero of an instantaneous dependence the hybrid
semantics never had (`frozen_discrete_walkthrough.md`). What makes the mixing
safe is the **embedding guarantee** ([§12.5][s12-5]), keyed on
**declared-`T` leaves** (row 33).

**Why.** A `Float64` observed at a declared-`T` leaf under a non-nominal
activation implies no `Dual` entered its computation: promotion is airtight,
and there is no lossy cast. Its true derivative along every seeded direction is
therefore zero, and embedding it as a zero-partial constant is exact.

Piecewise branches returning literal constants (`flow > 0 ? f(x) : 0.0`) are
legal as written, zero partials being the derivative of a locally-constant
branch. Which *invocation* carries partials is still chosen by seeding
([§14.10][s14-10]), never by typing: the declaration says which leaves *can*
carry them, the seed which directions do.

**The forgotten-`T` account, stated openly.** The whole-signature variant — a
continuous producer declared as though it were discrete — is unwritable by
construction: the tier mandate catches it in [Stratum](#g-stratum) A, before
any user code runs. What remains is per-leaf: an author writes `Float64` at a
leaf that really participates.

That bug **lurks, but is never silent**. No lossy `Dual → Float64` cast exists,
so the first `Dual` activation of that [component](#g-component) fails. It
fails at that activation's own lazy Stratum-C compile ([§12.4][s12-4]), not at
`build(world)`. The message carries the didactic hint ("if `F` participates in
differentiation, declare it `T`"), an observed `Dual` at a declared-pinned leaf
having exactly one honest cause.

The lurk is contained by policy rather than machinery: **the test suite builds
a `Dual` activation of every component**, which is the exhaustive set
[§12.4][s12-4] defines. An activation is a Stratum-C re-run, cheap enough to
make this unremarkable in CI. What the form buys in exchange is **reader
honesty**: participation is read off the declaration instead of reconstructed
from a framework rule carried in the reader's head, and a genuinely frozen leaf
can say so.

**The stores are walked; only the output side is evaluated.** The type derived
from a *continuous* leaf's `init_x` is walked: real leaves and `Real` type
parameters follow the activation scalar. `init_m` and a *discrete* leaf's
`init_x` pin wholesale, mirroring the discrete-producer rule. The asymmetry is
the register criterion stated above under the by-value declarations, not an
inconsistency: `init_*` declare *by value*, and [§7.1][s7-1] admits no pinned
state leaf for a `T` to record a choice about. Declared `Float64` initial
values embed as zero-partial constants under non-nominal activations, which is
the rule for `Float64` condition leaves ([§14.3][s14-3]) applied to the
defaults those conditions overlay.

Walking a continuous leaf's `init_x` presupposes the closed leaf vocabulary
[§7.1][s7-1] fixes: scalars and `SArray`s at the common eltype. On the discrete
tier, `init_x` keeps the full type freedom [§7.3][s7-3] allows. Stratum A
therefore checks the continuous vocabulary ([§12.1][s12-1]) and reports a
failure in the didactic register:
- "`init_x` field `gear_count::Int` is not a continuous state — integers,
  `Bool`s and enums belong in `init_m`";
- "`init_x` field `q_nb::RQuat` is not a state leaf — declare the `SVector{4}`
  backing and cast where rotation semantics are wanted ([§7.1][s7-1])".

#### `events(::C)`

`events` declares an ordered, named collection of [guard](#g-guard)/handler
pairs, spelled `Event(guard, handler)` with no detection keyword. Detection
policy is declared by the guard's return type instead: a `Bool` guard makes the
event [boundary-detected](#g-boundary-detected), checked for edges at step
boundaries only, with no root-finding. A guard returning the nominal scalar
makes it [localized](#g-localized), the crossing instant bracketed by
root-finding over probe sweeps ([§8.4][s8-4]). Order is semantics: declaration
order ([§5.3][s5-3]), priority with re-decision ([§8.6][s8-6]). Nothing here is
inferrable.

#### No stage tags anywhere

Which stage produces which [port](#g-port) stays invisible in the
[contract](#g-contract), preserving [§4.2][s4-2]: moving a port between stages
is non-breaking for consumers. Membership is *derived* instead, with no
chicken-and-egg. Stage-1 functions (`h_x`) structurally receive no inputs, so
the build [probes](#g-probe) them first, observes their contract ports, assigns
the remainder to stage 2, builds the graph, and probes the stage-2 chain in
topological order with real upstream values. The "decoder takes no inputs"
property is exactly what makes the derivation well-founded. Stage names carry
no [tier](#g-tier) information at all (row 173): what declares a leaf's tier is
the completeness rule below.

#### Custom structs as port types

A custom struct is a first-class port type — `contact = GearContact{T}` —
under the scoping [§7.2][s7-2] establishes: parametric in its real-scalar
leaves, constructors inferring the scalar, no [pinned](#g-walked) fields on the
continuous path. A participating struct leaf is declared with the scalar in its
parameter position, `GearContact{T}` recursively for nested parameters. A
struct with a hardcoded `Float64` field offers no such position, so it can only
be declared bare: a pinned leaf, honestly spelled. Any `Dual`-carrying
construction then detonates inside the stage with an `InexactError` naming the
offending constructor — the [§7.2][s7-2] CI invariant reached through the
declaration layer with no extra machinery.

#### Completeness of the declaration set

Four rules the build checks in [Stratum](#g-stratum) A ([§12.1][s12-1]), stated
here because they are properties of the declarations, not of the wiring.

**A store needs its update.** `init_x` with neither an `f` nor a `g` method is
a build error: continuous state with no [flow](#g-flow), or a discrete store
nothing updates. The framework will not silently supply `ẋ = 0`, which is a
model, not a default; and an unupdated discrete store is a parameter in
disguise, parameters being plain struct fields. The didactic register says
exactly that. `init_m` carries no such obligation: modes are written by
handlers, and a [component](#g-component) may legitimately declare modes no
event of its own transitions.

**An event needs both halves.** An `events` entry whose [guard](#g-guard) or
handler has no method for the component type is a build error, caught by method
lookup at declaration-reading time rather than as a `MethodError` at the first
firing. An event that fires only in a corner of the envelope would otherwise
hide the omission indefinitely.

**[Tier](#g-tier) is declared by the update law.** For a **stateful** leaf, `f`
marks continuous and `g` marks discrete: the flow/jump pair carries the tier,
the stage names being shared (row 173). The remaining tier-implying
declarations must agree. `init_m` and `events` are continuous-only, the event
system being continuous-side only ([§5.2][s5-2], [§3.2][s3-2],
[§14.1][s14-1]). `workspace`'s arity splits the tiers (`(::C, ::Type{T})`
versus `(::C)`), and so do the arities of `output_types` and `input_types`
(rows 166–167). Disagreement is `DeclarationOnWrongTier` ([Appendix C][sC]),
reported as the offending declaration with the tier the leaf's other
declarations announce; it covers declaring both `f` and `g`, and a `g` beside a
two-argument `output_types`.

A **stateless** leaf declares no store and no update law, so its tier is
decided by its [contract](#g-contract) arities: `output_types`, mandatory hence
always the decider, with `input_types` agreeing where declared. The arity is no
mere marker — it *is* the tier's semantics (rows 166–167): the two-argument
forms declare [cells](#g-cell) and tolerances at the
[activation](#g-activation) scalar, walking with it, where the plain forms
declare the [pinned](#g-walked) discrete world. A stateless `h_xu` component is
tier-transparent library material ([§13.7][s13-7]). Members of both families,
or of neither, are the [§11.5][s11-5] class errors.

**The root of a build is an [assembly](#g-assembly).** Root [slots](#g-slot)
are the root's input [faces](#g-face) declared through `input_connections`
([§6.1][s6-1], [§9.3][s9-3]), and only assemblies declare
[boundary](#g-boundary) connections ([§11.6][s11-6]). A primitive root
therefore has no root slots — its faces are just its own [port](#g-port) names
— and every input it declares is an unconnected-input error. Exercising a leaf
alone is what the [component test rig](#g-component-test-rig)
([§13.7][s13-7]) is for: it supplies the one-child assembly.

### 11.3 Visibility: the contract is the interface

**Rule.** Visibility is decided by *where the value goes*:

- a field declared in `output_types` is public;
- a field returned in `w`, the optional second slot of an output stage's
  return, is private *by construction*;
- a field returned in `y` — a stage's own published signals — and declared
  nowhere is a build error;
- a component with no `output_types()` method has no outputs.

That is the same move as class-by-declaration-shape.
[Ports](#g-port) in the [contract](#g-contract) are connectable, GUI-listed,
[snapshot](#g-snapshot)-carried and log-exported. The table is public
throughout, every [cell](#g-cell) a declared port or an auto-published one, so
nothing anywhere needs a presentation filter. Private intermediates are not
filtered — they are simply not there. `w` is handed from stage to consumer as a
value (the one-hop law, [§5.2][s5-2]), with no cell, no name in a contract and
nothing for a wire, a listing or a log to reach.

The inspection path for an intermediate is **promotion**: one line in
`output_types` makes it public, checked and visible everywhere at once
(row 165). FlightCore is the precedent, where an intermediate was inspected by
putting it in the `Model` output and no other way. Publicity is never implicit:
even the minimal [component](#g-component) writes
`output_types(::LowPassFilter, ::Type{T}) where {T <: Real} = (x = T,)`, one
line, in exchange for "public" always meaning someone wrote it down.

- **Conformance**: a declared port must be produced, by exactly one stage or by
  **auto-publication**. Auto-publication covers declared names matching state
  or mode fields that no stage produces ([§5.3][s5-3]). Stage membership is
  derived over `output_types` alone ([§12.1][s12-1]). Declared-but-unproduced
  and produced-by-two-stages are build errors. A declared port matching neither
  a stage product nor a state field errors with both lists in hand: "not
  produced by any stage and not a state field". A *returned port field*
  declared nowhere is a build error at [probe](#g-probe), with
  [did-you-mean](#g-did-you-mean) — the offending name plus the list-in-hand it
  should have matched — against `output_types`. That is the return-side
  analogue of [§11.4][s11-4] walkthrough 1 (rows 34, 55). The forgotten-branch
  walkthrough holds: a declared `P` missing from the taken branch's return
  fails at probe; missing from an *untaken* branch, it fails loudly at that
  branch's first execution via the always-on check.
- **`w`'s regime is probe-observed, and that is sound precisely because `w` is
  not a cell.** A cell needs a fixed type per [activation](#g-activation) (a
  re-run of Stratum C at a given scalar type), which only a declaration can
  supply. A value flowing between two functions in one fused pass has no type
  contract to violate, and mixed branches are handled exactly by promotion. So
  the probe takes `w` as it finds it: it checks that the second return
  [slot](#g-slot) is a `NamedTuple` at all, and it checks the *consumer's*
  reads against the observed field set. A destructured name that is not there
  fails inside the framing diagnostic ([§13.2][s13-2]) with did-you-mean from
  the actual fields. That is weaker than a declaration-backed message: it can
  say "`f` of `Foo` reads `w.q_dny`; the producing stage returned `q_dyn`" but
  cannot say which spelling was intended. It is located, name-shaped, and costs
  no declaration.
- **Branch-shape rule**: stage returns must have the same `NamedTuple` shape on
  every branch. Julia's type-stability discipline already demands that for
  performance; the framework merely makes it a stated rule with a good error.
  `w` is inside the rule: a stage whose `w` changes shape between branches is
  as broken as one whose `y` does.
- **The always-on check covers `w` at the nominal activation only.** Beside the
  `y` test ([§12.5][s12-5]) sits one baked `isa` against the type the *nominal
  probe observed*. It folds to nothing while the stage is stable, and it
  converts the unintended-branch-divergence class — which otherwise announces
  itself only as an allocation in somebody's benchmark — into a loud located
  field-naming error at the divergent branch's first execution. The blame text
  says what it honestly knows: "expected what the nominal probe observed". This
  complements the canary ([§7.5][s7-5]): the canary detects, this localizes.
  **Non-nominal activations run no `w` check at all**, and deliberately. There
  is no declaration to anchor a branch-independent expectation, and no store
  whose typing the check would be protecting. A strict probe-observed
  expectation would also fire on the legal constant-branch idiom — the
  one-probe-point argument, which is exactly what kills observation authority
  for cells. Correctness needs no [guard](#g-guard) there: a `Float64` in `w`
  under a `Dual` activation is an honest zero-partial constant by the embedding
  guarantee ([§12.5][s12-5]), and its downstream promotion is exact. Walking
  the nominal observation to synthesize non-nominal expectations was rejected
  (row 165).
- **[Schema authority](#g-schema-authority) is total over the table**
  (declarations define structure; evaluation only checks conformance): every
  *cell* traces to an authored declaration, the always-on check's expected type
  for `y` is fully declaration-derived, and return typos cannot silently define
  new cells. Protection against silently dropped partials rests on the
  embedding guarantee — promotion is airtight, so an observed `Float64` is a
  true constant, [§12.5][s12-5]. Probe-observed expected types remain rejected
  *for cells* (rows 34, 55, 165), and that rejection does not reach `w`, which
  declares nothing and types nothing.
- **What this rules out** (rows 16, 34, 55, 165): the `unlisted` flag
  ([§4.2][s4-2]) and its satellite-function representation; identity
  publication by default ([§7.4][s7-4] step 4); **probe-observed private
  cells**; the `Private(T)` fallback; and the opt-in variant with a
  `Float64`-under-`Dual` diagnostic.

### 11.4 Failure walkthroughs (the error-locality grounding)

The five mistakes that decided declaration-vs-inference, with their failure sites
under this layer. Each was traced under inference-by-evaluation too, and in every
case the failure surfaced inside *correct* code, later, or never; row 32 carries
the traces.

1. **Typo'd wire** (`:throtle`): build error at the connection, "no input
   `throtle`; did you mean `throttle`?"
2. **Forgotten wire** (`fuel_available`, read only by a [guard](#g-guard)): [§6.1][s6-1]
   unconnected-input error at build.
3. **Forgotten branch field** (`P` returned by one branch only): [probe](#g-probe) or
   first-execution error naming the declared [port](#g-port).
4. **Type mismatch** (a `Float64` fraction wired into a `Bool` input): wiring-time
   error naming both endpoints and both [faces](#g-face).
5. **Typo'd return field**, in its two [registers](#g-register). A typo'd *port*
   (`P_shft = ...` for a declared `P_shaft`) keeps the full strength of the
   declaration: a probe error with [did-you-mean](#g-did-you-mean) (the offending name plus the
   list-in-hand it should have matched) against `output_types`, plus the
   unproduced-`P_shaft` error with both the stage-product and state-field lists
   in hand. A typo *inside* `w` (`q_dny` where the consumer reads `q_dyn`) has
   no declaration to be checked against. It surfaces one hop later, at the
   consumer, as the framing diagnostic ([§13.2][s13-2]). That diagnostic
   carries the producing stage's observed field set: "`f` of `Foo` reads
   `w.q_dyn`; the stage returned `q_dny`". It is weaker than a
   declaration-backed error, since the framework cannot know which of the two
   spellings was meant. Even so, it is located at the pair of lines that
   disagree, and it is still name-shaped. That is the price of the private
   channel, paid where no interface is at stake. The remedy for an intermediate
   worth stronger checking is the one [§11.3][s11-3] names: declare it an
   output.

### 11.5 Assembly declaration: type-based, class by declaration shape

**Rule.** An [assembly](#g-assembly) is a plain struct: fields whose type is
`<: AbstractComponent` are its children, and all other fields are inert
parameters.

Field names are path segments. Substitutability and variants use ordinary
parametric fields — exactly today's `Cessna172X{K, A}` shape. Alongside the
struct come the well-known declarations: `child_connections(::A)`, mandatory
even when empty, plus `input_connections(::A)`, `output_connections(::A)` and
`sample_times(::A)`.

#### Container children

**Rule.** A field whose type is a `Tuple` or `NamedTuple` with *every* element
`<: AbstractComponent` contributes its elements as
[container children](#g-container-children).

They are path-named `"field/1"…"field/N"` (tuples) or `"field/key"`
(NamedTuples), declaration order governing layout. Containers are **transparent
grouping, not assemblies**: no [contract](#g-contract), no `child_connections`,
no [rate scope](#g-rate-scope), no existence beyond the path segment. The
elements are children *of the parent*, whose `child_connections`/
`input_connections`/`output_connections`/`sample_times` address them by element
name. Anything wanting its own wiring or [faces](#g-face) declares itself an
assembly.

The payoff is parametric composition. `struct Formation{NT <: NamedTuple};
aircraft::NT; … end` holds any [roster](#g-roster) — size, names, mixed
aircraft types — per instantiation, the declaration bodies generating wires by
comprehension over the keys. That is the arity-via-computed-contracts pattern
[§6.2][s6-2] uses for `SumJunction{W, N}`, here at structure scale. The swarm
worlds ([§14.9][s14-9]) consume it directly, and so does
[mounting](#g-mounting), the relocation of a whole problem or tap set with
[`at`](#g-at)`("aircraft/red", problem)`.

The edges of the container form are fixed by rule:

- A container mixing [component](#g-component) and non-component elements is a
  build error in this section's [did-you-mean](#g-did-you-mean) family — the
  offending name plus the list-in-hand it should have matched. All-component
  elements are children; zero-component elements are inert parameter data.
- Containers of containers are rejected in the first cut, deeper grouping being
  what assemblies are for.
- Empty containers are legal, contributing zero children: parametric code then
  needs no special case.
- Abstract element types follow the same concreteness discipline as plain
  fields — directly concrete, or concrete through type-parameter bounds. That
  is the [generic holding](#g-generic-holding) — a parent holding a child
  through a non-concrete field type — that [§11.8][s11-8] allows.

`sample_times` needs no rule change. Element names are immediate child names,
hence legal keys, and the bare field name is sugar for a uniform declaration
across all elements.

#### The builder is rejected

The builder — `Assembly()` plus `add!`/`connect!` — is rejected (row 39).

Its one real advantage, programmatic generation, survives intact in the
type-based form: a declaration is an ordinary function body, and loops and
comprehensions build the returned tuple.

#### `Group`: the on-the-fly assembly

The *immutable* version of "grouping components by plain calls" needs no
builder. It is already expressible under this section's rules as a single
library component (the starting inventory, [§13.7][s13-7]): the
container-children rule makes a `NamedTuple` field contribute its elements as
path-named children, and declarations are ordinary functions of the *instance*,
free to read its fields:

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

One type, defined once; every ad-hoc topology is a *value* of it. The type
parameters still carry the children's concrete types, so [Stratum](#g-stratum)
C specialization is unchanged (the strata being the build's three phases:
structure, schedule, activation). So is the [executor](#g-executor), the
compiled execution form of the schedule ([§12.7][s12-7]). Wiring validation,
did-you-mean errors and the two-producer check all run at build against the
instance exactly as for a named assembly.

What is given up relative to a named type is exactly what named types are
*for*: dispatching domain code on `::Cessna172X`, a reusable identity for the
topology. The exploratory and programmatic composition `Group` serves does not
want it anyway.

The reach of the builder rejection is fixed by row 184: it targets mutable
recipes, not type-based *semantics*. `Group` is the library's anonymous
assembly form beside the named types, shipped the way Julia ships anonymous
functions alongside named ones, and it serves the model assembler with a
library addition and zero new declaration rules.

#### Class by declaration shape

**No `AbstractAssembly`; one root `AbstractComponent`** (row 39).

**Why.** The domain hierarchies — `AbstractAircraft`, the engine families —
have to carry both classes: a [slot](#g-slot) declared `E <: AbstractEngine`
must accept a primitive `PistonEngine` and a composite turbofan assembly alike.
And class is implementation detail behind the contract ([§11.3][s11-3]).

[Class](#g-class) — a component's primitive-vs-assembly status — is declared
instead by *which* well-known declarations a type defines. `child_connections`
is the marker, mandatory even when empty (the `LowPassFilter` precedent), and
defining it makes an **assembly**. Any leaf declaration makes a **primitive**:
`init_x`/`init_m`, `workspace`, `input_types`/`output_types`, `events`, or any
stage, `f`, `g` or `project` method.

The rule is total. A `<: AbstractComponent` type declaring neither family has
no class to read, and is a build error naming both families rather than a
silence that fails later and elsewhere. That error sharpens into a did-you-mean
when the type has component-typed fields ("holds components but declares no
`child_connections`"). `child_connections` plus `init_x` or a stage on one type
is a build error as well: assemblies have no state of their own —
no-atomic-assemblies at declaration time ([§8.5][s8-5]).

Reading which declarations exist is reading declarations — the same move as
visibility-by-declaration-site ([§11.3][s11-3]), not the banned
inference-by-evaluation ([§11.1][s11-1]).

#### Contract signature shape follows the class

Class also **mandates the shape of the contract signatures** rather than merely
being read from them (rows 166, 167). **Both** contract declarations follow the
[tier](#g-tier): on a continuous leaf, `input_types` and `output_types` must
take the two-argument form `input_types(::C, ::Type{T}) where {T <: Real}` and
`output_types(::C, ::Type{T}) where {T <: Real}`; on a discrete leaf, both must
take the plain one-argument form.

Either violation — a continuous declaration missing the `T`-form, a discrete
declaration carrying one — is `TierSignatureMismatch` ([Appendix C][sC]). The
diagnostic reports the component path, the declaration at fault, the tier its
other declarations announce, and the form found versus the form mandated. The
check is Stratum A and collected: declaration shape is read, nothing is
evaluated.

The tier fact is therefore spelled in the signature *and* fixed by the class,
the two kept in agreement by a check rather than by convention. That is what
makes the whole-signature forgotten-`T` bug (the worst case, row 79)
unwritable.

### 11.6 Paths, wiring and faces

**Paths are slash-separated strings** — `"systems/ldg/left/trn"` — relative to the
[assembly](#g-assembly) being declared, no leading slash; one canonical form, shared verbatim by
declarations, error messages, [device](#g-device)/[trace](#g-trace) addressing ([§9.3][s9-3]) and the HDF5 log
tree. [Container children](#g-container-children) ([§11.5][s11-5]) add index and key segments — `"aircraft/2"`,
`"aircraft/red"` — ordinary segments, resolved against the container field. Instance navigation,
tuples of symbols and dotted paths were all rejected (row 40); a path-tracking
proxy remains addable sugar. One fact from that adjudication is load-bearing
downstream: symmetric immutable siblings are `===`-identical, so a path is
unrecoverable from an instance — which is why the helpers ([§11.8][s11-8]) name the child
by path.

**`child_connections(::A)`** is an ordered collection of `"src/port" => "dst/port"`
pairs, strictly child-port → child-port; the rules ([§6.1][s6-1]) apply (one wire per input,
deep paths through concretely-typed fields only, stopping at a generic child's
[faces](#g-face)). The assembly's **[boundary](#g-boundary)** is declared by two further methods, one per
direction. **`input_connections(::A)`** is an ordered collection of pairs, face
name => internal endpoint path — or a tuple of paths for an input face routed to
several internal endpoints (fan-out through the boundary)
(`"trn" => ("systems/ldg/left/trn", …)`). **`output_connections(::A)`** runs the
other way, internal source path => face name
(`"aircraft/pose" => "view_pose"`), so that its pairs, like every other pair in
the three declarations, read along the flow.

**Face names are arbitrary strings with two build-checked invariants.** The
first is that a face name contains no `/` (reserved for structural paths). The
second is uniqueness across the union of the two boundary declarations' face
names. Every other naming choice (separators, grouping prefixes like
`"pilot.throttle_axis"`) is author convention, not framework law. The
`input_passthrough` helper's defaults ([§11.8][s11-8]) document the house style
without legislating it.

The two-notation rule this rests on is directional — structure vs. contract, not
read vs. write. **Slash is structure**: endpoint paths walking real children and
ports; the inspection [register](#g-register)'s [snapshot](#g-snapshot) and log
addressing. **Face names are opaque contract tokens.** The
[periphery](#g-periphery)'s write side (input devices, mappings, the trace, the
GUI write path) speaks face names exclusively ([§9.3][s9-3]). The read side
speaks them wherever it wants meaning that outlives the build: integration
bindings (`get_face`, [§9.2][s9-2]) and load-bearing service reads
([§14.4][s14-4]).
The three declarations return pairs of strings rather than NamedTuples (row 46).

One invariant spans all three declarations: every pair's arrow points the way the
signal flows — the left side is a producer or entry point, the right side a
consumer — and every right side is fed exactly once. **Direction is therefore
declared by the method**, not inferred: the resolved endpoints only *cross-check*
it, and an entry whose endpoint resolves to a port of the wrong direction is a
build error naming the method, the entry and the resolved port's actual direction.
A mixed entry is not expressible: the single list that made that error class
possible does not exist; two entries producing the same output face
remain the ordinary two-producers error. Face *types and [tiers](#g-tier)* are derived from
the internal endpoints — the [blessed](#g-blessed) derivation-from-declarations ([§11.2][s11-2]) — and the
derivation is forced, not merely convenient (row 41): an assembly is
tier-neutral, exporting continuous-sourced and discrete-sourced ports side by
side, and a face's [cells](#g-cell) follow the producer's own declaration
([§11.5][s11-5]), evaluated at the [activation](#g-activation) scalar on the
continuous tier and [pinned](#g-walked) on the discrete.
Three alternative spellings are rejected (rows 41, 170): routing values under
the leaf names `input_types`/`output_types`, leaf-style *typed* faces with face
wires inside `child_connections`, and routing-as-wires with derived types and no
face list. Publicity is never implicit ([§11.3][s11-3]).

**Root [slots](#g-slot) fall out with no vocabulary**: at every non-root level an input face
declared through `input_connections` is fed by the parent's wire; at the root there
is no parent, and the face *is* the [write surface](#g-write-surface), the set
of faces a writer's batch entries may reach ([§9.3][s9-3]). The whole-tree
obligation model ([§6.1][s6-1]) states the complementary error rule. An
assembly never declares its external connections — those live in the parent
that instantiates it, exactly as a leaf's do.

**A [worked](#g-worked) assembly.** The IMU ([§15.5][s15-5]), spelled in full — a mixed-tier assembly
exercising paths, faces and sample times together:

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

sample_times(::IMU) = (sampler = Relative(1), errors = Relative(1))
```

Two spellings worth reading closely: `input_passthrough` enumerates the child's
**input** faces and nothing else ([§11.8][s11-8]), which is why the pass-through of the
integrals' kinematic-truth inputs (`q_eb`, `r_eb_e`, `ω_eb_b`, `a_ib_b`,
`α_ib_b`, [§15.5][s15-5]) is a bare splat with nothing to say about direction;
and the measured-increment face sources `errors/sample_meas`, the error
model's *output* port, not the `errors/sample` input the sampler already
feeds — listing `errors/sample` in `output_connections` would fail the direction
cross-check, and listing it in `input_connections` while it is wired is the
two-producers error of [§11.8][s11-8].

Three facts the example carries: the assembly is tier-neutral — every face's
type and tier derive from its internal endpoint, and a `sample_times` key on
`integrals`, the continuous child, would be a [§11.7][s11-7] build error; the two
discrete children default to `Relative(1)` anyway, so this `sample_times`
declaration is declaratory — their absolute rate arrives from the enclosing
scope at deployment ([§11.7][s11-7]); and the latch-back wire ([§15.5][s15-5]), where the integrals consume
the sampler's published latch, joins `child_connections` as one more ordinary pair.

### 11.7 Rate scopes

`sample_times(::A) = (nav = Relative(5), gnss = Absolute(Hz(10)))` — child name =>
`Relative` or `Absolute` declaration ([§8.5][s8-5]'s two registers: relative entries
composing affinely down the tree, absolute entries anchoring; all compiled to one
`(D, Φ)` pair per discrete [component](#g-component)). The wrappers are the whole value
vocabulary — a bare integer or bare quantity is a declaration error. The
declaration is optional, and so is any given key: an unlisted discrete child
defaults to `Relative(1)`, so only multiplied, phased or anchored children need
appear. Keys are **immediate child names only** — a deep key would edit another
type's design from outside, and the composition rule guarantees you never need
to. Container elements ([§11.5][s11-5]) are immediate children, so `"aircraft/red"` is a
legal key; the bare field name applies one declaration to every element. A
`sample_times` key on a continuous child is a build error (the Δt-on-continuous
error at declaration time, [§8.5][s8-5]). `Δt_base`, `h` and `n` appear in no
declaration — they are deployment decisions fixed at `Simulation` construction
(the three sources for `Δt_base`, [§12.1][s12-1]). The declaration belongs to the
[assembly](#g-assembly) type, not to the child instance: a sample time is a design ratio
or a modeled [device](#g-device)'s intrinsic rate ([§8.5][s8-5]), never a per-instance value.
The FlightCore-`Subsampled`-style instance wrapper is rejected in row 42.

### 11.8 Computed connections and generic boundaries

`input_connections` and `output_connections` are ordinary functions evaluated at
build against the concrete
instance, so they may *compute* entries from child [contracts](#g-contract) — derivation from
declarations, [§11.2][s11-2]-blessed. The framework helper, sketched:

```julia
# the two shapes of `declaration_error` used below:
declaration_error(path::AbstractString, why::Symbol)      # e.g. :both_given
declaration_error(path::AbstractString, unknown, legal)   # did-you-mean against the legal set

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

The child is named by path and never passed as an instance, because the `===`
problem ([§11.6][s11-6]) makes a path unrecoverable from an instance. A
[face](#g-face) name containing dots is a legal final path segment on the
internal-endpoint side, precisely because slash is the only structural
separator. Computed entries mix freely with hand-written ones in either
declaration. `resolve` and `input_faces` are build-pipeline primitives needed
anyway, and `input_passthrough` is a thin composition. That is what keeps the
helper sugar rather than machinery. There is no `rename` hook,
because the [boundary](#g-boundary) declarations are ordinary code (map over
the pairs). Normative signatures for both primitives are in [§13.3][s13-3].
Every error stays first-class: an `except` face the [assembly](#g-assembly)
then fails to wire is an ordinary unconnected input; a face both wired and
passed through is a two-producers error; `except`/`only` naming a nonexistent
face errors with the child's face list in hand; a `prefix = ""` collision is
caught by the build's uniqueness check like any hand-written duplicate. The
effective face list is plain printable data — the inspectable contract of this
instantiation. What computation does *not* do is
auto-bubble: the author wrote down "every input face of this child that I don't
feed, I expose under this prefix" — explicit at the type level, evaluated at
build.

**The name carries the direction.** `input_passthrough` reads
`input_faces(child)` and `except`/`only` filter *face names* within that set;
the helper exists for the pass-through case, where an assembly hands a child's
unfed requirements up one level. Computed *output* re-export — a sibling
`output_passthrough` splatted into `output_connections`, with the predicate
selection [§13.7][s13-7]
records for `except`/`only` — is a **[guarded addition](#g-guarded-addition)**: plausibly wanted by the
conventional-surface work ([§9.2][s9-2], [§16][s16]) and by test rigs, cheap to add, and not
adopted here, because every output face in the worked assemblies is an explicit
pair and no consumer has yet demonstrated the computed form. It could not be a
keyword on one helper in any case: after the boundary split a single call cannot
emit entries into two different declarations.

**One authored list, two declarations.** The `World` example's two-entry
`except` understates the real shape. Every level of a realistic tree is a
generic [seam](#g-seam), and an assembly that feeds some of a child's input faces while
passing the rest up must name the fed ones in `except` — at C172X scale, four
seams and roughly ten names at the innermost one, restating in each `except`
tuple the wire list sitting in the same assembly's `child_connections`. That is
"structure kept in two artifacts" ([§11.1][s11-1]; row 39), the shape this
design refuses
elsewhere. It needs no vocabulary: declaration bodies are ordinary code
([§11.5][s11-5]), so
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
[projections](#g-projection) of the authored list, so the drift class is removed rather than
detected. Every misspelling stays loud: a mistyped destination is an
unknown-face error with the child's face list in hand, whether the wire or
the `except` entry meets it first. One honest asymmetry: a pair *omitted*
from the list is not an error but a structural change — the face leaves the
`except` set and joins the input face surface, ultimately a root [slot](#g-slot) for
conditions to cover
([§14.6][s14-6]).
What the idiom preserves, and the helper below surrenders, is that the feed
statement exists to be reviewed: an omission is legible in one authored
artifact, not defined away as the complement of the wire list.

**The line not to cross** is deriving `except` from `child_connections` itself — a
helper spelled `except = fed(s, "aero")`, reading the assembly's own wire
list. That is auto-bubbling under another name (rows 43, 145). The single source
must be **authored data, never inferred structure**.

**[Generic holding](#g-generic-holding) = imposed contract.** A parent holding a child generically
constrains it exactly through the faces its wires and boundary connections reference: build a
`World` whose concrete aircraft lacks a referenced face and the error names the
`World` entry — build-time structural typing, no new vocabulary (a formal
required-faces declaration on domain abstract types remains possible sugar).
Scalar faces make partial scripting compose: a guidance [scenario component](#g-scenario-component) wires
`mode_req` and `EAS_ref` while the remaining faces stay exported for GUI or
defaults — impossible with a bundled face ([§4.3][s4-3] write-side corollary).

---

## 12. The build pipeline

The build consumes a root [component](#g-component) instance and produces the runnable artifact:
resolved wires, typed [signal table](#g-signal-table), evaluation [schedule](#g-schedule), absolute rate divisors,
flat state layout, root [slots](#g-slot). [§11][s11] states what is declared and what must hold;
this section states *when* each fact is checked, against what, and with which
failure. The [§11.4][s11-4] walkthroughs plus the error rules ([§6.1][s6-1]) are its acceptance tests.
Error-*reporting* policy is settled in [§13.1][s13-1]: declarative checking passes
collect, user-code evaluation fails fast, strata are barriers — the only partial
results carried past failures are violation lists from pure checks.

### 12.1 Three strata

Three ordering constraints are forced by settled decisions: [face](#g-face)
derivation is **bottom-up** (an [assembly](#g-assembly)'s
[boundary](#g-boundary) connections evaluate against child
[contracts](#g-contract), [§11.8][s11-8]); the unconnected-input obligation
check and cross-level two-producers detection are **global** — decidable only
at the root, after every assembly's wires and faces are in hand ([§6.1][s6-1]);
and stage membership is **derived by probing** the stage-1 functions
([§11.2][s11-2]), so evaluation interleaves with graph construction at exactly
one [blessed](#g-blessed) spot. The pipeline is therefore inherently
heterogeneous, organized as three [strata](#g-stratum).

#### Stratum A — structure

Stratum A is pure declaration reading: no user stage code executes in it. The
`input_connections`/`output_connections`/`input_passthrough` bodies are
declaration code ([§11.8][s11-8]).

The stratum is a tree walk from the root instance, in this order:

1. [Components](#g-component) are collected by path.
2. Each component's [class](#g-class) — its primitive-vs-assembly status — is
   read off declaration shape ([§11.5][s11-5]).
3. Leaf contracts are collected: `input_types`, `output_types`, `init_*`
   values, `events`.
4. Face derivation runs bottom-up.
5. Global wiring resolution then runs, resolving wires to absolute leaf
   terminals.

Resolution runs these checks:

- one-writer-per-input;
- the typo [did-you-mean](#g-did-you-mean) — the offending name plus the
  list-in-hand it should have matched — against the destination's input list;
- the two wiring type clauses ([§6.1][s6-1], [§11.2][s11-2]), stated below;
- the whole-tree obligation check;
- the closed leaf vocabulary ([§7.1][s7-1]), checked on every `init_x` because
  the walk in [§11.2][s11-2] rests on it.

Root [slots](#g-slot) fall out here too, as the root's input faces.

**The bound check** is the first type clause, and it applies at nominal faces:
the producer's declaration at `Float64` must be `<:` the entry at `Float64`.
Equality is the concrete degenerate case. Abstract-at-root is detected here.

**The walk-compatibility clause** is the second, and it applies to continuous
consumers only. It is decided by evaluating both declarations at a marker
scalar and comparing per leaf, and its diagnostic is
`WalkingFaceAtFrozenEntry`. It stays inside this stratum's charter because both
sides are declaration functions of `T`: declarations are evaluated, no user
stage code runs.

Stratum A also checks the declaration-completeness rules ([§11.2][s11-2]): a
store without its update, an event missing a [guard](#g-guard) or handler
method, a leaf mixing [tier](#g-tier) families, a contract signature whose form
contradicts the leaf's tier ([§11.5][s11-5]), and a primitive at the root.

`sample_times` validation is Stratum A's too, and it has two parts. The first
is per-entry validity against the constraints of [§8.5][s8-5]: wrapper-typed
values, `K ≥ 1`, `0 ≤ Φ < K`, `T > 0`, `0 ≤ τ < T`, and keys naming discrete or
scope children. Those violations are collected with path attribution. The
second is compilation into **`(anchor, m, c)` triples**. A triple carries a
discrete component's divisor and [phase](#g-phase) in the [tick](#g-tick) units
of its [anchor](#g-anchor) — the exact `(T, τ)` pair an `Absolute` entry
establishes.

The compilation is a fold down the tree, one rule per case:

| the fold meets | the triple it produces |
|---|---|
| the root scope | `(A₀, 1, 0)`, anchor 0 being the base grid itself: `(T, τ) = (Δt_base, 0)` |
| `Relative(K, φ)` under a scope at `(a, mₛ, cₛ)` | `(a, K·mₛ, cₛ + φ·mₛ)` |
| `Absolute(q, τ)` under any scope | **severs and re-seeds**: a fresh anchor `Aₖ = (period(q), τ)`, its subtree continuing at `(Aₖ, 1, 0)` |

Anchor 0 is symbolic until deployment. The `Relative` case is the affine law
([§8.5][s8-5]) in anchor-tick units. The canonical residue (`c < m`) holds
within each anchor's subtree by the same induction.

Everything except binding `Δt_base` — deployment's — happens in Stratum A.
Final divisors for anchored entries genuinely cannot exist until `Δt_base`
binds.

#### Stratum B — schedule

Stratum B is the single evaluation-feeds-structure step. It computes the
[schedule](#g-schedule):

- [Workspace](#g-workspace) — component-declared mutable scratch arriving as
  the `ws` bundle field — is allocated at the probing scalar. That is sound
  this early because the allocator reads only the instance and the scalar
  (row 77), so there is no layout dependence.
- Stage-1 [probes](#g-probe) run at `Float64`, on `init_x`/`init_m` values.
  They are well-founded, the no-[feedthrough](#g-feedthrough) stage taking no
  inputs.
- [Ports](#g-port) are classified over `output_types` alone: stage-1,
  auto-published, and the stage-2 remainder. A stage's `w` classifies nothing,
  being no part of the contract ([§11.3][s11-3]).
- The feedthrough graph is built from the wires carrying stage-2 ports, and a
  topological order over it follows. [§5.5][s5-5] cycle rejection applies.

The output is structural: names only, `T`-independent, branch-protected by the
branch-shape rule plus the always-on check ([§12.5][s12-5]).

#### Stratum C — activation, parametric in `T`

An [activation](#g-activation) is a re-run of Stratum C at a given scalar type.
The stratum holds everything type-shaped:

- The producers' output declarations are **evaluated** at the activation's `T`
  to type the [cells](#g-cell). That is the literal semantics
  ([§11.2][s11-2]): a continuous producer's two-argument declaration is called
  at `T`, and a discrete producer's plain one is read once and
  [pinned](#g-walked).
- The `init_x`-derived state type is [walked](#g-walked) by the leaf-walk rule
  ([§11.2][s11-2]).
- The probe chain runs in topological order, threading each stage's `w` to its
  one-hop consumers ([§5.2][s5-2], [§12.3][s12-3]), and observed is compared
  against declared.
- The flat `x` [buffer](#g-buffer) and the table are laid out.

The nominal `Float64` activation runs at build. Other activations re-run *only
this stratum* ([§12.4][s12-4]).

#### Deployment binding

Deployment binding sits after all three strata, at `Simulation` construction.
It binds `Δt_base`, `h`, `n`, `t_end`, the algorithm, `localization_tol`,
`localization_budget` and `firing_budget`, runs harmonic-grid validation, and
instantiates the tick schedule. Nothing in A–C depends on it.

`Δt_base` has exactly one of three sources, cross-validated:

- the explicit keyword, a `Rational`, `Period` or `Hz` value, from which `n` is
  derived as `Δt_base/h` and validated an integer ≥ 1;
- the `n·h` product when only `n` is given, today's rule, with the default
  `n = 1`;
- **derivation**, permitted only when every discrete component is anchored —
  that is, with anchor 0 unpopulated.

**Why.** Under that restriction `Δt_base` is pure bookkeeping that no
component's period depends on. An unanchored component's period is
`m·Δt_base`, and deriving with one present would let an anchor edit anywhere in
the tree silently rescale it — action at a distance.

**Rule.** If any unanchored component exists, deployment must declare
`Δt_base`. The refusal is constructive — the suggestion message
([§12.2][s12-2]).

Admissibility is exact GCD arithmetic over the **constraint pool**:

| quantity | value |
|---|---|
| the constraint pool | every anchor's period and every nonzero offset |
| an admissible `Δt_base` | one dividing every pool entry, equivalently one dividing `gcd(pool)` |
| the admissible set | `gcd(pool)/k`, for integer `k ≥ 1` |
| the derived `Δt_base` | `gcd(pool)` itself, the coarsest admissible value |
| per anchor `k` | `Dₖ = Tₖ/Δt_base`, `Φₖ = τₖ/Δt_base` |
| per component | `D = m·Dₖ`, `Φ = Φₖ + c·Dₖ`, `Δt = D·Δt_base` |

Resolution is therefore one division pair per anchor and one multiply-add per
component. `Dₖ` and `Φₖ` must both come out exact integers; otherwise a
`DeploymentInvalid`, naming the anchor with its declaring scope and key from
the provenance column. The per-component triples are the
[bound schedule](#g-bound-schedule), the printable artifact deployment binding
produces ([§12.2][s12-2]).

Deployment validation is collected like its declarative siblings
([§13.1][s13-1]). Collected and reported as `DeploymentInvalid`
([Appendix C][sC] — parameter, value, the violated constraint):

- a nonpositive `h`;
- an `n < 1`;
- a harmonic-grid violation;
- a non-dividing anchor period or offset;
- a declared `Δt_base` disagreeing with a declared `n`;
- an algorithm the [stepper seam](#g-seam) does not know;
- a nonpositive `localization_tol`;
- a `localization_budget` or a `firing_budget` that is not an integer ≥ 1.

The event parameters validate on their own terms only. They are
grid-independent, so they take no part in the harmonic-grid check
([§8.4][s8-4], [§8.6][s8-6]).

### 12.2 The `Build` artifact

`build(world) → Build` is a standalone entry point. `Simulation(world; …)`
(the spelling, [§15.4][s15-4]) is the convenience that calls it and adds
deployment binding, [buffers](#g-buffer) and the stopped-sim services.

The constructor is two entry points, not one: `Simulation(build::Build; …)`
accepts the artifact directly, and `Simulation(world; …)` is *defined as*
`Simulation(build(world); …)`. The build CI checked, the build an acceptance
test targeted, the build a [face](#g-face)-provenance table was printed from:
that artifact is the one deployed, never an assumed-equal reconstruction.

**Why.** Computed [boundary](#g-boundary)-connection bodies are ordinary user
code re-evaluated on every build, so equality between two builds of the same
world is an assumption the factorization removes.

Deployment binding still happens only at `Simulation` construction, whichever
entry point runs.

**The `Build` is immutable and may back any number of `Simulation`s,
concurrently** — true by construction once buffers are single-owner ([§12.4][s12-4]): each
`Simulation` materializes its own from the shared layouts, so nothing writable
is shared. The one mutable thing on the artifact is the lazily populated
[activation](#g-activation) cache, whose insertion [§12.4][s12-4] makes torn-state-free.
The `Build` is the
inspectable [contract](#g-contract) of the instantiation [§11.8][s11-8] gestures at — wire list, face
table, [schedule](#g-schedule), root [slots](#g-slot) as plain printable data. CI checks a model by
calling `build`; the acceptance tests target `build` errors directly;
`attach!` validates [device](#g-device) [bindings](#g-binding) against it. Build living only inside the
`Simulation` constructor was rejected (row 49).

**The schedule the `Build` carries is anchor-relative; the `Simulation` binds
it.** From [Stratum](#g-stratum) A the artifact gains two printable tables: the **[anchor](#g-anchor)
table** — each anchor's exact `(T, τ)` rationals with the declaring scope's
path and key — and the **[component](#g-component) table** of `(anchor, m, c)` triples with
their declaration provenance, the `Relative`/`Absolute` chain down the tree.
For a fully relative model the only anchor is `A₀` and the triples *are* the
final base-[tick](#g-tick) `(D, Φ)` pairs; when anchors exist, final divisors cannot
live here — they do not exist until `Δt_base` binds, and the same `Build`
already backs many `Simulation`s with different deployment parameters.
Binding ([§12.1][s12-1]) produces the **[bound schedule](#g-bound-schedule)**, a named printable artifact
on the `Simulation`: per discrete component, `(D, Φ, Δt)` with the anchor and
provenance columns carried through — the single source of truth for `Δt`
([§8.5][s8-5]), the substrate of the grid diagnostics below, and the table that
answers "when does what run, and what coincides with what". Its `show`-form
is the **hyperperiod chart**: the pattern repeats with period `lcm(Dᵢ)` base
ticks — the gate is pure modulo arithmetic, so one hyperperiod is the
complete truth, not a sample — rendered as a tick chart over
`k = 0 … lcm(Dᵢ) − 1`, with a guard for absurd hyperperiods.

**Example.** The model worked in [§8.5][s8-5] has three discrete components
under two scopes: `sample_times` declares `fcs = Relative(1)` and
`gnss = Absolute(Hz(50))` at the root, `inner = Relative(1)` and
`outer = Relative(5, 2)` under `fcs`. Deploy it at `Δt_base = 2 ms`. The
`Absolute` entry seeds the anchor `A₁ = (1//50, 0)` and the rest of the tree
stays on anchor 0, so the three components carry these values through binding:

| | `inner` | `outer` | `gnss` |
|---|---|---|---|
| declaration | `Relative(1)` under `fcs` | `Relative(5, 2)` under `fcs` | `Absolute(Hz(50))` at the root |
| anchor | `A₀` | `A₀` | `A₁` |
| triple `(m, c)` | `(1, 0)` | `(5, 2)` | `(1, 0)` |
| bound `(D, Φ)` | `(1, 0)` | `(5, 2)` | `(10, 0)` |
| bound `Δt` | 2 ms | 10 ms | 20 ms |

`gnss` is the entry whose divisor could not exist before `Δt_base` bound:
`D = m·D₁`, with `D₁ = T₁/Δt_base = (1//50)/(1//500) = 10`.

**Grid diagnostics print from the pool, exactly.** The refusal path's
suggestion message and the derivation path's info line share one substrate:
the coarsest admissible `Δt_base` with the admissible set `gcd(pool)/k`, and
per-entry attribution. Attribution has two forms.

**Leave-one-out refinement factors** are the first form:
`r_p = gcd(pool ∖ p)/gcd(pool)`, an integer ≥ 1 read as "how much coarser the
grid would be without this entry". Every `r_p > 1` is listed rather than one
culprit crowned, joint responsibility being the honest answer.

**Prime attribution** is the second form: each prime power of `1/Δt_base` is
traced to the pool entries whose denominators supply it. That pinpoints what
an edit changed.

When an offset is a driver, the message adds the nearest non-refining
alternatives — the admissible offsets on the grid the rest of the pool
supports. That turns the diagnostic into a repair.

Blame is computed against the actual pool. A simple-fraction-of-its-period
test stays authoring guidance and never becomes the engine's (row 187).

The derivation path — the one place refinement happens silently — always prints
the derived value with its drivers, and carries the one advisory:
`GridUtilization` ([Appendix C][sC]), reporting `min_i Dᵢ` — base ticks between
the fastest component's ticks — as "grid is N× finer than the fastest declared
work" with the drivers named, information rather than scolding, since a scope
deliberately declared finer than its fastest member to buy stagger room
([§8.5][s8-5]) legitimately inflates the metric.

### 12.3 Probing and input synthesis

**[Probe](#g-probe)-everything scope.** The nominal [activation](#g-activation) probes every user function —
stages, `f`, `g`, [guards](#g-guard), handlers, `project` — once, at the initial state,
with real values, checking shape/type conformance and discarding results (all
are pure; the cost is one evaluation each). "Fails loudly at build time where
possible" ([§11.1][s11-1]) decides this: a malformed `f` return must not wait for the
first integrator step. Probes see only the initial state's branch — the
marginal coverage is earliness, not completeness; the always-on check ([§12.5][s12-5])
remains the completeness backstop.

**Probe argument sourcing.** `x`/`m` come from `init_*` declarations
(declared by value); `y_x` from the stage-1 probes' *returns* (an
auto-published name is a framework write, never a probe product, so it is
absent from the hand-down — [§5.2][s5-2]); wired inputs from
upstream products, real values available because the stage-2 chain is probed in
topological order. **`w` threads through the probe in stage order.** It follows
the one-hop law of [§5.2][s5-2]. If the [component](#g-component) defines an
`h_xu`, the `h_x` probe's second return [slot](#g-slot) enters that stage's
probe [bundle](#g-bundle) (the NamedTuple of zero-copy views a component
function receives). Otherwise it enters the downstream bundles. The `w` that
survives the last output stage enters the `f`, guard and handler probes. So
every consumer is probed against the same value it will receive at run time.
Three checks ride the same pass. The first is the return's shape: a stage
returning something other than a `NamedTuple` or a `NamedTuple` pair fails
here, `nothing` in either slot included. The second checks the consumer's `w`
reads against the observed field set — the [did-you-mean](#g-did-you-mean) of
[§11.3][s11-3], delivered through the [§13.2][s13-2] framing diagnostic. The
third is the dead-stage rule: a stage returning bare `(;)`, neither
[ports](#g-port) nor `w`, is `DeadStage`, fail-fast.
The [bundle law](#g-bundle)'s two remaining fields ([§5.2][s5-2]): `t` is
probe-scoped `0.0` — deployment binds no clock and `t₀` post-dates even
deployment ([§14.5][s14-5]), so like `Δt` below it is a fabricated, probe-scoped
value; `ws` comes from invoking the component's `workspace` allocator at the
probing scalar, which reads only the instance and the scalar (row 77),
deriving nothing from layouts, so it runs before the [Stratum](#g-stratum) B probes that
need it. Exactly one kind of terminal has no producer: **root
slots**. The build synthesizes their values via `probe_value(::Type)`:
framework methods for `Real` (`zero(T)`), `Bool` (`false`), enums (first
instance), ultimate fallback the zero-argument constructor `T()` — which is
where well-behaved constrained types already put their valid default (`RQuat()`
= identity; the `@kwdef` convention supplies it broadly). `probe_value` is
**overridable**: a type whose valid default is not reachable that way declares
its own method, which is also the [seam](#g-seam) a [walked](#g-walked) type uses to state a
constrained default. No method → build error, in the didactic register: it
names the [face](#g-face) and the type, and asks for one of the two fixes ("no
`probe_value` for `Ranged{Float64, -1, 1}` at face `pilot.elevator_axis` —
define `probe_value(::Type{Ranged{Float64, -1, 1}})` or a zero-argument
constructor"). Synthesis never meets an abstract type:
root slots are concrete by the tight-bound rule ([§11.2][s11-2]; the slot type
is the consuming entry evaluated at `Float64`), and [abstract
entries](#g-abstract-entry) only
occur on component-fed inputs, which the probe sources from upstream products.
Physically silly values are acceptable by
construction: the probe checks types, and return types that depend on input
*values* are type instabilities (banned by the branch-shape rule); the [§4.3][s4-3]
write-side granularity corollary keeps root slots predominantly scalar, so the
surface is small. Rejected (row 51): inputs declared by value à la `init_x`, NaN
poison values, and init-service values.

**Probe values are strictly probe-scoped.** Everything the probe writes is
garbage once the build finishes; probe values never double as initial slot
values — that would smuggle in the default semantics rejected above. The
same doctrine covers the clock: `Δt` in seconds does not exist until
`Simulation` binds `Δt_base` (deployment post-dates the build), so
discrete-[tier](#g-tier) probes supply a placeholder period (`1.0`) in the bundle —
a fabricated, probe-scoped value like any synthesized input; the probe
checks types, not physics.
`Simulation` must not reach its first [boundary](#g-boundary) with uninitialized slots;
enforcement is the pre-write `UninitializedSlots` check carried by every
complete-world application — `init!`, trim setup, trim commit
([§14.6][s14-6]).

**The author's side of that bargain.** Silly values are acceptable *because* the
author is obliged to accept them: **stage code must be total over type-valid
inputs** — every probed user function (stages, `f`, `g`, guards, handlers,
`project`) evaluates without throwing on any input satisfying its declared
types. The domain is type-validity, not the probe's particular synthesized
values: the branch-shape rule already bans value-dependent return types, so
types are the only domain the framework can speak of, and the probe is the
enforcement moment, not the reason (row 142). Two consequence sites, the same throw at both: at
build it is a `UserCodeFraming`-wrapped build failure ([§13.1][s13-1]) whose diagnostic
points at code that is "correct" on every trajectory it has ever seen; at
runtime it is a `StepError` and the run ends `errored` ([§13.4][s13-4]) — exceptions from
model code are always abnormal ([§13.5][s13-5]). Three habits of shipped code have
sanctioned spellings. A *plausibility* check meaning "stop the run" — a strut
throwing on a touchdown overload — is a published `Bool` output face plus
`stop_on` ([§13.5][s13-5]), machinery already there. A *self-consistency* assert — an
author checking that their own contact algebra cancels a velocity component to a
hard tolerance — is a regression test about that algebra, and its home is the
test suite; it is also the most probe-fragile of the three, since a
near-degenerate synthesized geometry can keep the cancellation algebraically
exact while missing an absolute tolerance in floating point. And there is the
*defensive exhaustiveness* branch: an `else error("unrecognized surface type")`
over a closed enum, or a coefficient constructor asserting an ordering of its
arguments when that constructor runs per step inside a stage. Such a branch is
not banned validation but **mislocated** validation. Totality over a closed
enum means handling every instance, and an `else error` is an admission that
the function is partial. Parameter validation belongs where user-controlled
data enters: the constructors of parameter and instance values, which run
before the build, where asserts are perfectly legitimate. It never belongs
inside a stage, on probe-fed data.

### 12.4 Activations: executable sets, laziness, caching

An **[activation](#g-activation) at `T`** re-runs [Stratum](#g-stratum) C with a
different scalar:

- producer-fed [cells](#g-cell) re-typed by *evaluating* the producing
  [component](#g-component)'s output declaration at `T` ([§11.2][s11-2]) — a
  continuous producer's two-argument declaration called at the scalar, a
  discrete producer's plain declaration pinning;
- root-[slot](#g-slot) cells re-typed by *evaluating* the consuming
  `input_types` entry at `T`, which [§11.2][s11-2] reads permissively — a `T`
  entry follows the activation, a `Float64` entry stays frozen;
- the state type by the walk over `init_x`'s, with table and state
  [buffers](#g-buffer) re-laid-out;
- [workspace](#g-workspace) allocators re-invoked at `T`, not introduced — the
  first invocation precedes the Stratum B [probes](#g-probe)
  ([§12.1][s12-1]/[§12.3][s12-3]), and a
  [continuous component](#g-continuous-component)'s scratch carries the
  activation's scalar ([§7.3][s7-3]);
- the probe chain re-run.

Structure and [schedule](#g-schedule) are `T`-independent by construction.

**Each activation probes exactly the function set it can execute.** A `Dual`
activation (linearization, gradient trim) evaluates the model at a frozen
instant: discrete stages are gated off holding `Float64` values (the [§11.2][s11-2]
frozen-constant semantics), [guards](#g-guard) and handlers never run (event localization
is `Float64` [sweeps](#g-sweep) by design, [§8.4][s8-4]). Only the continuous output stages
(`h_x`/`h_xu`) and `f`
ever see a `Dual` — so only they are probed. Probing the discrete stages, `g`, or guards at `Dual`
would check code against a number type it cannot receive. One rule, no
special cases; the [§5.6][s5-6] tracer activation follows it identically. "Tracer
activation" names the *global* set-tracer (row 12) — a whole-model run at the
tracer scalar, an activation like any other. The cycle classifier ([§5.6][s5-6])
is the other variant (row 12), the schedule-free per-member local
[trace](#g-trace), which runs in
Stratum B's failure path and is not an activation at all.

**Lazy, with an opt-in exhaustive mode.** Non-nominal activations run at first
request, not at build: the dominant cost is compiling the continuous chain a
second time at `Dual`, pure waste for interactive fly-around use. The price,
stated openly: `build` succeeding does **not** certify the model linearizable —
a [pinned](#g-walked) `Float64` ([§7.2][s7-2]), whether hidden in a constructor or written into an
output declaration at a leaf that really participates (the per-leaf
forgotten-`T`, [§11.2][s11-2]), lurks until the first `Dual` activation
detonates it
at the probe, naming the offending constructor or leaf. The repository's test suite
pins the invariant instead, as policy rather than advice (row 166):
**every component gets a `Dual` activation built in CI** —
`build(world; activations = (Float64, ProbeDual))`
(or a `check` entry) runs the exhaustive set, catching both genericity
violations and forgotten-`T` leaves at PR time, at the cost of a Stratum-C
re-run per component. The same keyword is the recommended idiom for the
parallel-sweep register ([§9.1][s9-1]): pre-materialize the activations the
sweep will need and
the shared `Build` is a fully immutable artifact, with no synchronization on any
path. [`ProbeDual`](#g-probedual) is the framework's exported canonical probe
scalar — `const ProbeDual = ForwardDiff.Dual{ProbeTag, Float64, 1}` — because
an activation is keyed by a *concrete* scalar type and the bare `Dual`
`UnionAll` cannot key one, be [walked](#g-walked) to, or answer `zero(T)`. Its width is
arbitrary: what CI pins is genericity, not any particular Jacobian, so one
canonical width suffices even though [§14.10][s14-10] chunks at whatever widths it needs.

**Caching is implementation detail, not semantics.** An activation is a pure
function of the build and the concrete scalar type — so the cache is the
`Build`'s, and it holds (layouts, compiled plans, validated-flag) keyed by that
type: immutable once constructed, hence freely shareable. **Buffers are never
cached**, because every buffer set has exactly one owner. The `Simulation` owns
its nominal activation's buffers — materialized from the cached layouts at
construction, what the loop's zero-allocation stepping runs on — and every
service invocation owns the scratch set it instantiates from those same layouts;
[§14.8][s14-8] states this for `trim!`, and it is the general rule, not a trim-local one.
Compiled code is cached by Julia itself, process-wide; what the framework cache
saves is the expensive part — probe re-runs, layout construction, and Julia's
compilation of the `Dual` chain — which is what actually amortizes in
activation-reusing loops (the envelope-grid gain-schedule case: hundreds of
trim-then-linearize points paying those costs once). What does not amortize is
the per-point allocation of a working store set, O(model size) and trivial
against the solve it feeds: the zero-allocation invariant ([§7.5][s7-5]) is
scoped to the stepping loop, and the services were always allocation-tolerant.
Nothing
numerical is ever cached.
Note `Dual{Tag,V,N}` carries the partial count: a different seeding width is a
different scalar type, hence a separate entry and a separate Julia compile.
**Lazy materialization is torn-state-free**, normatively: concurrent first
requests for the same activation must never expose partially populated cache
state. The mechanism is unspecified — a guard around insertion suffices, paid at
service time and never on the hot path — and since an activation is a pure
function of build and scalar, the worst benign race is duplicated work; torn
state is excluded by [contract](#g-contract), not by luck.

### 12.5 The always-on conformance check

The [probe](#g-probe) validates each function *once*, on the initial state's branch; the
schema-authority bargain's second clause ("at first execution otherwise",
[§11.1][s11-1]) is discharged by leaving the probe's comparison permanently in place.
At the point where the [executor](#g-executor) (the compiled execution form of
the schedule) stores a stage return into the table, it holds the complete
expected return type at this [activation](#g-activation). That type comprises
the declared types of the names *this stage* produces, as fixed by
[Stratum](#g-stratum) B's stage classification: one concrete `NamedTuple` type
per ([component](#g-component), stage). It is drawn from `output_types(c, T)`
on a continuous producer and from `output_types(c)` on a discrete one
([§11.2][s11-2]). Auto-published names belong to no stage's expected type;
the framework writes those [cells](#g-cell) itself. The executor canonicalizes the
observed return to that type's field order by a type-level reorder
(`NamedTuple{names(Expected)}(y2)`) and performs a single
type test against that type (conceptually `y2 isa Expected`). Type-stable conformant
code: the compiler proves the return type, decides the test at compile time,
and deletes it — zero instructions. Branch-divergent code: the test survives
as a runtime check, nanoseconds on conformant branches, a loud located error
on the divergent one at its first execution. (Type-unstable-but-conformant
code pays nanoseconds on top of the dynamic dispatch it already bought.)

**The names are the pairing; field order carries no semantics.** `Expected`'s
order is an internal fact — derived from `output_types`, stage-filtered,
auto-published names removed, an order no single declaration shows the author.
The author never reproduces it: a return spelling the right names at the right
types conforms in any order. `(; P = M*ω, M_shaft = M)` and
`(; M_shaft = M, P = M*ω)` are the same return. This is the general rule at
every author↔framework `NamedTuple` [seam](#g-seam) ([§14.7][s14-7] states it
for the trim problem's decisions and residuals). It is also what downstream
consumption already assumes: the scatter writes each returned field into its
own *named* cell ([§4.3][s4-3]). Order-sensitivity in the check would therefore
be incidental strictness rather than protection. The canonicalizing reorder
costs nothing: it is a compile-time permutation of an already-typed value,
register shuffling SROA deletes. It folds exactly where the test folds, so the
economics (row 53) hold — one whole-type test, never per-field checks. The
canary ([§7.5][s7-5]) verifies the fold empirically rather than by assertion.
What is an error is a key-set mismatch or a per-field type mismatch, reported
by the [payload](#g-payload) below. A permutation is not an error at all —
which is equally why that diff never has to express one.

**Exact match at nominal; embed-accept at declared-`T` leaves.** At
the nominal activation — the only one that ever runs in real time — the check
is an exact type match, no convert-on-write, one baked `isa` that
folds away (row 53). The error can afford to be
didactic: "field `M_shaft`: expected `Float64`, got `Int64` — return
`zero(x.ω)`, not `0`". Under a non-nominal activation (a re-run of Stratum C at
a given scalar type) the two leaf kinds the declaration ([§11.2][s11-2]) distinguishes
are checked differently. A **declared-`T` leaf** — the author wrote `T` there —
accepts exactly two types: the activation scalar or `Float64`. The activation
scalar is the fast path, the baked `isa`. A `Float64` the executor **embeds**
as a zero-partial constant (`convert` through the leaf). Struct-valued
[ports](#g-port) use the standard cross-eltype constructor, a missing one
failing loudly with both types named. Nothing else is accepted. A
**declared-[pinned](#g-walked) leaf** — the author wrote a concrete type,
`Float64` at the head of the list — takes the nominal-style exact check at
*every* activation, its declaration having said the leaf never carries
partials. An observed `Dual` there is the per-leaf forgotten-`T` error,
that being the one honest cause. The didactic hint is attached: "if `F`
participates in differentiation, declare it `T`". The embedding is exact, not
lenient. Promotion is airtight and there is no lossy `Dual → Float64` cast, so
a `Float64` observed at a declared-`T` leaf means no `Dual` entered its
computation. Its true derivative along every seeded direction is zero, which
is precisely what the embedded constant says. This scopes the blanket
convert-on-write rejection to the nominal check (row 53). The bug that
rejection guards against — silently zeroed partials — cannot arise from honest
code, because accidental `Float64`s from `Dual` operands are impossible
(`MethodError` at the operation site). The residual is **deliberate stripping**
(`ForwardDiff.value`): a stated intent to discard partials, producing a silent
zero in the Jacobian. That is the stop-gradient idiom, occasionally legitimate
— deliberately frozen couplings, opaque non-Julia wrappers. Applied
mid-expression it is equally invisible to a strict exact-match rule, so the
leniency costs nothing. What it need not be is invisible to the schema:
**the declared-pinned leaf is the schema-visible freeze** — an
author who means to strip declares the leaf `Float64` and strips inside the
stage, and the check above holds the freeze to its word at every activation.
Stripping mid-expression at a leaf still declared `T` remains legal and remains
unseen, as the sharp tool it is.

**`w` is checked at the nominal activation, and nowhere else.** A stage that
returns the `(y, w)` pair ([§5.2][s5-2]) gets a second baked `isa` beside the first,
against the type the **nominal probe observed** — no declaration exists to draw
an expectation from, and none is wanted, `w` being a value in flight rather
than a cell to type. It folds exactly as its sibling does when the stage is
stable. When it does not fold it earns its nanoseconds: the branch that quietly
returns a `w` of a different shape is otherwise invisible until it shows up as
an allocation in someone's benchmark. Here it is a located error naming the
offending field at that branch's first execution. The message cites its
authority honestly — expected *what the nominal probe observed*, not what
anybody declared. Under **non-nominal activations the `w` test is
absent**: with no declaration there is no branch-independent anchor, there is
no cell whose typing the check would protect, and a probe-observed expectation
would reject the constant-branch idiom that [§11.2][s11-2] keeps legal on the
output side. Nothing is lost in correctness — a `Float64` arriving in `w` under
a `Dual` activation is a true zero-partial constant by the embedding guarantee
above, and promotion at its first use with a `Dual` operand is exact. The
reasoning is [§11.3][s11-3]'s, recorded in row 165.

**Uniform across all probed functions.** `f` checks against `X`'s own shape at
the activation's `T` ([§7.1][s7-1]: a scalar leaf expects a `T`, an `SArray`
leaf the same `SArray` at `T`). Its predicate is "every field scatters into its
field's block at `T`", which is what makes derivative completeness structural
rather than a matter of author discipline. [Guards](#g-guard) check against
their probe-derived [predicate](#g-predicate) form (below); `g` against its
leaf's `x` shape; handlers against the [§5.2][s5-2] return law, key by key.
`project` checks against `X`'s own shape at `T`, **complete**, since its result
is written back to the [buffer](#g-buffer) wholesale at both of the
[schedule](#g-schedule) positions ([§5.3][s5-3]) and a
[projection](#g-projection) with a mode-dependent branch first executes its
second branch at run time. That is the same predicate as a handler's `x` key.

**Handler returns, key by key.** The returned NamedTuple's key set is checked
first: an unknown key, or a key naming a store the component does not declare,
is a build error with [did-you-mean](#g-did-you-mean) against `{x, m}` narrowed to the stores that
exist — the [bundle law](#g-bundle)'s classification running in the return direction. Then,
per present key: `x` must be **complete** against the state field set (and
conformant at `T` like any state value), while `m` may be **partial**, checked
against a names-subset-with-matching-types predicate — still a type-level
computation that folds when inferred. An absent key is not an error and not a
no-op to diagnose: it is the handler saying it does not touch that store. The
completeness asymmetry is storage-shaped: `x` must be complete because it lives
in a flat buffer written back wholesale, while `m` may be partial because
`m` lives in per-field stores where a partial merge is the natural write.

**Guards have two admissible forms** ([§2.1][s2-1]), so their check is form-aware
rather than a flat `isa Bool`: a `Bool`-form guard's probed return is `Bool`, a
sign-form guard's is the nominal scalar. Guards run only at the nominal
activation (row 52), so no parametrized-leaf case arises here. Any other probed
return type is a build error naming both admissible forms. There is nothing
further to check: the probed form *is* the detection policy ([§8.4][s8-4], row 179),
so no form/policy mismatch can be declared.

**Failure payload:** component path, function, field-level diff (missing /
unexpected / per-field expected-vs-observed), simulation time. Deliberately
absent: the source branch (values carry no provenance; the diff identifies
it). The always-on input [trace](#g-trace) makes every such failure **reproducible by
[replay](#g-replay)** — the error names the [boundary](#g-boundary) to replay to
(`to_boundary`, [§10.7][s10-7]). At run time the failure
travels as a species of `StepError` through the single catch site ([§13.4][s13-4]),
which adds the loop-level nonfinite-state check as its divergence sibling.

### 12.6 Stopped-sim services as Stratum-C clients

Sketched here because it grounds the strata; the services themselves are [§14][s14].
The C172 trim problem (`c172.jl`: `TrimState`, `TrimParameters`,
`θ_constraint`, the `ẋ`-reading cost) transfers near-verbatim:

- **Trim** is a write-condition → [sweep](#g-sweep) → read loop on an [activation](#g-activation) — by
  default the `Dual` activation, decision variables seeded for exact residual
  Jacobians ([§14.7][s14-7]); the derivative-free fallback runs the same loop on the
  nominal `Float64` activation (a re-run of Stratum C at a given scalar type)
  with no new activation needed, and the always-on checks ride along either
  way. Decision variables stay opaque to the framework (only the assignment's
  *output* is framework vocabulary). `assign!` inverts from in-place mutation +
  self-invoked `f_ode!` to a pure function returning a [condition](#g-condition) value (state
  by path, modes, [slots](#g-slot) by [face](#g-face)) that the service writes and evaluates. Domain
  math — the pitch constraint, `Kinematics.Initializer`, per-residual scalings
  and the equilibrium-subset choice — survives aircraft-side, with one
  respelling: the initializer's `atmosphere::Model` argument becomes a
  [field handle](#g-field-handle) ([§4.4][s4-4]), built at value level by the atmosphere's
  [value-level constructor](#g-value-level-constructor) or held directly as a rig slot value ([§14.1][s14-1], [§14.9][s14-9]).
- **Linearization** is a `Dual` activation plus seeded sweeps: gather/scatter
  over the canonical layout replaces the hand-written
  `get_x_ss`/`assign_x_ss!` layer (the deletion discharged, [§7.1][s7-1]); root slots are
  the input surface; frozen discrete outputs are constants with zero partials,
  which is exactly "linearize with the discrete state held" ([§11.2][s11-2]). Gradient-based trim —
  decision variables seeded through the `T`-generic assignment math — is the
  default ([§14.7][s14-7]).
- The generic service loop (vectorization, optimizer setup, bounds packing,
  solved-condition write-back including root slots and the [trace header](#g-trace-header)'s
  slot capture) replaces today's per-aircraft NLopt plumbing. A failed trim
  leaves the simulation's stores untouched — an improvement over today's
  warn-but-assign `f_init!`.

### 12.7 The compiled executor

The [schedule](#g-schedule) exists in two representations at two lifecycle stages. In the
`Build` it is plain printable data ([§12.2][s12-2]) — paths, stage names, order — the
authoring and diagnostic form. At `Simulation` construction, and per
[activation](#g-activation) (a re-run of Stratum C at a given scalar type), that data is
compiled into the execution form: **a concretely-typed tuple of entries over
statically typed [cell](#g-cell) storage, traversed by a compile-time-unrolled walk**.
This is a forced move, not a preference: the zero-allocation invariant ([§7.5][s7-5]),
the fold-away conformance test ([§12.5][s12-5]) and the zero runtime graph logic ([§5.1][s5-1])
are reachable only under full specialization (row 86). An entry carries what
selects code — [component](#g-component) type, stage — in type parameters, and what is plain
data — [tick](#g-tick) divisor and [phase](#g-phase), the [bundle](#g-bundle)'s `Δt`, layout offsets — in fields;
gating compiles to `(idx − Φ) % D == 0` inside the specialized *[boundary](#g-boundary)* body,
the interior bodies holding no discrete entries to test ([§8.5][s8-5]).

**Cells are stored per element type, not per cell.** The [signal table](#g-signal-table) is one
contiguous [buffer](#g-buffer) per element type — the construction pointed at
signals rather than state ([§7.1][s7-1]) — and a cell address is a build-time
offset into it, carried
in an entry *field* with the [port](#g-port) type as the address's own parameter; gathers
reconstruct and scatters flatten through the same leaf walk, so the closed
vocabulary earns its keep twice. This is the entry rule above paying rent: two
instances of one component type then differ only in field values, share an
entry type, and compile to **one** body — where a store enumerating every cell
in its own type, addressed by index in the type domain, compiles one body per
instance and grows the store type with the model. The choice was measured
rather than argued (row 162, `prototypes/cellstore_bench`).

**[Phase bodies](#g-measurement-seam) are the outer decomposition, and they are
semantically forced.** The [boundary sweep](#g-sweep)'s `h_x` block is order-free
by definition (the no-[feedthrough](#g-feedthrough) stage reads no `u`). The
`h_xu` block gates in the [due](#g-due) discrete stages — those whose components
this boundary admits by their compiled `(D, Φ)` pair. It is the only
topologically ordered one. The `f` block — the [RHS](#g-flow) body the stepper
calls per stage evaluation — and the `g` block are order-free with disjoint
writes.
[Guards](#g-guard) and handlers are their own small callables inside the
[§8.6][s8-6] iteration.

**Each sweep block compiles in two arities off one entry list**, along the
interior/boundary split that [§8.5][s8-5] fixes. The zero-arg
`sweep_hx()`/`sweep_hxu()` are the interior variants, over continuous entries
only. That is what makes `@ballocated(sweep_hxu()) == 0` a well-defined
measurement *of the interior path*, rather than of whichever tick phase the
simulation happens to be sitting in. The `sweep_hx(tick)`/`sweep_hxu(tick)`
forms are the boundary variants, gating their discrete entries by
`(idx − Φ) % D` against the passed index, symmetric with `ticks(tick)`. `rhs`
takes no index (row 147). One gate serves all three tick-sensitive blocks —
due-ness is per component, per boundary, never per stage — and `t*`'s empty
due set is **arity selection, not an index trick** (rows 147, 185), so the
`t*` iteration runs the zero-arg arities, whose compiled bodies contain no
discrete entries ([§8.5][s8-5]).

Across passes these bodies communicate only through the stores and the table,
so the [seams](#g-seam) between passes cost nothing — no values cross them.
**Within** a pass one kind of value does: a stage's `w` ([§5.2][s5-2]) is
handed to its one-hop consumers as an ordinary SSA argument, across block and
chunk seams alike, never through storage. That is what makes the fusion a
design constraint rather than an optimization: a step's sweep and its `f`
block compile into one pass, and an event round's sweep, its guards and its
fired handlers into another ([§8.6][s8-6]), because that is the scope over
which a private intermediate is fresh by construction. Two doors this
structure opens for free, recorded not committed: deterministic parallel
evaluation of the order-free blocks (disjoint writes, and no floating-point
reductions to reorder — [§6.2][s6-2] made every sum an ordered junction
entry), and finer recompilation granularity (editing a discrete component
invalidates the boundary body, not the RHS body — literal under the two-arity
split, discrete entries existing only in the boundary variants).

**[Chunking](#g-chunking) bounds the compile cost.** Within a large block the tuple splits
into chunks behind non-inlined but statically-typed function barriers.
Inside a chunk everything the design relies on survives — static dispatch,
inlining, view SROA, check folding, zero allocation; at the seams only
cross-entry fusion is lost, which a table-mediated dataflow barely had.
Chunk size is the implementation's *only* representation freedom (fully
fused and chunk-of-one are its endpoints), and it converts the compile cost
from superlinear in the largest body to linear in entry count.

Measured anchors, taken 2026-07 over synthetic ~15-op bodies on Apple Silicon,
with the last two rows extrapolated to a C172X-scale model — roughly 200–400
entries, larger bodies ([§15.4][s15-4]):

| case | activation | compile time |
|---|---|---|
| 400-entry sweep, fused | `Float64` | ~0.8 s |
| 400-entry sweep, chunked | `Float64` | ~0.34 s |
| 400-entry sweep, chunked | 8-partial `Dual` | ~9 s |
| C172X-scale model, extrapolated | nominal | seconds |
| C172X-scale model, extrapolated | `Dual` | tens of seconds, before mitigation |

The fused curve is visibly superlinear. An 8-partial `Dual` activation
multiplies instruction count ~20×, and its chunked curve is linear —
instruction-bound rather than structure-bound. Re-measurement on the real
vehicle skeleton is a [§16][s16] migration item.

**The mitigation ladder**, in order: activations are lazy ([§12.4][s12-4] — a session
that never linearizes never compiles `Dual`); non-nominal activations may
compile at reduced optimizer level (their sweeps run inside service loops
where microseconds are irrelevant — a one-line per-module policy); and
activations bake into package images via ordinary precompile workloads — an
aircraft package exercising build-plus-one-sweep per activation turns TTFX
from a session tax into a CI artifact.

**[Views](#g-view) are spelled rebuild-per-call.** Every entry constructs its bundle (the
NamedTuple of zero-copy views a component function receives) at its own
position; there is no framework-maintained hoisting and therefore no
cache-invalidation obligation. Hoisting belongs to the code generator: CSE
merges repeated loads exactly where no intervening store invalidates them —
which is precisely the staleness rule — and the sweep-varying bundle fields
(`u`, `y_x`) are per-call by topological necessity either way ([§7.1][s7-1]).

**Construction is type-opaque; only the [executor](#g-executor) specializes.** Schedule
tuples are built from untyped buffers and splatted once. Generic tuple
utilities — range indexing, long `ntuple` closures, naive recursion — are
inference traps at schedule length (a 400-entry heterogeneous tuple can send
generic `getindex` inference into combinatorial collapse), so the compiled
tuple's type has exactly one consumer: the unrolled walk.

**The phase bodies are the [§7.5][s7-5] [measurement seam](#g-measurement-seam).**
`phase_bodies(sim)` returns the compiled bodies of the nominal activation as
named callables bound over the simulation's own buffers. The four blocks:

- `rhs` — the `f` block.
- `sweep_hx` — in both arities.
- `sweep_hxu` — in both arities.
- `ticks` — takes the tick index its entries gate on.

Returned with them are the per-event guards and handlers and the per-component
`project` callables, keyed by the model's own [roster](#g-roster).

The four-body roster is fixed and total: the accessor returns all of it always,
whatever the model happens to declare. A model with no discrete components, no
events or no continuous state at all still gets every body. The empty ones are
legal, compile to no-ops, and their `@ballocated` assertion passes vacuously.
That is the point, because consumers then iterate the roster uniformly, with no
existence checks and no per-model branching in the measurement code.

One promise, in the diagnostic register ([§13.5][s13-5]): **these are the bodies
the loop runs** — not re-derivations. That is what makes the measurement
honest, and why each callable carries the real in-loop argument types by
construction. Those types are the thing a hand-built standalone test cannot
reproduce, and row 116 records why per-component tests cannot discharge the
invariant.

CI is warm-then-assert over the roster: one call compiles, then
`@ballocated(body()) == 0`. It asserts at per-body granularity, each sweep
arity in its own right — the interior call bare, the boundary call at a due
index. So a documented [§7.5][s7-5] tolerance loosens exactly one assertion.
This is the successor of the migration suite's
`@ballocated f_ode!`/`f_step!`/`f_periodic!` idiom and the seam the
[§16][s16] FlightCore comparison measures through.

Publication is not a phase body — the carve-out ([§7.5][s7-5]) made structural:
what the accessor exposes is exactly what the invariant claims is zero.
Invoking bodies in isolation mutates the simulation's buffers outside any
[frame](#g-frame) sequence (a tick entry advances discrete state with no clock
advance), leaving them valid but off-trajectory. A session that wants to
continue meaningfully re-runs `init!`.

---

# Part IV — Failure and services

## 13. Error discipline

[§11.4][s11-4] fixed what must be caught and where; [§12][s12] fixed when each fact is checked.
This section fixes how failures are *reported* — the reporting policy, the
diagnostic representation, the runtime failure story, and the seam between "the
model reached a terminal state" and "the run should end". Two of FlightCore's
paid-for lessons ground it: the compact-backtrace discipline (parameterized
model types make rendered output unreadable) and the `SimulationTermination`
machinery, which [§13.5][s13-5] replaces.

### 13.1 Reporting policy: collect the checks, fail the evaluations fast

The build's failure sites split into two populations, and that split settles the
fail-fast vs. compiler-style question: each population takes the reporting
policy that fits it.

- **Declarative checks over collected structure** — unconnected inputs,
  two-producers, wire typos and type mismatches, [face](#g-face)-name uniqueness,
  `output_types`/state-field consistency, `sample_times` validation. Each is a
  pass over a list; the whole-tree obligation check literally computes *the
  set of* inputs whose obligation chain never terminates. Reporting every
  violation is the natural output of the pass; truncating to the first would be
  extra work. These failures also cluster in practice: a freshly written
  [assembly](#g-assembly) has five unwired inputs; a renamed [port](#g-port)
  breaks three wires. **These passes collect:** each returns its full violation
  list.
- **User-code evaluation** — the [boundary](#g-boundary)-connection bodies, run
  in [Stratum](#g-stratum) A (one of the build's three phases: structure,
  schedule, activation); the stage-1 [probes](#g-probe) in B; the probe chain in
  C. When user code throws there is no meaningful rest-of-collection. A failed
  `input_connections` leaves the parent's face derivation undefined, and a
  failed stage-2 probe starves every downstream probe of its wired inputs, since
  [probe values](#g-probe-value) flow topologically ([§12.3][s12-3]).
  **The first user-code exception aborts the phase** (row 57).

Strata are barriers. A stratum that produced any error-severity diagnostic, of
either kind, throws before the next stratum begins; probing against unresolved
wiring is meaningless. The only partial results ever carried past a failure are
violation lists from pure checking passes. So none of the three strata
([§12][s12]) needs the machinery for carrying partial internal results across a
failure that it would otherwise need. That machinery is the cost that kept this
decision open, and it never materializes.

**No cascade suppression within a stratum** — a deliberate simplification
(row 57). A wire typo'd as `:throtle` produces both a
[did-you-mean](#g-did-you-mean) error (the offending name plus the list-in-hand
it should have matched) and an unconnected-input error for the intended
`throttle`; both are reported. They render adjacently (diagnostics sort by
path), and the pairing is self-explanatory.

### 13.2 Diagnostics: structured values, one carrier exception

Both reporting policies move the same thing. What a collecting pass returns and
what a fail-fast site produces is in either case a *diagnostic*, so the shape of
that value decides what an acceptance test can assert and what a user reads.

**Rule.** A diagnostic is a plain value from a small closed set of
[kinds](#g-kind). **[Appendix C][sC]** enumerates that set normatively: kind
name, [payload](#g-payload) fields, owning section, severity. That table is the
artifact the [§11.4][s11-4] acceptance tests and the error-message work are
written against. Each kind carries its own structured payload: endpoint paths,
[face](#g-face) names, expected/observed types, a severity, and the
*list-in-hand* a [did-you-mean](#g-did-you-mean) needs (the offending name plus
the list it should have matched).

Checking passes return diagnostics; the [stratum](#g-stratum) barrier (a stratum
is one of the build's three phases: structure, schedule, activation) throws a
single `BuildError` wrapping the collection. `showerror` renders that
`BuildError` compiler-style, grouped by kind and sorted by path.

```julia
# a diagnostic value: its kind is its identity, its payload is plain data
struct WireTypeMismatch      # one kind of the closed Appendix C set
    …                        # payload fields per Appendix C: paths and names
                             # as Strings, expected/observed port types as types
end

# the carrier: one exception, thrown once at the stratum barrier
struct BuildError <: Exception
    …                        # the collection it wraps
end
```

A user-code exception is wrapped in a framing diagnostic — [component](#g-component)
path, which function, the [probe](#g-probe) context including synthesized inputs
— with the original exception as `cause`. The didactic frame therefore renders
first and the raw throw second.

One class is recognized rather than merely framed. A `FieldError` carries its
type and field as data. Matched against the [bundle](#g-bundle)'s own NamedTuple
type — the NamedTuple of zero-copy views a component function receives — it
becomes the bundle-law did-you-mean ([§5.2][s5-2]), carrying the legal set and
the undeclared-store / wrong-[tier](#g-tier) / illegal-for-this-function
classification. Nothing is recovered by reading message text.

The [§11.4][s11-4] walkthroughs as acceptance tests target diagnostics: tests match on
kind plus payload fields, never on message text. Messages become pure
presentation.

Two rendering rules are doctrine, not style:

- **Strings, never instances.** Diagnostics carry paths and names as strings,
  never component instances and never model types — the `compact_backtrace`
  lesson. Expected/observed *[port](#g-port)* types are the payload exception, and they
  are small: `Float64` vs. `Bool`, a NamedTuple field diff.
- **The didactic [register](#g-register) is policy.** Every diagnostic states the fix or the
  lists-in-hand, not just the violation: "return `zero(x.ω)`, not `0`"; "no
  input `throtle`; did you mean `throttle`?"; the child's face list alongside
  the unknown `except` entry.

**Two [warning streams](#g-warning-streams), scoped separately.** The *build* diagnostic stream is
the one the build kinds ([Appendix C][sC]) ride: warnings there carry warning severity,
render with the collection, and never trigger the throw. That stream's warning set
is **currently empty**, the unconnected-output warning having been rejected as
its sole candidate ([§6.1][s6-1], row 84). Better an empty, trusted stream than a
noisy one; a warnings-as-errors CI switch is addable, not built.

The *runtime* status/log stream is a different channel, and it is not empty. It
is per-occurrence, carried by the per-writer [diagnostic cells](#g-diagnostic-cell) ([§9.8][s9-8]) — the
single-writer cell each writer owns for diagnostics and heartbeat. The cells are
where the rate limit lives, as a structural bound (a bounded ring plus per-kind
suppressed counts, drained at frame top) rather than a policy layered over the
stream. So "rate-limited wherever its source can repeat" holds of every kind
below without any kind arranging it. The stream surfaces through the published
[framework status](#g-framework-status) ([§9.2][s9-2]) alongside the [§8.7][s8-7] pacer diagnostics and the
[§10.2][s10-2] liveness heartbeats, which ride in the same [cells](#g-cell). It is never
collected, since there is no collection to join. Nothing in the argument
(row 84) applies to it: that decision is about what the *build* warns on.

A *service* warning (`TrimCommitEvents`, `TrimCommitResiduals` —
[Appendix C][sC]) is neither stream. It is a synchronous per-call annotation,
emitted once at a stopped-sim service call's return beside the value it returns,
its payload duplicated as plain report fields. There is no carrier cell, no
collection, no rate limit to arrange. The committed runtime warnings, in one
place:

- **[chattering](#g-chattering) / localization-budget exhaustion** ([§8.4][s8-4]) — a [localized](#g-localized) event
  whose bracketing budget runs out at a [boundary](#g-boundary);
- **firing-budget exhaustion** ([§8.6][s8-6]) — an event that has spent its
  `firing_budget` at a boundary, its further edges there dropped;
- **forgiven-debt re-anchor** ([§8.7][s8-7]) — the pacer abandoning accumulated debt
  and re-anchoring its [schedule](#g-schedule);
- **the write-surface and entry violations** ([§9.3][s9-3]), all at staging — the
  [drain](#g-drain) checks nothing: `OutOfClaimEntry` (an
  enumerated surface's binding drift — no position in
  the attach-compiled schema), `ClaimedFaceEntry` (a harness write to
  a face claimed in the run's frozen partition, naming the incumbent; also
  fired by the stopped-sim attach renormalizing a pending [batch](#g-batch)) and
  `EntryTypeMismatch` (a value unconvertible to its
  [slot](#g-slot)'s declared type, rejected at staging for every writer);
- **a tolerated [device](#g-device)-side datum failure** ([§9.6][s9-6], [§13.4][s13-4]):
  `MalformedDatum` — emitted by the author's loop via `report!(handle, …)`
  into the device's own cell ([§9.8][s9-8]), the `InputMappingError` successor;
- **staging discarded during [replay](#g-replay)** ([§10.7][s10-7]): `ReplayDiscardedStaging` —
  a live batch found staged while the [trace](#g-trace) feeds the drain;
- **thread-budget tightness** ([§10.2][s10-2]) — once per `run!`, against the frozen
  [roster](#g-roster);
- **device join timeout** ([§10.4][s10-4]) — a device task exceeding the shutdown
  join timeout, abandoned by name rather than hanging `run!`;
- **device crash** ([§10.4][s10-4], [§13.4][s13-4]) — a device task's failure caught by the
  framework wrapper, the sim continuing with the device absent;
- **unbounded run** ([Appendix B][sB]) — no finite `t_end`, no `stop_on` faces,
  `pace = Inf` at run start.

### 13.3 Build primitives: `resolve` and the face-list accessors

The `input_passthrough` sketch ([§11.8][s11-8]) calls two of the primitives
below without defining them. All three are normative in the forms given here:

- `resolve(asm, path::String) → AbstractComponent` — the getfield walk along
  `/`-segments.
- `input_faces(c)` / `output_faces(c) → Vector{String}` — the stringified keys of
  a leaf's `input_types` (the key set is `T`-independent); the entries of
  `input_connections(c)` / `output_connections(c)` for an
  [assembly](#g-assembly). Declaration order is preserved: deterministic
  printouts, stable diagnostics.
- `resolve_terminal(asm, path) → (component, name)` — splits a terminal path's
  final segment and resolves the prefix through `resolve`. The split is
  unambiguous because [face](#g-face) names may contain dots, never slashes
  ([§11.6][s11-6]).

**Enforcing the generic-[boundary](#g-boundary) rule is the one non-obvious duty
in `resolve`** ([§6.1][s6-1]). The walk follows *declared field types* alongside
instances, and a segment that traverses **past** a generically-held field — one
whose declared type is non-concrete — is a diagnostic even though the concrete
instance in hand would resolve it. Resolving *to* a generic child is
[port](#g-port)-level access and legal. An unknown segment errors with the
sibling field list in hand.

**The duty is [register](#g-register)-scoped.** The load-bearing/diagnostic line
(row 83) is carried into resolution. Client policy rides on one primitive, the
same arrangement as the two application registers over one plan
([§14.4][s14-4]).

| register | who resolves under it | what the walk enforces |
|---|---|---|
| **structural** | wiring resolution, in [Stratum](#g-stratum) A (one of the build's three phases: structure, schedule, activation) | the strict rule above, verbatim |
| **load-bearing** | [condition](#g-condition) entries (the path-addressed sparse overlay that sets a build's state), trim `reads`, [taps](#g-taps) ([§14.3][s14-3], [§14.7][s14-7], [§14.10][s14-10]) | strict, evaluated **at the authoring or mount level** |
| **diagnostic** | [device](#g-device) read [bindings](#g-binding), GUI panels, [snapshot](#g-snapshot) and log inspection ([§9.2][s9-2], [§9.7][s9-7]) | the instance walk |

Each register's treatment has its own warrant. The structural register is the
one the law ([§6.1][s6-1]) lives in, so it applies that law verbatim. The
load-bearing register evaluates at the authoring or mount level for two
reasons: the locality law is an authoring-level law, absolute paths being a
compiled derivative ([§14.2][s14-2]); and a mount prefix is checked by the
mount itself, which validates the problem's declared type against the mount
point's [contract](#g-contract) ([§14.9][s14-9]). So this register checks the
authored path below that prefix. The diagnostic register walks instances
instead. A generic [seam](#g-seam) is not an error for a client that never
claimed substitutability: "what is in *this* build" is the inspection
register's defining question. Drift still stays loud — an unknown path is an
attach-time `ReadBindingUnresolved` with [did-you-mean](#g-did-you-mean) (the
offending name plus the list-in-hand it should have matched).

Which register a client resolves under is internal framework fact, never
user-facing API — the same status as the two `apply!` registers
([§14.4][s14-4]).

`resolve_terminal` is first-class because five clients share it across the three
registers: wiring resolution (structural); condition addressing
([§14.3][s14-3]) and tap resolution ([§14.10][s14-10]) (load-bearing);
device-binding validation ([§9.2][s9-2]) and snapshot inspection (diagnostic).
The result is one splitter, one did-you-mean site.

### 13.4 Runtime failures: one catch site, an execution cursor

**Where caught.** The loop wraps each execution of the [boundary](#g-boundary) macro-sequence
(integrate → project → event iteration → [ticks](#g-tick) → publication) in a single
`try`, never per stage or per [component](#g-component) (row 59). Framing information does
not need to be *caught* into existence. The [executor](#g-executor) (the compiled execution
form of the [schedule](#g-schedule)) maintains an **[execution cursor](#g-execution-cursor)**, a plain mutable
field in the loop state recording where in the compiled schedule execution is.
The cursor records three facts. The first is the component path, as a schedule
index. The second is which function is running: `h_x`, `h_xu`, `f`, `g`, a
[guard](#g-guard), a handler, or `project`. The third is the boundary phase:
integration stage *k*, event round *r*, localization [probe](#g-probe) at trial time, or
tick. Maintaining it costs one cheap store per dispatch on a single-tasked
executor: no allocation, no exception frames. And it covers every user-code
surface uniformly, including the forgettable ones: [RHS](#g-flow) evaluations at
interior RK stage points, guard evaluations inside ITP/Brent probes,
environment closures.

```julia
# the cursor: one mutable field of the loop state, overwritten per dispatch
mutable struct ExecutionCursor
    …    # what it records, per the prose above: component path (schedule
         # index), which function, and the boundary phase
end
```

**How handled.** The catch site wraps the original exception in `StepError`, the
runtime counterpart of `BuildError`. A `StepError` carries four things: the
cursor's [frame](#g-frame), the boundary time, the **frame-entry boundary index**, and
the original exception as `cause`. The frame-entry boundary index is the
[replay](#g-replay) pointer: the frame-top boundary at which the failing frame began.
That frame top is a grid boundary or [boundary zero](#g-boundary-zero) (the initialization
boundary: the ordinary macro-sequence with an empty integrate), and it is always
a legal replay halt ([§10.7][s10-7]). A `StepError` is rendered with compact frames
per the doctrine ([§13.2][s13-2]).

Conformance failure ([§12.5][s12-5]) needs no separate path. It is thrown as its
typed diagnostic at the table-write point, and it arrives at the same catch
site. There it is a species of `StepError` carrying the field-diff [payload](#g-payload).

Reproducibility holds by construction. Staged inputs are drained and recorded to
the [trace](#g-trace) at the frame top, *before* the boundary executes. So the failing
boundary's inputs are already in the trace when it fails. The error names the
frame-entry boundary `k` to replay to; `replay!(sim2, trc; to_boundary = k)`
halts exactly at that frame top. `step!` then re-executes the failing frame
under instrumentation, [localized](#g-localized) boundaries included ([§10.7][s10-7]).

```julia
replay!(sim2, trc; to_boundary = k)   # halt at the frame top the error names
step!(sim2; frames = 1)               # re-execute the failing frame, instrumented
```

**The one exception never wrapped.** An `InterruptException` is not model code
failing; it is the operator's stop command ([§10.4][s10-4]). So the catch site
discriminates it and routes it to the stop path. The run takes the ordinary
graceful tail and ends `stopped`, never `errored` under a `StepError`. With the
boundary masking in force ([§10.4][s10-4]) the branch is unreachable in practice: the
interrupt is deferred to a frame-top or wait unmask point, and never raises
inside the guarded sequence. It is kept defensively, because the cost of being
wrong about that is a terminally errored session in place of a clean stop.

**Disposition.** The `Simulation` ends in a terminal status, `stopped` or
`errored`, with the exception retrievable. A synchronous [unattended run](#g-unattended-run) (a run
with empty staging and no snapshot readers) rethrows after the shutdown tail
completes, so CI fails honestly. An interactive session logs the rendered error
and surfaces the status through the control plane and GUI.

**The nonfinite check.** Divergence is not termination. Dynamics that blow up
(ground penetration, an unstable gain) produce NaNs that defeat guards. NaN
comparisons are false, so no declared condition will catch them. A loop-level
`isfinite` [sweep](#g-sweep) over `x` at boundaries fails fast as a `StepError` species,
naming the offending component's state block and the boundary. It catches
diverging models generally, not just post-terminal ones.

*Placement is the whole value.* The sweep is the boundary's **first act** —
immediately after integrate returns, before `project` and before the boundary
sweep. Run there, `NonfiniteState` names the component whose own block
diverged. Run later, the NaN has already propagated: it reaches an innocent
downstream component through the ordinary signal path and surfaces as that
component's lookup-table `DomainError`, or as an `InexactError` in its
conversion. That is the error-locality inversion ([§11.4][s11-4]), designed out of the
build [tier](#g-tier) and quietly reintroduced at runtime. One `isfinite` pass over a
flat [buffer](#g-buffer) is cheap enough that placement, not cost, decides.

*Scope: `ẋ` does not participate* (row 157). A nonfinite derivative
contaminates its own state block's step result within that very step, so the
`x` check at the next boundary is the same detection with identical component
attribution. And `ẋ` buffers are integrator scratch: written per stage,
meaningful only inside a step, and not boundary-consistent in the sense the
check is stated over.

**Domain separation.** [Device](#g-device)-side user code fails in the device's own
domain: loop bodies and mappings run on the device task ([§9.4][s9-4], [§9.6][s9-6]). It
fails in two classes. A genuine bug takes the per-device crash path (liveness
heartbeat, `DeviceCrash`) while the sim keeps running. An unmappable datum is
not a failure at all; the loop body tolerates and reports it (`MalformedDatum`,
[§9.6][s9-6]). The two failure domains never mix — exactly what the
no-shared-mutable-model decision bought.

### 13.5 Termination is a state, not an exception

FlightCore's `SimulationTermination` idiom — model code throws, the loop
catches and logs it as informational — has **no counterpart here** (row 60).
The discipline: **exceptions from model code are always abnormal**. Graceful
termination is model *state*, reaching the loop through declared machinery:

- **Detection** is ordinary [guard](#g-guard)/handler/mode machinery. Declare
  the [predicate](#g-predicate) as a sign-form event — hence
  [localized](#g-localized), the crossing instant bracketed by root-finding
  over probe sweeps — if the stop should be localized. Touchdown overload is
  precisely a zero-crossing: the boundary is localized to the crossing, the
  handler sets `m.crashed`, and the [snapshot](#g-snapshot) at the crossing
  instant carries the touchdown state.
- **Publication** is an ordinary `Bool` output [face](#g-face), exported to
  the root. Within concretely-declared structure, deep wires gather the
  condition at its owning boundary in one visible block ([§6.1][s6-1]): `Ldg`
  ORs its three legs through a junction (the ownership idiom, [§6.2][s6-2];
  the library, [§13.7][s13-7]) and exports one `damaged` face; intermediate
  [assemblies](#g-assembly) are untouched. Each *generic* [seam](#g-seam)
  costs one output connection entry. That hop is the substitutability
  [contract](#g-contract) doing its job, not plumbing (the imposed contract,
  [§11.8][s11-8]).
- **Policy** binds at deployment: `Simulation(world; ..., stop_on = (...))`
  names root-exported `Bool` output faces. They are OR-combined, validated
  against the `Build`, and recorded in the [run metadata](#g-run-metadata) —
  the [trace header](#g-trace-header)'s deployment block ([§9.5][s9-5]). After
  *every* published boundary — grid, `t*` ([§8.4][s8-4]) and
  [boundary zero](#g-boundary-zero) ([§14.5][s14-5]) alike — the loop reads the
  named faces in the snapshot it just published. The first `true` initiates
  [§10.4][s10-4] shutdown with *this* snapshot as the final one: the terminal
  snapshot is the terminal state, no roll-back, nothing [§10.4][s10-4] doesn't
  already do. That terminal snapshot's status carries the run's final cumulative
  diagnostic counters ([§9.8][s9-8]). `run!` therefore checks the boundary-zero
  snapshot before the first step: an authored [condition](#g-condition) (the
  path-addressed sparse overlay that sets a build's state) already terminal
  ends the run at `t₀` with that snapshot final, integrating nothing. The
  default is no stop faces and a run to `t_end` — `stop_on` is `t_end`'s
  model-declared sibling at the same declaration site.

**Both are `run!`-time overridable, with the constructor value as the
default.** `Simulation(world; t_end, stop_on)` sets the defaults for the
simulation; `run!(sim; t_end = …, stop_on = …)` binds them for **that run
only**. The `run!` argument wins where given, and the constructor's value
stands where it is not. Nothing about the `Simulation` is mutated, so the next
`run!` without arguments gets the constructor's policy again. The effective
values are what the run metadata records. `stop_on` face validation against
the `Build` runs at **both** binding sites, identically: an unknown or
non-`Bool` face fails at `run!` exactly as it fails at construction.

```julia
sim = Simulation(world; h = 0.02, t_end = 60, stop_on = (…))  # the Simulation's defaults

init!(sim, cond); run!(sim)                # t_end = 60, stop_on as constructed
init!(sim, cond); run!(sim; t_end = 10)    # this run only: t_end = 10, stop_on as constructed
init!(sim, cond); run!(sim; stop_on = ())  # this run only: no stop faces, t_end still 60
init!(sim, cond); run!(sim)                # the constructor's pair again: nothing was mutated
```

This is not the root-declared stop policy rejected below: the `run!` argument
moves binding one notch *later* along the same axis, more deployment-flavored
rather than less (rows 60 and 91). The `stopped → init! → run!` cycle and the
`step!` register ([§10.6][s10-6]) are precisely where one `Simulation` wants
different stopping policies on different runs. The honest cost is two homes for
one fact; the precedence rule above settles it.

**The termination record names the source.** Where run metadata carries the
effective *policy*, the run's termination record carries its *outcome*: which
of the four sources ended the run. Those sources are `t_end` reached, a named
`stop_on` face reading `true`, a control-plane stop (GUI button,
[device](#g-device) handle, calling code), and an
[operator interrupt](#g-operator-interrupt) ([§10.4][s10-4]). A `stopped`
simulation therefore answers "why did it stop?" without its consumer
reconstructing the answer from the clock. The interrupt is a tag on an ordinary
stop, not a diagnostic [kind](#g-kind) of its own ([Appendix C][sC] gains
nothing here): nothing failed.

Taught contract: **stop faces are sampled at completed boundaries; declare a
sign-form event if you need the stop localized.** Both stop-flag shapes work
without framework latching. A handler-set `m` flag is sticky by nature, and a
transient stage-2 Bool is caught because the loop reacts to the first `true`.
Compound stop logic composes in-model — a monitor [component](#g-component)
reading the relevant signals and outputting one Bool — the same move
[§10.5][s10-5] made for scripts.

Post-terminal dynamics are the model's job, and that is a feature. Today
`robot2d` *throws* when it falls, because it has no other way to say "my
dynamics are no longer meaningful". Here it declares the fall as an event,
switches to a frozen mode — mode-dependent `f`, machinery it already has — and
exports `fallen`. Wired, the sim ends at the fall; unwired, it integrates a
frozen robot — well-defined, unlike an uncaught throw. The discipline forces
models to have well-defined terminal states, which is better modeling.

Rejected mechanisms — predicate closures (`stop_when = snap -> ...`),
root-type-declared stop policy, [blessed](#g-blessed) terminal types and
`terminal` event flags, a [control-plane](#g-control-plane) capability for
[components](#g-component), and observation-by-path (`stop_on` naming a deep path
into any public output) — are litigated in row 60. A root-declared *default*,
overridable at the constructor, is the one variant on record for reopening,
should the constructor argument prove chronically forgotten ([§16][s16]).

The observation-by-path line leaves doctrine behind it. **Diagnostic
observation** — the log retaining the full table, GUI panels rendering a
component's [ports](#g-port), [replay](#g-replay) inspection — is human-facing,
has no effect on run semantics, and legitimately sees every public
[cell](#g-cell). **Load-bearing observation** — a read that changes what the run
*does* — must speak the [contract](#g-contract). `stop_on` is the one read that
changes what the run does, which is why it alone names root-exported faces;
output devices are the other half of the same doctrine, their reads being
diagnostic snapshot-path bindings ([§9.2][s9-2], [§15.4][s15-4]).

The wall-clock channel — GUI stop button, device handle, code — is orthogonal
and untouched: that is the [control plane](#g-control-plane)'s operator path. The sim-time,
model-detected channel specified here meets it at [§10.4][s10-4] and nowhere else.

### 13.6 Abnormal shutdown: one tail, two entries

**The [boundary](#g-boundary) is all-or-nothing outside the sim task.** That is
why a `StepError` cannot break [§10.4][s10-4]. [Sweeps](#g-sweep) write into
table [buffers](#g-buffer), and integration intermediates live in
framework-owned integrator buffers, never in a [component](#g-component)'s
`workspace` as [§7.3][s7-3] defines it. The only externally visible act is
[snapshot](#g-snapshot) publication at the very end of the sequence. A boundary
that throws has published nothing. The last *published* snapshot — a complete,
consistent boundary by construction — is still the newest thing any
[device](#g-device), logger or waiter has seen.

The abnormal path is therefore: **discard the failed boundary, promote the
previous snapshot to final, and rejoin the ordinary tail.** The protocol becomes
one tail with two entry points: graceful entry after a *completed* final
boundary, abnormal entry after a *discarded* one. Everything downstream of
"final snapshot" runs identically: sticky stopped, waiters woken through the
boundary-counter + `Condition` path, `unblock!`/close hooks, named joins with
timeout. Those waiters observe stopped rather than a new boundary, so no device
task hangs.

```
  graceful entry                      abnormal entry
  final boundary completes            failed boundary discarded
  final snapshot published            previous snapshot promoted to final
             |                                   |
             +-----------------+-----------------+
                               |
                       "final snapshot"
                               |
       sticky stopped → waiters woken → unblock!/close hooks
                     → named joins with timeout
```

This fills the seat [§10.4][s10-4] reserved for it: a loop-side failure runs
the same protocol from the catch path.

Tail hygiene: the hooks are user code too, so each is individually
caught-and-logged. Shutdown therefore runs to completion even if a device's hook
misbehaves. The join timeout already bounds a hook that hangs rather than
throws.

What is lost is quarantined: the state [stores](#g-store) may hold mid-boundary
values (a half-written `m`, integration intermediates). They are retained on the
errored `Simulation` for post-mortem inspection, but an errored sim is
terminally stopped, not resumable. The reproduction tool is [trace](#g-trace)
[replay](#g-replay), not resurrection, and the stopped-sim services enforce that
non-resumability by refusing an errored simulation outright
(`ServiceLifecycle`, [§14][s14]). The published record (snapshot chain, log,
trace) ends at the last consistent boundary. Nothing downstream of the sim ever
sees half a boundary.

### 13.7 Tooling consequences: provenance and the component library

Termination chains are the second structural customer of computed
[boundary](#g-boundary) connections, after generic-boundary
[contracts](#g-contract). Computed connections are therefore prominent in this
section, and two commitments follow: a library and an idiom.

**The `input_passthrough` helper family grows deliberately.** Predicate-based
selection — an `endswith`-style filter alongside `except`/`only` — is a natural
extension: still explicit at the declaration site, still evaluated at build,
still printable. That is the [blessed](#g-blessed) side of the auto-bubbling
line, where the author writes down the *rule* and the build evaluates it into
inspectable data.

**The `Build` printer owes [face](#g-face) provenance.** For every root face,
that means the resolved chain down to the producing terminal (`"crashed" →
aircraft/monitor/out ← systems/ldg/{left,right,nose}/damaged`). Once faces are
computed rather than hand-listed, "what does this face actually reach" is a
question the artifact must answer, not the reader. The same rendering serves the
wiring diagnostics, which already carry endpoint paths.

#### The standard component library

**A standard [component](#g-component) library makes good on the junction
promise ([§6.2][s6-2]).** That promise rests on explicit junctions being
*cheap*, and a junction hand-written per arity per type is not.

The starting inventory comes strictly from demonstrated need: wrench/scalar
summing junctions, the Bool gates the termination chains use, `UnitDelay` — the
spelling the second loop-breaking remedy ([§5.5][s5-5]) needs — and
`Constant{V}`, the source block. The library grows by migration demand only:
Simulink's library is a language; this is a toolbox.

One member is admitted by persona rather than migration demand: `Group`, the
on-the-fly [assembly](#g-assembly), whose declaration-layer treatment lives in
[§11.5][s11-5] (row 184). `Group` serves the model assembler, for whom topology
is data rather than a named type, and it needs no rule the declaration layer
does not already have.

Doctrine: **library blocks are ordinary components** — no framework privileges,
no special vocabulary. Ordinary status is what keeps [schema
authority](#g-schema-authority) total (declarations define structure;
evaluation only checks conformance). That same status makes the library a
permanent ergonomics [torture test](#g-torture-test): if a three-input OR gate
is painful to write under the declaration rules, the rules are wrong.

Arity comes from a type parameter: `Or{N}` builds `(in1 = Bool, …, inN = Bool)`
programmatically — a derivation [§11.2][s11-2] blesses, and an early validation
that the contract functions support parametric components.

[Tier](#g-tier)-transparency falls out of settled semantics. A stateless
continuous `h_xu` recomputes every [sweep](#g-sweep), so fed ZOH-held discrete
signals its output changes only at [ticks](#g-tick). No tier-neutral class is
needed.

`UnitDelay{V}` is a **discrete** leaf at `K = 1`: the tier's native `z⁻¹`
([§8.6][s8-6]), the shape `sat_out_0` hand-writes in [§15.2][s15-2], and it
needs no framework support.

```julia
# UnitDelay{V} — a discrete leaf at K = 1; port face names elided
init_x(::UnitDelay{V}) where {V} = (v = zero(V),)
h_x(::UnitDelay, (; x)) = (; … = x.v)     # publishes the stored value
g(::UnitDelay, (; u)) = (; v = …)         # stores the incoming one, from u
```

`UnitDelay`'s tier semantics are the point, and they must be stated wherever the
remedy is recommended. Inserting a `UnitDelay` into a *continuous* loop moves
that signal onto the discrete tier and inserts a `Δt_base`-scale ZOH into the
model's mathematics. That is a modeling decision, the delayed signal being
genuinely sampled rather than a transparent wire. The diagnostic
([§5.5][s5-5]) therefore says so rather than offering the remedy as free.

`Constant{V}` is the **source block**: no inputs, no state, and a stage-1 body
returning the value the instance holds.

```julia
# Constant{V} — a stateless continuous leaf; its value is instance data
output_types(::Constant{V}, ::Type{T}) where {V, T <: Real} = (out = V,)
h_x(c::Constant, _) = (; out = …)         # the value the instance holds
```

`Constant` is a stateless continuous leaf, so the tier-transparency argument
above already covers discrete consumers and no discrete variant is needed. The
declaration takes the [activation](#g-activation) scalar (a re-run of Stratum C
at a given scalar type) and ignores it. That is the point: the block's output
*is* its stored value, so the leaf is **deliberately [pinned](#g-walked)** at
that value's own type. A `Constant{Float64}` declares a `Float64` port and means
it, and the embedding ([§11.2][s11-2]) turns it into the zero-partial constant
it already was under any `Dual` activation. The honest pin is spelled rather
than inferred.

Two demonstrated needs admit `Constant` under the inventory's demonstrated-need
charter. The zero-contributor configurations ([§6.2][s6-2]) are one: a required
aggregate input has no physical contributor, and the zero total must be spelled
as a wire. The rig stub below is the other. The block's value is instance data,
like junction arity, not an overridable default: a configuration wanting an
externally settable source uses a root [slot](#g-slot) ([§9.3][s9-3]), which
keeps the block from drifting into a back-door input default. A migration-phase
deliverable.

#### The component test rig

**The [component test rig](#g-component-test-rig) is the library's companion
idiom.** It is a one-child assembly whose `input_connections` surface the
child's entire input face set — `input_passthrough(rig, "child")` verbatim
([§11.8][s11-8]). Any component can therefore be built and simulated in
isolation: every input becomes a root slot fed by ordinary conditions and
[devices](#g-device), and every output is observable in the
[snapshot](#g-snapshot) table.

One qualification comes from the root-slot rule ([§11.2][s11-2]). An *abstract*
input entry (`terrain = AbstractTerrainField`, [§4.4][s4-4]) cannot surface as a
root slot, because abstract-at-root is a build error. The rig therefore
satisfies that entry *inside* the rig: a concrete stub child (a
`SampleTerrainField` provider) wired to the face, with the concrete remainder
exposed via `input_passthrough(rig, "strut"; except = ("terrain",))`.

```julia
struct StrutRig <: AbstractComponent    # the rig: component under test + stub
    strut::Strut
    stub::Constant{SampleTerrainField}  # the test handle, held as instance data
end

child_connections(::StrutRig) = ("stub/out" => "strut/terrain",)
input_connections(rig::StrutRig) =
    (input_passthrough(rig, "strut"; except = ("terrain",))...,)
```

That stub child is typically just a `Constant` holding the test handle — the
source block's first shipped instance. Bespoke stubs remain ordinary components
wherever the double must compute something.

The rig adds zero new machinery: wiring and `except` already exist. It is the
substitutability contract doing its job. An [abstract entry](#g-abstract-entry)
(an `input_types` entry admitting any concrete producer face) declares that a
substitute must be chosen, and the rig choosing its test double explicitly, as
ordinary inspectable code, is precisely the isolation the rig exists to provide.

The rig is [`design_world`](#g-design_world)'s little sibling
([§14.9][s14-9]): where `design_world(ac)` mounts an aircraft in a minimal world
for trim and linearization, the rig mounts one component behind a root contract
for unit tests and open-loop probing. The machinery is deliberately ordinary end
to end, with no framework support. That makes the rig, like the library blocks,
a standing ergonomics test of the declaration rules.

---

## 14. Stopped-sim services

[§12.6][s12-6] previewed the services as [Stratum](#g-stratum)-C clients:
initialization, trim, linearization and [capture](#g-capture) (reading the
current stores and slots back as a condition). Everything they share reduces to
one artifact: the **[condition](#g-condition) value**, the datum that says "set
this build to this state."
[§14.1][s14-1]–[§14.4][s14-4] settle its representation, composition and application;
[§14.5][s14-5]–[§14.6][s14-6] the [boundary](#g-boundary)-zero sequence and [slot totality](#g-slot-totality)
(the requirement that an application cover every root slot); [§14.7][s14-7]–[§14.9][s14-9] the
trim service in full; [§14.10][s14-10] linearization and `capture`.

**Lifecycle preconditions.** Every service requires a non-running simulation:
while a run exists the loop owns the [stores](#g-store) between [drains](#g-drain)
(the frame-top swap that publishes staged device inputs into the root slots),
and a service reading or writing them would race it. Pause is no exception, by
the doctrine that freezes the [roster](#g-roster): pause is a control-plane state
*inside* a run, so a prohibition that holds mid-run holds while paused
([§9.3][s9-3], [§10.1][s10-1]).

Within the stopped-sim states, legality follows each service's inputs:

| service | `built` | `initialized` | `stopped` | its inputs |
|---|---|---|---|---|
| `capture` | error | legal | legal | committed, boundary-consistent stores |
| `init!` | legal | legal | legal | authored conditions |
| `trim!` | legal | legal | legal | authored conditions; the scratch world is [`override`](#g-override)`(baseline, condition(guess))` ([§14.8][s14-8]), never the sim's stores |
| `linearize`, operating point defaulted to `capture(sim)` | error | legal | legal | inherits `capture`'s precondition |
| `linearize`, explicit `about` ([§14.10][s14-10]) | legal | legal | legal | inherits `init!`'s legality — legal wherever `init!` is |

**`errored` is terminal for all four** (rows 59, 108). Post-mortem inspection
of an errored sim's stores, log and [trace](#g-trace) stays available as a
diagnostic read; it may not become a condition value.

A violation is `ServiceLifecycle` ([Appendix C][sC] — the
operation, the current status, the legal statuses), the same kind
`attach!`/`detach!` raise while `running` ([§9.3][s9-3]): one [register](#g-register) for "this
operation is illegal in the current lifecycle state," distinct from
`MissingInit`, which names a missing prior step.

### 14.1 Conditions are path-addressed overlays on the declared defaults

A [condition](#g-condition) may specify state fields (`x`, on either [tier](#g-tier)) and
modes (`m`, [continuous components](#g-continuous-component) only, [§3.2][s3-2]); both are
addressed by [§11.6][s11-6] slash path plus field name. It may also specify root
input [slots](#g-slot), addressed by [face](#g-face). Never outputs, which are derived data.
Never [workspace](#g-workspace) (component-declared mutable scratch arriving as the `ws`
bundle field). Entries are validated in the [§13.1][s13-1] collecting [register](#g-register): full
list, violations collected, one `BuildError`.

**The overlay base is always the declared defaults.** Every [store](#g-store) has a
declared initial value (declaration-by-initial-value, [§11.2][s11-2]), so conditions
are naturally sparse. Applying one means "fresh run from the `init_*` defaults,
with these overrides" (row 63). Warm restart needs no second semantics. A
`capture` service reads the current stores **and root slots** back *as a
condition value* (capture → tweak → apply). Slot coverage is what makes the
captured condition total, hence re-applicable under [§14.6][s14-6]. That gather is
the one the [trace header](#g-trace-header) already needs: one mechanism, two uses.

**Doctrine.** Addressing conditions by path does not reopen the
observation-by-path rejection ([§13.5][s13-5]). That was *runtime* coupling: a
root-authored predicate reaching through generic [seams](#g-seam) the root does not
own, breaking on substitution. A condition is a *design-time statement about a
concrete build*, authored in the same register as `child_connections` (which
also speaks paths, about children its author owns). The composition law
([§14.2][s14-2]) makes the parallel exact.

**Pre-[sweep](#g-sweep) doctrine.** Condition writes precede the first sweep by
definition. A would-be init value that depends on swept outputs is therefore
either analytically known to the caller or an equilibrium constraint. The first
case is trim's `α_filt = α_a`: α is a *decision variable*, so the value is known
above, not computable below. The second is a job for the trim service, not for
init.

"Caller-computable" reaches past closed-form knowledge to **environment
queries**. A condition needing one constructs the same handle the sweep will
produce, and then calls the same query function the consuming component calls.
One route to that handle is the [value-level constructor](#g-value-level-constructor) (the plain
exported function building a field handle from input values, [§4.4][s4-4]), applied
to the same values the [`baseline`](#g-baseline) writes into the environment
[component](#g-component)'s slots. The other applies in a rig where the handle
itself is a root slot value: the condition simply holds the value the `baseline`
wrote there. One implementation of the field math, evaluated one level up: no
pre-sweep, no new mechanism. Where closed-form enforcement of a target is not
wanted at all, the second escape already covers the case: promote the eliminated
state coordinates to decision variables and enforce the targets as residuals on
swept outputs. That route needs no environment access at condition time
whatsoever.

### 14.2 Fragment composition: locality without schema

Init knowledge is [component](#g-component)-local — the engine knows
`n_eng → ω = n_eng·ω_rated`, and nothing above it should have to. Making that
locality a schema entry was rejected (row 64). That spelling would declare
`initialize(::C, spec)` — today's `f_init!` reborn declaratively — and add an
[assembly](#g-assembly)-level rule routing sub-specs to children.

What preserves the locality is an idiom, not schema: **[fragment](#g-fragment) functions**,
ordinary functions shipped beside the component and dispatched on the component:

```julia
condition(eng::PistonEngine; n_eng) =
    fragment(x = (ω = n_eng * eng.ω_rated,), m = (phase = Phase.running,))
```

Fragments are composed by *pull* from the structure's owner:

```julia
condition(sys::C172XSystems; n_eng, α_a, β_a) = merge(
    at("pwp/engine", condition(sys.pwp.engine; n_eng)),
    at("aero",       fragment(x = (α_filt = α_a, β_filt = β_a))))
```

Dispatch selects variant-specific methods, so the c172s/c172x actuation
split costs no upstream edits. The three combinators are constructors of an
**inert, lazy tree**: no path arithmetic happens at composition.

```julia
struct Fragment{X,M,S}  x::X; m::M; slots::S  end          #self-vocabulary payloads; no paths
struct Scoped{N}  prefix::String; node::N  end             #at(prefix, node): stores, never applies
struct Merged{T<:Tuple}  nodes::T  end                     #merge(ns...): collects; order = diagnostics only
```

Every node is isbits except the interned literal prefixes, so **rebuilding
the tree per trim iteration is stack-only construction** — the zero-alloc
property of today's `assign!` loop, preserved structurally.

`fragment`'s payloads speak only about the component at the authoring point;
addressing children is exclusively `at`'s job (one way to say everything). A
`slots` payload names faces *of the authoring level's [contract](#g-contract)*.
Resolution walks the export chain to the root slot and errors if the face never
surfaces. An internally-wired input has no slot, and writing it would be
meaningless because the first sweep overwrites it. Unexported stays unpokeable
for init exactly as it does for the GUI ([§9.7][s9-7], [§15.4][s15-4]).

**The locality law** here is the one [§6.1][s6-1] states for connections, now in
its third instance — child connections, computed [boundary](#g-boundary)
connections, conditions. Each level speaks its own fields, its declared
children's names, and its own faces; delegation runs by dispatch at every
genericity [seam](#g-seam); deep `at` paths are legitimate exactly where deep
connections are, within an owned concrete subtree. Absolute paths exist only in
the flattened entry list, a *compiled derivative* of the composition, as slot
offsets are of `child_connections`. Substituting a component invalidates
precisely the fragments its owner shipped, nothing else. The enforcement status
carries over from [§6.1][s6-1] as well: the law is convention. Ownership is a
fact about who maintains the code and the build cannot see it, so the law is
available and idiomatic rather than machine-checked. `fragment`/`at`/`merge`
are [§13.7][s13-7] standard-library material — ordinary artifacts, no privileges.

A [merge](#g-merge) collision is two entries on one leaf. Collisions are errors
at resolution, and the error reports *both* provenance chains. The message
names the layering combinator: "`merge` is collision-intolerant by design — use
`override(base, patch)` to layer." Last-writer-wins was rejected (row 65). That
rejection is exactly where this `merge` parts company with `Base.merge`, which
is last-wins on NamedTuples; the two share a name and not a semantics, and
dispatch on the node types keeps them apart mechanically. The mixed call is
closed explicitly. A `merge(::Fragment, ::NamedTuple)` — or any other blend of a
condition node with a bare NamedTuple — is an **error method**, defined so the
call cannot fall through to `Base.merge`'s last-wins semantics on a payload that
looks plausible. Its message is directive: wrap the NamedTuple in
`fragment(...)` (or `at(prefix, fragment(...))`) and merge nodes with nodes. The
rejection carries a [kind](#g-kind) like every other: `ConditionNodeMisuse`
([Appendix C][sC]), carrying the offending argument's type and the node kinds in
hand. It is raised at composition time, before any resolution pass or provenance
chain exists. That is why it is its own kind and not a `ConditionResolution`
sub-kind ([§14.3][s14-3]). The explicit, *ordered* layering spelling —
`override` — belongs with the use case [slot totality](#g-slot-totality)
produces ([§14.6][s14-6]).

### 14.3 Resolution: flatten, validate, compile once

Resolution takes the root node plus a `Build`. Flattening is the only place
path strings are ever concatenated: a trivial recursion with a path
accumulator that also records each entry's **tree position**, its
`getfield`/`getindex` step tuple.

The collecting pass then checks each flat entry:

- the path resolves, with [did-you-mean](#g-did-you-mean) (the offending name
  plus the list-in-hand it should have matched) over children;
- the field is declared in the target's `init_x`/`init_m`;
- the value type is convertible to the declared leaf type;
- [slot](#g-slot) [faces](#g-face) reach root slots;
- no `(path, store, field)` is duplicated.

The `Build` supplies two lookup families. **Schema** is the evaluated
declarations — may you write this field, at what leaf type — and it is the
authority. **Layout** is the destination: `x` backing ranges, store indices for
`m` and for discrete `x`, and slot indices from the
[activation](#g-activation) (a re-run of Stratum C at a given scalar type).
Layout also carries the face chains from [Stratum](#g-stratum) A (one of the
build's three phases: structure, schedule, activation).

A valid list compiles to a plan. Per leaf, the plan holds a `Getter{P}`
[lens](#g-lens) (the compiled navigation step of a condition entry), a
destination offset, and a converter. The lens is the position tuple lifted to a
type parameter, so navigation of the fixed tree type is type-stable.

**Rule.** The converter is baked now, selected per leaf from that leaf's type
in the resolved shape.

**Why.** Leaf types are shape facts, carried by the tree type along with the
full nesting and every field name ([§14.4][s14-4]). Selection therefore
consults no runtime fact and stays a resolution-time bake.

There are two cases:

| leaf in the resolved shape | converter baked |
|---|---|
| already at the activation's scalar type — decision-descended | the type's ordinary `convert`/constructor methods *at that eltype* |
| a plain `Float64` leaf against a non-nominal activation's scratch — a held constant | the `Float64 → Dual` zero-partial embedding |

**A leaf already at the activation's scalar type** is decision-descended: under
a `Dual`-seeded evaluation of a type-stable `trim_condition(d)`, every
decision-dependent leaf is `Dual`-typed in the shape ([§14.7][s14-7]). Such a
leaf takes the type's ordinary `convert`/constructor methods *at that eltype*:
an authored `RQuat` of `Dual`s → the `SVector{4}` state leaf at `Dual`,
partials flowing through untouched. That untouched flow of partials is what
makes the seeded decisions reach the [sweep](#g-sweep) at all. At the nominal
activation the same rule is the ordinary `Float64` conversion: an authored
`RQuat` value → the `SVector{4}` state leaf it initializes.

**A plain `Float64` leaf against a non-nominal activation's scratch** is a held
constant, and takes the `Float64 → Dual` zero-partial embedding. That embedding
is semantically exact in that case, and in no other: "held at the operating
point" *is* zero partials. Zero partials are the whole of a linearization
operating-point [condition](#g-condition), which is authored decision-free
([§14.10][s14-10]).

The selection is a one-time [boundary](#g-boundary) decision, and it leaves the
nominal exact-match doctrine for table [cells](#g-cell) ([§12.5][s12-5])
untouched. Converters run here and in `capture`'s gather ([§14.10][s14-10]) —
the write paths — never on state [views](#g-view) ([§7.1][s7-1]).

Overlay partiality for the `m` and discrete-`x` stores is baked the same way:
the writer holds `merge(init_m_defaults, overlay)` with the base resolved at
compile time (the fork, [§14.1][s14-1]).

### 14.4 Two application registers over one plan

**The paradigm-change tax feared at execution does not materialize.** All string
work, validation and addressing are functions of the *shape* of the
[condition](#g-condition) (the path-addressed sparse overlay that sets a build's
state), and every hot path holds the shape fixed while varying values. Execution
is therefore resolve-once/execute-many, with two [registers](#g-register) over
one plan:

- **Specialized `apply!`** — for services that iterate: trim's per-evaluation
  write, linearization's seeding. It unrolls stores through the baked lenses and
  converters, the same machine operations as today's in-place writes: zero-alloc,
  no strings, no dispatch. The per-iteration shape check is the mechanism
  ([§12.5][s12-5]) transferred. The tree type is proven by dispatch, and it
  carries the full nesting, every field name and leaf type. A `===`
  [sweep](#g-sweep) over the interned path literals closes the remainder. Those
  pointer compares fold to nothing in the all-literal case, and shape drift is a
  structured error, not silent corruption. Cost: Julia codegen of ~10–50 ms
  *once per condition shape*. That cost is noise against the model's own
  first-sweep warmup (seconds), and against the 10³–10⁴ optimizer evaluations
  the codegen amortizes over.
- **Dynamic walk** — for one-shot init. It executes the same validated entry
  list by runtime dispatch per write: microseconds total, allocation permitted,
  since the stopped-sim path was never under the zero-alloc regime
  ([§7.5][s7-5]). It needs no per-shape codegen: fifty structurally different
  scripted conditions cost fifty walks, not fifty compiles.

**Rule.** Which register a service uses is internal, never user-facing API.

#### The read-selector family

**The read-[selector](#g-selector) family is closed.** Its members are
`get_state(path, field[, i])`, `get_deriv(path, field[, i])`,
`get_output(path, field[, i])`, `get_slot(face)` and `get_face(name)` — one
address space for every reader of the model.

The names carry a deliberate `get_` prefix. A selector is a *deferred read*: a
value describing the read the compiled gather will perform. The prefix names
that action, and it keeps five short common nouns out of the namespace user
declarations share with domain code.

There is no selector for a [component](#g-component)'s private intermediates,
and there cannot be one: they are values in flight, not [cells](#g-cell)
anything could address ([§5.2][s5-2]). So a reader that wants one is asking the
producing component to declare it an output ([§11.3][s11-3]). `get_face`
addresses a root-exported output [face](#g-face) — the *integration* register
([§9.2][s9-2]).

**Rule.** A selector resolves against a source, before any client policy
applies.

The table selectors — `get_output`, `get_slot`, `get_face` — resolve against a
*table source*: a [boundary](#g-boundary) [snapshot](#g-snapshot), or the
scratch tables a service evaluation instantiates ([§14.8][s14-8]). The store
selectors — `get_state`, `get_deriv` — resolve only against live stores. The
table/store axis separates table-borne values from store-borne ones, not
snapshots from services.

Only stopped-sim service evaluations, `capture`, and post-run inspection of the
live stores (the [replay](#g-replay)-to-inspect, [§9.2][s9-2]) ever hold live
stores: the snapshot deliberately carries no state stores ([§9.2][s9-2]), and
`ẋ` [buffers](#g-buffer) are integrator scratch, not boundary-consistent
objects outside a service evaluation. A snapshot-bound reader naming a store
selector is therefore a resolution error at attach (`ReadBindingUnresolved`),
raised in the didactic register. The honest remedy ([§9.2][s9-2]) is to declare
the field public and read the [auto-published
port](#g-auto-published-port) (published by the framework from the state or
mode store).

Client policy rides on top — the registers (row 83) restated as a resolver
property:

- **Load-bearing services speak the [contract](#g-contract).** Trim's `reads`
  and linearization's [taps](#g-taps) (the three selector lists declaring what
  linearization seeds and reports) name
  `get_state`/`get_deriv`/`get_output`/`get_slot`/`get_face`, within the scopes
  the locality law ([§6.1][s6-1]) and [fragment](#g-fragment) scoping
  ([§14.2][s14-2]) own. `get_face` is the set's [seam](#g-seam)-crossing member:
  it resolves through export chains exactly as [mounting](#g-mounting)
  (relocating a whole problem or tap set with `at(prefix, …)`) resolves
  [slot](#g-slot) faces ([§14.9][s14-9]), the read side mirroring the write
  side. So an equilibrium equation reaching behind a generically-held child
  binds the curated face register instead of a path the locality law forbids. A
  service evaluation needing a private intermediate has one remedy, and it is
  the same at every register: the component exports it ([§14.7][s14-7]).
- **Diagnostic readers admit the whole family, within the source rule.**
  Output-[device](#g-device) bindings, GUI panels and log inspection take deep
  paths and `get_face` names alike. The store selectors reach only the
  diagnostic clients that actually hold stores (`capture`, post-run
  inspection): a snapshot-bound reader is barred from them by source, not by
  client.
- **`stop_on` is not a family client.** It names root-exported `Bool` output
  faces, period ([§13.5][s13-5], row 60): termination is run policy against the
  root contract, and no path selector reaches `stop_on`.

The five selectors, their sources, and their clients:

| selector | resolves against | load-bearing services | diagnostic readers |
|---|---|---|---|
| `get_state(path, field[, i])` | live stores | named in the contract | only clients that hold stores |
| `get_deriv(path, field[, i])` | live stores | named in the contract | only clients that hold stores |
| `get_output(path, field[, i])` | a table source | named in the contract | admitted |
| `get_slot(face)` | a table source | named in the contract | admitted |
| `get_face(name)` | a table source | named in the contract | admitted |

**Compiled readers are the gather twin** over this family and the layout
tables. Trim's cost read (`ẋ` and output fields), linearization's Jacobian
gather, and `capture`'s full-store readback are one primitive run in reverse:
one machinery, both directions, in the `Build`'s client kit.

The per-iteration ledger for trim is user fragment math (stack-only, the domain
computations unchanged from today) + leaf stores + folded shape check + sweep.
The sweep dominates, exactly as `f_ode!` does today. `apply!` ends at
established stores. Making the model *coherent* is [boundary
zero](#g-boundary-zero) (the initialization boundary: the ordinary
macro-sequence with an empty integrate), [§14.5][s14-5].

### 14.5 Boundary zero: an ordinary boundary with authored incoming transitions

`apply!` establishes the stores at `t₀`. The [trace header](#g-trace-header)
captures those stores, together with the [slot](#g-slot) values, *before anything
below runs* (the capture placement, [§9.5][s9-5]). A post-sequence capture would
hand [replay](#g-replay) already-transitioned state. The init service then
completes the [§8.6][s8-6] macro-sequence with an empty integrate:
project → [[sweep](#g-sweep) → [guards](#g-guard) → handlers]\* →
[due](#g-due) `g` updates → first [snapshot](#g-snapshot). The parity with an
ordinary boundary is exact, not approximate. Piece by piece:

- **Project runs.** Authored `x` can sit off-manifold: a hand-assembled
  quaternion ulps off unit norm, or a [condition](#g-condition) (the
  path-addressed sparse overlay that sets a build's state) writing part of a
  constrained block against fresh defaults. [Projection](#g-projection) after
  condition writes is the same position it holds after any other `x` mutation.
  And it costs nothing when the state is already clean.
- **The sweep runs with every `Φ = 0` [tick](#g-tick) due.** `t₀` is a grid point
  of every phase-free divisor, so those discrete output stages are gated in and
  publish from the authored `x`. They can do nothing else: no earlier tick exists
  for a ZOH to hold. An offset [component](#g-component) ([§8.5][s8-5]) is *not*
  due — its first tick is at `Φ·Δt_base`. Until then that component's
  [cells](#g-cell) hold the values the [activation](#g-activation)'s
  [probe](#g-probe) populated ([§12.3][s12-3]). The `t₀` snapshot carries a fully
  populated table either way.
- **Events run.** A condition can land a guard [predicate](#g-predicate) in
  [holding](#g-edge-semantics) territory (firing on not-holding → holding
  transitions, never bare sign changes): an authored stall flag, a strut authored
  into contact. The event then fires visibly at `t₀` rather than one step later.
  The [prior](#g-prior) rule grounds that timing ([§8.6][s8-6]):
  [boundary zero](#g-boundary-zero) establishes every guard prior — the event's
  stored predicate sample from the previous boundary — as not-holding.
  Suppressing those firings was rejected (row 67), on the
  [stage-on-interaction](#g-stage-on-interaction) lesson of [§9.7][s9-7] (widgets
  stage on edit or activation, never per render pass): insurance that masks
  invariant violations is anti-diagnostic. The header records the *resolved
  pre-sequence* stores and slots ([§9.5][s9-5]), so replay re-executes
  [boundary](#g-boundary) zero from the same starting point and whatever fires at
  `t₀` fires again identically ([§10.7][s10-7]). The firings are recomputed,
  never recorded. (A `stop_on` [face](#g-face) already `true` is a different
  category: nothing *fires* — the face simply reads `true` in the published `t₀`
  snapshot and the loop reacts, [§13.5][s13-5].)
- **Due `g` updates run.** This follows from an interval-alignment fact that
  is easy to mis-picture and is hereby a taught [contract](#g-contract), sibling
  to the boundary-sampling line ([§15.5][s15-5]): **a boundary's `g` is the
  *outgoing* transition** — at tick `t_k` it consumes the completed boundary's
  samples and produces `x_{k+1}`, the value the next tick reads. The transition
  that carried `x` *into* `t_k` ran at `t_{k-1}`. Boundary zero is missing its
  incoming transitions on *both* tiers, and both are replaced by authorship:

  | [tier](#g-tier) | `t₋₁` | `t₀` (boundary zero) | `t₁` |
  |---|---|---|---|
  | discrete | the `g` that would have produced a discrete leaf's `x(0)` never ran; the condition authored `x(0)` | `g` consumes the `t₀` samples and produces `x(1)` | the gated stages read `x(1)` |
  | continuous | the integration over `[t_{-1}, t_0]` that would have produced a continuous leaf's `x(0)` never ran; the condition authored `x(0)` | the authored `x(0)` is the initial condition of the outgoing integrate, $t_0 \to t_0 + h$ | |

  The outgoing work all runs, and `t₀`'s `g` has its only opportunity: `x(1)`
  must sit in the store before `t₁`'s gated stages read it. An accumulator
  $x_{k+1} = x_k + \Delta t \, e_k$ authored with $x_0 = 0$ under nonzero
  $e(t_0)$ would otherwise first integrate $e(t_1)$, putting the whole
  sampled-data lattice one period late (row 67). The authored `x(0)` needs no
  protection: it is published in the `t₀` snapshot regardless. The
  continuous-tier analogue of `g`-at-`t₀` is not the empty incoming integrate but
  that first *outgoing* one, and both authored values are the published initial
  conditions of their outgoing transitions.
- **`t₀` is an init-service argument** (default `0.0`), never a condition
  entry: time is not a store of any component. The
  [harmonic grid](#g-harmonic-grid) (every discrete period an integer multiple of
  `Δt_base`) anchors at whatever `t₀` boundary zero runs at. Both init-service
  entry points carry the argument, with the same default:
  `init!(sim, condition; t0)` and `trim!`'s commit
  (`trim!(sim, problem; baseline, t0, backend)`, [§14.8][s14-8]). Conditions are
  time-free, and `capture` returns condition and time separately for
  resume-at-time — the returned `t` passed back as `t0`.
- **Trim is untouched by all of this.** Optimizer iterations are raw
  write → sweep → read cycles on the activation — no boundaries, no events,
  no `g`. Only the committed solution executes boundary zero.
- **A guard firing at commit is a wanted failure signal.** Today's hand-written
  trim asserts (`!stall`, no weight-on-wheels, `ω > ω_idle`) become the model's
  own event logic, surfaced through the ordinary machinery instead of `@assert`.
  And the signal has a channel: boundary zero reports the fired set on the
  `TrimReport` and raises `TrimCommitEvents` when that set is non-empty
  ([§14.8][s14-8], [Appendix C][sC]). A handler that fires at commit moves the
  committed stores off the solved point, and saying nothing would be
  warn-but-assign relocated.
- **A commit-fired handler is not the only mover, and the second one is
  unconditional.** Boundary zero's *first* act is `project`, so the committed
  `x` is `project(x*)`, not the solver's `x*` — an attitude quaternion
  renormalized by a few ulps is the canonical case. That move is legitimate,
  wanted, and usually invisible in the residuals; but the point the stores sit
  at is no longer the point the verdict was read at.
- **The doctrine is the same for both movers, and so is the remedy.** The
  `TrimReport` carries the committed-state residuals beside the solved-point
  ones, and a converged solve whose committed-state residuals fail the box test
  raises `TrimCommitResiduals` ([§14.8][s14-8]) — the move made visible rather
  than silent.

### 14.6 Slot totality: the missing-value error and the `override` combinator

[Slots](#g-slot) are the one initialized datum without declared defaults — the
bare-types decision ([§9.3][s9-3]), upheld here (row 68). So a
slot's only source before [boundary zero](#g-boundary-zero) (the initialization boundary: the
ordinary macro-sequence with an empty integrate) is the condition, and three
consequences follow.

**Totality is a precondition of starting, checked by the service.** A
condition value is legitimately partial ([fragments](#g-fragment) compose; trim iterations
write subsets; capture-then-tweak patches leaves) — "every root slot
covered" is not a property of conditions but of *every application that
establishes a complete world over virgin stores*. That principle, not an
enumeration, names the sites: `init!`, trim's setup application to freshly
allocated scratch stores, and trim's commit through [boundary](#g-boundary) zero
([§14.8][s14-8]) — one class, one mechanism, one kind. Each compares the resolved
plan's slot coverage against the `Build`'s `input_faces` before writing
anything; a shortfall is one collected, declaration-ordered diagnostic
(`UninitializedSlots`, a [§13.2][s13-2] kind) naming every uncovered face.
Coverage is a *plan-level fact* — both operands are resolution-time data —
so the check is one comparison and runs before any evaluation, not merely
before any write. Pre-write means all-or-nothing: a rejected init leaves
the sim exactly as it was, the same posture as failed trim.

**The [probe](#g-probe)-value barrier is structural.** The `probe_value`
synthesis ([§12.3][s12-3]; zero/false/first-enum/`T()`) exists so build-time
probes can exercise code
with fabricated inputs; a fabricated zero is a fine probe input and a
terrible flight condition (a silently zeroed `mixture` kills the engine and
sends the user debugging aerodynamics). The services path simply contains
no call to it: a slot gets a condition value or the application errors —
no third branch. [Replay](#g-replay) likewise never synthesizes: the [trace header](#g-trace-header)
records every slot value, and with totality enforced its slot capture is
complete by construction (the requirement discharged, [§9.3][s9-3]).

**[Baselines](#g-baseline) are aircraft-shipped [condition](#g-condition) functions, layered by
`override`.** Nobody hand-writes ~20 slot values per script; today
`SystemsInitializer`'s `@kwdef` defaults carry that load, and their
successor is ordinary user math in one authoritative home —
`ready_for_taxi(ac)`, `cold_and_dark(ac)` — returning full-coverage
conditions. But "baseline plus tweaks" collides with the duplicate-leaf error
([§14.2][s14-2]) by design: the collision *is* the intent. Hence the
fourth node kind, **`override(base, patch)`** — ordered and asymmetric
where `merge` is symmetric and collision-intolerant. At resolution a leaf
present in both takes the patch's value, with provenance recording both
sources ("patch overrode base's `throttle`"); collisions *within* one
layer remain errors; variadic layering
(`override(campaign, aircraft, todays_case)`) composes. Trim uses it on
day one: the committed condition is `override(baseline, solution)` — the
solver's handful of values over full coverage (row 68).

### 14.7 The trim problem: NamedTuple decisions, declared reads, named residuals

The aircraft author ships one value: what the solver may vary, what those
decisions make of the model, what to read back after each evaluation, and which
equations the readings must satisfy.

**Rule.** The field set is normative. The lift ([§14.9][s14-9]) is
field-by-field, so this list is closed:

- `guess`, `lower`, `upper` — same-named all-`Float64` NamedTuples;
- `condition` — the condition-valued function over decisions;
- `reads` — the declared read set;
- `residuals` — the residual function;
- `tolerances` — an all-`Float64` NamedTuple, same-named as the residual
  function's return.

`tolerances` is carried *in the problem* because a relocated problem must carry
its own convergence test; `at` passes it through untouched.

[Worked](#g-worked), the C172 cruise case reduces to its three-equation core.
The real problem is the same shape with the full 7-variable search, and the
θ-constraint elimination survives inside `trim_condition`:

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

The rest of the section takes what the author ships one piece at a time,
against today's `c172.jl`:

- **Decision variables, initial guess and box bounds are plain, same-*named*,
  all-`Float64` NamedTuples.** The `AbstractTrimState{N}`/`FieldVector`
  supertype dies; its only job was vectorization, and vectorization is the
  service's. The service packs and unpacks by field order, that order being the
  `guess` NamedTuple's own. [§12.5][s12-5] states the rule this
  [seam](#g-seam) runs on: **the names are the pairing, order carries no
  semantics**. `lower` and `upper` are checked at setup for key-set equality
  with `guess` and `Float64` fields, then canonicalized to `guess`'s field
  order by the same type-level reorder (`NamedTuple{keys(guess)}(lower)`). A
  permuted bound spelling is therefore a non-event rather than `α`'s bound
  silently applied to `throttle`; a key-set or field-type mismatch is
  `TrimProblemInvalid`. Guess, bounds and the returned solution share one
  spelling, and `Base.merge(guess, (throttle = 0.3,))` is free warm-start
  tweaking. An author who wants a documented `@kwdef` struct keeps it privately
  and converts.
- **`TrimParameters` stays a plain user struct** the framework never sees;
  the assignment is the pure `trim_condition(ac, params, d)` [fragment](#g-fragment)-tree
  function ([§14.2][s14-2]), applied per iteration by the compiled plan ([§14.4][s14-4]).
- **The read side is declared, then compiled**: `reads(name = get_state(path,
  field) | get_deriv(path, field) | get_output(path, field) | get_slot([face](#g-face)) |
  get_face(name), ...)` — the load-bearing set ([§14.4][s14-4]). `get_state` and
  `get_deriv` address a declared state field and its derivative (validated
  against `init_x`), `get_output` a declared output [port](#g-port) (validated against
  `output_types`), `get_slot` and `get_face` a root input and output face
  (validated against the root face lists). The path [selectors](#g-selector) (the closed
  family of deferred reads resolving against a source) reach only through the
  locality scopes ([§6.1][s6-1]). An equilibrium equation crossing a generic seam reads
  a face. A [component](#g-component)'s private intermediates are not addressable at all ([§5.2][s5-2],
  [§14.4][s14-4]); a trim evaluation needing one is a signal the component should
  export it. A derivative wanted across a [contract](#g-contract) [boundary](#g-boundary) takes the same
  remedy: publish it as an ordinary output port computed in `h_xu` ([§7.4][s7-4] step
  2's one-line binding, made contract). That leaves `get_deriv` scoped to owned
  concrete subtrees. The compiled reader (the gather twin, [§14.4][s14-4]) fills a
  stack-only NamedTuple per evaluation.
- **The user supplies a residual *system*, not a scalar cost.** It is authored
  as a NamedTuple of named equations and packed to the solver's vector by field
  order — `tolerances`' field order here, the declared side again, with the
  residual return canonicalized to it. The decisions rule holds symmetrically
  on both ends of the seam: names pair, order never does.
- **The FlightCore formulation's core is correct and survives verbatim as user
  math.** That core is analytic elimination: `θ_constraint` substituting the
  pitch constraint, filter and actuator equilibria imposed by construction, and
  the minimal 7-variable search.
- **What changes is the numerics.** Trim is a square root-find, and
  FlightCore's derivative-free scalar minimization over $\|r\|^2$ against a
  hand-scaled absolute `stopval` was the rational choice only because Jacobians
  through the mutating `f_ode!` chain and the assignment math were out of reach
  (row 69).
- **Nonlinear least squares with exact AD Jacobians is the default.** The
  `Dual` [activation](#g-activation) seeds the decision variables through the
  `T`-generic assignment, [sweep](#g-sweep) and `f`. The seeds survive the condition
  write boundary because [§14.3][s14-3] selects the baked converter per leaf from the
  shape: a decision-descended leaf is `Dual`-typed there and takes the
  structural conversion, while the zero-partial embedding stays on the held
  `Float64` leaves. What [§12.6][s12-6] leaves open as an option is here the
  *default*: nonlinear least squares on $r(d)$ with exact AD Jacobians, in the
  trust-region/Levenberg–Marquardt register. Convergence is quadratic (~5–15
  evaluations), the tolerances are per-residual and physical, and failure
  reports name the unbalanced equations with magnitudes. The convergence
  verdict itself is service-owned and backend-independent ([§14.8][s14-8]).
- **Non-squareness degrades gracefully.** Redundant actuation becomes weighted
  or minimum-norm least squares; infeasible demands converge to a nonzero
  residual identifying the impossible balance. At the solution,
  $\partial r / \partial d$ is free flight-physics data (control effectiveness)
  cross-checking linearization.
- **The derivative-free scalar path survives as the fallback.** The service
  squares *and normalizes* the residuals against the tolerances —
  $\sum (r_i/\mathit{tol}_i)^2$ at `stopval = 1` ([§14.8][s14-8]). That
  normalization is where the hand-scaled absolute threshold is repaired,
  leaving today's algorithm as the degenerate case.
- **[Recorded, not built](#g-recorded-not-built)** (a worked-out extension
  deliberately left unimplemented, its seams named): closed-loop sampled-data
  trim and on-ground static equilibrium, each simply another problem value over
  the same service. Closed-loop trim appends $g(x) - x = 0$ residuals via a
  nondestructive scratch evaluation of `g`, structurally impossible under
  FlightCore's mutating `f_disc!`. On-ground static equilibrium is the other:
  strut compressions and attitude against gear forces.

The residual signature:

```julia
residuals(reads::NamedTuple, d::NamedTuple) → NamedTuple
```

**Rule.** What the solver varies is passed; what is fixed per problem is closed
over.

The gathered reads and the decision NamedTuple arrive as arguments; `d` is the
one value that *cannot* be closed over. `TrimParameters` stays behind the
closure, exactly as `condition` already holds it (the framework never sees it).
Being user-shaped, that record is also where any environment handles the
condition math needs conventionally ride (the [value-level constructor](#g-value-level-constructor),
[§4.4][s4-4]; the pre-sweep doctrine, [§14.1][s14-1]). The problem *receives* the
environment, and never writes it ([§14.9][s14-9]).

The returned NamedTuple's names are the equation names the report and the
failure messages use. The service packs residuals and tolerances by field
order: `tolerances`' order is canonical for the r-side, and each residual
return is reordered to it (`NamedTuple{keys(tolerances)}(r)`) before packing.
An equation list spelled in a different order inside the lambda therefore pairs
correctly and costs nothing. Names and types are checked at setup: the guess
evaluation the service performs anyway observes the residual key set. A
`tolerances` key-set mismatch — or any field-type disagreement, on either side
of the seam — is `TrimProblemInvalid` ([Appendix C][sC]), with the offending
field and the names or types in hand. Order is never a mismatch.

### 14.8 The trim service: solver seam, scratch stores, commit and report

A `TrimProblem` is an inert value until a service runs it. `trim!` is that
service: it drives a solver it holds behind one method, works on scratch stores
of its own, commits a converged solution the way `init!` would, and returns a
report.

#### The backend seam

`trim!(sim, problem; baseline, t0 = 0.0, backend = LevenbergMarquardt())`. The
default is an in-house dense Levenberg–Marquardt. For decision dimensions ~10
with exact Jacobians, the core is ~100 lines: a damping loop, a small linear
solve, a convergence test. That is the [§8.2][s8-2] stepper precedent exactly
(tiny needed core vs. heavy dependency). The per-residual physical tolerances
sharpen the case ([§14.7][s14-7]): they are a convergence test no external
package spells natively. That is precisely why the *service*, not the backend,
applies it.

The backend [contract](#g-contract) is a **pinned signature**, value-passed —
one required method per backend. That one method is the [seam](#g-seam):

```julia
solve(backend, eval!, d0, lower, upper, tol) -> (; d, status, nevals, niters)
```

`eval!(r, J, d)` is in-place and always fills `r`, the residual vector packed
in `tolerances`' field order. It fills `J` **iff `J !== nothing`**. The
Jacobian is requested by argument, so a Jacobian-free backend simply always
passes `nothing`, and there is still exactly one evaluation method to
implement.

`d0`, `lower` and `upper` are packed `Vector{Float64}` in `guess`'s field
order, with `±Inf` meaning unbounded. The declared side is canonical on both
ends, and the service has canonicalized before packing ([§14.7][s14-7]). A
backend that ignores bounds therefore ignores two vectors, not a missing
argument.

`tol` is a `Vector{Float64}` in `tolerances`' field order. It is data the
backend *may* stop on — that is the service's per-register translation,
below — and decisive of nothing.

Returned: `d`, the solution; `status::Symbol`, from a deliberately **open**
set; and `nevals`/`niters`, diagnostic counts. The status is recorded verbatim
in the report because the verdict is the service's (row 158). The name `solve`
is subject to the [§16][s16] naming audit like every other API spelling. The
backend sees vectors and never names, so the solution it returns unpacks by the
same order it was given.

#### The convergence verdict is the service's, uniformly

`converged` means `all(abs.(rᵢ) .≤ tolᵢ)`: the per-residual box test
([§14.7][s14-7]) in its own physical units, evaluated **by the service at the
backend's returned point**. That is one residual evaluation, noise against the
solve that produced it. That verdict, and nothing else, gates the commit and
fills `TrimReport.converged`. The backend's returned `status` is recorded in
the report as diagnostic data and is authoritative over nothing.

#### The tolerance translation, per register

The tolerance translation is the service's too, per register. In the
least-squares register the tolerances *are* the stopping criterion: they feed
the per-residual test directly, LM's damping loop testing exactly what the
service will re-test.

The derivative-free scalar fallback is `NLoptBackend(:LN_BOBYQA)` in a package
extension: it passes `nothing` for the Jacobian, keeps today's algorithm one
keyword away, and leaves the framework core carrying zero optimizer
dependencies. For that fallback the service squares *and normalizes*: the
objective minimized is $\sum_i (r_i/\mathit{tol}_i)^2$ with `stopval = 1`. That
objective is dimensionless where FlightCore's threshold was hand-scaled and
absolute ([§14.7][s14-7]), and a well-scaled valley where a raw $\|r\|^2$ sums
forces against moments.

The two criteria cannot disagree in the dangerous direction:

$$\sum_i (r_i/\mathit{tol}_i)^2 \le 1
\quad\Longrightarrow\quad (r_i/\mathit{tol}_i)^2 \le 1 \;\text{ for every } i$$

so the `stopval` sphere is *inscribed* in the tolerance box, and a fallback
stopping at `stopval` necessarily passes the service's box test. The converse
disagreement — a backend stopping early and reporting an optimistic status —
is caught by the re-check, which remains the single authority. What is *not*
claimed is point identity: different backends may land on different solutions,
an algorithmic difference and a legitimate one. What is eliminated is
per-backend meanings of `converged`.

#### Box bounds and saturated decisions

Box bounds are honored by step projection. A decision variable saturated *at
the solution* is flagged in the report: "converged with `elevator` at its upper
bound". That is the classic CG-limit diagnostic, today inferable only from
mysterious residuals.

#### Scratch stores, stated without type luck

Every `trim!` invocation instantiates a fresh working [store](#g-store) set:
`x` backing, `m` and discrete-`x` stores, [slot](#g-slot) and
[signal tables](#g-signal-table), derivative [buffer](#g-buffer). The set is
built from the [activation](#g-activation)'s *layout* (a re-run of Stratum C at
a given scalar type). The layout is the reusable compiled artifact; the buffers
are per-invocation and die with the call (stopped-sim allocation,
[§7.5][s7-5]).

The `Dual` backend's buffers being un-aliasable by type is defense in depth,
not the mechanism: a `Float64` backend (NLopt) gets equally fresh buffers. The
invariant is backend-independent: **the simulation's authoritative stores have
exactly one writer, the commit through [boundary zero](#g-boundary-zero)** (the
initialization boundary: the ordinary macro-sequence with an empty integrate).

Setup applies `override(baseline, condition(guess))` to the scratch set once,
and that application's full coverage is *checked here* — the comparison of the
resolved plan against the `Build`'s `input_faces` ([§14.6][s14-6]), one
plan-level comparison before the first evaluation. [Sweeps](#g-sweep) therefore
see a complete world. Raw instantiation is sound exactly because of that check:
every slot is written before any read. An incomplete `baseline` is one
declaration-ordered `UninitializedSlots` at setup rather than a whole solve
against undefined [cells](#g-cell).

The commit applies the same composite over the same `baseline`
(`override(baseline, condition(d*))`, [§14.9][s14-9]), so its coverage is
setup's. Commit's totality check is therefore structurally unfailable through
the trim path, and it stands as the shared `init!`-[boundary](#g-boundary)
defense. A converged solve is always committable, so `TrimReport` carries no
committed flag and the no-throw doctrine needs no exception.

Iterations rewrite only the problem's write-set via the compiled plan; an LM
evaluation is one Dual-seeded sweep yielding `r` (value parts) and `J`
(partials) together. No convergence — the service's box test failing at the
returned point, whatever status the backend attached to it — means no commit.
No commit means the sim is bit-for-bit untouched, including "never
initialized". Today's warn-but-assign is structurally impossible.

The same structure covers an interrupt. Ctrl-C during a long solve unwinds an
ordinary Julia call operating on per-invocation scratch stores: no commit has
happened, and the simulation is bit-for-bit untouched exactly as a
non-converged solve leaves it. The services need no counterpart to the loop's
boundary masking ([§10.4][s10-4]).

#### The commit, in full

The committed solution is applied as an `init!` in every respect:
`override(baseline, solution)` through boundary zero, with the pre-write
[slot totality](#g-slot-totality) check ([§14.6][s14-6]), the sequence
([§14.5][s14-5]) and [guards](#g-guard) at commit.

The [harmonic grid](#g-harmonic-grid) (every discrete period an integer
multiple of `Δt_base`) is anchored at `trim!`'s own `t0` argument, default
`0.0`. That is `init!`'s default too: one rule for both init-service entry
points. The [recorders](#g-recorders) are cleared exactly as [§10.6][s10-6]
states for `init!`: the [trace](#g-trace), the log, and any
[batches](#g-batch) still in [staging cells](#g-staging-cell) (where a device's
pending write batch waits between drains).

A fresh recording starting at its own anchor is the unattended register's
natural shape. Fly-then-retrim keeps continuity explicitly: `capture` returns
`(condition, t)` separately for exactly this ([§14.1][s14-1]), so the resumed
spelling is `trim!(sim, problem; baseline = c, t0 = t)`.

#### The report, not an exception

`trim!` returns a structured `TrimReport`:

```julia
# what trim! returns: a value to read, never an exception to catch
struct TrimReport
    converged   # the service's box test at the backend's returned point
    …           # the remaining fields, enumerated below; this section
                # does not name them
end
```

Field by field:

- the `converged` flag — the service's own box test at the returned point,
  never a backend's opinion;
- the solution NamedTuple — guess-shaped, hence warm-startable;
- the **solved-point residuals** with their tolerances — the very numbers the
  verdict is read off, gathered at the backend's returned point;
- the **committed-state residuals** — the same residuals re-gathered from the
  boundary-zero world after the commit;
- the backend's returned status together with its iteration/evaluation
  counts — diagnostic throughout: informative about *how* the solve went,
  decisive about nothing;
- the saturated-bounds list;
- the commit's fired events — component paths and event names, empty when
  boundary zero ran quiet ([§14.5][s14-5]).

The committed-state residuals are nearly free: that boundary's sweep has
already run, so the residuals' declared reads need only gather from it. The
offset caveat carries over from [§14.5][s14-5]: an offset
[component](#g-component) is not [due](#g-due) at boundary zero, so a residual
reading its [port](#g-port) reads the [probe](#g-probe)-populated cell, not a
commit-refreshed one. Those committed-state residuals are the numbers
describing the state the simulation is actually *in*, which is the point a
`capture`-defaulted `linearize` reads. A non-empty fired-event set also raises
`TrimCommitEvents` ([Appendix C][sC]): the committed stores then sit at the
post-handler point, not the reported solution, and a `capture`-defaulted
`linearize` ([§14.10][s14-10]) reads that point.

The two residual sets are what make the moved point auditable. A converged
solve whose *committed-state* residuals violate the box test raises
`TrimCommitResiduals` ([Appendix C][sC]), naming the offending residuals with
their committed values and tolerances. The move — `project`, or a commit-fired
handler ([§14.5][s14-5]) — is surfaced rather than left silent. The verdict
itself is not re-litigated: it gated the commit, at the solved point, and the
numbers (row 150) stand as reported.

Non-convergence never throws: it is an expected *outcome* (envelope-sweep data:
hitting the infeasible edge is information), per the
exceptions-are-broken-machinery line ([§13][s13]). A malformed problem is a
different case: a `BuildError`-class failure at setup, `TrimProblemInvalid`
([Appendix C][sC]). The malformed cases: a guess/bounds key-set or
field-type disagreement, an unknown `reads` [selector](#g-selector), a
`tolerances`/residual key-set mismatch observed at the setup guess evaluation.
The error carries the offending field with the names or types in hand,
collected, mirroring linearization's `TapResolution`. A permuted spelling is
none of these ([§14.7][s14-7]).

An *empty* problem is none of them either:
[`TrimProblem`](#g-trimproblem)`(guess = (;), …)` is legal, not
`TrimProblemInvalid`. With zero decision variables the solver is bypassed
outright — nothing to pack, no seeded activation, no backend call. The service
simply evaluates the residuals once at the [baseline](#g-baseline), the
ordinary box test deciding `converged` and the commit running as usual. The
degenerate problem is the "is this operating point an equilibrium?" probe:
evaluate this condition's equations and report, useful in its own right and
free.

#### The AD obligation, scoped

The default formulation requires `Dual` genericity of exactly the continuous
output-stage chains and `f`, plus the user's assignment and residual math. The
discrete [tier](#g-tier)'s stages and `g`, and the event system's guards and
handlers, never see a `Dual`: they are frozen constants with zero partials,
semantically exact ([§11.2][s11-2]).

This is *not a new obligation*. It is the same activation linearization is
defined on. AD-readiness is also a build-checked property: the Dual probe
detonates [pinned](#g-walked) intermediates with a culprit-naming
`InexactError`, and `build(world; activations)` puts it in CI. The robustness
comes from enforcement, not hope.

C172 migration audit (one afternoon, one genuine item):

- `Interpolations.jl` tables (propeller coefficient maps, engine maps) must
  accept generic scalars. They do; but prefer cubic knots over linear where
  partials matter, since linear knots make Jacobian entries piecewise-constant.
- In-model saturations (actuator limits, idle/FRC clamps) zero Jacobian columns
  when active. LM damping tolerates the rank deficiency and the report names
  it, and cruise trim leaves those saturations inactive.
- The landing gear is never evaluated off-zero airborne.
- `norm`-at-zero guards are already in place (e5efb3a).

Fallback per problem: one `backend =` keyword.

### 14.9 Mounting: problems as relocatable values

**What a `TrimProblem` is.** Not a condition, but an **implicitly specified
condition**: `condition` is a condition-*valued function* over the decision
space, `reads`/`residuals` are the equations that pin the free variables
down, `guess`/`bounds` say where to search. Solving makes the implicit
condition explicit. The commit is then literally an init:
`override(baseline, condition(d*))` through [boundary zero](#g-boundary-zero) (the initialization
boundary: the ordinary macro-sequence with an empty integrate). The services
unify as clients of one condition algebra: `init!` applies an explicit
condition, `capture` produces one, `trim!` searches a family for the member
satisfying its equations.

**`at` lifts to problems in five lines.** Every field of a problem is
either condition-producing (path-relative) or path-free. The rule that residual
math sees only the gathered NamedTuple ([§14.7][s14-7]) pays off here:

```julia
at(prefix::String, p::TrimProblem) = TrimProblem(
    guess      = p.guess,                      #path-free: pass through
    lower      = p.lower,
    upper      = p.upper,
    condition  = d -> at(prefix, p.condition(d)),  #post-compose: wrap each returned tree
    reads      = at(prefix, p.reads),              #inert selector data: same Scoped node
    residuals  = p.residuals,                      #path-free: pass through
    tolerances = p.tolerances)
```

Resolution then needs nothing new. The flattening accumulator of [§14.3][s14-3]
enters the `Scoped` wrapper and prefixes every entry
(`"vehicle/dynamics"` → `"wing/vehicle/dynamics"`). [Slot](#g-slot) entries authored
in the aircraft's [face](#g-face) vocabulary resolve through the export chain *from
the mount point* (`throttle` at `"wing"` → root slot `"wing.throttle"`).
An unexported face fails resolution by name, and correctly so. An internally
wired input (a [scenario component](#g-scenario-component) driving the wingman's throttle) is
untrimmable from outside, and the build says so. The service compiles the
scoped condition and reads and runs the identical loop — it never knows
where its paths are mounted.

**A problem never authors the environment.** The environment is the world's
and the `baseline`'s business: a sibling [component](#g-component)'s slots in a full world,
a handle-valued root slot in a thin rig. The problem *receives* its handles
through the user parameter record ([§14.7][s14-7]), and only queries them.

**Why.** The reason is the resolution rule just stated: a condition entry
naming a wired input fails by name, correctly. A problem writing an
environment face would therefore be applicable only to those rigs where that
face happens to be unconnected, and the relocatability this section exists to
guarantee would be lost.

**The world wrapper dissolves.** Today's `f_init!(::Model{<:SimpleWorld})`
(initialize environment, then call the aircraft's trim) has no successor
method: the environment, the other aircraft and all slots are covered by
the `baseline` condition ([§14.6][s14-6]), applied once at setup. The commit is
`override(baseline, at(mount, condition(d*)))`. Method nesting became value
layering.

**"Aircraft as root" is a thin world.** By default the aircraft is not
literally the root: its environment inputs ([§4.4][s4-4] function-valued
signals) are wired from provider components. Design tasks therefore use a
shipped rig, `design_world(ac)` = aircraft + `SimpleAtmosphere(wind = NoWind())` +
`HorizontalTerrain`. That rig is today's ad-hoc models inside `linearize`
promoted to a named artifact. One register: the "root" case is the shallowest
world, the trim problem mounts at `"aircraft"` like anywhere else. Leaving an
environment face *unconnected* is legal by construction, though. The face
becomes an ordinary root slot holding the handle **value**, written by the
`baseline` like any other slot. That is the test-rig register: the
function-valued sibling of a [constant source](#g-constant-source) (a library component publishing
a value its instance holds), zero ceremony for a frozen environment. For
design tasks the shipped rig stays `design_world(ac)`. That keeps the
environment's tunables in the slot vocabulary that conditions, `capture`,
linearization's input surface and the [trace header](#g-trace-header) already speak.

**Swarm doctrine.** The service solves *one problem at a time*. Sequential
independent trims (trim lead, commit, trim wing against the committed
world) cover weak/one-way coupling. A joint trim is user-side value
composition: concatenate decision NamedTuples under prefixed names, merge
the scoped condition trees, stack the residuals. If joint trims become
routine, a `product(p₁ => "lead", p₂ => "wing")` helper belongs in the
[§13.7][s13-7] library. That helper is [recorded, not built](#g-recorded-not-built) (a worked-out
extension deliberately left unimplemented, its seams named).

### 14.10 Linearization: tap selectors, one seeded pass, a pure query

**The tap declaration.** Today's per-aircraft `XStateSpace`/`UStateSpace`/
`YStateSpace` structs, plus the `get_*_ss`/`assign_*_ss!` shuttle methods, run
to ~150 lines of bookkeeping per variant. All of it becomes three
[selector](#g-selector) lists (the closed family of deferred reads resolving
against a source). Three members of the read-selector family supply them:
`get_state`, `get_slot` and `get_output` ([§14.4][s14-4]). The lists carry the
optional [component](#g-component) index, so a vector leaf yields *named
scalars*. The NamedTuple key is the label control design slices by:

```julia
x = (p = get_state("vehicle/dynamics", :ω_eb_b, 1),
     θ = get_state("vehicle/kinematics", :θ_nb), …)
u = (throttle_cmd = get_slot("throttle"), …)
y = (EAS = get_output("vehicle/airflow", :EAS), …)
```

The three lists are validated at resolution against
`init_x`/[faces](#g-face)/`output_types`, with
[did-you-mean](#g-did-you-mean) errors (the offending name plus the
list-in-hand it should have matched). They compile to offsets once, and
relocate whole via `at(prefix, taps)`. The shuttle layer's successor is that
compiled writer/reader pair, and the promised `get_x_ss` deletion
([§7.1][s7-1]) is discharged.

**The evaluation.** Each invocation instantiates its own scratch store set —
the trim service's mechanism verbatim ([§14.8][s14-8]) — and applies the
operating-point condition. It then runs **one** Dual evaluation, seeded with
one direction per `x`-tap and per `u`-tap entry (chunked internally). Value
parts give `ẋ₀` and `y₀`; partials give `A`, `B`, `C`, `D` simultaneously,
exact to machine precision:

```
  x-taps ─┐                                         ┌─ value parts → ẋ₀, y₀
          ├─ seed directions → one Dual evaluation ─┤
  u-taps ─┘          (chunked internally)           └─ partials    → A, B, C, D
```

That single pass replaces four `FiniteDiff` jacobians, their step-size
heuristics and ~4n perturbed evaluations.

**What the pass seeds, and what it holds frozen.** Unseeded states sit constant
at the operating point, and so do unseeded [slots](#g-slot): the condition
apply embeds their `Float64` values as zero-partial constants. A root-slot
[cell](#g-cell) follows the [activation](#g-activation) scalar (a re-run of
Stratum C at a given scalar type) by *evaluating* its consuming `input_types`
entry at that scalar ([§11.2][s11-2]). The discrete [tier](#g-tier) is frozen
with zero partials — precisely "linearize with the discrete state held"
([§11.2][s11-2]). Differentiation participation is a per-invocation
*seeding* fact for every slot the schema leaves seedable, one register for `x`
and slots alike. One declared exception is visible in the schema: a slot whose
entry is declared `Float64` is **declaredly unseedable**, its cell frozen at
every activation. Selecting it as a `B`-matrix tap is therefore rejected at tap
resolution with the offending entry in hand, rather than silently yielding a
zero column (row 167). Under fan-out the rejection names the **pinning
consumer**, not the face alone: a slot is unseedable whenever any one of its
consumers demands frozen, which is the fan-out meet ([§11.2][s11-2]; row 168).
The author's next move — promote that leaf to a tolerant entry, or route the
tap around it — depends on knowing which leaf froze the slot. Seeded and
frozen, side by side here:

| leaf | in the one pass | what fixes it |
|---|---|---|
| a state leaf named by an `x` tap | seeded, one direction | per-invocation seeding |
| any other state leaf | sits constant at the operating point | per-invocation seeding |
| a root slot named by a `u` tap | seeded, one direction | per-invocation seeding |
| any other seedable root slot | constant, its `Float64` value embedded as a zero-partial constant | per-invocation seeding |
| a root slot whose entry is declared `Float64` | frozen at every activation, and rejected as a `B`-matrix tap | the schema: a declaration |
| any discrete-tier leaf | frozen, zero partials | the tier |

**A pure query, and the shape of `capture`.** Linearization is the first
service with no commit and no [boundary zero](#g-boundary-zero) (the
initialization boundary: the ordinary macro-sequence with an empty integrate).
It works on scratch [buffers](#g-buffer) only, and nothing it computes becomes
authoritative. Today's restore-the-trim dance — re-`assign!` after
`FiniteDiff` dirtied the model — has no successor. The default operating point
is the sim's current committed state, taken through
`capture(sim) → (condition, t)`. That gather covers stores *and root slots* in
full, slot totality ([§14.6][s14-6]) making slot coverage mandatory for
capture → apply. After a `trim!` commit, `linearize(sim, taps)` is about the
trim point with nothing re-specified. An `about = <condition>` keyword
linearizes anywhere else without touching the sim.

**The returned object and `LinearizedSS`.** `linearize` returns labeled data:
`(ẋ₀, x₀, u₀, y₀, A, B, C, D)`, carrying the label sets of the
[taps](#g-taps) (the three selector lists declaring what linearization seeds
and reports). On that data, `subsystem`/`delete_vars` survive as pure
label-indexed matrix slicing, with no model involvement. The `c172x_ctl` LQR
pipeline consumes it with cosmetic changes. `LinearizedSS` the *component*
survives separately, as an ordinary
[continuous component](#g-continuous-component) in the migrated library:
`init_x` = the state vector, labeled faces, the affine update in `h_xu`/`f`.
It has no privileges, and its schema is everyone else's.

**Recorded guidance.** Linearization taps should select minimal-coordinate
mechanizations. Perturbing Euler-angle states is meaningful where seeding
quaternion components steps off the unit-norm manifold. This is why today's
code linearizes the `{NED}` variant. `design_world(ac)` rigs that variant,
promoting implicit practice to stated rule. The coordinate choice belongs to
the tap author, not the framework.

**The sampled-data Dual activation is
[recorded, not built](#g-recorded-not-built)** (a worked-out extension
deliberately left unimplemented, its seams named). The frozen-exact doctrine is
consumer-scoped, not a capability wall: today's services differentiate the
continuous dynamics with the discrete state held, and for that a frozen
discrete output — a ZOH constant with zero partials — is the exact answer. The
type system enforces it ([§11.2][s11-2]). Stated once, because the question
recurs: **the frozen discrete cell is not an AD limitation on the signal path;
it is the true zero of an instantaneous dependence that the hybrid semantics
never had**. The dataflow through a discrete component is temporal, not
instantaneous, and AD follows actual dataflow
(`frozen_discrete_walkthrough.md` works the three-component chain through).

Differentiating "through" the discrete side means differentiating a *different
object*: the sampled-data step map $\Phi : (x_k, \mathrm{slots}) \to x_{k+1}$,
taken over the model's *whole* state, continuous and discrete leaves alike. One
evaluation of $\Phi$ integrates one period, then runs the [due](#g-due)
[ticks](#g-tick). The extension is additive along existing [seams](#g-seam):

- **[Walked](#g-walked)-leaf parametrization** of the discrete tier's
  real-scalar state leaves; counters and enums stay [pinned](#g-walked), like
  `m`.
- **Opt-in participation** on discrete components, frozen-exact staying the
  default. A participating component opts in through an explicit trait, and
  that trait **brings the two-argument `T`-form of `output_types` with it**: it
  flips the leaf's mandated declaration shape from the plain form to the
  continuous one ([§11.2][s11-2], [§11.5][s11-5]). Participation therefore
  stays authored per leaf on that tier too. The hinge is recorded here so the
  two forms stay compatible — graceful migration, no flag day.
- **One new activation** ([§12.4][s12-4]): "continuous chain + `f` + discrete
  `h_x`/`h_xu` + `g`".
- **Forward sensitivities** through the in-house RK steppers, for free — a
  payoff of owning the loop ([§8.1][s8-1]).

The honest [boundary](#g-boundary): $\Phi$ is differentiable only where the
event pattern is locally constant. Exactness across a firing needs saltation
corrections. The scope is therefore event-quiescent operating points, which
trim points already are — [guards](#g-guard) at commit see to that
([§14.5][s14-5]). The scope comes with a loud diagnostic if an event fires
inside a differentiated step. Two consumers wait. The first is the closed-loop
trim door ([§14.7][s14-7]), whose $g(x) - x = 0$ residuals currently imply the
derivative-free fallback, since frozen `g` has no Jacobian columns. The second
is exact discrete-time linearization of the full loop — digital design on the
exact discretized plant instead of continuous linearization + Tustin.

**Declarative non-participation: what the schema states, and what stays
recorded.** **Both halves of this door have a spelling.** The output half
(row 166): a continuous producer's declaration is per-leaf, so "this
[port](#g-port) is frozen under differentiation" has a spelling. Declare the
leaf `Float64`, and strip with `ForwardDiff.value` inside the stage
([§11.2][s11-2], [§12.5][s12-5]). An opaque wrapper (an FMU, a C aerodynamic
table) and a deliberately severed coupling can therefore both say so in the
schema, instead of showing up in Jacobians as unexplained zero rows. The
conformance check holds them to it at every activation.

The input half (row 167): a consumer's entries are per-leaf too, so a `Float64`
entry declares "never hand me partials". That is the AD-incompatible
component's own statement, enforced at the wire ([§6.1][s6-1]). At a root slot,
such an entry *is* the forbid-seeding marker itself. That marker carries
semantics rather than mere protection: it types the slot cell at every
activation. An unseeded slot is therefore a *choice* where a `Float64`-entry
slot is a *declaration*, and tap resolution rejects the latter with the
offending entry in hand instead of returning a silent zero column.

What stays recorded is only the remaining **tooling** over that visibility.
Pinned-face validation by the tap declaration: selecting a declared-frozen
output = warning. A [feedthrough](#g-feedthrough)-graph lint: a frozen output
fed by participating inputs names the severed coupling. Both are additive when
a consumer shows up, no flag day; until then the declared pins and the visible
zero rows suffice.

---

# Part V — Grounding

## 15. Case studies

### 15.1 `Vehicle` today → this framework

This case study is the grounding exercise that validated [§5][s5]. Current
`Vehicle.f_ode!` (`aircraftbase.jl:142-170`) is a hand-woven instance of the machinery
specified here:

| Today (convention) | This design (checked structure) |
|---|---|
| `kinematics.u .= dynamics.x` — velocity extracted directly from the state vector because `f_ode!(dynamics)` can't run yet | `dyn`'s stage-1 output, scheduled first by construction; the artificial loop in `VehicleDynamics` dissolves |
| Hand-ordered `f_ode!` body (kinematics → airdata → systems → route five `dynamics.u` assignments → dynamics last) | Build-time topological sort; wrong wiring = build error naming the cycle or dangling [port](#g-port) |
| Velocity state duplicated in `dynamics.x` and `kinematics.u`, kept in sync by hand | One state, one owner; consumers wire to `dyn.vel` |
| `get_wr_b`/`get_mp_b`/`get_hr_b` generated tree-walk sums | [Summing junctions](#g-summing-junction) at ownership [boundaries](#g-boundary), one explicit wire per contributor, exported totals ([§6.2][s6-2]) |
| `f_step!` quaternion renorm + engine-phase/stall-latch checks | `project` hook + [boundary-detected](#g-boundary-detected) events with defined semantics |
| `Aircraft.f_ode!` runs avionics before the vehicle → continuous avionics reads one-stage-stale `vehicle.y` (implicit delay) | Avionics scheduled inside the [sweep](#g-sweep), after the stage-1 outputs avionics consumes — no delay. Or avionics declared periodic, sampling post-step by stated semantics |
| `atmosphere`/`terrain` threaded as arguments through every signature | Field-handle signals through ordinary ports ([§4.4][s4-4]) |

Two of those entries carry detail a cell cannot hold. The artificial loop in
`VehicleDynamics` is the pairing of a state-only velocity output with
[feedthrough](#g-feedthrough) accelerations. The hand sync of the duplicated velocity
state reaches into initialization, where `f_init!` carries the line
`dynamics.x .= kinematics.u  #essential`.

The same exercise surfaced a migration cost: today's monolithic `KinData` splits in
two, because its parts genuinely have different dependencies.

- `pose` — stage 1: `q_eb`, `r_eb_e`, `ϕ_λ_h`, ...
- `kin_vel` — stage 2: `v_eb_n`, `v_gnd`, `χ`, `γ`, ...

The recurring trade, stated once: the framework asks authors to write down structure
previously kept in their heads, and pays them back by never letting that structure
silently rot.

The genuine [algebraic loop](#g-algebraic-loop) in the domain — α̇-dependent
aerodynamics — is already broken in the current C172 model by a filter state, exactly
the explicit break [§5.5][s5-5] prescribes. That precedent is evidence that the
reject-loops policy matches domain practice rather than fighting it.

### 15.2 Torture tests for the §5.2 interfaces: `PistonEngine` and the FCS PID cascade

Three exercises, each starting from code that exists today. Two [components](#g-component) were
transliterated to validate the decoder interfaces before adoption: `PistonEngine`
on the continuous side, `PID` and the C172X FCS on the discrete one. A third
exercise takes the supervisor sitting one level above those compensators. Each
is read first as what today's code does, then as what this design makes of it.

#### `PistonEngine`: the continuous side

The current engine (piston.jl:310-449) carries a mode enum with three [flow](#g-flow)
regimes, four table lookups, two embedded continuous PI compensators, boolean
transitions and an argument-threaded `fuel_available`. The points below place
each of those features under the decoder interfaces.

- The compensator paths (`idle`, `frc`) are pure functions of the engine's own
  state `ω`. Their complete PI laws — outputs *and* state derivatives —
  therefore evaluate in `h_x`. (The alternative factoring, compensators as
  child components of an engine [assembly](#g-assembly), also [schedules](#g-schedule) cleanly from the
  core's stage-1 [ports](#g-port).)
- `h_xu` runs the lookup chain and the mode branch once. `f` is a three-field
  copy (`ω̇`, `ẋ_idle`, `ẋ_frc`). Under the orthodox split, `f` would reproduce
  essentially the whole `f_ode!` body — four lookups and the mode branch — ×4
  RK stages per step (row 15).
- `f_step!`'s transitions become [boundary-detected](#g-boundary-detected) events with mixed [predicate](#g-predicate)/threshold [guards](#g-guard)
  ([§2.1][s2-1]).
- `fuel_available` becomes an ordinary port. It is state-derived at the fuel
  system, hence a stage-1 port — no loop.
- Forced publications: none — everything `f` reads was already in `PistonEngineY`.

#### `PID` and the C172X FCS: the discrete side

`PID` (control.jl:431-471) and the C172X FCS around it are the discrete side's
representative.

- The current update entangles outputs and next state by construction. The
  spelling is `y_i = x_i`: this [tick](#g-tick)'s integral-path output *is* the updated
  integrator state.
- Under [§5.3][s5-3] the law runs once in `h_xu`, publishing paths, saturation and the
  updated states; `g` is a three-field copy.
- Under the orthodox split, `g(x, u, t)` would reproduce the entire law per
  compensator per tick (row 15).

**Discovered latent delay.** The FCS chains anti-windup: outer compensators
take `sat_ext` from the inner LQR's `sat_out` (c172x_ctl.jl:332,345,...). Wired
naively, that chain is a *genuine* tick-domain [algebraic loop](#g-algebraic-loop), and the build
correctly rejects it:

```
outer.output → inner.input → inner.sat_out → outer.sat_ext →
outer.int_halted → outer.y_i → outer.output
```

Today's code escapes the loop only through hand-managed call order: the outer
loops read the LQR's `sat_out` *before* the LQR updates, silently consuming
the **previous tick's** value. That is a unit delay existing nowhere in the
code, only in statement ordering.

Under this design the fix is one visible wire: connect `outer.sat_ext` to the
inner compensator's stage-1 port for the previous saturation, `sat_out_0`.
That port is an `x` field declared in the LQR's output [contract](#g-contract), hence
auto-published at stage-1 position ([§11.3][s11-3]). The delay becomes an explicit
property of the wiring. The loop and its fix are formalism-independent: the
framework's contribution is refusing to let the ambiguity through, and stage
1's contribution is having the delayed value already on a port.

Both components passed without blockers, with zero publications forced beyond
current practice. That result is the empirical basis for the claim ([§5.3][s5-3]) that
derivative/output overlap is the domain norm and that the decoder matches the
codebase's grain.

#### The supervisor slice: scheduled gains and bumpless engage

One level above the compensators, today's `c172x_ctl.jl` runs on two idioms
the [stores](#g-store)-and-[views](#g-view) rules deliberately remove. The first is per-tick gain
scheduling by mutation: `assign!` writes `Ref`-[cell](#g-cell) parameters from
EAS/altitude lookups every 50 Hz tick, LQR matrix sets included. The second is
mode-transition resets: `f_init!` plus a bumpless-transfer latch, hand-ordered
*before* the same tick's `f_periodic!`. Both survive as dataflow.

*Scheduled gains are inputs.* A scheduler component owns the lookup tables as
inert parameters, reads the scheduling variables as inputs, and publishes one
gain [bundle](#g-bundle) per compensator; compensators consume gains as `u`. What mutation
hid, ports expose. Gain trajectories become observable in log, [trace](#g-trace) and
[replay](#g-replay); the `Ref` writes were invisible to all three. The [feedthrough](#g-feedthrough)
graph carries the dependency. Linearization holds unseeded gain [slots](#g-slot)
constant with zero special-casing ([§14.10][s14-10]). One-shot design-time gains
(`robot2d`'s controller synthesis at init) are construction-time parameters or
stopped-sim service outputs — not a runtime write path.

*Resets are same-tick inputs, consumed in the output stage.* The supervisor
publishes `engage` and the latch value from its own feedthrough stage; the
compensator, topologically after the supervisor, honors them **this tick**:

```julia
h_xu(c::PI, (; x, u)) = (; u_cmd = u.engage ? u.u_latch : c.k_p*u.e + x.x_i)
g(c::PI, (; x, u, Δt)) = (; x_i = u.engage ? u.u_latch - c.k_p*u.e
                                           : x.x_i + c.k_i*Δt*u.e)
```

Honoring the reset only in `g` is legal, and it means something else. The
state still lands correctly at the next tick. But the *output at the
engagement tick* was already published from the stale state during the [sweep](#g-sweep),
and under ZOH the plant integrates a full step under that stale command. That
one-tick-late command is exactly the bump that bumpless transfer exists to
remove. No diagnostic can catch the bump: both spellings are meaningful
designs.

The update stage cannot rescue its own [boundary](#g-boundary) — republish-from-`x⁺` is
rejected (row 67) — so the output stage is the *only* same-tick path. Today's
hand-ordering (`f_init!` before `f_periodic!` in one call) is that contract
enforced manually. [Appendix A][sA] carries the contract as the same-tick reset
entry, and the bumpless-engage answer ([§9.7][s9-7]) presupposes exactly this
spelling — engage semantics live in the FCS.

One relative lives outside the FCS. The landing gear's level-triggered
cross-component reset (`!wow` re-initializing the friction regulator every
step) becomes an edge-triggered event owned by the regulator, a semantic
tightening recorded in the migration mapping ([§16][s16]). There the respelling is
not a stylistic one: the continuous [tier](#g-tier) admits no input spelling at all,
because only handlers write `x` ([§3.1][s3-1]). The event is therefore necessity
rather than taste, and the reimplemented `PIVector`'s optional reset [face](#g-face)
([§16][s16]) is sugar over exactly that event. [Appendix A][sA] carries the
continuous-reset contract too.

### 15.3 Torture test for the §9 staging shapes: filter, joystick and GUI

This case study is the exercise that selected per-[device](#g-device) [cells](#g-cell)
([§9.4][s9-4]) and produced the [§9.7][s9-7] [contracts](#g-contract). Setup (the
`sketch_io.jl` listing): a first-order filter with root inputs `u_cmd` and `τ`;
a fictitious 100 Hz single-axis joystick streaming a slow ramp onto `u_cmd`
(complete writer); a 60 Hz GUI with sliders for both [slots](#g-slot) (sparse
writer); 50 Hz [boundaries](#g-boundary); pace 1. The interference on `u_cmd` is
the point.

Three candidate staging shapes were on the table: **per-slot cells**, a shared
**[batch](#g-batch) stack**, and **per-device cells**. The user-level listing
came out identical across all three — ergonomics cannot discriminate between
them; behavior under a concrete interleaving did.

**Slot exclusivity rules out the very contest the setup builds.** Under slot
exclusivity ([§9.3][s9-3]) the contested-`u_cmd` scenario cannot arise — a
second writer on `u_cmd` is an attach-time error. What the test settles is
therefore the cell *shapes*: atomicity, [coalescing](#g-coalescing), pause
behavior and the peek rule. Its conflict-precedence comparison and the
active-widget stage-every-pass contract ([§9.7][s9-7]) describe a contested-slot
world the design does not have. The findings below are read under that scope.

- **Drag against the stream** — the user grabs the `u_cmd` slider while the
  joystick streams. Under per-slot cells and the batch stack, each
  [frame](#g-frame)'s conflict resolves by last-store/last-push wall-clock order
  (row 24). With 16.7 ms renders against 10 ms polls, the applied input
  alternates between drag value and ramp on the cadence beat, the filter visibly
  wobbles, and the pattern differs run to run. The [trace](#g-trace) replays any
  given run exactly; the behavior is still a timing artifact. Under per-device
  cells the GUI stages in every drag frame: ≥ one render per 20 ms frame, plus
  the active-widget contract. The GUI therefore wins every [drain](#g-drain) (the
  frame-top swap that publishes staged device inputs into the root slots) by
  attachment order. That win is a clean, deterministic override for exactly the
  grab duration. Same user code, qualitatively different physics.
- **Edits while paused.** Under per-slot cells, the user's `u_cmd` edit is
  overwritten by the still-polling joystick ~10 ms later (row 24); the knob
  visibly snaps back and the edit never applies. Under the batch stack, the edit
  is buried under newer pushes, and the pending chain grows at the polling rate
  — ~10³ nodes per 10 s pause — with every [peek](#g-peek) walking that chain
  (row 24). Under per-device cells, the `u_cmd` edit holds in the GUI's own cell
  across the pause, the knob keeping the edit by the [§9.7][s9-7] peek rule. That
  edit merges with the `τ` edit — the sparse-accumulation case — and applies at
  the un-pause drain, for one deterministic frame before the joystick reclaims
  the slot. That one-frame application is the honest semantics of one-shot
  editing a streamed input. The uncontested `τ` edit works under all three
  shapes.
- **Corrections the exercise forced.** The sparse-writer lost-write hazard is
  specific to one-cell-per-device layouts: per-slot cells cannot lose
  independent-slot writes; the CAS merge is per-device cells' antidote, not a
  general need. And the batch stack's conflict order is temporal, not an
  attachment-order policy (row 24).
- **Discoveries.** Two: the active-widget contract, and the
  [port](#g-port)-resolution answer to panel reuse ([§9.7][s9-7]) — prompted by
  asking how the filter's panel survives the filter becoming an embedded
  [component](#g-component) with `u_cmd` driven by another component. That
  embedded-filter case is the `Cessna172Xv0` → `Xv1` throttle situation.

### 15.4 The interactive C172X demo: the periphery under load

The full-fidelity successor to [§15.3][s15-3], against the real deployment:
`generic_simulation()` (`FlightApps/demos/c172_demos.jl`) —
`SimpleWorld(Cessna172Xv1, SimpleAtmosphere, HorizontalTerrain)`, GUI, joysticks,
an XPlane12 output [device](#g-device), ground/trim init, paced run, post-run plots. Method:
FlightCore's mechanisms are reference *behavior*, not requirements — the question
is whether the new machinery expresses the experience (move stick, plane banks),
never how to reproduce `assign_input!`. The interactive surface is *not* one
thing: pilot commands cluster under a prefix; environment knobs stay with their
components' panels. Inventory of the complete interactive surface, with each
item's home:

- **Streamed commands** (`throttle_axis`, `elevator/aileron/rudder_axis`): today
  written by joystick mappings after shaping *and* by GUI sliders on the same
  fields. Every dual-writer field in the demo is this pattern — a stream shadowed
  by a mirror, where simultaneous live writing is a bug. This adjudicated [slot](#g-slot)
  exclusivity ([§9.3][s9-3]): [claim](#g-claim)/disable covers every case found; none needs two
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
  preference) and command semantics (FCS design); the [face](#g-face) [contract](#g-contract) splits it —
  conditioning upstream as mapping data, semantics in-model ([§9.4][s9-4]).
- **Mode engage** (`mode_req` + setpoint capture from current measurements — the
  GUI handler does `u.EAS_ref = EAS` read from `vehicle.y`): the one place the
  GUI composes writes from model state. Resolved under *Frame anatomies* below.
- **Vehicle-direct and environment tunables** (engine start/stop/mixture, payload
  masses, terrain surface enum, sea-level T/p, wind NED): ordinary [component](#g-component)
  inputs exported to root faces; the GUI writes them under its [greedy claim](#g-greedy-claim)
  (the unclaimed complement, computed by the framework instead of returned) via
  [§9.7][s9-7]; no machinery.
- **The Xv1 actuator sliders**: FlightCore's dead sliders; resolved read-only by
  [§9.7][s9-7]. No action.
- **Outbound** (XPlane12: control-surface angles, nose-wheel steering, prop
  speed/phase, pose, `t`): a [snapshot](#g-snapshot)-consuming device, pure `map_output` on the
  device task ([§9.2][s9-2]). No friction found.
- **Init/trim, pause/pace, post-run plots**: stopped-sim services ([§14][s14]), control
  plane ([§10.1][s10-1]), log/[trace](#g-trace) ([§9.2][s9-2], [§9.5][s9-5]).

#### Architectures examined here and rejected

The [§9][s9], [§10][s10] [periphery](#g-periphery) decisions were
forced by this cast: [devices](#g-device) as [components](#g-component) (a `T16000M`
component wrapping SDL), a root-level `PilotInterface` cockpit component, and
bundled command [faces](#g-face) (`pilot_inputs` as one struct [port](#g-port)) —
all three litigated in row 45. What each leaves behind is the design's own
answer:

- the *knowledge* half of a device model — its semantics — is expressible as an
  ordinary in-model [component](#g-component) wherever wanted; only the wall-clock
  pump stays outside;
- the cockpit component's claimed jobs are covered where they belong: struct
  [assembly](#g-assembly) in-model downstream of scalar faces, curves as mapping
  data, widget arbitration by [§9.7][s9-7] plus exclusivity, and the stateful
  residue (accumulators, capture-on-engage) in the avionics;
- the routing convenience a command bundle bought in FlightCore's
  argument-threading world is provided by the namespace prefix and
  `input_passthrough` ([§11.8][s11-8]), with the struct reappearing legitimately
  downstream, assembled in-model by a single producer.

#### Surface walkthrough

The demo line by line:

- `SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_LOWS15))` —
  pure value construction; no `Model` wrapper (its jobs move into the build).
  `HorizontalTerrain`'s elevation is a plain field (parameter), its surface type
  an input [port](#g-port): the parameter/port split FlightCore kept implicit in
  `U()`-vs-field convention is now the declaration itself. The aircraft's
  `input_connections` block carries the `pilot.*` [face](#g-face) group in one place, deep routes
  spanning avionics *and* systems — today's mapping writes flaps/brakes directly
  into `act`, bypassing avionics; that bypass becomes a declared route.
- `Simulation(world; algorithm = RK4(), h = 0.02, n = 1, t_end = 1000)` — `n`
  binds `Δt_base = n·h` ([§8.5][s8-5]; default 1: base [tick](#g-tick) every step). The entire
  build pipeline runs here: [class](#g-class) resolution, path validation, face derivation
  (computed [boundary](#g-boundary) connections expanded, printable), two-producers/unconnected checks,
  topological sort, [probe](#g-probe) passes, rate compilation, flat layout, [slot](#g-slot) table.
- `init!(sim, ready_for_taxi(ac); t0 = 0.0)` — stopped-sim services ([§14][s14];
  trim is its own service, `trim!(sim, problem; baseline, ...)`, whose commit
  runs the same boundary): they write `(x, m)`, **establish every root
  slot's initial value**, and capture the [trace header](#g-trace-header). Slot initialization
  decisively belongs here, not in declarations: the trim service writes slot
  values it *solved for* (throttle, elevator) — not declaration constants.
- `attach!(sim, XPlane12Control(...), binding)` — output [device](#g-device): [claims](#g-claim) nothing,
  consumes [snapshots](#g-snapshot) via [§10.3][s10-3], pure `map_output` on its task. Its [binding](#g-binding) names
  snapshot paths, **validated at attach against the actual [contract](#g-contract)** — an
  aircraft substitution that breaks the binding fails at attach, not with silent
  garbage UDP (a new, cheap [§9.2][s9-2] obligation).
- `attach!(sim, joystick, T16000MBinding())` — the binding is a declarative
  table: axis/button → face name + conditioning params
  (`stick_y = (face = "aircraft.pilot.elevator_axis", expo = 1.0, deadzone =
  0.05)`, `button_3 = (face = "aircraft.pilot.flaps_up_count", as = :count)`).
  At attach: faces resolved against the root contract (typo → [did-you-mean](#g-did-you-mean)),
  claim set registered (second joystick on the same faces errors here). The
  Gladiator variant is the same table with different keys, zero shaping code —
  the duplication smell structurally gone.
- `run!(sim; gui = true, pace = 1)` — a [greedy claim](#g-greedy-claim) over every
  unclaimed face and liveness with zero configuration, both settled at run start
  against the [frozen roster](#g-roster) ([§9.3][s9-3]): axis mirrors read-only
  (claimed, with provenance), mode/setpoint/mixture/payload/environment widgets
  live, actuator sliders read-only ([component](#g-component)-fed). The `gui`
  flag's attachment lasts exactly this run ([§10.6][s10-6]).
  Unplug the joystick → its task exits. The mirrors stay read-only with the
  death in their provenance ("claimed by `T16000M` — task dead"), and the axes
  hold their last-drained values; those two behaviors are the accepted orphan
  anomaly ([§9.3][s9-3]). Recovery is between runs: stop, `detach!`, then
  `init!` for a fresh trajectory, or `replay!`-to-end + `run!` to continue the
  interrupted one ([§10.7][s10-7]).
  Post-run: `TimeSeries` over retained snapshots; the [trace](#g-trace) can
  re-drive a fresh `Simulation(world)` bit-identically. That replay is also the
  state-trajectory inspector (row 38 paying its way).

#### Frame anatomies

One [frame](#g-frame) each:

- *Stick motion*: [device](#g-device) task polls, conditioning helper applies binding params,
  complete [batch](#g-batch) overwrites the [cell](#g-cell) (inter-frame polls coalesce, ZOH-correct);
  [drain](#g-drain) applies + [traces](#g-trace); avionics [tick](#g-tick) reads the [slot](#g-slot) fresh; worst-case
  stick-to-physics latency = poll interval + frame, now by stated semantics.
- *Flaps click*: button [peeks](#g-peek) counter `k` (own-pending-else-[snapshot](#g-snapshot)), stages
  level `k+1` on [activation](#g-activation); drain applies; avionics compares slot counter to its
  `x` counter, moves the detent, stores. Multi-click in one window counts via
  own-pending-first peek; repeated staging idempotent ([§9.7][s9-7]).
- *Mode engage*: the GUI stages `mode_req`, plus optionally peek-captured
  setpoint slots. **Bumpless-engage semantics live in the FCS already**: the
  current `ControlLaws` latches each controller's reference from the present
  command vector on mode transitions. So the fork dissolved: semantic capture is
  aircraft design. That arrangement is the status quo, and it is uniform across
  writers — a script engages sanely staging one value.
  The GUI [peek](#g-peek)-batch ([§9.7][s9-7]) therefore
  survives as display/slot-sync sugar only. Residual check for migration:
  order-sensitivity of latch vs. sync-write on the same [boundary](#g-boundary)
  (believed none — both derive from the same measurements).
- *Wind slider*: sparse CAS-merge, the uncontested-`τ` case ([§15.3][s15-3]), live in the
  real cast.
- *Pause/un-pause*: [control plane](#g-control-plane); GUI edits hold in its cell (peek displays),
  joystick cell coalesces bounded; un-pause drain applies both (disjoint slots —
  exclusivity makes the contested question unaskable), pacer re-anchors.
- *Window close*: [§10.4][s10-4] verbatim — complete boundary, final snapshot, sticky
  stopped, wake waits, unblock hooks, named-timeout joins.

Remaining open (feeding [§16][s16]): the `q_sf` home (thin mapping entry vs.
avionics-internal derivation — aircraft design, not framework design).

### 15.5 The strapdown IMU: integrate-and-dump across the tier boundary

The strongest challenge mounted against the [§3][s3] class split, and its resolution.
The general question first: why two leaf classes at all — why not one all-in-one
primitive carrying continuous state, modes *and* discrete state, with `f`, events *and* `g`, purely
continuous or discrete [components](#g-component) falling out by whichever facets an author
declares? ([Class](#g-class) is already read off declaration shape, [§11.5][s11-5] — the question is
whether the two declaration sets should be exclusive.)

#### Why the merge buys nothing

The split is between *time bases*, not state
classes — the continuous primitive is already hybrid (`m`, [guards](#g-guard),
handlers, [§3.1][s3-1]); what separates the classes is [sweep](#g-sweep)-driven
versus [tick](#g-tick)-driven execution. And the settled rules force a merged
[component](#g-component)'s two halves to communicate exactly as two siblings
do: [one home per datum](#g-one-home-per-datum) ([§5.2][s5-2]), `f` sees only
the continuous state and `g` only the discrete one, and `x⁺` is decoded only at
the owner's next tick (`g` runs last). That deferred decode is what makes
ticks→events structurally impossible and what terminates the
[boundary](#g-boundary) iteration ([§8.6][s8-6]). Cross-[tier](#g-tier)
influence inside the merged class still routes through published table
[cells](#g-cell). The all-in-one component is therefore an
[assembly](#g-assembly) of two primitives in a trench coat, buying no
expressiveness and incurring costs of its own that row 56 enumerates. One of
those costs is not bookkeeping. The sampling [seam](#g-seam) — ZOH and the
`z⁻¹` delay — is the most bug-prone boundary in a flight-control stack. A
monolith swallows it; the split keeps it a visible wire.

#### The counterexample

The source is a pre-design FlightCore sketch, `navsensors.jl`, whose operative
content and companion derivation note are recorded here in full. In that sketch
a strapdown IMU integrates raw increments continuously — `ẋ.ϑ_c = ω_ic_c`,
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
cumulative [stores](#g-store), term by term. The reset is periodic, not
condition-triggered, so events are the wrong [tier](#g-tier); and the reset is a
discrete-tier write into continuous state, exactly the operation this design
forbids (`g` writes only its own `x`; handlers are the sole resetters of
continuous state, and they are [guard](#g-guard)-driven).
Integrate-and-dump falls squarely into the crack between the classes:
tightly-coupled continuous and periodic dynamics in one physical [device](#g-device).

#### The idiom: integrate-and-difference

The reset is eliminable by algebra, not
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
  [RHS](#g-flow) evaluation, RK stages included, applies the current cumulative attitude,
  exactly as the sketch applies the current `q_c_cc`.

#### Exactness condition, stated once

Interval-relative integrals factor into
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

The `IMU` [assembly](#g-assembly) wires the four integral [ports](#g-port) across, holds the error model as
a discrete sibling consuming `sample`, and leaves the sampler at `K = 1` in its
own scope — the parent sets the IMU's rate ([§11.7][s11-7]). `Δt` in the stage [bundle](#g-bundle)
(the NamedTuple of zero-copy views a component function receives) is the [§8.5][s8-5]
single source of truth, put there for exactly this kind of discretized law.
(Initialization consistency — the sampler's `x` must equal the initial
integrals or the `t₀` sample is wrong — holds by default at zeros/identity, and
[boundary zero](#g-boundary-zero) discharges the rest: its [due](#g-due) `g` latches `x ← integrals(t₀)` for
every subsequent sample, so only the `t₀` sample itself depends on the authored
`x` — a [condition](#g-condition)-authoring obligation under trim, [§14.5][s14-5].)

#### Why `u.V` is fresh — the line that would silently zero

The sculling line is
correct only because a due [tick](#g-tick) samples the *completed* [boundary](#g-boundary): if `u.V` still
held the previous boundary's decode it would equal `x.V` exactly (that is the
value `g` latched), and sculling would vanish without an error anywhere. The
guarantee is the [§8.6][s8-6] macro-sequence, not a scheduling accident: integrate →
project → [sweep](#g-sweep), with the due sampler's stages gated *into* that sweep ([§8.5][s8-5]) and
the integrals arriving at stage-1 position (auto-published state, [§5.3][s5-3]) — before
any stage-2 function runs, regardless of topological placement. The rest of the
timeline closes consistently: the sampler's `h_xu` decodes `x` (the `t_{k-1}` latch)
*before* `g` runs — the `z⁻¹` semantics — and after event [quiescence](#g-quiescence) `g` latches
the `t_k` values for the next tick; same-boundary events re-run the gated stages
in their re-sweeps, so `g` and external readers see the settled boundary.

#### Author-knowledge note

User observation, recorded as a documentation
obligation: the clean implementation leans on the author *knowing* that "sampling
at `t_k`" means post-integration, post-[projection](#g-projection), stage-1-fresh state. That
knowledge must be part of the framework's taught [contract](#g-contract) — [§8.5][s8-5]/[§8.6][s8-6] semantics
stated in [component](#g-component)-author documentation with this IMU as the [worked](#g-worked) example — not
internal lore. The failure mode of not knowing it is instructive: an author who
distrusts the [sweep](#g-sweep) order adds a defensive one-[tick](#g-tick) delay or re-derives the
integrals in the sampler, silently degrading the model.

**When the coupling is genuinely two-way: the latch-back wire.** The IMU's
coupling is one-directional (integrals → sampler). If the [flow](#g-flow) itself needed the
interval-relative value — integrator saturation within the sampling interval,
say — the latch becomes a wire back: the sampler publishes the sample-instant
values from its *[feedthrough](#g-feedthrough)* stage (`h_xu` reads `u`, so the latch [port](#g-port) carries
the current tick's values, ZOH until the next; an `h_x`-published latch would be
one period stale), and the continuous `f` computes `x − u.latch`. Both cross-wires
consume the other side's ports and the [schedule](#g-schedule) stays acyclic (integrals stage 1 →
sampler `h_xu`; sampler `h_xu` → the integrals' `f`-edge, [§5.4][s5-4]). The "reset"
becomes a visible [tier](#g-tier)-crossing feedback loop — which is what it always was,
physically.

**Verdict.** The strongest counterexample landed on the two-class taxonomy with
*less* code than the fused original — same thirteen integral scalars, same math,
minus the reset block — and three structural gains. The sampling [seam](#g-seam) became a
wire. The sketch's incidental violations became visible structure: the
`CircularBuffer` mutated inside the component struct (constants) moves to the
consumer's `x` or falls out of the log; the parent-called `f_disc!(errors)`
becomes a discrete sibling, making the truth/corrupted sample pair separately
loggable. And linearization got sane: under a `Dual` [activation](#g-activation) the discrete tier
is held ([§11.2][s11-2]), and "integrators that never reset" *is* the cumulative
formulation — the framework's rules pushed the model into the only form its own
linearization semantics could coherently handle. Residual escape hatch, recorded
unbuilt: if interval-relative dynamics ever neither factor algebraically nor
tolerate the latch-back wire, the [guarded addition](#g-guarded-addition) is a **tick-triggered handler**
on [continuous components](#g-continuous-component) (periodic events). Nothing surveyed needs it (row 56).

---

## 16. Open axes

Three axes are still to be settled: the migration of FlightPhysics and
FlightApps, the GUI panel authoring API, and log and [trace](#g-trace)
persistence.

#### Migration

What follows is an outline for FlightPhysics/FlightApps, not a specification.
The table carries one row per item: the item, the disposition recorded for it,
the section owning the machinery it touches, and the governing decision row. A
dash means the outline names the item and records nothing further. Items whose
disposition exceeds a cell are expanded below the table.

| item | disposition | section | decision row |
|---|---|---|---|
| The [walked](#g-walked)-leaf parametrization pass | the `Ranged` rewrite targets the walk rule wherever `Ranged` survives, at ports and parameters | [§11.2][s11-2] | — |
| The `KinData`-style output splits | — | — | — |
| The contributor survey feeding the aggregation chains | mechanical to extract from today's trait implementations | [§6.2][s6-2] | — |
| Comparison criteria against FlightCore's demonstrated strengths | three strengths to compare against: zero-alloc stepping, flexibility, interactive operation | [§12.7][s12-7] | — |
| The [component](#g-component) library's starting inventory | — | [§13.7][s13-7] | — |
| The conventional exported aircraft surface for generic periphery consumers | pose and velocity faces with wrapper types, the periphery-facing half of the `KinData` successor | [§9.2][s9-2] | — |
| The supervisor seam | three respellings — gain ports and schedulers, mode-transition latches, the gear's reset — the last of which lands on the *library* side | [§15.2][s15-2] | rows 89, 139 and 141 |
| The steering contract re-factoring | `AbstractSteering` moves from "give me the angle" to `(engaged, ψ_cmd)` | [§5.4][s5-4] | — |
| Splitting `Strut` | the residual remedy, recorded and not taken | — | — |
| The state-declaration conversion to the closed vocabulary | each `RQuat` state field becomes its `SVector{4}` backing, each `Ranged` state field a plain scalar | [§7.1][s7-1] | — |
| The exported-name surface | decided deliberately rather than by accident, by a full-surface audit under the four-register naming convention | [§14.2][s14-2] | row 144 |
| The [executor](#g-executor) compile-cost re-measurement | runs on the real vehicle skeleton, early — before the executor's shape hardens | [§12.7][s12-7] | — |
| *Residual*: the `q_sf` home | aircraft design, so it belongs on this list | [§15.4][s15-4] | — |
| *Residual*: a root-declared overridable `stop_on` default | reopen only if the constructor argument proves chronically forgotten | [§13.5][s13-5] | — |

**The parametrization pass.** `Ranged` survives at [ports](#g-port) and
parameters, and there the rewrite targets the walk rule ([§11.2][s11-2]):
constructor discipline admitting the walked scalar with the value parameters
left alone, plus a `probe_value` method. State fields are not among those
survival sites; the state-declaration conversion below turns each `Ranged`
state field into a plain scalar.

**Comparison criteria.** FlightCore's demonstrated strengths are three:
zero-alloc stepping, flexibility, interactive operation. Zero-alloc stepping is
measured through the `phase_bodies` [seam](#g-seam) ([§12.7][s12-7]),
apples-to-apples with today's `@ballocated f_ode!` suites.

**The conventional exported aircraft surface.** Generic
[periphery](#g-periphery) consumers read the integration
[register](#g-register) ([§9.2][s9-2]). What that surface exports is pose and
velocity [faces](#g-face) with wrapper types — `VelocityData`, field meaning
defined at the type — the periphery-facing half of the `KinData` successor.

**The supervisor seam.** The supervisor sitting above the compensators
([§15.2][s15-2]) contributes three respellings. Compensator gains become input
ports fed by scheduler components (~7 for the C172X). Every mode-transition
latch is respelled as a same-[tick](#g-tick) reset. The gear's level-triggered
reset becomes an edge event — and that last one lands on the *library* side.

On the library side, the reimplemented `PIVector` gains a **flag-gated reset
face**. `PIVector(; reset = true)` adds a `Bool` input face plus the event; the
default omits both. Declarations are ordinary functions of the instance
([§11.5][s11-5]), which is what makes this the honest version of Simulink's
checkbox. One fixed policy governs the face: a rising edge resets to the
declared `init_x` values. The implementation is internal — an ordinary
[guard](#g-guard)/handler event, the continuous-reset [contract](#g-contract) in
its [worked](#g-worked) instance ([Appendix A][sA]).

Falling-edge consumers wire a NOT gate (the Bool gates, [§13.7][s13-7]).
Level-pinning and reset-to-an-external-value, which is tracking, are different
blocks rather than options on this one (row 141).

The gear then wires `strut.wow → frc.reset`. That is the **touchdown** edge —
the not-[holding](#g-edge-semantics) → holding semantics ([§2.1][s2-1]) — and it
gives fresh regulator state per contact episode. The liftoff edge (`!wow`) was
rejected (row 141). [Boundary-detected](#g-boundary-detected) policy (checked
for edges at step boundaries only, no root-finding) suffices, because the
regulator's input ramps from zero at touchdown, so localization buys nothing. A
sim initialized on ground fires the reset at [boundary zero](#g-boundary-zero)
(the initialization boundary: the ordinary macro-sequence with an empty
integrate). It fires harmlessly there: declared inits are zero, and
[boundary](#g-boundary)-zero [priors](#g-prior) are not-holding
([§14.5][s14-5]).

The engine's two `PIVector` instances — `PistonEngine`'s `idle` and `frc` —
migrate **unchanged, flag off**. They are verified reset-free in today's code,
where windup across unused phases is already handled by the saturation bounds
and `int_halted`. Their `f_init!` gain writes become construction-time
parameters, as `Contact`'s do (row 89). The PI *law* is shared as plain pure
functions called by the block's stages, the laws-as-plain-functions pattern
(row 139). `sat_ext` poses the same always-on-vs-flag-gated face question, to
be decided at reimplementation time on the same axis.

**The steering contract re-factoring.** This is the middle rung
([§5.4][s5-4]), worked on the shipped instance. `AbstractSteering` moves from
"give me the angle" to `(engaged, ψ_cmd)`, with the castoring fallback
`ψ_sw = engaged ? ψ_cmd : ψ_v` computed inside `Strut`. That move deletes the
strut → steering → strut artificial loop that stage-2 conservatism would
otherwise manufacture. The `VehicleDynamics` instance standing beside it
([§15.1][s15-1]) needs no such move: it dissolves under the two-stage split
alone.

**Splitting `Strut`.** The residual remedy is to split `Strut`, its shared
geometry crossing the new boundary as one `StrutGeometry` [bundle](#g-bundle)
port. It is recorded and not taken. The call is an aircraft-library one — a
component's own contract — recorded here rather than in framework vocabulary.

**The state-declaration conversion.** State declarations move to the closed
vocabulary ([§7.1][s7-1]). Each `RQuat` state field becomes its `SVector{4}`
backing, with the explicit `normalization = false` cast at its use sites; the
4-wide rate is already what today's `Attitude.dt` delivers. Each `Ranged` state
field becomes a plain scalar, its clamp respelled as dynamics or
[projection](#g-projection), never as construction.

**The exported-name surface.** This surface is to be decided deliberately
rather than by accident. `condition`, `fragment`, `at`, `capture` and the
`merge` overload ([§14.2][s14-2]) are generic names sharing a namespace with
FlightPhysics domain code. `merge` in particular is a piracy surface, and its
mixed-argument methods must stay error methods. For the readers, the `get_`
prefix of the [selector](#g-selector) family already settles the question
([§14.4][s14-4]). Whether the [condition](#g-condition) algebra ships behind a
submodule is the packaging question.

The audit is a full-surface [sweep](#g-sweep) (per user, 2026-08-01). Every API
method name is either specific enough to export, or gets renamed, or is left
unexported — and for extension-only surface, *unexported* is the preferred
disposition. Extension-only surface has three parts:

- the declaration and stage family of the import list ([§11.1][s11-1]), the
  larger half of the question: it sits on every component file's first line and
  is settled there;
- the [binding](#g-binding) interface `claims`/`reads` ([§9.6][s9-6]) and the
  side traits `is_input`/`is_output`/`is_greedy`, with `map_input`/`map_output`
  outside the question as loop-idiom conventions the framework never calls;
- the [device](#g-device) contract
  `init!`/`loop`/`shutdown!`/`unblock!`/`needs_calling_task`, which authors
  extend by `import` or qualified name, `Base.show`-style, rather than call
  every day.

The audit's criterion is the **four-register naming convention** (row 144):

1. **Declarations**, which the author defines and the framework calls, are bare
   nouns or `init_*`/`_types`: `child_connections`,
   `input_connections`/`output_connections`, `events`, `input_types`,
   `workspace`, the stage letters, and `claims(b)` from the binding interface
   ([§9.6][s9-6]).
2. **Value selectors**, called against `reads` and [snapshots](#g-snapshot),
   carry `get_` ([§14.4][s14-4]).
3. **Lifecycle and mutating actions** are verbs, with `!` when they mutate.
4. **Build primitives** ([§13.3][s13-3]) are plain verbs.

A name in the wrong register is a rename candidate on that ground alone.

The convention also has a **semantic axis**: right register, wrong noun.
`input_passthrough` ([§11.8][s11-8], row 171) and the binding methods
`claims`/`reads` ([§9.6][s9-6], row 146) are what settle it — bare-noun
declarations name the *consequence* a declaration has rather than its
*content*. `exports` is that axis's retired exemplar (row 170). The
`*_connections` family names content deliberately, for authoring transparency;
that is a recorded choice, not register drift.

Five items are flagged for the sweep and deliberately not settled here:

- `input_faces`/`output_faces` — noun accessors punning on the `_types`
  declarations, mitigated by being framework-facing.
- `workspace` — a declaration whose bare noun reads as an accessor; every
  candidate replacement is clunkier, so lean keep.
- `loop` — the device contract ([§9.6][s9-6]): a mutating task body spelled as
  a bare noun among its verb-`!` siblings `init!`/`shutdown!`/`unblock!`. With
  `run!` taken and the "loop body" prose entrenched, it needs the audit's
  whole-surface view.
- The bare-noun accessor family `trace(sim)`, `latest(sim)`, `binding(handle)`,
  `phase_bodies(sim)` — value selectors outside register (2)'s `get_` rule.
  `trace` is the sharpest of them: the constructor kill-switch `trace = false`
  and the post-run accessor `trace(sim)` are one name in two senses, the
  overload pattern rows 122 and 144 retire.
- Whether register (1) needs an explicit exemption for predicate traits
  (`is_greedy`, `needs_calling_task`).

All five are boundary cases row 144's convention does not settle, and they
are not defects of its list.

#### GUI panel authoring API

The semantics are settled ([§9.7][s9-7]): derived liveness, first-class
read-only rendering, own-pending-else-snapshot [peek](#g-peek),
[stage-on-interaction](#g-stage-on-interaction), orphan display. What is
deferred to migration is the calling convention — context contents, port
naming, child composition. That convention is to be co-designed against the GUI
library under the four constraints ([§9.7][s9-7]).

#### Log and trace persistence

The in-memory artifacts are settled; nothing on-disk is. Three facts stand on
the in-memory side:

- the log is the retained boundary snapshots ([§9.2][s9-2]);
- the input trace is always on and device-tagged, carrying its header of
  initial [stores](#g-store) and [slot](#g-slot) values ([§9.5][s9-5],
  [§14.5][s14-5], [§14.6][s14-6]);
- the primary/derived rule holds: the log is recomputable from the trace, never
  the reverse.

The on-disk questions are deferred to migration, where the consumers exist to
ground the choices:

- the HDF5 export scope — the whole snapshot log, or selected subtrees;
- field-handle summarization over retained snapshots, the successor to the
  `getproperty` navigation of `TimeSeries`, which is today's post-processing
  entry point;
- the trace file format, which doubles as the reproducibility carrier: the
  [replay](#g-replay) pointers ([§13.4][s13-4]) name positions in it.

---

## Appendix A. Taught contracts: the author-facing index

The build pipeline enforces structure — declarations, wiring, types,
conformance. A residue of *semantic* facts is unenforceable by any check.
Knowing them is what makes component and periphery code come out right.
Not knowing them produces defensive delays, duplicated math or mistimed
samples with no diagnostic firing anywhere; the author-knowledge note
([§15.5][s15-5]) is the archetype.
This appendix is an **index, not a second home**: one recall line per
contract, with the normative statement staying in the owning section. That is
one home per datum, applied to the document itself.

For component authors:

- **The stage funnel** ([§5.2][s5-2]). Stage name ⊇ bundle ⊇ destructured reads:
  the stage name fixes the maximal legal view set, the component's
  declarations narrow it to the bundle, and the signature's destructuring
  narrows the bundle to actual reads — so "no `u` in the name" *is* the
  no-feedthrough property. The teaching line: stage 1 publishes what you
  know from state alone; stage 2 adds what needs inputs; your dynamics
  read your own published results instead of recomputing them.
- **One home per datum** ([§5.2][s5-2], [§4.3][s4-3]). The signal table holds *produced*
  signals only, never transported ones: buffer for continuous `x`, stores for
  discrete `x` and for `m`, table for signals — no store mirrors another.
- **The value-level constructor** ([§4.4][s4-4]). A field-emitting component ships
  the map (component, input values) → handle as a plain exported function,
  and its output stage merely calls it: the condition math ([§14.1][s14-1]) must
  be able to produce the sweep's exact handle outside any sweep, and only the
  component's author can write that function without re-creating the
  drift class.
- **Boundary sampling** ([§8.5][s8-5]/[§8.6][s8-6]; worked example [§15.5][s15-5]). "Sampling at
  `t_k`" means post-integration, post-projection, stage-1-fresh state: a
  due tick's gated stages run inside the boundary sweep and sample the
  *completed* boundary. Distrusting that guarantee — a defensive one-tick
  delay, a re-derivation inside the sampler — silently degrades the model.
- **Interval alignment** ([§14.5][s14-5]). A boundary's `g` is the *outgoing*
  transition: at tick `t_k` it consumes the completed boundary's samples
  and produces `x_{k+1}` — the value the component's *next* tick decodes
  (the sampled-data `z⁻¹` delay, by construction). Hence `g` runs at
  boundary zero: that run is the `t₀` sample's only chance.
- **Same-tick reset consumption** ([§15.2][s15-2]) — *discrete tier*. A commanded reset
  of a discrete component's `x` is an input. For same-tick output semantics
  the *output stage* consumes that input — overriding the state-derived path —
  and `g` stores the matching `x⁺`. A reset honored only in `g` reaches the
  outputs one tick late: the plant integrates a full step under the stale
  command. Both spellings are legal; they mean different things. The
  continuous tier has no such choice — next entry.
- **A continuous component's state reset is an event** ([§3.1][s3-1], [§8.6][s8-6], [§15.2][s15-2]).
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
  stale-output hazard to manage: [§8.6][s8-6] re-sweeps outputs to quiescence after
  handlers, so a continuous edge-reset is same-boundary by construction.
- **Guard predicates, edges and priors** ([§2.1][s2-1], [§8.6][s8-6]). A guard defines
  a predicate — a `Bool` form, or a sign value `σ` with
  positive = holding — and the form chosen *is* the detection policy: `Bool`
  boundary-detected, sign localized. Events fire on not-holding → holding *edges* against per-event
  priors (the previous boundary's quiescent sample): a predicate that
  keeps holding fires once, at the boundary where it first held. Boundary
  zero sets every prior to not-holding, so a predicate already holding in
  the authored state fires at `t₀`. The opposite crossing direction is a
  second event with the negated guard.
- **Handler-phase visibility** ([§5.3][s5-3], [§8.6][s8-6]). A handler executes
  against exactly the world its guard fired on: own `y`, foreign `u` and own
  `x`/`m` are all the firing round's sweep, and a component fires at most one
  event per round — later own events are re-decided against the next round's
  sweep, one round per causal link, within and across components alike. The
  signal table is written only by sweeps.
- **Stage totality** ([§12.3][s12-3]; [§13.4][s13-4], [§13.5][s13-5]). Stage code is total over
  type-valid inputs: the probe evaluates every user function against values
  chosen for their types alone, and a value-level throw is a build failure
  there and a `StepError` at runtime. Physical plausibility is a published
  `Bool` and `stop_on`; self-consistency asserts belong in tests; parameter
  validation belongs at instance construction, not inside a stage.
- **Stop-face sampling** ([§13.5][s13-5]). Stop faces are read in completed-boundary
  snapshots; declare a sign-form (localized) event if the stop needs localizing.

For periphery authors and consumers:

- **Levels, never deltas** ([§9.4][s9-4]). Staged input values are levels
  (`press_count = 17`, never `presses += 1`) — idempotent under
  coalescing; button edges ride as monotonic counters. Cross-datum state
  (press counters, edge detection) lives in the device struct, maintained
  by the loop, arriving *inside* the datum — `map_input` is pure ([§9.6][s9-6]).
- **The device loop idioms** ([§9.6][s9-6], [§10.4][s10-4]). Loop on `running(handle)`;
  make every blocking call interruptible (an `unblock!` override, or
  timeouts); voluntary exit is returning. Three canonical shapes:
  timer-poll (sleep, poll, stage), source-driven (block on your socket;
  `unblock!` closes it), boundary-driven (`wait_next_snapshot`, gather,
  send). A forgotten predicate check surfaces as `DeviceJoinTimeout` with
  your device's name; a stall as a stale heartbeat.
- **`shutdown!` closes only what is open** ([§9.6][s9-6], [§10.4][s10-4]). The framework
  runs `shutdown!` on every exit path, your own `init!`'s failure included:
  a throw half-way through acquisition hands the half-built device straight
  back to you, so guard each release (`isopen`, a `nothing` handle) rather than
  assuming initialization completed. The converse is a burden you do *not*
  carry: `init!` owes no cleanup of its own.
- **Binding traits are declarations, mappings are your own idiom** ([§9.6][s9-6]).
  Keep `is_input`/`is_output`/`is_greedy` trivial — a literal, or a flag read
  off a field fixed at the constructor call — because the framework calls them
  once, at attach, and cross-checks each against the enumeration method it
  implies. `map_input`/`map_output` are the other kind of thing: conventions of
  the loop idiom, called only by your own `loop`, never by the framework — the
  names are worth keeping for readers, and nothing enforces them.
- **[Bad datum](#g-bad-datum) vs. bug** ([§9.6][s9-6], [§13.4][s13-4]). Catch what your parser can throw,
  `report!(handle, MalformedDatum(cause))`, stage nothing, continue; let
  everything else propagate — the wrapper makes it `DeviceCrash`.
  Tolerating everything hides bugs as "device attached, nothing happens";
  tolerating nothing kills a live link on its first truncated datagram.
- **Derived liveness** ([§9.7][s9-7]). A widget is live iff its port's feed chain
  terminates in a root slot inside the GUI's own claim in the run's frozen
  partition; there is no per-port marking, and unexported ports are
  unpokeable.
- **The two observation registers** ([§9.2][s9-2], [§13.5][s13-5]). A deep snapshot path is
  the *inspection* register: it sees everything and promises nothing
  across builds. An exported output face is the *integration* register:
  curated, writer-independent meaning — the only shield against silent
  semantic drift. Bind faces in anything meant to outlive the current
  build. The store selectors (`get_state`/`get_deriv`) belong to neither:
  they read live stores, never snapshots (the source rule, [§14.4][s14-4]).

---

## Appendix B. API synopsis: the entry points

The user-facing surface on one page — same rule as [Appendix A][sA]: an index, not a
second home, with each signature normative only where its owning section settles
it. The author-side declaration surface first, then the operator surface by
lifecycle:

**Authoring** — what a component or assembly defines ([§11.2][s11-2], [§11.5][s11-5]–[§11.7][s11-7]):

- Continuous leaf: `init_x`/`init_m` (by value), `workspace(::C, ::Type{T})`
  (by allocation), `input_types(::C, ::Type{T})` and
  `output_types(::C, ::Type{T})` (by type),
  `events` — stages `h_x`, `h_xu`, `f`, guard/handler pairs
  (`Event(guard, handler)`; detection policy comes from the guard's return
  type, [§8.4][s8-4]), `project`.
- Discrete leaf: `init_x`, `workspace(::C)`,
  `input_types`/`output_types` — stages `h_x`, `h_xu`, `g`.
- Assembly: `child_connections` (mandatory — the class marker),
  `input_connections`, `output_connections`, `sample_times`.
- Shipped conditions: `condition(::C; kw)` fragment functions ([§14.2][s14-2]).

Bundle contents by function family (the maximal legal sets, [§5.2][s5-2] — signatures
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

Table footnotes, from the bundle law ([§5.2][s5-2]) — the sets above are maximal, and
each field is present only if it exists for the component: `u` iff the function
family may see inputs **and** the component declares `input_types`; `y` iff the
component produces any table cell (`output_types` ∪
auto-published); `x`/`m`/`ws` iff declared; `y_x` iff the stage-1
*return* is non-empty (auto-published names excluded — [§5.2][s5-2], row 169);
`w` iff the stage handing it down returned one (the one-hop law, [§5.2][s5-2]);
`Δt` on the discrete tier only. Returns: a stage returns a NamedTuple of
port values, or the pair `(y, w)` adding its private intermediates
([§4.3][s4-3], [§5.2][s5-2]); `f` returns the layout image of `X` ([§7.1][s7-1]); a **handler
returns `(; x, m)` with each key present iff that store exists and the handler
updates it** (the return law, [§5.2][s5-2] — no padding, `x` complete, `m` partial).

**Build.**

- `build(world) → Build` — standalone; the inspectable contract artifact:
  wire list, face table with provenance, schedule, root slots ([§12.2][s12-2]).
  `build(world; activations = (Float64, ProbeDual))` additionally pins
  activation invariants for CI (`ProbeDual` the exported canonical concrete
  probe scalar, [§12.4][s12-4]), and pre-materializes activations so a parallel
  sweep shares a fully immutable `Build` ([§9.1][s9-1], [§12.4][s12-4]).
- `resolve(asm, path) → AbstractComponent` — the getfield walk along `/`
  segments, enforcing the generic-boundary rule ([§6.1][s6-1]) at the primitive ([§13.3][s13-3]).
- `input_faces(c)` / `output_faces(c) → Vector{String}` — declaration-ordered
  face names ([§13.3][s13-3]).
- `input_passthrough(asm, path; prefix, sep, except, only)` — the declaration-site
  helper for computed input connections ([§11.8][s11-8]).

**Deployment.**

- `Simulation(world; algorithm = RK4(), h, n = 1, Δt_base = nothing,
  t_end = Inf,
  stop_on = (), localization_tol = 1e-6, localization_budget = 8,
  firing_budget = 4,
  trace = true, log = true, log_every = 1, log_max = 65536)` —
  wraps the build (`Simulation(world; ...) = Simulation(build(world); ...)`;
  the `Build` overload takes the same deployment keywords and deploys an
  inspected artifact directly, [§12.2][s12-2]).

  | keyword | default | meaning | owning section |
  |---|---|---|---|
  | `algorithm` | `RK4()` | the stepper; `RK4` is the default | [§8.2][s8-2] |
  | `h` | — | required: a domain rate is not a framework default | [§8.2][s8-2] |
  | `n` | `1` | absent the `Δt_base` keyword, the `n·h` product is the base tick period (the default path); given it, `n` is instead derived and validated an integer ≥ 1 | [§12.1][s12-1] |
  | `Δt_base` | `nothing` | the base tick period as a `Rational`, `Period` or `Hz` value; one of three binding sources | [§12.1][s12-1] |
  | `t_end` | `Inf` | the run's end time — a **default**, overridable per run at `run!` | [§13.5][s13-5], [§10.6][s10-6] |
  | `stop_on` | `()` | root-exported `Bool` output faces, OR-combined — a **default**, overridable per run at `run!` | [§13.5][s13-5], [§10.6][s10-6] |
  | `localization_tol` | `1e-6` | the root-finder's relative bracket-width convergence test (`localization_tol · h`) | [§8.4][s8-4] |
  | `localization_budget` | `8` | the per-frame localization allowance | [§8.4][s8-4] |
  | `firing_budget` | `4` | the per-event, per-boundary firing allowance of the event iteration, an integer ≥ 1 | [§8.6][s8-6] |
  | `trace` | `true` | the input trace's plain kill switch | [§9.5][s9-5] |
  | `log` | `true` | the snapshot log's plain kill switch | [§9.2][s9-2] |
  | `log_every` | `1` | the log's keep-every-kth decimation | [§9.2][s9-2] |
  | `log_max` | `65536` | the maximum number of retained snapshots, finite by default with `Inf` the opt-out | [§9.2][s9-2] |

  `Δt_base` binds from exactly one of three sources ([§12.1][s12-1]): the
  `Δt_base` keyword — a `Rational`, `Period` or `Hz` value, `n` then derived
  and validated an integer ≥ 1 — the `n·h` product when the keyword is absent
  (the default path), or, in a fully anchored model omitting both, derivation
  from the constraint pool at the coarsest admissible value, printed with its
  drivers ([§12.2][s12-2]).

  `t_end = Inf` is the honest interactive default — open-ended in time but
  bounded in memory, `log_max` being what keeps such a session from growing
  without limit ([§9.2][s9-2]). A run with no finite `t_end`, no `stop_on`
  faces and `pace = Inf` warns at start, an unbounded unattended run being
  almost always an oversight. A run ends at the first grid boundary reaching
  or exceeding `t_end`, whole frames only ([§10.4][s10-4]). The `stop_on` faces
  are recorded in run metadata — the trace header's deployment block
  ([§9.5][s9-5], [§13.5][s13-5]; walkthrough [§15.4][s15-4]).

  An event that exhausts `firing_budget` at a boundary loses its further edges
  there, under a `FiringBudget` warning ([§8.6][s8-6]). `localization_tol`,
  `localization_budget` and `firing_budget` are all three
  trajectory-determining like their siblings, hence validated with them
  (`DeploymentInvalid`) and recorded in the deployment block, where replay
  compares them ([§9.5][s9-5], [§10.7][s10-7]).

  Recording: `log_every` is admissible on the derived artifact only, never on
  the trace ([§9.2][s9-2], [§9.5][s9-5], row 29). When the log fills, the
  retention stride doubles, so the whole run stays covered at coarsening
  density, the boundary-zero and terminal snapshots being retained
  unconditionally and outside the bound ([§9.2][s9-2]). Both
  are view policies, not trajectory-determining: neither enters the deployment
  block, and replay neither records nor compares them.
- `attach!(sim, dev::AbstractDevice, binding::AbstractBinding; should_abort = false)`
  — the roots are mandatory and the signature is the gate.

  | argument | default | meaning | owning section |
  |---|---|---|---|
  | `dev::AbstractDevice` | — | the device instance; the root type is mandatory | [§9.6][s9-6] |
  | `binding::AbstractBinding` | — | the binding value; the root type is mandatory | [§9.6][s9-6] |
  | `should_abort` | `false` | the per-attachment failure policy: set, the device's departure also requests a sim stop; clear, the run continues with the device absent and its claims held to run end | [§9.6][s9-6], [§10.4][s10-4] |

  A departure is the loop body returning, a crash, or a failed `init!`.

  **Sides are declared** by the binding's Bool traits, and each declared side
  carries its own methods:

  | trait / method | default | meaning | owning section |
  |---|---|---|---|
  | `is_input(b)` | `false` on `AbstractBinding` | declares the input side | [§9.6][s9-6] |
  | `is_output(b)` | `false` on `AbstractBinding` | declares the output side | [§9.6][s9-6] |
  | `is_greedy(b)` | — | `true` switches the claim's *source* | [§9.3][s9-3], [§9.6][s9-6] |
  | `claims(b)` / `map_input(datum, b)` | — | the input side: the enumerated face set *is* the claim — what the device may write, not what it will | [§9.4][s9-4] |
  | `reads(b)` / `map_output(nt, b)` | — | the output side | [§14.4][s14-4], [§9.2][s9-2] |

  The conformance check runs at attach, pairing each trait against its method: error fallbacks
  for a declared side whose `claims`/`reads` was never written,
  `which`-against-the-fallback for a method defined under a false trait, both
  `BindingContractMismatch` ([§9.6][s9-6]). A claim is registered with
  exclusivity enforced, and the staged shape and normalization shim are
  compiled ([§9.4][s9-4]); the `reads` selectors are validated and compiled to
  one gather ([§14.4][s14-4], [§9.2][s9-2]). Under `is_greedy(b) = true` the
  framework computes the unclaimed complement at attach instead of calling
  `claims`, everything downstream being identical, an empty remainder legal and
  reported (`EmptyGreedyClaim`), and `is_greedy` without `is_input` an error
  ([§9.3][s9-3], [§9.6][s9-6]). `TableBinding` is the shipped data-driven
  binding, the standard GUI binding the shipped greedy one ([§9.6][s9-6]).

  `attach!` is a stopped-sim operation — legal in `built`, `initialized` and
  `stopped`, an error while `running` (`ServiceLifecycle`; the roster freeze,
  [§9.3][s9-3]). Admission checks identity (`AlreadyAttached` — one roster
  entry per instance, rebinding = `detach!` + `attach!`), calling-task affinity
  (`CallerTaskConflict` — at most one holder) and claims (`ClaimConflict`),
  [§9.3][s9-3]. It registers only: the task appears at the next `run!`.
- `detach!(sim, device)` — removes the roster entry and releases the
  device's claims; stopped-sim only, like `attach!`. A loop body's
  voluntary exit or crash mid-run does *not* detach: the task dies, the
  claims persist to run end ([§9.3][s9-3], [§9.6][s9-6], [§10.4][s10-4]).
- The device contract — `MyDevice <: AbstractDevice` plus `init!(dev)` /
  `loop(dev, handle)` /
  `shutdown!(dev)` / optional `unblock!(dev)` / optional trait
  `needs_calling_task(dev) = false`, a trait the task topology admits at most
  one of per roster and whose device runs its loop body inline on the calling
  task ([§9.1][s9-1]). Around those functions: per-run `init!`
  on the calling task — bracketed, so a throw there is `shutdown!` plus
  `DeviceCrash` by name and the device is dead from boundary zero
  ([§10.4][s10-4]) — the author-owned task body inside the framework's
  try/catch/finally wrapper, voluntary exit = return ([§9.6][s9-6], [§10.4][s10-4]).
- The device handle — one type, capabilities not taxonomy: `running`,
  `latest`, `wait_next_snapshot` ([§10.3][s10-3]), `stage!`, `binding`, `gather`,
  `report!` ([§9.6][s9-6]).

**Condition algebra** ([§14.1][s14-1]–[§14.6][s14-6]).

- `fragment(; x, m, slots)` — self-vocabulary payloads at the authoring
  level; `slots` names faces of that level's contract.
- `at(prefix, node)` — scoping; stores, never applies. Also lifts whole
  `TrimProblem`s and linearization tap sets ([§14.9][s14-9], [§14.10][s14-10]).
- `merge(nodes...)` — symmetric collection; duplicate leaves error with dual
  provenance; mixing a node with a bare NamedTuple is an error method, not
  `Base.merge` ([§14.2][s14-2]).
- `override(base, patches...)` — ordered layering; patch wins, provenance
  keeps both ([§14.6][s14-6]).
- `condition(comp; kw)` — the shipped fragment-function idiom; aircraft
  baselines (`ready_for_taxi(ac)`, `cold_and_dark(ac)`) are its
  full-coverage instances.

**Stopped-sim services** ([§14][s14]).

- `init!(sim, condition; t0 = 0.0)` — slot totality checked pre-write
  ([§14.6][s14-6]), then boundary zero: project → sweep → events → due `g` updates →
  header + first snapshot ([§14.5][s14-5]).
- `trim!(sim, problem; baseline, t0 = 0.0, backend) → TrimReport` —
  nonlinear least squares on the packed residuals with exact Dual
  Jacobians, against the problem's own `tolerances`
  (`residuals(reads, d) → NamedTuple`, packed in `tolerances`' field order as
  decisions pack in `guess`'s — names pair, order is the declared side's,
  within the problem [§14.7][s14-7] closes at seven fields); setup and commit
  both carry the slot-totality
  check ([§14.6][s14-6]); commit = `init!` with `override(baseline, solution)` —
  boundary zero anchored at `t0`, recordings cleared ([§10.6][s10-6]); resume-at-
  time = `capture`'s returned `t` as `t0`; `converged` = the service's
  per-residual box test at the backend's returned point, backend-independent
  and the commit's gate, with the backend's status and counts recorded
  diagnostically; the backend seam a pinned one-method signature,
  `solve(backend, eval!, d0, lower, upper, tol) → (; d, status, nevals, niters)`
  — in-place `eval!(r, J, d)` filling `J` only when it is not `nothing`,
  packed vectors in the declared orders, `status` an open `Symbol` recorded
  verbatim; non-convergence reports, never
  throws ([§14.7][s14-7], [§14.8][s14-8]).
- `capture(sim) → (condition, t)` — full-store gather including root slots;
  warm restart = capture → tweak → apply ([§14.1][s14-1], [§14.10][s14-10]).
- `linearize(sim, taps) → labeled (ẋ₀, x₀, u₀, y₀, A, B, C, D)` — pure query, one
  seeded Dual pass on scratch; operating point defaults to `capture(sim)`;
  taps = `get_state`/`get_slot`/`get_output` selector lists with control-design
  labels ([§14.10][s14-10]).

**Running.**

- `run!(sim; gui = false, pace = 1, margin = 0.002, t_end = <ctor value>,
  stop_on = <ctor value>)` — `run!` blocks until the run ends; deviceless it is
  fully synchronous on the calling task; `init!` required first
  ([§10.6][s10-6]). Paced and unpaced runs are bit-identical ([§8.7][s8-7]).

  | keyword | default | meaning | owning section |
  |---|---|---|---|
  | `gui` | `false` | **run-scoped attachment**: at run entry it attaches the standard GUI device under the standard greedy binding, with `should_abort = true`, **iff no GUI is already rostered** | [§10.4][s10-4], [§9.6][s9-6], [§9.7][s9-7] |
  | `pace` | `1` | the run's pacing rate | [§8.7][s8-7] |
  | `margin` | `0.002` | the single pacing knob, in seconds | [§8.7][s8-7] |
  | `t_end` | the constructor's value | overrides that default **for this run only** | [§13.5][s13-5] |
  | `stop_on` | the constructor's value | overrides that default **for this run only**, validated against the `Build` here exactly as at construction | [§13.5][s13-5] |

  The GUI is an ordinary rostered device rendered on the calling task
  ([§9.6][s9-6], [§9.7][s9-7]). Because the flag attaches only if no GUI is
  already rostered, a hand-attached GUI makes it a no-op rather than an
  admission error; and because the run's shutdown tail detaches that GUI again
  ([§10.4][s10-4]), the error path included, nothing the flag
  did survives the run; a persistent GUI session is spelled `attach!`/`detach!`
  by hand. Placement follows the roster, not the flag: a rostered GUI moves the
  loop to a spawned task for as long as it is rostered ([§9.1][s9-1],
  [§10.6][s10-6]); sugar never activates by default.

  `margin` defaults to 2 ms, the sleep primitive's granularity plus its
  measured overshoot, with `0` / 2 ms / `∞` spanning the design space
  ([§8.7][s8-7]). The effective `t_end`/`stop_on` pair is recorded in the run
  metadata ([§13.5][s13-5]).
- `step!(sim; frames = 1) → frames_advanced` — synchronous partial advance
  through the ordinary frame sequence, bit-identical to the same frames under
  `run!`; `t_plus = <duration>` is the mutually-exclusive duration spelling
  (whole frames until the boundary time covers that duration); returns the
  frames *actually* advanced, fewer than requested when
  `t_end` or a `stop_on` face ended the run inside the call. Between calls the
  simulation reports `initialized`; `run!` may follow and continues from the
  current boundary; a stepping session is deviceless — write via `stage!`,
  read via `latest` ([§10.6][s10-6]).
- `stage!(sim, "face" => value, ...)` — task-free staging from the calling
  task into the harness register ([§9.3][s9-3]; surface = the
  currently-unclaimed faces): traced, drained last at the next frame top,
  surface-checked exactly
  as the GUI's writes (the harness cell, [§10.6][s10-6]; legal under `run!` and
  `step!` alike).
- `latest(sim) → snapshot` — the current published snapshot, the same
  immutable value device handles read ([§9.2][s9-2]); the assertion/inspection
  accessor of the harness and REPL registers ([§10.6][s10-6]).
- `phase_bodies(sim) → named callables` — the compiled phase bodies of the
  nominal activation, bound over the simulation's own buffers: the four
  blocks (`rhs`, `sweep_hx`, `sweep_hxu`, `ticks` — the sweeps in both
  arities, zero-arg interior and tick-indexed boundary; `ticks` takes the tick
  index) plus per-event guards/handlers and per-component `project`, keyed
  by the model's roster. The [§7.5][s7-5] allocation seam: warm, then
  `@ballocated(body()) == 0` per body; diagnostic register, the one promise
  being identity with what the loop runs; isolated invocation leaves buffers
  valid but off-trajectory — re-run `init!` to continue ([§12.7][s12-7]).
- Control plane — pause/un-pause, pace and `margin` changes, stop on a
  separate atomic surface, never staged ([§10.1][s10-1]; pacing sits outside the
  semantics, so pace and `margin` are both safe to change live).
- Termination — model state via `stop_on` faces read at every published
  boundary ([§13.5][s13-5]); shutdown completes a boundary, publishes the final
  snapshot, then joins ([§10.4][s10-4]).
- Post-run — the log is retained snapshots; `trace(sim) → trc` retrieves the
  always-on input trace, and `replay!(sim2, trc; to_boundary = k)` re-drives
  a fresh `Simulation(world)` bit-identically through the ordinary loop
  (boundary zero from the trace header, drain fed by frame ordinal), ending
  `initialized` — inspect via `latest`/live stores, advance via `step!`,
  continue via `run!`; the state-trajectory inspector and the `StepError`
  reproduction tool ([§9.2][s9-2], [§9.5][s9-5], [§10.7][s10-7]; on-disk persistence deferred,
  [§16][s16]).

---

## Appendix C. The diagnostic kind set

The kinds below are the closed set row 58 commits to, made normative. **Tests
match on kind plus payload fields, never on message text** ([§13.2][s13-2]). So
the tables below — not any message — are the acceptance-test contract. Adding a
kind is a decision-log entry. Every row's payload is *in addition to* what
[§13.2][s13-2] requires of all diagnostics: paths and names as strings, never
instances; the list-in-hand wherever a did-you-mean renders; the didactic
register (state the fix). Owning sections stay the normative home of each rule;
this appendix is an index of the values, in the manner of Appendices A and B.

Severities, in the vocabulary [§13][s13] fixes:

- **build (collected)** — a diagnostic from a declarative pass, collected with its
  siblings and thrown as one `BuildError` at the stratum barrier ([§13.1][s13-1]);
- **build (fail-fast)** — raised while *user code* runs (a boundary-connection body, a
  probe); the first one aborts the phase ([§13.1][s13-1]);
- **service** — raised by a stopped-sim service, or by
  `attach!`/`Simulation`/`run!` validating against the `Build`; collected into one
  carrier wherever the owning section says so (the register, [§14.1][s14-1]; the
  pre-write check, [§14.6][s14-6]), a single throw at the call otherwise;
- **runtime** — fail-fast during a boundary, reaching the single catch site ([§13.4][s13-4])
  as a species of `StepError`;
- **warning (runtime)** — the per-occurrence runtime stream of [§13.2][s13-2],
  carried by the per-writer diagnostic cells ([§9.8][s9-8]) and rate-limited by them:
  every kind in this severity is bounded per writer per boundary (a ring of
  sixteen retained values, the excess becoming per-kind suppressed counts).
  The per-row qualifiers below record where that bound is load-bearing — a
  source that can repeat within a frame — and where the source itself fires
  once.
- **warning (service)** — raised by a stopped-sim service call that
  *completed*: emitted at the call site through the standard logging backend,
  beside the returned value, never thrown, part of no collection; no rate
  limit — each kind fires at most once per call, and its payload is drawn from
  the report the call returns ([§14.5][s14-5], [§14.8][s14-8]).
- **warning (build)** — this severity exists and its set is currently empty
  (row 84).

**Declaration and wiring** (Stratum A):

| kind | payload | owner | severity |
|---|---|---|---|
| `UnknownPort` | the wire end (`source`/`destination`), that end's path, the unknown port name, that end's port list (did-you-mean) | [§6.1][s6-1], [§11.4][s11-4] w1 | build (collected) |
| `UnconnectedInput` | leaf path, input name, declared entry type, the obligation chain's last level | [§6.1][s6-1], [§11.4][s11-4] w2 | build (collected) |
| `TwoProducers` | destination terminal, both producer terminals with provenance (sibling wire / ancestor deep route / boundary connection entry) | [§6.1][s6-1], [§11.8][s11-8] | build (collected) |
| `WireTypeMismatch` | both endpoint paths, both face names, declared entry type, producer face type | [§6.1][s6-1], [§11.2][s11-2], [§11.4][s11-4] w4 | build (collected) |
| `WalkingFaceAtFrozenEntry` | consumer path and entry name, producer path and face name, the offending leaf, both declared leaf types; both remedies in the message ("declare the entry `T` if the consumer promotes; feed it from a non-walking source if the freeze is genuine") | [§6.1][s6-1], [§11.2][s11-2] | build (collected) |
| `PathResolution` | path, offending segment, sibling field list; for a traversal past a generically-held field, that field's declared type | [§6.1][s6-1], [§13.3][s13-3] | build (collected) |
| `AbstractAtRoot` | face name, consuming leaf path, the abstract entry; remedy hint (wire a concrete producer — in a rig, a stub child, [§13.7][s13-7]) | [§11.2][s11-2] | build (collected) |
| `RootSlotTypeConflict` | face name, the consuming paths, their conflicting concrete declarations at nominal (a tolerance difference is not a conflict — the meet, [§11.2][s11-2]) | [§11.2][s11-2] | build (collected) |
| `IllegalStateLeaf` | component path, `init_x` field name, leaf type, the closed vocabulary (scalar / `SArray` at the common eltype) | [§7.1][s7-1], [§11.2][s11-2] | build (collected) |
| `StoreWithoutUpdate` | component path, the `init_x` store, the missing update (neither `f` nor `g` has a method); shadowing note when the parent module defines its own `f`/`g` ([§11.1][s11-1]) | [§11.2][s11-2] | build (collected) |
| `EventHalfMissing` | component path, event name, which half, the function that has no method | [§11.2][s11-2] | build (collected) |
| `PrimitiveAtRoot` | root path, component type | [§11.2][s11-2] | build (collected) |
| `ClassUnreadable` | component path, type, declarations found, both family lists; did-you-mean when the type holds component-typed fields; shadowing note when the parent module defines same-named declaration functions ([§11.1][s11-1]) | [§11.5][s11-5] | build (collected) |
| `ClassMixed` | component path, the `child_connections` declaration and the offending leaf declarations | [§11.5][s11-5] | build (collected) |
| `ContainerMixed` | container field path, offending element keys/indices, their types | [§11.5][s11-5] | build (collected) |
| `DeclarationOnWrongTier` | component path, the offending declaration (`f`/`g`, `events`, `init_m`, or a `workspace`/`output_types` arity — stage names carry no tier since row 173), the tier the leaf's other declarations announce | [§5.2][s5-2], [§11.2][s11-2], [§11.5][s11-5] | build (collected) |
| `TierSignatureMismatch` | component path, the declaration at fault (`input_types` or `output_types`), the leaf's tier, the signature form found versus the form mandated (two-argument `(::C, ::Type{T})` on the continuous tier, plain `(::C)` on the discrete); stateful leaves only — on a stateless leaf `output_types`' arity *is* the tier ([§11.2][s11-2]), so there is nothing to mismatch | [§11.2][s11-2], [§11.5][s11-5] | build (collected) |
| `FaceNameIllegal` | assembly path, face name, the violated invariant (contains `/`) | [§11.6][s11-6] | build (collected) |
| `FaceNameCollision` | assembly path, face name, both entries' provenance (hand-written / computed) | [§11.6][s11-6] | build (collected) |
| `FaceDirectionConflict` | assembly path, the declaring method, the offending entry, the resolved port's actual direction | [§11.6][s11-6] | build (collected) |
| `UnknownFaceSelection` | child path, reason (unknown names / both `except` and `only` given), the offending names, the child's face list | [§11.8][s11-8] | build (collected) |
| `RatesViolation` | assembly path, offending key, reason (deep key / unknown child / `K` on a continuous child) | [§8.5][s8-5], [§11.7][s11-7] | build (collected) |
| `MissingProbeValue` | face name, type | [§12.3][s12-3] | build (collected) |

**Schedule and contract conformance** (Strata B and C):

| kind | payload | owner | severity |
|---|---|---|---|
| `AlgebraicCycle` | the SCC's member terminals in slash form, the wires among them, optional classification (`real`/`artificial`) with the member whose hop died | [§5.5][s5-5], [§5.6][s5-6] | build (collected) |
| `ProducedByTwoStages` | component path, port name, both stage names | [§4.3][s4-3], [§11.3][s11-3] | build (collected) |
| `DeclaredNotProduced` | component path, declared name, the stage-product list and the state-field list | [§11.3][s11-3] | build (collected) |
| `UndeclaredReturnField` | component path, stage, returned field name, candidates (`output_types`) | [§11.3][s11-3], [§11.4][s11-4] w5 | build (fail-fast) |
| `DeadStage` | component path, stage — a stage method returning bare `(;)`, producing neither ports nor `w` | [§5.2][s5-2], [§12.3][s12-3] | build (fail-fast) at probe |
| `ConformanceFailure` | component path, function, field-level diff (missing / unexpected / per-field expected-vs-observed — order-insensitive, the return having been canonicalized first), simulation time | [§12.5][s12-5] | build (fail-fast) at probe; **runtime** as a `StepError` species |
| `GuardForm` | component path, event name, observed probe return type, both admissible forms | [§12.5][s12-5] | build (fail-fast) |
| `BundleFieldError` | component path, function family, requested field, the legal field set, classification (undeclared store / wrong-tier fact / illegal for this function family) | [§5.2][s5-2], [§13.2][s13-2] | build (fail-fast) at probe; **runtime** thereafter |
| `HandlerReturnKey` | component path, event name, offending key, the legal set `{x, m}` narrowed to the stores that exist | [§5.2][s5-2], [§12.5][s12-5] | build (fail-fast) |
| `UserCodeFraming` | component path, which function, the probe context including synthesized inputs; the original exception as `cause` | [§13.2][s13-2] | build (fail-fast) |

**Deployment, periphery and services:**

| kind | payload | owner | severity |
|---|---|---|---|
| `MissingInit` | the simulation's status, the entry point called (`run!`/`step!`) | [§10.6][s10-6] | service |
| `ServiceLifecycle` | the operation (`attach!`/`detach!`/`init!`/`trim!`/`capture`/`linearize`), the current status, the legal statuses | [§9.3][s9-3], [§14][s14] | service |
| `StopFaceInvalid` | face name, reason (unknown / not root-exported / not `Bool`), the root output-face list; the binding site (constructor or `run!`) | [§13.5][s13-5] | service |
| `DeploymentInvalid` | the deployment parameter (`h`, `n`, `Δt_base`, algorithm, `localization_tol`, `localization_budget`, `firing_budget`, the harmonic-grid relation, a non-dividing anchor period or offset — the anchor named with its declaring scope and key), the value in hand, the violated constraint | [§12.1][s12-1] | service (collected) |
| `AttachUnknownFace` | device id, binding entry, face name, the root input-face list | [§9.3][s9-3] | service |
| `AlreadyAttached` | the device id of the existing roster entry, its binding | [§9.3][s9-3] | service |
| `CallerTaskConflict` | both device ids — the rostered `needs_calling_task` holder and the candidate | [§9.1][s9-1], [§9.3][s9-3] | service |
| `ClaimConflict` | face name, claiming device id, incumbent device id | [§9.3][s9-3] | service |
| `EmptyGreedyClaim` | the greedy device's id and its binding — the computed complement was empty, every root input face being claimed already | [§9.3][s9-3], [§9.6][s9-6] | warning (service) |
| `BindingContractMismatch` | the binding type, the trait and the method at fault, and the direction: a declared side whose enumeration method is missing (`is_input`/`is_output` true, the root's error fallback reached), or a `claims`/`reads` method defined under a false trait (detected by `which` against the fallback); `is_greedy` without `is_input`, and `claims` defined on a greedy binding, report here too | [§9.6][s9-6] | service |
| `ReadBindingUnresolved` | device id, the selector, path and field, candidates; a `reason` distinguishing an unresolved path from a store selector in a snapshot binding (the source rule, [§14.4][s14-4]) | [§9.2][s9-2], [§14.4][s14-4] | service |
| `ConditionResolution` | entry path, store, field, offending value type and declared leaf type, provenance chain; sub-kinds: unknown path, undeclared field, unconvertible value, unexported slot face | [§14.2][s14-2], [§14.3][s14-3] | service (collected) |
| `DuplicateConditionLeaf` | the leaf `(path, store, field)`, both provenance chains, the `override` advice | [§14.2][s14-2] | service (collected) |
| `ConditionNodeMisuse` | the offending argument's type, the node kinds in hand | [§14.2][s14-2] | service |
| `UninitializedSlots` | every uncovered root face, in declaration order | [§14.6][s14-6] | service (collected), pre-write |
| `TapResolution` | tap set (`x`/`u`/`y`), selector kind, path, field, optional index, candidates; for a declaredly-unseedable slot, the pinning consumer's path and its `input_types` entry | [§14.10][s14-10] | service (collected) |
| `TrimProblemInvalid` | the offending `TrimProblem` field, the names or types in hand (a key-set or field-type mismatch; never a field-order difference) | [§14.7][s14-7], [§14.8][s14-8] | service (collected) |
| `TrimCommitEvents` | the events fired at boundary zero: component paths and event names; the same list rides the `TrimReport` | [§14.5][s14-5], [§14.8][s14-8] | warning (service) |
| `TrimCommitResiduals` | the offending residual names with committed-state values and tolerances — a converged solve whose committed-state residuals violate the box test | [§14.5][s14-5], [§14.8][s14-8] | warning (service) |
| `GridUtilization` | the derived `Δt_base`, its driver entries with provenance and refinement factors, and `min_i Dᵢ` — the grid rendered as "N× finer than the fastest declared work" | [§12.1][s12-1], [§12.2][s12-2] | warning (service), at deployment binding (derivation path only) |
| `ReplayHeaderMismatch` | the mismatch, discriminated: a store or slot (component path, store, expected vs. found layout/type) or a deployment parameter (`Δt_base`/`h`/`n`/algorithm/`localization_tol`/`localization_budget`/`firing_budget`, recorded vs. bound value); the build's and the trace's provenance | [§9.5][s9-5], [§10.7][s10-7] | service |
| `ReplayUnknownFace` | face name, frame ordinal, the trace's device tag, the root input-face list | [§10.7][s10-7] | service (collected) |

**Runtime:**

| kind | payload | owner | severity |
|---|---|---|---|
| `StepError` | the carrier: cursor frame (component path, function, boundary phase — RK stage, event round, localization probe, tick), boundary time, frame-entry boundary index (replay pointer), original exception as `cause` | [§13.4][s13-4] | runtime |
| `NonfiniteState` | component path, the offending state block, boundary time and index | [§13.4][s13-4] | runtime |
| `ChatteringBudget` | component path, event name, boundary time, the exhausted `localization_budget` and the frame's localization count | [§8.4][s8-4] | warning (runtime) |
| `FiringBudget` | component path, event name, boundary time, the exhausted `firing_budget` and the boundary's firing count | [§8.6][s8-6] | warning (runtime) |
| `DebtReanchor` | forgiven debt, the new schedule anchor, boundary time | [§8.7][s8-7] | warning (runtime) |
| `ClaimedFaceEntry` | face name, the incumbent (claiming) device id, the discarded value; the site (staging, or a stopped-sim attach's renormalization). Harness-register only — a device's out-of-surface entry is `OutOfClaimEntry` | [§9.3][s9-3], [§9.4][s9-4] | warning (runtime) |
| `OutOfClaimEntry` | device id, face name, the discarded value, the device's claim set; the incumbent's device id when the face is claimed elsewhere | [§9.3][s9-3] | warning (runtime) |
| `ThreadBudget` | thread count, device-task count | [§10.2][s10-2] | warning (runtime), at `run!` |
| `DeviceJoinTimeout` | device id, the join timeout, boundary time and index at shutdown | [§10.4][s10-4] | warning (runtime) |
| `DeviceCrash` | device id, the original exception as `cause`, whether `should_abort` was set; also the init-time failure, reported pre-spawn from the initialization bracket after its `shutdown!` | [§10.4][s10-4], [§9.6][s9-6], [§13.4][s13-4] | warning (runtime) |
| `ReplayDiscardedStaging` | device id, the discarded batch's face names, frame ordinal | [§10.7][s10-7] | warning (runtime), repeating source — rate-limited per writer ([§9.8][s9-8]) |
| `MalformedDatum` | device id, the cause exception; emitted by the author's loop body via `report!(handle, …)` | [§9.6][s9-6], [§13.4][s13-4] | warning (runtime), repeating source — rate-limited per writer ([§9.8][s9-8]) |
| `EntryTypeMismatch` | writer id, face name, the offending value's type, the slot's declared type, the discarded value | [§9.4][s9-4] | warning (runtime) |
| `UnboundedRun` | the effective `t_end`, `stop_on` set and `pace`; the remedy names both, and — interactively — the operator interrupt as the sanctioned escape from the configuration warned about ([§10.4][s10-4]) | [Appendix B][sB], [§13.5][s13-5] | warning (runtime), at run start |

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

<a id="g-abstract-entry"></a>**abstract entry** — an `input_types` entry whose declared type is abstract,
stating **structural substitutability**: any concrete producer face below the
bound wires to it (the field handles, [§4.4][s4-4], are the demonstrated client). Never
needed for eltype genericity, and illegal where the face surfaces as a root
slot (`AbstractAtRoot`) ([§11.2][s11-2]).

<a id="g-assembly"></a>**assembly** — pure composition: component-typed fields as children, plus
`child_connections` (mandatory, the class marker), `input_connections`,
`output_connections` and `sample_times`, with no
dynamics of its own; flattened away for scheduling, retained as the
navigation hierarchy and as declaration-level rate scopes ([§3.3][s3-3], [§11.5][s11-5]).

<a id="g-auto-published-port"></a>**auto-published port** — a declared output that matches a state or mode field
by name and type and that no stage produces: the framework publishes it from
the store at stage-1 position (`h_x`, either tier; the match is against
`init_x`, plus `init_m` on the continuous tier). Contract-driven — row 16
rejected blanket identity publication of state — a framework write, never a
probe product, and excluded from the stage-1 hand-down `y_x` ([§5.3][s5-3], [§11.3][s11-3],
[§5.2][s5-2], row 169).

<a id="g-class"></a>**class** — a component's primitive-vs-assembly status, read off *which*
well-known declarations its type defines: `child_connections` ⇒ assembly, any leaf
declaration ⇒ primitive, neither ⇒ `ClassUnreadable` ([§11.5][s11-5]). Not to be
confused with *tier* (continuous vs. discrete, [§D.4][sD-4]) or with a diagnostic
*kind* ([§D.9][sD-9]).

<a id="g-component"></a>**component** — the unit of modeling: a leaf (continuous or periodic discrete
primitive) or an assembly of components; "primitive" and "leaf" are used
interchangeably for the non-assembly classes ([§3][s3]).

<a id="g-container-children"></a>**container children** — a `Tuple`/`NamedTuple` field whose elements are all
components, contributing them as children path-named `"field/1"` or
`"field/key"`. Transparent grouping, not an assembly: no contract, no
`child_connections`, no rate scope ([§11.5][s11-5]).

<a id="g-continuous-component"></a>**continuous component** — the hybrid primitive: continuous state `x`, modes
`m`, flow `f`, two output stages, events (guards + handlers) and optional
`project`; any facet may be empty, so a state-free instance is an FSM ([§3.1][s3-1]).

<a id="g-contract"></a>**contract** — a component's declared interface: `input_types` (its
requirements, read permissively — what each entry *allows* to arrive) and
`output_types` (its public ports, read literally — what each cell *carries*).
Both take the two-argument `T`-form on the continuous tier and the plain form on
the discrete one. Declared in
`output_types` = public, returned in `w` = private by construction, returned in
`y` and declared nowhere = build error ([§11.2][s11-2], [§11.3][s11-3]).

<a id="g-declaration-inventory"></a>**declaration inventory** — the closed set of well-known functions a component
or assembly defines — `init_x`/`init_m`, `workspace`,
`input_types`/`output_types`, `events`, the stages, `f`/`g`/
`project`, and `child_connections`/`input_connections`/`output_connections`/`sample_times` — each declared in a stated
register of authority: by value, by type, by allocation ([§11.2][s11-2]).

<a id="g-function-family"></a>**function family** — which bundle fields a given function may legally
receive: `h_x`/`h_xu`/`f`/`g`/guard/handler/`project` (the `h_*` sets being
per-tier, [§5.2][s5-2]), with the comment block ([§5.2][s5-2]) stating each
family's maximal legal set and
`BundleFieldError` classifying a read as illegal for the family ([§5.2][s5-2]).
Not a diagnostic *kind* ([§D.9][sD-9]).

<a id="g-generic-holding"></a>**generic holding** — a parent holding a child through a non-concrete field
type; the child is opaque below its faces, and the wires and boundary connections
referencing those faces *are* the imposed contract, checked per instantiation
([§11.8][s11-8], [§6.1][s6-1]).

<a id="g-hybrid-causal-system"></a>**hybrid causal system** — what the framework simulates: continuous flow with
algebraic outputs, multi-rate periodic discrete dynamics, zero-crossing
events, post-step manifold projection, and externally injected inputs ([§2][s2]).

<a id="g-the-letters"></a>**the letters** — `f` the continuous flow, `g` the discrete update, `h_*` the
output stages (suffix = dependence class, not argument list), `x` the state on
either tier and `m` the continuous-only mode store, `u` wired inputs, `y` own
published signals, `ws` the
workspace. Bare `h` means the integration step size only ([§5.3][s5-3], [§8][s8]);
bare `z` means only the shift operator `z⁻¹` (row 173).

<a id="g-periodic-discrete-component"></a>**periodic discrete component** — a leaf with state `x`, update `g`
at a declared rate, and two output stages whose cells hold zero-order between
ticks; it has no `m` store, and its `x` reaches others only through signals
([§3.2][s3-2]).

<a id="g-rate-scope"></a>**rate scope** — an assembly's `sample_times` declaration: immediate child
name ⇒ `Relative` or `Absolute` declaration against the enclosing scope,
relative entries composing affinely down the tree, absolute entries
anchoring; all compiled to one `(D, Φ)` pair per discrete component
([§11.7][s11-7], [§8.5][s8-5]).

<a id="g-schema-authority"></a>**schema authority** — the principle that declarations *define* structure and
evaluation only *checks* conformance against them, never the reverse; types by
declaration, values by execution, conformance by comparison ([§11.1][s11-1]).

<a id="g-stage-function"></a>**stage function / two-stage outputs** — every component provides exactly two
output stages, `h_x` (no `u` in the bundle, hence structurally no
feedthrough) and `h_xu`; feedthrough is thereby declared by signature,
with no dependency annotations anywhere ([§5.2][s5-2]).

<a id="g-workspace"></a>**workspace** — component-declared mutable scratch, declared *by allocation*
(`workspace(::C, ::Type{T})` continuous, `workspace(::C)` discrete), arriving
as the `ws` bundle field; excluded from state semantics, never a condition
target, and never inspected or mutated by the framework — contents at call
entry are unspecified ([§7.3][s7-3]).

### D.2 Signals and data homes

<a id="g-buffer"></a>**buffer** — the framework-owned contiguous `Vector{T}` backing all continuous
state, laid out at build time; authoritative, with typed state values as
ephemeral reconstructions of it ([§7.1][s7-1]). The integration intermediates ([§13.6][s13-6]) live
in framework-owned integrator buffers, never in a component's workspace.

<a id="g-bundle"></a>**bundle** — the single `NamedTuple` of zero-copy views a component function
receives beside the component itself. Under the bundle law a name is present
**iff** the corresponding store or fact exists for that component; undeclared
stores are absent, never `nothing`-filled ([§5.2][s5-2]).

<a id="g-cell"></a>**cell** — one concretely-typed entry of the signal table, one per output port
of the flattened model, written by its producing
stage and read by every gatherer ([§4.1][s4-1]). Every cell is public, private
intermediates never being cells ([§11.3][s11-3]). Bare "cell" is only this — see
*staging cell* ([§D.6][sD-6]) and *store*.

<a id="g-constant-source"></a>**constant source** — an ordinary library component with no inputs and no
state, publishing a value its instance holds (`Constant{V}`); the spelling for
an aggregate input with zero contributors and for the rig stub feeding an
abstract face. Its value is instance data, never a default ([§13.7][s13-7], [§6.2][s6-2]).

<a id="g-entry"></a>**entry** — never used bare: the spec's compounds are table entry (a cell),
input entry ([§11.2][s11-2]), executor entry ([§12.7][s12-7]), roster entry ([§9.3][s9-3]), batch entry
([§9.3][s9-3]) and condition entry ([§14.3][s14-3]), each a different thing.

<a id="g-face"></a>**face** — the name a port wears on its component's boundary: for a leaf the
port's own name, for an assembly a name declared in `input_connections` or
`output_connections`, aliasing an interior port. An opaque token with
two build-checked invariants (no `/`, unique within the assembly), its type
derived from its internal endpoint and its direction declared by the method
that names it. The periphery's write side
speaks face names only; the read side speaks them wherever it wants contract
rather than structure ([§11.6][s11-6], [§9.2][s9-2], [§14.4][s14-4]).

<a id="g-feedthrough"></a>**feedthrough** — an instantaneous input→output dependence. **Structural
feedthrough** is this design's version: fixed by which stage produces a port
rather than annotated, with every stage-2 output conservatively presumed
dependent on every wired input ([§5.3][s5-3]).

<a id="g-field-handle"></a>**field handle / function-valued signal** — an immutable query object carried
on an ordinary port (`ISAField`, `TerrainField`) that consumers evaluate at
arguments of their own choosing; bulk data rides as build-time-frozen
references, never as mutable caches ([§4.4][s4-4]).

<a id="g-immutable-value-semantics"></a>**immutable value semantics** — the signal rule, stated precisely as
immutability *plus frozen references* (`isbits` is the common case, not the
rule): no aliasing, safe concurrent reads, and a definite per-cell freshness
tied to the producer's schedule position ([§4.1][s4-1]).

<a id="g-one-home-per-datum"></a>**one home per datum** — buffer for continuous `x`, stores for discrete `x` and for `m`, table for
produced signals; no store mirrors another, and the table never holds
transported data ([§5.2][s5-2], [§7.1][s7-1]).

<a id="g-port"></a>**port** — the addressable unit of the model: one declared name, one cell, one
root slot, one staged write, one device claim, one trace address, one GUI
liveness verdict. Wiring is port-granular, and which stage computes a port is
invisible outside the component ([§4.2][s4-2], [§4.3][s4-3]).

<a id="g-signal-table"></a>**signal table** — the framework-owned collection of cells holding every
produced signal of the flattened model; consumers gather views from it, and
its consistency is a boundary property, transiently integrator scratch within
a step ([§4.1][s4-1], [§8.3][s8-3]).

<a id="g-slot"></a>**slot** — a root input slot: the root assembly's own input face,
produced by no component, constant within a frame, and the only thing the
periphery may write ([§9.3][s9-3], [§11.6][s11-6]).

<a id="g-staging-cell"></a>**staging cell** — the per-device atomic holding place where a device's pending
write batch waits between drains; mutated frame by frame, hence outside the
table's publish-once discipline ([§9.4][s9-4]). Not a table cell ([§4.1][s4-1]).

<a id="g-store"></a>**store** — the typed home of `m` and of a discrete leaf's `x`: overwritten by the framework when a
handler or update returns a new value, never arithmetic-touched, snapshot-free
to copy. Never called a cell — root slots, by contrast, *are* source cells of
the table ([§7.3][s7-3], [§4.1][s4-1], [§9.2][s9-2]).

<a id="g-summing-junction"></a>**summing junction** — an ordinary library component performing N-to-1
aggregation through explicit wires (`SumJunction{W, N}` or a named
site-specific variant); there is no framework aggregation mechanism, and fold
order is the junction's positional input order ([§6.2][s6-2]).

<a id="g-value-level-constructor"></a>**value-level constructor** — the plain exported function (component, input
values) → field handle that every field-emitting component is obliged to
provide, its own swept output stage being a one-line call to it; the device by
which [§14.1][s14-1] condition math queries the environment before any sweep exists
([§4.4][s4-4]).

<a id="g-view"></a>**view** — a zero-copy reconstruction of a store handed to a function through
its bundle; it materializes in the caller's frame for the duration of the call
and is value-identical on re-materialization within a sweep ([§7.1][s7-1], [§5.2][s5-2]).

<a id="g-w"></a>**`w` (private intermediates)** — the optional second slot of an output stage's
return, `(y, w)`: an `isbits`-leaf NamedTuple of values the component computes
for its own later use. Private by construction — no cell, no contract entry,
nothing to wire, list, filter or address — travelling exactly one hop by the law
([§5.2][s5-2]), SSA-passed within a fused pass, probe-observed in type and
checked at the nominal activation only ([§11.3][s11-3], [§12.5][s12-5]). The inspection path
for one is promotion to a declared output.

### D.3 Evaluation and scheduling

<a id="g-algebraic-loop"></a>**algebraic loop** — a genuine cycle in the instantaneous dependency graph: a
build error naming the path in canonical slash form, broken by the author with
inserted dynamics, an explicit `UnitDelay` or restructuring ([§5.5][s5-5]). Not to be
confused with an **artificial loop**, port-level acyclic but stage-level
cyclic, whose remedy is a ladder — the two-stage split, contract re-factoring,
and as residual a component split ([§5.4][s5-4]).

<a id="g-flow"></a>**flow / RHS** — `f`, the continuous derivative function. Evaluating the RHS
means running the whole sweep, since `f` reads the fresh table: there is no
incremental `f`-only re-evaluation ([§3.1][s3-1], [§5.3][s5-3]).

<a id="g-frame"></a>**frame** — one iteration of the loop — drain, integrate, boundary sequence,
publication — the unit `step!` counts and the trace's ordinal key ([§9.1][s9-1]).
Distinct from a *boundary* ([§D.4][sD-4]), and from the kinematic reference frames of
the aircraft domain, which always appear compounded ("the b frame").

<a id="g-projection"></a>**projection** — the optional per-component hook `x ← project(x)`, run in the
only two schedule positions between a state write and its decode (after
integration, after a handler's `x`-reset); the cheap end of geometric
integration's projection methods ([§2][s2], [§5.3][s5-3]).

<a id="g-schedule"></a>**schedule** — the static evaluation order computed once at build time from
wiring edges plus intra-component feedthrough: all stage-1 functions in any
order, stage 2 in topological order, then `f`. The hot loop runs a flat list of
`(component, stage)` entries, with zero runtime graph logic ([§5.1][s5-1]).

<a id="g-sweep"></a>**sweep** — one execution of that schedule against the current state, in one of
two statically distinct variants compiled from the same entry list: the
**interior sweep** walks continuous entries only — what RK stage evaluations and
localization guard probes run, so discrete cells hold ZOH mid-step by
construction — and the **boundary sweep** walks the full list, with the
boundary's due discrete entries gated in by counter modulo. Mid-step sweeps are
integrator scratch; the boundary sweep restores table consistency, and the event
phase re-runs whole boundary sweeps, against that boundary's fixed due set,
until quiescence ([§5.3][s5-3], [§8.3][s8-3], [§8.5][s8-5], [§8.6][s8-6]).

### D.4 Time and events

<a id="g-anchor"></a>**anchor** — the exact `(T, τ)` pair an `Absolute` entry establishes: period
and offset in rational seconds, severing its subtree from the enclosing
scope's grid; anchor 0 is the base grid itself. Anchors join the deployment
constraint pool; relative declarations below one compose against it exactly
as against the root grid ([§8.5][s8-5], [§12.1][s12-1]).

<a id="g-bound-schedule"></a>**bound schedule** — the named printable artifact on the `Simulation`
produced by deployment binding: per discrete component, `(D, Φ, Δt)` with
anchor and provenance columns — the single source of truth for `Δt` and the
substrate of the grid diagnostics and the hyperperiod chart ([§12.2][s12-2], [§8.5][s8-5]).

<a id="g-boundary"></a>**boundary** — a published consistency point: where the [§8.6][s8-6] macro-sequence
completes and a snapshot goes out. Every grid point is a boundary, but `t*`
and boundary zero are boundaries that are not frame tops ([§8.4][s8-4]). *Boundary
zero* ([§14.5][s14-5], [§D.8][sD-8]) is a hyponym — it is an ordinary boundary whose incoming
transitions are authored rather than computed.

<a id="g-boundary-detected"></a>**boundary-detected** — the detection policy a `Bool`-returning guard declares:
guards are checked for
not-holding → holding edges against their priors at step boundaries only, with
no root-finding and no step rejection, the handler firing at the end of the
step in which the edge was observed. Exact, not approximate, for guards over
`u`/`m` alone — those predicates are piecewise frame-constant ([§8.4][s8-4]).

<a id="g-chattering"></a>**chattering / localization budget** — the bounded per-frame localization
allowance, `localization_budget`, a `Simulation` deployment keyword defaulting
to 8; exhaustion *degrades* rather than throws — localization stops for the
rest of the frame and further crossings fire at the next boundary, under a
`ChatteringBudget` warning naming the event ([§8.4][s8-4]).

<a id="g-dt_base"></a>**`Δt_base`** — the base tick period, an integer multiple `n·h` of the
continuous step, bound at `Simulation` construction from one of three
sources: explicit keyword, `n·h`, or — fully anchored models only —
derivation from the constraint pool; every discrete component's period is an
integer multiple of it ([§8.5][s8-5], [§12.1][s12-1]).

<a id="g-due"></a>**due** — a discrete component is due at a boundary when its compiled `(D, Φ)`
pair admits that boundary's tick index (`(idx − Φ) % D == 0`); due components'
output stages are gated into the *boundary* sweep (never the interior one) and
their `g` updates run after quiescence. The due set is a property of the
boundary, fixed for its whole event iteration: the gate's image of the frame
index at a frame top, empty at `t*`, the `Φ = 0` set at boundary zero
([§8.5][s8-5], [§8.6][s8-6]).

<a id="g-edge-semantics"></a>**edge semantics / holding** — an event fires on a not-holding → holding
transition of its predicate, never on a bare sign change; the opposite
crossing direction is declared as a second event with the negated guard ([§2.1][s2-1],
[§8.6][s8-6]).

<a id="g-firing-budget"></a>**firing budget** — the rule bounding the boundary event iteration: each
declared event fires at most `firing_budget` times per boundary (a
`Simulation` deployment keyword, an integer ≥ 1 defaulting to 4), eligibility
being a not-holding → holding edge on the event's last-observed sample. An
event re-enabled within the boundary therefore fires *at* that boundary;
exhaustion drops its further edges there under a `FiringBudget` warning
([§8.6][s8-6]).

<a id="g-guard"></a>**guard** — the declared function defining an event's predicate, evaluated
against the fresh boundary table and paired with a handler in an ordered,
named `events` collection; its detection policy is declared by its return
type — `Bool` boundary-detected, the nominal scalar localized ([§8.4][s8-4],
[§11.2][s11-2]).

<a id="g-harmonic-grid"></a>**harmonic grid** — the rule that every discrete period is an integer multiple
of `Δt_base` — and every anchor period and offset an integer multiple
likewise — itself an integer multiple of `h`, so ticks land only on step
boundaries; grid times are indexed from the frame count, never accumulated
([§8.5][s8-5], [§8.4][s8-4]).

<a id="g-input-epoch"></a>**input epoch** — a maximal span of constant `u`, delimited by the frame-top
drains ([§9.4][s9-4]). Within an epoch a guard changes only through the
trajectory; at a seam it can jump without crossing. Hence the **θ = 0
validation**: the first act of a triggered localization is a probe at `xₙ`
under the frame's own `u`, whose `σ₀` both supplies the left bracket value and
tells a *trajectory-caused* edge (root-find) from an *epoch-caused* one (no
in-frame crossing exists — discard the localization and let the event fire in
the boundary's ordinary iteration, no budget, no warning) ([§8.4][s8-4]).

<a id="g-interpolant"></a>**interpolant** — the lazily built cubic Hermite continuous extension over the
last completed step, from which localization probes read the states they sweep;
built only after the θ = 0 validation confirms an in-frame crossing, and
invalidated at `t*`, where the handlers have made it a lie ([§8.4][s8-4]).

<a id="g-localized"></a>**localized** — the detection policy a sign-form guard declares: the crossing
instant is
bracketed by derivative-free root-finding over probe sweeps of interpolated
states, to a bracket narrower than `localization_tol · h` (a deployment
keyword, default `1e-6`). Only the sign form can declare it — the `Bool` form
offers no root to bracket — and it runs
identically paced or unpaced ([§8.7][s8-7], [§8.4][s8-4]).

<a id="g-pacing"></a>**pacing / pacer debt** — the pacer inserts waits between completed frames and
never alters the boundary sequence; a frame exceeding its wall budget leaves
**debt** that later frames repay, with excess forgiven by re-anchor plus
warning ([§8.7][s8-7]).

<a id="g-phase"></a>**phase (`Φ`)** — a schedule's offset against its grid: in scope ticks for
`Relative(K, Φ)`, in rational seconds for `Absolute(q, τ)`, compiled to base
ticks with `0 ≤ Φ < D` by construction; the boundary gate is
`(idx − Φ) % D == 0`, and a phase shifts firing instants, never the period
([§8.5][s8-5]).

<a id="g-predicate"></a>**predicate** — what a guard defines: a `Bool`-valued form, or the sign of a
continuous function with positive = holding (writing the sign value `σ`,
holding = `σ ≥ 0`) ([§2.1][s2-1]). Not to be confused with the *condition*
([§14][s14]), the value that sets a build's state ([§D.8][sD-8]).

<a id="g-prior"></a>**prior** — the per-event stored sample of its predicate at the previous
boundary's quiescence, always an honest observation and never a manufactured
one; held in loop state and never in a state store; "newly fired"
is defined against it for the boundary's first round (later rounds test the
last-observed sample), and boundary zero establishes every prior as
not-holding ([§8.6][s8-6]).

<a id="g-quiescence"></a>**quiescence** — the fixed point of the boundary event phase: rounds of
[sweep → guards → handlers] iterate until a round fires nothing, after which
the priors are updated and due `g` updates run ([§8.6][s8-6]).

<a id="g-remainder-step"></a>**remainder step** — the integration from `t*` to the original grid target
after a localized event, with `h′` derived at use; guards are re-checked on it
under the localization budget ([§8.4][s8-4]).

<a id="g-t"></a>**`t*`** — the localized event time: the holding endpoint of the root-finder's
final bracket, structurally strictly later than `tₙ`. A full boundary runs
there, but no ticks are due and no staged inputs are drained ([§8.4][s8-4]).

<a id="g-tick"></a>**tick** — an instant at which a discrete component's stages and update run,
gated by counter modulo against the harmonic grid inside the boundary sweep;
different boundaries therefore run different subsets of the schedule ([§8.5][s8-5]).

<a id="g-tier"></a>**tier** — the continuous or discrete side of the hybrid formalism, read off a
leaf's declaration shape (`DeclarationOnWrongTier` names a violation) ([§11.2][s11-2],
[§11.5][s11-5]). Bare "tier" means only this: the genericity classes are *walked /
pinned / exempt* ([§D.5][sD-5]) and the detection policies *boundary-detected /
localized*.

### D.5 Build pipeline

<a id="g-activation"></a>**activation** — a re-run of Stratum C at a given scalar type `T`: cells
re-typed (producer-fed ones by evaluating the producer's output declaration at
`T`, root slots by evaluating the consuming `input_types` entry at `T`, the
state type by the leaf walk), buffers re-laid-out,
workspace allocators re-invoked, probe chain re-run. Structure and schedule are `T`-independent;
non-nominal activations are lazy, with an opt-in exhaustive set for CI ([§12.4][s12-4]).

<a id="g-always-on-conformance-check"></a>**always-on conformance check** — the probe's comparison left permanently in
place: one type test of a stage return against the declaration-derived
expected `NamedTuple` at the table-write point, preceded by a type-level
reorder to that type's field order (the names pair; order carries no
semantics), both folding to zero instructions
when the return type is proven ([§12.5][s12-5]).

<a id="g-build"></a>**`Build`** — the artifact `build(world)` produces: wire list, face table with
provenance, schedule and root slots as plain printable data — the inspectable
contract of the instantiation, and what `attach!`, `stop_on`, replay and
condition resolution all validate against ([§12.2][s12-2]).

<a id="g-chunking"></a>**chunking** — splitting a large phase body's entry tuple into statically
typed chunks behind non-inlined function barriers; the implementation's only
representation freedom, converting compile cost from superlinear in body size
to linear in entry count ([§12.7][s12-7]).

<a id="g-executable-set"></a>**executable set** — the function set an activation can actually run, hence
exactly what it probes: a `Dual` activation sees only the continuous output
stages and `f` — never the discrete stages, `g`, guards or handlers ([§12.4][s12-4]).

<a id="g-executor"></a>**executor** — the compiled execution form of the schedule: a concretely-typed
tuple of entries over statically typed cell storage, traversed by a
compile-time-unrolled walk, with code-selecting facts in type parameters and
plain data in fields ([§12.7][s12-7]).

<a id="g-leaf-walk"></a>**leaf walk** — the framework's derivation of per-activation types from a
declared nominal type: real leaves and `Real` type parameters follow the
activation scalar, everything else pins. It applies on the **state** side alone
(the type derived from a continuous leaf's `init_x`; `init_m` and a discrete
leaf's `init_x` pin wholesale). **Cells are not
walked**: an output cell comes from evaluating the producer's `output_types` at
the activation scalar (row 166) and a root-slot cell from evaluating the
consuming `input_types` entry at it (row 167), participation and tolerance
authored per leaf in both ([§11.2][s11-2]; applied in Stratum C,
[§12.1][s12-1]).

<a id="g-lens"></a>**lens (`Getter`)** — the compiled navigation step of a condition entry: its
tree position tuple lifted to a type parameter, giving type-stable access to
the authored value at apply time ([§14.3][s14-3]).

<a id="g-measurement-seam"></a>**measurement seam / phase bodies** — `phase_bodies(sim)` returns the compiled
bodies of the nominal activation bound over the simulation's own buffers
(`rhs`, `sweep_hx`, `sweep_hxu` — the sweeps in both arities, zero-arg interior
and tick-indexed boundary — `ticks`, plus per-event guards and handlers
and per-component `project`). Its one promise is identity with what the loop
runs, which is what makes the allocation assertions ([§7.5][s7-5]) honest ([§12.7][s12-7]).

<a id="g-nominal"></a>**nominal** — the `Float64` activation, and of a declaration its `Float64`
face (for a continuous producer's output declaration, its evaluation at
`Float64`); the only activation that runs in real time, and the one where the
conformance check demands exact type match ([§11.2][s11-2], [§12.4][s12-4], [§12.5][s12-5]).

<a id="g-probe"></a>**probe** — the build's single evaluation of a user function with real values,
checking shape and type conformance and discarding the result. Every user
function is probed once, at the initial state; probes see only that state's
branch ([§12.3][s12-3]).

<a id="g-probe-value"></a>**probe value / input synthesis** — `probe_value(::Type)` fabricates values
for the one kind of terminal with no producer, root slots
(`zero(T)`/`false`/first enum/`T()`, overridable). Strictly probe-scoped:
never an initial slot value, which [§14.6][s14-6] makes a structural barrier ([§12.3][s12-3]).

<a id="g-probedual"></a>**`ProbeDual`** — the framework's exported canonical concrete probe scalar
(`ForwardDiff.Dual{ProbeTag, Float64, 1}`), which keys the CI activation
pinning walked-leaf genericity; its width is arbitrary, since what CI pins is
genericity, not a particular Jacobian ([§12.4][s12-4]).

<a id="g-schema-vs-layout"></a>**schema vs. layout** — the two lookup families the `Build` supplies to
condition resolution: *schema* is the evaluated declarations (may you write
this field, at what leaf type — the authority), *layout* is where it
physically lives (buffer ranges, store and slot indices) ([§14.3][s14-3]).

<a id="g-stratum"></a>**stratum** — one of the build's three phases: A structure (pure declaration
reading), B schedule (the single evaluation-feeds-structure step), C
activation (everything type-shaped). Strata are barriers — a stratum that
produced any error-severity diagnostic throws before the next begins ([§12.1][s12-1],
[§13.1][s13-1]).

<a id="g-walked"></a>**walked / pinned / exempt** — the eltype-genericity classes: walked
payload/value types follow the activation scalar, pinned parameters and
definitions stay `Float64`, and the discrete side is exempt. Enforced by the
leaf walk on the state side and stated per leaf in a continuous leaf's contract
declarations on the cell side — `output_types` for what a producer's cells
carry, `input_types` for what a consumer's entries tolerate ([§7.2][s7-2], [§11.2][s11-2]).

### D.6 Runtime periphery

<a id="g-bad-datum"></a>**bad datum** — a datum unmappable for environmental reasons (truncated
datagram, malformed JSON, out-of-range field): tolerated *in the loop body* —
catch, stage nothing, `report!(handle, MalformedDatum(cause))`, continue —
while any other exception propagates and becomes `DeviceCrash`. The
classification is the device author's ([§9.6][s9-6]).

<a id="g-batch"></a>**batch** — a device's staged set of face ⇒ value writes, coalesced in its
staging cell and applied whole at the next drain ([§9.4][s9-4]). The word means only
this; error reporting *collects* ([§D.9][sD-9]).

<a id="g-binding"></a>**binding** — the value passed at `attach!` that makes a device
framework-legible: a subtype of `AbstractBinding` declaring its sides by the
Bool traits `is_input`/`is_output` (false by default on the root), with
`is_greedy` switching the input side's claim source from returned to computed
(the unclaimed complement, in place of `claims`). `claims` and `reads` carry
error fallbacks on the root, and attach cross-checks each trait against its
method in both directions (`BindingContractMismatch`); `map_input`/`map_output`
are loop-idiom conventions the framework never calls. Every input-side binding
stakes a claim; `TableBinding` is the shipped
data-driven one ([§9.6][s9-6], [§9.4][s9-4]).

<a id="g-boundary-counter"></a>**boundary counter** — the monotonic count of *published boundaries* carried
in the snapshot and mirrored in the loop state the wait predicate tests;
incremented after the `latest` release-store, so a waking waiter can never see
a stale snapshot ([§10.3][s10-3]).

<a id="g-calling-task"></a>**calling task** — the task that invoked `run!`. It runs the loop itself (the
unattended register) unless a `needs_calling_task` device is rostered, in
which case it runs that device's loop body inline and the loop moves to a
spawned task ([§9.1][s9-1]).

<a id="g-claim"></a>**claim** — the set of faces a device *may* write, registered at attach —
either returned by its binding's `claims` or computed as the unclaimed
complement under `is_greedy` — and released at detach; claiming an
already-claimed face is an attach-time error (`ClaimConflict`), and a broad
claim costs GUI liveness ([§9.3][s9-3]).

<a id="g-coalescing"></a>**coalescing** — the CAS merge keeping one pending batch per device:
untouched faces survive, re-staged faces take the newest level (the per-face
ZOH). Its outbound mirror is newest-wins snapshot delivery ([§9.4][s9-4], [§10.3][s10-3]).

<a id="g-control-plane"></a>**control plane** — the separate few-word atomic surface carrying pause,
un-pause, pace, `margin` and stop, consulted at frame top and inside the
loop's wait and pause states; structurally not staging, since a paused loop
drains nothing ([§10.1][s10-1]).

<a id="g-derived-liveness"></a>**derived liveness** — the rule that a GUI widget is live iff its port's feed
chain terminates in a root slot inside the GUI's own claim in the run's frozen
surface partition; baked once at run start, with no per-port "GUI-controlled"
marking anywhere ([§9.7][s9-7]).

<a id="g-device"></a>**device** — any attached participant in the periphery: a subtype of
`AbstractDevice` under one authoring
contract (`init!`/`loop`/`shutdown!`, optional `unblock!` and
`needs_calling_task`) and one handle; input-only and output-only are
degenerate uses, and the GUI is an ordinary device ([§9.6][s9-6]).

<a id="g-diagnostic-cell"></a>**diagnostic cell** — the single-writer cell each rostered device and the loop
itself owns for runtime diagnostics and liveness: a bounded ring (capacity 16)
of diagnostic values plus per-kind suppressed counts — the bound being the
rate limit itself — and an atomic heartbeat timestamp, taken by the loop with
`atomicswap` at the frame-top drain and frozen into the published status
([§9.8][s9-8]).

<a id="g-drain"></a>**drain** — the single point at the top of each frame where the loop takes
each staging cell by `atomicswap` and applies it through the attach-compiled
scatter, in attachment order; never at a `t*` boundary, and under the roster
freeze it performs no checks at all ([§9.1][s9-1], [§9.4][s9-4]). The diagnostic
cells are taken at the same point ([§9.8][s9-8]).

<a id="g-framework-status"></a>**framework status** — the concrete frozen value each snapshot carries beside
the signal table: the pacer diagnostics ([§8.7][s8-7]) plus, per writer, this
boundary's drained diagnostics (`recent`), the counts the ring refused
(`suppressed`), the loop's cumulative per-writer × per-kind counters copied in
(`totals`) and the liveness timestamp ([§9.8][s9-8], [§9.2][s9-2]).

<a id="g-greedy-claim"></a>**greedy claim** — the claim a binding declaring `is_greedy` receives: the
unclaimed complement computed by the framework at attach instead of returned
by `claims`, ordinary in every respect afterwards; an empty remainder is legal
and reported (`EmptyGreedyClaim`), and the shipped GUI binding is the shipped
instance ([§9.3][s9-3], [§9.6][s9-6]).

<a id="g-harness-cell"></a>**harness cell** — the always-present staging cell of the harness register,
written by `stage!(sim, "face" => value, …)` from the calling task itself:
ordinary batches, traced and surface-checked, drained last by convention
([§10.6][s10-6], [§9.3][s9-3]).

<a id="g-harness-register"></a>**harness register** — the framework-owned write path of the calling task —
`stage!(sim, …)` and its cell — and the design's sole *derived* surface: the
unclaimed complement, the faces no rostered device claims, recomputed at every
stopped-sim roster change; a write to a claimed face is `ClaimedFaceEntry`
([§9.3][s9-3], [§10.6][s10-6]).

<a id="g-latest"></a>**`latest`** — the `@atomic` reference a published snapshot is release-stored
into and readers acquire-load; `latest(sim)` hands the calling task the same
immutable value a device handle gets ([§9.2][s9-2]).

<a id="g-next-snapshot-wait"></a>**next-snapshot wait** — `wait_next_snapshot(handle)`: the boundary counter
plus one `Threads.Condition` under the canonical predicate loop
(`counter > last_seen && running`), newest-wins, no queues, no per-frame reset
([§10.3][s10-3]).

<a id="g-operator-interrupt"></a>**operator interrupt** — Ctrl-C in an interactive session, read as a
control-plane stop rather than a failure: delivery is masked across the boundary
macro-sequence (`disable_sigint`) and raised at a frame-top or wait unmask
point, so the run takes the ordinary graceful tail and ends `stopped`. A second
one collapses the device joins; outside the REPL, SIGINT still kills the process
([§10.4][s10-4], [§10.1][s10-1]).

<a id="g-orphaned-claims"></a>**orphaned claims** — the claims of a device whose task died mid-run. Death is
not detach: the roster entry and claims persist to run end, the slots hold
their last-drained values, and the GUI renders the fact in the widget's
provenance; recovery is between runs ([§9.3][s9-3]).

<a id="g-peek"></a>**peek** — the GUI display rule: a widget shows its own pending write if any,
else the snapshot value. Own-cell only, which is what makes multi-click
counting and paused editing correct ([§9.7][s9-7]).

<a id="g-periphery"></a>**periphery** — everything outside the loop that exchanges data with it — GUI,
input devices, network I/O, logging — together with the concurrency model
binding them: staged writes inbound, snapshot reads outbound, control on its
own surface ([§9][s9], [§10][s10]).

<a id="g-roster"></a>**roster** — the list of attached device entries (binding, claims, stable
device id, attachment order): a plain immutable value the loop reads once at
`run!`, since `attach!`/`detach!` are stopped-sim operations ([§9.3][s9-3]).

<a id="g-scenario-component"></a>**scenario component** — the home of a sim-time script under the mid-run
mutation doctrine: an ordinary periodic discrete component executed
synchronously in the loop, deterministic paced or unpaced and replayed by
recomputation. The clock is the criterion — wall-clock interactions are
devices ([§10.5][s10-5]).

<a id="g-selector"></a>**selector (read-selector family)** — the closed set of deferred reads
`get_state`/`get_deriv`/`get_output`/`get_slot`/`get_face`, each
resolving against a source (table sources — a boundary snapshot or a service
evaluation's scratch tables — vs. live stores) before any client policy
applies ([§14.4][s14-4]).

<a id="g-should_abort"></a>**`should_abort`** — the per-attachment failure policy, an `attach!` keyword
defaulting to `false`: set, a device's departure — loop body returning, crash,
or a failed `init!` — also requests a control-plane stop; clear, the run
continues with the device absent and its claims held to run end. An attachment
fact, never a device property, the same device being advisory in one deployment
and load-bearing in another; the shipped GUI attaches with `true` ([§9.6][s9-6],
[§10.4][s10-4]).

<a id="g-snapshot"></a>**snapshot** — the immutable per-boundary publication: boundary-consistent
signal table (root slots included), `t`, boundary index and
framework status. It deliberately carries no state stores — the state
trajectory is derived data ([§9.2][s9-2]).

<a id="g-stage-on-interaction"></a>**stage-on-interaction** — the GUI staging contract: value widgets stage the
new level on edit, edge widgets on activation as a level computed from the
peek; held buttons do not re-stage, and no widget stages per render pass
([§9.7][s9-7]).

<a id="g-unattended-run"></a>**unattended run** — a run with empty staging and no snapshot readers: the
same loop, fully synchronous on the calling task, rethrowing after the
shutdown tail so CI fails honestly ([§9.1][s9-1], [§13.4][s13-4]).

<a id="g-write-surface"></a>**write surface** — the set of faces a writer's batch entries may reach: a
device's claim set, whether returned by `claims` or computed under
`is_greedy` ([§9.6][s9-6]), and for the harness register the derived
unclaimed complement. Static per run and enforced entirely at
staging — `OutOfClaimEntry` for a device, `ClaimedFaceEntry` for the harness
([§9.3][s9-3]).

### D.7 Recording and replay

<a id="g-decimation"></a>**decimation** — the log's keep-every-kth retention policy (`log_every`),
admissible on the log alone because it is derived data; every boundary still
runs, publishes to live readers and enters the trace. Bounded by `log_max`, the
maximum number of retained snapshot references (default finite, `Inf` the
opt-out): when the log fills, the effective stride doubles — *progressive
re-decimation*, so coverage stays global at `log_every · 2^k` instead of
collapsing to a rolling window — with the boundary-zero and terminal snapshots
retained unconditionally and outside the bound. A view policy throughout, never
trajectory-determining ([§9.2][s9-2]).

<a id="g-frame-ordinal"></a>**frame ordinal** — the trace's key: replay applies the recording's batches
for frame *k* at frame *k*, exact because the frame sequence is itself
deterministic under replay ([§10.7][s10-7], [§9.1][s9-1]).

<a id="g-log"></a>**log** — the retained sequence of published snapshots (the same objects, no
copies), with a plain kill switch and `log_every` decimation; derived data,
recomputable from the trace by replay ([§9.2][s9-2]).

<a id="g-recorders"></a>**recorders** — the trace and the log jointly, cleared together at `init!` and
at a trim commit so they restart with the run they record ([§10.6][s10-6], [§14.8][s14-8]).

<a id="g-replay"></a>**replay** — the ordinary loop with exactly two substitutions: boundary zero
from the trace header, and a drain reading the trace by frame ordinal. It
re-records, ends `initialized`, and validates the header — stores, slot
faces and deployment block — up front, applying the header's `t₀`
([§10.7][s10-7]).

<a id="g-run-metadata"></a>**run metadata** — the trace header's deployment block: `t₀`, `Δt_base`,
`h`, `n`, the algorithm identifier, `localization_tol`, `localization_budget`,
`firing_budget` and the
effective `t_end`/`stop_on`
pair ([§9.5][s9-5], [§13.5][s13-5]).

<a id="g-trace"></a>**trace** — the primary record of a session: the sequence of drained,
device-tagged batches per frame, plus its header. On by default, because the
log is recomputable from the trace and never the reverse ([§9.5][s9-5]).

<a id="g-trace-header"></a>**trace header** — the trace's preamble: the resolved initial stores
`(x, m)`, the initial root-slot values, each writer's face-name →
position schema, and the deployment block — captured after `apply!` and the
slot writes, before the boundary-zero sequence runs ([§9.5][s9-5], [§14.5][s14-5]).

<a id="g-trace-record"></a>**trace record** — the retained form of a drained batch, uniform for every
writer: (position ⇒ value) pairs for the non-`nothing` entries, converted at
the drain against the header's schema, so trace size tracks information
rather than surface width and consumers meet one format and one replay path
([§9.5][s9-5], row 176).

<a id="g-what-if-register"></a>**what-if register** — replaying a trace against the same structure with
changed parameters: deterministic re-driving of the recorded inputs through a
modified model, promising determinism but never bit-identical reproduction
([§10.7][s10-7]).

### D.8 Stopped-sim services and the condition algebra

<a id="g-at"></a>**`at` / `Scoped`** — the scoping combinator: `at(prefix, node)` stores a
prefix beside a condition node and applies nothing, path concatenation
happening once at resolution. It also lifts whole `TrimProblem`s and
linearization tap sets ([§14.2][s14-2], [§14.9][s14-9]).

<a id="g-baseline"></a>**baseline** — an aircraft-shipped, full-coverage condition function
(`ready_for_taxi(ac)`, `cold_and_dark(ac)`) layered under tweaks by
`override`, and the `baseline` keyword `init!`/`trim!` take ([§14.6][s14-6]). Not to be
confused with an event *prior* ([§D.4][sD-4]).

<a id="g-boundary-zero"></a>**boundary zero** — the initialization boundary: the ordinary macro-sequence
with an empty integrate — project → [sweep → guards → handlers]\* → due `g`
updates → header and first snapshot — run at `t₀` once `apply!` has
established the stores ([§14.5][s14-5]).

<a id="g-capture"></a>**capture** — the service reading the current committed stores *and* root
slots back as a condition value, returning `(condition, t)`; the gather twin
of `apply!`, and what makes warm restart need no second semantics ([§14.1][s14-1],
[§14.10][s14-10]).

<a id="g-component-test-rig"></a>**component test rig** — a one-child assembly exporting the child's entire
input face set, so any component can be built and simulated in isolation; an
abstract entry is satisfied *inside* the rig by a concrete **stub child**
wired to that face ([§13.7][s13-7]).

<a id="g-condition"></a>**condition** — the datum that says "set this build to this state": a
path-addressed sparse overlay on the declared defaults, covering `x` and `m`
fields plus root slots by face — never outputs, never workspace ([§14.1][s14-1]).
[§14][s14] owns the word; a guard defines a *predicate* ([§D.4][sD-4]).

<a id="g-design_world"></a>**`design_world`** — the shipped thin world (aircraft +
`SimpleAtmosphere(wind = NoWind())` + `HorizontalTerrain`) that mounts an
aircraft for trim and linearization; "aircraft as root" is the shallowest
world, not a special case ([§14.9][s14-9]).

<a id="g-fragment"></a>**fragment** — the leaf node of the condition algebra: `fragment(; x, m,
slots)` payloads speaking only about the component at the authoring point
(**self-vocabulary**), with addressing left entirely to `at` ([§14.2][s14-2]).

<a id="g-fragment-tree"></a>**fragment tree** — the inert, lazy composition of `Fragment`/`Scoped`/
`Merged`/override nodes; isbits but for the interned prefixes, so rebuilding
it per trim iteration is stack-only construction ([§14.2][s14-2]).

<a id="g-merge"></a>**merge** — the symmetric, collision-intolerant combinator over condition
nodes: a duplicate leaf is an error naming both provenance chains, and mixing
a node with a bare `NamedTuple` is an error method, not `Base.merge` ([§14.2][s14-2]).

<a id="g-mounting"></a>**mounting** — relocating a whole problem or tap set with `at(prefix, …)`:
every field is either condition-producing (path-relative, post-composed) or
path-free, so the service never knows where its paths sit ([§14.9][s14-9]).

<a id="g-override"></a>**override** — the ordered, asymmetric layering combinator: on a shared leaf
the patch wins and provenance keeps both sources, while collisions *within* a
layer remain errors; variadic ([§14.6][s14-6]).

<a id="g-service-lifecycle"></a>**service lifecycle** — the `Simulation` states `built` / `initialized` /
`running` / `stopped` / `errored` ([§10.6][s10-6]) and each service's legality against
them; a violation is `ServiceLifecycle`, and `errored` is terminal for all
four services ([§14][s14]).

<a id="g-slot-totality"></a>**slot totality** — the pre-write requirement that an application establishing
a complete world over virgin stores — `init!`, trim setup, trim commit —
cover every root slot. Conditions themselves are legitimately partial; a
shortfall is `UninitializedSlots`, collected and declaration-ordered, leaving
the simulation untouched ([§14.6][s14-6]).

<a id="g-taps"></a>**taps** — the three selector lists (`x`, `u`, `y`) declaring what
linearization seeds and reports, with an optional component index so a vector
leaf yields named scalars; validated at resolution (`TapResolution`) and
relocatable via `at` ([§14.10][s14-10]).

<a id="g-trimproblem"></a>**`TrimProblem`** — the closed seven-field value
`guess`/`lower`/`upper`/`condition`/`reads`/`residuals`/`tolerances`: an
*implicitly specified* condition, solved as a square root-find over named
residuals and committed as an `init!` of `override(baseline, solution)`
([§14.7][s14-7], [§14.8][s14-8]).

### D.9 Error discipline and diagnostics

<a id="g-carrier-exception"></a>**carrier exception** — the single exception a set of diagnostics travels in:
`BuildError` thrown at a stratum barrier, `StepError` at the runtime catch
site. Diagnostics themselves are plain values ([§13.2][s13-2], [§13.4][s13-4]).

<a id="g-collect-the-checks-fail-the-evaluations-fast"></a>**collect the checks, fail the evaluations fast** — the reporting policy:
declarative passes over collected structure return their full violation list,
while the first user-code exception aborts the phase; strata are barriers, and
the site column spells the collecting case "build (collected)" ([§13.1][s13-1],
[Appendix C][sC]).

<a id="g-did-you-mean"></a>**did-you-mean** — the required shape of any name-shaped failure: the
offending name plus the list-in-hand it should have come from, carried as
payload rather than baked into message text ([§13.2][s13-2]).

<a id="g-error-locality"></a>**error locality** — the property the declaration layer buys: a mistake fails
at the site of the mistake, not later and inside correct code. The five
walkthroughs ([§11.4][s11-4]) are its grounding cases and the acceptance tests
([§11.4][s11-4]).

<a id="g-execution-cursor"></a>**execution cursor** — the plain mutable field recording where in the compiled
schedule execution is (component path, function, boundary phase); one cheap
store per dispatch, so a runtime failure gets its frame without exception
frames in the hot path ([§13.4][s13-4]).

<a id="g-feedthrough-tracer"></a>**feedthrough tracer** — the set-propagation instrument (global value-blind,
or local primal-carrying at sampled states) used to classify a rejected cycle
as real or artificial; diagnostic only, never an input to scheduling ([§5.6][s5-6]).

<a id="g-kind"></a>**kind** — a diagnostic's identity in the closed set enumerated normatively in
[Appendix C][sC], with payload fields, owning section and severity; tests match on
kind plus payload, never on message text ([§13.2][s13-2]). Not a component *class*
([§D.1][sD-1]) or a *function family* ([§D.1][sD-1]).

<a id="g-payload"></a>**payload** — the structured data a diagnostic carries beside its kind: paths
and names as strings (never instances or model types), expected/observed port
types, the list-in-hand, the severity ([§13.2][s13-2], [Appendix C][sC]).

<a id="g-stop_on"></a>**`stop_on` / termination is a state** — graceful termination is model state,
never an exception: detection is ordinary event machinery, publication an
ordinary root-exported `Bool` output face, and `stop_on` the deployment policy
naming the faces the loop reads after every published boundary ([§13.5][s13-5]).

<a id="g-warning-streams"></a>**warning streams** — two, scoped separately: the *build* stream, whose
warning set is deliberately empty, and the *runtime* stream — per-occurrence,
carried by the per-writer diagnostic cells that structurally rate-limit it
([§9.8][s9-8]), surfaced through published
framework status, with its committed inventory listed in [§13.2][s13-2].

### D.10 Meta-vocabulary

<a id="g-blessed"></a>**blessed** — the spec's marker for a practice it explicitly sanctions where a
neighboring one is forbidden: derivation from other declarations ([§11.2][s11-2]), the
one spot where evaluation feeds structure ([§12.1][s12-1]), the workspace-plus-snapshot
idiom for zero-allocation ticks ([§7.3][s7-3]).

<a id="g-the-freeze"></a>**the freeze** — the roster freeze: `attach!`/`detach!` are stopped-sim
operations, so the roster, its claims and the run's partition of the root face
set into write surfaces are static, inspectable facts of each run ([§9.3][s9-3], rows
106–107).

<a id="g-guarded-addition"></a>**guarded addition** — a capability the design admits but does not build,
weighed against Flight.jl's fundamental strengths and recorded with its shape
so adoption stays additive ([§1][s1]; e.g. field-addressed staging,
[§4.3][s4-3]; mid-run reader attach, [§9.3][s9-3]).

<a id="g-normative"></a>**normative / index, not a second home** — the spec is the normative statement
of the design, and its appendices are indices: each recall line's normative
statement stays in the owning section ([Appendix A][sA], and the same rule for
Appendices B, C and D). The design directory's walkthrough explainers are
non-normative companions by their own preambles.

<a id="g-recorded-not-built"></a>**recorded, not built** — the disposition of a worked-out extension
deliberately left unimplemented, with its seams named so adoption is additive
(the closed-loop trim, [§14.7][s14-7]; the sampled-data `Dual` activation and
declarative non-participation, [§14.10][s14-10]).

<a id="g-register"></a>**register** — the spec's word for a mode or idiom in which something is done,
always compounded: the didactic register ([§13.2][s13-2]), the inspection and
integration registers ([§9.2][s9-2]), the by-allocation register ([§11.2][s11-2]), the
harness, unattended and what-if registers ([§10.6][s10-6], [§10.7][s10-7]). Reserved for this
sense — the recording artifacts are the *recorders* ([§D.7][sD-7]).

<a id="g-row"></a>**row** — a numbered entry of `framework_decisions.md`, cited throughout as
"row N": one settled decision with the alternatives weighed against it. Row
numbers are stable and never reused, and each row states its *current*
position, a superseded one being demoted into its rejected-alternatives column
rather than rewritten ([§1][s1]).

<a id="g-seam"></a>**seam** — a narrow, named interface kept deliberately thin so what sits
behind it can be replaced or measured: the stepper seam ([§8.2][s8-2]), the backend
seam ([§14.8][s14-8]), the measurement seam ([§12.7][s12-7]), the phase-body seams of the
compiled executor ([§12.7][s12-7]).

<a id="g-torture-test"></a>**torture test** — an existing, maximally awkward artifact transliterated
against a proposed mechanism to validate it before adoption: `PistonEngine`
and the FCS cascade against [§5.2][s5-2] ([§15.2][s15-2]), filter/joystick/GUI against the
[§9][s9] staging shapes ([§15.3][s15-3]), the strapdown IMU against the leaf split ([§15.5][s15-5]); the
standard component library is the standing ergonomics one ([§13.7][s13-7]).

<a id="g-worked"></a>**worked (example)** — a full spelling of a mechanism against a real artifact,
carried in the spec rather than left to the reader: the worked assembly of
[§11.6][s11-6], the worked C172 cruise problem of [§14.7][s14-7], and the IMU
([§15.5][s15-5]) as the boundary-sampling example [Appendix A][sA] points at.

<!-- citation link definitions — generated by tools/linkify.jl; do not edit -->
[s1]: #1-purpose-and-method
[s10]: #10-runtime-periphery-lifecycle-and-orchestration
[s10-1]: #101-control-plane
[s10-2]: #102-loop-scheduling-wait-primitive-yields-thread-budget
[s10-3]: #103-the-next-snapshot-wait
[s10-4]: #104-shutdown-protocol
[s10-5]: #105-scripts-and-the-mid-run-mutation-doctrine
[s10-6]: #106-run-lifecycle-and-partial-advance
[s10-7]: #107-replay-the-trace-re-drives-the-ordinary-loop
[s11]: #11-the-declaration-layer-components-and-assemblies
[s11-1]: #111-position-a-declarative-trait-layer--plain-julia-no-macros
[s11-2]: #112-the-declaration-inventory
[s11-3]: #113-visibility-the-contract-is-the-interface
[s11-4]: #114-failure-walkthroughs-the-error-locality-grounding
[s11-5]: #115-assembly-declaration-type-based-class-by-declaration-shape
[s11-6]: #116-paths-wiring-and-faces
[s11-7]: #117-rate-scopes
[s11-8]: #118-computed-connections-and-generic-boundaries
[s12]: #12-the-build-pipeline
[s12-1]: #121-three-strata
[s12-2]: #122-the-build-artifact
[s12-3]: #123-probing-and-input-synthesis
[s12-4]: #124-activations-executable-sets-laziness-caching
[s12-5]: #125-the-always-on-conformance-check
[s12-6]: #126-stopped-sim-services-as-stratum-c-clients
[s12-7]: #127-the-compiled-executor
[s13]: #13-error-discipline
[s13-1]: #131-reporting-policy-collect-the-checks-fail-the-evaluations-fast
[s13-2]: #132-diagnostics-structured-values-one-carrier-exception
[s13-3]: #133-build-primitives-resolve-and-the-face-list-accessors
[s13-4]: #134-runtime-failures-one-catch-site-an-execution-cursor
[s13-5]: #135-termination-is-a-state-not-an-exception
[s13-6]: #136-abnormal-shutdown-one-tail-two-entries
[s13-7]: #137-tooling-consequences-provenance-and-the-component-library
[s14]: #14-stopped-sim-services
[s14-1]: #141-conditions-are-path-addressed-overlays-on-the-declared-defaults
[s14-10]: #1410-linearization-tap-selectors-one-seeded-pass-a-pure-query
[s14-2]: #142-fragment-composition-locality-without-schema
[s14-3]: #143-resolution-flatten-validate-compile-once
[s14-4]: #144-two-application-registers-over-one-plan
[s14-5]: #145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions
[s14-6]: #146-slot-totality-the-missing-value-error-and-the-override-combinator
[s14-7]: #147-the-trim-problem-namedtuple-decisions-declared-reads-named-residuals
[s14-8]: #148-the-trim-service-solver-seam-scratch-stores-commit-and-report
[s14-9]: #149-mounting-problems-as-relocatable-values
[s15-1]: #151-vehicle-today--this-framework
[s15-2]: #152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade
[s15-3]: #153-torture-test-for-the-9-staging-shapes-filter-joystick-and-gui
[s15-4]: #154-the-interactive-c172x-demo-the-periphery-under-load
[s15-5]: #155-the-strapdown-imu-integrate-and-dump-across-the-tier-boundary
[s16]: #16-open-axes
[s2]: #2-formalism
[s2-1]: #21-events-two-detection-policies
[s2-2]: #22-exclusions-deliberate
[s3]: #3-component-taxonomy
[s3-1]: #31-continuous-component-the-hybrid-primitive
[s3-2]: #32-periodic-discrete-component
[s3-3]: #33-assembly
[s4-1]: #41-immutable-value-semantics
[s4-2]: #42-consumers-see-ports-not-stages
[s4-3]: #43-table-mechanics-and-port-granularity
[s4-4]: #44-function-valued-signals-environment-access
[s5]: #5-evaluation-order-and-feedthrough
[s5-1]: #51-the-scheduling-problem
[s5-2]: #52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws
[s5-3]: #53-structural-feedthrough-stage-roles-schedule-and-step-boundaries
[s5-4]: #54-artificial-loops-and-the-escape-hatch
[s5-5]: #55-algebraic-loop-policy-reject-at-build-time
[s5-6]: #56-diagnostics-feedthrough-tracing
[s6-1]: #61-connections-and-hierarchy
[s6-2]: #62-aggregation-explicit-summing-junctions
[s7]: #7-state-and-data-representation
[s7-1]: #71-continuous-state-structured-immutable-flat-backing
[s7-2]: #72-numeric-genericity-eltype
[s7-3]: #73-discrete-state-modes-and-workspace
[s7-4]: #74-the-fused-evaluation-lineage-prior-art-and-how-we-got-here
[s7-5]: #75-allocation-policy-a-scoped-invariant
[s8]: #8-time-and-execution
[s8-1]: #81-loop-ownership-the-framework-owns-the-simulation-loop
[s8-2]: #82-the-stepper-seam
[s8-3]: #83-signal-table-consistency-is-a-boundary-property
[s8-4]: #84-localization-mechanics
[s8-5]: #85-multi-rate-tick-scheduling
[s8-6]: #86-event-iteration-at-boundaries-to-quiescence-budgeted
[s8-7]: #87-real-time-pacing
[s9]: #9-runtime-periphery-the-data-plane
[s9-1]: #91-no-shared-mutable-model-staged-writes-snapshot-reads
[s9-2]: #92-outbound-snapshot-publication
[s9-3]: #93-inbound-root-input-slots-claims-and-the-frozen-roster
[s9-4]: #94-inbound-per-device-staging-representation-and-the-drain
[s9-5]: #95-inbound-the-input-trace
[s9-6]: #96-devices-one-authoring-contract-no-taxonomy
[s9-7]: #97-the-gui-write-path-port-resolution-peek-staging-contract
[s9-8]: #98-diagnostics-and-liveness-the-per-writer-cell
[sA]: #appendix-a-taught-contracts-the-author-facing-index
[sB]: #appendix-b-api-synopsis-the-entry-points
[sC]: #appendix-c-the-diagnostic-kind-set
[sD-1]: #d1-component-model-and-declaration-layer
[sD-4]: #d4-time-and-events
[sD-5]: #d5-build-pipeline
[sD-6]: #d6-runtime-periphery
[sD-7]: #d7-recording-and-replay
[sD-8]: #d8-stopped-sim-services-and-the-condition-algebra
[sD-9]: #d9-error-discipline-and-diagnostics
