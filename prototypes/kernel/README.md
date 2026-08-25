# Kernel prototype

The walking skeleton for the framework in `docs/notes/design/framework_spec.md`,
built to keepable standards and grown one increment at a time. Increment 1 (the
cell-store representation bench) lives in `../cellstore_bench` and stays frozen
there — D-162 cites its numbers.

**Increment 2 — the continuous tier walks.** A model builds, schedules itself
from its own feedthrough structure, integrates, and does it without allocating.

**Increment 3 — the discrete tier, at one rate.** Both tiers now share the
model: tier is read off the declaration shape, discrete state and modes live in
their own stores, `g` runs at boundaries after the output stages, and the ZOH
holds mid-step because the interior sweep has no discrete entries to gate.
Every discrete component sits at `D = 1, Φ = 0`; the grid arrives in
increment 5, after assemblies.

**Increment 4 — hierarchy and assemblies.** `build` takes the root component
instance. Class is read off declaration shape, children are struct fields and
container elements, and the three connection declarations resolve against
slash-separated paths into flat primitives with one producer per input.
Assemblies are virtual for execution: what reaches the executor is exactly what
increment 3 already ran. The whole-tree obligation model replaces the silent
unwired face — an input fed by nothing is a build error, and the root
assembly's own input faces are the slots.

**Increment 5 — the multi-rate grid.** The two-register `sample_times`
declaration folds to one `(anchor, m, c)` triple per component in Stratum A;
deployment binding — now a real seam, `Simulation(build(root); h, n, Δt_base)`
— resolves the triples to `(D, Φ, Δt)` by exact GCD arithmetic and compiles the
executor behind it, one deployment-free `Build` backing many `Simulation`s. The
boundary sweep gates discrete entries by `(idx − Φ) % D`, and a step boundary
that is not a base tick runs the zero-arg bodies: an empty due set is arity
selection, never a sentinel index failing every gate (D-185).

**Increment 6 — events at boundaries.** The declaration surface — ordered named
`events` of `Event(guard, handler)` pairs, positional `project` — and the whole
§10.6 machinery behind it: detection policy read off the guard's probed return
type, the handler return law held key by key at the probe, and the boundary
macro-sequence in its final form, integrate → project → [sweep → guards →
handlers] iterated to quiescence under `firing_budget` → all due `g` updates.
The three per-event registers are plain vectors on the `Simulation`, and
boundary zero derives fires-at-`t₀` from all-not-holding priors.

**Increment 7 — localization.** The runtime now consults the policy the probe
recorded: a sign-form guard's crossing is bracketed by trial evaluations over
the cubic Hermite interpolant and fired at a `t*` boundary strictly inside the
frame. The frame loop runs §10.4's chain — arrival sweep and trigger, the
θ = 0 validation with its epoch discriminator, the lazily-paid ẋₙ₊₁, ITP
root-finding, the `t*` boundary as the off-tick §10.6 sequence, the remainder
step targeting the indexed frame top — iterated under `localization_budget`,
with `localization_tol` and `localization_budget` joining `firing_budget` as
validated deployment keywords. The stand-ins table empties.

**Increment 8 — the stepper seam's second backend.** §10.2's seam becomes
dispatch: `step!`, `startpoint` and `dense!` on a stepper type, RK4's scratch
moving off the `Simulation` into its own struct, and Heun joining as
`method = Heun` — a deployment keyword defaulting to `RK4`. The frame loop's
reach into RK4's work registers is gone: localization reads the seam's
retained pair and dense output and never a scratch index, and the fitted
convergence orders — 4 and 2 through one unchanged deployment surface and one
unchanged localization machinery — are the interchangeability proof.

**Increment 9 — the data plane's core exchange.** §11.1's planes 1–3 with the
roster empty: the harness register's staging cell — the positional
`Union{Nothing, Tᵢ}` batch compiled over the root face set, `stage!`'s shim
with every check at staging, CAS merge under newest-wins-per-face — the
frame-top drain by one `atomicswap` and a compiled scatter, and snapshot
publication, one buffer copy per frame-top boundary behind a release/acquire
`@atomic latest` read by `latest(sim)`. `run!`'s frame takes its §11.1
anatomy — drain, integrate, boundary sequence, publication — and the epoch
tests' slot writes move onto the machinery the θ = 0 discriminator was built
expecting. Two stand-ins enter the table.

**Increment 10 — the roster and claims.** §11.3 wholesale over increment 9's
cell machinery: `attach!`/`detach!` as stopped-sim operations behind the run
freeze, the three-part admission in spec order, both claim sources — the
enumeration called once and validated, the unclaimed complement computed at
the attach instant and never recomputed — per-device staging cells compiled
over claim sets exactly as the harness register's is over the complement that
remains, that register recompiled at every roster change with its pending
batch renormalized through the new schema, and the drain over the roster in
attachment order, harness last. From §11.6, exactly the attach-point slice:
the two periphery roots, the declared sides, the enumeration contract with
its error-throwing fallback, and the bidirectional conformance check. One
stand-in enters the table: device staging without handles.

    julia --project=. test/runtests.jl

## What is real here

| piece | spec | file |
| --- | --- | --- |
| leaf walk: flatten / reconstruct / the activation retype | §7.1, §7.2 | `src/leaves.jl` |
| declaration layer, both tiers' disjoint name families and arities, the bundle law, `probe_value`, `AbstractComponent`, the three connection declarations, the two rate registers (`Relative`/`Absolute` over `Period`/`Hz`, plain data carriers) with `sample_times`, and the event surface — `Event(guard, handler)` with no detection keyword, the ordered named `events` collection, positional `project`, the guard/handler name set and the §2.1 predicate | §5.2, §8.2, §8.6, §8.7, §9.3, §2.1, D-179, D-185, D-195 | `src/declare.jl` |
| class by declaration shape (with `events` and `project` in the leaf family), children and containers, paths and the reach rule, endpoint and face resolution, the flatten pass with the obligation model and the sample-time fold to `(anchor, m, c)` triples | §8.5, §8.6, §6.1, §9.1 | `src/assembly.jl` |
| per-eltype cell stores and the store bundle | §9.7, D-162 | `src/store.jl` |
| entries, chunked unrolled walk, the interior/boundary split, the `(idx − Φ) % D` gate on discrete boundary entries, and the event machinery — event/projection entries, the register vectors (§10.6's three plus §10.4's policy mask and σ samples, the guard walk feeding both) and the guard/fire/project walks with handler-store latching | §9.7, §10.5, §10.6, §10.4, §5.3 | `src/executor.jl` |
| tier classification, probe, feedthrough graph, layout, embed-accept, the `Build` artifact, the activation seam (nominal at build, `activation(b, T)` as a cached Stratum-C re-run, frozen products carried across, the eager `activations` keyword), deployment binding (three cross-validated `Δt_base` sources, the GCD constraint pool, per-anchor division pairs), per-deployment entry compilation, and the event probe — both-halves at declaration reading, policy off the guard's return type onto `Build.policies`, the handler return law key by key, `project` held complete | §8.2, §9.3, §5.3, §9.1, §9.2, §9.4, §10.4, D-166, D-179 | `src/build.jl` |
| `Simulation` with its deployment constructor, bound schedule and the `method`/`firing_budget`/`localization_tol`/`localization_budget` keywords, the seam's framework side (the never-entered-empty short-circuit), the boundary macro-sequence in its final form (project → iterated event phase → due updates) at ticks and off ticks, the §10.6 iteration with its three registers and the `FiringBudget` degradation, `phase_bodies`, the path- and face-addressed table accessors | §10.2, §10.3, §10.5, §10.6, §10.4, §9.7, §11.3 | `src/sim.jl` |
| the stepper seam's backend side: the three-clause contract as dispatch — `step!` by arbitrary `h`, dense output via the retained `startpoint` pair and `dense!` (the shared cubic Hermite both backends inherit), one-step methods only — with fixed-step RK4 (the default) and Heun, each owning exactly its own scratch | §10.2, D-017 | `src/stepper.jl` |
| the frame loop: arrival sweep and trigger, the θ = 0 validation with its epoch discriminator, the lazily-paid arrival derivative, bracketed trial evaluations over the seam's dense output under ITP, `t*` boundaries and the remainder step against the indexed grid, `localization_budget` counting `t*` boundaries per frame, and the `ChatteringBudget` degradation | §10.4, §10.2, D-018, D-133 | `src/localize.jl` |
| the data plane's core exchange: the compiled writer — one schema/shim/merge/scatter unit over any write surface, with every check at staging and the out-of-schema warning discriminated by writer — staging cells, the frame-top drain via `atomicswap` behind stopped-sim-compiled thunks, `run!`'s frame anatomy (drain → integrate → boundary sequence → publication), and publication — the boundary-consistent whole-table snapshot with the state stores excluded, `capture`'s buffer copy, the release/acquire `@atomic latest` pair, and `latest(sim)` as the inspection register | §11.1, §11.3, §11.4, §11.2, §12.6 | `src/dataplane.jl`, `src/sim.jl` |
| the roster and claims: `AbstractDevice`/`AbstractBinding` with the declared sides and the false root defaults, the enumeration contract and the bidirectional conformance check (`which` against the fallback), `attach!`/`detach!` behind the §11.3 freeze with the three-part admission in spec order, both claim sources with `EmptyGreedyClaim` for the staked empty remainder, monotonic never-reused device ids, per-device cells over claim sets, the harness register as the derived remainder — recompiled at every roster change, its pending batch renormalized (`ClaimedFaceEntry` at the seam), emptied outright by a rostered greedy (D-192) — and the attachment-order drain, harness last | §11.3, §11.4, §11.6 | `src/roster.jl`, `src/sim.jl` |
| the coverage set: `Plant`, `Gain`, `Sum`, `DiscreteIntegrator`, `TickCounter`, `Smoother`, `WorkGain`, `ModedSource`, `ZOH`, `Ramp`, the event set (`Trigger`, `Follower`, `Sawtooth`, `Rotor`, `Chatterer`, `TwoShot`, `Preempted`), the localization set (`Stamper`, `GatedStamper`, `Bouncer`, `Relaxer`), plus `Group`, the named `SampledLoop`/`Vehicle`/`MultiRate` assemblies, and the periphery set — `Pad`/`Panel` devices (mutable, `===` being identity) and the `Enumerated`/`Greedy` bindings, one per claim source | §5.2, §6.2, §7.3, §8.5, §10.4, §10.5, §10.6, §11.3, §11.6 | `src/library.jl` |

The properties the tests pin down, each of which is a spec claim rather than a
programming convenience:

- **The schedule is derived, not authored.** Stage-1 (`h_x` continuous, `h_s`
  discrete) ports carry no
  input dependence, so consuming one adds no edge to the feedthrough graph. The
  reference model's feedback path closes through the plant's `h_x` port and
  schedules as `sum → ctl → plant`; rewiring the same loop through the plant's
  stage-2 port makes it a genuine algebraic loop and the build rejects it by
  name.
- **Derivative completeness is structural.** `f` returns a value shaped like
  `init_x`; the probe checks the shape once and the runtime scatter into the
  flat `ẋ` block is then safe by construction. A forgotten field is impossible,
  not merely detectable.
- **The phase-body roster is fixed and total.** `phase_bodies(sim)` always
  returns all four bodies. A model with no discrete components still gets
  `ticks`, empty — it compiles to a no-op, and its `@ballocated` assertion
  passes vacuously, which is the point: consumers iterate the roster with no
  per-model branching.
- **Tier is read off the declaration shape, never announced.** For a stateful
  leaf the whole name family carries it — `init_x`/`h_x`/`h_xu`/`f` continuous
  against `init_s`/`h_s`/`h_su`/`g` discrete, the two disjoint (D-195), with
  the update law the decider; for a stateless one the contract arity does.
  Every other tier-implying declaration must agree, and disagreement names the
  offending one — `g` beside a two-argument `output_types`, `init_m` on a
  discrete leaf, both arities of one contract, and the wrong-letter case the
  split families restore: a continuous stage name on a leaf whose update law is
  `g`.
- **The ZOH is not implemented.** It is the absence of any way to change a
  discrete cell mid-step: the interior sweep is compiled from continuous
  entries alone, so the hot path carries no gating test, and a discrete cell
  cannot move across a step because nothing in that walk writes it.
- **A path reaches exactly as far as the declaration knows.** The reach rule is
  about the declaring type's knowledge, not the instance's: the same
  `"inner/sum/a"` against the same `SampledLoop` value resolves when the field
  is declared `::SampledLoop` and is a build error when it is declared `::L`,
  because deep-routing *past* a generically held child hard-codes one
  implementation. Resolving *to* one is face-level access and legal (§13.3): a
  route through concretely-declared structure may end at a generic child's
  face, and a single hop always may — which is what lets a container's
  elements be addressed by their parent's own declarations.
- **Being fed is a whole-tree obligation, not a per-declaration one.** An unfed
  child input inside one assembly is merely awaiting a claim from above — a
  sibling wire, an ancestor's deep route, or an `input_connections` entry
  handing it up a level. The error fires at the root, for the chain that never
  terminates; the one legitimate terminus is a root input face, which is a
  slot. The one-producer rule spans levels the same way, so an ancestor's route
  onto an input a sub-assembly already wires is caught where the two claims
  meet, with both entries named.
- **A face is derived, never declared.** It owns no cell: it resolves — through
  as many levels of re-export as there are — to its ultimate internal endpoint,
  and takes that endpoint's type and tier. `Vehicle`'s `y` and `cmd` are one
  assembly's faces sourced from the two tiers, and at a `Dual` activation `y`
  walks while `cmd` stays pinned, with nothing anywhere declaring either.
- **An anchored divisor does not exist until `Δt_base` binds.** The `Build`
  carries `(anchor, m, c)` triples against exact rational anchors, and the same
  `Build` lands `MultiRate`'s gnss at `D = 10` or `D = 5` depending on the
  deployment keywords — with nothing writable shared between the two
  `Simulation`s, because each materializes its own stores and buffers. The
  worked example compiles to exactly §9.2's table.
- **Due-ness is one subtraction and one remainder, at boundaries only.** The
  interior walk is untouched by increment 5 — no index reaches it — and a step
  boundary that is not a base tick runs the zero-arg bodies, which *are* the
  boundary walk with every discrete entry gated out. The multi-rate sampled
  loop only matches its exact discretization at `Δt_ctl = D·Δt_base = 4h` if
  the gate admits the controller at exactly its own ticks, the hold spans the
  sub-ticks, and the bundle's `Δt` is the compiled schedule's — one test, three
  ways to fail it.
- **A cell holds what the build probe populated until a sweep first writes it.**
  Entry compilation seeds every cell from the probe products, so an offset
  component's pre-first-tick reads and a frozen component's pinned cells are
  the same story (§10.5, §9.3) rather than a lucky zero.
- **The nominal activation runs at build; any other is a Stratum-C re-run.**
  `activation(b, T)` materializes at first request and caches on the `Build`
  (§9.4); structure and schedule are computed once and shared. A frozen
  component sits outside the non-nominal executable set: its stages are never
  probed there, and its complete products — cells included — are carried
  across from the nominal activation, holding what a tick at `t₀⁻` computed
  from real upstream values rather than anything synthesized. That is what
  makes a discrete stage reading `t` lawful (it is `Float64` whenever the
  stage actually runs), and why a pinned-leaf lurk detonates at the `Dual`
  activation — or eagerly, via `build(root; activations = (Float64, D8))`,
  the CI idiom — rather than at `build`.
- **The guard's form is the declared policy.** `Event(guard, handler)` carries
  no detection keyword: the probe runs every guard at the nominal activation
  and its return type decides — `Bool` boundary-detected, the nominal scalar
  localized, anything else an error naming both admissible forms (D-179). The
  policy lands on `Build.policies` and compiles into the event set's mask,
  which is what the frame loop's trigger check reads — the mixed predicate
  rides the gate idiom `(gate) ? σ : -one(σ)`, whose return type keeps it
  localized. Guards and handlers are outside every non-nominal activation's
  executable set (D-052), so at a `Dual` deployment the event phase compiles
  to the bare sweep, a mode never moves, and the frame loop takes the bare
  step with no arrival machinery at all.
- **A cascade is logically simultaneous, not one link per step.** The event
  phase iterates [sweep → guards → handlers] to quiescence within one boundary,
  so a supervisor → follower → follower chain settles in three rounds at any
  `h` — the latency an integrator parameter must never buy. Eligibility inside
  the boundary is an edge like any other, read against the last-observed
  sample: a sticky predicate fires once, and firing is at most one event per
  component per round, first eligible in declaration order — priority with
  re-decision. The one register subtlety is D-191's: an eligible-but-blocked
  event's sample is *not* overwritten, so its edge stands into the next round,
  while a premise the earlier transition falsified re-decides against the
  post-transition sweep and simply does not fire.
- **Budget exhaustion degrades; it does not throw.** A toggling pair spends
  `firing_budget` per event at one boundary under a `FiringBudget` warning (at
  most once per event per boundary), the rest of the model iterates untouched,
  and the run proceeds — the quiescent samples become honest priors, so the
  chatterer presents no fresh edge afterward. The registers are detection
  bookkeeping in plain vectors, in no state store; boundary zero establishes
  every prior as not-holding, so a predicate holding in the authored state
  fires at `t₀`, and a warm restart fires it again — derived, not asserted.
- **Projection runs between a state write and its decode.** At every boundary
  `project` runs after the integrate and before the sweep — boundary zero
  included, so an off-manifold `init_x` lands on the manifold before the first
  cell is published — and each firing is `handler → project`. The probe holds
  its return complete against the state shape, which is what makes the
  wholesale buffer write-back safe by construction.
- **Due updates run after quiescence.** Ticks stay outside the iteration: at a
  wrap boundary the handler resets the state, the re-sweep publishes the
  post-transition value, and only then does a sampled consumer's `g` read it —
  one test, and the wrong order costs the full pre-reset value, orders of
  magnitude outside its tolerance.
- **`t*` is an observation, and the stamp proves it.** A handler reads `t` from
  its bundle, so a stamping handler is the direct observable: a localized
  firing lands within `localization_tol` of the true crossing — exactly, on a
  trajectory the cubic Hermite reproduces; at `O(h⁴)` on the rotor — and at
  the *holding endpoint* of the final bracket, never before the crossing. The
  trajectory-level counterpart is the discarding reset: `Bouncer`'s period is
  exact under localization where boundary-resolution firings would accumulate
  a full overshoot per reset, while `Sawtooth`'s carrying handler shows the
  complement — its boundary states are firing-resolution-invariant, so the
  boundary recursion stays its exact reference.
- **The θ = 0 validation discriminates the edge's cause, before any
  interpolant cost.** A slot write between runs is the frame-top epoch seam:
  the re-measured σ₀ holds under the frame's own `u`, so the edge is the
  drain's, no in-frame crossing exists, and the event falls through to fire at
  the frame top *exactly* — the indexed grid time, no budget spent, nothing
  warned. The gate idiom's `Bool` flip is the same story. And a crossing
  landing on the grid point itself degenerates the other way: every interior
  trial is not-holding, the localization is discarded, and the firing is
  bitwise the boundary-detected one.
- **Multiple crossings fire earliest-first, and ties need no decision.** Each
  triggered event is root-found, the boundary fires at the minimum `t*`, and
  the later crossing simply does not hold there — it re-triggers on the
  remainder against a fresh interpolant. A tie is one localization: both edges
  stand at the shared `t*` and fire inside that one boundary's iteration,
  which a budget of 1 proves by not warning.
- **The budget counts `t*` boundaries, and exhaustion degrades.** Per segment,
  root-finding runs are structurally bounded at one per declared event; the
  segment count is what chattering inflates, so `localization_budget` prices
  exactly that. The relaxer spends 8 localizations in one frame, warns
  `ChatteringBudget` once — naming the event and the count — and its next
  crossing fires at the frame top: boundary granularity for that frame alone,
  the remainder having already completed.
- **Ticks are never due at `t*`.** The `t*` boundary is the off-tick §10.6
  sequence — projection, the full iteration, priors settled — and nothing
  else: a sampled consumer beside a localized event accumulates only frame-top
  samples, where a spurious tick would have added the mid-frame value.
  Discrete cells ZOH-hold through the trials for the increment-3 reason:
  the interior sweep the trials run has no discrete entries to gate.
- **The integration method is a deployment binding, and nothing outside its
  own struct knows which one ran.** `method = Heun` swaps the backend with the
  loop, the localization machinery and every model unchanged, and the
  convergence orders fitted through that one surface — 4 and 2 against the
  same matrix-exponential reference — pin each backend as itself, where a
  shared loose tolerance would let a mislabeled one hide. The dense-output
  obligation rides the seam too: a localized stamp under Heun lands at the
  discrete solution's O(h²) through the same Hermite trials, and the
  frame-top stamps stay bitwise, having never depended on the method.
- **A staged write is pending, never applied.** `stage!` from any task touches
  no live slot; the batch lands at the top of the next frame `run!` advances —
  and nowhere earlier. A batch staged while stopped waits through boundary
  zero (`init!` runs on the un-drained slot — the contrast with `set_slot!`
  before `init!`, which makes a holding predicate fire at `t₀`), and the
  frame's outcome is a pure function of the drained batch: the staged
  trajectory is bitwise the directly-written one.
- **Merge is the only coalescing policy, and every check runs at staging.**
  Newest wins per face and untouched faces survive — a sparse `b` batch cannot
  clobber a pending `a` (§15.3's flaps/gear hazard), a re-staged face takes
  the newest level. The checks run on the writer's side: a face with no
  position in the schema is discarded under `OutOfClaimEntry`, an
  unconvertible value under `EntryTypeMismatch`, the rest of the batch stands,
  and the drain is pure application — the shim having already converted to the
  activation's slot types, so a `Dual` deployment stages plain reals. The
  empty drain is allocation-free: a quiet loop's frame top costs nothing.
- **Nothing reachable from a published snapshot is ever written again.**
  Publication is one buffer copy and one release-store after the boundary
  sequence completes — boundary zero included, off-tick frame tops included —
  and the snapshot is the whole table, root slots riding along as the source
  cells they are, state stores deliberately absent (§11.2). The run moves on;
  the snapshot holds. A concurrent reader acquire-loading `latest` sees only
  coherent worlds — the in-lockstep pair `g2 = 2·g1` holds bitwise in every
  observed snapshot, and `t` never decreases — and a task staging against
  running frames loses nothing to the CAS merge: the last staged level is what
  stands.
- **One writer per slot at any time, structurally.** Claims are disjoint by
  admission, so cross-writer races on a slot cannot arise and no drain-order
  arbitration policy exists to need. The admission's order is itself pinned:
  the same instance re-attached under an overlapping claim is
  `AlreadyAttached`, never a self-`ClaimConflict`, because identity runs before
  claims — which is what keeps `ClaimConflict` always naming two *distinct*
  devices. Identity is the instance (`===`): the library's stub devices are
  mutable structs precisely so two same-named `Pad`s are two devices.
- **One claim mechanism, two claim sources, exhausted at the attach point.**
  A returned claim is the enumeration, validated face by face
  (`AttachUnknownFace`) — the empty enumeration an honest may-write-nothing
  degenerate, not a back door. A computed claim is the unclaimed complement at
  the attach instant, never recomputed, so attaching the greedy claimant last
  is the idiom and a second greedy stakes the empty remainder under
  `EmptyGreedyClaim`. Downstream nothing tells the sources apart: the greedy
  entry's out-of-claim staging is the same `OutOfClaimEntry` as anyone's.
- **The harness register's surface is derived, never claimed.** It is the
  unclaimed complement, recomputed at every roster change: a harness write to
  a claimed face is `ClaimedFaceEntry` naming the incumbent, a rostered greedy
  empties the surface outright and every harness `stage!` in such a session is
  rejected that way (D-192), and detach regrows it. The recompilation's one
  seam is renormalization: a batch staged before a stopped-sim `attach!` is
  reshaped through the new schema — the newly claimed face discarded with the
  incumbent and the site named, the rest surviving to the drain — so the run
  always starts with cells matching the run's schemas.
- **The roster is frozen per run, and the drain stays free.** `attach!` and
  `detach!` against a running loop are `ServiceLifecycle`, from any task, and
  the freeze lifts with the run; device ids are monotonic per `Simulation`,
  never reused, and a rejected attach consumes none. The per-writer drain
  thunks are compiled at the same stopped-sim points, so a populated roster's
  empty drain still allocates nothing — and the frame's outcome stays a pure
  function of the drained batches, whoever staged them: the device-staged
  trajectory is bitwise the directly-written one.

**The store bundle is in, and with it the *plural* in "per-eltype stores".**
The bench that settled the representation (D-162) measured exactly one buffer;
`StoreBundle` now holds one `CellStore` per leaf element type, keyed by the
eltype's name. Selection is static — an address carries its port type, whose
leaf eltypes name the buffers at compile time, with one cursor per eltype in
address *fields* — so a deliberately pinned `Float64` leaf (D-166) keeps a
buffer of its own beside the `Dual` one instead of being flattened into it as
a zero-partial, which is what increment 2 had to refuse by name, and a
mixed-leaf cell (an `Int` tag beside `T` leaves, or a pin inside a declared
struct) simply spans several buffers, its leaves each living with their own
kind. The `K = 1` homogeneous cell is a case, not a special case; the C2M
point in `../cellstore_bench` (2026-08-20) re-confirmed D-162's flat compile
curve for exactly this address shape. The bundle is one concrete type per model, keeping `Chunk`'s
store parameter single: chunk-type count, not model size, is what bounds the
compile curve D-162 blessed.

Correctness is checked against analytically integrated references — the
continuous loop by matrix exponential, the sampled loop by its exact ZOH
discretization `q[k+1] = Ad q[k] + Bd s[k]`, localized event times against
analytic crossing instants — with a tolerance, never `==` (D-163), except for
the frame-top stamps, which are asserted bitwise against the indexed grid time
because that is the claim. The sampled-data reference is the sharp one: it
only matches if the hold is real *and* output stages run before updates within
a boundary.

## What is deliberately absent

**D-187's grid diagnostics:** the bound schedule here is plain data on the
`Simulation` (path, `D`, `Φ`, `Δt`), not the named printable artifact — no
hyperperiod-chart `show`-form, no anchor/provenance columns carried through, no
leave-one-out refinement factors or prime attribution, no nearest-non-refining
offset suggestions, no `GridUtilization` advisory. The refusal messages name
the offending anchor and the pool's GCD, and stop there.

Beyond that: computed connections and the §8.8 helpers (`input_passthrough` and
the generic-holding sugar), visibility enforcement (§8.3), did-you-mean
suggestion lists (an error names the offender plainly), auto-published ports,
§9.5's always-on conformance check (the return laws are checked once, at the
probe), §13.2's diagnostic framing (build errors here are a plain `BuildError`
with a good message, not the structured carrier, and runtime degradations are a
plain `@warn` carrying the spec'd payload, not Appendix C's named advisory
values), partial advance and the per-run overrides (§12.6 — `run!(sim, t_end)`
takes `t_end` to the nearest step boundary), and the runtime periphery beyond
increments 9–10's data plane:

- **the device task machinery (§11.6) and §11.1's task topology:**
  `init!`/`loop`/`shutdown!`/`unblock!`, the framework wrapper with its
  `DeviceCrash` discrimination, handles and their primitives, `should_abort`,
  `map_input`/`map_output` and the shipped `TableBinding`, the bad-datum
  channel, device death and orphaned claims — a rostered device here is an
  identity holding a claim and a staging cell, with no task behind it, and
  `run!` blocks its caller: with no device tasks to spawn, that is §11.1's
  calling-task register, not a deviation. The output side is absent wholesale —
  a binding declaring `is_output` is rejected at attach, `reads` and the
  compiled gather unbuilt;
- **the snapshot's framework status** (§11.2, §11.8) — the snapshot carries
  the table, `t` and the frame ordinal, no status value — along with §11.8's
  diagnostics and liveness machinery wholesale;
- **the log** (§11.2): retention, `log_every`/`log_max`, progressive
  re-decimation, the kept endpoints — readers hold snapshots or lose them;
- **output-device bindings and the selectors** (§11.2, §14.4), the §11.5
  input trace, and the §11.7 GUI write path;
- **all of §12 beyond `latest`'s inspection register:** control plane, loop
  scheduling, the next-snapshot wait, shutdown, real-time pacing (§10.7) and
  replay (§12.7). The §11.3 freeze is enforced by a single running flag —
  §12.6's status machine is absent, so `built`, `initialized` and `stopped`
  are indistinguishable here and every not-running point admits
  `attach!`/`detach!`.

## Stand-ins: where the prototype's shape is not the spec's

**Rule: nothing deviates silently.** Every construct a reader could mistake
for the design's is accounted for in exactly one of three places: the "what is
real" table, the absence list above, or a row here naming the spec shape it
replaces. The table is transactional — an increment that introduces a stand-in
adds its row in the same commit, and one that retires a stand-in deletes it. A
deviation in none of the three places is a defect in this file, not a liberty
the prototype has. `prototypes/` sits outside the design tools' checked
rosters, so no battery enforces this section: the diff review is the
enforcement.

| spec shape | stand-in here | retirement |
| --- | --- | --- |
| one published snapshot per boundary, `t*` boundaries included (§11.2, §10.6) | publication at frame tops and boundary zero only; a `t*` boundary fires unpublished | when a consumer needs the per-boundary record — the log or the trace increment |
| slot initial values owned by the init/trim services (§11.3, §14.6); the running write path is the drain alone | `set_slot!` writes a root slot directly at any stopped-sim point | the §14 services increment |
| a device stages through the handle `attach!` returned, `stage!(handle, batch)` on its own task (§11.6) | `stage!(sim, dev, pairs...)` addresses the rostered device's cell directly, and `attach!` returns the device id | the device-contract increment, when handles exist |

(Every stand-in introduced through increment 5 was retired on 2026-08-20, and
increment 6's one row — localized guards detecting at boundary resolution —
was retired by increment 7 on 2026-08-25. The first two rows above entered
with increment 9, the third with increment 10.)

## Authoring caveat found while building this

Declarations written in a **local scope** never reach the framework. Inside a
`let`, a function body or a `@testset`, `h_x(::MyComp, (; x)) = …` does not add
a method to the global `h_x`: it binds a *new local function* of that name.
Calls inside the block see it and look fine; the global generic function the
build dispatches on is untouched, so `build` sees a component that declares
nothing at all. Ordinary top-level definitions — the
script, module and notebook-cell cases — are unaffected. The files under
`test/` define their malformed test components at top level for exactly this
reason.

This is the local-scope sibling of the `using Flight` trap §8.1 already
documents, and §8.1's shadowing check does not reach it: there is no
parent-module binding to compare against, because the shadow is a local binding
that disappears with its block. Ratified as D-164: a component that declares
nothing and defines no stage is a build error — an inert component cannot be
intentional — spec'd as the inert-component check in §8.1's stage register,
raised as `DeadStage` at the probe (§9.3). Increment 4 catches the case where
the trap lands, one stratum earlier and under the other rule: a component that
declares nothing has no *class* to read either (§8.5), so the build names both
families rather than reaching the probe at all. `DeadStage` itself is not built.

The periphery's trait and enumeration methods (`is_input`, `is_greedy`,
`claims`, `needs_calling_task`) are the same class of declaration and hit the
same trap — a trait "declared" inside a `@testset` is a local function the
conformance check never sees, presenting a binding that declares nothing. The
malformed bindings in `test/test_roster.jl` sit at top level for exactly this
reason.
