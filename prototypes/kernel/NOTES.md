# Kernel prototype — increment history and test-property notes

The optional companion to `README.md`: the per-increment build narrative and
the property-by-property record of what the tests pin down. Read it when
modifying an existing test or wondering why one asserts what it does; nothing
here is needed to orient a coding session — that is the README's job. The
accretion convention follows the content: an increment adds its narrative
paragraph and its property bullets *here*, keeping the README constant-size.

## The increments

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

**Increment 11 — the log.** Publication moves to its §11.2 cadence — every
boundary, `t*` boundaries included, each snapshot released before integration
resumes — and the log dissolves into it: a vector of retained snapshot
references, the same objects with zero extra copies, behind the
`log`/`log_every`/`log_max` keywords with progressive re-decimation — the
stride doubling at each fill, one release per retained append, compaction once
per generation — and the two endpoints held unconditionally outside the bound,
the terminal one standing in as the latest published boundary while §12.4's
run end is absent. The keywords are view policies, never
trajectory-determining, and increment 9's publication stand-in retires.

**Increment 12 — the device contract: tasks, the handle, and the run's end.**
§11.6's authoring contract in full — `init!`/`loop`/`shutdown!`/`unblock!`
with no-op resource defaults and an error-throwing `loop` fallback, the
framework wrapper with its `DeviceCrash` catch, `shutdown!` guaranteed on
every exit path and `should_abort` (now an `attach!` keyword) consulted at
departure — behind the handle `attach!` constructs and returns: the
capability carrier, read/stage/control and never the `Simulation`. `run!`
takes §11.1's topology — the §12.4 pre-spawn init bracket in attachment
order, one run-scoped task per live entry, the loop moving to a spawned task
when a `needs_calling_task` holder is rostered while the calling task runs
that body inline through the identical wrapper — and §12.4's tail closes
every run: stop observed at frame top, final snapshot, sticky stopped, waits
woken, `unblock!`, the join under `join_timeout` (the deployment keyword
D-198 promoted, one shared deadline) with `DeviceJoinTimeout` abandonment.
Beneath it, exactly the minimal §12 slices: §12.1's stop word alone, §12.3's
counter-plus-condition wait with the normative publication order and the
boundary ordinal riding in the snapshot, §12.2's per-frame yield with devices
rostered. The device-staging stand-in retires; the log's terminal endpoint is
now §12.4's run-end snapshot in fact.

**Increment 13 — the binding's runtime half.** §11.6's mapping conventions
with their framework-owned instances: `TableBinding` — the table in the type,
the claim derived from it, deadzone/expo conditioning in the generic
`map_input`, the shared pure helper with an owner — behind the handle's
`binding` capability; the bad-datum channel — `report!(handle,
MalformedDatum(cause))` into the per-device diagnostic cell, §11.8's ring of
sixteen plus suppressed counts under the sentinel-swap drain at frame top and
run end, its presentation standing in as loop-side warns and the per-run
totals behind `diagnostics(sim)`; and the output side — §14.4's table
selectors, `reads` fixed as a labeled NamedTuple resolved at attach and
compiled to the one gather `gather(handle, snap)` runs, `map_output`
receiving exactly its NamedTuple, and the conformance check completed over
both (trait, method) pairs. Two stand-ins enter the table.

**Increment 14 — the diagnostic channel and the framework status.** §11.8
completed over §13.2's runtime stream: the closed kind set as types over the
fixed-shape isbits `KindCounts`, one cell per writer — the harness register
gaining its own beside the devices' and the loop's — the heartbeat field
stored by the handle primitives and judged against §12.2's 2 s threshold,
`DeviceCrash` written through the cell (the wrapper's catch on the device
task, the init bracket's on the calling task, §12.4's entry-addressed
report), the budget degradations on the loop's cell, the frame-top fold into
loop-owned accounts, and the `FrameworkStatus` every snapshot carries — the
delta riding exactly one snapshot, totals monotone since the run began,
`task_state` off the run's task registry (D-193). The two §11.8 presentation
stand-ins retire; three rows enter: the harness cell (a spec-enumeration
gap), the tail window (`DeviceJoinTimeout` and the post-final-frame-top
residue present synchronously, no snapshot following them), and the status's
per-publication vector against the spec's inline-allocation claim.
`diagnostics(sim)` is deleted, its account subsumed by the status; the
warning-matching tests convert to the spec's own doctrine — kind plus
payload fields, never message text.

## The properties the tests pin down

Each of these is a spec claim rather than a programming convenience:

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
- **Every boundary publishes, and the log is the publications themselves.** A
  localized firing's `t*` boundary releases its own snapshot before
  integration resumes — the log holds a snapshot whose `t` is bitwise the
  stamped `t*`, carrying the settled post-transition table while the boundary
  before it still shows the armed value — and `logged(sim)[end] ===
  latest(sim)`: retention is reference-keeping over what publication already
  built, zero extra copies (D-023). What a localizing frame allocates is
  exactly that publication, the framework-side §7.5 carve-out, and nothing of
  the localization machinery's own.
- **Retention is a view policy under a continuously-held bound.** Deployments
  differing only in `log`/`log_every`/`log_max` produce bitwise-identical
  trajectories. The middle never exceeds `log_max` at any boundary — one
  release pairs with each retained append — and when the log fills the stride
  doubles: coverage stays global at `log_every · 2^k`, the middles landing on
  consecutive multiples of the effective stride with the whole run spanned,
  while the boundary-zero and latest snapshots outlive any policy, outside
  the bound. The off switch retains nothing; publication is upstream of it,
  so `latest` still sees every boundary.
- **The handle is the device's whole world.** It carries read, stage and
  control — never the `Simulation` — and a device staging through it from its
  own spawned task changes nothing the drain can see: the frame's outcome
  stays a pure function of the drained batches, and a device-requested stop
  truncates the run at a frame top with the batch it staged still pending,
  applied at the next run's first drain and nowhere earlier.
- **The bracket holds on every exit path.** Voluntary return, a crash, and a
  failed `init!` all run `shutdown!` exactly once — the crash warned as
  `DeviceCrash` with the run continuing, the failed `init!` spawning no task
  and never reaching `loop` — and in every case the claims persist to run
  end: death is not detach, and the harness register still cannot take a dead
  device's face. `should_abort` splits the run's disposition uniformly: a
  flagged crash stops the run, a flagged `init!` failure leaves the stop
  already pending and the run advances zero frames.
- **The run's end is a protocol, not a return.** The stop word is observed at
  frame top only, the final snapshot goes out before the sticky status is
  set — so a body observing `running(handle)` false can still read a complete
  final world — the waits wake, `unblock!` turns a parked blocking call into
  a clean exit inside the cap, and a body that ignores the predicate is
  abandoned under `join_timeout` with `DeviceJoinTimeout` naming it, `run!`
  returning regardless. Abandonment is not a kill: the straggler's wrapper
  still runs `shutdown!` when it finally returns. The word clears at the next
  run's top, which owes nothing to the last one's stop.
- **The calling task is pinned by trait, and the loop is the movable piece.**
  A `needs_calling_task` device's body runs on exactly the task that invoked
  `run!`, inside the same wrapper as any spawned loop, while the frame loop
  runs spawned — and the trajectory is bitwise the deviceless one, topology
  being scheduling, never semantics. A device with no `loop` method crashes
  by name at its first spawn instead of idling as attached-but-inert.
- **A waiter never wakes onto a stale snapshot.** The release-store of
  `latest` precedes the counter increment under the condition's lock — the
  §12.3 normative order — so the boundary ordinals a waiting consumer
  observes only ever advance, each snapshot indexes itself, and the stop wake
  hands the waiter the final world rather than nothing. While stopped the
  wait returns at once; only a running loop can park it.
- **Conditioning runs upstream of the face, and the face never knows.** A
  face carries post-conditioning semantics (§11.4's GUI-parity test):
  `TableBinding`'s deadzone/expo run in `map_input` on the device task, so
  the slot receives the conditioned level and the model output is bitwise the
  directly-staged one — while an entry declaring no parameter passes its
  value through untouched, which is what carries a throttle's `[0, 1]` level
  or a press counter without ever meeting the axis convention. The claim is
  the table: what the binding may write is exactly what its entries name.
- **A bad datum and a bug have two fates, and the bound is the rate limit.**
  The loop's tolerance idiom — catch its own parser error, stage nothing,
  `report!`, continue — keeps the link alive with the stream's good datums
  landing, while an unclassified throw is the wrapper's `DeviceCrash`; the
  classification is the author's, because only they know their parser. A
  twenty-report flood costs sixteen retained values — earliest-in-frame, the
  ones with diagnostic content — plus one integer, and nothing is lost by
  not looking: the status's totals carry the full account (since increment
  14; the per-run totals behind `diagnostics(sim)` before it), the tail
  sweeping the cells one last time at the run's end.
- **Binding drift fails at attach, never on the wire.** Every `reads`
  selector resolves against the build at the attach point — the unknown
  path, the unknown slot, the input face offered to `get_face` — so
  `map_output` receives exactly the compiled gather's labeled NamedTuple or
  the device never rosters, and a rejected attach consumes nothing. The two
  read registers agree by construction: an exported face aliases its
  producer's cell, so `get_face(:y)` and the deep `get_output` are bitwise
  equal in every observed snapshot. An output-only binding stakes no claim,
  and the harness register keeps every face.
- **The delta rides exactly one snapshot; the totals ride every one.** A
  drained occurrence appears in `recent` in precisely the first snapshot
  published after its frame-top fold — `[0, 1, 0, 0, 0, 0]` across a
  six-boundary log — while `totals` is monotone from that boundary on. That
  is §11.8's any-cadence legibility made a vector equality: a 60 Hz reader
  sees each occurrence once, an occasional sampler still reads the complete
  account, and decimation loses *which* boundary, never *how many*.
- **A diagnostic is accounted exactly once, wherever timing lands it.** A
  device-task report races the run's last frame top, so the tests assert the
  exclusive-or: in the terminal status's totals, or presented by the
  run's-end sweep — never both, never neither (`accounted` in `utils.jl`).
  The deterministic corners pin the two pure cases: a failed `init!` reports
  pre-spawn and always makes the first fold; a `should_abort` init failure
  runs zero frames and can only be swept.
- **Dead-from-boundary-zero needs no machinery.** The failed-`init!` device's
  record shows the never-stored heartbeat (`0.0`, stale against any clock)
  and `task_state === :none` — §12.4's rule that the cell's silence *is* the
  marking, asserted off the published status alone.
- **Staging diagnostics wait like staged entries.** A rejection at stopped-sim
  staging sits in the writer's cell — readable there, `(@atomic cell.batch)`
  in the tests — and surfaces in the next run's first status, the same
  lifecycle as the batch it was rejected from; the attach renormalization's
  `ClaimedFaceEntry` carries `site = :renormalization` to tell it from an
  ordinary staging rejection.

## Stand-in retirement history

Every stand-in introduced through increment 5 was retired on 2026-08-20;
increment 6's one row — localized guards detecting at boundary resolution —
was retired by increment 7, and increment 9's publication row by increment 11,
both on 2026-08-25. The `set_slot!` row in the README's table entered with
increment 9; increment 10's device-staging row was retired by increment 12 on
2026-08-26, `attach!` returning the handle, and the two §11.8 presentation
rows entered with increment 13. Increment 14 retired both presentation rows
on 2026-08-26 and entered its own three: the harness register's cell, the
tail window, and the status's per-publication allocation.
