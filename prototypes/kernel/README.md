# Kernel prototype

The walking skeleton for the framework in `docs/notes/design/framework_spec.md`,
built to keepable standards and grown one increment at a time. Increment 1 (the
cell-store representation bench) lives in `../cellstore_bench` and stays frozen
there — D-162 cites its numbers.

Increments 2–16 are built: the two tiers, hierarchy and assemblies, the
multi-rate grid, events with localization, the stepper seam (RK4/Heun), the
data plane's core exchange, the roster and claims, the log, the device
contract with its tasks and the run's end, the binding's runtime half, the
diagnostic channel with its published framework status, and the run
lifecycle with its termination policy and record.
The per-increment narrative and the property-by-property record of what the
tests pin down live in `NOTES.md` — optional reading, for when you are
modifying an existing test or wondering why one asserts what it does; new
increments add their paragraph and property bullets there.

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
| the data plane's core exchange: the compiled writer — one schema/shim/merge/scatter unit over any write surface, with every check at staging and the out-of-schema warning discriminated by writer — staging cells, the frame-top drain via `atomicswap` behind stopped-sim-compiled thunks, `run!`'s frame anatomy (drain → integrate → boundary sequence → publication), and publication — the boundary-consistent whole-table snapshot with the state stores excluded and the framework status riding beside the table, `capture`'s buffer copy, the release/acquire `@atomic latest` pair, and `latest(sim)` as the inspection register | §11.1, §11.3, §11.4, §11.2, §12.6 | `src/dataplane.jl`, `src/sim.jl` |
| the roster and claims: `AbstractDevice`/`AbstractBinding` with the declared sides and the false root defaults, the enumeration contract and the bidirectional conformance check (`which` against the fallback), `attach!`/`detach!` behind the §11.3 freeze with the three-part admission in spec order, both claim sources with `EmptyGreedyClaim` for the staked empty remainder, monotonic never-reused device ids, per-device cells over claim sets, the harness register as the derived remainder — recompiled at every roster change, its pending batch renormalized (`ClaimedFaceEntry` at the seam), emptied outright by a rostered greedy (D-192) — and the attachment-order drain, harness last | §11.3, §11.4, §11.6 | `src/roster.jl`, `src/sim.jl` |
| the log: publication at every boundary (`t*` included, before integration resumes), retention as a vector of the published references themselves under `log`/`log_every`/`log_max`, progressive re-decimation (the stride doubling at each fill, one release per retained append, once-per-generation compaction, the middle at consecutive multiples of `log_every · 2^k`), the two endpoints outside the bound — the terminal one being §12.4's run-end snapshot, published before the sticky status — and `logged(sim)` as the stopped-sim reader behind the §11.3 gate | §11.2, §12.4, D-137, D-023, D-038 | `src/dataplane.jl`, `src/sim.jl` |
| the device contract and its tasks: the four contract functions (no-op resource defaults, the error-throwing `loop` fallback), the handle as the capability carrier — `stage!`, `latest`, `wait_next_snapshot`, `running`, `stop!`, `binding`, `gather`, `report!`, never the `Simulation` — the wrapper (crash → `DeviceCrash`, `shutdown!` on every exit path, `should_abort` at departure), §11.1's topology (run-scoped spawn per live entry, the movable loop, the inline calling-task body, topology derived after initialization), the §12.4 init bracket and tail (1)–(5) under the shared `join_timeout` deadline, and the §12 slices beneath them — the §12.1 stop word, §12.3's wait (counter mirrored under the condition's lock, release-store before increment, monotonic across runs, the ordinal in the snapshot) and §12.2's per-frame yield | §11.6, §11.1, §12.4, §12.1, §12.2, §12.3, D-198 | `src/devices.jl`, `src/sim.jl` |
| the binding's runtime half: `TableBinding` — the table in the type, the validating constructor, `claims` derived from the entries — with its generic `map_input` (sparse over the datum, an unknown channel a deliberate throw) and the conditioning helper (deadzone rescale, expo blend `(1−e)·a + e·a³`, endpoints fixed, pass-through where no parameter is declared), the loop-idiom conventions the framework never calls, §14.4's table selectors as deferred-read values (`get_output`/`get_slot`/`get_face`, whole cells), `reads` fixed as a labeled NamedTuple, resolved at attach (`ReadBindingUnresolved`, the two root registers discriminated) and compiled to one gather — the compiled scatter's mirror — behind `gather(handle, snap)`, the conformance check completed over both (trait, method) pairs with the error-throwing `reads` fallback, and the output-only attach staking no claim | §11.6, §11.4, §11.2, §14.4 | `src/bindings.jl`, `src/roster.jl`, `src/devices.jl` |
| the diagnostic channel and the framework status: the closed kind set as types — `MalformedDatum` (`report!`'s one author-facing kind, the payload the cause, attribution the cell's), the three staging kinds written at staging on the writer's task, the two budget degradations on the loop's own cell, `DeviceCrash` from the wrapper and the init bracket — over the fixed-shape isbits `KindCounts`; one cell per writer (each rostered device's, the harness register's, the loop's own) with the ring of `DIAG_RING` = 16 plus per-kind suppressed counts (earliest-in-frame retained), the CAS append mirroring the staging cell's, and the heartbeat field the handle primitives store (`stale` against the 2 s threshold, initial `0.0` = dead-from-boundary-zero); the frame-top fold into per-writer accounts, the published `FrameworkStatus` every snapshot carries — `recent` riding exactly one snapshot, `totals` monotone since the run began, `heartbeat` and `task_state` (D-193, off the run's task registry) fresh per publication — and the run's-end sweep behind the tail | §11.8, §13.2, §11.2, §12.2, §12.4, §11.6 | `src/dataplane.jl`, `src/roster.jl`, `src/devices.jl`, `src/sim.jl`, `src/localize.jl` |
| the run lifecycle and termination: the five-state machine behind `lifecycle(sim)` — mandatory `init!` as the one door into `initialized`, opening the trajectory wholesale (log, stop word, termination record and every staged batch clear with it), the `:running` state as the §11.3 freeze spanning the tail, terminal `stopped`/`errored` — the §13.5 policy pair `t_end`/`stop_on` as constructor defaults with `run!`'s per-run overrides, validated identically at both sites (stop faces are root-exported Bool output faces, OR-combined, sampled at every publication — a `t*` hit makes that snapshot final and abandons the frame's remainder), partial advance (`step!` by `frames`/`t_plus` over the same frame loop `run!` drives, deviceless, returning the count actually advanced), the termination record behind `termination(sim)` naming the source (`:t_end`, `:stop_on` with the face, `:control_stop`, `:error` with the cause), and §13.6's abnormal entry — the failed boundary discarded by publication-last construction, the previous snapshot promoted, stores retained for post-mortem, the ordinary tail, the synchronous rethrow | §12.6, §13.5, §13.6, §12.1, §13.4 | `src/sim.jl`, `src/devices.jl`, `src/localize.jl` |
| the coverage set: `Plant`, `Gain`, `Sum`, `DiscreteIntegrator`, `TickCounter`, `Smoother`, `WorkGain`, `ModedSource`, `ZOH`, `Ramp`, the event set (`Trigger`, `Follower`, `Sawtooth`, `Rotor`, `Chatterer`, `TwoShot`, `Preempted`), the localization set (`Stamper`, `GatedStamper`, `Bouncer`, `Relaxer`), the termination set (`Overload`, the §13.5 touchdown archetype; `Exploder`, §13.6's specimen), plus `Group`, the named `SampledLoop`/`Vehicle`/`MultiRate` assemblies, and the periphery set — `Pad`/`Panel` devices (mutable, `===` being identity, immediate-return loop bodies), the `Enumerated`/`Greedy` bindings, one per claim source, and `Readout`, the output side's coverage binding with the identity `map_output` | §5.2, §6.2, §7.3, §8.5, §10.4, §10.5, §10.6, §11.3, §11.6, §11.2 | `src/library.jl` |

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
probe), §13.2's build-side diagnostic framing (build errors here are a plain
`BuildError` with a good message, not the structured carrier — the *runtime*
warning stream does carry Appendix C's kinds as typed values, through the
§11.8 cells), and the runtime periphery beyond increments 9–16's data plane,
device tasks, diagnostic channel and lifecycle:

- **§11.8's remainder:** the pacer diagnostics are absent with §10.7, and the
  kinds whose sources are absent — `DebtReanchor`, `ThreadBudget`,
  `ReplayDiscardedStaging`, `UnboundedRun` — are absent with them; the
  presentation-side maxlog-25 renderer (count-only display past 25 cumulative
  occurrences per writer × kind) is presentation policy nothing here renders,
  the status being read raw by the tests. The snapshot still carries no
  run-status value: the lifecycle is read off the simulation
  (`lifecycle(sim)`), the framework status stays the diagnostic account. Two unguarded edges remain: staging through a handle
  whose device was detached lands in an orphaned cell and is silently lost
  (handles are run-scoped task equipment; a guard would put a roster scan
  back into `stage!`), and an `InterruptException` in a device loop reports
  as `DeviceCrash`, the discrimination belonging to the absent operator
  interrupt;
- **the selector family beyond the binding client (§14.4):** the store
  selectors `get_state`/`get_deriv`, whose holders — trim's reads,
  linearization's taps, `capture` — are §14's and absent with it; sub-port
  field and index addressing (`get_output(path, field, i)`: a selector here
  reads whole cells, as every reader of this table does); and did-you-mean
  candidate lists at resolution, per the general rule above;
- **the §11.5 input trace and the §11.7 GUI write path**;
- **§12 beyond increments 12 and 16's slices** (the §12.1 stop word, §12.2's
  yield, §12.3's wait, §12.4's bracket and tail, §12.6's lifecycle and
  partial advance): pause and the rest of the control plane's surface,
  real-time pacing (§10.7) with §12.2's thread-budget warning, the operator
  interrupt and its sigint masking (the runs here are test-driven, and the
  interrupt's termination-record tag is absent with it), and replay (§12.7)
  with `replay!`'s alternative entry into boundary zero. §13.4's `StepError`
  wrap and execution cursor are absent with the trace — the errored record
  retains the raw cause — and so is the nonfinite sweep; the unbounded
  register is absent with `UnboundedRun`, `run!` requiring a finite `t_end`
  from one of its two binding sites. Every terminal or not-yet-run state
  admits `attach!`/`detach!`, `errored` included: the freeze is exactly the
  lifecycle's `:running`.

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
| slot initial values owned by the init/trim services (§11.3, §14.6); the running write path is the drain alone | `set_slot!` writes a root slot directly at any stopped-sim point | the §14 services increment |
| the per-writer status rides inline in the snapshot's one per-boundary allocation — zero additional heap allocation on a quiet frame (§11.8) | a `Vector` of per-writer records built at each publication, the small extra allocation the simple shape costs | an allocation-tightening pass (an `NTuple` status type fixed per run) |

(The retirement history of past stand-ins is in `NOTES.md`.)

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

The periphery's trait and enumeration methods (`is_input`, `is_output`,
`is_greedy`, `claims`, `reads`, `needs_calling_task`), the device contract's
four functions (`init!`, `loop`, `shutdown!`, `unblock!`) and the mapping
conventions (`map_input`, `map_output`) are the same class of declaration
and hit the same trap — a trait "declared" inside a `@testset` is a local
function the conformance check never sees, and a `loop` "declared" there
leaves the global fallback in place, crashing the device by name. The
malformed bindings in `test/test_roster.jl` and `test/test_bindings.jl` and
the devices in `test/test_devices.jl`, `test/test_bindings.jl` and
`test/test_diagnostics.jl` sit at top level for exactly this reason.
