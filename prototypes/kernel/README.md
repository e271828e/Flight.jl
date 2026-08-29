# Kernel prototype

The walking skeleton for the framework in `docs/notes/design/framework_spec.md`,
built to keepable standards and grown one increment at a time — increments 2–21
so far (increment 1, the cell-store bench, is frozen in `../cellstore_bench`;
D-162 cites its numbers). This file is the orientation; read it first and
alone. On demand:

- `MAP.md` — what each source file implements, piece by piece, with spec
  citations: the long form of the table below, the long absence list, and the
  long form of the authoring trap.
- `NOTES.md` — the per-increment narrative, the property-by-property record of
  what the tests pin down, and the stand-in retirement history. Read it when
  modifying an existing test or wondering why one asserts what it does; new
  increments add their paragraph and property bullets there.

    julia --project=. test/runtests.jl

## What is real here

| file | implements | spec |
| --- | --- | --- |
| `src/leaves.jl` | the leaf walk: flatten / reconstruct / the activation retype | §7.1, §7.2 |
| `src/diagnostics.jl` | the diagnostic kinds, `severity`/`path`/`message`, the `BuildError` carrier and its compiler-style rendering, `logged`, `InternalInvariant` | §13.1, §13.2, Appendix C, D-214, D-215 |
| `src/declare.jl` | the declaration layer: both tiers' name families and arities, the bundle law, `probe_value`, the connection declarations beside `transparent_container`, the rate registers with `sample_times`, the event surface | §5.2, §8.2, §8.5–§8.7, §9.3, D-211 |
| `src/assembly.jl` | class by declaration shape; children and containers (bare-key transparency and its three-arm collision family); paths, §6.1's one-level rule, endpoint and face resolution, the root's face invariants; the flatten pass with its two-sided face graph and the sample-time fold; §13.3's `resolve`/`resolve_terminal`/face-list primitives and §8.8's `input_passthrough`/`output_passthrough` | §6.1, §8.5–§8.8, §9.1, §9.2, §13.3, D-207–D-212 |
| `src/store.jl` | per-eltype cell stores, the `StoreBundle`, gather/scatter, `_cell_key`, the `Clock` | §9.7, D-162 |
| `src/executor.jl` | entries, the chunked unrolled walk, the interior/boundary split, the `(idx − Φ) % D` gate and boundary zero's `ESTABLISH` beside it, the event set with its registers and the guard/fire/project walks | §9.7, §10.4–§10.6, §14.5, D-205 |
| `src/build.jl` | tier classification, the probe, the feedthrough graph, the layout, embed-accept, the `Build` and its activations, deployment binding, and `compile` → `Executor{T}` — one activation's buffer set with the bodies closed over it, one owner per set, `evaluate!`/`_round!`/`apply!` on it | §8.2, §9.1–§9.4, §9.7, §10.4, D-166, D-208, D-210 |
| `src/readers.jl` | the closed read-selector family, `reads`, the internal `_compile_reads` → `Reader{T}`, `gather` as `apply!`'s twin over an executor; activation identity on readers and plans as an internal invariant | §14.4, §14.7, §14.10 |
| `src/sim.jl` | `Simulation` (owning its `exec`), the deployment keywords, the boundary macro-sequence and event phase, `init!`, `run!`/`step!`, `attach!`/`detach!`, staging/drain/publication, the lifecycle and termination record, the accessors | §10.2–§10.6, §11.1–§11.4, §11.8, §12.1–§12.6, §13.5, §13.6, §14.5, §14.6, D-203 |
| `src/stepper.jl` | the seam's backend side: RK4 and Heun, the retained `startpoint`, dense output | §10.2, D-017 |
| `src/localize.jl` | the frame loop: arrival sweep, θ = 0 validation, ITP bracketing, `t*` boundaries, the localization budget | §10.4, D-018, D-133 |
| `src/dataplane.jl` | the compiled writer and staging cells, the drain, snapshots and the log with re-decimation, the typed diagnostic kinds and cells, the framework status | §11.1–§11.4, §11.8, §12.6, D-137 |
| `src/roster.jl` | device/binding traits and conformance, the roster, both claim sources, the harness register | §11.3, §11.4, §11.6 |
| `src/bindings.jl` | `TableBinding`, `map_input` and the conditioning helper, binding reads resolved at attach (`ReadBindingUnresolved`, the source rule) | §11.2, §11.6, §14.4 |
| `src/devices.jl` | the device contract, the handle, the task wrapper, the init bracket and the tail under `join_timeout` | §11.1, §11.6, §12.1–§12.4, D-198 |
| `src/conditions.jl` | the condition algebra, one collecting pass behind both application registers — `resolve_condition` (values) and `compile_plan` (`Getter{P}` lenses, `SpecializedPlan`, `ConditionShapeDrift`) — root-input totality, `capture` | §14.1–§14.6, §13.1, D-063–D-068, D-204, D-205 |
| `src/trim.jl` | `TrimProblem`, the `solve` seam with `LevenbergMarquardt`, `trim!` over D-213's two-half scratch world, `TrimReport`, the `Trim*` kinds | §14.7, §14.8, D-213 |
| `src/library.jl` | the coverage component set, `Group` and the named assemblies, the devices and bindings, the `condition` fragment-function idiom, `Pendulum` — user material | — |

Correctness is checked against analytically integrated references with a
tolerance, never `==` (D-163) — except the frame-top stamps, asserted bitwise
against the indexed grid time because that is the claim.

## What is deliberately absent

The long form, with reasons, is in `MAP.md`. In brief:

- **The Appendix C kinds whose mechanism is absent** — an absence gets no
  struct, so no `ThreadBudget`, `DeadStage`, `BundleFieldError`,
  `UserCodeFraming`, `UnboundedRun` or replay kind is defined here, and
  `TapResolution` is raised by the read register alone, never by §14.10's
  absent tap register. Absent with them: did-you-mean **ranking** (the
  candidate list a site holds is carried and rendered; nothing orders it by
  edit distance, and a mistyped *path* gets no list at all) and §11.8's
  maxlog renderer.
- **§9.5's always-on conformance check** (the return laws are checked once,
  at the probe); **§8.3 visibility**; auto-published ports; §13.3's
  load-bearing generic-holding check.
- **§8.8 beyond the helper pair** (the feed-list idiom, generic-holding sugar,
  required-faces declarations); **D-187's grid diagnostics** (the bound
  schedule is plain data; refusals name the anchor and the pool's GCD).
- **§14**: `linearize` (§14.10), mounting (§14.9), the NLopt fallback and the
  nominal-activation loop it would run on; sub-port-field addressing; index
  addressing in the binding register.
- **§11.5's input trace**, **§11.7's GUI write path**, §10.7 pacing and its
  diagnostics, the §11.8 remainder (`DebtReanchor`, `ThreadBudget`,
  `ReplayDiscardedStaging`, `UnboundedRun`, the maxlog renderer).
- **§12 beyond its built slices**: pause and the control plane's surface, the
  operator interrupt, replay (§12.7); §13.4's `StepError` wrap, cursor and
  nonfinite sweep. `run!` requires a finite `t_end`; every non-running state
  admits `attach!`/`detach!`.

## Stand-ins: where the prototype's shape is not the spec's

**Rule: nothing deviates silently.** Every construct a reader could mistake
for the design's is in exactly one of three places: the table above, the
absence list, or a row here naming the spec shape it replaces. Transactional:
the commit introducing a stand-in adds its row, the one retiring it deletes
it. No tooling enforces this (`prototypes/` is outside the design tools'
rosters); the diff review is the enforcement.

| spec shape | stand-in here | retirement |
| --- | --- | --- |
| the per-writer status rides inline in the snapshot's one per-boundary allocation — zero additional heap allocation on a quiet frame (§11.8) | a `Vector` of per-writer records built at each publication, the small extra allocation the simple shape costs | an allocation-tightening pass (an `NTuple` status type fixed per run) |

## Authoring caveats

**Declarations in a local scope never reach the framework.** Inside a `let`,
a function body or a `@testset`, `h_x(::MyComp, (; x)) = …` binds a *new local
function*, not a method of the global `h_x`, so the build sees a component
that declares nothing. Test fixtures — components, devices, bindings, traits —
live at top level for this reason (long form, and its D-164 ratification, in
`MAP.md`).

Traps hit more than once while building, for whoever builds next:

- every `src/` and `test/` file is included into `Main`, so a test fixture
  named like a framework type silently clobbers it — `rg "struct <Name>"
  test/` before naming one;
- `===` has no curried form (`all(===(x), v)` fails — use a lambda); a
  `where`-clause method's `.sig` is a `UnionAll` (`Base.unwrap_unionall`
  before `.parameters`);
- a local named `events` inside `compile` shadows the `events(c)` accessor;
- assert a store's type on the `Ref` — `(v[ci]::Base.RefValue{S})[]` — never
  after `[]`, which boxes 16 bytes;
- `Symbol(::Type)` printing depends on the printing module; key buffers with
  `_cell_key`.
