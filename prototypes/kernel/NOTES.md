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
assembly's own input faces are the root inputs.

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
tests' root-input writes move onto the machinery the θ = 0 discriminator was
built expecting. Two stand-ins enter the table.

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

**Increment 15 — the values-plus-mask batch.** D-202 applied: the staged
batch becomes a pair of positional tuples — a concrete values tuple and a
`Bool` touched-mask — replacing the `Union{Nothing,T}` carrier, whose
covariant narrow batch types made the frame-top scatter dispatch per
touched-face combination (a boxed call per populated drain, a mid-run
compile per never-seen pattern) and whose merge fell off Base's small-tuple
`map` path past 32 faces. The writer gains its `blank` batch — placeholders
drawn from the layout's probe values under an all-clear mask — the shim sets
the mask where it used to fill `nothing`, merge and scatter unroll by
generation over the writer's one concrete type, and the drain is
allocation-free populated as well as empty, at any width, whatever the
batch touches. No stand-ins enter or retire.

**Increment 16 — the run lifecycle and termination.** §12.6's five-state
machine behind `lifecycle(sim)`, §13.5's stop policy and termination record
behind `termination(sim)`, §12.6's partial advance, and §13.6's abnormal
entry, in one increment. The single `plane.running` flag retires: the §11.3
freeze becomes the lifecycle's `:running` — spanning `run!`, tail included —
while `ctl.stopped` stays the device-facing sticky word at tail step (1).
`init!` becomes mandatory and the one door into `initialized`, opening the
trajectory wholesale: the log, the stop word, the record *and every staged
batch* clear (§12.6's no-stale-batch rule — the pre-run register is now
`init!` → `stage!` → `run!`, which restructured the increment-9-era
stage-before-`init!` tests). `run!` trades its positional `t_end` for the
§13.5 keyword pair — constructor defaults, per-run overrides, both sites
validating identically, a finite bound owed from one of them (the unbounded
register stays absent with `UnboundedRun`). Stop faces are sampled at every
publication, `t*` included: a `t*` hit abandons the frame's remainder and
leaves the clock at `t*` — a real bug caught by the increment's tests,
`frame!` having stamped the frame top unconditionally. `step!(frames`/
`t_plus)` shares the one frame loop with `run!` — bit-identity by
construction — and returns the frames actually advanced; a stepping session
is deviceless and reports `initialized` between calls, which is what
absorbed the suite's run-chaining idiom. A loop-side throw leaves terminal
`errored` — the cause retained raw, §13.4's `StepError` wrap and cursor
absent — with the previous snapshot already final by publication-last
construction, and `run!` rethrows after the tail. Two decisions worth
flagging: a stop word pending at `step!`'s entry stops the session as
`:control_stop` rather than being cleared (`run!` clears at its top, `init!`
with the trajectory — the spec is silent on the stepped case); and
test_bindings' `Poller` became confirm-then-stop — stage, observe the write
applied in a snapshot, then stop — because a departing device cannot
otherwise guarantee its batch drained before the run ends. New coverage:
`Overload` (§13.5's touchdown archetype) and `Exploder` (§13.6's specimen).
No stand-ins enter or retire.

**Increment 17 — D-203: the typed termination record and the recorded tail.**
The spec amendment applied the day it landed (commit 81c3cba7): the flat
`Termination` with its `nothing`-padded per-source fields becomes
`TerminationRecord{T}` — the final boundary time in the deployment's own
scalar with `nothing` the out-of-band no-boundary arm (D-202's anti-sentinel
argument retiring the old `NaN`), a `TerminationSource` in §13.2's
kind-is-identity convention (`EndTimeReached`, `ModelRequestedStop(face)`,
`ControlRequestedStop(issuer)`, `LoopError(exception)`), and the tail residue
per writer (`ResidueRecord`: the final ring plus suppressed counts). The
§12.1 stop word becomes an issuer word — `Union{Nothing,Symbol,String}`
under a first-wins CAS from `nothing`, `:code` from `stop!(sim)`, the
device's name from `stop!(handle)`, the `:interrupt` arm absent with its
machinery — so every control stop names its initiator, the init bracket's
`should_abort` included. `DeviceJoinTimeout` joins the built kind set
(`who`/`timeout`/`t`/`boundary`, per Appendix C), written to the loop's own
cell at tail step (5) instead of warned inline; the run's-end sweep now
*returns* what it takes, and both advance entries assemble the record once,
in their outermost `finally`, after the sweep — the terminal lifecycle
release-store following the record write on the errored path too, which the
old code ordered the other way. `_advance!` returns bare sources, the
record-per-return-site construction dissolving. Presentation is unchanged in
shape (the sweep still warns per entry, now as the record's renderer), so
the `accounted` either-status-or-sweep invariant held without edits. No
stand-ins enter or retire.

**Increment 18 — the condition algebra, and `init!` as a real service.** §14's
first end-to-end slice: the four node kinds (`Fragment`/`Scoped`/`Combined`
plus §14.6's `Override`), their constructors `fragment`/`at`/`combine`/
`override` — top-level functions, none extending a `Base` generic, `combine`
being D-204's rename away from `Base.merge`'s last-wins contract — resolution
against a `Build`, §14.6's totality check, and §14.4's dynamic-walk
application register, all in `src/conditions.jl`. Composition is inert: `at`
stores a prefix and flattening at resolution is the only place path strings
join, so the deep tree a misaddressed fragment function builds constructs
fine and fails with the build in hand. The collecting pass (§13.1) checks
path, declared field, convertibility, input-face reachability and leaf
uniqueness, and throws once with the full list; the plan it compiles bakes
`xbuf` offsets, cell addresses and the `merge(defaults, overlay)` store
values — the *other* merge, §14.1's fork, and the reason D-204 wanted the
combinator's name back.

One design point inside the algebra is worth naming because it is not the
literal reading of §14.3's `(path, store, field)` duplicate rule: **an input
entry's leaf is the root input it resolves to**, not the face it was written
with. The export chain is therefore walked during flattening, at the same
moment the paths join, and the resolved face becomes the entry's key. That is
what makes §14.6's central use case work — a full-coverage baseline authored
at the root, layered under a patch a component's own `condition` method ships
against its own local face — while keeping two spellings of one root input a
`combine` collision, which they plainly are.

Three integration facts are worth naming. **`init!` now resets the declared
state before applying** — D-063's "fresh run from the `init_*` defaults, with
these overrides", which the previous `init!` never did: the establishment of
`xbuf` and the `s`/`m` stores was factored out of `compile` as
`establish_defaults!`, and both sites call it, so a warm restart no longer
overlays the last trajectory's endpoint. Two existing assertions moved with
the reset and are the fix made visible: a re-`init!`ed `Trigger` fires from
`count = 0` again (1, not 2), and the chatterer's warm restart re-exhausts its
budget from `flips = 0` (8, not 16).

**D-205 re-adjudicated the cells the same day**, in the register D-203 took for
increment 17. Review of the first cut found a second staleness class: a warm
boundary zero left an offset component's cell holding the previous
trajectory's value, where §10.5 and `build.jl`'s own comment promised the
build probe's. The first fix was symmetric with the state reset — an
`establish_cells!` re-seeding the probe products, `init!` and `compile`
sharing it. The design then went the other way and removed the class outright
(D-205): **at boundary zero every discrete output stage runs, due
or not**, publishing from the authored `s` and the `t₀` table, while the `g`
updates keep the `Φ` gate. So `establish_cells!` was retired the day it was
written, and the prototype is the better for it — the virgin-table
precondition is gone, every published cell at `t₀` is boundary-zero-derived at
the deployment activation, and §14.6's barrier now covers the whole table
rather than the root inputs alone. The declared world has **three**
establishment registers, not four: `init!` resets the state homes, boundary
zero derives the cells, totality supplies the root inputs. Root inputs are
pointedly not re-seeded in any of them — re-seeding would put `probe_value`
back on the services path,
which is precisely what §14.6 removed.

The mechanism is one dispatch. `Establish` is a singleton marker the boundary
walk takes in place of a tick index; `run_at!(::Gated, …, ::Establish)` admits
the entry outright where `run_at!(::Gated, …, ::Int)` tests `(idx − Φ) % D`.
The walk's index parameter lost its `::Int` annotation from `Chunk` down so
the one compiled tuple serves both, and since every caller passes a concrete
`Int` or `ESTABLISH` each specializes exactly as before — the §7.5 gates
confirm the frame loop's path is unchanged. `boundary_zero!` in `sim.jl` is
the ordinary macro-sequence with `event_phase!(sim, ESTABLISH)` and
`bodies.ticks(0)`. Frozen components need no exemption logic: at a non-nominal
activation they have no compiled entries at all (§9.4's executable set), so
there is nothing for the wide gate to admit and their pinned cells stand.
**The `Simulation` gained a `build` field** — §14.3 says resolution takes the
root node plus a `Build`, and the schema (declared fields, leaf types) lives
nowhere else. **Root-input totality is structural**, so the probe-value barrier
bit immediately: every test whose model has root inputs now authors them,
which is most of the churn in this increment's diff. The `condition(comp; kw)`
idiom ships in `src/library.jl` as user material — `Plant`, `DiscreteIntegrator`,
and `SampledLoop`/`Vehicle` composing them by pull — and no framework code
knows the name exists. The `set_slot!` stand-in retires here; the two tests
that genuinely needed a *mid-trajectory* direct write (the drain's
counterfactual in test_dataplane and test_roster) keep it as a test-local
`poke!` in `test/utils.jl`, which is test equipment and not framework API.

**Increment 19 — D-206 to D-208, in two commits.** The first is vocabulary
alone: "slot" / "root slot" become **root input** across identifiers, error
messages, comments and both markdown files, with the four API spellings
D-206 fixes — `fragment(; x, s, m, inputs)`, `get_input(face)`,
`UninitializedInputs`, and `RootInputTypeConflict`, which this prototype never
had. Six English other-sense uses survive on purpose (a "single-slot
resource", a "released slot" of a vector, a "scheduling slot"), as does the
retired `set_slot!` stand-in's own name in the retirement history below: the
rename runs forward, never over the log.

The second commit is the semantics. **The reach rule collapses to one level**
(D-207): a wiring endpoint names an immediate child — plus the key segment
where the child is a container element — and one of its faces, and anything
deeper is a build error naming the entry and the offending path.
`resolve_terminal` loses its walk loop and `generically_held` goes with it:
under one-level routing an endpoint stops before any field it could traverse
past, so the generic-holding diagnostic has nothing left to police in the
structural register (§13.3), and the concrete/generic distinction leaves it
entirely — the same declarations now build under either holder, which is
D-207's point. One library spelling and `test_structure.jl`'s reach-rule and
double-feed fixtures were deep and are re-spelled per level; the one that
carries a lesson is `Vehicle`, which used to source `"loop/plant/power"` and
now re-exports `power` from `SampledLoop`'s own boundary, the re-export entry
per level being exactly the ceremony D-207 costs.

**The face graph becomes total, and the `Flat` retains both sides** (§9.2):
`faces` is now `out_faces`, and `in_faces` beside it records every input face
at every level with the producer it routes to. The derivation waits for the
obligation pass, which is the moment it becomes well defined: a face and the
leaf entries behind it are claimed together by the one route above them, so
the consumers of a face share a producer. A primitive's own entries complete
the record, so the graph is uniform in the level's class.

**Any component may be the root** (D-208). A primitive root flattens to the
single leaf at the root path `""`, its `input_types` keys pushed as root
inputs and claimed by the same `("", face)` pseudo-producer an assembly root's
faces get — after which every downstream site was already right: the stage-2
schedule skips the empty producer path, the probe reads the synthesized value
from the layout, `input_addr` finds the root input's own cell, and totality
covers it like any other. Abstract-at-root needed no clause either, the type
derivation being the one it always was.

**Per-level `at` addressing falls out of the retained graph** (D-207's second
ruling): `_root_input` looks the authored face up at `(path, face)` in
`in_faces` rather than in the addressed component's flat entry list, so a
prefix stopping at an *assembly* resolves through that level's face to the
root input behind it, and the component-fed face keeps its own refusal at
every level. State payloads at an assembly prefix stay refused — assemblies
own no state, and only the `inputs` payload gained a level to resolve at.

**The third commit closes the two silences the conformance build exposed**
(D-210). §8.6's face-name uniqueness invariant now reads the root by the
root's *class*: a primitive root's face set is its `input_types` and
`output_types` keys together, checked where `flatten` meets the root, because
that is where the two declarations first share the periphery's address space.
The repro is worth keeping in mind — `input_types = (u = T,)` beside
`output_types = (u = T, v = T)` had `cell_layout` place the root input's cell
over the output port's, so the authored `5.0` read back as the stage's `10.0`,
boundary-zero-ordered and silent. Non-root primitives are deliberately
untouched: below the root an input face aliases its producer's cell and places
nothing, so there is no address to collide over. The second ruling is one level
down: an `input_connections` entry routes to at least one internal endpoint, so
the empty tuple is a declaration error wherever it appears. It used to
contribute no consumer, leave no row in the face graph, and let a condition
addressing that face read "declares no input face" — byte-identical to a bare
typo's diagnostic. The check sits in the walk over `input_connections`, where
the fan-out is already in hand, and it makes the `isempty(consumers)` guard in
the `in_faces` derivation dead: every route now pushes its row unconditionally.
`_root_input_type`'s "routes to no input" throw stays as an internal invariant
with nothing left that can reach it.

**Increment 20, part 1 — name-transparent containers** (D-211). The
declaration layer gains one optional entry, `transparent_container(c)`,
defaulting to `nothing`: a component may name at most one of its container
fields, and that field's elements are then contributed under their bare keys —
`"key"`, `"1"` — in place of `"field/key"` and `"field/1"`. Naming is the only
thing it changes, so the whole of it lives in `children`: the same fields in
the same order, one branch on `name === tf` deciding the segment. What pays for
it is one declaration check and a collision family in three arms. The
declaration must name a container field of the type, validated *after* the
field walk so that a mixed container still reports as one. Then no two children
may share a name — a check written generally over the collected list rather
than over the transparent case, because bare keys only make the collision
*reachable*, they do not define it, and the general form catches the
pathological cases (a bare key against a sibling field, against another
container's composite name) in one place. The family's other two arms sit where
the bare key is formed, because no child name collides in either and the
general check can see neither: an element keyed with its own field's name,
which `sample_times`' field-name sugar already spells, and one keyed with a
sibling *container* field's name **where that field contributes children**,
which collides with no child yet shadows the `"field/key"` grammar reaching
them — the review's find, spec'd as the §8.5 clause D-212 records: the elements
stay in the flat list and become unreachable in the structural register, the
exact-match-first lookup in `_one_level` answering with the bare child and the
one-level rejection then blaming it rather than the shadowed container. The
qualifier is load-bearing and was itself adjudicated: an *empty* container
reaches nothing, so there is no grammar to shadow, and reserving its name anyway
would refuse over empty *inert* data — an empty `Tuple` is both at once, the
value cannot tell them apart, and the walk here already treats them as one case.
That leaves legality per-instantiation, which is the framework's norm rather
than an exception to it: every wire is validated against the instance too. All
three arms name both parties by provenance. `resolve_terminal`'s two-segment lookahead is untouched: it still serves
the undeclared containers, whose key segment D-211 leaves exactly as it was.

The one place the rename reached past naming was the rate declaration's
field-name sugar, which `_rate_entry` had implemented as a *path-segment*
prefix match (`startswith(seg, "children/")`). §8.5 keys that sugar on the
field, and the two spellings coincide only while the field segment rides the
child name — under transparency the prefix match would have silently stopped
applying to `(children = Relative(2),)`. `_children` now returns the
contributing field beside each child and `_rate_entry` matches on it, which is
the spec's rule and a strict generalization of the old one: for an undeclared
container, field-keying and prefix-keying are the same lookup. `_check_sample_times`
admits the field name as a key on the same basis. What that leaves is the one
ambiguity §8.5 names — a transparent element keyed with its own field's name —
and it is refused where the bare key is formed.

`Group` is the library case and the reason for the decision: it declares
`children` transparent, drops the 4-arg positional convenience constructor for
`Group(children; wires, inputs, outputs, rates)` — a bare `Pair` normalized to
a one-entry tuple, since a single wire should not have to be written with a
trailing comma — and its declarations then read exactly like a named
assembly's. The test sweep is the visible half: 335 `children/` occurrences
across `src/` (19) and `test/` (316) became bare keys, including the 16
`var"children/ctl"` rate spellings that a slash-bearing key had forced, the nested read paths
(`"children/f/children/a"` → `"f/a"`) and the error-message assertions, two of
which were re-pinned against a quoted form rather than left matching a
one-letter substring.

**Increment 20, part 2 — §13.3's primitives and §8.8's passthrough pair.** The
declaration surface gains its computed half. `resolve` was a reserved name the
condition algebra had taken, so the rename came first: `resolve(node, b::Build)`
is `resolve_condition`, internal and behaviour-preserving, and §13.3's path
primitive takes the name it is spec'd under. `resolve(asm, path)` and the public
`resolve_terminal(asm, path)` are then *entry-less* forms of the walk wiring
resolution already ran: the one-level check was factored out of the internal
`resolve_terminal` into `_one_level`, parameterized by how many segments follow
the child — one face name for a terminal path, none for a bare child path — so
there is a single implementation of the rule and the container lookahead serves
both. Two accessors were already there and needed only their §13.3 form:
`input_faces`/`output_faces` returned a `Tuple` for a primitive (`String.` over
a NamedTuple's keys broadcasts to a tuple) and a `Vector` for an assembly; both
now return `Vector{String}`, which is what the section specifies and what the
helpers' `setdiff` wants.

The helpers themselves are the thin composition §8.8 promises: `input_faces` of
the resolved child, a shared filter, and a `Tuple` comprehension. What is worth
noticing is how much of their error surface is *not* theirs. Exclusivity and the
unknown-face list are checked in the helper because nothing downstream could
name the offender; everything else is left where it already lives — a
`prefix = ""` collision meets §8.6's face-name uniqueness, a face both wired and
passed through meets §6.1's one-producer rule with both claimants named, and an
`except`ed face the assembly then fails to wire is an ordinary unconnected
input. The tests pin exactly that division. `declaration_error`'s two shapes in
the spec sketch are §13.2's structured carrier, which the README records as
absent: these throw a plain `BuildError`, the house framing. One default was
amended after the increment closed (the §8.8 doc ruling in the same commit as
this paragraph's): `prefix` defaults to `child_path` with its slash folded into
`sep`, so the blessed container-element path (`"units/1"`) labels its faces
legally (`"units.1.…"`) by default instead of minting a slash-bearing face name
for §8.6's check to refuse; an explicit `prefix` is used verbatim.

**Increment 21, part 1 — the compile product gets a name.** `compile` returned
a NamedTuple of fresh buffers with the phase bodies closed over them; that
product is the spec's *executor* (§9.7) — the compiled execution form of the
schedule over its own buffers, at one scalar — and it is now `Executor{T}`,
declared beside `compile` in `build.jl`. The `Simulation` traded ten fields for
one: `store`, `xbuf`, `ẋbuf`, `clock`, `bodies`, `events`, `sstores` and
`mstores` all live on the executor it owns, and `layout` and `flat` were
duplicates of `exec.act.layout` and `build.flat`. No `getproperty` forwarding
was added, deliberately: every call site now spells `sim.exec.xbuf`, and that
is what keeps ownership visible once the services arrive and instantiate
executors of their own (§9.2: "every buffer set has exactly one owner"). The
three evaluation entry points moved with the buffers — `evaluate!`, `_round!`
in its three arities and the condition register's `apply!` all take an
`Executor` now, with one-line `Simulation` methods delegating — so a service
can evaluate a scratch world with no `Simulation` to hang it on.

Nothing else changed, and nothing was meant to: the extraction is
behaviour-preserving by construction, and §7.5's gates are its acceptance test.
`step!`, `evaluate!`, the four phase bodies in both arities, `boundary!`,
`offtick_boundary!`, `frame!` and `drain!` all still measure zero — the extra
`getfield` through `sim.exec` folds away exactly as the seam's other
indirections do. One trap surfaced on the way, worth recording because it is
invisible in the diff: naming the local `evset` was forced, because a local
named `events` inside `compile` shadows the `events(c)` declaration accessor
the same function calls, and the model then fails to build at all.

**Increment 21, part 2 — the read side, and `capture`.** §14.4's read-selector
family was already half-built: increment 15 gave the output binding
`get_output`/`get_input`/`get_face` as deferred-read values, because a binding
needed them. The two store members close the family, and closing it moved it:
the five selectors, `reads(…)` as the labeled set in one type, and the compiled
`Reader` now live in `src/readers.jl`, which is included *above* the data
plane — the binding register is a client of the family, not its owner, and its
`_resolve_read` methods dispatch on selector types that therefore have to exist
by then. What readers.jl needs from the condition algebra in the other
direction — the `x`-offset walk, the "is this path a level of the build at all"
predicate — it calls at resolution time, long after conditions.jl has been read,
which is the ordinary late-binding this file layout already relies on
(`init!` calls `resolve_condition` the same way).

The reader is `apply!`'s **gather twin** in the literal sense: one entry per
selector, its leaf type baked into the entry's own type, and a `map` over the
entry tuple that the compiler unrolls into the four reads a model can offer —
an `xbuf` offset, the `ẋbuf` offset beside it (`ẋ` has `x`'s shape at the
activation scalar, so the derivative of a field sits at the field's own
offset), a discrete store's index and field, and a cell address. One trap is
worth recording. The `s` stores are held in a `Vector{Any}`, so the read needs
a type assertion to stay inferable — and the assertion has to go on the
*reference*: `ex.sstores[ci][]::S` infers correctly and still allocates 16
bytes per read, because the `[]` on an `Any` is a dynamic call that boxes its
return before the assertion sees it. `(ex.sstores[ci]::Base.RefValue{S})[]`
measures zero.

Closing the family also made §14.4's **source rule** enforceable for the first
time: the store selectors resolve only against live stores, and a binding reads
a published snapshot, which carries none — so `get_state` in a binding's
`reads` is now `ReadBindingUnresolved` at attach, with the honest remedy
named, where before it could not be spelled at all. `capture` sits at the end
of `conditions.jl` rather than in readers.jl: it reads stores directly rather
than through selectors, it builds a condition tree out of `fragment`/`at`/
`combine`, and it is the gather twin of the `apply!` it now follows in the
file. Its round trip has one caveat that belongs to boundary zero rather than
to capture, and the test fixture is built around it: re-applying a captured
condition runs the §14.5 sequence, so `project` and any guard already holding
fire again — and so does a *due* `g`, the outgoing transition of §14.5's
table. A bit-for-bit round trip is a claim about the establishment, not about
a `g` that would run a second time, so the fixture holds the integrator's
input at zero and the claim stands undiluted.

Two names had to move. `Reader` was taken by a coverage binding in
`test_roster.jl` (`struct Reader <: AbstractBinding`), and since every file
here is included into `Main`, the test's definition silently replaced the
framework type before `test_readers.jl` ran; the binding is now `Unwritten`,
which is what it was for. And the collected refusals carry no kind of
their own: Appendix C names none for an unresolved read as such — its spellings
are `TrimProblemInvalid` where trim owns the setup and `TapResolution` where
linearization does — and a read set is only ever compiled inside a client, so
the standalone entry point is internal (`_compile_reads`) and its collecting
half is factored apart from the throw exactly so trim can fold the same list
into its kind.

**Increment 21, part 3 — the specialized register.** §14.4's two application
registers finally both exist, and the second one is what an iterating service
spends. The dynamic walk bakes *values*, so a plan is good for exactly one
tree; the specialized register bakes *lenses*, so a plan compiled from a tree's
shape applies to every later tree of that shape. Making that real took one new
fact in the flattening pass — each entry's **tree position**, the
`getfield`/`getindex` step tuple from the root node down to the authored value,
which for an `override` leaf is the winning layer's — and then `Getter{P}`, the
position lifted to a type parameter with a generated navigation, exactly as the
glossary describes it.

The interesting part of the increment was not the lens but the **factoring**.
§14.3's checks had one implementation and had to keep having one, so
`resolve_condition`'s body split in two: `_resolve_entries` runs the collecting
pass and hands back the survivors as `Resolved` values — the authored entry, the
destination the activation's layout supplies, the destination leaf type, and the
value through that type's `convert` — and the two registers then bake what they
need from the same list. The dynamic one takes the value; the specialized one
takes the position and the leaf type. Nothing else differs, which is exactly the
claim §14.4 makes about the pair.

Both of §14.3's converter arms turned out to be the same call. The destination
leaf type at the activation *is* the converter, so `convert(L, lens(tree))`
covers a `Dual` into a `Dual` leaf (partials untouched) and a plain `Float64`
into a `Dual` leaf (the zero-partial embedding) with no branch anywhere — and
the one case no converter covers, a decision variable authored into a leaf that
is pinned `Float64` at the seeded activation, was already refused by the
existing `_unconvertible` path, which now carries a clause saying *why*: a
decision variable descends into neither a frozen discrete `s` nor a pinned leaf.
That clause fires on a type test rather than a message guess — the value's leaf
eltypes contain the activation scalar and the destination's do not.

The shape check folds in two halves, per §9.5's mechanism. The tree type is a
plan type parameter, so a tree of another shape simply does not match
`apply!`'s method and lands on a fallback that raises `ConditionShapeDrift`
naming both types; a prefix, being a runtime `String` field, cannot ride in the
type, so the plan records the prefixes it was compiled from and `apply!` sweeps
them with `===` before any write. Two things make that sweep honest here.
`===` on `String` is *content* equality, not pointer identity — `egal` is
specialized for it, and `f(i) = string("ge", "ar/", i)` gives `f(1) === f(1)`
at two different addresses — so what the sweep tests is that the prefix still
spells the same path, and a prefix a service *computes* afresh per iteration
(`at("gear/$i", …)`) passes it exactly when it should. Only the all-literal
case additionally folds the compare away at compile time; elsewhere the cost is
a length check and a short `memcmp`. And the sweep running *first* is what
makes a refused application leave the executor bit-for-bit as it found it,
which the drift tests assert on both halves.

The store write is where the composite matters. A service compiles its plan
from the whole tree it builds, `override(baseline, condition(d))`, never from
the patch alone — the `s`/`m` write is one whole `merge(defaults, overlay)`
value with the base baked, so a plan over the patch alone would silently reset
every field the baseline authored and the patch did not. No coverage component
declares a two-field `s`, which is what a single-field store hides behind
`override`'s own last-wins, so the test fixture `Ledger` exists to make the
property visible. The allocation trap Stage 1 recorded reappears verbatim on
the write side and takes the same remedy: the store type is baked into the
plan entry's own type and the assertion goes on the *reference*,
`(ex.sstores[ci]::Base.RefValue{S})[] = v`.

**Increment 21, part 4 — the trim service.** With both application registers,
the compiled reader and `Executor{T}` in place, `trim!` turned out to be mostly
the loop *between* them plus the vectorization at its edges. What it genuinely
adds is the discipline §14.8 is about: the simulation's stores have exactly one
writer, the commit through boundary zero. Every evaluation runs on executors
the invocation instantiates from the build's cached layouts and drops on
return, so "no convergence, no commit, the sim bit-for-bit untouched" is
structural rather than promised — the tests assert it on a never-initialized
simulation (still `built`, `run!` still refusing) and on an initialized one
(every buffer equal to its pre-call copy), and neither assertion needed any
code to defend it.

D-213 is the part with real content, and it is invisible unless the fixture is
chosen for it. The seeded activation freezes the discrete tier (§9.4), so
nothing there can derive a discrete output cell from an authored `s`; without
the ruling those cells hold the build probe's synthesized values for the whole
solve. The fixture is a pendulum whose torque comes from a `DiscreteIntegrator`
whose `acc` the *baseline* authors: the answer is `asin(acc/(g/l))` exactly when
the scratch world's frozen cell holds the authored output, and `0` when it holds
the probe's. Both halves fell out of machinery that already existed — the
dynamic walk for the nominal set, `_round!(ex, ESTABLISH)` for the establishment
round, and a value copy of the frozen components' output cells into the seeded
set — which is the ruling's own argument for itself. Its other half, a decision
variable authored *into* that frozen `s`, is refused by the seeded compile's own
`ConditionResolution` and deliberately not folded into `TrimProblemInvalid`: the
condition is what does not resolve, and the collected read-set violations *are*
folded in, so the file now shows both dispositions side by side.

The backend seam behaved exactly as §14.8 predicts. `LevenbergMarquardt`'s
`solve` is forty-odd lines — a damping loop, a small dense solve, the
per-residual box test as its own stopping rule — and it converges the one-step
linear problem in 3 iterations / 4 evaluations and the nonlinear one in 4 / 5,
which is §14.7's "quadratic, ~5–15 evaluations" on a two-equation model. The
one decision worth recording is the *acceptance* test, which measures
`‖r ./ tol‖` rather than raw `norm(r)`. §14.8 puts the normalization on the
derivative-free objective, but its reason is not the backend's: a raw norm sums
forces against moments, and on a problem whose residuals differ by orders of
magnitude in physical scale *and* whose box is active it can reject every
projected step and stall short of a point the box test would accept. In the
least-squares register the section's own sentence settles it — "the tolerances
*are* the stopping criterion … LM's damping loop testing exactly what the
service will re-test" — so descent test and stopping rule are read in the same
units, and `_scaled_norm` is a fold rather than `norm(r ./ tol)` so it costs no
temporary. The division is what makes a non-positive tolerance a *setup* error
rather than a numerical curiosity: `tol = 0` sends the acceptance test to
`Inf`/`NaN`, every trial step is rejected, and the solve returns `:stalled` at
the guess — so the collecting pass now requires each tolerance finite and
`> 0`, beside the `lower ≤ upper` clause and for the same reason, that a bound
no arithmetic can honor is the problem being malformed. Neither iteration count
moved: still 3/4 and 4/5.

One trap cost an afternoon and is worth recording because it is invisible in
the diff. The store bundle keys its per-eltype buffers by `Symbol(L)`, and
`show` for a type abbreviates the module prefix when the type is *visible from
the printing module* — which is not the same module inside a `@generated` body
as it is in `compile`. Every scalar the prototype had used until now was tagged
by something in `Base` (`Float64`, `Dual{Nothing,…}`), where the two spellings
coincide. `TrimTag` is the service's own, declared where the prototype declares
everything, and the bundle was built with a field named
`ForwardDiff.Dual{TrimTag, Float64, 1}` while the generated gather looked up
`ForwardDiff.Dual{Main.TrimTag, Float64, 1}` — a `FieldError` from inside
`compile`, on the first `Dual` activation any service ever asks for. The fix is
one shared `_cell_key(L)` printing with no module context at all, so both sites
agree on the fully-qualified spelling; the three tests that named the old
spelling were updated with it. Any user tag would have hit this, so it is a
defect the increment found rather than one it caused.

A guess outside its own box is **projected at the pack site**, `clamp.(d0,
lower, upper)`, not handed to the backend as given. Step projection alone does
not cover it: `solve` has two returns that never take a step — the
already-within-tolerance return at the top of the first iteration, and a stall
at iteration one — and through either of them an unprojected guess comes back
as the solution, is committed as converged, and is reported with `saturated`
empty because a point outside both bounds sits at neither. Projecting the one
point the backend starts from makes every point it sees and returns lie in the
box, which is what §14.8's bound treatment is *for*, and it costs one line. The
degenerate box `lower == upper` falls out of the same clamp: the decision is
pinned, and `_saturated` testing `lower` first reports it as `:lower`.

`nevals`/`niters` in the report are the *backend's*, verbatim — the service's
own verdict evaluation at the returned point is not counted, because §14.8 calls
those counts "the backend's returned status together with its
iteration/evaluation counts" and the verdict evaluation is noise against the
solve by the same paragraph's argument.

One thing the no-throw doctrine deliberately does not cover: a residual lambda
that throws at the *committed* point — after `init!` has run — escapes `trim!`
as that exception with the simulation committed, because the doctrine is about
non-convergence being an outcome rather than about catching broken user
machinery, and the service does not catch it.

**Increment 22 — the diagnostic carrier, in four stages.** §13.2's kind values
replace the message strings the prototype had raised since increment 2. Stage 1
adds `src/diagnostics.jl`: `Diagnostic` as the root of the closed set, each kind
an immutable `@kwdef` struct whose fields are Appendix C's payload column, three
methods apiece — `severity(d)`, a property of the *kind* under D-214 (`:error`
by default, the warning kinds overriding, never stored per occurrence),
`path(d)`, the renderer's sort key, and `message(d)`, which has no default so a
kind added without one fails loudly at its first rendering rather than printing
a stub — and the one `BuildError` carrier over a `Vector{Diagnostic}`, whose
`showerror` renders compiler-style: a lone diagnostic on one line, a collection
under a count line, grouped by kind in first-appearance order and stably sorted
by `path` within a group. Where Appendix C's payload column reads "reason (a /
b / c)" the kind carries a `Symbol` `reason` — one struct per kind, never one
per reason, which is what lets `d.reason === :inverted_box` be a test's whole
assertion. `dataplane.jl`'s eight runtime kinds are re-parented under
`Diagnostic` so `severity` covers them too. A `BuildError(::String)` constructor
wrapping a `LegacyMessage` kind rode through stages 2–3 so the suite never went
red, with its own README stand-in row; stage 4 deleted both.

Stages 2–4 convert the sites: the build side (`assembly.jl`, `build.jl`), the
services and periphery (`sim.jl`, `roster.jl`, `bindings.jl`, `devices.jl`,
`declare.jl`), and the three collecting registers — `conditions.jl`'s
`_report_violations`, `readers.jl`'s `_compile_reads`, `trim.jl`'s
`_report_trim!` — whose `Vector{String}` becomes the `Vector{Diagnostic}` the
barrier hands to one `BuildError`. Trim's setup throw now carries two kinds at
once, its own `TrimProblemInvalid`s beside the read set's `TapResolution`s,
which is what a carrier over a heterogeneous vector buys and what the string
splice was imitating. The `logged(d)` rendering — the kind name, then the
message, the carrier's line without the carrier — is the `logged` policy's
form, used by the two trim warnings and by `attach!`'s `EmptyGreedyClaim`.

**What stage 2 could not make collect.** §13.1 asks for one barrier per
stratum, and Stratum A's passes are interleaved with user code:
`input_connections`, `output_connections`, `sample_times` and
`transparent_container` all run mid-walk. The *checking* parts now collect; the
*resolution* parts — `classify`, `resolve_source`/`resolve_dest`/
`resolve_terminal`, the one-level and direction rules — still abort at the
first violation and drop what the walk had gathered, because a resolver that
kept going would need a sentinel return the whole wiring API does not have.
That is its own increment. The split is worth naming because it is invisible in
the diff: a model with two bad wire *ends* reports both, and one with a bad
wire and an unresolvable path reports the path alone.

**`InternalInvariant` is not a diagnostic** (D-215). The activation-identity
refusals and the handful of `error(...)` assertions raise their own exception
type, carrying a message and no payload, because they name no failure a user
can fix: a `Reader` compiled at `Float64` meeting a `Dual` executor is a
framework invariant broken, not a model that is wrong. Keeping them outside the
kind set keeps the acceptance-test contract — match on kind and payload —
honest, and the four tests that read their text do so on purpose.

**One rendering testset, marked.** `message(d)` is presentation, and no test
outside `test_diagnostics.jl`'s `@testset "rendering"` reads it; everywhere
else the assertion is the kind plus the payload fields carrying the same fact
(`"5 violations"` becomes `length(e.diagnostics) == 5`,
`startswith(e.msg, "TrimProblemInvalid:")` becomes `d isa TrimProblemInvalid`,
a provenance chain becomes the `provenance` field compared as a value). The
directive prose the old assertions matched — "use `override(base, patch)`",
"stopped → init! → run!" — is rendering and went with it; the fact each was
pinning survives as payload. The exception exists because the didactic register
is a claim of its own: state the fix, show the list in hand, lead with the kind
name. That is what the rendering testset pins, and it is the only place in the
suite that may.


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
- **The allocation invariant is width-independent.** Chunking bounds compile
  cost, never the runtime claim: a phase body of 18 one-entry chunks
  (`chunk_size = 1`, a chunk count no other fixture approaches) walks both
  variants at zero allocation, so the §7.5 canary covers the outer walk over
  the chunk tuple itself, not just the entry walks within one chunk.
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
- **A wiring endpoint reaches exactly one level.** The same `"inner/sum/a"`
  against the same `SampledLoop` value is a build error whether the field is
  declared `::SampledLoop` or `::L` — the declaring type's knowledge stopped
  being the question with D-207 — and `"inner/ref"`, the sub-assembly's own
  face, builds under either holder to the same trajectory. A container's key
  segment rides along inside the one level (`"units/1/e"`), which is what lets
  a container's elements be addressed by their parent's own declarations.
- **Being fed is a whole-tree obligation, not a per-declaration one.** An unfed
  child input inside one assembly is merely awaiting a claim from above — a
  sibling wire, or an `input_connections` entry handing it up a level, which
  an ancestor's route then feeds through the face. The error fires at the root,
  for the chain that never terminates; the one legitimate terminus is the root
  component's own input face. The one-producer rule spans levels the same way,
  so an ancestor's route onto an input a sub-assembly already wires is caught
  where the two claims meet, with both entries named.
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
  interpolant cost.** A root-input write between runs is the frame-top epoch
  seam: the re-measured σ₀ holds under the frame's own `u`, so the edge is the
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
  no live root input; the batch lands at the top of the next frame `run!`
  advances — and nowhere earlier. A batch staged while stopped waits through
  boundary zero (`init!` runs on the un-drained root input — the contrast with
  `set_slot!` before `init!`, which makes a holding predicate fire at `t₀`),
  and the frame's outcome is a pure function of the drained batch: the staged
  trajectory is bitwise the directly-written one.
- **Merge is the only coalescing policy, and every check runs at staging.**
  Newest wins per face and untouched faces survive — a sparse `b` batch cannot
  clobber a pending `a` (§15.3's flaps/gear hazard), a re-staged face takes
  the newest level. The checks run on the writer's side: a face with no
  position in the schema is discarded under `OutOfClaimEntry`, an
  unconvertible value under `EntryTypeMismatch`, the rest of the batch stands,
  and the drain is pure application — the shim having already converted to the
  activation's root-input types, so a `Dual` deployment stages plain reals. The
  empty drain is allocation-free: a quiet loop's frame top costs nothing.
- **Nothing reachable from a published snapshot is ever written again.**
  Publication is one buffer copy and one release-store after the boundary
  sequence completes — boundary zero included, off-tick frame tops included —
  and the snapshot is the whole table, root inputs riding along as the source
  cells they are, state stores deliberately absent (§11.2). The run moves on;
  the snapshot holds. A concurrent reader acquire-loading `latest` sees only
  coherent worlds — the in-lockstep pair `g2 = 2·g1` holds bitwise in every
  observed snapshot, and `t` never decreases — and a task staging against
  running frames loses nothing to the CAS merge: the last staged level is what
  stands.
- **One writer per root input at any time, structurally.** Claims are disjoint
  by admission, so cross-writer races on one cannot arise and no drain-order
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
  the root input receives the conditioned level and the model output is
  bitwise the directly-staged one — while an entry declaring no parameter
  passes its value through untouched, which is what carries a throttle's
  `[0, 1]` level or a press counter without ever meeting the axis convention.
  The claim is the table: what the binding may write is exactly what its
  entries name.
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
  path, the unknown root input, the input face offered to `get_face` — so
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
- **The drain is free, populated or empty, at any width.** The batch's
  values-plus-mask pair is one concrete isbits layout per writer (D-202), so
  the scatter is a single specialization compiled at the stopped-sim point:
  a never-before-drained sparsity pattern allocates exactly what the warmed
  one does — nothing — for the harness register and device writers alike,
  and a 34-face surface (past the 32-wide threshold where Base's tuple `map`
  leaves its inlined path, which the generated merge never touches) drains
  as free as a two-face one.
- **Termination is a state with a typed record, never an exception (D-203).**
  Every ended run answers "why did it stop?" and "how did the stop go?"
  through `termination(sim)`: `EndTimeReached` with the final frame top,
  `ModelRequestedStop` carrying the first holding face in declaration order
  and the boundary time of the very snapshot that held it,
  `ControlRequestedStop` carrying its issuer — `:code` for calling code, the
  device's name when a handle or its `should_abort` spoke, first CAS wins —
  and `LoopError` with the cause retained. The record is `nothing` unless
  the lifecycle is terminal, `init!` clears it with the trajectory, and its
  residue holds what landed past the final account: the tests pin the
  `DeviceJoinTimeout` entry (who/timeout/t/boundary) on the loop's record
  for an abandoned join, the pre-spawn `DeviceCrash` on a `should_abort`
  init failure's zero-frame run, and the empty vector on a clean tail.
- **A stop face ends the run at its own boundary — `t*` included.** The
  boundary-detected face ends the run at the sweep that saw it (frame 4 for
  a ramp crossing 0.35 at h = 0.1, never t_end's frame); the localized one
  ends it at the crossing's `t*` within tol of the analytic instant, the
  clock left at `t*`, the frame's remainder never integrated, and the log's
  terminal endpoint *is* the `t*` snapshot; an authored condition already
  terminal ends the run at t₀ with boundary zero final and zero frames
  advanced.
- **`init!` is the one door, and it opens the trajectory wholesale.**
  `run!`/`step!` before it refuse naming it; `stopped → init! → run!` is the
  re-run cycle; and the §12.6 clearing is pinned beside the §11.4 wait in
  one testset — a batch staged *before* `init!` is discarded with the
  trajectory it predates, the same batch staged *after* waits for frame 1's
  drain and never reaches boundary zero. The freeze on the far side of that
  door — `init!` and `run!` both refusing while the lifecycle is `:running` —
  is pinned against a run the test holds open at *both* ends: it spins for
  `:running`, probes, and only then stages the root input that trips the stop
  face, so the run's end is that staged write's consequence rather than a frame
  count the probes have to outrace. The same testset is therefore the
  end-to-end witness for a staged root input causing a model-requested stop:
  `stage!` mid-run on a deviceless simulation reaches the whole root-input
  surface through the harness register (§11.3), the next frame's drain applies
  it, and the boundary that follows publishes the holding face.
- **A stepped frame is a run frame, bitwise.** The two advance entries share
  one frame loop, so the equality is by construction and the tests assert it
  with `===`; the count actually advanced is the return value, so t_end and
  stop-face truncation are detected without inspecting the clock.
- **§13.6 costs nothing because publication is a boundary's last act.** The
  failing frame published nothing, so the "promoted" final snapshot is
  simply the newest published one: the record and `latest` agree on it, the
  log ends at it, the device's bracket closed through the ordinary tail, the
  mid-boundary stores stay readable, and `errored` refuses `run!`, `step!`
  and `init!` alike, each by name.
- **Composition is inert, and that is testable.** A `combine`/`at` tree over
  a path that resolves against nothing *constructs*: the nodes are inspected
  for their stored prefixes, unconcatenated, and the fragments are `isbits`;
  only `resolve` throws. The same test pins the stack-only property the trim
  loop will need before trim exists.
- **`combine` collides, `override` layers, and both say who wrote what.** A
  duplicate leaf reports both provenance chains — the `at` prefixes and the
  payload position, verbatim — and names `override(base, patch)` as the
  remedy; `override` gives the shared leaf to the patch, passes untouched
  leaves through, composes variadically, and still refuses a collision
  *within* a layer. The dual provenance is asserted through a violation on
  the overridden leaf, which is where the chain surfaces.
- **Layering is keyed on the root input, not on the spelling.** A root-level
  baseline under a component-local patch on the same root input resolves to the
  patch's value, with the baseline's others untouched — §14.6's central
  use case, which a literal `(path, store, field)` key would reject with the
  very advice it was following. The same two spellings under `combine` still
  collide, so the directive that branch prints is advice that works.
- **Boundary zero publishes every output stage, due or not (D-205).** Two
  components at `Relative(2, 1)` — neither due at `t₀` — both publish there:
  the ZOH from the condition's root-input value, the integrator from its
  authored `s`, and the probe's synthesized `0.0` reaches no published cell,
  live table or `t₀` snapshot. Dueness still governs the updates: the authored
  `s` survives boundary zero untouched and the first sample the integrator
  *consumes* is its own `Φ·Δt_base` tick's, pinned analytically as one period
  times one sample. Cold and warm `init!` agree exactly, the authored
  condition determining the `t₀` table outright — which is what retired the
  virgin-table precondition rather than guarding it.
  `test_multirate.jl`'s offset testset carries the same rule at its own
  fixture: the assertions were already the evaluated values, so only the claim
  behind them was re-pointed.
- **The service entry point speaks the algebra's diagnostics.** A bare
  NamedTuple handed to `init!` where a condition belongs gets
  `ConditionNodeMisuse` and the wrap-it directive, never a `MethodError`.
- **The misuse kind is a composition-time fact.** Blending a node with a bare
  NamedTuple throws `ConditionNodeMisuse` in both argument orders, at every
  arity, in `at` and `override` too, and with no build anywhere in sight —
  which is exactly why it is its own kind rather than a `ConditionResolution`
  sub-kind.
- **Resolution collects.** Five different violations in one condition —
  unknown path, undeclared field, unconvertible value, an input face wired
  internally, one root input written twice — come back as one throw naming
  all
  five. The undeclared-field message discriminates outputs and workspace by
  name, because "conditions never specify outputs" is the rule the author is
  most likely to test.
- **An input face is resolved through the export chain, not guessed.** A face
  authored at the component that declares it lands on the root input its
  obligation chain ends at; an internally wired input is refused, because the
  first sweep would overwrite it and unexported stays unpokeable.
- **Totality is pre-write, and "pre-write" is asserted by inspection.** A
  condition short one root input is rejected with every uncovered face named in
  declaration order, and the simulation is then checked to be bit-for-bit what
  it was: same `x`, same stores, same root-input cells, same lifecycle, the
  same
  published snapshot *object*. `init!(sim)` with no condition stays legal
  exactly where the build has no root inputs.
- **The overlay base is the declared defaults, always.** A sparse condition
  leaves unnamed fields at their `init_*` values, and a second `init!` after a
  full run restores them bitwise — `SVector(0.0, 0.0)` and `(state = :armed,
  count = 0)` again, not the trajectory's endpoint. This is the D-063
  behaviour the previous `init!` silently lacked.
- **An authored state fires its guard at t₀.** Boundary zero establishes every
  prior as not-holding, so a condition landing a predicate in holding
  territory fires visibly at `t₀` — authored, not staged, and with zero frames
  advanced.
- **The fragment-function idiom composes by pull across two levels.** The
  vehicle scopes the loop's fragment under `at("loop", …)`, the loop scopes
  the plant's under `at("plant", …)`, and nothing anywhere writes
  `"loop/plant"`: the deep path exists only in the flattened entry list.
- **One-level routing is class-blind and holder-blind.** A route through the
  sub-assembly's own face builds identically whether the field is declared
  concretely or generically, and both holders are rejected identically one
  segment further in — the endpoint stops before any field the old rule could
  have policed. The message names the child it reaches past, so the repair
  (declare the face on that child's boundary) reads off the error.
- **The cross-level one-producer rule survives the reach rule that used to
  demonstrate it.** An ancestor can no longer route deep onto an input a
  sub-assembly already wires; it feeds that input *through the face*, and the
  two claims still meet with both entries named — `child_connections` at the
  sub-assembly and at the root component.
- **A primitive root is the whole model.** `build(Plant())` flattens to one
  leaf at the root path, its `input_types` keys the root inputs, and the
  deployed bare leaf integrates `ẋ = A x + B u` against the matrix-exponential
  reference under a root-input-held `u` — no wrapper assembly, no vocabulary
  the assembly root does not already use. Totality reaches it identically, and
  the leaf's own `condition` fragment composes at the root with no `at` prefix.
- **An `at` prefix stopping at an assembly resolves to the same plan the root
  spelling produces.** `at("loop", fragment(inputs = (ref = …)))`
  and `fragment(inputs = (in = …))` compile entry for entry, because the face
  graph the `Build` retains has a name to follow at every level. The two
  refusals that flank it are unchanged in kind: a component-fed face reaches
  no root input, and a state payload at an assembly prefix has nothing to
  write.
- **Face uniqueness reads the root's class, not its family (D-210).** A
  primitive root declaring one key in both contracts is rejected as the same
  build error a duplicate assembly face name is — the shape that used to build
  and place two cells at one address, the root input's over the port's. The
  identical leaf one level down still builds, which is the rule's other half:
  below the root an input face places nothing, so there is nothing to forbid.
- **A face feeding nothing declares nothing (D-210).** An `input_connections`
  entry routing to the empty tuple is refused at the level that declares it,
  root or not, with the offending entry named. The condition misdiagnosis it
  used to produce is *unreachable* rather than reworded: the test authors the
  dead face, addresses it with a condition, and gets the build's refusal — the
  resolver never sees it.
- **Transparency is a naming change and nothing else (D-211).** The same two
  children, in the same declaration order, wired and read through bare keys on
  a transparent container and through `"units/1"` on an undeclared one — the
  flat list, the wiring endpoints, the read path and the exported face all land
  identically. The default is `nothing`, and `TupleRoster` sitting beside
  `TransparentRoster` in the same file is the control.
- **The collision family has three arms, and each names both parties.** The
  duplicate-*name* check is the general one, written over the collected list
  rather than over the transparent case: a bare key against a sibling component
  field, and a bare key against another container's composite name, are one
  error there. The other two arms exist because no child name collides at all,
  so that check can see neither — an element keyed with its own field's name,
  which `sample_times`' field-name sugar already spells, and one keyed with a
  sibling *container* field's name **where that field contributes children**,
  which shadows the `"field/key"` grammar reaching them. Both are refused where
  the bare key is formed, in the same naming-both-parties form the general arm
  uses. The qualifier is pinned by the same parametric type one instantiation
  apart: populated, the bare key is refused; empty, it stands and both the
  wiring register and the reads resolve it to the bare child — an empty
  container reaches nothing, so nothing is shadowed, and the value cannot tell
  it from empty inert data anyway. The declaration itself is checked too: a
  component field and an absent name are refused alike.
- **Bare rate keys drive the grid the composite ones did.** The same two
  children under the same two `sample_times` entries compile to the same
  `(D, Φ)` pairs whether the container is transparent (`a`, `b`) or not
  (`kids/a`, `kids/b`), and the field-name sugar keys on the *field*, so
  `(children = Relative(2, 1),)` still applies one declaration to every element
  of a `Group`.
- **A computed boundary is the authored one.** The assembly whose two boundary
  declarations are `input_passthrough`/`output_passthrough` calls and the twin
  that writes both out by hand produce equal declaration tuples, equal root
  inputs, equal exported faces, equal flat lists, and equal port values at the
  exported face after `init!`. The helper is sugar, and the test is what says
  so: nothing downstream can tell which of the two it was handed. The twin is
  stateless and never runs — a trajectory would pin nothing the boundary-zero
  read does not.
- **The helpers own two refusals and no more.** `except` and `only` together,
  and a filter naming a face the child does not have (with the child's face
  list in hand), are the helper's errors, because nothing downstream could name
  the offender. Everything else is left where the rule already lives:
  `prefix = ""` colliding with a hand-written face is §8.6's face-name
  uniqueness, and a face both wired and passed through is §6.1's one-producer
  rule naming `child_connections` and `input_connections` as the two claimants.
- **`child_path` is a child path, one level, whatever names the child.** A
  deeper `child_path` meets the same rejection a wiring endpoint does, and a
  child living in a name-transparent container is addressed by its bare key —
  the helper never learns which kind of field held it.

- **A simulation owns one executor, and it is the one the loop runs.**
  `phase_bodies(sim) === sim.exec.bodies`, and `sim.exec.act` is
  `activation(sim.build, Float64)` itself: the bodies §7.5's measurements time
  are the executor's own rather than a re-derivation, and the executor carries
  the activation it was materialized from. `evaluate!(sim)` and
  `evaluate!(sim.exec)` leave the same derivative buffer, and the executor form
  allocates nothing — the `Simulation` method is a spelling, not a layer.

- **Five selectors, one address space, both activations.** Each member reads
  exactly what it names off an executor — the whole `x` leaf out of `xbuf`, the
  discrete `s` out of the component's own store, `f`'s output out of `ẋbuf`, a
  port's cell, a root input's cell, and an exported face's producer cell, the
  face reading identically to the port it aliases — and `i` indexes the read
  value. At `D8` the leaf types are the activation's throughout, while the
  frozen discrete store stays pinned `Float64` (§9.4, D-166): the same read set
  compiles at both scalars and the test runs it at both.
- **The gather twin costs nothing.** `gather(reader, executor)` allocates zero
  and infers a concrete NamedTuple, which is what makes a per-iteration read
  free against the sweep it follows. The empty read set gathers `(;)` — a
  service with nothing to read pays a no-op, not a special case.
- **One call, every violation.** A misspelled path, an undeclared port, a
  `get_deriv` on a discrete `s` and an unknown root face come back as one
  refusal naming all four, each by its own label and its selector as authored;
  a second call collects the assembly path, a root input read as a face, an
  index on a scalar leaf and an undeclared state field. §13.1's register, one
  register over from the condition algebra's.
- **A read set is a type, not a NamedTuple.** The bare spelling is refused with
  the directive `combine` gives for its own misuse, and a non-selector value in
  a `reads(…)` call is refused where it is written.
- **The source rule is enforced where the source is known.** A store selector
  in a binding's `reads` is `ReadBindingUnresolved` at attach, naming the
  remedy (declare the field public, read the published port); so is an indexed
  selector, the binding register reading whole cells. Both rejections leave the
  roster untouched, like every other attach refusal.
- **`capture` is the door back out, and it is total.** The captured pair
  re-establishes the world it was read from — `xbuf`, every `s` and `m` store,
  every root input cell and the clock — through an ordinary `init!` on a twin,
  at boundary zero and again after a trajectory; its coverage is the build's
  root inputs exactly, so nothing lies under it and §14.6's check passes by
  construction. The condition is time-free and `t` rides beside it, which is
  what `t₀ = t` spends. Legality is the §14 table's: `initialized` and
  `stopped`, with `built` refused for want of committed stores and `running`
  refused by the §11.3 freeze, all as one `ServiceLifecycle`.

- **One plan, one shape, many trees.** A plan compiled from a condition tree
  lands, from a *second* tree of that shape with other values everywhere, the
  same four homes the dynamic walk lands from that second tree — flat buffer,
  discrete store, mode store and both root-input cells compared whole. The plan
  holds lenses, not the values it was shown, and its type carries the tree type
  it was compiled from.
- **The specialized write is free.** `apply!(ex, plan, tree)` allocates zero
  over the four-home fixture — the prefix sweep, the flat-buffer write, both
  stores as whole values and both root-input scatters — with the tree handed in
  already built. Building it is the caller's cost, and it is *not* free here: an
  `at` node holds a `String`, so any tree carrying one is not isbits, and this
  fixture's construction measures 912 bytes — §14.2's "rebuilding the tree per
  trim iteration is stack-only construction" does not hold for a tree with `at`
  prefixes as things stand, which is recorded here rather than papered over.
  The register is measured where its own work is.
- **The store merge is the composite's.** A baseline authoring one field of a
  store and a patch authoring the other both survive the write, because the
  plan is compiled from the composite tree and the merge base is the declared
  defaults. A plan over the patch alone would reset the baseline's field to its
  default, which is the trap the composite exists to avoid.
- **Shape drift is structured, and nothing is written.** A tree of another type
  reaches the fallback and names both types; a tree of the right type with a
  different `at` prefix at the same position reaches the `===` sweep and names
  the position and both strings. Either way the executor is bit-for-bit what it
  was before the call.
- **The prefix sweep compares content, not pointers.** A prefix built afresh
  per call — a `String` at a different address, asserted to be one — passes the
  sweep of a plan compiled from the literal, because `===` on `String` is
  content equality. That is what a service rebuilding its tree per evaluation
  needs, and it is a stronger claim than "equal literals are one object".
- **The converter is the destination leaf type, and it is baked.** A plain
  `Float64` leaf into a seeded activation's `x` reads back with zero partials —
  the embedding that is exact for a value held at the operating point — and a
  leaf already at the activation's scalar reads back with its partials intact.
  A decision variable authored into a discrete `s`, frozen `Float64` at every
  activation, is refused at resolution with the clause that says why, and the
  nominal activation's own refusals are unchanged.
- **A compiled product belongs to one activation, by dispatch.** A `Float64`
  `Reader`, `ConditionPlan` and `SpecializedPlan` each refuse a `Dual` executor
  naming both scalars, and the executor is untouched. The refusal is an
  `InternalInvariant` and no diagnostic kind: §14.4 makes the pairing a
  framework invariant the services uphold, so reaching it is an internal
  assertion firing, outside the kind set on purpose (D-215).
  The scalar rides in the product's own type, so the pairing costs a method
  signature rather than a runtime test — the §7.5 gates measure zero unchanged.

- **A trim problem solves, commits, and the commit is an `init!`.** The
  one-step linear problem lands `u = (g/l)·sin θ` to its tolerance in four
  evaluations, and afterwards the simulation is `initialized` at the `t₀` it
  was given, with the authored attitude in its state and the solved torque in
  the root-input cell. A twin `init!`ed by hand with `override(baseline,
  condition(report.solution))` at the same anchor has identical buffers,
  stores, root inputs and clock — the commit is not *like* an init, it is one.
- **The nonlinear problem converges, and the box picks the branch.**
  `−(g/l)·sin θ + u = 0` has two roots in `[0, π]`; the declared box admits
  one, and the solve lands on it in four iterations, with the `Dual` decision
  leaf and the held `Float64` leaf meeting inside one fragment payload —
  §14.3's two converter arms in a single write.
- **Names pair, order never does.** The same two-decision problem with both
  bound spellings permuted *and* the tolerances permuted with them returns the
  same solution bit for bit; only the residual NamedTuple's key order follows
  the tolerances it was canonicalized to.
- **No convergence, no commit, and nothing moves.** An infeasible balance
  reports `converged == false` with the backend's status recorded and
  `committed_residuals === nothing` — the absence of a commit, not a flag. On a
  never-initialized simulation the lifecycle stays `built` and `run!` still
  refuses; on an initialized one every buffer equals its pre-call copy.
- **A saturated decision is named at the solution.** An actuator bound below
  the equilibrium torque, under a force balance declared to a tolerance the
  bound still meets, converges with `saturated == [(:u, :upper)]` and the
  bound's own value committed. The same problem unbounded names nothing.
- **The box is honored at every point the backend returns.** A guess of `100.0`
  under `[-1, 1]` comes back as `1.0`, saturated at `:upper`, with `1.0`
  committed: the guess is projected at the pack site, so the returns that never
  step — the already-within-tolerance one, a stall at iteration one — cannot
  hand back an out-of-box point as converged with `saturated` empty. A
  degenerate box (`lower == upper`) pins the decision to that one value, named
  `:lower` because the saturation check tests `lower` first.
- **The empty problem is the equilibrium probe.** `guess = (;)` bypasses the
  solver outright — `status == :bypassed`, one evaluation, no iterations, no
  seeded activation — and the nominal half's establishment round *is* the
  evaluation. On an equilibrium baseline it converges and commits; on one that
  is not, it answers no by the ordinary box test and leaves the sim `built`.
- **The setup diagnostic collects, in three observable stages.** A bounds
  key-set mismatch, an `Int` guess field and an unresolvable selector come back
  as one throw of three diagnostics — two `TrimProblemInvalid`s and the read
  set's own `TapResolution`, spliced in as the value it is; an inverted box (`lower` above
  `upper` on one decision) is collected there too, with both values named,
  because no projection can honor it, and so is a non-positive tolerance —
  zero and negative in one problem come back as one throw naming both, the
  acceptance test's divisor being the reason. The residual/tolerances key-set
  disagreement is the check only an evaluation can make (§14.7 says so), so it
  is reported from there, in its own collected throw of the same kind — and
  from *both* evaluations that can observe it: a lambda branching on the scalar
  it is handed answers one key set at the nominal guess and another at the
  first seeded point, and the seeded re-check turns what would be a bare
  `ErrorException` from the reorder into the same named refusal. Every refusal
  leaves the simulation untouched.
- **An incomplete baseline is `UninitializedInputs` at setup.** Trim's
  application to its own scratch stores is one of §14.6's three sites, and
  coverage is a plan-level fact: the check runs before any evaluation, names
  the uncovered face, and names `trim!` as the operation.
- **The scratch world's frozen cells are the authored world's (D-213).** With
  the pendulum's torque arriving from a discrete producer whose `s` the
  baseline authors, the solution is `asin(acc/(g/l))` — reachable only if the
  seeded world's frozen cell holds what that `s` publishes; the probe's
  synthesized zero would put it at `0`. Changing the authored `acc` moves the
  solution with it, which is what makes the copy load-bearing rather than
  incidental. And a decision variable authored *into* that frozen `s` is
  refused by the seeded compile's own `ConditionResolution`, with the
  frozen/pinned clause, before anything is written.
- **The commit is an ordinary boundary zero, and its movers are reported.** A
  guard the solved attitude already holds fires at `t₀` — derived, not asserted
  — and the report carries `[("trig", :fire)]` beside a `TrimCommitEvents`
  warning. A handler that only writes modes moves no residual, so the
  committed-state numbers are still the solved ones.
- **A residual mover shows in the committed-state residuals.** A handler that
  resets the solved attitude raises both warnings, in order, and
  `committed_residuals.torque` is the residual at the *reset* point, not the
  solved one. The verdict is not re-litigated: it gated the commit at the
  solved point, and both sets of numbers stand as reported.
- **A warm restart is `capture` handed back.** `run!`, then `(c, t) =
  capture(sim)`, then `trim!(sim, problem; baseline = c, t₀ = t)` commits with
  the clock at `t` — the resumed spelling §14.8 names, with continuity explicit
  rather than implied.
- **The per-iteration write and read are free.** At the service's own seeded
  scalar, `apply!(ex, plan, tree)` and `gather(reader, ex)` both measure zero
  allocations, and one seeded pass yields the residual value and its Jacobian
  column together (`value` and `partials` of the same gathered leaf). The
  residual lambda is the user's and is not gated; building the tree is not
  free either, for the `at`-prefix reason increment 21's part 3 records.
- **`trim!` is a stopped-sim service on a nominal deployment.** A
  `Simulation{D8}` is refused outright — the commit runs through boundary zero
  on the simulation's own stores, and those are the nominal world's — a
  non-problem value gets a directive rather than a `MethodError`, and `running`
  is the §11.3 freeze, refused as `ServiceLifecycle` exactly as `init!` and
  `capture` are.

- **A collecting register returns kind values, and the barrier throws once.**
  The condition resolver's five violations over one tree come back as five
  diagnostics in one `BuildError` — four `ConditionResolution`s discriminated
  by `reason` and one `DuplicateConditionLeaf` — and the count is read off
  `length(e.diagnostics)` rather than out of a rendered sentence. The read
  register's four come back as four `TapResolution`s in the read set's own
  declaration order, each carrying its label and the selector as authored.
  Every payload the old message text carried is now a field: the candidate
  list in hand, the tier that has no such store, the role that made a field
  ineligible, the seeded activation on the frozen-leaf refusal.
- **Provenance is payload, not prose.** A `combine` collision names both chains
  as the two entries of `d.provenance`, and an `override` patch's chain records
  the layer it overrode inside the one string the flattening built — asserted
  as values, which is what makes them a property rather than a wording.
- **Rendering is pinned in exactly one testset.** The carrier over two kinds ×
  two paths renders the kind names leading, the groups in first-appearance
  order, the paths sorted within a group and the count line above; a lone
  diagnostic renders on one line with no count. Beside it, one did-you-mean
  render showing the candidates the site held (carried, never ranked) and one
  remedy render showing the list in hand, and `logged(d)` for a warning kind.
  Nothing else in the suite may match rendered diagnostic text (§13.2).

## Stand-in retirement history

Every stand-in introduced through increment 5 was retired on 2026-08-20;
increment 6's one row — localized guards detecting at boundary resolution —
was retired by increment 7, and increment 9's publication row by increment 11,
both on 2026-08-25. The `set_slot!` row in the README's table entered with
increment 9; increment 10's device-staging row was retired by increment 12 on
2026-08-26, `attach!` returning the handle, and the two §11.8 presentation
rows entered with increment 13. Increment 14 retired both presentation rows
on 2026-08-26 and entered its own three: the harness register's cell, the
tail window, and the status's per-publication allocation. The harness-cell and
tail-window rows were retired the same day by the D-200/D-201 spec amendment,
leaving `set_slot!` and the status allocation. Increment 18 retired
`set_slot!` on 2026-08-27, the §14 services increment it was booked against:
root-input initial values now arrive through `init!`'s condition (stopped) or
`stage!` (running), and the function is gone. What survives it is `poke!` in
`test/utils.jl`, used at two sites where the *counterfactual to the drain* is
the point — a mid-trajectory direct write, which neither framework path can
express — and test equipment, not a stand-in. Increment 22's `LegacyMessage`
row entered with its stage 1 and was retired by its stage 4 on 2026-08-29, the
same increment, the string constructor having existed only so the suite stayed
green while the sites converted file by file. The status allocation is the only
row left.
