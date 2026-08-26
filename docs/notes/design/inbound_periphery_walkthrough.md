# The inbound periphery, from the ground up

*A companion explainer, not normative text. The ground truth is
`framework_spec.md` [§11.3][s11-3] (slots, claims, the roster freeze), [§11.4][s11-4] (staging and the drain), [§11.7][s11-7] (GUI write path),
[§12.6][s12-6] (harness cell) and decision [D-044][d-044], [D-093][d-093], [D-096][d-096], [D-106][d-106], [D-174][d-174], [D-176][d-176], [D-177][d-177].
Written
2026-07-31 after the round-3 write-surface settlement; rewritten 2026-08-01
for the roster freeze ([D-106][d-106]–[D-107][d-107]) and 2026-08-11 for the surface
unification ([D-174][d-174]–[D-177][d-177]: the GUI as an ordinary enumerated writer, the
harness register as the sole derived surface, uniform sparse trace records,
declared binding traits). If this document and the spec ever
disagree, the spec wins.*

Everything in [§11.3][s11-3]–[§11.4][s11-4] answers one question: **how do things outside the loop
feed values into a simulation that owns its data exclusively?** FlightCore's
answer was a lock around the live model. This design's answer is a small
pipeline of immutable handoffs, configured entirely while the simulation is
stopped and frozen for the duration of each run. Each concept below names one
station of that pipeline; they are introduced in dependency order, each with a
code shape, then one frame runs end to end, and the write-surface rule closes
the story.

## 1. Slots and faces: *what* can be written

The root assembly exports input faces ([§8.6][s8-6]) — named, typed inputs that no
component inside the model produces. Say our aircraft's root exports four:

```julia
input_faces(world)  # ["throttle", "elevator", "flaps", "brake"] — all Float64 here
```

A **root input slot** is the storage backing one such face: a source cell of
the signal table. "Slot" is reserved for exactly these. Conceptually:

```julia
# inside the Simulation, part of the signal table the loop owns:
slots = (throttle = 0.62, elevator = -0.03, flaps = 0.0, brake = 0.0)
```

Two rules give slots their character. They are **constants within a frame**:
the sweep reads them like any other signal, and nothing may change them
between frame top and frame end. And they are **the only thing the periphery
may write** — the entire outside world (joysticks, network peers, GUI, your
REPL) influences the model *only* by proposing new slot values. The write
side addresses them by face *name* (`"throttle"`), never by structural path —
the root contract is the vocabulary. (The read side is different: snapshot
consumers address table cells by path or face, [§11.2][s11-2]/[§14.4][s14-4].)

## 2. The roster, attaching, claiming: *who* may write

The **roster** is the registry of attached devices, and the first thing to
know about it is *when* it can change: **only while the simulation is
stopped** (`built`, `initialized` or `stopped` — never `running`, and pause
is inside a run). `attach!` and `detach!` are configuration operations in
the same family as `init!` and trim; a run begins by reading the roster
once, as a plain immutable value, and nothing about it moves until the run
ends. One entry per device:

```julia
struct RosterEntry
    id::DeviceId              # stable across runs — trace provenance
    cell::StagingCell         # this device's mailbox (§3)
    claims::NTuple{N,String}  # face names this device owns
    schema::Schema            # face name → position, compiled at attach (§3)
end
```

**Attaching** (`attach!(sim, device, binding)`) is: validate the binding's
face names against the root contract (unknown face → `AttachUnknownFace`),
admit through [§11.3][s11-3]'s three-part ordered check — **identity** (`AlreadyAttached`),
**affinity** (at most one `needs_calling_task` holder, `CallerTaskConflict`),
**claims** (face exclusivity, `ClaimConflict`) — ordered so a failing later
check always names two *distinct* devices, then compile the staging shape
([section 3](#3-batches-and-staging-cells-how-a-write-is-proposed)), add the entry.
A **claim** is exclusive ownership of a face: one writer per slot at any
time, and `detach!` releases the claims. The claim set comes from the
**binding**, by one of two routes — *returned*, the declarative table that
also drives the mapping, or *computed*, when the binding declares
`is_greedy`:

```julia
joystick_binding = (stick_y  = (face = "elevator", expo = 0.6),
                    throttle = (face = "throttle",))
# ⇒ claims = ("elevator", "throttle")           # returned: enumerated by the table

is_greedy(::StandardGUIBinding) = true
# ⇒ claims = every root input face not already claimed, computed at attach
```

Past that point nothing can tell the two apart: the complement is disjoint
from every incumbent claim by construction, so it validates for exclusivity
like any other claim and gets the same schema, cell, drain and trace
([§11.6][s11-6]). A second greedy binding attached afterwards finds nothing left
and stakes the empty claim — legal, and reported as `EmptyGreedyClaim`.

Because the roster is frozen per run, so is the **partition of the face set**
it induces: every face is either claimed by exactly one device or belongs to
the harness remainder ([section 5](#5-registers-modes-of-use-not-more-machinery)), and that partition is a static,
printable fact of the run — who writes what, decided before the first frame.

One consequence is deliberate and worth stating early: **device death is not
detach**. A task that crashes, returns voluntarily, or loses its hardware
mid-run simply stops filling its cell; the [§12.2][s12-2] heartbeat reports the death
by name, and the entry — claims included — persists to run end. The orphaned
slots hold their last-drained values, and the GUI renders the fact where the
user is looking ("claimed by `T16000M` — task dead", [§11.7][s11-7]). Recovery is
between runs: stop, `detach!`, then `init!` for a fresh trajectory or
`replay!`-to-end + `run!` to continue the interrupted one ([§12.7][s12-7]).

## 3. Batches and staging cells: *how* a write is proposed

Nobody outside the loop ever assigns a slot. Instead, a writer produces a
**batch** and deposits it in its **staging cell**: a one-element atomic
mailbox, one per attached device, written by exactly one device task:

```julia
mutable struct StagingCell
    @atomic pending::Union{Nothing, Batch}    # the latest not-yet-applied batch
end
```

The batch's shape is **fixed at attach**: a pair of positional tuples over
the writer's surface — a concretely typed values tuple, and a parallel
`Bool` touched-mask whose set positions mean *staged this time*. Untouched
positions carry placeholders and are never read. Authors never build it by
hand — `map_input` returns sparse face ⇒ value pairs, and `stage!`
normalizes them through the entry's attach-compiled schema:

```julia
map_input(datum, binding)      # ⇒ (throttle = 0.70, elevator = -0.05)
# stage! normalizes against the schema ("elevator" ⇒ 1, "throttle" ⇒ 2):
values = (-0.05, 0.70)         # positional over the claim set
mask   = (true, true)          # false where the batch touched nothing
```

**Staging** is depositing that batch into your own cell. It happens at any
wall-clock moment, on the device's own task, and touches nothing the loop is
using. Staging three times between frames leaves one batch: the incoming
batch **CAS-merges** into the pending one, positionally — a clear mask
position keeps the pending value, a set one wins as newest. Merge is the only
policy because it is always correct: for a *complete* writer (a joystick: full
write-set every poll) merge and overwrite are provably the same operation,
while for a *sparse* writer (the GUI: only what the user touched) overwrite
would silently lose the untouched pending edits. Staged values are **levels,
never deltas** (`press_count = 17`, never `presses += 1`): levels are
idempotent and survive coalescing.

## 4. The drain: *when* proposals become slot values

At the top of each frame — and only there — the loop takes each cell's
contents atomically and applies it through the entry's attach-compiled
**scatter** (position → slot store, statically typed, unmasked positions
skipped — the mirror of [§11.2][s11-2]'s output gather), in attachment order
(the harness cell last, [section 5](#5-registers-modes-of-use-not-more-machinery)):

```julia
for entry in roster                                    # frame top, loop task
    batch = @atomicswap entry.cell.pending = nothing   # indivisible take
    batch === nothing && continue                      # nothing staged: fine
    scatter!(slots, entry, batch)                      # no checks — see §6
    record_in_trace!(frame, entry.id, batch)
end
```

Note what is *absent*: the drain validates nothing. Every check ran at
staging ([section 6](#6-the-write-surface-rule-d-044-d-106-and-d-174)), so the drain is pure application — and since the roster is a
fixed value at `run!`, the whole thing is compilable: the cells and their
scatters form a known tuple the frame function can specialize on, with no
name resolved and no dynamic dispatch at frame top.

**Draining** is that swap-and-apply. It is the single point where the
periphery's proposals become the frame's slot constants, and everything after
it — the sweep, the boundaries, the published snapshot — is a pure function
of the drained batches. That purity is what makes the **input trace** (the
per-frame sequence of device-tagged drained batches, plus the header's
initial state, slot values and per-writer schemas) a complete record: replay
feeds the same batches to the same drain and gets a bit-identical trajectory.
One retention detail ([D-176][d-176]): **every** drained batch is converted at the
drain into sparse (position ⇒ value) pairs against the writer's schema, so
trace size tracks information rather than surface width and every consumer
meets one record format; replay converts back once, up front, off the loop
([§12.7][s12-7]).

## 5. "Registers": modes of use, not more machinery

When the spec says "register" here it means a *sanctioned usage mode* of the
one mechanism above (the word as in "speech register", not a hardware
register). There are two:

- **The claimed register**: enumerate faces in your binding — or declare
  `is_greedy` and let the framework enumerate the remainder for you — claim
  them at attach, own them exclusively. The joystick lives here; so does the
  GUI ([§11.7][s11-7]), which is an ordinary claimed writer whose claim happens to
  be computed, and so do several interactive front ends at once when their
  claims are explicit and disjoint (a web console on the autopilot faces, a
  local GUI on the stick faces). The only thing still limited to one holder
  per roster is `needs_calling_task`.
- **The harness register**: the framework-owned task-free entry point
  `stage!(sim, "face" => value, ...)` and its always-present cell
  ([§12.6][s12-6]). This is the design's one **derived** surface — the unclaimed
  remainder of the run's partition, the complement of the union of all
  claims, computed rather than staked and every bit as static for the run as
  a claim set. It has its own compiled positional shape and scatter over that
  set, recompiled at each stopped-sim `attach!`/`detach!`, and nobody
  declares residency in it: it belongs to whoever is holding the REPL.

The one sequencing rule: the harness cell drains **last**. With every surface
disjoint the order is unobservable, so this is a convention that keeps traces
reading the same way, not an arbitration policy.

## 6. The write-surface rule (D-044, D-106 and D-174)

**Every writer has a write surface, and staging enforces it**: a batch entry
reaches a slot iff the named face is inside the writer's surface; anything
else is rejected in `stage!`'s normalization, on the writer's own task, with
a runtime warning. A device's surface is its claim set, however that set was
acquired; the harness register's is the derived remainder — and under the
roster freeze both are static per run, which is why staging can be the
*only* enforcement point: there is no fact left that only the drain could
know. Two warning kinds cover the violations:

- `OutOfClaimEntry` — a device-surface violation: the face has no
  position in the writer's schema. Writer id, face, discarded value, the
  claim set (plus the incumbent's id when the face is claimed elsewhere).
  "Your peer drifted from your binding." The GUI's out-of-surface writes are
  this kind too, its claim being an ordinary one — though its widgets cannot
  produce them, since liveness renders out-of-claim slots read-only.
- `ClaimedFaceEntry` — a harness-register violation, and only that:
  `stage!` named a face some device claims in this run's partition. Face,
  incumbent device id, discarded value, the site.

Worked illustration: attach a UDP telecommand peer whose binding enumerates
`("brake",)`; its remote peer sends `{"flaps": 0.5}` — a face it never
enumerated, currently unclaimed. `map_input` dutifully produces
`(flaps = 0.5,)`, and `stage!`'s normalization rejects it on the spot:
`flaps` has no position in the peer's schema → `OutOfClaimEntry`, attributed
to the device whose mapping drifted, before the value ever nears the loop.

The GUI is therefore not an exception: one device contract ([§11.6][s11-6]), one
claim mechanism with two claim *sources*, one staging rule, one checkless
drain. Opportunistic writing by autonomous devices does not exist — a device
that wants a face enumerates it, and greediness is an explicit declaration —
so cross-writer races on one slot
structurally cannot arise: every claim is exclusive whatever its source, the
unclaimed remainder admits only the harness, and drain order is a diagnostic
fact, not an arbitration policy.

One seam remains, and it lives at the attach, not the drain: `stage!` works
while stopped, so a pending harness batch may predate a stopped-sim
`attach!` that reshapes the unclaimed set. The attach renormalizes it —
reshape to the new schema, discard entries on newly-claimed faces with
`ClaimedFaceEntry` — so every run starts with cells matching the run's
schemas. The enumeration's other enforcement point is `attach!` itself
(`AttachUnknownFace`, `ClaimConflict`); the frozen roster is what extends
its authority, unchanged, to every frame of the run.

## 7. One frame, end to end

Joystick attached while stopped (claims `throttle`, `elevator`); GUI
attached next under the greedy binding (claim computed at that instant:
`flaps`, `brake`); `run!` reads the roster, bakes widget liveness
([§11.7][s11-7]) and specializes the drain; slots as above.

1. *Between frames*: the joystick task polls at its own rate, runs
   `map_input` (deadzone, expo — pure, on the device task), stages
   `(throttle = 0.70, elevator = -0.05)` — normalized to its two-position
   tuple; three polls happen before the next frame, so the cell holds only
   the newest. You drag the flaps slider; the GUI stages `flaps = 1.0` into
   *its* cell, at the flaps position of its own claim's shape. The throttle
   slider renders read-only — its face is claimed by the joystick and so
   falls outside the GUI's own claim, a fact baked at run start.
2. *Frame top*: drain. Both cells swap and scatter — no checks, both
   surfaces were enforced at staging. Slots are now
   `(throttle = 0.70, elevator = -0.05, flaps = 1.0, brake = 0.0)`, frozen
   for the frame. Both batches enter the trace in the one record format —
   the joystick's as `(throttle ⇒ 0.70, elevator ⇒ -0.05)`, the GUI's as
   `(flaps ⇒ 1.0)`.
3. The sweep runs against those constants; the boundary completes; the
   snapshot publishes; the GUI's next render reads the snapshot and shows
   flaps at 1.0 — its own write, round-tripped, no model poking anywhere.

<!-- citation link definitions — generated by tools/linkify.jl; do not edit -->
[d-044]: framework_decisions.md#d-044--slot-exclusivity-and-the-write-surface-rule
[d-093]: framework_decisions.md#d-093--spawn-device-tasks-per-run-not-per-attach
[d-096]: framework_decisions.md#d-096--define-the-harness-register-stagelateststep-duration
[d-106]: framework_decisions.md#d-106--freeze-the-device-roster-for-the-duration-of-a-run
[d-107]: framework_decisions.md#d-107--match-trace-record-density-to-batch-density-not-surface-width
[d-174]: framework_decisions.md#d-174--re-class-the-gui-as-an-ordinary-enumerated-writer
[d-176]: framework_decisions.md#d-176--unify-trace-retention-on-one-sparse-record-format
[d-177]: framework_decisions.md#d-177--re-found-the-periphery-on-mandatory-roots-plus-declared-traits
[s11-2]: framework_spec.md#112-outbound-snapshot-publication
[s11-3]: framework_spec.md#113-inbound-root-input-slots-claims-and-the-frozen-roster
[s11-4]: framework_spec.md#114-inbound-per-device-staging-representation-and-the-drain
[s11-6]: framework_spec.md#116-devices-one-authoring-contract-no-taxonomy
[s11-7]: framework_spec.md#117-the-gui-write-path-port-resolution-peek-staging-contract
[s12-2]: framework_spec.md#122-loop-scheduling-wait-primitive-yields-thread-budget
[s12-6]: framework_spec.md#126-run-lifecycle-and-partial-advance
[s12-7]: framework_spec.md#127-replay-the-trace-re-drives-the-ordinary-loop
[s14-4]: framework_spec.md#144-two-application-registers-over-one-plan
[s8-6]: framework_spec.md#86-paths-wiring-and-faces
