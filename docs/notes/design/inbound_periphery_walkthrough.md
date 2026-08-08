# The inbound periphery, from the ground up

*A companion explainer, not normative text. The ground truth is
`framework_spec.md` [§9.3](framework_spec.md#93-inbound-root-input-slots-claims-and-the-frozen-roster) (slots, claims, the roster freeze), [§9.4](framework_spec.md#94-inbound-per-device-staging-representation-and-the-drain) (staging and the drain), [§9.7](framework_spec.md#97-the-gui-write-path-port-resolution-peek-staging-contract) (GUI write path),
[§10.6](framework_spec.md#106-run-lifecycle-and-partial-advance) (harness cell) and decision rows 44, 93, 96, 106, 107. Written
2026-07-31 after the round-3 write-surface settlement; rewritten 2026-08-01
for the roster freeze (rows 106–107). If this document and the spec ever
disagree, the spec wins.*

Everything in [§9.3](framework_spec.md#93-inbound-root-input-slots-claims-and-the-frozen-roster)–[§9.4](framework_spec.md#94-inbound-per-device-staging-representation-and-the-drain) answers one question: **how do things outside the loop
feed values into a simulation that owns its data exclusively?** FlightCore's
answer was a lock around the live model. This design's answer is a small
pipeline of immutable handoffs, configured entirely while the simulation is
stopped and frozen for the duration of each run. Each concept below names one
station of that pipeline; they are introduced in dependency order, each with a
code shape, then one frame runs end to end, and the write-surface rule closes
the story.

## 1. Slots and faces: *what* can be written

The root assembly exports input faces ([§11.6](framework_spec.md#116-paths-wiring-and-exports)) — named, typed inputs that no
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
consumers address table cells by path or face, [§9.2](framework_spec.md#92-outbound-snapshot-publication)/[§14.4](framework_spec.md#144-two-application-registers-over-one-plan).)

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
admit through [§9.3](framework_spec.md#93-inbound-root-input-slots-claims-and-the-frozen-roster)'s four-part ordered check — **identity** (`AlreadyAttached`),
**affinity** (at most one `needs_calling_task` holder, `CallerTaskConflict`),
**interactivity** (at most one `interactive` device, `InteractiveConflict`),
**claims** (face exclusivity, `ClaimConflict`) — ordered so a failing later
check always names two *distinct* devices, then compile the staging shape
([§3](#3-batches-and-staging-cells-how-a-write-is-proposed)), add the entry.
A **claim** is exclusive ownership of a face: one writer per slot at any
time, and `detach!` releases the claims. The claim set is derived from the
**binding** — the declarative table that also drives the mapping:

```julia
joystick_binding = (stick_y  = (face = "elevator", expo = 0.6),
                    throttle = (face = "throttle",))
# ⇒ claims = ("elevator", "throttle")
```

Because the roster is frozen per run, so is the **partition of the face set**
it induces: every face is either claimed by exactly one device or belongs to
the shared interactive remainder ([§5](#5-registers-modes-of-use-not-more-machinery)), and that partition is a static,
printable fact of the run — who writes what, decided before the first frame.

One consequence is deliberate and worth stating early: **device death is not
detach**. A task that crashes, returns voluntarily, or loses its hardware
mid-run simply stops filling its cell; the [§10.2](framework_spec.md#102-loop-scheduling-wait-primitive-yields-thread-budget) heartbeat reports the death
by name, and the entry — claims included — persists to run end. The orphaned
slots hold their last-drained values, and the GUI renders the fact where the
user is looking ("claimed by `T16000M` — task dead", [§9.7](framework_spec.md#97-the-gui-write-path-port-resolution-peek-staging-contract)). Recovery is
between runs: stop, `detach!`, then `init!` for a fresh trajectory or
`replay!`-to-end + `run!` to continue the interrupted one ([§10.7](framework_spec.md#107-replay-the-trace-re-drives-the-ordinary-loop)).

## 3. Batches and staging cells: *how* a write is proposed

Nobody outside the loop ever assigns a slot. Instead, a writer produces a
**batch** and deposits it in its **staging cell**: a one-element atomic
mailbox, one per attached device, written by exactly one device task:

```julia
mutable struct StagingCell
    @atomic pending::Union{Nothing, Batch}    # the latest not-yet-applied batch
end
```

The batch's shape is **fixed at attach**: a positional tuple over the
writer's surface, `Union{Nothing, T}` per face, `nothing` meaning *not
touched this time*. Authors never build it by hand — `map_input` returns
sparse face ⇒ value pairs, and `stage!` normalizes them through the entry's
attach-compiled schema:

```julia
map_input(datum, binding)      # ⇒ (throttle = 0.70, elevator = -0.05)
# stage! normalizes against the schema ("elevator" ⇒ 1, "throttle" ⇒ 2):
(-0.05, 0.70)                  # positional over the claim set; nothing = untouched
```

**Staging** is depositing that tuple into your own cell. It happens at any
wall-clock moment, on the device's own task, and touches nothing the loop is
using. Staging three times between frames leaves one batch: the incoming
tuple **CAS-merges** into the pending one, positionally — `nothing` keeps
the pending value, anything else wins as newest. Merge is the only policy
because it is always correct: for a *complete* writer (a joystick: full
write-set every poll) merge and overwrite are provably the same operation,
while for a *sparse* writer (the GUI: only what the user touched) overwrite
would silently lose the untouched pending edits. Staged values are **levels,
never deltas** (`press_count = 17`, never `presses += 1`): levels are
idempotent and survive coalescing.

## 4. The drain: *when* proposals become slot values

At the top of each frame — and only there — the loop takes each cell's
contents atomically and applies it through the entry's attach-compiled
**scatter** (position → slot store, statically typed, `nothing` skips — the
mirror of [§9.2](framework_spec.md#92-outbound-snapshot-publication)'s output gather), in attachment order (the harness cell
last, [§5](#5-registers-modes-of-use-not-more-machinery)):

```julia
for entry in roster                                    # frame top, loop task
    batch = @atomicswap entry.cell.pending = nothing   # indivisible take
    batch === nothing && continue                      # nothing staged: fine
    scatter!(slots, entry, batch)                      # no checks — see §6
    record_in_trace!(frame, entry.id, batch)
end
```

Note what is *absent*: the drain validates nothing. Every check ran at
staging ([§6](#6-the-write-surface-rule-rows-44-and-106)), so the drain is pure application — and since the roster is a
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
One retention detail (row 107): an enumerated writer's batch enters the trace
verbatim — claim-narrow, dense by nature — while the interactive register's
wide, mostly-`nothing` tuple is converted on retention to sparse
(position ⇒ value) pairs, so trace size tracks information, not surface
width; replay converts back once, up front, off the loop ([§10.7](framework_spec.md#107-replay-the-trace-re-drives-the-ordinary-loop)).

## 5. "Registers": modes of use, not more machinery

When the spec says "register" here it means a *sanctioned usage mode* of the
one mechanism above (the word as in "speech register", not a hardware
register). There are two:

- **The claimed register**: enumerate faces in your binding, claim them at
  attach, own them exclusively. The joystick. Autonomous devices live here.
- **The interactive register**: the GUI ([§9.7](framework_spec.md#97-the-gui-write-path-port-resolution-peek-staging-contract)) and the harness/REPL's
  task-free entry point `stage!(sim, "face" => value, ...)` ([§10.6](framework_spec.md#106-run-lifecycle-and-partial-advance)). Its
  writers share the **unclaimed remainder** of the run's partition — the
  complement of the union of all claims, computed rather than staked, and
  every bit as static for the run as a claim set. It has its own compiled
  positional shape and scatter, over the unclaimed set, recompiled at each
  stopped-sim `attach!`/`detach!`. Residency here is declared, not granted: a
  binding elects the derived surface with the `interactive(b)` marker
  ([§9.6](framework_spec.md#96-devices-one-authoring-contract-no-taxonomy)), and at most one interactive device sits on a roster beside the
  always-present harness cell — a second candidate is an `InteractiveConflict`
  at attach ([§9.3](framework_spec.md#93-inbound-root-input-slots-claims-and-the-frozen-roster)). A second interactive front end therefore has a
  spelling; the shipped GUI is just the resident.

The one sequencing rule: among interactive writers the harness cell drains
**last** — the explicit hand of code beats a widget interaction. Sequencing
within one register, not cross-device arbitration; in practice the two
rarely coexist, since during a GUI run the calling task renders.

## 6. The write-surface rule (rows 44 and 106)

**Every writer has a write surface, and staging enforces it**: a batch entry
reaches a slot iff the named face is inside the writer's surface; anything
else is rejected in `stage!`'s normalization, on the writer's own task, with
a runtime warning. A surface arises in one of two ways — enumerated (a
device's claim set) or derived (the interactive remainder) — and under the
roster freeze both are static per run, which is why staging can be the
*only* enforcement point: there is no fact left that only the drain could
know. Two warning kinds cover the violations:

- `OutOfClaimEntry` — an enumerated-surface violation: the face has no
  position in the writer's schema. Writer id, face, discarded value, the
  claim set (plus the incumbent's id when the face is claimed elsewhere).
  "Your peer drifted from your binding."
- `ClaimedFaceEntry` — an interactive-register violation: the GUI or
  `stage!` named a face some device claims in this run's partition. Face,
  incumbent device id, discarded value, which interactive writer.

Worked illustration: attach a UDP telecommand peer whose binding enumerates
`("brake",)`; its remote peer sends `{"flaps": 0.5}` — a face it never
enumerated, currently unclaimed. `map_input` dutifully produces
`(flaps = 0.5,)`, and `stage!`'s normalization rejects it on the spot:
`flaps` has no position in the peer's schema → `OutOfClaimEntry`, attributed
to the device whose mapping drifted, before the value ever nears the loop.

The GUI is therefore not an exception: one device contract ([§9.6](framework_spec.md#96-devices-one-authoring-contract-no-taxonomy)), two
*binding* sides (enumerated vs. derived), one staging rule, one checkless
drain. Opportunistic writing by autonomous devices does not exist — a device
that wants a face enumerates it — so cross-writer races on one slot
structurally cannot arise: claimed faces have exclusivity, unclaimed faces
admit only the interactive register, and drain order is a diagnostic fact,
not an arbitration policy.

One seam remains, and it lives at the attach, not the drain: `stage!` works
while stopped, so a pending interactive batch may predate a stopped-sim
`attach!` that reshapes the unclaimed set. The attach renormalizes it —
reshape to the new schema, discard entries on newly-claimed faces with
`ClaimedFaceEntry` — so every run starts with cells matching the run's
schemas. The enumeration's other enforcement point is `attach!` itself
(`AttachUnknownFace`, `ClaimConflict`); the frozen roster is what extends
its authority, unchanged, to every frame of the run.

## 7. One frame, end to end

Joystick attached while stopped (claims `throttle`, `elevator`); GUI
attached (claims nothing); `run!` reads the roster, bakes widget liveness
([§9.7](framework_spec.md#97-the-gui-write-path-port-resolution-peek-staging-contract)) and specializes the drain; slots as above.

1. *Between frames*: the joystick task polls at its own rate, runs
   `map_input` (deadzone, expo — pure, on the device task), stages
   `(throttle = 0.70, elevator = -0.05)` — normalized to its two-position
   tuple; three polls happen before the next frame, so the cell holds only
   the newest. You drag the flaps slider; the GUI stages `flaps = 1.0` into
   *its* cell, at the flaps position of the interactive shape. The throttle
   slider renders read-only — its face is claimed, a fact baked at run
   start.
2. *Frame top*: drain. Both cells swap and scatter — no checks, both
   surfaces were enforced at staging. Slots are now
   `(throttle = 0.70, elevator = -0.05, flaps = 1.0, brake = 0.0)`, frozen
   for the frame. The joystick's batch enters the trace verbatim; the GUI's
   is retained as the sparse pair `(flaps ⇒ 1.0)`.
3. The sweep runs against those constants; the boundary completes; the
   snapshot publishes; the GUI's next render reads the snapshot and shows
   flaps at 1.0 — its own write, round-tripped, no model poking anywhere.
