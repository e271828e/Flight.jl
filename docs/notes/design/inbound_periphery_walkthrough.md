# The inbound periphery, from the ground up

*A companion explainer, not normative text. The ground truth is
`framework_spec.md` §12.3 (staging, claims, drain), §12.5 (GUI write path),
§12.11 (harness register) and decision rows 44, 93, 96. Written 2026-07-31,
after the round-3 write-surface settlement; if this document and the spec ever
disagree, the spec wins.*

Everything in §12.3 answers one question: **how do things outside the loop
feed values into a simulation that owns its data exclusively?** FlightCore's
answer was a lock around the live model. This design's answer is a small
pipeline of immutable handoffs. Each concept below names one station of that
pipeline; they are introduced in dependency order, each with a code shape,
then one frame runs end to end, and the write-surface rule closes the story.

## 1. Slots and faces: *what* can be written

The root assembly exports input faces (§13.6) — named, typed inputs that no
component inside the model produces. Say our aircraft's root exports four:

```julia
input_faces(build)  # ("throttle", "elevator", "flaps", "brake") — all Float64 here
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
consumers address table cells by path or face, §12.2/§16.4.)

## 2. Batches and staging cells: *how* a write is proposed

Nobody outside the loop ever assigns a slot. Instead, a writer produces a
**batch** — an immutable set of face ⇒ value pairs:

```julia
batch = (throttle = 0.70, elevator = -0.05)   # a joystick poll, conditioned
```

and deposits it in a **staging cell**: a one-element atomic mailbox, one per
attached device, written by exactly one device task:

```julia
mutable struct StagingCell
    @atomic pending::Union{Nothing, Batch}    # the latest not-yet-applied batch
end
```

**Staging** is depositing a batch into your own cell. It happens at any
wall-clock moment, on the device's own task, and touches nothing the loop is
using. Staging three times between frames leaves one batch: a *complete*
writer (a joystick: full write-set every poll) overwrites its cell —
coalescing to newest is the same zero-order-hold decision the loop makes
everywhere — while a *sparse* writer (the GUI: only what the user touched)
CAS-merges into its own cell, so an untouched slider never clobbers anything.
Staged values are **levels, never deltas** (`press_count = 17`, never
`presses += 1`): levels are idempotent and survive coalescing.

## 3. The roster, attaching, claiming: *who* may write

The **roster** is the framework's registry of currently attached devices — an
immutable array, republished by one atomic reference operation whenever it
changes. One entry per device:

```julia
struct RosterEntry
    id::DeviceId              # stable across sessions — trace provenance
    cell::StagingCell         # this device's mailbox
    claims::NTuple{N,String}  # face names this device owns
    seq::Int                  # attachment order
end
```

**Attaching** (`attach!(sim, device, binding)`) is: validate the binding's
face names against the root contract (unknown face → `AttachUnknownFace`),
register the **claims**, allocate the cell, republish the roster. A **claim**
is exclusive ownership of a face: one writer per slot at any time, a second
claimant is a `ClaimConflict` at attach, and detaching releases the claims.
The claim set is derived from the **binding** — the declarative table that
also drives the mapping:

```julia
joystick_binding = (stick_y  = (face = "elevator", expo = 0.6),
                    throttle = (face = "throttle",))
# ⇒ claims = ("elevator", "throttle")
```

Note what a claim is *for*: it is not access control on staging (any task can
put anything in its own cell — it's a private mailbox), it is an
**arbitration promise** checked where writes become real: at the drain.

## 4. The drain: *when* proposals become slot values

At the top of each frame — and only there — the loop acquire-loads the
current roster once and, in attachment order (the harness cell last, §6
below), atomically takes each cell's contents:

```julia
for entry in roster                                    # frame top, loop task
    batch = @atomicswap entry.cell.pending = nothing   # indivisible take
    batch === nothing && continue                      # nothing staged: fine
    for (face, value) in pairs(batch)
        face_in_surface(face, entry, roster) ?         # the write-surface rule, §7
            (slots[face] = value) : warn_and_discard(...)
    end
    record_in_trace!(frame, entry.id => batch)
end
```

**Draining** is that swap-and-apply. It is the single point where the
periphery's proposals become the frame's slot constants, and everything after
it — the sweep, the boundaries, the published snapshot — is a pure function
of the drained batches. That purity is what makes the **input trace** (the
per-frame sequence of device-tagged drained batches, plus the header's
initial state and slot values) a complete record: replay feeds the same
batches to the same drain and gets a bit-identical trajectory.

## 5. "Registers": modes of use, not more machinery

When the spec says "register" here it means a *sanctioned usage mode* of the
one mechanism above (the word as in "speech register", not a hardware
register). There are two:

- **The claimed register**: enumerate faces in your binding, claim them at
  attach, own them exclusively. The joystick. Autonomous devices live here.
- **The interactive register**: the GUI (§12.5) and the harness/REPL's
  task-free entry point `stage!(sim, "face" => value, ...)` (§12.11). Its
  writers stage to currently-unclaimed faces without owning them.

## 6. One frame, end to end

Joystick attached (claims `throttle`, `elevator`); GUI attached (claims
nothing); slots as above.

1. *Between frames*: the joystick task polls at its own rate, runs
   `map_input` (deadzone, expo — pure, on the device task), stages
   `(throttle = 0.70, elevator = -0.05)`; three polls happen before the next
   frame, so the cell holds only the newest. You drag the flaps slider; the
   GUI stages `(flaps = 1.0)` into *its* cell. The throttle slider renders
   read-only — its face is claimed.
2. *Frame top*: drain. Joystick batch: both faces are its claims → applied.
   GUI batch: `flaps` is unclaimed → applied. Slots are now
   `(throttle = 0.70, elevator = -0.05, flaps = 1.0, brake = 0.0)`, frozen
   for the frame. Both batches enter the trace, tagged.
3. The sweep runs against those constants; the boundary completes; the
   snapshot publishes; the GUI's next render reads the snapshot and shows
   flaps at 1.0 — its own write, round-tripped, no model poking anywhere.

## 7. The write-surface rule (the settled arbitration story, row 44)

**Every writer has a write surface, and the drain enforces it**: a batch
entry is applied iff the named face is inside the writer's surface at drain
time; anything else is discarded with a runtime warning. A surface arises in
one of two ways:

- **Enumerated** — an ordinary device's surface is its claim set. Static for
  the attachment, exclusively its own (claims are disjoint by construction),
  and binding-bounded even where nobody else is involved.
- **Derived** — the interactive register's surface is the currently-unclaimed
  face set: the complement of the union of all claims, computed rather than
  staked.

The GUI is therefore not an exception: one device kind (§12.4), two *binding*
kinds (enumerated vs. derived), one drain rule. Opportunistic writing by
autonomous devices does not exist — a device that wants a face enumerates
it — so cross-writer races on one slot structurally cannot arise: claimed
faces have exclusivity, unclaimed faces admit only the interactive register.
Drain order is a diagnostic fact, not an arbitration policy; the one
sequencing rule is that the harness cell drains **last** (the explicit hand
of code beats a widget interaction — sequencing within one register, and in
practice the two rarely coexist, since during a GUI run the calling task
renders).

The two surface kinds come with a **guarantee asymmetry**. An enumerated
surface is stable: nothing can take `throttle` from the joystick
mid-attachment. The derived surface is shared and precarious: any face in it
can vanish the moment a device claims it, and a batch staged before the claim
and drained after it is discarded — no one's bug, the surface moved.

Worked illustration: attach a UDP telecommand peer whose binding enumerates
`("brake",)`; its remote peer sends `{"flaps": 0.5}` — a face it never
enumerated, currently unclaimed. `map_input` dutifully produces
`(flaps = 0.5,)`, and the drain discards it: `flaps` is outside the writer's
enumerated surface. The discard warns with the drifted mapping named —

- `OutOfClaimEntry` — an enumerated-surface violation: writer id, face,
  discarded value, the claim set (plus the incumbent's id when the face is
  claimed elsewhere). "Your peer drifted from your binding."
- `StaleInteractiveEntry` — a derived-surface violation: the GUI or `stage!`
  staged to a face a device has claimed since (or holds). "The surface moved
  between staging and drain."

The enumeration's other enforcement point is `attach!` itself
(`AttachUnknownFace`, `ClaimConflict`); the drain's surface check is what
extends the roster's authority to every frame thereafter.
