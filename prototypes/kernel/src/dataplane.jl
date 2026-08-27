# The data plane (§11): staging cells with their frame-top drain — planes 1
# and 2 of §11.1's replacement for the shared mutable model — the diagnostic
# channel with the published framework status (§11.8, §13.2), snapshot
# publication, plane 3, and the log riding behind it (§11.2). The cells'
# owners — the roster's device entries and the harness register beside them —
# live in roster.jl, the device tasks that stage into them in devices.jl;
# what stands further out in the spec — the trace, the pacer diagnostics —
# is deliberately absent (README).
#
# This file holds the types and the pure mechanics; the `Simulation`-facing
# surface — `stage!`, `drain!`, `publish!`, `latest` — lives in sim.jl, beside
# the loop that calls it.

"""
The staging cell (§11.4): one atomic reference holding the pending batch — CAS
merge on the writer's side, `atomicswap` at the drain. The batch rides boxed in
a `Ref` so every handoff is one pointer-width atomic operation (§11.1) whatever
the surface's size, and the CAS compares batch *identity*: a failed replace
means the loaded batch was intercepted — drained, or merged past by another
stager — and the retry re-reads exactly what is pending now. The failure case
is precisely correct: an intercepted batch is already applied (or already
carries the other stager's entries), and must not be re-staged.
"""
mutable struct StagingCell{B}
    @atomic pending::Union{Nothing,Base.RefValue{B}}
end

# --- the diagnostic channel (§11.8, §13.2, Appendix C) -------------------------

"""
The prototype's closed kind set for the runtime warning stream (§13.2,
Appendix C): each kind is a Julia type, its identity, and its payload is
plain data — paths and names as strings and symbols, never component
instances; the declared/observed *port* types are the payload exception, and
they are small. These are the kinds whose sources the prototype has built;
the four whose features are absent — `DebtReanchor`, `ThreadBudget`,
`ReplayDiscardedStaging`, `UnboundedRun` — are absent with them (README).
Writer attribution is never a payload field: the channel is per-writer, so
the cell supplies it (§11.8, §12.4: no call passes a device id).
`DeviceJoinTimeout`'s `who` is not that attribution — it is the payload's
subject per Appendix C, the abandoned device; the *writer* is the loop,
whose own cell carries it (§12.4, D-203).

`MalformedDatum` is the author's kind: a datum that could not be mapped for
environmental reasons — a truncated datagram, malformed JSON, an
out-of-range field — reported by the author's own loop body via
`report!(handle, MalformedDatum(cause))` after catching its own parser error
(§11.6: catch, stage nothing, report, continue). The classification is the
author's — only they know their parser — and any *other* exception
propagates to the wrapper as the `DeviceCrash` it is. The rest are the
framework's: the three staging kinds (§11.4 — the write-surface and entry
violations, every check at staging), the two budget degradations (§10.4,
§10.6, on the loop's own cell), and the device crash (§12.4, from the
wrapper on the device's task, or from the calling task at the init bracket).
"""
struct MalformedDatum
    cause::Any
end

"§11.4's out-of-schema entry: no position in the writer's compiled surface."
struct OutOfClaimEntry
    face::Symbol
    value::Any                      # the discarded value
    surface::Vector{Symbol}         # the writer's claim set, the list-in-hand
    incumbent::Union{Nothing,String} # who claims the face, when someone does
end

"§11.3's harness write to a claimed face, naming the incumbent; `site` is `:staging` or `:renormalization`."
struct ClaimedFaceEntry
    face::Symbol
    incumbent::String
    value::Any
    site::Symbol
end

"§11.4's unconvertible value, rejected at staging for every writer."
struct EntryTypeMismatch
    face::Symbol
    value::Any
    declared::Any                   # the root input's declared type
end

"§10.4's localization-budget exhaustion: further crossings this frame fire at boundary resolution."
struct ChatteringBudget
    path::String
    event::Symbol
    t::Float64
    budget::Int
    count::Int
end

"§10.6's firing-budget exhaustion: the event's further edges at this boundary are lost."
struct FiringBudget
    path::String
    event::Symbol
    t::Float64
    budget::Int
    count::Int
end

"§12.4's device failure — the wrapper's catch, or the init bracket's; `abort` is the attachment's `should_abort`."
struct DeviceCrash
    cause::Any
    abort::Bool
end

"""
§12.4(5)'s join-timeout abandonment, written to the loop's own cell at the
tail and collected by the run's-end sweep into the termination record
(D-203). `who` names the abandoned device — Appendix C's payload, beside the
cap and the final snapshot's boundary time and index at shutdown.
"""
struct DeviceJoinTimeout
    who::String
    timeout::Float64
    t::Float64
    boundary::Int
end

"The closed set as a union: what a ring holds, and what `_report!` admits."
const DiagValue = Union{MalformedDatum,OutOfClaimEntry,ClaimedFaceEntry,
                        EntryTypeMismatch,ChatteringBudget,FiringBudget,
                        DeviceCrash,DeviceJoinTimeout}

"""
The per-kind counter record (§11.8): a **fixed-shape isbits record, never a
`Dict`** — licensed by the closed kind set, which makes the counter layout a
type rather than a lookup. One field per kind, in the kinds' declaration
order; the same shape serves a batch's suppressed counts and the loop's
cumulative totals, and `+` is the fold between them.
"""
struct KindCounts
    malformed::Int
    out_of_claim::Int
    claimed_face::Int
    type_mismatch::Int
    chattering::Int
    firing::Int
    crash::Int
    join_timeout::Int
end
KindCounts() = KindCounts(0, 0, 0, 0, 0, 0, 0, 0)

_kind(::MalformedDatum)   = :malformed
_kind(::OutOfClaimEntry)  = :out_of_claim
_kind(::ClaimedFaceEntry) = :claimed_face
_kind(::EntryTypeMismatch) = :type_mismatch
_kind(::ChatteringBudget) = :chattering
_kind(::FiringBudget)     = :firing
_kind(::DeviceCrash)      = :crash
_kind(::DeviceJoinTimeout) = :join_timeout

_bump(c::KindCounts, k::Symbol) =
    KindCounts((getfield(c, f) + (f === k) for f in fieldnames(KindCounts))...)
Base.:+(a::KindCounts, b::KindCounts) =
    KindCounts((getfield(a, f) + getfield(b, f) for f in fieldnames(KindCounts))...)
_total(c::KindCounts) = sum(f -> getfield(c, f), fieldnames(KindCounts))

"The ring's capacity (§11.8): the channel's normative bound, which *is* the rate limit."
const DIAG_RING = 16

"""
One frame's bounded accumulation (§11.8): a small ring of diagnostic values,
capacity `DIAG_RING`, plus the per-kind counts of what the ring could not
hold. The drop policy is earliest-in-frame retained, excess becomes counts:
the first occurrences are the ones with diagnostic content, the hundredth is
noise the count already reports. Immutable, because the CAS handoff below
needs value semantics: an append builds a new batch, and a drained batch is
frozen — never written again — with the writer allocating afresh at its next
emission, on its own task (§11.8's GC-over-reuse trade).
"""
struct DiagBatch
    ring::Vector{DiagValue}
    suppressed::KindCounts
end

"The shared empty sentinel (§11.8): what the drain swaps *in*, and what a quiet frame gets back."
const EMPTY_DIAG = DiagBatch(DiagValue[], KindCounts())

"""
The diagnostic cell (§11.8): one per writer — each rostered device's, the
harness register's, the loop's own — single-writer, the same ownership
argument as the staging cells: no locking, no arbitration, no new primitive.
The CAS mirrors `_stage!`'s: a failed replace means the loaded batch was
intercepted by the drain, and the retry re-reads what is pending now — the
intercepted entries are already taken, and the entry in hand merges against
the fresh (sentinel) state.

The liveness heartbeat rides in the same cell (§11.8, §12.2), as the atomic
timestamp the device task stores on every loop pass from inside the handle
primitives and the loop acquire-loads where it publishes. It is a field,
always present, never a diagnostic kind; its initial `0.0` is what makes a
never-heartbeated cell — a failed `init!`'s — read stale against the §12.2
threshold from the first frame, with no marking machinery (§12.4).
"""
mutable struct DiagCell
    @atomic batch::DiagBatch
    @atomic heartbeat::Float64
end
DiagCell(b::DiagBatch) = DiagCell(b, 0.0)

# The writer's side (§11.8), on the writer's own task: append under the bound,
# or count past it by kind. Reached through `report!(handle, …)` (devices.jl)
# and from the framework's own emission sites.
function _report!(cell::DiagCell, d::DiagValue)
    while true
        cur = @atomic cell.batch
        next = length(cur.ring) < DIAG_RING ?
            DiagBatch(push!(copy(cur.ring), d), cur.suppressed) :
            DiagBatch(cur.ring, _bump(cur.suppressed, _kind(d)))
        (; success) = @atomicreplace cell.batch cur => next
        success && return nothing
    end
end

# The indivisible take (§11.8): one `atomicswap` per cell at frame top, at the
# same point and under the same argument as the staging drain. A quiet frame
# swaps the sentinel in and gets the sentinel back — no allocation, and no
# load-only code path that goes untested on healthy runs.
_take!(cell::DiagCell) = @atomicswap cell.batch = EMPTY_DIAG

# The heartbeat's two sides (§11.8, §12.2): the store from inside the handle
# primitives, the acquire-load where the loop publishes.
_beat!(cell::DiagCell) = (@atomic :release cell.heartbeat = time(); nothing)
_heartbeat(cell::DiagCell) = @atomic :acquire cell.heartbeat

"The shared empty delta (§11.8): a quiet writer's `recent`, one object for all of them."
const EMPTY_RECENT = DiagValue[]

"""
The loop's private account for one writer (§11.8): everything on the loop's
side of the cell, owned by the loop task alone. `recent` and `suppressed` are
the pending delta — what the frame-top drain took and no snapshot has carried
yet; the first publication after the drain takes them, so a reader at any
cadence sees each occurrence exactly once in `recent`. `totals` is the
cumulative per-kind record since the run began, reset at each run's top and
*copied* into every status (`KindCounts` being isbits, the read is the copy):
delta plus total is what makes the status legible at any reading cadence —
nothing is lost by not looking.
"""
mutable struct WriterAccount
    recent::Vector{DiagValue}
    suppressed::KindCounts
    totals::KindCounts
end
WriterAccount() = WriterAccount(EMPTY_RECENT, KindCounts(), KindCounts())

function _reset!(a::WriterAccount)
    a.recent = EMPTY_RECENT
    a.suppressed = KindCounts()
    a.totals = KindCounts()
    nothing
end

# The frame-top fold (§11.8): the indivisible take, then — the batch
# exclusively the loop's — retained values into the pending delta and every
# occurrence, retained and suppressed alike, into the totals. The quiet path
# is the sentinel coming back: no allocation, nothing touched.
function _fold!(a::WriterAccount, cell::DiagCell)
    batch = _take!(cell)
    batch === EMPTY_DIAG && return nothing
    counts = batch.suppressed
    for d in batch.ring
        counts = _bump(counts, _kind(d))
    end
    a.recent = a.recent === EMPTY_RECENT ? copy(batch.ring) : append!(a.recent, batch.ring)
    a.suppressed = a.suppressed + batch.suppressed
    a.totals = a.totals + counts
    nothing
end

"""
One writer's record in the published framework status (§11.8): `recent` — the
ring the frame's drain took, at most `DIAG_RING` entries, riding exactly one
snapshot; `suppressed` — the per-kind counts the ring refused alongside them;
`totals` — the cumulative account since the run began, monotone across
snapshots (so log decimation loses *which* boundary an occurrence fell on,
never *how many*); and, for a rostered device, the liveness `heartbeat`
(§12.2) beside the `task_state` the loop reads off the run's `Task` handle at
publication (D-193) — `:none` when no task exists (a failed `init!`, or a
stopped sim), `:running`, `:done`, or `:failed`. The harness register's and
the loop's own records carry `nothing` for both: no task of their own to be
alive or dead.
"""
struct WriterStatus
    who::String
    recent::Vector{DiagValue}
    suppressed::KindCounts
    totals::KindCounts
    heartbeat::Union{Nothing,Float64}
    task_state::Union{Nothing,Symbol}
end

"""
The published framework status (§11.8, §11.2): a concrete frozen value, never
a window onto live bookkeeping — per-writer records in the drain's own order,
each rostered device in attachment order, then the harness register, then the
loop itself. Built fresh in `publish!` and frozen into the snapshot; the
binding rule holds because everything here is either immutable, a copy, or a
vector the account has released.
"""
struct FrameworkStatus
    writers::Vector{WriterStatus}
end

"§12.2's staleness threshold, in seconds of wall clock: advisory, a display rule, never a kill trigger."
const STALE_S = 2.0

"""
    stale(w::WriterStatus; now = time())

§12.2's liveness read: a device whose heartbeat is more than `STALE_S` behind
wall clock — deliberately loose, tolerating a device legitimately parked in a
blocking read between rare data. A starved, blocked or crashed device shows
as a stale heartbeat with a name on it, not as mysteriously frozen physics;
the never-heartbeated cell's `0.0` reads stale unconditionally. `false` for
the harness's and the loop's records, which have no heartbeat to judge.
"""
stale(w::WriterStatus; now::Float64 = time()) =
    w.heartbeat !== nothing && now - w.heartbeat > STALE_S

_task_state(::Nothing) = :none
_task_state(t::Task) = istaskfailed(t) ? :failed : istaskdone(t) ? :done : :running

"""
The batch (§11.4, D-202): a pair of positional tuples over one writer's
surface — a values tuple, concretely typed `Tuple{T₁, …, Tₙ}`, and a parallel
`Bool` touched-mask. A set mask position means *staged this time*; a clear one
means *not touched*, never "reset" — its value is the placeholder the writer
was compiled with, dead behind the mask guard. The concrete layout is the
point (D-202): a `Union{Nothing,Tᵢ}` carrier reads more directly, but tuple
covariance makes each touched-face combination its own concrete type, and the
frame-top scatter would specialize per pattern — mid-run compilation and a
boxed dispatch where the drain promises pure application.
"""
struct Batch{V<:Tuple,M<:Tuple}
    vals::V                  # position → the staged (or placeholder) value
    mask::M                  # position → touched this time (NTuple{n,Bool})
end

"""
The compiled writer (§11.4): one write surface's staged representation, fixed
at a stopped-sim point — a device's compiled at attach against its claim set,
the harness register's at deployment binding and at every roster change
against its *derived* surface, the unclaimed complement (§11.3). The batch is
the values-plus-mask pair above (D-202), isbits with one concrete layout per
writer. The face-name → position schema, the root-input types, the per-position
cell addresses — the compiled scatter's data — and the blank batch the shim
starts from (placeholders drawn from the layout's probe values, mask all clear)
live here beside the staging cell they describe; a roster entry carries one of
these per device (roster.jl), and the harness register carries the one whose
shape the framework derives rather than receives.
"""
struct Writer{B<:Batch,A<:Tuple}
    faces::Vector{Symbol}    # position → face name: the schema
    types::Vector{Any}       # position → the root input's declared type Tᵢ
    addrs::A                 # position → the root input's cell address
    blank::B                 # placeholders under an all-clear mask
    cell::StagingCell{B}
end

_port_type(::CellAddr{P,K}) where {P,K} = P

function Writer(layout::Layout, faces::Vector{Symbol})
    addrs = Tuple(layout.addr[("", f)] for f in faces)
    types = Any[_port_type(a) for a in addrs]
    probes = Dict(f => v for (f, v) in layout.root_inputs)
    vals = Tuple(convert(types[i], probes[faces[i]]) for i in eachindex(faces))
    blank = Batch(convert(Tuple{types...}, vals), ntuple(_ -> false, length(faces)))
    Writer{typeof(blank),typeof(addrs)}(faces, types, addrs, blank,
                                        StagingCell{typeof(blank)}(nothing))
end

"""
§11.4's shim, run at staging on the writer's own side: the author's face ⇒
value pairs become a batch — name → position, convert to the root input's declared
type, set the mask. Every diagnostic site follows the compilation here, to
staging: face validity, surface membership and value convertibility are all
static facts of the run, so a face with no position in the schema is
discarded under an `OutOfClaimEntry`/`ClaimedFaceEntry` and an unconvertible
value under `EntryTypeMismatch`, while the rest of the batch stands. Nothing
remains at the drain. The diagnostics are written into the writer's own cell
(§11.8) — the device's for a handle staging, the harness register's
otherwise — on the staging task, surfacing in the status at the next frame
top's drain; a batch staged while stopped waits in the cell exactly as its
entries wait in the staging cell.

The out-of-schema kind discriminates by writer (§11.3, Appendix C): a
device's entry is always `OutOfClaimEntry` — naming the incumbent when the
face is claimed elsewhere — while the harness register's is `ClaimedFaceEntry`
naming the incumbent when a rostered claim covers the face, and
`OutOfClaimEntry` only when the face names no root input at all. `claimedby` is
the exclusivity index the roster maintains; `device` identifies a device
writer and `nothing` the harness; `site` distinguishes ordinary staging from
an attach's renormalization in the `ClaimedFaceEntry` payload. Returns
`nothing` for a batch with no surviving entry.
"""
function _normalize(w::Writer, pairs, claimedby::Dict{Symbol,String},
                    cell::DiagCell; device::Union{Nothing,String} = nothing,
                    site::Symbol = :staging)
    vals = Any[w.blank.vals...]
    mask = fill(false, length(w.faces))
    for (face, v) in pairs
        s = Symbol(face)
        i = findfirst(==(s), w.faces)
        if i === nothing
            incumbent = get(claimedby, s, nothing)
            if device !== nothing
                _report!(cell, OutOfClaimEntry(s, v, w.faces, incumbent))
            elseif incumbent !== nothing
                _report!(cell, ClaimedFaceEntry(s, incumbent, v, site))
            else
                _report!(cell, OutOfClaimEntry(s, v, w.faces, nothing))
            end
            continue
        end
        vals[i] = try
            convert(w.types[i], v)
        catch
            _report!(cell, EntryTypeMismatch(s, v, w.types[i]))
            continue
        end
        mask[i] = true
    end
    any(mask) || return nothing
    Batch(convert(typeof(w.blank.vals), (vals...,)), (mask...,))
end

_facelist(faces) = isempty(faces) ? "empty" : "{$(join(faces, ", "))}"

# The one coalescing policy (§11.4): merge, newest wins per face. Untouched
# faces survive; re-staged faces take the newest level — the per-face ZOH.
# Positional and mask-driven, unrolled by generation: both arguments share the
# writer's one concrete batch type, and the unroll leans on no small-tuple
# heuristic, so width does not degrade it (D-202).
@generated function _merge(pending::Batch{V,M}, incoming::Batch{V,M}) where {V,M}
    vals = [:(incoming.mask[$i] ? incoming.vals[$i] : pending.vals[$i])
            for i in 1:fieldcount(M)]
    mask = [:(pending.mask[$i] | incoming.mask[$i]) for i in 1:fieldcount(M)]
    :(Batch{V,M}(($(vals...),), ($(mask...),)))
end

# The CAS merge loop (§11.4), on the writer's task.
function _stage!(w::Writer{B}, batch::B) where {B}
    cell = w.cell
    while true
        pending = @atomic cell.pending
        merged = pending === nothing ? batch : _merge(pending[], batch)
        (; success) = @atomicreplace cell.pending pending => Base.RefValue{B}(merged)
        success && return nothing
    end
end

# The drain's application (§11.4): the compiled scatter, position → the root
# input's cell, statically typed, masked-off positions skipped — pure
# application, no checks. One specialization per writer, compiled at the
# stopped-sim point that compiled the writer, whatever the batch touches (D-202).
@generated function _apply!(store, addrs::Tuple, batch::Batch)
    stmts = [:(batch.mask[$i] && scatter!(store, addrs[$i], batch.vals[$i]))
             for i in 1:fieldcount(fieldtype(batch, :mask))]
    quote
        $(stmts...)
        nothing
    end
end

# One cell's drain (§11.4): the indivisible take by `atomicswap`, then the
# compiled scatter. Reached at frame top through the stopped-sim-compiled
# thunks (roster.jl's `_drain_thunk`), inside which both arguments are
# concrete.
function _drain!(store, w::Writer)
    ref = @atomicswap w.cell.pending = nothing
    ref === nothing && return nothing
    _apply!(store, w.addrs, ref[])
    nothing
end

# --- publication (§11.2) -------------------------------------------------------

"""
A snapshot (§11.2): the boundary-consistent signal table — the *whole* table,
root inputs riding along as the source cells they are — with `t`, the frame
ordinal and the §12.3 boundary ordinal, the counter-home rule: the index
rides *in* the snapshot, so any holder of one — the log, a post-run
inspector — indexes it without consulting the loop. The framework status
rides beside the table (§11.8): the per-writer diagnostic account, frozen.
Built in private memory by copying the cell buffers, then handed out through
one release-store; nothing reachable from a published snapshot is ever
written again, which is what makes the lock-free read sound. The state
stores (`x`, `s`, `m`) are deliberately not carried (§11.2). Read it with
`port(snap, path, name)`, addressed exactly as the live table.
"""
struct Snapshot{T,S<:StoreBundle}
    t::T
    frame::Int
    boundary::Int     # the §12.3 published-boundary ordinal; boundary zero = 0
    store::S
    layout::Layout    # build-frozen addressing, shared, never copied
    status::FrameworkStatus
end

port(s::Snapshot, path::String, name::Symbol) = gather(s.store, s.layout.addr[(path, name)])

# One boundary's capture: fresh buffers, one allocation per boundary — the
# framework side of §7.5's scope, which carved publication and logging out.
capture(b::StoreBundle) = StoreBundle(map(cs -> CellStore(copy(cs.buf)), b.stores))

"""
§11.2's `@atomic latest` reference in its own mutable holder, so the
`Simulation` itself stays immutable: release-store at publication,
acquire-load at `latest(sim)`. The exchange is wait-free in both directions —
a wedged reader cannot delay publication, and the loop cannot tear a reader's
view.
"""
mutable struct Published
    @atomic latest::Union{Nothing,Snapshot}
end

# --- the log (§11.2) -----------------------------------------------------------

"""
The log (§11.2): logging dissolves into publication. What is retained are the
published snapshots themselves — the same objects, zero extra copies (D-023
rejected preallocated buffers; unlogged snapshots die young) — under three
policies validated at deployment binding: the plain switch, the keep-every-kth
stride `log_every`, and the retention bound `log_max` (`Inf` stored as
`typemax(Int)`). Normative are the guarantees, not the mechanism below: the
bound is respected *continuously*, the retained count never exceeding
`log_max`; coverage is global at the effective stride `log_every · 2^k` after
`k` generations; and the two endpoints are kept unconditionally.

The endpoints ride outside the bounded middle — two extra references, never
counted against it. `first` is the boundary-zero snapshot (§14.5); `last` is
the latest published boundary — which at any stopped moment is the tail's
final snapshot, published before the sticky status was set (§12.4), the
terminal snapshot the spec retains unconditionally.

The middle re-decimates progressively (D-137: a rolling window was rejected —
what coarsens is density, never extent), amortized: `snaps` holds one
generation's retained snapshots, `nothing` marking a released slot until the
once-per-generation compaction. The arithmetic rests on one invariant — at
each generation's start the compacted vector holds the boundaries at ordinals
`stride · (1..max)`, index by index — so the entries the doubled stride
abandons are exactly the odd indices, and `cursor` walks them.
"""
mutable struct SnapshotLog
    enabled::Bool
    every::Int                      # log_every: the authored stride
    max::Int                        # log_max: what bounds `live`, never the endpoints
    stride::Int                     # the effective stride, every · 2^generation
    nb::Int                         # published-boundary ordinal; boundary zero = 0
    first::Union{Nothing,Snapshot}  # the boundary-zero endpoint (§14.5)
    last::Union{Nothing,Snapshot}   # the terminal endpoint: the latest published boundary
    snaps::Vector{Union{Nothing,Snapshot}}   # the bounded middle; `nothing` = released
    live::Int                       # retained middles: the count `log_max` bounds
    cursor::Int                     # the thinning cursor over odd indices; 0 = inactive
end

SnapshotLog(enabled::Bool, every::Int, max::Int) =
    SnapshotLog(enabled, every, max, every, 0, nothing, nothing,
                Union{Nothing,Snapshot}[], 0, 0)

# A warm restart is a new trajectory (§10.6's register reset, carried through):
# the log starts over, its boundary zero a fresh first endpoint.
function _reset!(L::SnapshotLog)
    L.stride = L.every
    L.nb = 0
    L.first = nothing
    L.last = nothing
    empty!(L.snaps)
    L.live = 0
    L.cursor = 0
    nothing
end

"""
One published boundary enters the log (§11.2), on the loop task, right behind
the release-store: boundary zero lands in `first`, a stride multiple is
retained into the middle, and every snapshot re-points `last` — one field
store, which is all the terminal endpoint costs. Off, the switch retains
nothing at all: retention is what it gates, publication being upstream of it.
"""
function log!(L::SnapshotLog, snap::Snapshot)
    L.enabled || return nothing
    if L.nb == 0
        L.first = snap
    elseif L.nb % L.stride == 0
        _retain!(L, snap, L.nb)
    end
    L.last = snap
    L.nb += 1
    nothing
end

# One retained append (§11.2, D-137), the amortized re-decimation in full. At
# the fill point — middle at `max`, no thinning under way — the stride doubles
# *immediately*, and the candidate in hand re-tests against it, half the old
# stride's multiples being no longer retained; the appends that follow continue
# the surviving even multiples seamlessly. While thinning is active each
# retained append first releases one predecessor at the cursor, so the bound
# holds continuously and a generation's thinning completes exactly when its
# refill does; compaction then runs, once per generation, restoring the
# index-by-ordinal invariant for the next fill.
function _retain!(L::SnapshotLog, snap::Snapshot, nb::Int)
    if L.cursor == 0 && L.live == L.max
        L.stride *= 2
        L.cursor = 1
        nb % L.stride == 0 || return nothing
    end
    if L.cursor > 0
        L.snaps[L.cursor] = nothing
        L.live -= 1
        L.cursor += 2
        L.cursor > L.max && (filter!(!isnothing, L.snaps); L.cursor = 0)
    end
    push!(L.snaps, snap)
    L.live += 1
    nothing
end
