# The data plane (§11): staging cells with their frame-top drain — planes 1
# and 2 of §11.1's replacement for the shared mutable model — snapshot
# publication, plane 3, and the log riding behind it (§11.2). The cells'
# owners — the roster's device entries and the harness register beside them —
# live in roster.jl, the device tasks that stage into them in devices.jl;
# what stands further out in the spec — the trace, the snapshot's framework
# status — is deliberately absent (README).
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

"""
The compiled writer (§11.4): one write surface's staged representation, fixed
at a stopped-sim point — a device's compiled at attach against its claim set,
the harness register's at deployment binding and at every roster change
against its *derived* surface, the unclaimed complement (§11.3). The batch is
a positional tuple over the surface, `Union{Nothing, Tᵢ}` per face, isbits and
pointer-free, with `nothing` meaning *not touched this time*, never "reset".
The face-name → position schema, the slot types and the per-position cell
addresses — the compiled scatter's data — live here beside the staging cell
they describe; a roster entry carries one of these per device (roster.jl), and
the harness register carries the one whose shape the framework derives rather
than receives.
"""
struct Writer{B<:Tuple,A<:Tuple}
    faces::Vector{Symbol}    # position → face name: the schema
    types::Vector{Any}       # position → the slot's declared type Tᵢ
    addrs::A                 # position → the slot's cell address
    cell::StagingCell{B}
end

_port_type(::CellAddr{P,K}) where {P,K} = P

function Writer(layout::Layout, faces::Vector{Symbol})
    addrs = Tuple(layout.addr[("", f)] for f in faces)
    types = Any[_port_type(a) for a in addrs]
    B = Tuple{(Union{Nothing,T} for T in types)...}
    Writer{B,typeof(addrs)}(faces, types, addrs, StagingCell{B}(nothing))
end

"""
§11.4's shim, run at staging on the writer's own side: the author's face ⇒
value pairs become a batch — name → position, convert to the slot's declared
type, fill `nothing`. Every diagnostic site follows the compilation here, to
staging: face validity, surface membership and value convertibility are all
static facts of the run, so a face with no position in the schema is discarded
with a warning and an unconvertible value with `EntryTypeMismatch`, while the
rest of the batch stands. Nothing remains at the drain.

The out-of-schema warning discriminates by writer (§11.3, Appendix C): a
device's entry is always `OutOfClaimEntry` — naming the incumbent when the
face is claimed elsewhere — while the harness register's is `ClaimedFaceEntry`
naming the incumbent when a rostered claim covers the face, and
`OutOfClaimEntry` only when the face names no root slot at all. `claimedby` is
the exclusivity index the roster maintains; `device` identifies a device
writer and `nothing` the harness; `site` distinguishes ordinary staging from
an attach's renormalization in the `ClaimedFaceEntry` payload. Returns
`nothing` for a batch with no surviving entry.
"""
function _normalize(w::Writer{B}, pairs, claimedby::Dict{Symbol,String};
                    device::Union{Nothing,String} = nothing,
                    site::String = "at staging") where {B}
    vals = Any[nothing for _ in w.faces]
    for (face, v) in pairs
        s = Symbol(face)
        i = findfirst(==(s), w.faces)
        if i === nothing
            incumbent = get(claimedby, s, nothing)
            if device !== nothing
                held = incumbent === nothing ? "" : " — the face is claimed by $incumbent"
                @warn "OutOfClaimEntry: `$face` is outside $device's claim, whose set is " *
                      "$(_facelist(w.faces))$held; the entry is discarded (§11.3)"
            elseif incumbent !== nothing
                @warn "ClaimedFaceEntry: `$face` is claimed by $incumbent, and the harness " *
                      "register speaks only for unclaimed faces; the entry is discarded " *
                      "$site (§11.3)"
            else
                @warn "OutOfClaimEntry: `$face` names no face on the staging surface — " *
                      "the surface is $(_facelist(w.faces)); the entry is discarded (§11.4)"
            end
            continue
        end
        vals[i] = try
            convert(w.types[i], v)
        catch
            @warn "EntryTypeMismatch: `$face` cannot take $(repr(v))::$(typeof(v)) — " *
                  "its slot declares $(w.types[i]); the entry is discarded (§11.4)"
            continue
        end
    end
    all(isnothing, vals) ? nothing : convert(B, (vals...,))
end

_facelist(faces) = isempty(faces) ? "empty" : "{$(join(faces, ", "))}"

# The one coalescing policy (§11.4): merge, newest wins per face. Untouched
# faces survive; re-staged faces take the newest level — the per-face ZOH.
# Positional, so it compiles straight-line and union-splits. Plain `::Tuple`
# arguments, because a sparse batch's *runtime* type is the narrow tuple of
# what it touched — covariantly inside `B`, but two batches rarely share it.
_merge(pending::Tuple, incoming::Tuple) =
    map((p, i) -> i === nothing ? p : i, pending, incoming)

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

# The drain's application (§11.4): the compiled scatter, position → slot cell,
# statically typed, `nothing` skipped — pure application, no checks.
@generated function _apply!(store, addrs::Tuple, batch::Tuple)
    stmts = [:(batch[$i] === nothing || scatter!(store, addrs[$i], batch[$i]))
             for i in 1:fieldcount(batch)]
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

# --- the diagnostic cell (§11.8) -----------------------------------------------

"""
The one diagnostic value of the prototype's device channel: a datum that could
not be mapped for environmental reasons — a truncated datagram, malformed
JSON, an out-of-range field — reported by the author's own loop body via
`report!(handle, MalformedDatum(cause))` after catching its own parser error
(§11.6: catch, stage nothing, report, continue). The classification is the
author's — only they know their parser — and any *other* exception
propagates to the wrapper as the `DeviceCrash` it is. The payload is the
cause (Appendix C: the exception, or any diagnostic-content value); device
attribution is the cell's, the channel being per-writer.
"""
struct MalformedDatum
    cause::Any
end

"The ring's capacity (§11.8): the channel's normative bound, which *is* the rate limit."
const DIAG_RING = 16

"""
One frame's bounded accumulation (§11.8): a small ring of diagnostic values,
capacity `DIAG_RING`, plus a count of what the ring could not hold — the
per-kind record collapsed to one integer, the closed kind set having exactly
one member here (`MalformedDatum`; the staging kinds warn synchronously as a
stand-in, README). The drop policy is earliest-in-frame retained, excess
becomes counts: the first occurrences are the ones with diagnostic content,
the hundredth is noise the count already reports. Immutable, because the CAS
handoff below needs value semantics: an append builds a new batch, and a
drained batch is frozen — never written again — with the writer allocating
afresh at its next emission, on its own task (§11.8's GC-over-reuse trade).
"""
struct DiagBatch
    ring::Vector{MalformedDatum}
    suppressed::Int
end

"The shared empty sentinel (§11.8): what the drain swaps *in*, and what a quiet frame gets back."
const EMPTY_DIAG = DiagBatch(MalformedDatum[], 0)

"""
The diagnostic cell (§11.8): one per rostered device, single writer — the
device's own task, the same ownership argument as the staging cells: no
locking, no arbitration, no new primitive. The CAS mirrors `_stage!`'s: a
failed replace means the loaded batch was intercepted by the drain, and the
retry re-reads what is pending now — the intercepted entries are already
taken, and the entry in hand merges against the fresh (sentinel) state.
"""
mutable struct DiagCell
    @atomic batch::DiagBatch
end

# The writer's side (§11.8), on the device task: append under the bound, or
# count past it. Reached through `report!(handle, …)` (devices.jl).
function _report!(cell::DiagCell, d::MalformedDatum)
    while true
        cur = @atomic cell.batch
        next = length(cur.ring) < DIAG_RING ?
            DiagBatch(push!(copy(cur.ring), d), cur.suppressed) :
            DiagBatch(cur.ring, cur.suppressed + 1)
        (; success) = @atomicreplace cell.batch cur => next
        success && return nothing
    end
end

# The indivisible take (§11.8): one `atomicswap` per cell at frame top, at the
# same point and under the same argument as the staging drain. A quiet frame
# swaps the sentinel in and gets the sentinel back — no allocation, and no
# load-only code path that goes untested on healthy runs.
_take!(cell::DiagCell) = @atomicswap cell.batch = EMPTY_DIAG

"""
The loop's per-run totals for one writer (§11.8): cumulative per-kind counts
since the run began — one field, the kind set's one member — owned by the
loop, reset at each run's top, and read while stopped by `diagnostics(sim)`
(the stand-in for the published framework status the snapshot does not carry,
README). Mutable in its own holder so the roster entry stays immutable.
"""
mutable struct DiagTotals
    malformed::Int
end

# --- publication (§11.2) -------------------------------------------------------

"""
A snapshot (§11.2): the boundary-consistent signal table — the *whole* table,
root slots riding along as the source cells they are — with `t`, the frame
ordinal and the §12.3 boundary ordinal, the counter-home rule: the index
rides *in* the snapshot, so any holder of one — the log, a post-run
inspector — indexes it without consulting the loop. Built in private memory
by copying the cell buffers, then handed out through one release-store;
nothing reachable from a published snapshot is ever written again, which is
what makes the lock-free read sound. The state stores (`x`, `s`, `m`) are
deliberately not carried (§11.2) and the framework status is absent here
(README). Read it with `port(snap, path, name)`, addressed exactly as the
live table.
"""
struct Snapshot{T,S<:StoreBundle}
    t::T
    frame::Int
    boundary::Int     # the §12.3 published-boundary ordinal; boundary zero = 0
    store::S
    layout::Layout    # build-frozen addressing, shared, never copied
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
