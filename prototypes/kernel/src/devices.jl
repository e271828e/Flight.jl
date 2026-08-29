# The device contract and its task machinery (§11.6, §11.1) with the §12.4
# lifecycle slice: the handle every attached device receives, the authoring
# contract's four functions, the framework wrapper around the author-owned
# loop body, the pre-spawn init bracket and the tail. The control surface here
# is §12.1's stop word, §12.3's counter-plus-condition wait and §12.6's
# lifecycle with §13.5's termination record beside it: pause, pacing and the
# operator interrupt are absent (README). A device
# failure reports as `DeviceCrash` into the device's own diagnostic cell
# (§11.8, §12.4); what the tail alone produces — the join timeout, and
# whatever landed past the final frame top — is folded into the termination
# record and presented through the logging backend, the last snapshot having
# already been published (D-201, D-203).
#
# This file holds the types, the contract and the per-run task mechanics;
# `run!`'s frame anatomy and the Simulation-typed surface (`attach!`,
# `stop!(sim)`) live in sim.jl. The helpers here take `sim` untyped for
# include order only — they run once per run, never inside a frame.

"""
§13.5's termination sources (D-203): the diagnostic convention applied to the
run's outcome — each kind is its identity, its payload plain data — without
joining Appendix C's diagnostic set, the record being outcome, not warning.
`EndTimeReached` carries nothing: the record's own `t` is the fact, and the
configured bound lives with the policy. `ModelRequestedStop` carries the
first named `stop_on` face observed holding, in declaration order.
`ControlRequestedStop` carries its issuer — `:code` from `stop!(sim)`, or the
requesting device's name from `stop!(handle)`; the spec's `:interrupt` arm is
absent with the interrupt's machinery (README). `LoopError` is §13.6's
abnormal entry, `exception` the retained cause (raw: §13.4's `StepError` wrap
and its cursor are absent, README).
"""
abstract type TerminationSource end

struct EndTimeReached <: TerminationSource end

struct ModelRequestedStop <: TerminationSource
    face::Symbol
end

struct ControlRequestedStop <: TerminationSource
    issuer::Union{Symbol,String}
end

struct LoopError <: TerminationSource
    exception::Any
end

"""
One writer's share of the tail residue (§11.8, D-203): what the run's-end
sweep took past the final frame top — the final ring, at most `DIAG_RING`
entries, and the per-kind counts the ring refused. A quiet writer contributes
no record.
"""
struct ResidueRecord
    writer::String
    recent::Vector{DiagValue}
    suppressed::KindCounts
end

"""
§13.5's termination record: the run's *outcome*, where the deployment carries
its policy — so a stopped simulation answers "why did it stop?", and "how did
the stop go?", without its consumer reconstructing either from the clock or
the log stream (D-203). `t` is the final snapshot's boundary time in the
deployment's own scalar (§7.2), `nothing` when no boundary ever ran; `source`
is the typed source above; `residue` is what the run's-end sweep collected —
recorded here and presented through the logging backend, never published
(D-201, D-203). `init!` clears the record with the trajectory.
"""
struct TerminationRecord{T}
    t::Union{Nothing,T}
    source::TerminationSource
    residue::Vector{ResidueRecord}
end

"""
The control surface — §12.1's stop and §12.6's lifecycle. `stop_issuer` is
the control-plane stop word, and it carries its issuer (D-203): set by
`stop!` from any task through a compare-and-swap from `nothing` — the first
writer wins, so the recorded issuer is the request that actually initiated
the tail — consulted for non-empty by the loop at frame top, cleared at the
top of the next run (which is what lets an init bracket's `should_abort`
failure leave a stop *already pending* at the run's start, §12.4) and at
`init!` (a fresh trajectory owes nothing to the last one's stop). `stopped`
is §12.4(1)'s sticky status: set only after the final snapshot is published,
and read by `running(handle)`, which is how loop bodies observe the run's
end.

`lifecycle` is §12.6's five-state machine: `:built`, `:initialized`,
`:running`, and terminally `:stopped` or `:errored` (§13.6). `:running` is
the §11.3 freeze, and it deliberately spans the whole of `run!` — tail
included, the terminal state landing in the outermost `finally` — while
`stopped` flips at tail step (1), so device loops exit while the joins are
still ahead. `termination` is the §13.5 record beside it, written by the
run's own task before the terminal state's release-store; readers reach both
through `lifecycle(sim)`/`termination(sim)`, never the raw fields.

`cond` and `counter` are §12.3's two artifacts: the counter counts *published
boundaries* — grid, `t*`, boundary zero — mirrored under the lock right
behind each release-store of `latest`, in that normative order, so a waiter
observing `counter > last_seen` can never wake onto a stale snapshot. The
counter is monotonic across runs and never re-armed: its absolute value is
nowhere normative, and monotonicity keeps the predicate sound with no
per-run reset.
"""
mutable struct Control
    @atomic stop_issuer::Union{Nothing,Symbol,String}
    @atomic stopped::Bool
    cond::Threads.Condition
    counter::Int
    @atomic lifecycle::Symbol
    termination::Union{Nothing,TerminationRecord}
end
Control() = Control(nothing, true, Threads.Condition(), 0, :built, nothing)

# The stop word's one write path (§12.1, D-203): first CAS from empty wins —
# the same arbitration as the loop reacting to the first holding stop face —
# and a later issuer is dropped, the tail already having its initiator.
_request_stop!(ctl::Control, issuer::Union{Symbol,String}) =
    (@atomicreplace ctl.stop_issuer nothing => issuer; nothing)

"""
The §11.3 freeze, keyed on the lifecycle (§12.6): `attach!`, `detach!` and the
stopped-sim readers are refused exactly while `run!` or `step!` holds the
simulation `:running` — which spans the tail, so a roster change cannot race
the joins. Every other state admits them, `:errored` included: post-mortem
inspection of a terminally stopped simulation is legitimate (§13.6).
"""
assert_stopped(ctl::Control, op::Symbol) =
    (@atomic ctl.lifecycle) === :running ?
    throw(BuildError(ServiceLifecycle(op = op, status = :running))) : nothing

"""
The handle (§11.6): the one object every attached device receives, carrying
exactly the primitive capabilities — read (`latest`, `wait_next_snapshot`),
stage (`stage!`) and control access (`running`, `stop!`) — and deliberately
*not* the `Simulation`: what a device may touch is what the handle holds.
`attach!` constructs it and returns it, and the same handle is what the
wrapper passes to `loop(dev, handle)` on the device's task. It also carries
the attachment's binding, read back by `binding(handle)`: the loop's own
`map_input`/`map_output` calls take it from the handle instead of the device
carrying its configuration (§11.6). `last_seen` is the §12.3 waiter's private
register, refreshed at spawn so a run's first wait observes that run's
boundaries. One unguarded edge (README): staging through a handle whose
device was detached lands in an orphaned cell and is silently lost — handles
are run-scoped task equipment, and guarding would put a roster scan back
into `stage!`.
"""
mutable struct DeviceHandle
    const id::Int
    const who::String
    const b::AbstractBinding
    const writer::Writer
    const plane::DataPlane
    const ctl::Control
    const published::Published
    const diag::DiagCell
    const gatherer::Union{Nothing,ReadGather}   # the compiled reads; nothing without an output side
    last_seen::Int
end

# --- the authoring contract (§11.6) --------------------------------------------

"""
The authoring contract: four functions, one optional, one trait
(`needs_calling_task`, roster.jl). `init!` and `shutdown!` are the per-run
resource bracket, with no-op defaults — a stub device holds nothing — and
`shutdown!` must tolerate a partially initialized device: it is guaranteed
on every exit path, the failed-`init!` bracket included. `unblock!` is
§12.4(3)'s optional hook, default no-op: an override makes the device's own
blocking call return. `loop` is the one mandatory function — the
author-owned task body, looping on `running(handle)` with interruptible
blocking. A device whose loop was never written would otherwise present as
"attached, nothing happens", the hidden-bug class §11.6 tolerates nowhere, so
`attach!` refuses it by name (`check_device`, roster.jl) before it is ever
rostered — `DeviceContractMismatch`, service, fail-fast (Appendix C). This
fallback exists only as the comparison target `which` needs to detect a
device with no method of its own; it is unreachable, `attach!` having already
refused, and its own throw says so.
"""
init!(::AbstractDevice) = nothing
shutdown!(::AbstractDevice) = nothing
unblock!(::AbstractDevice) = nothing
loop(dev::AbstractDevice, handle) =
    throw(InternalInvariant("unreachable: attach! refuses a device with no loop method"))

# --- the handle primitives (§11.6) ---------------------------------------------

"""
    running(handle)

§12.4's exit predicate: true from a run's start until tail step (1) sets the
sticky stopped status — after the final snapshot, so a body observing false
may still read a complete final world. The author's loop obligation is to
check it between blocking points (§11.6).

Every handle primitive that a loop pass touches — this one, `latest`,
`stage!`, `wait_next_snapshot`, `gather`, `report!` — stores the liveness
heartbeat on its way through (§11.8, §12.2): the framework observes activity
without owning the loop body, and there is no separate liveness channel to
remember to feed.
"""
running(h::DeviceHandle) = (_beat!(h.diag); !(@atomic h.ctl.stopped))

"""
    stop!(handle)

Request a control-plane stop (§12.1): sets the stop word from any task, the
device's name riding as its issuer into the termination record (D-203); the
loop observes it at the next frame top, completes that boundary, publishes,
and enters the tail (§12.4). Idempotent — a second request loses the CAS and
changes nothing — and inert while already stopped.
"""
stop!(h::DeviceHandle) = _request_stop!(h.ctl, h.who)

"""
    binding(handle)

The attachment's binding (§11.6), for the author-owned loop idiom —
`stage!(handle, map_input(datum, binding(handle))...)`. The framework never
calls `map_input`/`map_output` itself: they are conventions of the loop
body, and the handle carrying the binding is what keeps the device struct
free of its per-deployment configuration (the binding stays an `attach!`
argument, never a device field).
"""
binding(h::DeviceHandle) = h.b

"""
    latest(handle)

The handle's primitive read (§11.6): acquire-load the most recently published
snapshot — exactly `latest(sim)`, through the capability the handle carries.
"""
latest(h::DeviceHandle) = (_beat!(h.diag); @atomic :acquire h.published.latest)

"""
    stage!(handle, "face" => value, ...)

The handle's stage primitive (§11.6, §11.4): the same compiled shim, CAS
merge and per-face newest-wins as every writer's, against this attachment's
claim set — every check at staging, an out-of-claim face discarded under
`OutOfClaimEntry` naming the incumbent, the rest of the batch standing. From
any task, at any wall-clock moment; the batch lands at the top of the next
frame `run!` advances.
"""
function stage!(h::DeviceHandle, pairs::Pair...)
    _beat!(h.diag)
    batch = _normalize(h.writer, pairs, h.plane.claimedby, h.diag; device = h.who)
    batch === nothing || _stage!(h.writer, batch)
    nothing
end

"""
    gather(handle, snap)

The output side's read (§11.2, §11.6): run the attachment's compiled gather —
`reads(b)`, resolved at attach — over a snapshot, returning the labeled
NamedTuple `map_output` receives. The loop idiom is
`send(dev.socket, map_output(gather(handle, snap), binding(handle)))`, on the
device's own task, against the snapshot §12.3's wait handed it: the compiled
addresses read the frozen store, so no name is resolved per datum and nothing
here touches the running loop. On a handle whose binding declares no output
side the call is a contract misuse, and throws by name.
"""
function gather(h::DeviceHandle, s::Snapshot)
    _beat!(h.diag)
    h.gatherer === nothing && throw(BuildError(
        DeviceContractMismatch(device = h.who, reason = :no_output_side)))
    _gather(h.gatherer, s)
end

"""
    report!(handle, MalformedDatum(cause))

The bad-datum channel (§11.6, §11.8): the single-writer entry point into the
device's own diagnostic cell, from the author's loop body after catching its
own parser error — catch, stage nothing, report, continue. The tolerance is
bounded by the cell (`DIAG_RING` retained values per frame, the excess a
per-kind count), so a peer flooding malformed datagrams costs sixteen
retained values and an integer increment, whatever it does, and no writer can
starve another — the cells are disjoint. This is not a general
user-diagnostics channel: `MalformedDatum` is the one thing the *author* may
put in it, and any exception the loop does *not* classify as a bad datum
propagates to the wrapper as the `DeviceCrash` it is. The loop drains the
cell at frame top, folding it into the published framework status —
device-attributed, delta plus totals (§11.8) — and sweeps it once more at the
run's end for whatever landed past the last frame top.
"""
report!(h::DeviceHandle, d::MalformedDatum) = (_beat!(h.diag); _report!(h.diag, d))

"""
    wait_next_snapshot(handle)

§12.3's next-snapshot wait: block until a boundary this waiter has not seen
is published, or until the run stops, and return the latest snapshot — at
least the boundary whose publication woke the wait, newest-wins beyond it (a
slow consumer skips boundaries and always receives the current world). The
predicate loop is the canonical idiom: the condition carries no facts, and
the two that matter — the counter and the stopped status — are re-tested at
every wake, which is what makes shutdown work: tail step (2) wakes every
waiter and the predicate routes it out. After a stop return the author's
loop re-checks `running(handle)`, exactly as after any blocking call.
"""
function wait_next_snapshot(h::DeviceHandle)
    _beat!(h.diag)
    ctl = h.ctl
    lock(ctl.cond)
    try
        while ctl.counter <= h.last_seen && !(@atomic ctl.stopped)
            wait(ctl.cond)
        end
        h.last_seen = ctl.counter
    finally
        unlock(ctl.cond)
    end
    latest(h)
end

# --- the wrapper and the run's bracket (§11.6, §12.4) --------------------------

# The guarded release: `shutdown!` is guaranteed on every exit path, and a
# throw out of it must not wreck the bracket or the tail around it (§11.6).
function _shutdown!(e::RosterEntry)
    try
        shutdown!(e.dev)
    catch err
        @warn "shutdown! of $(_who(e)) threw; its resources may leak (§11.6)" #=
            =# exception = (err, catch_backtrace())
    end
    nothing
end

"""
The framework wrapper (§11.6): the author owns the loop body, the framework
owns the bracket. Every exit path — voluntary return, a crash, the
stop-drained predicate — runs `shutdown!` and then consults the attachment's
`should_abort`; a crash is caught here and reported as `DeviceCrash` into
the device's own diagnostic cell (§11.8, §12.4(6)) — on the device's task,
the cell's writer — the run continuing with the device's task absent and its
claims held to run end. A `needs_calling_task` device runs this identical
wrapper inline on the calling task: the invocation site is its only
difference (§11.1). Death is marked nowhere beyond the record: the task has
ended, `task_state` says so at the next publication, and the heartbeat goes
stale (§12.2).
"""
function _wrap(e::RosterEntry)
    try
        loop(e.dev, e.handle)
    catch err
        _report!(e.diag, DeviceCrash(err, e.should_abort))
    finally
        _shutdown!(e)
        e.should_abort && stop!(e.handle)
    end
    nothing
end

"""
§12.4's pre-spawn initialization bracket, a step of the shutdown protocol
taken at the top of a run: `init!` once per roster entry, in attachment
order, on the calling task, each call in its own bracket. A throw goes
straight back to `shutdown!` — which is why `shutdown!` owes tolerance of a
partially initialized device — is reported as the ordinary `DeviceCrash`,
and spawns no task: the device is dead from the run's first frame, its
claims persisting to run end and the orphaned root inputs holding their values.
With `should_abort` set the failure requests a stop, already pending when
the loop would start: the run advances zero frames and ends through the same
tail, every remaining entry still getting its `init!`/`shutdown!` pair
uniformly. Returns the live entries, from which §11.1's topology is derived
— derived *after* initialization, never from the roster alone.
"""
function _init_devices!(sim)
    live = RosterEntry[]
    for e in sim.plane.roster
        ok = try
            init!(e.dev)
            true
        catch err
            _shutdown!(e)
            _report!(e.diag, DeviceCrash(err, e.should_abort))   # addressed by the entry:
            e.should_abort && stop!(e.handle)                    # no task holds a handle yet (§12.4)
            false
        end
        ok && push!(live, e)
    end
    live
end

# Spawn the wrapper, one run-scoped task per entry (§11.1): inside `run!`
# only, never at `attach!` — a task exists only once the run it serves does.
# The §12.3 registers are refreshed first, on the calling task.
function _spawn!(entries::Vector{RosterEntry})
    for e in entries
        e.handle.last_seen = e.handle.ctl.counter
    end
    [Threads.@spawn _wrap(e) for e in entries]
end

# Tail steps (1)'s close and (2) (§12.4): the final snapshot is whatever the
# loop last published; the sticky status flips only after it, and the notify
# under the lock wakes every §12.3 waiter, whose predicate routes it out.
# Idempotent, and run even when the loop leaves by a throw — the §13.6 catch
# path is absent (README), but device tasks must never be left parked.
function _finish!(sim)
    ctl = sim.control
    @atomic ctl.stopped = true
    lock(ctl.cond)
    try
        notify(ctl.cond)
    finally
        unlock(ctl.cond)
    end
    nothing
end

"""
Tail steps (3)–(5) (§12.4), on the calling task, after `_finish!` has set
the sticky status and woken the waits: `unblock!` per spawned entry — the
override's own blocking call returns; a throw out of the hook is warned and
the tail proceeds — then the join under one shared `join_timeout` deadline,
D-198's one patience for the whole tail. A task exceeding what remains of
the deadline is abandoned rather than left to hang `run!`, reported by name
under `DeviceJoinTimeout` into the loop's own cell (D-203): the terminal
snapshot precedes the join by construction, so the run's-end sweep — not a
drain — is what collects it, into the termination record and the logging
backend. The calling-task device sits outside the join: nothing can abandon
the task `run!` stands on.
"""
function _tail!(sim, entries::Vector{RosterEntry}, tasks::Vector{Task})
    for e in entries
        try
            unblock!(e.dev)
        catch err
            @warn "unblock! of $(_who(e)) threw; its task can now exit only " *
                  "through the join timeout (§12.4)" exception = (err, catch_backtrace())
        end
    end
    deadline = time() + sim.join_timeout
    for (e, t) in zip(entries, tasks)
        remaining = deadline - time()
        joined = istaskdone(t) ||
            (remaining > 0 &&
             timedwait(() -> istaskdone(t), remaining; pollint = min(0.01, remaining)) === :ok)
        if !joined
            s = latest(sim)                  # after init!, never nothing (§14.5)
            _report!(sim.loop_diag, DeviceJoinTimeout(_who(e), sim.join_timeout,
                                                      Float64(s.t), s.boundary))
        end
    end
    nothing
end

"""
The run's last sweep (§11.8), in `run!`'s outermost `finally`: one more take
per cell, folding what landed after the final frame top — a crash caught on
the way out, a report from a device's exit path, the tail's own
`DeviceJoinTimeout` — into the accounts, so the per-run record is complete
and nothing leaks into the next run's status. The terminal snapshot is
already out, so no snapshot can carry this residue: it is collected into the
returned `ResidueRecord`s — the termination record's third field — and
presented through the logging backend, the record's renderer (D-201, D-203).
The terminal status's account is therefore complete up to its own frame top,
and the tail's remainder is loud *and* recorded, still never published.
"""
function _sweep_tail!(sim)
    plane = sim.plane
    residue = ResidueRecord[]
    for e in plane.roster
        _fold!(e.acct, e.diag)
        _residue!(residue, _who(e), e.acct)
    end
    _fold!(plane.harness_acct, plane.harness_diag)
    _residue!(residue, "harness", plane.harness_acct)
    _fold!(sim.loop_acct, sim.loop_diag)
    _residue!(residue, "loop", sim.loop_acct)
    residue
end

# One writer's take: a quiet account contributes no record; a noisy one hands
# its pending vector over — the account re-arms the shared empty, so the
# record's vector is never written again — and is rendered entry by entry.
function _residue!(out::Vector{ResidueRecord}, who::String, a::WriterAccount)
    isempty(a.recent) && _total(a.suppressed) == 0 && return nothing
    r = ResidueRecord(who, a.recent, a.suppressed)
    push!(out, r)
    for d in r.recent
        @warn "$(nameof(typeof(d))) from $who, past the final snapshot's account: $d (§11.8)"
    end
    s = _total(r.suppressed)
    s > 0 && @warn "$s more suppressed occurrence(s) from $who past the final " *
                   "snapshot's account (§11.8)"
    a.recent = EMPTY_RECENT
    a.suppressed = KindCounts()
    nothing
end
